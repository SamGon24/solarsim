// =============================================================================
//  SolarSim - ThreadPool.cpp
// =============================================================================
#include "ThreadPool.hpp"

#include <algorithm>

namespace solsim {

unsigned ThreadPool::defaultThreadCount() {
    const unsigned hc = std::thread::hardware_concurrency();
    return hc == 0 ? 1u : hc;
}

ThreadPool::ThreadPool(unsigned threads) {
    if (threads == 0) threads = defaultThreadCount();
    const unsigned extra = threads > 1 ? threads - 1 : 0;
    m_workers.reserve(extra);
    for (unsigned i = 0; i < extra; ++i) m_workers.emplace_back([this, i] { workerLoop(i); });
}

ThreadPool::~ThreadPool() {
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_stop = true;
    }
    m_wake.notify_all();
    for (auto& t : m_workers)
        if (t.joinable()) t.join();
}

void ThreadPool::parallelFor(std::size_t count, const std::function<void(std::size_t, std::size_t)>& fn) {
    if (count == 0) return;
    if (m_workers.empty() || count < 64) {
        fn(0, count);
        return;
    }
    // Over-decompose a little for load balancing.
    const std::size_t chunks = std::min<std::size_t>(count, size() * 4);
    {
        std::lock_guard<std::mutex> lock(m_mutex);
        m_job = &fn;
        m_count = count;
        m_chunks = chunks;
        m_nextChunk = 0;
        m_finishedChunks = 0;
        ++m_generation;
    }
    m_wake.notify_all();

    // The calling thread works too.
    for (;;) {
        std::size_t c;
        {
            std::lock_guard<std::mutex> lock(m_mutex);
            if (m_nextChunk >= m_chunks) break;
            c = m_nextChunk++;
        }
        const std::size_t b = c * count / chunks;
        const std::size_t e = (c + 1) * count / chunks;
        fn(b, e);
        std::lock_guard<std::mutex> lock(m_mutex);
        ++m_finishedChunks;
    }
    std::unique_lock<std::mutex> lock(m_mutex);
    m_done.wait(lock, [this] { return m_finishedChunks == m_chunks; });
    m_job = nullptr;
}

void ThreadPool::workerLoop(unsigned /*index*/) {
    unsigned long long seenGeneration = 0;
    for (;;) {
        std::unique_lock<std::mutex> lock(m_mutex);
        m_wake.wait(lock, [&] { return m_stop || (m_generation != seenGeneration && m_job != nullptr); });
        if (m_stop) return;
        seenGeneration = m_generation;
        while (m_job && m_nextChunk < m_chunks) {
            const std::size_t c = m_nextChunk++;
            const auto* job = m_job;
            const std::size_t count = m_count, chunks = m_chunks;
            lock.unlock();
            const std::size_t b = c * count / chunks;
            const std::size_t e = (c + 1) * count / chunks;
            (*job)(b, e);
            lock.lock();
            ++m_finishedChunks;
            if (m_finishedChunks == m_chunks) m_done.notify_all();
        }
    }
}

} // namespace solsim
