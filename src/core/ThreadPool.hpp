// =============================================================================
//  SolarSim - ThreadPool.hpp
//  A small persistent worker pool specialised for data-parallel loops
//  (parallelFor). Workers sleep on a condition variable between jobs, so the
//  per-call overhead is a couple of microseconds - acceptable for force
//  evaluations over thousands of test particles.
// =============================================================================
#pragma once

#include <condition_variable>
#include <cstddef>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace solsim {

class ThreadPool {
public:
    // threads == 0 -> use std::thread::hardware_concurrency()
    explicit ThreadPool(unsigned threads = 0);
    ~ThreadPool();

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

    unsigned size() const { return static_cast<unsigned>(m_workers.size()) + 1; } // + calling thread

    // Calls fn(begin, end) on disjoint chunks covering [0, count). Blocks until
    // all chunks are done. The calling thread participates.
    void parallelFor(std::size_t count, const std::function<void(std::size_t, std::size_t)>& fn);

    static unsigned defaultThreadCount();

private:
    void workerLoop(unsigned index);

    std::vector<std::thread> m_workers;
    std::mutex m_mutex;
    std::condition_variable m_wake;
    std::condition_variable m_done;
    const std::function<void(std::size_t, std::size_t)>* m_job = nullptr;
    std::size_t m_count = 0;
    std::size_t m_chunks = 0;
    std::size_t m_nextChunk = 0;
    std::size_t m_finishedChunks = 0;
    unsigned long long m_generation = 0;
    bool m_stop = false;
};

} // namespace solsim
