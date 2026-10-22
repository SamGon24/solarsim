// =============================================================================
//  SolarSim - Terminal.cpp
// =============================================================================
#include "Terminal.hpp"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <thread>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <conio.h>
#include <direct.h>
#include <io.h>
#include <windows.h>
#ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
#endif
#else
#include <cerrno>
#include <sys/ioctl.h>
#include <sys/select.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <termios.h>
#include <unistd.h>
#endif

namespace solsim::platform {

namespace {
std::atomic<bool> g_interrupted{false};
bool g_vtEnabled = true;
} // namespace

bool interruptRequested() { return g_interrupted.load(); }

void writeStdout(const std::string& data) {
    std::size_t off = 0;
    while (off < data.size()) {
        const std::size_t n = std::fwrite(data.data() + off, 1, data.size() - off, stdout);
        if (n == 0) break;
        off += n;
    }
    std::fflush(stdout);
}

void sleepMs(int ms) {
    if (ms > 0) std::this_thread::sleep_for(std::chrono::milliseconds(ms));
}

bool supportsTrueColor() {
    const char* ct = std::getenv("COLORTERM");
    if (ct && (std::strstr(ct, "truecolor") || std::strstr(ct, "24bit"))) return true;
#if defined(_WIN32)
    return g_vtEnabled; // Windows 10 1703+ console supports 24-bit colour
#else
    const char* term = std::getenv("TERM");
    if (term && (std::strstr(term, "kitty") || std::strstr(term, "alacritty") || std::strstr(term, "wezterm") ||
                 std::strstr(term, "direct")))
        return true;
    const char* prog = std::getenv("TERM_PROGRAM");
    if (prog && (std::strstr(prog, "iTerm") || std::strstr(prog, "vscode") || std::strstr(prog, "WezTerm")))
        return true;
    return false;
#endif
}

#if defined(_WIN32)
// =============================================================================
//  Windows implementation
// =============================================================================
namespace {
DWORD g_originalInMode = 0;
DWORD g_originalOutMode = 0;
bool g_haveOriginal = false;

BOOL WINAPI consoleHandler(DWORD type) {
    if (type == CTRL_C_EVENT || type == CTRL_BREAK_EVENT || type == CTRL_CLOSE_EVENT) {
        g_interrupted = true;
        return TRUE;
    }
    return FALSE;
}
} // namespace

bool initConsole() {
    SetConsoleOutputCP(CP_UTF8);
    HANDLE out = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD mode = 0;
    if (out == INVALID_HANDLE_VALUE || !GetConsoleMode(out, &mode)) {
        g_vtEnabled = false;
        return false;
    }
    g_originalOutMode = mode;
    g_vtEnabled = SetConsoleMode(out, mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING) != 0;
    return g_vtEnabled;
}

bool stdoutIsTerminal() { return _isatty(_fileno(stdout)) != 0; }
bool stdinIsTerminal() { return _isatty(_fileno(stdin)) != 0; }
bool stderrIsTerminal() { return _isatty(_fileno(stderr)) != 0; }

TerminalSize terminalSize() {
    TerminalSize ts;
    CONSOLE_SCREEN_BUFFER_INFO info;
    if (GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info)) {
        ts.columns = info.srWindow.Right - info.srWindow.Left + 1;
        ts.rows = info.srWindow.Bottom - info.srWindow.Top + 1;
    }
    return ts;
}

void installInterruptHandler() { SetConsoleCtrlHandler(consoleHandler, TRUE); }

RawTerminal::RawTerminal() {
    HANDLE in = GetStdHandle(STD_INPUT_HANDLE);
    DWORD mode = 0;
    if (GetConsoleMode(in, &mode)) {
        g_originalInMode = mode;
        g_haveOriginal = true;
        SetConsoleMode(in, mode & ~(ENABLE_LINE_INPUT | ENABLE_ECHO_INPUT));
    }
    writeStdout("\x1b[?1049h\x1b[?25l\x1b[2J");
    m_active = true;
}

RawTerminal::~RawTerminal() {
    writeStdout("\x1b[0m\x1b[?25h\x1b[?1049l");
    if (g_haveOriginal) SetConsoleMode(GetStdHandle(STD_INPUT_HANDLE), g_originalInMode);
}

KeyEvent RawTerminal::poll() {
    KeyEvent ev;
    if (!_kbhit()) return ev;
    const int c = _getch();
    if (c == 0 || c == 224) {
        const int c2 = _getch();
        switch (c2) {
            case 72: ev.key = Key::Up; break;
            case 80: ev.key = Key::Down; break;
            case 75: ev.key = Key::Left; break;
            case 77: ev.key = Key::Right; break;
            case 73: ev.key = Key::PageUp; break;
            case 81: ev.key = Key::PageDown; break;
            case 71: ev.key = Key::Home; break;
            case 79: ev.key = Key::End; break;
            default: ev.key = Key::None; break;
        }
        return ev;
    }
    switch (c) {
        case 27: ev.key = Key::Escape; break;
        case 13: ev.key = Key::Enter; break;
        case 9: ev.key = Key::Tab; break;
        case 8: ev.key = Key::Backspace; break;
        case 3: g_interrupted = true; break;
        default:
            ev.key = Key::Char;
            ev.ch = static_cast<char>(c);
    }
    return ev;
}

bool makeDirectories(const std::string& path) {
    std::string cur;
    for (std::size_t i = 0; i < path.size(); ++i) {
        const char c = path[i];
        cur.push_back(c);
        if ((c == '/' || c == '\\') && cur.size() > 1 && cur[cur.size() - 2] != ':') _mkdir(cur.c_str());
    }
    if (!path.empty()) _mkdir(path.c_str());
    DWORD attr = GetFileAttributesA(path.c_str());
    return attr != INVALID_FILE_ATTRIBUTES && (attr & FILE_ATTRIBUTE_DIRECTORY);
}

#else
// =============================================================================
//  POSIX implementation
// =============================================================================
namespace {
struct termios g_original;
bool g_haveOriginal = false;

void onSignal(int) { g_interrupted = true; }
} // namespace

bool initConsole() {
    const char* term = std::getenv("TERM");
    g_vtEnabled = !(term && std::strcmp(term, "dumb") == 0);
    return g_vtEnabled;
}

bool stdoutIsTerminal() { return isatty(STDOUT_FILENO) != 0; }
bool stdinIsTerminal() { return isatty(STDIN_FILENO) != 0; }
bool stderrIsTerminal() { return isatty(STDERR_FILENO) != 0; }

TerminalSize terminalSize() {
    TerminalSize ts;
    struct winsize w;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0 && w.ws_col > 0 && w.ws_row > 0) {
        ts.columns = w.ws_col;
        ts.rows = w.ws_row;
        return ts;
    }
    if (const char* c = std::getenv("COLUMNS")) ts.columns = std::max(20, std::atoi(c));
    if (const char* l = std::getenv("LINES")) ts.rows = std::max(10, std::atoi(l));
    return ts;
}

void installInterruptHandler() {
    std::signal(SIGINT, onSignal);
    std::signal(SIGTERM, onSignal);
}

RawTerminal::RawTerminal() {
    if (isatty(STDIN_FILENO) && tcgetattr(STDIN_FILENO, &g_original) == 0) {
        g_haveOriginal = true;
        struct termios raw = g_original;
        raw.c_lflag &= static_cast<tcflag_t>(~(ICANON | ECHO));
        raw.c_cc[VMIN] = 0;
        raw.c_cc[VTIME] = 0;
        tcsetattr(STDIN_FILENO, TCSANOW, &raw);
    }
    writeStdout("\x1b[?1049h\x1b[?25l\x1b[2J");
    m_active = true;
}

RawTerminal::~RawTerminal() {
    writeStdout("\x1b[0m\x1b[?25h\x1b[?1049l");
    if (g_haveOriginal) tcsetattr(STDIN_FILENO, TCSANOW, &g_original);
}

namespace {
bool readByte(unsigned char& c, int timeoutMs) {
    fd_set set;
    FD_ZERO(&set);
    FD_SET(STDIN_FILENO, &set);
    struct timeval tv;
    tv.tv_sec = 0;
    tv.tv_usec = timeoutMs * 1000;
    if (select(STDIN_FILENO + 1, &set, nullptr, nullptr, &tv) <= 0) return false;
    return read(STDIN_FILENO, &c, 1) == 1;
}
} // namespace

KeyEvent RawTerminal::poll() {
    KeyEvent ev;
    unsigned char c = 0;
    if (!readByte(c, 0)) return ev;
    if (c == 27) {
        unsigned char c1 = 0, c2 = 0;
        if (!readByte(c1, 15)) { ev.key = Key::Escape; return ev; }
        if (c1 != '[' && c1 != 'O') { ev.key = Key::Escape; return ev; }
        if (!readByte(c2, 15)) { ev.key = Key::Escape; return ev; }
        switch (c2) {
            case 'A': ev.key = Key::Up; return ev;
            case 'B': ev.key = Key::Down; return ev;
            case 'C': ev.key = Key::Right; return ev;
            case 'D': ev.key = Key::Left; return ev;
            case 'H': ev.key = Key::Home; return ev;
            case 'F': ev.key = Key::End; return ev;
            default: break;
        }
        if (c2 >= '0' && c2 <= '9') {
            unsigned char t = 0;
            readByte(t, 15); // trailing '~'
            if (c2 == '5') ev.key = Key::PageUp;
            else if (c2 == '6') ev.key = Key::PageDown;
            else if (c2 == '1' || c2 == '7') ev.key = Key::Home;
            else if (c2 == '4' || c2 == '8') ev.key = Key::End;
        }
        return ev;
    }
    if (c == '\n' || c == '\r') ev.key = Key::Enter;
    else if (c == '\t') ev.key = Key::Tab;
    else if (c == 127 || c == 8) ev.key = Key::Backspace;
    else if (c == 3) g_interrupted = true;
    else {
        ev.key = Key::Char;
        ev.ch = static_cast<char>(c);
    }
    return ev;
}

bool makeDirectories(const std::string& path) {
    std::string cur;
    for (std::size_t i = 0; i < path.size(); ++i) {
        cur.push_back(path[i]);
        if (path[i] == '/' && cur.size() > 1) mkdir(cur.c_str(), 0755);
    }
    if (!path.empty()) mkdir(path.c_str(), 0755);
    struct stat st;
    return stat(path.c_str(), &st) == 0 && S_ISDIR(st.st_mode);
}
#endif

} // namespace solsim::platform
