// =============================================================================
//  SolarSim - Terminal.hpp
//  Thin platform layer for the interactive terminal viewer:
//    * raw (non-canonical, no echo) keyboard input with escape-sequence decoding
//    * terminal size query
//    * enabling ANSI/VT processing and UTF-8 output on Windows 10+
//    * Ctrl+C handling that always restores the terminal
//  POSIX (Linux/macOS/BSD) uses termios + select(); Windows uses the console API.
// =============================================================================
#pragma once

#include <string>

namespace solsim::platform {

enum class Key {
    None, Char, Up, Down, Left, Right, PageUp, PageDown, Home, End, Escape, Enter, Tab, Backspace
};

struct KeyEvent {
    Key key = Key::None;
    char ch = 0; // for Key::Char
};

struct TerminalSize {
    int columns = 100;
    int rows = 32;
};

// Called once at program start: UTF-8 console + VT sequences on Windows.
// Returns true if ANSI escape sequences are believed to work.
bool initConsole();

bool stdoutIsTerminal();
bool stdinIsTerminal();
bool stderrIsTerminal();
TerminalSize terminalSize();
// Heuristic: does the terminal advertise 24-bit colour (COLORTERM=truecolor)?
bool supportsTrueColor();

// RAII guard: raw mode + alternate screen + hidden cursor while alive.
class RawTerminal {
public:
    RawTerminal();
    ~RawTerminal();
    RawTerminal(const RawTerminal&) = delete;
    RawTerminal& operator=(const RawTerminal&) = delete;

    bool active() const { return m_active; }
    // Non-blocking: returns Key::None when no key is pending.
    KeyEvent poll();

private:
    bool m_active = false;
};

// Set by the SIGINT / console control handler; the main loop should exit.
bool interruptRequested();
void installInterruptHandler();

// Writes the whole buffer to stdout (handles partial writes).
void writeStdout(const std::string& data);

// Sleeps for the given number of milliseconds.
void sleepMs(int ms);

// Creates a directory (and parents). Returns true if it exists afterwards.
bool makeDirectories(const std::string& path);

} // namespace solsim::platform
