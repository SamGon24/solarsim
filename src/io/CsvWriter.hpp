// =============================================================================
//  SolarSim - CsvWriter.hpp
//  RFC 4180 CSV output (quotes fields containing separators/quotes/newlines).
// =============================================================================
#pragma once

#include <fstream>
#include <string>
#include <vector>

namespace solsim::io {

class CsvWriter {
public:
    CsvWriter() = default;
    explicit CsvWriter(const std::string& path) { open(path); }

    bool open(const std::string& path);
    bool isOpen() const { return m_out.is_open() && m_out.good(); }
    void header(const std::vector<std::string>& cols) { row(cols); }
    void row(const std::vector<std::string>& cells);
    std::size_t rows() const { return m_rows; }
    void close() { m_out.close(); }

    static std::string escape(const std::string& field);
    static std::string num(double v, int precision = 10);

private:
    std::ofstream m_out;
    std::size_t m_rows = 0;
};

} // namespace solsim::io
