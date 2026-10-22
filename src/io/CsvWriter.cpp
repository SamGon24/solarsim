// =============================================================================
//  SolarSim - CsvWriter.cpp
// =============================================================================
#include "CsvWriter.hpp"

#include "../core/StringUtil.hpp"

namespace solsim::io {

bool CsvWriter::open(const std::string& path) {
    m_out.open(path, std::ios::binary | std::ios::trunc);
    m_rows = 0;
    return isOpen();
}

std::string CsvWriter::escape(const std::string& f) {
    if (f.find_first_of(",\"\r\n") == std::string::npos) return f;
    return "\"" + str::replaceAll(f, "\"", "\"\"") + "\"";
}

std::string CsvWriter::num(double v, int precision) { return str::format("%.*g", precision, v); }

void CsvWriter::row(const std::vector<std::string>& cells) {
    for (std::size_t i = 0; i < cells.size(); ++i) {
        if (i) m_out << ',';
        m_out << escape(cells[i]);
    }
    m_out << "\r\n";
    ++m_rows;
}

} // namespace solsim::io
