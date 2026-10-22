// =============================================================================
//  SolarSim - Table.hpp
//  Pretty printer for column aligned text tables (box drawing or plain ASCII).
// =============================================================================
#pragma once

#include <ostream>
#include <string>
#include <vector>

namespace solsim {

class TextTable {
public:
    enum class Align { Left, Right, Center };
    enum class Style { Unicode, Ascii, Markdown };

    explicit TextTable(std::vector<std::string> headers);

    void setAlign(std::size_t column, Align a);
    void setAllAlign(Align a);
    void addRow(std::vector<std::string> cells);
    void addSeparator();
    void setTitle(const std::string& title) { m_title = title; }

    void print(std::ostream& os, Style style = Style::Unicode) const;
    std::string str(Style style = Style::Unicode) const;

    std::size_t rowCount() const { return m_rows.size(); }

    // Global default (set from --ascii); Unicode box drawing otherwise.
    static void setDefaultStyle(Style s);
    static Style defaultStyle();

private:
    std::string cell(const std::string& text, std::size_t width, Align a) const;

    std::string m_title;
    std::vector<std::string> m_headers;
    std::vector<Align> m_align;
    std::vector<std::vector<std::string>> m_rows; // empty row == separator
};

} // namespace solsim
