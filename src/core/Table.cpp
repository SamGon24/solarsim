// =============================================================================
//  SolarSim - Table.cpp
// =============================================================================
#include "Table.hpp"

#include "StringUtil.hpp"

#include <algorithm>
#include <sstream>

namespace solsim {

namespace {
TextTable::Style g_defaultStyle = TextTable::Style::Unicode;
}

void TextTable::setDefaultStyle(Style s) { g_defaultStyle = s; }
TextTable::Style TextTable::defaultStyle() { return g_defaultStyle; }

TextTable::TextTable(std::vector<std::string> headers)
    : m_headers(std::move(headers)), m_align(m_headers.size(), Align::Right) {
    if (!m_align.empty()) m_align[0] = Align::Left;
}

void TextTable::setAlign(std::size_t column, Align a) {
    if (column < m_align.size()) m_align[column] = a;
}

void TextTable::setAllAlign(Align a) { std::fill(m_align.begin(), m_align.end(), a); }

void TextTable::addRow(std::vector<std::string> cells) {
    cells.resize(m_headers.size());
    m_rows.push_back(std::move(cells));
}

void TextTable::addSeparator() { m_rows.emplace_back(); }

std::string TextTable::cell(const std::string& text, std::size_t width, Align a) const {
    const std::size_t len = solsim::str::utf8Length(text);
    if (len >= width) return text;
    const std::size_t pad = width - len;
    switch (a) {
        case Align::Left: return text + std::string(pad, ' ');
        case Align::Right: return std::string(pad, ' ') + text;
        case Align::Center: return std::string(pad / 2, ' ') + text + std::string(pad - pad / 2, ' ');
    }
    return text;
}

void TextTable::print(std::ostream& os, Style style) const {
    const std::size_t nc = m_headers.size();
    std::vector<std::size_t> w(nc, 0);
    for (std::size_t c = 0; c < nc; ++c) w[c] = solsim::str::utf8Length(m_headers[c]);
    for (const auto& r : m_rows)
        for (std::size_t c = 0; c < r.size() && c < nc; ++c) w[c] = std::max(w[c], solsim::str::utf8Length(r[c]));

    const bool uni = style == Style::Unicode;
    const bool md = style == Style::Markdown;
    const std::string V = uni ? "\xE2\x94\x82" : "|";
    const std::string H = uni ? "\xE2\x94\x80" : "-";

    auto rule = [&](const char* l, const char* m, const char* r) {
        if (md) return;
        std::string line = uni ? l : "+";
        for (std::size_t c = 0; c < nc; ++c) {
            for (std::size_t k = 0; k < w[c] + 2; ++k) line += H;
            line += (c + 1 < nc) ? (uni ? m : "+") : (uni ? r : "+");
        }
        os << line << "\n";
    };

    if (!m_title.empty()) os << m_title << "\n";
    rule("\xE2\x94\x8C", "\xE2\x94\xAC", "\xE2\x94\x90");
    os << V;
    for (std::size_t c = 0; c < nc; ++c) os << ' ' << cell(m_headers[c], w[c], Align::Center) << ' ' << V;
    os << "\n";
    if (md) {
        os << "|";
        for (std::size_t c = 0; c < nc; ++c) {
            std::string dashes(w[c] + 2, '-');
            if (m_align[c] == Align::Right) dashes.back() = ':';
            os << dashes << "|";
        }
        os << "\n";
    } else {
        rule("\xE2\x94\x9C", "\xE2\x94\xBC", "\xE2\x94\xA4");
    }
    for (const auto& r : m_rows) {
        if (r.empty()) {
            rule("\xE2\x94\x9C", "\xE2\x94\xBC", "\xE2\x94\xA4");
            continue;
        }
        os << V;
        for (std::size_t c = 0; c < nc; ++c) os << ' ' << cell(r[c], w[c], m_align[c]) << ' ' << V;
        os << "\n";
    }
    rule("\xE2\x94\x94", "\xE2\x94\xB4", "\xE2\x94\x98");
}

std::string TextTable::str(Style style) const {
    std::ostringstream os;
    print(os, style);
    return os.str();
}

} // namespace solsim
