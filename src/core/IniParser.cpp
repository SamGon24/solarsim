// =============================================================================
//  SolarSim - IniParser.cpp
// =============================================================================
#include "IniParser.hpp"

#include "StringUtil.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace solsim {

bool IniSection::has(const std::string& key) const {
    for (const auto& kv : entries)
        if (str::iequals(kv.first, key)) return true;
    return false;
}

std::string IniSection::get(const std::string& key, const std::string& fallback) const {
    for (const auto& kv : entries)
        if (str::iequals(kv.first, key)) return kv.second;
    return fallback;
}

bool IniSection::getDoubleChecked(const std::string& key, double& out) const {
    if (!has(key)) return false;
    return str::parseDouble(get(key), out);
}

double IniSection::getDouble(const std::string& key, double fallback) const {
    if (!has(key)) return fallback;
    double v = 0.0;
    if (!str::parseDouble(get(key), v))
        throw std::runtime_error("section [" + name + "] (line " + std::to_string(line) + "): key '" + key +
                                 "' is not a number: '" + get(key) + "'");
    return v;
}

bool IniDocument::loadFile(const std::string& path, std::string& error) {
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        error = "cannot open '" + path + "'";
        return false;
    }
    std::ostringstream ss;
    ss << in.rdbuf();
    return loadString(ss.str(), error, path);
}

bool IniDocument::loadString(const std::string& text, std::string& error, const std::string& sourceName) {
    m_sections.clear();
    m_source = sourceName;
    std::istringstream in(text);
    std::string raw;
    int lineNo = 0;
    IniSection* current = nullptr;
    while (std::getline(in, raw)) {
        ++lineNo;
        // Strip UTF-8 BOM on first line.
        if (lineNo == 1 && raw.size() >= 3 && static_cast<unsigned char>(raw[0]) == 0xEF &&
            static_cast<unsigned char>(raw[1]) == 0xBB && static_cast<unsigned char>(raw[2]) == 0xBF)
            raw = raw.substr(3);
        std::string line = str::trim(raw);
        if (line.empty() || line[0] == ';' || line[0] == '#') continue;
        if (line.front() == '[') {
            if (line.back() != ']') {
                error = sourceName + ":" + std::to_string(lineNo) + ": unterminated section header";
                return false;
            }
            IniSection sec;
            sec.name = str::toLower(str::trim(line.substr(1, line.size() - 2)));
            sec.line = lineNo;
            m_sections.push_back(sec);
            current = &m_sections.back();
            continue;
        }
        const auto eq = line.find('=');
        if (eq == std::string::npos) {
            error = sourceName + ":" + std::to_string(lineNo) + ": expected 'key = value'";
            return false;
        }
        if (!current) {
            error = sourceName + ":" + std::to_string(lineNo) + ": key outside of any [section]";
            return false;
        }
        std::string key = str::toLower(str::trim(line.substr(0, eq)));
        std::string value = str::trim(line.substr(eq + 1));
        // Inline comments must be preceded by whitespace (keeps "a;b" and "#RRGGBB" values intact).
        for (std::size_t i = 1; i < value.size(); ++i) {
            if ((value[i] == ';' || value[i] == '#') && (value[i - 1] == ' ' || value[i - 1] == '\t')) {
                value = value.substr(0, i);
                break;
            }
        }
        value = str::trim(value);
        if (value.size() >= 2 && value.front() == '"' && value.back() == '"') value = value.substr(1, value.size() - 2);
        current->entries.emplace_back(key, value);
    }
    return true;
}

std::vector<const IniSection*> IniDocument::sectionsNamed(const std::string& name) const {
    std::vector<const IniSection*> out;
    for (const auto& s : m_sections)
        if (str::iequals(s.name, name)) out.push_back(&s);
    return out;
}

} // namespace solsim
