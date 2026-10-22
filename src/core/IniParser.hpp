// =============================================================================
//  SolarSim - IniParser.hpp
//  Parser for a simple INI dialect used by body definition files:
//
//      ; comment              # also a comment
//      [body]                 <- sections may repeat; order is preserved
//      name = Halley
//      e    = 0.96714         ; trailing comments allowed after whitespace
//
// =============================================================================
#pragma once

#include <string>
#include <utility>
#include <vector>

namespace solsim {

struct IniSection {
    std::string name;
    int line = 0;
    std::vector<std::pair<std::string, std::string>> entries;

    bool has(const std::string& key) const;
    std::string get(const std::string& key, const std::string& fallback = "") const;
    // Throws std::runtime_error with file/line context when the value is not numeric.
    double getDouble(const std::string& key, double fallback) const;
    bool getDoubleChecked(const std::string& key, double& out) const;
};

class IniDocument {
public:
    // Returns false and fills `error` on I/O or syntax errors.
    bool loadFile(const std::string& path, std::string& error);
    bool loadString(const std::string& text, std::string& error, const std::string& sourceName = "<string>");

    const std::vector<IniSection>& sections() const { return m_sections; }
    std::vector<const IniSection*> sectionsNamed(const std::string& name) const;
    const std::string& source() const { return m_source; }

private:
    std::vector<IniSection> m_sections;
    std::string m_source;
};

} // namespace solsim
