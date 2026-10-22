// =============================================================================
//  SolarSim - ArgParser.cpp
// =============================================================================
#include "ArgParser.hpp"

#include "StringUtil.hpp"

#include <algorithm>
#include <sstream>

namespace solsim {

ArgParser::ArgParser(std::string programName) : m_program(std::move(programName)) {}

void ArgParser::addCommand(const std::string& name, const std::string& summary, const std::string& usage) {
    m_commands.push_back({name, summary, usage});
}

void ArgParser::addOption(const std::string& name, const std::string& valueName, const std::string& help,
                          const std::string& group, const std::string& defaultText, char shortName) {
    OptionSpec s;
    s.name = name;
    s.valueName = valueName.empty() ? "VALUE" : valueName;
    s.help = help;
    s.group = group;
    s.defaultText = defaultText;
    s.shortName = shortName;
    s.isFlag = false;
    m_options.push_back(s);
}

void ArgParser::addFlag(const std::string& name, const std::string& help, const std::string& group,
                        char shortName) {
    OptionSpec s;
    s.name = name;
    s.help = help;
    s.group = group;
    s.shortName = shortName;
    s.isFlag = true;
    m_options.push_back(s);
}

const ArgParser::OptionSpec* ArgParser::find(const std::string& name) const {
    for (const auto& o : m_options)
        if (o.name == name) return &o;
    return nullptr;
}

const ArgParser::OptionSpec* ArgParser::findShort(char c) const {
    for (const auto& o : m_options)
        if (o.shortName == c) return &o;
    return nullptr;
}

bool ArgParser::commandKnown(const std::string& name) const {
    return std::any_of(m_commands.begin(), m_commands.end(),
                       [&](const CommandSpec& c) { return c.name == name; });
}

std::string ArgParser::suggest(const std::string& unknown) const {
    std::size_t best = 1000;
    std::string bestName;
    for (const auto& o : m_options) {
        const std::size_t d = str::editDistance(unknown, o.name);
        if (d < best) { best = d; bestName = o.name; }
    }
    if (!bestName.empty() && best <= std::max<std::size_t>(2, unknown.size() / 3))
        return " (did you mean --" + bestName + "?)";
    return "";
}

void ArgParser::parse(int argc, const char* const* argv) {
    m_command.clear();
    m_positionals.clear();
    m_values.clear();

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i] ? argv[i] : "";
        if (arg == "--") {
            for (int j = i + 1; j < argc; ++j) m_positionals.emplace_back(argv[j]);
            break;
        }
        if (str::startsWith(arg, "--")) {
            std::string name = arg.substr(2);
            std::string value;
            bool hasInline = false;
            const auto eq = name.find('=');
            if (eq != std::string::npos) {
                value = name.substr(eq + 1);
                name = name.substr(0, eq);
                hasInline = true;
            }
            const OptionSpec* spec = find(name);
            if (!spec) throw UsageError("unknown option --" + name + suggest(name));
            if (spec->isFlag) {
                if (hasInline) {
                    bool b = false;
                    if (!str::parseBool(value, b))
                        throw UsageError("flag --" + name + " expects true/false, got '" + value + "'");
                    if (b) m_values[name] = "1"; else m_values.erase(name);
                } else {
                    m_values[name] = "1";
                }
                continue;
            }
            if (!hasInline) {
                if (i + 1 >= argc) throw UsageError("option --" + name + " requires a value <" + spec->valueName + ">");
                value = argv[++i];
            }
            m_values[name] = value;
            continue;
        }
        if (arg.size() >= 2 && arg[0] == '-' && !(arg[1] >= '0' && arg[1] <= '9') && arg[1] != '.') {
            // Short option(s): -h, -v, or -o value
            for (std::size_t k = 1; k < arg.size(); ++k) {
                const OptionSpec* spec = findShort(arg[k]);
                if (!spec) throw UsageError(std::string("unknown option -") + arg[k]);
                if (spec->isFlag) {
                    m_values[spec->name] = "1";
                } else {
                    std::string value = arg.substr(k + 1);
                    if (value.empty()) {
                        if (i + 1 >= argc) throw UsageError("option -" + std::string(1, arg[k]) + " requires a value");
                        value = argv[++i];
                    }
                    m_values[spec->name] = value;
                    break;
                }
            }
            continue;
        }
        if (m_command.empty()) m_command = arg;
        else m_positionals.push_back(arg);
    }
}

bool ArgParser::has(const std::string& name) const { return m_values.count(name) != 0; }

std::string ArgParser::getString(const std::string& name, const std::string& fallback) const {
    const auto it = m_values.find(name);
    return it == m_values.end() ? fallback : it->second;
}

double ArgParser::getDouble(const std::string& name, double fallback) const {
    const auto it = m_values.find(name);
    if (it == m_values.end()) return fallback;
    double v = 0.0;
    if (!str::parseDouble(it->second, v))
        throw UsageError("option --" + name + " expects a number, got '" + it->second + "'");
    return v;
}

long long ArgParser::getInt(const std::string& name, long long fallback) const {
    const auto it = m_values.find(name);
    if (it == m_values.end()) return fallback;
    long long v = 0;
    if (!str::parseInt(it->second, v))
        throw UsageError("option --" + name + " expects an integer, got '" + it->second + "'");
    return v;
}

std::vector<std::string> ArgParser::getList(const std::string& name) const {
    std::vector<std::string> out;
    const auto it = m_values.find(name);
    if (it == m_values.end()) return out;
    for (const auto& p : str::split(it->second, ',')) {
        const std::string t = str::trim(p);
        if (!t.empty()) out.push_back(t);
    }
    return out;
}

std::string ArgParser::helpText() const {
    std::ostringstream os;
    os << "Usage: " << m_program << " <command> [options]\n\nCommands:\n";
    std::size_t w = 0;
    for (const auto& c : m_commands) w = std::max(w, c.name.size());
    for (const auto& c : m_commands) os << "  " << str::padRight(c.name, w + 2) << c.summary << "\n";

    std::vector<std::string> groups;
    for (const auto& o : m_options)
        if (std::find(groups.begin(), groups.end(), o.group) == groups.end()) groups.push_back(o.group);

    for (const auto& g : groups) {
        os << "\n" << g << ":\n";
        for (const auto& o : m_options) {
            if (o.group != g) continue;
            std::string left = "  ";
            left += o.shortName ? std::string("-") + o.shortName + ", " : "    ";
            left += "--" + o.name;
            if (!o.isFlag) left += " <" + o.valueName + ">";
            os << str::padRight(left, 34) << o.help;
            if (!o.defaultText.empty()) os << " [default: " << o.defaultText << "]";
            os << "\n";
        }
    }
    os << "\nRun '" << m_program << " help <command>' for command specific usage.\n";
    return os.str();
}

std::string ArgParser::commandHelp(const std::string& command) const {
    for (const auto& c : m_commands) {
        if (c.name != command) continue;
        std::ostringstream os;
        os << m_program << " " << c.name << " - " << c.summary << "\n";
        if (!c.usage.empty()) os << "\n" << c.usage << "\n";
        return os.str();
    }
    return "Unknown command '" + command + "'.\n";
}

} // namespace solsim
