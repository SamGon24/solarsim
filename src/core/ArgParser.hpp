// =============================================================================
//  SolarSim - ArgParser.hpp
//  Command line parser supporting:
//      program <command> [positionals...] [--opt value] [--opt=value] [--flag] [-x]
//  Options are declared up front (with help text and groups) so the parser can
//  generate --help output and suggest corrections for misspelled options.
// =============================================================================
#pragma once

#include <map>
#include <stdexcept>
#include <string>
#include <vector>

namespace solsim {

class UsageError : public std::runtime_error {
public:
    explicit UsageError(const std::string& what) : std::runtime_error(what) {}
};

class ArgParser {
public:
    struct OptionSpec {
        std::string name;        // long name without leading dashes
        char shortName = 0;      // optional single-letter alias
        std::string valueName;   // e.g. "DAYS"; empty for flags
        std::string help;
        std::string group;
        std::string defaultText; // shown in help only
        bool isFlag = false;
    };

    struct CommandSpec {
        std::string name;
        std::string summary;
        std::string usage;
    };

    explicit ArgParser(std::string programName);

    void addCommand(const std::string& name, const std::string& summary, const std::string& usage = "");
    void addOption(const std::string& name, const std::string& valueName, const std::string& help,
                   const std::string& group, const std::string& defaultText = "", char shortName = 0);
    void addFlag(const std::string& name, const std::string& help, const std::string& group,
                 char shortName = 0);

    // Throws UsageError on malformed input or unknown options.
    void parse(int argc, const char* const* argv);

    const std::string& command() const { return m_command; }
    const std::vector<std::string>& positionals() const { return m_positionals; }
    bool commandKnown(const std::string& name) const;

    bool has(const std::string& name) const;
    bool flag(const std::string& name) const { return has(name); }
    std::string getString(const std::string& name, const std::string& fallback = "") const;
    double getDouble(const std::string& name, double fallback) const;
    long long getInt(const std::string& name, long long fallback) const;
    std::vector<std::string> getList(const std::string& name) const; // comma separated

    std::string helpText() const;
    std::string commandHelp(const std::string& command) const;
    const std::vector<CommandSpec>& commands() const { return m_commands; }

private:
    const OptionSpec* find(const std::string& name) const;
    const OptionSpec* findShort(char c) const;
    std::string suggest(const std::string& unknown) const;

    std::string m_program;
    std::vector<OptionSpec> m_options;
    std::vector<CommandSpec> m_commands;
    std::string m_command;
    std::vector<std::string> m_positionals;
    std::map<std::string, std::string> m_values;
};

} // namespace solsim
