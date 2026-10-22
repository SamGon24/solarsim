// =============================================================================
//  SolarSim - SvgWriter.hpp
//  Minimal SVG 1.1 document builder for vector orbit diagrams.
// =============================================================================
#pragma once

#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace solsim::io {

class SvgWriter {
public:
    SvgWriter(double width, double height, const std::string& background = "#05070d");

    void comment(const std::string& text);
    void rect(double x, double y, double w, double h, const std::string& fill, const std::string& stroke = "none",
              double strokeWidth = 0.0);
    void circle(double cx, double cy, double r, const std::string& fill, const std::string& stroke = "none",
                double strokeWidth = 0.0, double opacity = 1.0);
    void line(double x1, double y1, double x2, double y2, const std::string& stroke, double width,
              double opacity = 1.0, const std::string& dash = "");
    void polyline(const std::vector<std::pair<double, double>>& pts, const std::string& stroke, double width,
                  double opacity = 1.0, bool closed = false, const std::string& dash = "");
    void text(double x, double y, const std::string& content, const std::string& fill, double size,
              const std::string& anchor = "start", const std::string& weight = "normal");
    void beginGroup(const std::string& attributes = "");
    void endGroup();
    void raw(const std::string& markup);

    std::string str() const;
    bool save(const std::string& path, std::string& error) const;

    static std::string escape(const std::string& s);

private:
    double m_w, m_h;
    std::ostringstream m_body;
    std::string m_background;
};

} // namespace solsim::io
