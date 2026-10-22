// =============================================================================
//  SolarSim - SvgWriter.cpp
// =============================================================================
#include "SvgWriter.hpp"

#include "../core/StringUtil.hpp"

#include <fstream>

namespace solsim::io {

SvgWriter::SvgWriter(double width, double height, const std::string& background)
    : m_w(width), m_h(height), m_background(background) {}

std::string SvgWriter::escape(const std::string& s) {
    std::string out;
    for (char c : s) {
        switch (c) {
            case '&': out += "&amp;"; break;
            case '<': out += "&lt;"; break;
            case '>': out += "&gt;"; break;
            case '"': out += "&quot;"; break;
            default: out.push_back(c);
        }
    }
    return out;
}

void SvgWriter::comment(const std::string& text) { m_body << "<!-- " << solsim::str::replaceAll(text, "--", "- -") << " -->\n"; }

void SvgWriter::rect(double x, double y, double w, double h, const std::string& fill, const std::string& stroke,
                     double strokeWidth) {
    m_body << solsim::str::format("<rect x=\"%.2f\" y=\"%.2f\" width=\"%.2f\" height=\"%.2f\" fill=\"%s\" stroke=\"%s\" "
                          "stroke-width=\"%.2f\"/>\n",
                          x, y, w, h, fill.c_str(), stroke.c_str(), strokeWidth);
}

void SvgWriter::circle(double cx, double cy, double r, const std::string& fill, const std::string& stroke,
                       double strokeWidth, double opacity) {
    m_body << solsim::str::format("<circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" fill=\"%s\" stroke=\"%s\" stroke-width=\"%.2f\" "
                          "opacity=\"%.3f\"/>\n",
                          cx, cy, r, fill.c_str(), stroke.c_str(), strokeWidth, opacity);
}

void SvgWriter::line(double x1, double y1, double x2, double y2, const std::string& stroke, double width,
                     double opacity, const std::string& dash) {
    m_body << solsim::str::format("<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"%s\" stroke-width=\"%.2f\" "
                          "opacity=\"%.3f\"",
                          x1, y1, x2, y2, stroke.c_str(), width, opacity);
    if (!dash.empty()) m_body << " stroke-dasharray=\"" << dash << "\"";
    m_body << "/>\n";
}

void SvgWriter::polyline(const std::vector<std::pair<double, double>>& pts, const std::string& stroke, double width,
                         double opacity, bool closed, const std::string& dash) {
    if (pts.size() < 2) return;
    m_body << (closed ? "<polygon" : "<polyline") << " fill=\"none\" stroke=\"" << stroke << "\" stroke-width=\""
           << solsim::str::format("%.2f", width) << "\" opacity=\"" << solsim::str::format("%.3f", opacity)
           << "\" stroke-linejoin=\"round\"";
    if (!dash.empty()) m_body << " stroke-dasharray=\"" << dash << "\"";
    m_body << " points=\"";
    for (std::size_t i = 0; i < pts.size(); ++i) {
        if (i) m_body << ' ';
        m_body << solsim::str::format("%.2f,%.2f", pts[i].first, pts[i].second);
    }
    m_body << "\"/>\n";
}

void SvgWriter::text(double x, double y, const std::string& content, const std::string& fill, double size,
                     const std::string& anchor, const std::string& weight) {
    m_body << solsim::str::format("<text x=\"%.2f\" y=\"%.2f\" fill=\"%s\" font-size=\"%.1f\" text-anchor=\"%s\" "
                          "font-weight=\"%s\">",
                          x, y, fill.c_str(), size, anchor.c_str(), weight.c_str())
           << escape(content) << "</text>\n";
}

void SvgWriter::beginGroup(const std::string& attributes) { m_body << "<g " << attributes << ">\n"; }
void SvgWriter::endGroup() { m_body << "</g>\n"; }
void SvgWriter::raw(const std::string& markup) { m_body << markup << "\n"; }

std::string SvgWriter::str() const {
    std::ostringstream os;
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
       << solsim::str::format("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%.0f\" height=\"%.0f\" viewBox=\"0 0 %.0f %.0f\" "
                      "font-family=\"Helvetica, Arial, sans-serif\">\n",
                      m_w, m_h, m_w, m_h)
       << solsim::str::format("<rect width=\"100%%\" height=\"100%%\" fill=\"%s\"/>\n", m_background.c_str()) << m_body.str()
       << "</svg>\n";
    return os.str();
}

bool SvgWriter::save(const std::string& path, std::string& error) const {
    std::ofstream f(path, std::ios::binary);
    if (!f) {
        error = "cannot open '" + path + "' for writing";
        return false;
    }
    f << str();
    if (!f) {
        error = "write error on '" + path + "'";
        return false;
    }
    return true;
}

} // namespace solsim::io
