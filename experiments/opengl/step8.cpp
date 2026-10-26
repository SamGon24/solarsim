// =============================================================================
//  Step 8: HDR rendering, bloom and tone mapping.
//  Builds on step 7. The scene is drawn into an off-screen floating-point image
//  (multisampled) and resolved into a texture. Then:
//    bright pass -> separable Gaussian blur (ping-pong, half resolution)
//    -> composite scene + bloom -> exposure -> tone map -> gamma -> screen.
//  Controls:
//    left mouse drag  orbit        mouse wheel  zoom
//    Tab / 0-9        focus body   Z            reset camera
//    Space            pause        [ ]          slower / faster time
//    R                reverse time F            display flux <-> physical flux
//    M                radial scale (sqrt / log / linear)
//    T                trails on/off              S   stars on/off
//    H                show/hide the on-screen legend
//    , .              exposure down / up         Y   tone map (clamp/reinhard/aces)
//    B                bloom on/off               - =  bloom weaker / stronger
//    Esc or Q         quit
// =============================================================================
#define GL_SILENCE_DEPRECATION
#include <SDL.h>
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#else
#define GL_GLEXT_PROTOTYPES
#include <SDL_opengl.h>
#endif

#include "astro/JulianDate.hpp"
#include "physics/Simulation.hpp"
#include "physics/SystemBuilder.hpp"
#include "render/Camera.hpp"
#include "render/Color.hpp"
#include "render/DisplayMapping.hpp"
#include "render/Font.hpp"
#include "render/StarField.hpp"
#include "render/Trails.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <utility>
#include <vector>

using solsim::Vec3;
using solsim::render::Color;

// -----------------------------------------------------------------------------
//  Shaders: bodies (lit spheres), lines (trails) and stars (points)
// -----------------------------------------------------------------------------
static const char* kBodyVert = R"(#version 330 core
layout(location = 0) in vec3 aPos;
uniform mat4 uMVP;
uniform mat4 uModel;
out vec3 vN;
out vec3 vWorld;
void main() {
    vN = mat3(uModel) * aPos;
    vWorld = (uModel * vec4(aPos, 1.0)).xyz;
    gl_Position = uMVP * vec4(aPos, 1.0);
}
)";

static const char* kBodyFrag = R"(#version 330 core
in vec3 vN;
in vec3 vWorld;
uniform vec3 uColor;
uniform vec3 uLightPos;
uniform float uFlux;
uniform float uEmissive;
out vec4 FragColor;
void main() {
    vec3 lin;
    if (uEmissive > 0.5) {
        lin = uColor * 12.0;
    } else {
        vec3 N = normalize(vN);
        vec3 L = normalize(uLightPos - vWorld);
        float diffuse = max(dot(N, L), 0.0);
        lin = uColor * (0.02 + diffuse * uFlux * 2.0);
    }
    FragColor = vec4(lin, 1.0);
}
)";

static const char* kLineVert = R"(#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in float aAlpha;
uniform mat4 uViewProj;
out float vAlpha;
void main() {
    vAlpha = aAlpha;
    gl_Position = uViewProj * vec4(aPos, 1.0);
}
)";

static const char* kLineFrag = R"(#version 330 core
in float vAlpha;
uniform vec3 uColor;
uniform float uStrength;
out vec4 FragColor;
void main() {
    FragColor = vec4(uColor * 2.5, vAlpha * uStrength);
}
)";

// A star is a direction, not a position: w = 0 makes the translation part of
// any matrix irrelevant, so stars stay at infinity however far the camera moves.
// Such a point lands just beyond the far plane and would be clipped, so we pin
// its depth to the far plane with .xyww (z = w).
static const char* kStarVert = R"(#version 330 core
layout(location = 0) in vec3 aDir;
layout(location = 1) in vec3 aColor;
layout(location = 2) in float aSize;
uniform mat4 uViewProj;
uniform float uPointScale;
out vec3 vColor;
void main() {
    vColor = aColor;
    vec4 p = uViewProj * vec4(aDir, 0.0);
    gl_Position = p.xyww;
    gl_PointSize = max(1.0, aSize * uPointScale);
}
)";

static const char* kStarFrag = R"(#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
    float d = length(gl_PointCoord * 2.0 - 1.0);
    float falloff = smoothstep(1.0, 0.2, d);
    FragColor = vec4(vColor * falloff, 1.0);
}
)";

// Fullscreen pass: one oversized triangle generated from gl_VertexID covers the
// whole screen (vertices at (0,0), (2,0), (0,2) in UV space).
static const char* kPostVert = R"(#version 330 core
out vec2 vUV;
void main() {
    vec2 p = vec2(float((gl_VertexID << 1) & 2), float(gl_VertexID & 2));
    vUV = p;
    gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0);
}
)";

// Bright pass: keep only what is brighter than the threshold (with a soft knee so
// there is no hard edge). Reading a full-size texture while writing a half-size
// target also downsamples it for free through bilinear filtering.
static const char* kBrightFrag = R"(#version 330 core
in vec2 vUV;
uniform sampler2D uScene;
uniform float uExposure;
uniform float uThreshold;
uniform float uKnee;
out vec4 FragColor;
void main() {
    vec3 c = texture(uScene, vUV).rgb * uExposure;
    float br = max(c.r, max(c.g, c.b));
    float soft = clamp(br - uThreshold + uKnee, 0.0, 2.0 * uKnee);
    soft = soft * soft / (4.0 * uKnee + 1e-5);
    float contrib = max(soft, br - uThreshold) / max(br, 1e-5);
    FragColor = vec4(c * contrib, 1.0);
}
)";

// One direction of a separable Gaussian blur (9 taps). uDir is the step in UV
// units: (1/width, 0) for horizontal, (0, 1/height) for vertical.
static const char* kBlurFrag = R"(#version 330 core
in vec2 vUV;
uniform sampler2D uImage;
uniform vec2 uDir;
out vec4 FragColor;
const float kWeight[5] = float[](0.227027, 0.1945946, 0.1216216, 0.054054, 0.016216);
void main() {
    vec3 r = texture(uImage, vUV).rgb * kWeight[0];
    for (int i = 1; i < 5; ++i) {
        r += texture(uImage, vUV + uDir * float(i)).rgb * kWeight[i];
        r += texture(uImage, vUV - uDir * float(i)).rgb * kWeight[i];
    }
    FragColor = vec4(r, 1.0);
}
)";

// The only place where bloom is added and exposure, tone mapping and gamma happen.
static const char* kPostFrag = R"(#version 330 core
in vec2 vUV;
uniform sampler2D uScene;
uniform sampler2D uBloom;
uniform float uExposure;
uniform float uBloomStrength;
uniform int uTonemap;
out vec4 FragColor;

vec3 aces(vec3 x) {
    return clamp((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14), 0.0, 1.0);
}

void main() {
    vec3 c = texture(uScene, vUV).rgb * uExposure + texture(uBloom, vUV).rgb * uBloomStrength;
    if (uTonemap == 1) c = c / (1.0 + c);
    else if (uTonemap == 2) c = aces(c);
    else c = clamp(c, 0.0, 1.0);
    FragColor = vec4(pow(c, vec3(1.0 / 2.2)), 1.0);
}
)";

static const char* kTonemapNames[] = {"clamp", "reinhard", "aces"};

// Screen-space textured rectangle for the legend. There is no vertex buffer at
// all: the four corners are generated from gl_VertexID (0..3) and drawn as a
// triangle strip. uRect is x, y, width, height in pixels from the top-left.
static const char* kUiVert = R"(#version 330 core
uniform vec4 uRect;
uniform vec2 uViewport;
out vec2 vUV;
void main() {
    vec2 c = vec2(float(gl_VertexID & 1), float(gl_VertexID >> 1));
    vUV = c;
    vec2 px = uRect.xy + c * uRect.zw;
    gl_Position = vec4(px.x / uViewport.x * 2.0 - 1.0, 1.0 - px.y / uViewport.y * 2.0, 0.0, 1.0);
}
)";

static const char* kUiFrag = R"(#version 330 core
in vec2 vUV;
uniform sampler2D uTex;
out vec4 FragColor;
void main() {
    FragColor = texture(uTex, vUV);
}
)";

static GLuint compileShader(GLenum type, const char* src) {
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok = 0;
    glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char log[1024];
        glGetShaderInfoLog(s, sizeof log, nullptr, log);
        std::printf("shader compile error:\n%s\n", log);
    }
    return s;
}

static GLuint makeProgram(const char* vsSrc, const char* fsSrc) {
    GLuint vs = compileShader(GL_VERTEX_SHADER, vsSrc);
    GLuint fs = compileShader(GL_FRAGMENT_SHADER, fsSrc);
    GLuint p = glCreateProgram();
    glAttachShader(p, vs);
    glAttachShader(p, fs);
    glLinkProgram(p);
    GLint ok = 0;
    glGetProgramiv(p, GL_LINK_STATUS, &ok);
    if (!ok) {
        char log[1024];
        glGetProgramInfoLog(p, sizeof log, nullptr, log);
        std::printf("link error:\n%s\n", log);
    }
    glDeleteShader(vs);
    glDeleteShader(fs);
    return p;
}

// -----------------------------------------------------------------------------
//  Matrices (column-major: element (row r, column c) is m[c * 4 + r])
// -----------------------------------------------------------------------------
struct Mat4 { float m[16]; };

static Mat4 identity() {
    Mat4 r = {};
    r.m[0] = r.m[5] = r.m[10] = r.m[15] = 1.0f;
    return r;
}

static Mat4 mul(const Mat4& a, const Mat4& b) {
    Mat4 r = {};
    for (int c = 0; c < 4; ++c)
        for (int row = 0; row < 4; ++row) {
            float s = 0;
            for (int k = 0; k < 4; ++k) s += a.m[k * 4 + row] * b.m[c * 4 + k];
            r.m[c * 4 + row] = s;
        }
    return r;
}

static Mat4 translate(float x, float y, float z) {
    Mat4 r = identity();
    r.m[12] = x;
    r.m[13] = y;
    r.m[14] = z;
    return r;
}

static Mat4 scaleUniform(float s) {
    Mat4 r = identity();
    r.m[0] = r.m[5] = r.m[10] = s;
    return r;
}

static Mat4 viewRotation(const Vec3& right, const Vec3& up, const Vec3& forward) {
    Mat4 r = identity();
    r.m[0] = float(right.x);     r.m[4] = float(right.y);     r.m[8] = float(right.z);
    r.m[1] = float(up.x);        r.m[5] = float(up.y);        r.m[9] = float(up.z);
    r.m[2] = float(-forward.x);  r.m[6] = float(-forward.y);  r.m[10] = float(-forward.z);
    return r;
}

static Mat4 perspective(float fovY, float aspect, float nearZ, float farZ) {
    float f = 1.0f / std::tan(fovY * 0.5f);
    Mat4 r = {};
    r.m[0] = f / aspect;
    r.m[5] = f;
    r.m[10] = (farZ + nearZ) / (nearZ - farZ);
    r.m[11] = -1.0f;
    r.m[14] = 2.0f * farZ * nearZ / (nearZ - farZ);
    return r;
}

static void buildSphere(int stacks, int slices, std::vector<float>& verts, std::vector<unsigned int>& idx) {
    const float PI = 3.14159265358979f;
    for (int i = 0; i <= stacks; ++i) {
        float phi = PI * float(i) / float(stacks);
        float y = std::cos(phi);
        float r = std::sin(phi);
        for (int j = 0; j <= slices; ++j) {
            float theta = 2.0f * PI * float(j) / float(slices);
            verts.push_back(r * std::cos(theta));
            verts.push_back(y);
            verts.push_back(r * std::sin(theta));
        }
    }
    for (int i = 0; i < stacks; ++i)
        for (int j = 0; j < slices; ++j) {
            unsigned int a = i * (slices + 1) + j;
            unsigned int b = a + slices + 1;
            idx.push_back(a);
            idx.push_back(b);
            idx.push_back(a + 1);
            idx.push_back(a + 1);
            idx.push_back(b);
            idx.push_back(b + 1);
        }
}

// -----------------------------------------------------------------------------
//  HDR render targets
//    scene pass   -> msaaFbo    (4x multisampled RGBA16F + depth renderbuffers)
//    resolve blit -> resolveFbo (a normal RGBA16F texture we can sample)
//    bright pass  -> bloomFbo[0]; blur ping-pongs between bloomFbo[0] and [1]
//    post pass    -> screen
// -----------------------------------------------------------------------------
struct Targets {
    GLuint msaaFbo = 0, msaaColor = 0, msaaDepth = 0;
    GLuint resolveFbo = 0, resolveTex = 0;
    GLuint bloomFbo[2] = {0, 0}, bloomTex[2] = {0, 0};  // half resolution, ping-pong pair
    int w = 0, h = 0, bw = 0, bh = 0;
};

static void destroyTargets(Targets& t) {
    if (t.msaaFbo) glDeleteFramebuffers(1, &t.msaaFbo);
    if (t.msaaColor) glDeleteRenderbuffers(1, &t.msaaColor);
    if (t.msaaDepth) glDeleteRenderbuffers(1, &t.msaaDepth);
    if (t.resolveFbo) glDeleteFramebuffers(1, &t.resolveFbo);
    if (t.resolveTex) glDeleteTextures(1, &t.resolveTex);
    for (int i = 0; i < 2; ++i) {
        if (t.bloomFbo[i]) glDeleteFramebuffers(1, &t.bloomFbo[i]);
        if (t.bloomTex[i]) glDeleteTextures(1, &t.bloomTex[i]);
    }
    t = Targets{};
}

static bool createTargets(Targets& t, int w, int h, int samples) {
    destroyTargets(t);
    t.w = w;
    t.h = h;

    glGenFramebuffers(1, &t.msaaFbo);
    glBindFramebuffer(GL_FRAMEBUFFER, t.msaaFbo);
    glGenRenderbuffers(1, &t.msaaColor);
    glBindRenderbuffer(GL_RENDERBUFFER, t.msaaColor);
    glRenderbufferStorageMultisample(GL_RENDERBUFFER, samples, GL_RGBA16F, w, h);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, t.msaaColor);
    glGenRenderbuffers(1, &t.msaaDepth);
    glBindRenderbuffer(GL_RENDERBUFFER, t.msaaDepth);
    glRenderbufferStorageMultisample(GL_RENDERBUFFER, samples, GL_DEPTH_COMPONENT24, w, h);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, t.msaaDepth);
    const bool okMsaa = glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;

    glGenTextures(1, &t.resolveTex);
    glBindTexture(GL_TEXTURE_2D, t.resolveTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glGenFramebuffers(1, &t.resolveFbo);
    glBindFramebuffer(GL_FRAMEBUFFER, t.resolveFbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, t.resolveTex, 0);
    const bool okResolve = glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;

    t.bw = std::max(1, w / 2);
    t.bh = std::max(1, h / 2);
    bool okBloom = true;
    for (int i = 0; i < 2; ++i) {
        glGenTextures(1, &t.bloomTex[i]);
        glBindTexture(GL_TEXTURE_2D, t.bloomTex[i]);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, t.bw, t.bh, 0, GL_RGBA, GL_FLOAT, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glGenFramebuffers(1, &t.bloomFbo[i]);
        glBindFramebuffer(GL_FRAMEBUFFER, t.bloomFbo[i]);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, t.bloomTex[i], 0);
        okBloom = okBloom && glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE;
    }

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    return okMsaa && okResolve && okBloom;
}

// -----------------------------------------------------------------------------
//  Legend texture: text from SolarSim's 5x7 bitmap font rasterised on the CPU
//  into an RGBA image (semi-transparent panel + coloured glyphs), then uploaded
//  once as an OpenGL texture.
// -----------------------------------------------------------------------------
struct Texture {
    GLuint id = 0;
    int w = 0, h = 0;
};

static Texture makeLegendTexture(const std::vector<std::pair<std::string, std::string>>& rows) {
    using solsim::render::BitmapFont;
    const int pad = 7, gap = 14, rowH = BitmapFont::kGlyphH + 5;
    int keyW = 0, descW = 0, headW = 0;
    for (const auto& r : rows) {
        if (r.second.empty()) headW = std::max(headW, BitmapFont::textWidth(r.first));
        else {
            keyW = std::max(keyW, BitmapFont::textWidth(r.first));
            descW = std::max(descW, BitmapFont::textWidth(r.second));
        }
    }
    Texture t;
    t.w = std::max(headW, keyW + gap + descW) + 2 * pad;
    t.h = int(rows.size()) * rowH + 2 * pad;

    std::vector<unsigned char> px(size_t(t.w) * size_t(t.h) * 4);
    auto put = [&](int x, int y, unsigned char r, unsigned char g, unsigned char b, unsigned char a) {
        if (x < 0 || y < 0 || x >= t.w || y >= t.h) return;
        unsigned char* p = &px[(size_t(y) * size_t(t.w) + size_t(x)) * 4];
        p[0] = r; p[1] = g; p[2] = b; p[3] = a;
    };
    for (int y = 0; y < t.h; ++y)
        for (int x = 0; x < t.w; ++x) {
            const bool edge = x == 0 || y == 0 || x == t.w - 1 || y == t.h - 1;
            if (edge) put(x, y, 70, 90, 130, 230);
            else put(x, y, 8, 12, 24, 170);
        }

    auto text = [&](int x, int y, const std::string& s, unsigned char r, unsigned char g, unsigned char b) {
        for (char c : s) {
            for (int row = 0; row < BitmapFont::kGlyphH; ++row)
                for (int col = 0; col < BitmapFont::kGlyphW; ++col)
                    if (BitmapFont::pixel(c, col, row)) put(x + col, y + row, r, g, b, 255);
            x += BitmapFont::kAdvance;
        }
    };
    for (size_t i = 0; i < rows.size(); ++i) {
        const int y = pad + int(i) * rowH + 2;
        if (rows[i].second.empty()) text(pad, y, rows[i].first, 140, 200, 255);
        else {
            text(pad, y, rows[i].first, 255, 214, 102);
            text(pad + keyW + gap, y, rows[i].second, 222, 230, 242);
        }
    }

    glGenTextures(1, &t.id);
    glBindTexture(GL_TEXTURE_2D, t.id);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, t.w, t.h, 0, GL_RGBA, GL_UNSIGNED_BYTE, px.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    return t;
}

// -----------------------------------------------------------------------------
//  Helpers around the SolarSim classes
// -----------------------------------------------------------------------------
static double systemExtent(const solsim::render::DisplayMapping& map, const solsim::physics::NBodySystem& sys) {
    double outer = 1.0;
    for (std::size_t i = 1; i < sys.namedCount(); ++i) {
        if (sys.info(i).parent != 0) continue;
        outer = std::max(outer, map.position(sys, i).norm());
    }
    return outer;
}

static double frameDistance(const solsim::render::DisplayMapping& map, const solsim::physics::NBodySystem& sys,
                            int focus, double homeDistance) {
    if (focus <= 0) return homeDistance;
    const Vec3 c = map.position(sys, std::size_t(focus));
    double span = map.displayRadius(sys, std::size_t(focus)) * 14.0;
    for (std::size_t i = 0; i < sys.namedCount(); ++i)
        if (sys.info(i).parent == focus) span = std::max(span, 2.2 * (map.position(sys, i) - c).norm());
    return span;
}

int main(int argc, char** argv) {
    long maxFrames = -1;
    for (int i = 1; i + 1 < argc; ++i)
        if (std::strcmp(argv[i], "--frames") == 0) maxFrames = std::atol(argv[i + 1]);

    // ---- simulation -------------------------------------------------------------------
    const double epoch = solsim::astro::utcToTT(solsim::astro::nowJDUtc());
    solsim::physics::BuildOptions bo;
    bo.epochJD = epoch;
    solsim::physics::Simulation sim(solsim::physics::buildSystem(bo), epoch, solsim::physics::SimulationConfig{});
    const auto& sys = sim.system();

    solsim::render::DisplayMapping map;
    map.prepare(sys);
    solsim::render::TrailRecorder trails;
    trails.record(sys, sim.time());
    const solsim::render::StarField starField;

    solsim::render::Camera cam;
    cam.setFieldOfView(45.0 * 3.14159265358979 / 180.0);
    double homeDistance = systemExtent(map, sys) * 2.4;
    const double homeYaw = -100.0 * 3.14159265358979 / 180.0;
    const double homePitch = 28.0 * 3.14159265358979 / 180.0;
    auto goHome = [&]() {
        homeDistance = systemExtent(map, sys) * 2.4;
        cam.setTarget(Vec3{});
        cam.setAngles(homeYaw, homePitch);
        cam.setDistance(homeDistance);
    };
    goHome();

    // Advance the simulation in small chunks and record trail samples between them.
    auto advanceWithTrails = [&](double days) {
        const double chunk = std::max(0.25, 4.0 * sim.timeStep());
        const int sub = std::clamp(int(std::ceil(std::fabs(days) / chunk)), 1, 64);
        for (int k = 0; k < sub; ++k) {
            sim.advance(days / sub);
            trails.record(sys, sim.time());
        }
    };

    // ---- window and OpenGL context --------------------------------------------------
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        std::printf("SDL: %s\n", SDL_GetError());
        return 1;
    }
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_Window* win = SDL_CreateWindow("SolarSim GL", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1280, 720,
                                       SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_GLContext ctx = win ? SDL_GL_CreateContext(win) : nullptr;
    if (!ctx) {
        std::printf("GL: %s\n", SDL_GetError());
        return 1;
    }
    SDL_GL_SetSwapInterval(1);
    std::printf("GL_VERSION:  %s\nGL_RENDERER: %s\n", glGetString(GL_VERSION), glGetString(GL_RENDERER));

    // ---- programs and uniform locations ------------------------------------------------
    GLuint bodyProg = makeProgram(kBodyVert, kBodyFrag);
    GLint locMVP = glGetUniformLocation(bodyProg, "uMVP");
    GLint locModel = glGetUniformLocation(bodyProg, "uModel");
    GLint locColor = glGetUniformLocation(bodyProg, "uColor");
    GLint locLight = glGetUniformLocation(bodyProg, "uLightPos");
    GLint locFlux = glGetUniformLocation(bodyProg, "uFlux");
    GLint locEmissive = glGetUniformLocation(bodyProg, "uEmissive");

    GLuint lineProg = makeProgram(kLineVert, kLineFrag);
    GLint locLineVP = glGetUniformLocation(lineProg, "uViewProj");
    GLint locLineColor = glGetUniformLocation(lineProg, "uColor");
    GLint locLineStrength = glGetUniformLocation(lineProg, "uStrength");

    GLuint starProg = makeProgram(kStarVert, kStarFrag);
    GLint locStarVP = glGetUniformLocation(starProg, "uViewProj");
    GLint locStarScale = glGetUniformLocation(starProg, "uPointScale");

    GLuint uiProg = makeProgram(kUiVert, kUiFrag);
    GLint locUiRect = glGetUniformLocation(uiProg, "uRect");
    GLint locUiViewport = glGetUniformLocation(uiProg, "uViewport");

    GLuint postProg = makeProgram(kPostVert, kPostFrag);
    GLint locPostExposure = glGetUniformLocation(postProg, "uExposure");
    GLint locPostTonemap = glGetUniformLocation(postProg, "uTonemap");
    GLint locPostBloom = glGetUniformLocation(postProg, "uBloomStrength");
    glUseProgram(postProg);
    glUniform1i(glGetUniformLocation(postProg, "uScene"), 0);   // texture unit 0
    glUniform1i(glGetUniformLocation(postProg, "uBloom"), 1);   // texture unit 1

    GLuint brightProg = makeProgram(kPostVert, kBrightFrag);
    GLint locBrightExposure = glGetUniformLocation(brightProg, "uExposure");
    GLint locBrightThreshold = glGetUniformLocation(brightProg, "uThreshold");
    GLint locBrightKnee = glGetUniformLocation(brightProg, "uKnee");

    GLuint blurProg = makeProgram(kPostVert, kBlurFrag);
    GLint locBlurDir = glGetUniformLocation(blurProg, "uDir");

    // ---- sphere mesh --------------------------------------------------------------------
    std::vector<float> verts;
    std::vector<unsigned int> idx;
    buildSphere(48, 96, verts, idx);
    GLuint sphereVao = 0, sphereVbo = 0, sphereEbo = 0;
    glGenVertexArrays(1, &sphereVao);
    glGenBuffers(1, &sphereVbo);
    glGenBuffers(1, &sphereEbo);
    glBindVertexArray(sphereVao);
    glBindBuffer(GL_ARRAY_BUFFER, sphereVbo);
    glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(verts.size() * sizeof(float)), verts.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphereEbo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, GLsizeiptr(idx.size() * sizeof(unsigned int)), idx.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    // ---- star buffer: built once, never changes (static) ---------------------------------
    // Per star: direction (3), colour * brightness (3), point size in pixels (1).
    std::vector<float> starData;
    for (const auto& s : starField.stars()) {
        const float flux = 0.55f * std::pow(10.0f, -0.24f * (s.magnitude - 1.0f));
        const float size = 2.0f * (0.55f + 0.22f * std::max(0.0f, 2.5f - s.magnitude)) + 0.5f;
        const Color c = s.color * flux;
        starData.insert(starData.end(), {float(s.direction.x), float(s.direction.y), float(s.direction.z), c.r, c.g, c.b, size});
    }
    const GLsizei starCount = GLsizei(starData.size() / 7);
    GLuint starVao = 0, starVbo = 0;
    glGenVertexArrays(1, &starVao);
    glGenBuffers(1, &starVbo);
    glBindVertexArray(starVao);
    glBindBuffer(GL_ARRAY_BUFFER, starVbo);
    glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(starData.size() * sizeof(float)), starData.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 7 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(2, 1, GL_FLOAT, GL_FALSE, 7 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);
    glBindVertexArray(0);

    // ---- trail buffer: rewritten every frame (dynamic) ----------------------------------
    // Per vertex: position (3, camera-relative), alpha (1).
    GLuint trailVao = 0, trailVbo = 0;
    glGenVertexArrays(1, &trailVao);
    glGenBuffers(1, &trailVbo);
    glBindVertexArray(trailVao);
    glBindBuffer(GL_ARRAY_BUFFER, trailVbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
    std::vector<float> trailScratch;

    // ---- legend: one texture, drawn with an attribute-less quad ---------------------------
    const Texture legend = makeLegendTexture({
        {"SolarSim GL controls", ""},
        {"", ""},
        {"Mouse drag", "orbit"},
        {"Wheel", "zoom"},
        {"Tab", "next body"},
        {"0-9", "Sun .. Pluto"},
        {"Z", "reset camera"},
        {"", ""},
        {"Space", "pause"},
        {"[  ]", "time slower / faster"},
        {"R", "reverse time"},
        {"", ""},
        {"T", "trails"},
        {"S", "stars"},
        {"F", "flux: display / physical"},
        {"M", "radial scale"},
        {",  .", "exposure down / up"},
        {"Y", "tone map: clamp / reinhard / aces"},
        {"B", "bloom on / off"},
        {"-  =", "bloom weaker / stronger"},
        {"", ""},
        {"H", "show / hide legend"},
        {"Esc  Q", "quit"},
    });
    GLuint emptyVao = 0;
    glGenVertexArrays(1, &emptyVao);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_PROGRAM_POINT_SIZE);

    // ---- state ------------------------------------------------------------------------
    int focus = 0;
    bool paused = false;
    bool physicalFlux = false;
    bool showTrails = true;
    bool showStars = true;
    bool showLegend = true;
    float exposure = 1.0f;
    int tonemap = 1;
    bool bloomOn = true;
    float bloomStrength = 0.6f;
    const int blurIterations = 4;        // each iteration = one horizontal + one vertical pass
    const float bloomThreshold = 1.2f;   // in exposed HDR units
    const float bloomKnee = 0.5f;
    GLint maxSamples = 4;
    glGetIntegerv(GL_MAX_SAMPLES, &maxSamples);
    const int samples = std::max(0, std::min(4, int(maxSamples)));
    Targets targets;
    double warp = 10.0;
    double direction = 1.0;
    double fps = 60.0;
    double titleTimer = 0.0;
    long frame = 0;
    Uint64 last = SDL_GetPerformanceCounter();
    const double perfFreq = double(SDL_GetPerformanceFrequency());

    bool running = true;
    while (running) {
        // ---- events -------------------------------------------------------------------
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = false;
            else if (e.type == SDL_MOUSEMOTION && (e.motion.state & SDL_BUTTON_LMASK)) {
                cam.orbit(-e.motion.xrel * 0.005, e.motion.yrel * 0.005);
            } else if (e.type == SDL_MOUSEWHEEL) {
                cam.zoom(std::pow(0.9, double(e.wheel.y)));
            } else if (e.type == SDL_KEYDOWN && !e.key.repeat) {
                const SDL_Keycode k = e.key.keysym.sym;
                if (k == SDLK_ESCAPE || k == SDLK_q) running = false;
                else if (k == SDLK_SPACE) paused = !paused;
                else if (k == SDLK_RIGHTBRACKET) warp = std::min(warp * 2.0, 1e6);
                else if (k == SDLK_LEFTBRACKET) warp = std::max(warp / 2.0, 1e-4);
                else if (k == SDLK_r) direction = -direction;
                else if (k == SDLK_f) physicalFlux = !physicalFlux;
                else if (k == SDLK_t) showTrails = !showTrails;
                else if (k == SDLK_s) showStars = !showStars;
                else if (k == SDLK_h) showLegend = !showLegend;
                else if (k == SDLK_COMMA) exposure = std::max(exposure / 1.25f, 0.05f);
                else if (k == SDLK_PERIOD) exposure = std::min(exposure * 1.25f, 50.0f);
                else if (k == SDLK_y) tonemap = (tonemap + 1) % 3;
                else if (k == SDLK_b) bloomOn = !bloomOn;
                else if (k == SDLK_MINUS) bloomStrength = std::max(bloomStrength / 1.25f, 0.02f);
                else if (k == SDLK_EQUALS) bloomStrength = std::min(bloomStrength * 1.25f, 5.0f);
                else if (k == SDLK_z) { focus = 0; goHome(); }
                else if (k == SDLK_m) {
                    using solsim::render::ScaleMode;
                    map.setMode(map.mode() == ScaleMode::Sqrt ? ScaleMode::Log
                                : map.mode() == ScaleMode::Log ? ScaleMode::Linear : ScaleMode::Sqrt);
                    focus = 0;
                    goHome();
                } else if (k == SDLK_TAB || (k >= SDLK_0 && k <= SDLK_9)) {
                    if (k == SDLK_TAB) {
                        focus = (focus + 1) % int(sys.namedCount());
                    } else {
                        static const char* byDigit[] = {"Sun", "Mercury", "Venus", "Earth", "Mars",
                                                        "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"};
                        const int found = sys.indexOf(byDigit[k - SDLK_0]);
                        if (found >= 0) focus = found;
                    }
                    cam.setDistance(frameDistance(map, sys, focus, homeDistance));
                }
            }
        }

        // ---- time ---------------------------------------------------------------------
        const Uint64 now = SDL_GetPerformanceCounter();
        const double dt = std::min(double(now - last) / perfFreq, 0.1);
        last = now;
        fps += (1.0 / std::max(dt, 1e-4) - fps) * 0.05;
        if (!paused) advanceWithTrails(direction * warp * dt);

        // ---- camera --------------------------------------------------------------------
        map.prepare(sys);
        cam.setTarget(focus > 0 ? map.position(sys, std::size_t(focus)) : Vec3{});

        int w = 1, h = 1, winW = 1, winH = 1;
        SDL_GL_GetDrawableSize(win, &w, &h);
        SDL_GetWindowSize(win, &winW, &winH);
        cam.setViewport(w, h);

        // (Re)create the off-screen targets whenever the drawable size changes.
        if (targets.w != w || targets.h != h) {
            if (!createTargets(targets, w, h, samples)) {
                std::printf("framebuffer incomplete\n");
                return 1;
            }
        }

        // ---- scene pass: everything below draws into the HDR target -------------------------
        glBindFramebuffer(GL_FRAMEBUFFER, targets.msaaFbo);
        glViewport(0, 0, w, h);
        glClearColor(0.0004f, 0.0004f, 0.0015f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const Vec3 eye = cam.eye();
        const double dist = cam.distance();
        const float nearZ = float(std::max(1e-6, dist * 0.001));
        const float farZ = float(dist * 20.0 + 50.0);
        const Mat4 proj = perspective(float(cam.fieldOfView()), float(w) / float(h), nearZ, farZ);
        const Mat4 viewProj = mul(proj, viewRotation(cam.right(), cam.up(), cam.forward()));
        const double focalPx = 0.5 * h / std::tan(0.5 * cam.fieldOfView());

        // ---- 1. stars: first, additive, no depth ---------------------------------------
        if (showStars) {
            glDisable(GL_DEPTH_TEST);
            glEnable(GL_BLEND);
            glBlendFunc(GL_ONE, GL_ONE);
            glUseProgram(starProg);
            glUniformMatrix4fv(locStarVP, 1, GL_FALSE, viewProj.m);
            glUniform1f(locStarScale, float(w) / float(std::max(1, winW)));
            glBindVertexArray(starVao);
            glDrawArrays(GL_POINTS, 0, starCount);
            glDisable(GL_BLEND);
            glEnable(GL_DEPTH_TEST);
        }

        // ---- 2. bodies: opaque, depth tested ---------------------------------------------
        glUseProgram(bodyProg);
        glBindVertexArray(sphereVao);
        const Vec3 sunRel = map.position(sys, 0) - eye;
        glUniform3f(locLight, float(sunRel.x), float(sunRel.y), float(sunRel.z));

        for (std::size_t i = 0; i < sys.namedCount(); ++i) {
            const Vec3 p = map.position(sys, i);
            const Vec3 rel = p - eye;
            const double depth = rel.dot(cam.forward());
            if (depth <= nearZ) continue;
            const double radius = std::max(map.displayRadius(sys, i), 2.5 * depth / focalPx);

            double flux = 1.0;
            if (i > 0) {
                const double d = physicalFlux ? (sys.pos(i) - sys.pos(0)).norm() : p.norm();
                flux = 1.0 / std::max(d * d, 1e-6);
            }
            const Color col = Color::fromHex(sys.info(i).color);

            const Mat4 model = mul(translate(float(rel.x), float(rel.y), float(rel.z)), scaleUniform(float(radius)));
            const Mat4 mvp = mul(viewProj, model);
            glUniformMatrix4fv(locMVP, 1, GL_FALSE, mvp.m);
            glUniformMatrix4fv(locModel, 1, GL_FALSE, model.m);
            glUniform3f(locColor, col.r, col.g, col.b);
            glUniform1f(locFlux, float(flux));
            glUniform1f(locEmissive, i == 0 ? 1.0f : 0.0f);
            glDrawElements(GL_TRIANGLES, GLsizei(idx.size()), GL_UNSIGNED_INT, (void*)0);
        }

        // ---- 3. trails: last, blended, depth tested but not depth written -----------------
        if (showTrails) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glDepthMask(GL_FALSE);
            glUseProgram(lineProg);
            glUniformMatrix4fv(locLineVP, 1, GL_FALSE, viewProj.m);
            glUniform1f(locLineStrength, 0.85f);
            glBindVertexArray(trailVao);
            glBindBuffer(GL_ARRAY_BUFFER, trailVbo);

            for (std::size_t i = 1; i < sys.namedCount(); ++i) {
                const std::size_t n = trails.count(i);
                if (n < 2) continue;
                const int parent = trails.parent(i);
                trailScratch.clear();
                for (std::size_t k = 0; k < n; ++k) {
                    const Vec3 d = map.relativeToParent(sys, parent, trails.point(i, k)) - eye;
                    trailScratch.insert(trailScratch.end(), {float(d.x), float(d.y), float(d.z), float(k + 1) / float(n + 1)});
                }
                const Vec3 cur = map.position(sys, i) - eye;
                trailScratch.insert(trailScratch.end(), {float(cur.x), float(cur.y), float(cur.z), 1.0f});

                Color c = Color::fromHex(sys.info(i).color);
                c = c * (1.0f / std::max(0.05f, c.maxComponent()));
                glUniform3f(locLineColor, c.r, c.g, c.b);
                glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(trailScratch.size() * sizeof(float)), trailScratch.data(),
                             GL_DYNAMIC_DRAW);
                glDrawArrays(GL_LINE_STRIP, 0, GLsizei(trailScratch.size() / 4));
            }
            glDepthMask(GL_TRUE);
            glDisable(GL_BLEND);
        }

        // ---- resolve: multisampled renderbuffer -> plain texture -------------------------
        glBindFramebuffer(GL_READ_FRAMEBUFFER, targets.msaaFbo);
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, targets.resolveFbo);
        glBlitFramebuffer(0, 0, w, h, 0, 0, w, h, GL_COLOR_BUFFER_BIT, GL_NEAREST);

        // ---- bloom: bright pass, then separable blur, all at half resolution ---------------
        glDisable(GL_DEPTH_TEST);
        glBindVertexArray(emptyVao);
        if (bloomOn) {
            glViewport(0, 0, targets.bw, targets.bh);

            glBindFramebuffer(GL_FRAMEBUFFER, targets.bloomFbo[0]);
            glUseProgram(brightProg);
            glUniform1f(locBrightExposure, exposure);
            glUniform1f(locBrightThreshold, bloomThreshold);
            glUniform1f(locBrightKnee, bloomKnee);
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, targets.resolveTex);
            glDrawArrays(GL_TRIANGLES, 0, 3);

            glUseProgram(blurProg);
            for (int i = 0; i < blurIterations; ++i) {
                glBindFramebuffer(GL_FRAMEBUFFER, targets.bloomFbo[1]);      // horizontal: [0] -> [1]
                glBindTexture(GL_TEXTURE_2D, targets.bloomTex[0]);
                glUniform2f(locBlurDir, 1.0f / float(targets.bw), 0.0f);
                glDrawArrays(GL_TRIANGLES, 0, 3);

                glBindFramebuffer(GL_FRAMEBUFFER, targets.bloomFbo[0]);      // vertical: [1] -> [0]
                glBindTexture(GL_TEXTURE_2D, targets.bloomTex[1]);
                glUniform2f(locBlurDir, 0.0f, 1.0f / float(targets.bh));
                glDrawArrays(GL_TRIANGLES, 0, 3);
            }
        }

        // ---- post pass: scene + bloom, exposure, tone map, gamma -> the real screen --------
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glViewport(0, 0, w, h);
        glUseProgram(postProg);
        glUniform1f(locPostExposure, exposure);
        glUniform1f(locPostBloom, bloomOn ? bloomStrength : 0.0f);
        glUniform1i(locPostTonemap, tonemap);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, targets.bloomTex[0]);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, targets.resolveTex);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glEnable(GL_DEPTH_TEST);

        // ---- 4. legend: screen-space overlay, no depth, alpha blended ----------------------
        if (showLegend) {
            const float pixelScale = float(w) / float(std::max(1, winW));
            const float s = std::max(1.0f, std::round(2.0f * pixelScale));
            const float margin = 12.0f * pixelScale;
            glDisable(GL_DEPTH_TEST);
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
            glUseProgram(uiProg);
            glUniform4f(locUiRect, margin, margin, float(legend.w) * s, float(legend.h) * s);
            glUniform2f(locUiViewport, float(w), float(h));
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, legend.id);
            glBindVertexArray(emptyVao);
            glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
            glDisable(GL_BLEND);
            glEnable(GL_DEPTH_TEST);
        }

        SDL_GL_SwapWindow(win);

        // ---- window title as a tiny HUD -----------------------------------------------
        titleTimer += dt;
        if (titleTimer > 0.2) {
            titleTimer = 0.0;
            char title[256];
            std::snprintf(title, sizeof title, "SolarSim GL | %s UTC | %s%.3g d/s | focus %s | scale %s | %s flux | %.0f fps | exposure %.2f %s | bloom %.2f | H: controls",
                          solsim::astro::formatUTC(sim.time(), false).c_str(), direction < 0 ? "-" : "", warp,
                          sys.info(std::size_t(focus)).name.c_str(), solsim::render::scaleModeName(map.mode()),
                          physicalFlux ? "physical" : "display", fps, double(exposure), kTonemapNames[tonemap],
                          bloomOn ? double(bloomStrength) : 0.0);
            SDL_SetWindowTitle(win, paused ? (std::string(title) + " | PAUSED").c_str() : title);
        }
        if (maxFrames >= 0 && ++frame >= maxFrames) running = false;
    }

    glDeleteBuffers(1, &sphereVbo);
    glDeleteBuffers(1, &sphereEbo);
    glDeleteBuffers(1, &starVbo);
    glDeleteBuffers(1, &trailVbo);
    glDeleteVertexArrays(1, &sphereVao);
    glDeleteVertexArrays(1, &starVao);
    glDeleteVertexArrays(1, &trailVao);
    glDeleteTextures(1, &legend.id);
    glDeleteVertexArrays(1, &emptyVao);
    glDeleteProgram(uiProg);
    destroyTargets(targets);
    glDeleteProgram(postProg);
    glDeleteProgram(brightProg);
    glDeleteProgram(blurProg);
    glDeleteProgram(bodyProg);
    glDeleteProgram(lineProg);
    glDeleteProgram(starProg);
    SDL_GL_DeleteContext(ctx);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}