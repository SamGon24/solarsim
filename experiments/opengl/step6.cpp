// =============================================================================
//  Step 6: the real simulation drawn with OpenGL.
//  Reuses SolarSim's physics (Simulation), display mapping and orbit camera.
//  Controls:
//    left mouse drag  orbit        mouse wheel  zoom
//    Tab / 0-9        focus body   Z            reset camera
//    Space            pause        [ ]          slower / faster time
//    R                reverse time F            display flux <-> physical flux
//    M                radial scale (sqrt / log / linear)
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
#include "render/DisplayMapping.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

using solsim::Vec3;

// -----------------------------------------------------------------------------
//  Shaders
// -----------------------------------------------------------------------------
static const char* kVert = R"(#version 330 core
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

static const char* kFrag = R"(#version 330 core
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
        lin = uColor * 3.0;
    } else {
        vec3 N = normalize(vN);
        vec3 L = normalize(uLightPos - vWorld);
        float diffuse = max(dot(N, L), 0.0);
        lin = uColor * (0.02 + diffuse * uFlux * 2.0);
    }
    lin = lin / (1.0 + lin);
    FragColor = vec4(pow(lin, vec3(1.0 / 2.2)), 1.0);
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

// Rotation-only view matrix built from the camera basis. The camera position
// is NOT in this matrix: every object is positioned relative to the eye
// (camera-relative rendering), so large offsets never reach the GPU as floats.
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

// -----------------------------------------------------------------------------
//  Sphere mesh
// -----------------------------------------------------------------------------
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
//  Helpers around the SolarSim classes
// -----------------------------------------------------------------------------
// sRGB 0xRRGGBB -> linear floats (lighting maths is done in linear space).
static void linearColor(std::uint32_t c, float out[3]) {
    out[0] = std::pow(float((c >> 16) & 255) / 255.0f, 2.2f);
    out[1] = std::pow(float((c >> 8) & 255) / 255.0f, 2.2f);
    out[2] = std::pow(float(c & 255) / 255.0f, 2.2f);
}

// Radius of the planetary system in display units (farthest body orbiting the Sun).
static double systemExtent(const solsim::render::DisplayMapping& map, const solsim::physics::NBodySystem& sys) {
    double outer = 1.0;
    for (std::size_t i = 1; i < sys.namedCount(); ++i) {
        if (sys.info(i).parent != 0) continue;
        outer = std::max(outer, map.position(sys, i).norm());
    }
    return outer;
}

// Camera distance that frames a body and its satellites (same idea as `solarsim view`).
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
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4);
    SDL_Window* win = SDL_CreateWindow("SolarSim GL", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1280, 720,
                                       SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_GLContext ctx = win ? SDL_GL_CreateContext(win) : nullptr;
    if (!ctx) {
        std::printf("GL: %s\n", SDL_GetError());
        return 1;
    }
    SDL_GL_SetSwapInterval(1);
    std::printf("GL_VERSION:  %s\nGL_RENDERER: %s\n", glGetString(GL_VERSION), glGetString(GL_RENDERER));

    GLuint prog = makeProgram(kVert, kFrag);
    GLint locMVP = glGetUniformLocation(prog, "uMVP");
    GLint locModel = glGetUniformLocation(prog, "uModel");
    GLint locColor = glGetUniformLocation(prog, "uColor");
    GLint locLight = glGetUniformLocation(prog, "uLightPos");
    GLint locFlux = glGetUniformLocation(prog, "uFlux");
    GLint locEmissive = glGetUniformLocation(prog, "uEmissive");

    std::vector<float> verts;
    std::vector<unsigned int> idx;
    buildSphere(48, 96, verts, idx);

    GLuint vao = 0, vbo = 0, ebo = 0;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ebo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, GLsizeiptr(verts.size() * sizeof(float)), verts.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, GLsizeiptr(idx.size() * sizeof(unsigned int)), idx.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_MULTISAMPLE);

    // ---- state ------------------------------------------------------------------------
    int focus = 0;
    bool paused = false;
    bool physicalFlux = false;
    double warp = 10.0;      // simulated days per real second
    double direction = 1.0;  // +1 forward, -1 backward
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
        if (!paused) sim.advance(direction * warp * dt);

        // ---- camera follows the focused body ------------------------------------------
        map.prepare(sys);
        cam.setTarget(focus > 0 ? map.position(sys, std::size_t(focus)) : Vec3{});

        int w = 1, h = 1;
        SDL_GL_GetDrawableSize(win, &w, &h);
        cam.setViewport(w, h);
        glViewport(0, 0, w, h);
        glClearColor(0.004f, 0.004f, 0.012f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        const Vec3 eye = cam.eye();
        const double dist = cam.distance();
        const float nearZ = float(std::max(1e-6, dist * 0.001));
        const float farZ = float(dist * 20.0 + 50.0);
        const Mat4 proj = perspective(float(cam.fieldOfView()), float(w) / float(h), nearZ, farZ);
        const Mat4 viewProj = mul(proj, viewRotation(cam.right(), cam.up(), cam.forward()));
        const double focalPx = 0.5 * h / std::tan(0.5 * cam.fieldOfView());

        glUseProgram(prog);
        glBindVertexArray(vao);
        const Vec3 sunRel = map.position(sys, 0) - eye;
        glUniform3f(locLight, float(sunRel.x), float(sunRel.y), float(sunRel.z));

        // ---- draw every named body ----------------------------------------------------
        for (std::size_t i = 0; i < sys.namedCount(); ++i) {
            const Vec3 p = map.position(sys, i);
            const Vec3 rel = p - eye;                       // camera-relative, still in double
            const double depth = rel.dot(cam.forward());
            if (depth <= nearZ) continue;
            // Bodies are tiny at overview scale: never let them shrink below ~2.5 pixels.
            const double radius = std::max(map.displayRadius(sys, i), 2.5 * depth / focalPx);

            double flux = 1.0;
            if (i > 0) {
                const double d = physicalFlux ? (sys.pos(i) - sys.pos(0)).norm() : p.norm();
                flux = 1.0 / std::max(d * d, 1e-6);
            }
            float col[3];
            linearColor(sys.info(i).color, col);

            const Mat4 model = mul(translate(float(rel.x), float(rel.y), float(rel.z)), scaleUniform(float(radius)));
            const Mat4 mvp = mul(viewProj, model);
            glUniformMatrix4fv(locMVP, 1, GL_FALSE, mvp.m);
            glUniformMatrix4fv(locModel, 1, GL_FALSE, model.m);
            glUniform3f(locColor, col[0], col[1], col[2]);
            glUniform1f(locFlux, float(flux));
            glUniform1f(locEmissive, i == 0 ? 1.0f : 0.0f);
            glDrawElements(GL_TRIANGLES, GLsizei(idx.size()), GL_UNSIGNED_INT, (void*)0);
        }

        SDL_GL_SwapWindow(win);

        // ---- window title as a tiny HUD -----------------------------------------------
        titleTimer += dt;
        if (titleTimer > 0.2) {
            titleTimer = 0.0;
            char title[256];
            std::snprintf(title, sizeof title, "SolarSim GL | %s UTC | %s%.3g d/s | focus %s | scale %s | %s flux | %.0f fps",
                          solsim::astro::formatUTC(sim.time(), false).c_str(), direction < 0 ? "-" : "", warp,
                          sys.info(std::size_t(focus)).name.c_str(), solsim::render::scaleModeName(map.mode()),
                          physicalFlux ? "physical" : "display", fps);
            SDL_SetWindowTitle(win, paused ? (std::string(title) + " | PAUSED").c_str() : title);
        }
        if (maxFrames >= 0 && ++frame >= maxFrames) running = false;
    }

    glDeleteBuffers(1, &vbo);
    glDeleteBuffers(1, &ebo);
    glDeleteVertexArrays(1, &vao);
    glDeleteProgram(prog);
    SDL_GL_DeleteContext(ctx);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}