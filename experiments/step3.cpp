#define GL_SILENCE_DEPRECATION
#include <SDL.h>
#include <OpenGL/gl3.h>
#include <cstdio>
#include <cmath>

static const char* kVert = R"(#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aColor;
uniform mat4 uMVP;
out vec3 vColor;
void main() {
    vColor = aColor;
    gl_Position = uMVP * vec4(aPos, 1.0);
}
)";

static const char* kFrag = R"(#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main() {
    FragColor = vec4(vColor, 1.0);
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

struct Mat4 { float m[16]; };
struct V3 { float x, y, z; };

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

static Mat4 rotateY(float a) {
    Mat4 r = identity();
    r.m[0] = std::cos(a);  r.m[8] = std::sin(a);
    r.m[2] = -std::sin(a); r.m[10] = std::cos(a);
    return r;
}

static Mat4 rotateX(float a) {
    Mat4 r = identity();
    r.m[5] = std::cos(a);  r.m[9] = -std::sin(a);
    r.m[6] = std::sin(a);  r.m[10] = std::cos(a);
    return r;
}

static V3 sub(V3 a, V3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
static float dot(V3 a, V3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
static V3 cross(V3 a, V3 b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
static V3 normalize(V3 v) { float l = std::sqrt(dot(v, v)); return {v.x / l, v.y / l, v.z / l}; }

static Mat4 lookAt(V3 eye, V3 center, V3 up) {
    V3 f = normalize(sub(center, eye));
    V3 s = normalize(cross(f, up));
    V3 u = cross(s, f);
    Mat4 r = identity();
    r.m[0] = s.x;  r.m[4] = s.y;  r.m[8]  = s.z;
    r.m[1] = u.x;  r.m[5] = u.y;  r.m[9]  = u.z;
    r.m[2] = -f.x; r.m[6] = -f.y; r.m[10] = -f.z;
    r.m[12] = -dot(s, eye);
    r.m[13] = -dot(u, eye);
    r.m[14] = dot(f, eye);
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

int main() {
    if (SDL_Init(SDL_INIT_VIDEO) != 0) return 1;
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_Window* win = SDL_CreateWindow("Step 3", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                       1280, 720, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_GLContext ctx = SDL_GL_CreateContext(win);
    if (!ctx) return 1;
    SDL_GL_SetSwapInterval(1);

    GLuint prog = makeProgram(kVert, kFrag);
    GLint locMVP = glGetUniformLocation(prog, "uMVP");

    const float verts[] = {
        -1, -1, -1,  1.0f, 0.2f, 0.2f,
         1, -1, -1,  0.2f, 1.0f, 0.2f,
         1,  1, -1,  0.2f, 0.4f, 1.0f,
        -1,  1, -1,  1.0f, 1.0f, 0.2f,
        -1, -1,  1,  1.0f, 0.2f, 1.0f,
         1, -1,  1,  0.2f, 1.0f, 1.0f,
         1,  1,  1,  1.0f, 1.0f, 1.0f,
        -1,  1,  1,  0.6f, 0.6f, 0.6f,
    };
    const unsigned int idx[] = {
        0, 1, 2,  2, 3, 0,
        4, 5, 6,  6, 7, 4,
        0, 4, 7,  7, 3, 0,
        1, 5, 6,  6, 2, 1,
        0, 1, 5,  5, 4, 0,
        3, 2, 6,  6, 7, 3,
    };

    GLuint vao = 0, vbo = 0, ebo = 0;
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ebo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof verts, verts, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof idx, idx, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);

    glEnable(GL_DEPTH_TEST);

    bool running = true;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) running = false;
            if (e.type == SDL_KEYDOWN && e.key.keysym.sym == SDLK_ESCAPE) running = false;
        }
        int w, h;
        SDL_GL_GetDrawableSize(win, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.01f, 0.01f, 0.04f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        float t = SDL_GetTicks() / 1000.0f;
        Mat4 model = mul(rotateY(t), rotateX(t * 0.6f));
        Mat4 view = lookAt({0.0f, 1.5f, 5.0f}, {0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f});
        Mat4 proj = perspective(1.0f, float(w) / float(h), 0.1f, 100.0f);
        Mat4 mvp = mul(proj, mul(view, model));

        glUseProgram(prog);
        glUniformMatrix4fv(locMVP, 1, GL_FALSE, mvp.m);
        glBindVertexArray(vao);
        glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, (void*)0);

        SDL_GL_SwapWindow(win);
    }
    glDeleteBuffers(1, &vbo);
    glDeleteBuffers(1, &ebo);
    glDeleteVertexArrays(1, &vao);
    glDeleteProgram(prog);
    SDL_GL_DeleteContext(ctx);
    SDL_DestroyWindow(win);
    SDL_Quit();
}