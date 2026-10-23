#define GL_SILENCE_DEPRECATION
#include <SDL.h>
#include <OpenGL/gl3.h>
#include <cmath>
#include <cstdio>
#include <vector>

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
uniform float uEmissive;
out vec4 FragColor;
void main() {
    vec3 lin;
    if (uEmissive > 0.5) {
        lin = uColor;
    } else {
        vec3 N = normalize(vN);
        vec3 toLight = uLightPos - vWorld;
        float d2 = dot(toLight, toLight);
        vec3 L = toLight * inversesqrt(d2);
        float diffuse = max(dot(N, L), 0.0);
        float light = 0.03 + diffuse * 9.0 / d2;
        lin = uColor * light;
    }
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

struct Body {
    float orbitRadius;
    float speed;
    float size;
    float r, g, b;
};

int main() {
    if (SDL_Init(SDL_INIT_VIDEO) != 0) return 1;
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_CORE);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_FLAGS, SDL_GL_CONTEXT_FORWARD_COMPATIBLE_FLAG);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 24);
    SDL_Window* win = SDL_CreateWindow("Step 5", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                       1280, 720, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE | SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_GLContext ctx = SDL_GL_CreateContext(win);
    if (!ctx) return 1;
    SDL_GL_SetSwapInterval(1);

    GLuint prog = makeProgram(kVert, kFrag);
    GLint locMVP = glGetUniformLocation(prog, "uMVP");
    GLint locModel = glGetUniformLocation(prog, "uModel");
    GLint locColor = glGetUniformLocation(prog, "uColor");
    GLint locLight = glGetUniformLocation(prog, "uLightPos");
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
    glBufferData(GL_ARRAY_BUFFER, verts.size() * sizeof(float), verts.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, idx.size() * sizeof(unsigned int), idx.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    const Body planets[] = {
        {3.0f, 1.6f, 0.30f, 0.7f, 0.7f, 0.7f},
        {5.0f, 1.0f, 0.45f, 0.2f, 0.5f, 1.0f},
        {7.5f, 0.6f, 0.38f, 1.0f, 0.4f, 0.2f},
    };

    glEnable(GL_DEPTH_TEST);

    auto draw = [&](const Mat4& viewProj, const Mat4& model, float r, float g, float b, float emissive) {
        Mat4 mvp = mul(viewProj, model);
        glUniformMatrix4fv(locMVP, 1, GL_FALSE, mvp.m);
        glUniformMatrix4fv(locModel, 1, GL_FALSE, model.m);
        glUniform3f(locColor, r, g, b);
        glUniform1f(locEmissive, emissive);
        glDrawElements(GL_TRIANGLES, (GLsizei)idx.size(), GL_UNSIGNED_INT, (void*)0);
    };

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
        glClearColor(0.005f, 0.005f, 0.02f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        float t = SDL_GetTicks() / 1000.0f;
        V3 eye = {13.0f * std::sin(t * 0.15f), 7.0f, 13.0f * std::cos(t * 0.15f)};
        Mat4 view = lookAt(eye, {0.0f, 0.0f, 0.0f}, {0.0f, 1.0f, 0.0f});
        Mat4 proj = perspective(0.9f, float(w) / float(h), 0.1f, 200.0f);
        Mat4 viewProj = mul(proj, view);

        glUseProgram(prog);
        glUniform3f(locLight, 0.0f, 0.0f, 0.0f);
        glBindVertexArray(vao);

        draw(viewProj, scaleUniform(1.2f), 1.0f, 0.85f, 0.3f, 1.0f);

        int n = 0;
        for (const Body& p : planets) {
            float a = t * p.speed;
            float px = p.orbitRadius * std::cos(a);
            float pz = p.orbitRadius * std::sin(a);
            draw(viewProj, mul(translate(px, 0.0f, pz), scaleUniform(p.size)), p.r, p.g, p.b, 0.0f);

            if (n == 1) {
                float ma = t * 3.0f;
                float mx = 0.9f * std::cos(ma);
                float mz = 0.9f * std::sin(ma);
                Mat4 moon = mul(translate(px + mx, 0.0f, pz + mz), scaleUniform(0.12f));
                draw(viewProj, moon, 0.8f, 0.8f, 0.75f, 0.0f);
            }
            ++n;
        }

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