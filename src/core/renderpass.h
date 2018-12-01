/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/platform.h"
#include <GL/glew.h>

#ifdef __APPLE__
#include "SDL.h"
#include <OpenGL/gl.h>
#else
#ifdef _WIN32
#include <GL/gl.h>
#include "SDL.h"
#else
#include <GL/gl.h>
#include "SDL2/SDL.h"
#endif
#endif

#include <core/core.h>

TR_NAMESPACE_BEGIN

struct Config;
struct BSDF;

/**
 * GL Object structure.
 */
struct GLObject {
    GLuint vao{0};
    GLuint vbo{0};
    GLuint shaderID{0};
    int shaderIdx{0};

    int nVerts{0};
    std::vector<GLfloat> vertices;

    // uniforms for diffuse shader
    GLuint albedoUniform{0};
    v3f albedo{0};

    // uniforms for phong shader
    GLuint exponentUniform{0};
    GLuint rho_d_Uniform{0};
    GLuint rho_s_Uniform{0};
    float exponent{0};
    v3f rho_d{0};
    v3f rho_s{0};
};

/**
 * RenderPass structure.
 */
struct RenderPass {
    const Scene& scene;

    SDL_Window* window{nullptr};
    SDL_GLContext contextGL{0};

    int width, height;
    int nPixel;

    bool isSaved = false;

    std::string shadersFilePath;

    std::vector<GLObject> objects;

    // BSDFs shaders
    static const int N_SHADERS = 3;
    GLuint shaders[N_SHADERS];
    const char* shadersName[N_SHADERS] = {"diffuse", "phong", "emitter"};
    const int DIFFUSE_SHADER_IDX = 0;
    const int PHONG_SHADER_IDX = 1;
    const int EMITTER_SHADER_IDX = 2;

    // Matrices
    glm::mat4 modelMat;
    glm::mat4 normalMat;

    // Camera real-time
    CameraRT camera;

    // Camera
    glm::vec3 camPos;

    // Light
    glm::vec3 lightPos;
    glm::vec3 lightIntensity;

    // shader attributes
    GLuint posAttrib{0};
    GLuint normalAttrib{1};
    const int N_ATTR_PER_VERT{6}; // 3 for pos, 3 for normals

    explicit RenderPass(const Scene& scene) : scene(scene) { }
    virtual bool init(const Config& config);
    virtual void cleanUp();
    virtual void render();

    virtual void buildVBO(size_t objectIdx);
    virtual void buildVAO(size_t objectIdx);

    bool save(GLfloat* data);
    void updateCamera(SDL_Event& e);

    // Utils
    bool initOpenGL(int width, int height);
    GLuint compileProgram(GLuint vs, GLuint fs);
    GLuint compileShader(const char* shaderPath_, GLenum shaderType);
    GLuint compileShader_(const char* codePtr, GLenum shaderType);
    std::string readFile(const char* filePath);
    void assignShader(GLObject& obj, const tinyobj::shape_t& s, const std::vector<std::unique_ptr<BSDF>>& bsdfs);

    // For Linear->sRGB post-process shader
    GLuint postprocess_quadShader;
    GLuint postprocess_quadVAO;
    GLuint postprocess_quadVBO;

    GLuint postprocess_fboScreen;
    GLuint postprocess_textureColor;
    GLuint postprocess_textureDepth;

    GLuint uvAttrib{1};

    void initPostProcessShader();
    void renderPostProcessShader();
};

TR_NAMESPACE_END