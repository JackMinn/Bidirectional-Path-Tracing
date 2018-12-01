/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/renderpass.h"
#include "tiny_obj_loader.h"

TR_NAMESPACE_BEGIN

/**
 * Surface normal renderpass.
 */
struct NormalPass : RenderPass {
    GLuint shader{0};

    GLuint modelMatUniform{0};
    GLuint viewMatUniform{0};
    GLuint projectionMatUniform{0};
    GLuint normalMatUniform{0};

    explicit NormalPass(const Scene& scene) : RenderPass(scene) { }

    bool init(const Config& config) override {
        RenderPass::init(config);

        // Create shader
        GLuint vs = compileShader("normal.vs", GL_VERTEX_SHADER);
        GLuint fs = compileShader("normal.fs", GL_FRAGMENT_SHADER);
        shader = compileProgram(vs, fs);
        glDeleteShader(vs);
        glDeleteShader(fs);

        // Create uniforms
        modelMatUniform = GLuint(glGetUniformLocation(shader, "model"));
        viewMatUniform = GLuint(glGetUniformLocation(shader, "view"));
        projectionMatUniform = GLuint(glGetUniformLocation(shader, "projection"));
        normalMatUniform = GLuint(glGetUniformLocation(shader, "normalMat"));

        // Create vertex buffers
        objects.resize(scene.worldData.shapes.size());
        for (size_t i = 0; i < objects.size(); i++) {
            buildVBO(i);
            buildVAO(i);
        }

        return true;
    }

    void cleanUp() override {
        // Delete vertex buffers
        for (size_t i = 0; i < objects.size(); i++) {
            glDeleteBuffers(1, &objects[i].vbo);
            glDeleteVertexArrays(1, &objects[i].vao);
        }

        RenderPass::cleanUp();
    }

    void render() override {
        glBindFramebuffer(GL_FRAMEBUFFER, postprocess_fboScreen);
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);

        // Define shader to use
        glUseProgram(shader);

        // Update camera
        glm::mat4 model, view, projection;
        camera.Update();
        camera.GetMatricies(projection, view, model);

        // Pass uniforms
        glUniformMatrix4fv(modelMatUniform, 1, GL_FALSE, &(modelMat[0][0]));
        glUniformMatrix4fv(viewMatUniform, 1, GL_FALSE, &(view[0][0]));
        glUniformMatrix4fv(projectionMatUniform, 1, GL_FALSE, &(projection[0][0]));
        glUniformMatrix4fv(normalMatUniform, 1, GL_FALSE, &(normalMat[0][0]));

        // Draw
        for (auto& object : objects) {
            /**
             * 1) Bind vertex array of current object.
             * 2) Draw its triangles.
             * 3) Bind vertex array to 0.
             */
            // TODO: Add previous assignment code (if needed)
        }

        RenderPass::render();
    }

};

TR_NAMESPACE_END
