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
 * Simple direct illumination (No shadows) renderpass.
 */
struct SimplePass : RenderPass {

    explicit SimplePass(const Scene& scene) : RenderPass(scene) { }

    virtual bool init(const Config& config) override {
        RenderPass::init(config);

        // Create BSDF shaders
        for (int i = 0; i < N_SHADERS; i++) {
            const std::string name = shadersName[i];
            GLuint vs = compileShader("simple.vs", GL_VERTEX_SHADER);
            GLuint fs = compileShader((name + ".fs").c_str(), GL_FRAGMENT_SHADER);
            shaders[i] = compileProgram(vs, fs);
            glDeleteShader(vs);
            glDeleteShader(fs);
        }

        // Create vertex buffers
        const auto& shapes = scene.worldData.shapes;
        objects.resize(shapes.size());
        for (size_t i = 0; i < objects.size(); i++) {
            buildVBO(i);
            buildVAO(i);
            assignShader(objects[i], shapes[i], scene.bsdfs);
        }

        return true;
    }

    virtual void cleanUp() override {
        // Delete vertex buffers
        for (size_t i = 0; i < objects.size(); i++) {
            glDeleteBuffers(1, &objects[i].vbo);
            glDeleteVertexArrays(1, &objects[i].vao);
        }

        RenderPass::cleanUp();
    }

    virtual void render() override {
        glBindFramebuffer(GL_FRAMEBUFFER, postprocess_fboScreen);
        glClearColor(0.f, 0.f, 0.f, 1.f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);

        // Update camera
        glm::mat4 model, view, projection;
        camera.Update();
        camera.GetMatricies(projection, view, model);

        for (size_t i = 0; i < objects.size(); i++) {
            GLObject obj = objects[i];

            // Define shader to use
            glUseProgram(obj.shaderID);

            // Pass uniforms
            GLuint modelMatUniform = GLuint(glGetUniformLocation(obj.shaderID, "model"));
            GLuint viewMatUniform = GLuint(glGetUniformLocation(obj.shaderID, "view"));
            GLuint projectionMatUniform = GLuint(glGetUniformLocation(obj.shaderID, "projection"));
            GLuint normalMatUniform = GLuint(glGetUniformLocation(obj.shaderID, "normalMat"));
            glUniformMatrix4fv(modelMatUniform, 1, GL_FALSE, &(modelMat[0][0]));
            glUniformMatrix4fv(viewMatUniform, 1, GL_FALSE, &(view[0][0]));
            glUniformMatrix4fv(projectionMatUniform, 1, GL_FALSE, &(projection[0][0]));
            glUniformMatrix4fv(normalMatUniform, 1, GL_FALSE, &(normalMat[0][0]));

            GLuint camPosUniform = GLuint(glGetUniformLocation(obj.shaderID, "camPos"));
            glUniform3f(camPosUniform, camera.camera_position.x, camera.camera_position.y, camera.camera_position.z);

            // Pass light position & power via uniforms
            //   this->lightPos
            //   this->lightIntensity
            // TODO: Add previous assignment code (if needed)

            // Pass shader-specific parameters via uniforms
            //   obj.albedo
            //   obj.rho_d
            //   obj.rho_s
            //   obj.exponent
            // TODO: Add previous assignment code (if needed)

            // Draw
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
