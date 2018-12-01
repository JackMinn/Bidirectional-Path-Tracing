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
 * SSAO (Screen Space Ambient Occlusion) renderpass.
 */
struct SSAOPass : RenderPass {
    GLuint gbuffer;
    GLuint texturePosition;
    GLuint textureNormal;
    GLuint textureDepth;
    GLuint geometryShader;

    GLuint quadVBO;
    GLuint quadVAO;
    GLuint shaderSSAO;

    GLuint uvAttrib{1};

    explicit SSAOPass(const Scene& scene) : RenderPass(scene) { }

    virtual bool init(const Config& config) override {
        RenderPass::init(config);

        // Create vertex buffers
        const auto& shapes = scene.worldData.shapes;
        objects.resize(shapes.size());
        for (size_t i = 0; i < shapes.size(); i++) {
            const tinyobj::shape_t& s = shapes[i];
            GLObject& obj = objects[i];
            buildVBO(i);
            buildVAO(i);
        }

        // Create shader to build GBuffer
        GLuint vs = compileShader("geometry.vs", GL_VERTEX_SHADER);
        GLuint fs = compileShader("geometry.fs", GL_FRAGMENT_SHADER);
        geometryShader = compileProgram(vs, fs);
        glDeleteShader(vs);
        glDeleteShader(fs);

        // Create position texture (GBuffer)
        glGenTextures(1, &texturePosition);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texturePosition);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, config.width, config.height, 0, GL_RGB, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);

        // Create normal texture (GBuffer)
        glGenTextures(1, &textureNormal);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, textureNormal);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, config.width, config.height, 0, GL_RGB, GL_FLOAT, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);

        // Create depth texture (GBuffer)
        glGenTextures(1, &textureDepth);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, textureDepth);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, config.width, config.height, 0, GL_DEPTH_COMPONENT, GL_UNSIGNED_BYTE, NULL);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);

        // Create the GBuffer
        glGenFramebuffers(1, &gbuffer);
        glBindFramebuffer(GL_FRAMEBUFFER, gbuffer);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texturePosition, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, textureNormal, 0);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, textureDepth, 0);
        unsigned int attachments[2] = {GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1};
        glDrawBuffers(2, attachments);
        GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
        if (status != GL_FRAMEBUFFER_COMPLETE) return false;
        glBindFramebuffer(GL_FRAMEBUFFER, 0);

        // Create quad VBO
        glGenVertexArrays(1, &quadVAO);
        glBindVertexArray(quadVAO);
        GLfloat quadVertices[] = {
            //x, y     u, v
            -1, 1, 0, 1,
            -1, -1, 0, 0,
            1, -1, 1, 0,
            -1, 1, 0, 1,
            1, -1, 1, 0,
            1, 1, 1, 1
        };
        glGenBuffers(1, &quadVBO);
        glBindBuffer(GL_ARRAY_BUFFER, quadVBO);
        glBufferData(GL_ARRAY_BUFFER, 4 * 6 * sizeof(GLfloat), (GLvoid*) (&quadVertices[0]), GL_STATIC_DRAW);
        glEnableVertexAttribArray(posAttrib);
        glEnableVertexAttribArray(uvAttrib);
        glVertexAttribPointer(posAttrib, 2, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 4, (GLvoid*)(0 * sizeof(GLfloat)));
        glVertexAttribPointer(uvAttrib, 2, GL_FLOAT, GL_FALSE, sizeof(GLfloat) * 4, (GLvoid*)(2 * sizeof(GLfloat)));

        // Create SSAO shader
        {
            GLuint vs = compileShader("quad.vs", GL_VERTEX_SHADER);
            GLuint fs = compileShader("ssao.fs", GL_FRAGMENT_SHADER);
            shaderSSAO = compileProgram(vs, fs);
            glDeleteShader(vs);
            glDeleteShader(fs);
        }

        return true;
    }

    virtual void cleanUp() override {
        // Delete GBuffer
        glDeleteTextures(1, &texturePosition);
        glDeleteTextures(1, &textureNormal);
        glDeleteTextures(1, &textureDepth);
        glDeleteFramebuffers(1, &gbuffer);
        glDeleteProgram(geometryShader);

        // Delete SSAO shader
        glDeleteBuffers(1, &quadVBO);
        glDeleteVertexArrays(1, &quadVAO);
        glDeleteProgram(shaderSSAO);

        // Delete vertex buffers
        for (size_t i = 0; i < objects.size(); i++) {
            GLObject obj = objects[i];
            glDeleteBuffers(1, &obj.vbo);
            glDeleteVertexArrays(1, &obj.vao);
        }

        RenderPass::cleanUp();
    }

    virtual void render() override {
        // 1 - Geometry pass (GBuffer)
        // =======================================================================================

        /**
        * 1) Bind the GBuffer.
        */
        // TODO: Add previous assignment code (if needed)
        
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
        
        // Update camera
        glm::mat4 model, view, projection;
        camera.Update();
        camera.GetMatricies(projection, view, model);

        /**
        * 1) Use the shader for the geometry pass.
        * 2) Pass the necessary uniforms.
        * 3) Bind vertex array of current object.
        * 4) Draw its triangles.
        * 5) Unbind the vertex array.
        */
        // TODO: Add previous assignment code (if needed)
        
        // 2 - SSAO pass
        // =======================================================================================
        /**
        * 1) Bind the screen buffer (postprocess_fboScreen).
        */
        // TODO: Add previous assignment code (if needed)
        
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        glDisable(GL_DEPTH_TEST);

        /**
        * 1) Use the shader for the SSAO pass.
        * 2) Pass the necessary uniforms.
        * 3) Bind the textures for position and normal from the GBuffer.
        * 4) Bind vertex array of the quad representing the screen texture.
        * 5) Draw the quad.
        * 6) Unbind the vertex array.
        * 7) Unbind the textures.
        */
        // TODO: Add previous assignment code (if needed)
        
        RenderPass::render();
    }

};

TR_NAMESPACE_END
