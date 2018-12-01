/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <core/platform.h>
#include <core/core.h>
#include <core/integrator.h>
#include <core/renderpass.h>

TR_NAMESPACE_BEGIN

/**
 * Renderer structure (offline and real-time).
 */
struct Renderer {
    std::unique_ptr<Integrator> integrator;
    std::unique_ptr<RenderPass> renderpass;
    Scene scene;
    bool realTime;
    bool nogui;
    bool realTimeCameraFree;
    unsigned int previousTime = 0, currentTime = 0;
    const int frameDuration = 30;

    explicit Renderer(const Config& config);
    bool init(bool isRealTime, bool nogui);
    void render();
    void cleanUp();

	//this function here is just for unit testing, its redundant and exists in bdpt.h
	inline void splatToImagePlane(const v3f position, int& xPixel, int& yPixel) const {
		float aspectRatio = (float)integrator->rgb->width / integrator->rgb->height;
		auto worldToCamera = glm::lookAt(scene.config.camera.o, scene.config.camera.at, scene.config.camera.up);
		auto cameraToClip = glm::perspective(deg2rad * scene.config.camera.fov, aspectRatio, 1.f, 1000.f);
		auto NDCToScreen = glm::scale(glm::mat4(1.f), v3f(integrator->rgb->width, integrator->rgb->height, 1.0f)) * glm::scale(glm::mat4(1.f), v3f(0.5f, -0.5f, 1.f)) * glm::translate(glm::mat4(1.f), v3f(1.f, -1.f, 0.f));
		v4f uv = worldToCamera * v4f(position, 1.f);
		uv = cameraToClip * uv;
		uv /= uv.w;
		uv = NDCToScreen * uv;
		xPixel = static_cast<int>(uv.x);
		yPixel = static_cast<int>(uv.y);
	}
};

TR_NAMESPACE_END