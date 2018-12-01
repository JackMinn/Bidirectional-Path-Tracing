/*
	This file is part of TinyRender, an educative rendering system.

	All code in this file was written by Jack Minnetian.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Perfectly mirror reflectance model
 */
	struct MirrorBSDF : BSDF {

	MirrorBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
		const tinyobj::material_t& mat = scene.materials[matID];
		std::cout << "ive been made" << std::endl;

		components.push_back(EDeltaReflection);

		combinedType = 0;
		for (size_t i = 0; i < components.size(); ++i)
			combinedType |= components[i];
	}

	inline v3f reflect(const v3f& d) const {
		return v3f(-d.x, -d.y, d.z);
	}

	v3f eval(const SurfaceInteraction& i) const override {
		v3f val(0.f);

		//It is impossible to have the reflection direction in advance, so we return 0
		//We do not try to compare the current direction to the reflection direction due to floating point imprecision
		return val;
	}

	float pdf(const SurfaceInteraction& i) const override {
		float pdf = 0.f;
		
		//It is a delta distribution, so the pdf returns 0
		return pdf;
	}

	//I dont think all the light is reflected, but only part of it (fresnel equations). Which would imply that it is not energy conserving either. 
	v3f sample(SurfaceInteraction& i, const v2f& sample, float* pdf) const override {
		v3f val(0.f);
		*pdf = 1.f;
		
		i.wi = reflect(i.wo);

		//brdf is delta / cosTheta, and we need to return brdfCosTheta, so the cosThetas cancel and we return 1
		val = v3f(1.f); 

		return val;
	}

	std::string toString() const override { return "Mirror"; }
};

TR_NAMESPACE_END