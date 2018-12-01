/*
	This file is part of TinyRender, an educative rendering system.

	All code in this file was written by Jack Minnetian.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Glass reflectance model with transmission
 */
	struct GlassBSDF : BSDF {
	float m_IndexOfRefraction;
	v3f m_Transmittance;

	GlassBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
		const tinyobj::material_t& mat = scene.materials[matID];
		m_IndexOfRefraction = mat.ior;
		m_Transmittance = glm::make_vec3(mat.transmittance);

		std::cout << m_IndexOfRefraction << std::endl;
		std::cout << getLuminance(m_Transmittance) << std::endl;

		components.push_back(EDeltaReflection);
		components.push_back(EDeltaTransmission);

		combinedType = 0;
		for (size_t i = 0; i < components.size(); ++i)
			combinedType |= components[i];
	}

	inline v3f reflect(const v3f& d) const {
		return v3f(-d.x, -d.y, d.z);
	}

	float FresnelDielectric(float eta_i, float eta_t, float cos_i, float cos_t) const {
		float eta = eta_i / eta_t;
		float sin2_t = eta * eta * (std::fmax(0.f, 1.f - cos_i * cos_i));
		// Total Internal Reflection
		if (sin2_t >= 1.f) 
			return 1.f;

		float rParallel =		((eta_t * cos_i) - (eta_i * cos_t)) /
								((eta_t * cos_i) + (eta_i * cos_t));
		float rPerpendicular =  ((eta_i * cos_i) - (eta_t * cos_t)) /
								((eta_i * cos_i) + (eta_t * cos_t));

		return (rParallel * rParallel + rPerpendicular * rPerpendicular) * 0.5f;
	}

	v3f eval(const SurfaceInteraction& i) const override {
		v3f val(0.f);

		return val;
	}

	float pdf(const SurfaceInteraction& i) const override {
		float pdf = 0.f;
		
		return pdf;
	}

	v3f sample(SurfaceInteraction& i, const v2f& sample, float* pdf) const override {
		v3f val(0.f);
		*pdf = 1;
		// TODO: Implement this
		
		bool isEntering = i.wo.z > 0.f;
		//we assuming the other medium at the surface boundary is air with IOR of 1
		float eta_i = 1.f, eta_t = m_IndexOfRefraction; 
		if (!isEntering)
			std::swap(eta_i, eta_t);

		float eta = eta_i / eta_t;
		float sin2_i = std::fmax(0.f, 1.f - i.wo.z * i.wo.z);
		float sin2_t = eta * eta * sin2_i;
		float cos_t = std::sqrtf(std::fmax(0.f, 1.f - sin2_t));
		cos_t = isEntering ? -cos_t : cos_t;

		//Total internal reflection. Light cannot pass the surface boundary between the 2 mediums
		//from the current incident direction. If I am entering the more optically dense medium from 
		//the less optically dense medium then I can always find a transmission ray, however as I 
		//approach 90 degrees, less is transmitted and more is reflected based on fresnel equations 
		//if (sin2_t >= 1.f) 
		//	return v3f(0.f); //should we not bounce the total internal reflection? I understand that we assuming this is a transmission sample, then there the result is 0, but what if we assume its a reflection sample?

		float fresnel = FresnelDielectric(eta_i, eta_t, std::fabs(i.wo.z), std::fabs(cos_t));
		
		//This is effectively russian rouletting, since we assume a pdf of 1
		//If we have total internal reflection, fresnel is always 1, and so
		//we never do a transmittance sample but always a reflectance sample.
		if (sample.x < fresnel) {
			val = v3f(1.f); //might make this a material property after
			i.wi = reflect(i.wo);
			//*pdf = fresnel;
		}
		else {
			val = m_Transmittance;
			i.wi = v3f(eta * -i.wo.x, eta * -i.wo.y, cos_t);
			//*pdf = 1.f - fresnel;
		}

		return val;
	}

	std::string toString() const override { return "Glass"; }
};

TR_NAMESPACE_END