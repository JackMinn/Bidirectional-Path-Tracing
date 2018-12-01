/*
	This file is part of TinyRender, an educative rendering system.

	All code in this file was written by Jack Minnetian.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Dielectric Microfacet model using GGX 
 */
	struct MicrofacetBSDF : BSDF {
	float m_PerceptualRoughness;

	MicrofacetBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
		const tinyobj::material_t& mat = scene.materials[matID];
		m_PerceptualRoughness = mat.roughness; //roughness should be between 0 and 1, what would it mean to go higher than 1?

		components.push_back(BSDF::EDiffuseReflection);
		components.push_back(BSDF::EGlossyReflection);

		combinedType = 0;
		for (size_t i = 0; i < components.size(); ++i)
			combinedType |= components[i];
	}

	float FresnelSchlick(float R0, float VDotH)
	{
		return R0 + (1.f - R0) * std::powf(1.f - VDotH, 5.f);
	}

	//https://agraphicsguy.wordpress.com/2015/11/01/sampling-microfacet-brdf/
	inline v3f sampleGGX(const float& roughness, const v2f& sample) const {
		float theta = glm::atan(roughness * glm::sqrt(sample.x / (1 - sample.x)));
		float phi = M_PI * 2.f * sample.y;

		float cosTheta = std::cosf(theta);
		float sinTheta = std::sqrtf(std::fmax(1.f - cosTheta * cosTheta, 0));
		v3f v = v3f(sinTheta * std::cosf(phi), sinTheta * std::sinf(phi), cosTheta);
		return v;
	}

	//the pdf with respect to a measure of solid angle around the half vector, v is the half vector
	inline float pdfGGX(const float& roughness, const v3f& v) const {
		float roughnessSqr = roughness * roughness;
		float cosTheta = v.z;
		float cosThetaSqr = cosTheta * cosTheta;

		//If roughness = 1, the larger cosTheta, the higher chance of sampling it, so we must sample smaller theta more often
		float pdf = (roughnessSqr * cosTheta) / (M_PI * std::powf(((roughnessSqr - 1) * cosThetaSqr + 1), 2.f));
	}

	v3f eval(const SurfaceInteraction& i) const override {
		v3f val(0.f);

		return val;
	}

	float pdf(const SurfaceInteraction& i) const override {
		float pdf = 0.f;
		float academicRoughness = m_PerceptualRoughness * m_PerceptualRoughness;

		v3f halfVector = glm::normalize(i.wi + i.wo);
		pdf = pdfGGX(academicRoughness, halfVector);
		float halfVectorToIncidentDirection = 1.f / (4.f * glm::dot(halfVector, i.wi));
		pdf *= halfVectorToIncidentDirection;

		return pdf;
	}

	//Get a half vector, using that you can decide what the new direction is to bounce in
	v3f sample(SurfaceInteraction& i, const v2f& sample, float* pdf) const override {
		v3f val(0.f);
		float academicRoughness = m_PerceptualRoughness * m_PerceptualRoughness;
		v3f halfVector = sampleGGX(academicRoughness, sample); //this is in shading tangent space
		v3f wi = i.wo - halfVector;
		wi = glm::normalize(wi);

		v3f halfVectorTest = glm::normalize(wi + i.wo);
		//check that halfVectorTest and halfVector are equal to each other

		i.wi = wi;



		val = eval(i);

		return val;
	}

	std::string toString() const override { return "Dielectric Microfacet"; }
};

TR_NAMESPACE_END