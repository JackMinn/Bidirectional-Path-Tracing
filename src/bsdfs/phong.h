/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core/core.h"

TR_NAMESPACE_BEGIN

/**
 * Modified Phong reflectance model
 */
struct PhongBSDF : BSDF {

    std::unique_ptr<Texture < v3f>> specularReflectance;
    std::unique_ptr<Texture < v3f>> diffuseReflectance;
    std::unique_ptr<Texture < float>> exponent;
    float specularSamplingWeight;
    float scale;

    PhongBSDF(const WorldData& scene, const Config& config, const size_t& matID) : BSDF(scene, config, matID) {
        const tinyobj::material_t& mat = scene.materials[matID];

        if (mat.specular_texname.empty())
            specularReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.specular)));
        else
            specularReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.specular_texname));

        if (mat.diffuse_texname.empty())
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new ConstantTexture3f(glm::make_vec3(mat.diffuse)));
        else
            diffuseReflectance = std::unique_ptr<Texture<v3f>>(new BitmapTexture3f(config, mat.diffuse_texname));

        exponent = std::unique_ptr<Texture<float>>(new ConstantTexture1f(mat.shininess));

        //get scale value to ensure energy conservation
        v3f maxValue = specularReflectance->getMax() + diffuseReflectance->getMax();
        float actualMax = max(max(maxValue.x, maxValue.y), maxValue.z);
        scale = actualMax > 1.0f ? 0.99f * (1.0f / actualMax) : 1.0f;

        float dAvg = getLuminance(diffuseReflectance->getAverage() * scale);
        float sAvg = getLuminance(specularReflectance->getAverage() * scale);
        specularSamplingWeight = sAvg / (dAvg + sAvg);

        components.push_back(EGlossyReflection);
        components.push_back(EDiffuseReflection);

        combinedType = 0;
        for (unsigned int component : components)
            combinedType |= component;
    }

    inline v3f reflect(const v3f& d) const {
        return v3f(-d.x, -d.y, d.z);
    }

	v3f eval(const SurfaceInteraction& i) const override {
		v3f val(0.f);
		// TODO: Add previous assignment code (if needed)
		if (Frame::cosTheta(i.wi) >= 0.f && Frame::cosTheta(i.wo) >= 0.f)
		{
			val += diffuseReflectance->eval(worldData, i) * INV_PI;

			float exp = exponent->eval(worldData, i);
			float cosTheta = glm::fclamp(glm::dot(i.wi, reflect(i.wo)), 0.f, 1.f);
			val += specularReflectance->eval(worldData, i) * (exp + 2) * INV_TWOPI * std::powf(cosTheta, exp);

			val *= scale;
			val *= Frame::cosTheta(i.wi);
		}
		return val;
	}

	float pdf(const SurfaceInteraction& i) const override {
		float pdf = 0.f;
		// TODO: Implement this
		v3f reflectionLocalDir = reflect(i.wo);
		Frame reflectionSpace(reflectionLocalDir);
		v3f brdfSample = reflectionSpace.toLocal(i.wi);
		float exp = exponent->eval(worldData, i);

		pdf = Warp::squareToPhongLobePdf(brdfSample, exp);
		return pdf;
	}

	v3f sample(SurfaceInteraction& i, const v2f& _sample, float* pdf) const override {
		v3f val(0.f);
		// TODO: Implement this
		v3f reflectionLocalDir = reflect(i.wo);
		Frame reflectionSpace(reflectionLocalDir);

		float exp = exponent->eval(worldData, i);
		v3f brdfSample = Warp::squareToPhongLobe(_sample, exp);
		*pdf = Warp::squareToPhongLobePdf(brdfSample, exp);

		v3f brdfSampleLocal = reflectionSpace.toWorld(brdfSample);
		i.wi = brdfSampleLocal;
		val = eval(i);

		return val;
	}

    std::string toString() const override { return "Phong"; }
};

TR_NAMESPACE_END