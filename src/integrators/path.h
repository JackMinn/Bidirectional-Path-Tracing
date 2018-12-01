/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
 * Path tracer integrator
 */
struct PathTracerIntegrator : Integrator {
    explicit PathTracerIntegrator(const Scene& scene) : Integrator(scene) {
        m_isExplicit = scene.config.integratorSettings.pt.isExplicit;
        m_maxDepth = scene.config.integratorSettings.pt.maxDepth;
        m_rrDepth = scene.config.integratorSettings.pt.rrDepth;
        m_rrProb = scene.config.integratorSettings.pt.rrProb;
		m_emitterSamples = scene.config.integratorSettings.pt.emitterSamples;
		m_bsdfSamples = scene.config.integratorSettings.pt.bsdfSamples;

		std::cout << "Samples per pixel are: " << scene.config.spp << std::endl;
		std::cout << "DI emitter samples are: " << m_emitterSamples << std::endl;
		std::cout << "DI bsdf samples are: " << m_bsdfSamples << std::endl;
		std::cout << "Max depth is: " << m_maxDepth << std::endl;
    }

	static inline float balanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
		float f = nf * fPdf, g = ng * gPdf;
		return f / (f + g);
	}

	v3f recursiveImplicit(Sampler& sampler, SurfaceInteraction& hit, int depth) const {
		if (depth < m_maxDepth) {
			p2f sample = sampler.next2D();
			float pdf;
			v3f brdfCosTheta = getBSDF(hit)->sample(hit, sample, &pdf);

			SurfaceInteraction visibilityInteraction;
			v3f wiW = hit.frameNs.toWorld(hit.wi);
			bool intersects = scene.bvh->intersect(Ray(hit.p, wiW, Epsilon, std::numeric_limits<float>::infinity()), visibilityInteraction);
			if (intersects)
			{
				if (getEmission(visibilityInteraction) == v3f(0.f)) {
					v3f Li = recursiveImplicit(sampler, visibilityInteraction, depth + 1);
					return Li * brdfCosTheta * (1.0 / pdf);
				}
				else {
					v3f Li = getEmission(visibilityInteraction);
					//This next line is to handle the situation with the ceiling square light, without it we assume uniform emission profile
					Li = glm::dot(visibilityInteraction.frameNs.n, -wiW) > 0.f ? Li : v3f(0.f); 
					return Li * brdfCosTheta * (1.0 / pdf);
				}
			}
			else {
				return v3f(0.f);
			}
		}
		else {
			return v3f(0.f);
		}
	}

	v3f recursiveExplicit(Sampler& sampler, SurfaceInteraction& hit, int depth) const {
		v3f Lr(0.f);
		float rrSample = sampler.next();
		v3f indirectEstimator(0.f);
		v3f directEstimator(0.f);

		//Either the depth or the russian roulette allows us to enter this branch
		if (depth < m_maxDepth || (m_maxDepth == -1 && (depth < m_rrDepth || rrSample < m_rrProb))) 
		{
			//Indirect illumination estimator
			//-Do we need to cumulate the RR probability for when we succeed the probability but hit a light source and so must try again? Or is it just a constant value?
			//-Also, if we rejection sample, does that not change the pdf? Say all but 1% of directions are permitted, can we still use the full pdf? Or do we need to track the failed samples?
			{
				unsigned int nSamples = 0;
				float indRRProb = 0.95f;
				float cumRRProb = 1.f;
				SurfaceInteraction visibilityInteraction;
				bool intersects = false;
				float pdf;
				v3f brdfCosTheta;
				do
				{
					//if (sampler.next() < indRRProb) {
						visibilityInteraction = SurfaceInteraction();
						p2f sample = sampler.next2D();
						brdfCosTheta = getBSDF(hit)->sample(hit, sample, &pdf);

						v3f wiW = hit.frameNs.toWorld(hit.wi);
						intersects = scene.bvh->intersect(Ray(hit.p, wiW, Epsilon, std::numeric_limits<float>::infinity()), visibilityInteraction);
						nSamples++;
						//cumRRProb *= indRRProb;
						//cumRRProb = indRRProb;
					/*}
					else {
						break;
					}*/
				} while (getEmission(visibilityInteraction) != v3f(0.f) && sampler.next() < indRRProb); //maybe put the 2nd russian roulette in here, so the first attempt is guaranteed, but rerolling a miss is rng
				cumRRProb = nSamples > 1 ? indRRProb : 1.f;

				//for debugging purposes
				//nSamples = 1;
				if (intersects && brdfCosTheta != v3f(0.f))
				{
					if (getEmission(visibilityInteraction) == v3f(0.f) && brdfCosTheta != v3f(0.f)) {
						v3f Li = recursiveExplicit(sampler, visibilityInteraction, depth + 1);
						indirectEstimator = Li * brdfCosTheta * (1.f / pdf) * (1.f / nSamples) * (1.f / cumRRProb);
					}
				}
			}

			//Direct illumination estimator
			{
				v3f emitterEstimator(0.f);
				v3f bsdfEstimator(0.f);

				for (unsigned int i = 0; i < m_emitterSamples; ++i)
				{
					float emitterPdf, emitterAreaPdf;
					v3f normalOut, positionOut;
					size_t id = selectEmitter(sampler.next(), emitterPdf);
					const Emitter& emitter = getEmitterByID(id);
					sampleEmitterPosition(sampler, emitter, normalOut, positionOut, emitterAreaPdf);

					v3f wiW = glm::normalize(positionOut - hit.p);
					v3f wiLocal = hit.frameNs.toLocal(wiW);
					float lightDistanceSquared = glm::distance2(positionOut, hit.p);
					float cosOutgoing = glm::dot(-wiW, normalOut);
					if (cosOutgoing > 0.f && Frame::cosTheta(wiLocal) > 0.f) {
						Ray directIlluminationRay = Ray(hit.p, wiW, Epsilon, std::numeric_limits<float>::infinity());
						SurfaceInteraction visibilityInteraction = SurfaceInteraction();
						bool intersects = scene.bvh->intersect(directIlluminationRay, visibilityInteraction);
						if (intersects)
						{
							if (visibilityInteraction.shapeID == emitter.shapeID) {
								v3f Li = getEmission(visibilityInteraction);
								hit.wi = wiLocal;
								float areaToSolidAngle = cosOutgoing * (1.f / lightDistanceSquared);

								float bsdfPdf = getBSDF(hit)->pdf(hit);
								//bring all measures into solid angle, probability of sampling a direction gets larger the further away the area light is as there a less possible directions
								//because the subtended solid angle gets smaller with distance, limit case is infinity (aka a delta function, which integrates to 1 over a single dw)
								float weight = balanceHeuristic(m_emitterSamples, emitterAreaPdf * emitterPdf * (1.f / areaToSolidAngle), m_bsdfSamples, bsdfPdf);

								emitterEstimator += weight * Li * getBSDF(hit)->eval(hit) * (1.f / emitterAreaPdf) * (1.f / emitterPdf) * areaToSolidAngle;
							}
						}
					}
				}
				emitterEstimator = m_emitterSamples == 0 ? v3f(0.f) : (emitterEstimator / m_emitterSamples);

				for (unsigned int i = 0; i < m_bsdfSamples; ++i)
				{
					p2f sample = sampler.next2D();
					float bsdfPdf;
					v3f brdfCosTheta = getBSDF(hit)->sample(hit, sample, &bsdfPdf);

					if (brdfCosTheta != v3f(0.f))
					{
						SurfaceInteraction visibilityInteraction = SurfaceInteraction();
						v3f wiW = hit.frameNs.toWorld(hit.wi);
						bool intersects = scene.bvh->intersect(Ray(hit.p, wiW, Epsilon, std::numeric_limits<float>::infinity()), visibilityInteraction);
						if (intersects)
						{
							v3f Li = getEmission(visibilityInteraction);
							if (Li == v3f(0.f)) {
								continue; //we did not intersect an emitter, so there is nothing more to do
							}

							const Emitter& emitter = getEmitterByID(getEmitterIDByShapeID(visibilityInteraction.shapeID));
							float emitterPdf = getEmitterPdf(emitter);
							float emitterAreaPdf = 1.f / emitter.area; //measure is area

							float lightDistanceSquared = glm::distance2(visibilityInteraction.p, hit.p);
							float cosOutgoing = glm::dot(-wiW, visibilityInteraction.frameNg.n);
							if (cosOutgoing > 0.f)
							{
								float areaToSolidAngle = cosOutgoing * (1.f / lightDistanceSquared); //we need this to take emitterAreaPdf to measure of solid angle

								//bring all measures into solid angle
								float weight = balanceHeuristic(m_bsdfSamples, bsdfPdf, m_emitterSamples, emitterPdf * emitterAreaPdf * (1.f / areaToSolidAngle));

								bsdfEstimator += weight * Li * brdfCosTheta * (1.f / bsdfPdf);
							}
						}
					}
				}
				bsdfEstimator = m_bsdfSamples == 0 ? v3f(0.f) : (bsdfEstimator / m_bsdfSamples);

				directEstimator = emitterEstimator + bsdfEstimator;
			}

			Lr = directEstimator + indirectEstimator;

			//If we are using russian roulette, and we are past m_rrDepth, then we got here by russian roullette, so scale the sample accordingly
			if (m_maxDepth == -1 && !(depth < m_rrDepth) && (rrSample < m_rrProb)) {
				Lr *= 1.f / m_rrProb;
			}
		}
		
		return Lr;
	}

    v3f renderImplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
        v3f Li(0.f);

        // TODO: Implement this

		//In this case Li is radiance incident on the camera plane, which has no brdf or surface normal to consider
		if (getEmission(hit) != v3f(0.f))
		{
			Li = getEmission(hit);
		}
		else {
			Li = recursiveImplicit(sampler, hit, 0);
		}

        return Li;
    }

    v3f renderExplicit(const Ray& ray, Sampler& sampler, SurfaceInteraction& hit) const {
        v3f Li(0.f);

        // TODO: Implement this

		//In this case Li is radiance incident on the camera plane, which has no brdf or surface normal to consider
		if (getEmission(hit) != v3f(0.f))
		{
			Li = getEmission(hit);
		}
		else {
			Li = recursiveExplicit(sampler, hit, 0);
		}

        return Li;
    }


    v3f render(const Ray& ray, Sampler& sampler) const override {
        Ray r = ray;
        SurfaceInteraction hit;

        if (scene.bvh->intersect(r, hit)) {
            if (m_isExplicit)
                return this->renderExplicit(ray, sampler, hit);
            else
                return this->renderImplicit(ray, sampler, hit);
        }
        return v3f(0.0);
    }

    int m_maxDepth;     // Maximum number of bounces
    int m_rrDepth;      // When to start Russian roulette
    float m_rrProb;     // Russian roulette probability
    bool m_isExplicit;  // Implicit or explicit
	size_t m_emitterSamples;     // Number of emitter samples for direct illumination
	size_t m_bsdfSamples;        // Number of BSDF samples for direct illumination
};

TR_NAMESPACE_END
