/*
	This file is part of TinyRender, an educative rendering system.

	Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
	Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

/**
* Jack's addendum to solve a mystery.
*/

bool quadratic(double a, double b, double c, double& t0, double& t1) {
	double discriminant = b * b - 4 * a*c;

	if (discriminant > 0)
	{
		double sqrtDiscriminant = sqrt(discriminant);
		double inv2a = 1 / (2 * a);
		t0 = (-b + sqrtDiscriminant) * inv2a;
		t1 = (-b - sqrtDiscriminant) * inv2a;
		return true;
	}
	else if (discriminant == 0)
	{
		t0 = (-b + sqrt(discriminant)) / (2 * a);
		t1 = t0;
		return true;
	}

	return false;
}

bool raySphereIntersect(const Ray &ray_, v3f sphereCenter, float sphereRadius, float& outT) {
	v3f newOrigin(ray_.o - sphereCenter);

	double c = glm::dot(newOrigin, newOrigin) - (sphereRadius * sphereRadius);
	double b = glm::dot(newOrigin, ray_.d) * 2.0;
	double a = glm::dot(ray_.d, ray_.d);

	double r1;
	double r2;

	r1 = ray_.min_t;

	bool intersects = quadratic(a, b, c, r1, r2);

	double closestT;
	r1 = r1 < 0 ? std::numeric_limits<double>::infinity() : r1;
	r2 = r2 < 0 ? std::numeric_limits<double>::infinity() : r2;

	//If a hit point is out of ray bounds, then it shouldnt be considered in fmin, hence we set it to infinity.
	r1 = (r1 > ray_.min_t && r1 < ray_.max_t) ? r1 : std::numeric_limits<double>::infinity();
	r2 = (r2 > ray_.min_t && r2 < ray_.max_t) ? r2 : std::numeric_limits<double>::infinity();
	if (intersects)
	{
		closestT = fmin(r1, r2);
		if (closestT > ray_.min_t && closestT < ray_.max_t)
		{
			outT = closestT;
			return true;
		}
	}

	return false;
}

/**
 * Direct illumination integrator with MIS
 */
struct DirectIntegrator : Integrator {
	explicit DirectIntegrator(const Scene& scene) : Integrator(scene) {
		m_emitterSamples = scene.config.integratorSettings.di.emitterSamples;
		m_bsdfSamples = scene.config.integratorSettings.di.bsdfSamples;
		m_samplingStrategy = scene.config.integratorSettings.di.samplingStrategy;
	}

	static inline float balanceHeuristic(float nf, float fPdf, float ng, float gPdf) {
		float f = nf * fPdf, g = ng * gPdf;
		return f / (f + g);
	}

	void sampleSphereByCosineHemisphere(const p2f& sample,
		const v3f& n,
		const p3f& pShading,
		const v3f& emitterCenter,
		float emitterRadius,
		v3f& wiW,
		float& pdf) const {
		// TODO: Implement this
	}

	void sampleSphereByArea(const p2f& sample,
		const p3f& pShading,
		const v3f& emitterCenter,
		float emitterRadius,
		v3f& pos,
		v3f& ne,
		v3f& wiW,
		float& pdf) const {
		// TODO: Implement this
		ne = Warp::squareToUniformSphere(sample);
		pos = ne * emitterRadius + emitterCenter;
		wiW = normalize(pos - pShading);
		pdf = 1.f / (4 * M_PI * emitterRadius * emitterRadius);
	}

	void sampleSphereBySolidAngle(const p2f& sample,
		const p3f& pShading,
		const v3f& emitterCenter,
		float emitterRadius,
		v3f& wiW,
		float& pdf) const {
		// TODO: Implement this

		v3f centerDir = glm::normalize(emitterCenter - pShading);
		Frame subtendedSolidAngleFrame(centerDir);

		//From the points perspective, the sphere looks like a circle, and from the center of the circle to the edge, the subtended angle is thetaMax
		float sinThetaMaxSquared = emitterRadius * emitterRadius / glm::distance2(pShading, emitterCenter);
		float cosThetaMax = std::sqrtf(std::fmax(0.f, 1.f - sinThetaMaxSquared));

		//float thetaMax = glm::asin(emitterRadius / glm::distance(pShading, emitterCenter));
		//float cosThetaMax = glm::cos(thetaMax);

		//Sample cone uniformly based on cosThetaMax
		float cosTheta = (1.f - sample.x) + (sample.x * cosThetaMax);
		float phi = sample.y * M_PI * 2.0f;
		float sinTheta = std::sqrtf(std::fmax(1.f - (cosTheta * cosTheta), 0));
		v3f sampleDir(sinTheta * std::cosf(phi), sinTheta * std::sinf(phi), cosTheta);

		//Transform the local sample direction to one oriented around the subtended solid angle direction
		sampleDir = subtendedSolidAngleFrame.toWorld(sampleDir);
		wiW = sampleDir;

		//Use the pdf for uniform cone sampling (measure is solid angle)
		pdf = INV_TWOPI * (1.f / (1.f - cosThetaMax));
	}

	v3f renderArea(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO: Implement this

		SurfaceInteraction sInteraction;
		bool intersects = scene.bvh->intersect(ray, sInteraction);
		if (intersects)
		{
			v3f emittedRadiance = getEmission(sInteraction);
			if (emittedRadiance != v3f(0.0f)) {
				Lr = emittedRadiance;
			}
			else {
				for (unsigned int i = 0; i < m_emitterSamples; ++i)
				{
					float emitterPdf;
					size_t id = selectEmitter(sampler.next(), emitterPdf);
					const Emitter& emitter = getEmitterByID(id);
					v3f emitterCenter = scene.getShapeCenter(emitter.shapeID);
					float emitterRadius = scene.getShapeRadius(emitter.shapeID);

					p2f sample = sampler.next2D();
					v3f outPos;
					v3f outNe;
					v3f outwiW;
					float outPdf;

					sampleSphereByArea(sample, sInteraction.p, emitterCenter, emitterRadius, outPos, outNe, outwiW, outPdf);
					float lightDistanceSquared = glm::distance2(outPos, sInteraction.p);
					float cosOutgoing = glm::dot(-outwiW, outNe);
					v3f wiLocal = sInteraction.frameNs.toLocal(outwiW);
					if (cosOutgoing <= 0.f || Frame::cosTheta(wiLocal) <= 0.f) {
						continue;
					}

					Ray directIlluminationRay = Ray(sInteraction.p, outwiW, Epsilon, sqrt(lightDistanceSquared) - Epsilon);
					SurfaceInteraction visibilityInteraction;
					bool intersects2 = scene.bvh->intersect(directIlluminationRay, visibilityInteraction);

					if (!intersects2)
					{
						v3f emittedRadiance = emitter.radiance;
						sInteraction.wi = wiLocal;
						float areaToSolidAngle = cosOutgoing * (1.f / lightDistanceSquared); //this takes d_A to d_O
						//I could alternative turn emitterPdf from a density over area to 1 over solid angle via lightDistanceSquared/cosOutgoing, the 1.f / pdf will make it equal to the jacobian transformation above

						Lr += emittedRadiance * getBSDF(sInteraction)->eval(sInteraction) * (1.f / outPdf) * (1.f / emitterPdf) * areaToSolidAngle;
					}
				}
				Lr /= m_emitterSamples;
			}
		}

		return Lr;
	}

	//Do I need to consider the cosine of the sphere light normal and wi??
	v3f renderCosineHemisphere(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO: Implement this

		SurfaceInteraction sInteraction;
		bool intersects = scene.bvh->intersect(ray, sInteraction);
		if (intersects)
		{
			v3f emittedRadiance = getEmission(sInteraction);
			if (emittedRadiance != v3f(0.0f)) {
				Lr = emittedRadiance;
			}
			else {
				for (unsigned int i = 0; i < m_emitterSamples; ++i)
				{
					p2f sample = sampler.next2D();
					v3f localRayDir = Warp::squareToCosineHemisphere(sample);
					v3f worldRayDir = normalize(sInteraction.frameNs.toWorld(localRayDir));
					Ray directIlluminationRay = Ray(sInteraction.p, worldRayDir, Epsilon, std::numeric_limits<float>::infinity());

					SurfaceInteraction visibilityInteraction;
					bool intersects2 = scene.bvh->intersect(directIlluminationRay, visibilityInteraction);
					if (intersects2)
					{
						v3f emittedRadiance = getEmission(visibilityInteraction);
						sInteraction.wi = localRayDir;
						Lr += emittedRadiance * getBSDF(sInteraction)->eval(sInteraction) * (1.0f / Warp::squareToCosineHemispherePdf(localRayDir));
					}
				}
				Lr /= m_emitterSamples;
			}
		}

		return Lr;
	}

	v3f renderBSDF(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO: Implement this

		SurfaceInteraction sInteraction;
		bool intersects = scene.bvh->intersect(ray, sInteraction);
		if (intersects)
		{
			v3f emittedRadiance = getEmission(sInteraction);
			if (emittedRadiance != v3f(0.0f)) {
				Lr = emittedRadiance;
			}
			else {
				for (unsigned int i = 0; i < m_bsdfSamples; ++i)
				{
					p2f sample = sampler.next2D();
					float pdf;
					v3f brdfCosTheta = getBSDF(sInteraction)->sample(sInteraction, sample, &pdf);

					SurfaceInteraction visibilityInteraction;
					bool intersects2 = scene.bvh->intersect(Ray(sInteraction.p, sInteraction.frameNs.toWorld(sInteraction.wi), Epsilon, std::numeric_limits<float>::infinity()), visibilityInteraction);
					if (intersects2)
					{
						v3f emittedRadiance = getEmission(visibilityInteraction);
						Lr += emittedRadiance * brdfCosTheta * (1.0 / pdf);
					}
				}
				Lr /= m_bsdfSamples;
			}
		}

		return Lr;
	}

	//My shape and intensity is generally correct, but I have so much less noise, and I dont understand why??
	v3f renderSolidAngle(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO: Implement this

		SurfaceInteraction sInteraction;
		bool intersects = scene.bvh->intersect(ray, sInteraction);
		if (intersects)
		{
			v3f emittedRadiance = getEmission(sInteraction);
			if (emittedRadiance != v3f(0.0f)) {
				Lr = emittedRadiance;
			}
			else {
				for (unsigned int i = 0; i < m_emitterSamples; ++i)
				{
					float emitterPdf;
					size_t id = selectEmitter(sampler.next(), emitterPdf);
					const Emitter& emitter = getEmitterByID(id);
					v3f emitterCenter = scene.getShapeCenter(emitter.shapeID);
					float emitterRadius = scene.getShapeRadius(emitter.shapeID);

					p2f sample = sampler.next2D();
					v3f outwiW;
					float outPdf;

					sampleSphereBySolidAngle(sample, sInteraction.p, emitterCenter, emitterRadius, outwiW, outPdf);
					v3f wiLocal = sInteraction.frameNs.toLocal(outwiW);

					if (Frame::cosTheta(wiLocal) <= 0.f) {
						//Lr += v3f(1.f); //debugging
						continue;
					}

					Ray directIlluminationRay = Ray(sInteraction.p, outwiW, Epsilon, glm::distance(sInteraction.p, emitterCenter) + Epsilon);

					//SurfaceInteraction visibilityInteraction = SurfaceInteraction();
					SurfaceInteraction visibilityInteraction;
					bool intersects2 = scene.bvh->intersect(directIlluminationRay, visibilityInteraction);
					if (intersects2) {
						if (emitter.shapeID == visibilityInteraction.shapeID)
						{
							v3f emittedRadiance = emitter.radiance;
							sInteraction.wi = wiLocal;
							Lr += emittedRadiance * getBSDF(sInteraction)->eval(sInteraction) * (1.f / outPdf) * (1.f / emitterPdf);

							//debugging
							if (Frame::cosTheta(sInteraction.wo) < 0.f) {
								//Lr = v3f(0.f, 16.f, 0.f);
							}
						}
					}
					else {
						float outT;
						bool test = raySphereIntersect(directIlluminationRay, emitterCenter, emitterRadius, outT);
						if (test) {
							v3f emittedRadiance = emitter.radiance;
							sInteraction.wi = wiLocal;
							Lr += emittedRadiance * getBSDF(sInteraction)->eval(sInteraction) * (1.f / outPdf) * (1.f / emitterPdf);
						}
					}
				}
				Lr /= m_emitterSamples;

			}
		}

		return Lr;
	}

	v3f renderMIS(const Ray& ray, Sampler& sampler) const {
		v3f Lr(0.f);
		// TODO: Implement this

		v3f emitterEstimator(0.f);
		v3f bsdfEstimator(0.f);

		SurfaceInteraction sInteraction;
		bool intersects = scene.bvh->intersect(ray, sInteraction);
		if (intersects)
		{
			v3f emittedRadiance = getEmission(sInteraction);
			if (emittedRadiance != v3f(0.0f)) {
				Lr = emittedRadiance;
			}
			else {
				//Build estimator for emitter distribution
				for (unsigned int i = 0; i < m_emitterSamples; ++i)
				{
					float emitterPdf;
					size_t id = selectEmitter(sampler.next(), emitterPdf);
					const Emitter& emitter = getEmitterByID(id);
					v3f emitterCenter = scene.getShapeCenter(emitter.shapeID);
					float emitterRadius = scene.getShapeRadius(emitter.shapeID);

					p2f sample = sampler.next2D();
					v3f outwiW;
					float outPdf;

					sampleSphereBySolidAngle(sample, sInteraction.p, emitterCenter, emitterRadius, outwiW, outPdf);
					v3f wiLocal = sInteraction.frameNs.toLocal(outwiW);

					if (Frame::cosTheta(wiLocal) <= 0.f) {
						continue;
					}

					//Ray directIlluminationRay = Ray(sInteraction.p, outwiW, Epsilon, glm::distance(sInteraction.p, emitterCenter) - Epsilon);
					Ray directIlluminationRay = Ray(sInteraction.p, outwiW, Epsilon, std::numeric_limits<float>::infinity());
					SurfaceInteraction visibilityInteraction;
					bool intersects2 = scene.bvh->intersect(directIlluminationRay, visibilityInteraction);
					if (intersects2) {
						if (emitter.shapeID == visibilityInteraction.shapeID)
						{
							v3f emittedRadiance = emitter.radiance;
							sInteraction.wi = wiLocal;
							float bsdfPdf = getBSDF(sInteraction)->pdf(sInteraction);
							float weight = balanceHeuristic(m_emitterSamples, outPdf * emitterPdf, m_bsdfSamples, bsdfPdf);

							emitterEstimator += emittedRadiance * getBSDF(sInteraction)->eval(sInteraction) * weight * (1.f / outPdf) * (1.f / emitterPdf);
						}
					}
					else {
						float outT;
						bool test = raySphereIntersect(directIlluminationRay, emitterCenter, emitterRadius, outT);
						if (test) {
							v3f emittedRadiance = emitter.radiance;
							sInteraction.wi = wiLocal;
							float bsdfPdf = getBSDF(sInteraction)->pdf(sInteraction);
							float weight = balanceHeuristic(m_emitterSamples, outPdf * emitterPdf, m_bsdfSamples, bsdfPdf);

							emitterEstimator += emittedRadiance * getBSDF(sInteraction)->eval(sInteraction) * weight * (1.f / outPdf) * (1.f / emitterPdf);
						}
					}
				}
				emitterEstimator = m_emitterSamples == 0 ? v3f(0.f) : (emitterEstimator / m_emitterSamples);

				//Build estimator for bsdf distribution
				for (unsigned int i = 0; i < m_bsdfSamples; ++i)
				{
					p2f sample = sampler.next2D();
					float outPdf;
					v3f brdfCosTheta = getBSDF(sInteraction)->sample(sInteraction, sample, &outPdf);

					SurfaceInteraction visibilityInteraction;
					bool intersects2 = scene.bvh->intersect(Ray(sInteraction.p, sInteraction.frameNs.toWorld(sInteraction.wi), Epsilon, std::numeric_limits<float>::infinity()), visibilityInteraction);
					if (intersects2)
					{
						v3f emittedRadiance = getEmission(visibilityInteraction);
						if (emittedRadiance == v3f(0.f)) {
							continue; //we did not intersect an emitter, so there is nothing more to do
						}

						const Emitter& emitter = getEmitterByID(getEmitterIDByShapeID(visibilityInteraction.shapeID));
						v3f emitterCenter = scene.getShapeCenter(emitter.shapeID);
						float emitterRadius = scene.getShapeRadius(emitter.shapeID);

						//From the points perspective, the sphere looks like a circle, and from the center of the circle to the edge, the subtended angle is thetaMax
						float sinThetaMaxSquared = emitterRadius * emitterRadius / glm::distance2(sInteraction.p, emitterCenter);
						float cosThetaMax = std::sqrtf(std::fmax(0.f, 1.f - sinThetaMaxSquared));

						//Use the pdf for uniform cone sampling (measure is solid angle)
						float emitterSolidAnglePdf = INV_TWOPI * (1.f / (1.f - cosThetaMax));
						emitterSolidAnglePdf *= getEmitterPdf(emitter);

						float weight = balanceHeuristic(m_bsdfSamples, outPdf, m_emitterSamples, emitterSolidAnglePdf);

						bsdfEstimator += emittedRadiance * brdfCosTheta * weight * (1.f / outPdf);
					}
				}
				bsdfEstimator = m_bsdfSamples == 0 ? v3f(0.f) : (bsdfEstimator / m_bsdfSamples);

				//Combine the 2 estimators
				Lr += emitterEstimator + bsdfEstimator;
			}
		}

		return Lr;
	}

	v3f render(const Ray& ray, Sampler& sampler) const override {
		if (m_samplingStrategy == "mis")
			return this->renderMIS(ray, sampler);
		else if (m_samplingStrategy == "area")
			return this->renderArea(ray, sampler);
		else if (m_samplingStrategy == "solidAngle")
			return this->renderSolidAngle(ray, sampler);
		else if (m_samplingStrategy == "cosineHemisphere")
			return this->renderCosineHemisphere(ray, sampler);
		else if (m_samplingStrategy == "bsdf")
			return this->renderBSDF(ray, sampler);
		std::cout << "Error: wrong strategy" << std::endl;
		exit(EXIT_FAILURE);
	}

	size_t m_emitterSamples;     // Number of emitter samples
	size_t m_bsdfSamples;        // Number of BSDF samples
	string m_samplingStrategy;   // Sampling strategy to use
};

TR_NAMESPACE_END