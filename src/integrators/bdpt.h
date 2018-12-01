/*
	This file is part of TinyRender, an educative rendering system.

	All code in this file was written by Jack Minnetian.
*/

#pragma once

#include <core/platform.h>
#include <core/integrator.h>

extern std::unique_ptr<std::mutex[]> g_FrameBufferLocks;

TR_NAMESPACE_BEGIN

#define LIGHT_TRACING 0
#define PATH_TRACING 0
#define NO_RR 1


/**
 * Bidirectional path tracer integrator
 */
	struct PathVertex {
	//Intersection which generated the current vertex. 
	SurfaceInteraction hit;
	//This is the radiance or importance carried by this vertex in the path for light paths and eye paths respectively. 
	v3f throughput;
	float vcm;
	float vc;
	//The probability that the subpath is continued past the current vertex.
	float rr;

	PathVertex(SurfaceInteraction hit, v3f throughput, float vcm, float vc, float rr) : hit(hit), throughput(throughput), vcm(vcm), vc(vc), rr(rr) {}
};

struct BDPTIntegrator : Integrator {
	explicit BDPTIntegrator(const Scene& scene) : Integrator(scene) {
		m_rrDepth = scene.config.integratorSettings.pt.rrDepth;
		m_rrProb = scene.config.integratorSettings.pt.rrProb;

		std::cout << "Samples per pixel are: " << scene.config.spp << std::endl;
	}

	//Returns true if the eye subpath hit a light source, in which case no light subpath can or should be generated (light source aka emitters dont have bsdfs in this model)
	inline void eyeSubpathWalk(const std::vector<PathVertex>& lightSubpath, const Ray& ray, SurfaceInteraction hit, Sampler& sampler, std::vector<PathVertex>& eyeSubpath, v3f& Li) const {
		bool isPathPureSpecular = true;

		v3f cameraForward = glm::normalize(scene.config.camera.at - scene.config.camera.o);
		float cosCamera = glm::dot(cameraForward, ray.d);
		//This is the distance of the near plane at which the pdf for the area of a pixel is 1 
		float virtualNearPlaneDistance = (1.f / std::tanf(deg2rad * scene.config.camera.fov * 0.5f)) * rgb->height * 0.5f;
		float imagePointToCameraDist = virtualNearPlaneDistance / cosCamera;
		float imageAreaToCameraSolidAngle = imagePointToCameraDist * imagePointToCameraDist / cosCamera;

		//Pdf for sampling the given point in the pixel is 1, then we convert the measure to solid angle as the tech paper states (Eq. 31)
		float t1Pdf = 1.f * imageAreaToCameraSolidAngle;

		Ray wi(ray);
		v3f throughput(1.f);
		float vc = 0.f;
		float vcm = rgb->width * rgb->height * (1.f / t1Pdf);

		int depth = 1;
		float rrProbability = 1.f;
		//We can have infinite recursion atm by having 2 mirrors point to each other with perfect reflection
		//We likely want to code in a maximum bounce count regardless of bias concerns
		while (depth < m_rrDepth || (sampler.next() < rrProbability && !NO_RR)) {
			SurfaceInteraction hit = SurfaceInteraction();
			if (!scene.bvh->intersect(wi, hit))
				break;

			const float distSquared = hit.t * hit.t;
			const float absCosIn = std::fabs(hit.wo.z);

			vcm *= (distSquared / absCosIn); //forward geometric term for solid angle to area jacobian
			vc *= (1.f / absCosIn);          //same as above but the reverse geometric term in vc also has distSquared, so they cancel out and cosOut/cosIn is left

			if (getEmission(hit) != v3f(0.f))
			{
				const Emitter& emitter = getEmitterByID(getEmitterIDByShapeID(hit.shapeID));
				float emitterPdf = getEmitterPdf(emitter);
				if (depth > 1)
				{
					v3f contribution = emitter.radiance * throughput;

					float emitterPositionPdf_a = 1.f / (emitter.area * emitterPdf); //this is p_connect and p_trace
					float emitterDirectionPdf_w = Warp::squareToUniformHemispherePdf(hit.wo);

					//p_trace has a measure of emitter area (t-1 vertex is on the emitter itself, and p_trace is w.r.t area)
					float cameraWeight = emitterPositionPdf_a * vcm + (emitterPositionPdf_a * emitterDirectionPdf_w) * vc;
					float misWeight = 1.f / (1.f + cameraWeight);


#if !PATH_TRACING && !LIGHT_TRACING
					//We only use MIS weights if the path had non specular interactions (otherwise the path was deterministic)
					if (!isPathPureSpecular) {
						contribution *= misWeight;
					}
#endif
					//The following code is adding a 2nd sample of direct illumination (where depth = 2)
					//because at depth = 1, we did an explicit connection to the light and added it to Li
					//and now this provides another connection to the light and also adds it to Li, so we 
					//are effectively getting 2 samples of direct illumination into 1 sample. However, if 
					//the previous vertex was specular, then NEE doesnt work, and we need this contribution.
					//This issue only arises in path tracing. Otherwise the MIS weight properly combines
					//the 2 samples of direct illumination. Also we only add this contribution for pure
					//specular paths because checking the last vertex alone causes fireflies due to sampling
					//a hotspot with a lot probability (also only an issue with path tracing).
#if PATH_TRACING
					if (isPathPureSpecular) {
						Li += contribution;
					}
#else
					//When using BDPT, we can always add this contribution due to MISing the multiple
					//direct illumination samples. We also need this for transmission.
					Li += contribution;
#endif
				}
				else if (depth == 1) //This corresponds to an LE path
				{
					Li += getEmission(hit);
				}
				break;
			}

			//float rrProbability = (depth + 1) < m_rrDepth ? 1.f : m_rrProb;
			//float rrProbability = (depth + 1) < m_rrDepth ? 1.f : std::fmin(1.f, getLuminance(throughput));
			rrProbability = (depth + 1) < m_rrDepth ? 1.f : (getLuminance(throughput) < 0.01f ? 0.5f : 1.f);
#if NO_RR
			rrProbability = 1.f;
#endif
			//eyeSubpath.emplace_back(hit, throughput, vcm, vc, rrProbability);
			PathVertex eyeVertex(hit, throughput, vcm, vc, rrProbability);

			auto hitBSDFType = getBSDF(hit)->getType();
			if (!(hitBSDFType & BSDF::EDelta))
			{
				isPathPureSpecular = false;

				//Connect light (next event estimation)
				Li += connectToLight(eyeVertex, sampler);

				//Connect vertices
#if !LIGHT_TRACING && !PATH_TRACING
				for (int i = 0; i < lightSubpath.size(); i++) {
					Li += connectVertices(lightSubpath[i], eyeVertex);
				}
#endif
			}

			if (!ContinuePathRandomWalk(hit, sampler, rrProbability, eyeVertex, throughput, depth, vc, vcm, wi))
				break;
		}
	}

	//I need to give this more thought, I am not sure about the geometric term and pdf measures nor RR
	inline void lightSubpathWalk(Sampler& sampler, std::vector<PathVertex>& lightSubpath, const SurfaceInteraction& shadingPoint) const {
		//Generate a light sample to start the subpath
		float emitterPdf, emitterAreaPdf, emitterEmissionPdf;
		v3f normalOut, positionOut;
		size_t id = selectEmitter(sampler.next(), emitterPdf);
		const Emitter& emitter = getEmitterByID(id);
		sampleEmitterPosition(sampler, emitter, normalOut, positionOut, emitterAreaPdf);
		v3f emitterDirection = Warp::squareToUniformHemisphere(sampler.next2D()); //Condition pdf of the next vertex given the vertex on the light. Maybe cosine sampling is better?
		emitterEmissionPdf = Warp::squareToUniformHemispherePdf(emitterDirection) * emitterAreaPdf;
		emitterAreaPdf *= emitterPdf;
		emitterEmissionPdf *= emitterPdf;

		Frame lightFrame(normalOut);
		Ray wi = Ray(positionOut, lightFrame.toWorld(emitterDirection));

		v3f throughput = emitterDirection.z * emitter.getRadiance() * (1.f / emitterEmissionPdf);
		//Both P0 trace and P0 connect are emitterAreaPdf and P1 * P0trace = emitterEmissionPdf
		float vc = emitterDirection.z * (1.f / emitterEmissionPdf); //distance squared terms are cancelled out by geometric term in numerator and solid angle to area jacobian in denominator
		//vcm is just 1 / probability of going in the direction to v1 from the light, this is equal to 1 / emitterDirectionPdf, measure is solid angle at the emitter
		float vcm = emitterAreaPdf / emitterEmissionPdf;

		if (emitterDirection.z <= 0.f) {
			std::cout << "Emitter cosOut was lequal to 0, deal with this." << std::endl;
			return;
		}

		int depth = 1;
		//How does this while loop work for the first iteration when there is no lightSubpath.back()? what if m_rrDepth is 1? will this crash?
		//What happens if a light path hits a light again on one of its bounces? what do we do?
		float rrProbability = 1.f;
		while (depth < m_rrDepth || (sampler.next() < rrProbability && !NO_RR)) {
			SurfaceInteraction hit = SurfaceInteraction();
			if (!scene.bvh->intersect(wi, hit))
				break;

			const float distSquared = hit.t * hit.t;
			const float absCosIn = std::fabs(hit.wo.z);

			vcm *= (distSquared / absCosIn); //forward geometric term for solid angle to area jacobian
			vc *= (1.f / absCosIn);          //same as above but the reverse geometric term in vc also has distSquared, so they cancel out and cosOut/cosIn is left

			//float rrProbability = (depth+1) < m_rrDepth ? 1.f : m_rrProb;
			//float rrProbability = (depth + 1) < m_rrDepth ? 1.f : std::fmin(1.f, getLuminance(throughput));
			rrProbability = (depth + 1) < m_rrDepth ? 1.f : (getLuminance(throughput) < 0.01f ? 0.5f : 1.f);
#if NO_RR
			rrProbability = 1.f;
#endif
			auto hitBSDFType = getBSDF(hit)->getType();
			PathVertex lightVertex(hit, throughput, vcm, vc, rrProbability);

			if (!(hitBSDFType & BSDF::EDelta))
				connectToCamera(lightVertex);

			if (!ContinuePathRandomWalk(hit, sampler, rrProbability, lightVertex, throughput, depth, vc, vcm, wi))
				break;

			if (!(hitBSDFType & BSDF::EDelta))
				lightSubpath.push_back(lightVertex);
		}
	}

	v3f render(const Ray& ray, Sampler& sampler) const override {
		v3f Li(0.f);

		Ray r = ray;
		SurfaceInteraction hit;

		if (scene.bvh->intersect(r, hit)) {
			std::vector<PathVertex> eyeSubpath;
			std::vector<PathVertex> lightSubpath;

#if LIGHT_TRACING
			lightSubpathWalk(sampler, lightSubpath, hit);
			Li += getEmission(hit);
#elif PATH_TRACING
			eyeSubpathWalk(lightSubpath, ray, hit, sampler, eyeSubpath, Li);
#else 
			lightSubpathWalk(sampler, lightSubpath, hit);
			eyeSubpathWalk(lightSubpath, ray, hit, sampler, eyeSubpath, Li);
#endif
		}

		return Li;
	}

	bool ContinuePathRandomWalk(SurfaceInteraction& hit, Sampler& sampler, const float& rrProbability, PathVertex& pathVertex,
		v3f& throughput, int& depth, float& vc, float& vcm, Ray& wi) const
	{
		const BSDF* bsdf = getBSDF(hit);
		bool isDeltaBsdf = (bsdf->getType() & BSDF::EDelta);
		float bsdfPdf_w;
		v3f bsdfCosTheta = bsdf->sample(hit, sampler.next2D(), &bsdfPdf_w);
		bsdfPdf_w *= rrProbability;
		pathVertex.hit.wi = hit.wi;

		float absCosOut = std::fabs(hit.wi.z);
		if (bsdfCosTheta == v3f(0.f)) //Transmissive materials can have cosOut <= 0.f but positive bsdfCosTheta
			return false;

		//Be wary that from the light paths perspective, the cosTheta below is cosOut. This is correct, for the same reason that the throughput from the emitter has the cosOut term. 
		//cosOut is all that will remain throughout all the bounces after we do all the change of measures for the pdfs and factor in cosIn for the projected solid angle.
		//We have a surface point that emits light with a probability density of sampling that surface point whose measure is the area at that point,
		//and we choose an outgoing direction vector from a pdf with measure of solid angle at the emitting point, which leads us to another point in the scene who will receive this light.
		//This is exactly like the first case where we pick a random point on an emitter and a direction from it to a point in the scene with a pdf of solid angle measure.
		//We already have the area pdf of picking the shading point (its the pdf we have accumulated so far), and bsdfPdf is the probability of picking the direction.
		//Thus just like before, when all is said and done, only the outgoing cosine will remains of the terms. It is possible to see this by doing it by hand for the base case and the recursion is obvious.
		throughput *= bsdfCosTheta * (1.f / bsdfPdf_w); //geometric terms in the throughput and pdf solid angle area jacobian cancel out, thus we can ignore them

		depth++;

		//Given the new edge, find the probability of computing the old edge
		SurfaceInteraction reverseHit(hit);
		std::swap(reverseHit.wi, reverseHit.wo);
		//Delta bsdfs are deterministic and thus calling pdf on them will always return 0
		float previousVertexReversePdf_w = isDeltaBsdf ? bsdfPdf_w : bsdf->pdf(reverseHit) * rrProbability;

		if (isDeltaBsdf) {
			//vcm is 1/delta_distrib, so it is 0. For vc, there is a 1/delta_distrib 
			//on the outside of Eq.35 which means the vcm_prev term cant be evaluated
			//but previousVertexReversePdf_w is also a delta_distrib and they cancel out
			//leaving the geometric term (absCosOut) * vc_prev
			vc = (absCosOut / bsdfPdf_w) * (previousVertexReversePdf_w * vc); //Eq.54
			vcm = 0.f;														  //Eq.53	
		}
		else {
			vc = (absCosOut / bsdfPdf_w) * (vcm + previousVertexReversePdf_w * vc);
			vcm = 1.f / bsdfPdf_w;
		}

		//Generate the new ray for the next iteration
		wi = Ray(hit.p, hit.frameNs.toWorld(hit.wi));

		return true;
	}

	//This corresponds to Vertex Connection (t=1), Eqs 46 & 47 in ImplementingVCM_TechRep2012_rev2
	//I think the issue is going to lie in here if anywhere
	void connectToCamera(const PathVertex& lightVertex) const {
		v3f cameraForward = glm::normalize(scene.config.camera.at - scene.config.camera.o);
		v3f eyeToLightVert = lightVertex.hit.p - scene.config.camera.o;
		float invDistanceSquared = 1.f / glm::length2(eyeToLightVert);
		eyeToLightVert *= glm::sqrt(invDistanceSquared);

		int xPixel, yPixel;
		splatToImagePlane(lightVertex.hit.p, xPixel, yPixel);
		if (xPixel < 0 || yPixel < 0 || xPixel >= rgb->width || yPixel >= rgb->height) {
			return;
		}

		float cosCamera = glm::dot(cameraForward, eyeToLightVert);
		if (cosCamera <= 0.f)
			return;

		SurfaceInteraction bsdfEvalInteraction(lightVertex.hit);
		bsdfEvalInteraction.wi = bsdfEvalInteraction.frameNs.toLocal(-eyeToLightVert);
		v3f bsdfCosTheta = getBSDF(bsdfEvalInteraction)->eval(bsdfEvalInteraction);
		if (bsdfCosTheta == v3f(0.f) || bsdfEvalInteraction.wi.z <= 0.f)
			return;

		//visibility test
		if (this->visibilityQuery(scene.config.camera.o, lightVertex.hit.p))
			return;

		//This is the distance of the near plane at which the pdf for the area of a pixel is 1 
		float virtualNearPlaneDistance = (1.f / std::tanf(deg2rad * scene.config.camera.fov * 0.5f)) * rgb->height * 0.5f;
		float imagePointToCameraDist = virtualNearPlaneDistance / cosCamera;
		//Our MC estimator is for the measurement equation defined w.r.t pixel area, thus we need to change our pdf measure from scene surface area to image plane pixel area
		float imageAreaToCameraSolidAngle = imagePointToCameraDist * imagePointToCameraDist / cosCamera;
		float cameraSolidAngleToSurfaceArea = bsdfEvalInteraction.wi.z * invDistanceSquared;
		float imageAreaToSurfaceArea = imageAreaToCameraSolidAngle * cameraSolidAngleToSurfaceArea;
		float surfaceAreaToImageArea = 1.f / imageAreaToSurfaceArea;

		int nlightSubpath = rgb->width * rgb->height;
		//The pdf of the following estimate is 1. Our MC estimator is of the measurement equation, and so we need to divide by the pdf for the pixel area, but that is 1 so we dont write anything.
		v3f radiance = lightVertex.throughput * (bsdfCosTheta * (1.f / bsdfEvalInteraction.wi.z)) * (1.f / surfaceAreaToImageArea) * (1.f / nlightSubpath) * (1.f / scene.config.spp);

		//We now wish to compute the MIS weight, as given Eq 46 & 47
		float reversePdf_a = 1.f * imageAreaToSurfaceArea; //The pdf of picking a point in the pixel area is 1, and we convert the measure to scene surface area
		//reversePdf_a = imageAreaToCameraSolidAngle * bsdfEvalInteraction.wi.z * invDistanceSquared;

		SurfaceInteraction reverseHit(bsdfEvalInteraction);
		//reverseHit.wi = bsdfEvalInteraction.wo;
		//reverseHit.wo = bsdfEvalInteraction.wi;
		std::swap(reverseHit.wi, reverseHit.wo);
		float previousVertexReversePdf_w = getBSDF(reverseHit)->pdf(reverseHit) * lightVertex.rr;

		//We only consider spp when dividing the total contribution provided by the spp number of iterations (in the radiance computation)
		//The MIS weight is a weight used for a single estimate of the pixel. For a single estimate of a pixel the number of lightpaths is just the number of pixels. 
		//The total number of sampled subpaths using sampling techniques with t=1 like this one is nlightSubpath (not 1, like in the other MIS weights).
		//The majority of those samples contribute 0 (they only contribute to the pixel they intersect), but that doesnt change the fact that we took a sample from them
		//Ultimately, when sampling across all points in the scene, a vast majority will have a 0 contribution. To make up for this, while the probability of sampling one that does 
		//contribute is low, this means that in the MC estimate, when dividing by the pdf, we get a large contribution. It is thus important to note that every single light
		//path is a sample used in the MC estimator of the measurement equation for every pixel, it just happens to only contribute to 1 pixel hence the code.
		float lightWeight = (reversePdf_a / nlightSubpath) * (lightVertex.vcm + previousVertexReversePdf_w * lightVertex.vc); //NOTE: should this be + or multiply?? 
		float eyeWeight = 0.f; //eye weight is explicitly here for clarity sake, the compiler likely optimizes this out anyways
		float misWeight = 1.f / (lightWeight + 1.f + eyeWeight);

#if !LIGHT_TRACING && !PATH_TRACING
		radiance *= misWeight;
#endif

		//for now we will test light tracing only, so no MIS weights will be used
		int pixel = yPixel * rgb->width + xPixel;
		//rgb->data[pixel] += radiance; //this is subject to race conditions, I will single thread for now but multi-thread requires me to have a lock on this

		while (true) {
			std::mutex& lock = g_FrameBufferLocks[pixel];
			if (lock.try_lock()) {
				rgb->data[pixel] += radiance;
				lock.unlock();
				break;
			}
		}
	}

	//Implementing Eq.44 and Eq.45
	v3f connectToLight(const PathVertex& eyeVertex, Sampler& sampler) const {
		float emitterPdf;
		size_t emitterID = selectEmitter(sampler.next(), emitterPdf);
		const Emitter& emitter = getEmitterByID(emitterID);

		v3f emitterNormal, emitterPosition;
		float emitterPositionPdf;
		sampleEmitterPosition(sampler, emitter, emitterNormal, emitterPosition, emitterPositionPdf);
		v3f lightToEyeVertexDir = eyeVertex.hit.p - emitterPosition;
		float lightToEyeDistSqrd = glm::length2(lightToEyeVertexDir);
		lightToEyeVertexDir *= (1.f / glm::sqrt(lightToEyeDistSqrd));

		SurfaceInteraction bsdfEvalInteraction(eyeVertex.hit);
		bsdfEvalInteraction.wi = bsdfEvalInteraction.frameNs.toLocal(-lightToEyeVertexDir);

		float cosAtLight = glm::dot(emitterNormal, lightToEyeVertexDir);
		float cosAtEyeVertex = bsdfEvalInteraction.wi.z;

		if (cosAtLight <= 0.f || cosAtEyeVertex <= 0.f)
			return v3f(0.f);

		float emitterConnectPdf_a = emitterPdf * emitterPositionPdf;
		float emitterConnectPdf_w = emitterConnectPdf_a * lightToEyeDistSqrd / cosAtLight;
		float emitterDirectionPdf_w = Warp::squareToUniformHemispherePdf(lightToEyeVertexDir);

		v3f Li = getBSDF(bsdfEvalInteraction)->eval(bsdfEvalInteraction) * (1.f / emitterConnectPdf_w) * eyeVertex.throughput * emitter.radiance;

		if (Li == v3f(0.f))
			return v3f(0.f);

		//visibility test
		if (this->visibilityQuery(eyeVertex.hit.p, emitterPosition))
			return v3f(0.f);

		//We want the probability of going from the eyeVertex to the light if we were not doing next event estimation
		//Is wi and wo in the correct order for this surface interaction?
		float lightPathReversePdf_w = getBSDF(bsdfEvalInteraction)->pdf(bsdfEvalInteraction) * eyeVertex.rr; //the probability that the eyePath generated the lightVertex
		float lightWeight = lightPathReversePdf_w / emitterConnectPdf_w; //Eq.44 - they are both in solid angle measure, no need to change them both to area measure

		SurfaceInteraction bsdfReverseInteraction(bsdfEvalInteraction);
		//bsdfReverseInteraction.wo = bsdfEvalInteraction.wi;
		//bsdfReverseInteraction.wi = bsdfEvalInteraction.wo;
		std::swap(bsdfReverseInteraction.wi, bsdfReverseInteraction.wo);
		float eyePathPreviousVertexReversePdf_w = getBSDF(bsdfReverseInteraction)->pdf(bsdfReverseInteraction) * eyeVertex.rr; //the probability that the previous eyeVertex would be sampled if the current eyeVertex was in the lightpath

		//Now we want the probability of picking the current eyeVertex from the light subpath given the point on the light
		float eyePathCurrentVertexReversePdf_a = cosAtEyeVertex * (1.f / lightToEyeDistSqrd) * emitterDirectionPdf_w;

		float eyeWeight = eyePathCurrentVertexReversePdf_a * (eyeVertex.vcm + eyePathPreviousVertexReversePdf_w * eyeVertex.vc); //Eq.45

		float misWeight = 1.f / (lightWeight + 1.f + eyeWeight);

#if !PATH_TRACING && !LIGHT_TRACING
		Li *= misWeight;
#endif
		return Li;
	}

	//Implementing Eq.40 and Eq.41
	//The issue at the moment is the russian roulette
	v3f connectVertices(const PathVertex& lightVertex, const PathVertex& eyeVertex) const {
		v3f lightToEyeDir = eyeVertex.hit.p - lightVertex.hit.p;
		float invDistanceSquared = 1.f / glm::length2(lightToEyeDir);
		lightToEyeDir *= glm::sqrt(invDistanceSquared);

		SurfaceInteraction lightVertexHit(lightVertex.hit);
		lightVertexHit.wi = lightVertexHit.frameNs.toLocal(lightToEyeDir);
		SurfaceInteraction eyeVertexHit(eyeVertex.hit);
		eyeVertexHit.wi = eyeVertexHit.frameNs.toLocal(-lightToEyeDir);

		float cosAtLightVertex = lightVertexHit.wi.z;
		float cosAtEyeVertex = eyeVertexHit.wi.z;

		if (cosAtLightVertex <= 0.f || cosAtEyeVertex <= 0.f)
			return v3f(0.f);

		//visibility test
		if (this->visibilityQuery(eyeVertex.hit.p, lightVertex.hit.p))
			return v3f(0.f);

		//There are no pdfs to be considered here, the connection between the 2 points was a deterministic connection
		v3f Li = getBSDF(lightVertexHit)->eval(lightVertexHit) * getBSDF(eyeVertexHit)->eval(eyeVertexHit);
		Li *= lightVertex.throughput * eyeVertex.throughput * invDistanceSquared; //the 2 cosines of the geomtric terms are in the bsdfeval

		//Now we compute the MIS weight

		//Given the lightPath, what is the probability of generating the connection edge
		float eyePathReversePdf_w = getBSDF(lightVertexHit)->pdf(lightVertexHit) * lightVertex.rr;
		//lightVertexHit.wo = lightVertexHit.wi;
		//lightVertexHit.wi = lightVertex.hit.wo;
		std::swap(lightVertexHit.wi, lightVertexHit.wo);
		//Given the connection edge, what is the probability of generating the previous vertex in the light path
		float lightPathPreviousVertexReversePdf_w = getBSDF(lightVertexHit)->pdf(lightVertexHit) * lightVertex.rr;

		//Same as above
		float lightPathReversePdf_w = getBSDF(eyeVertexHit)->pdf(eyeVertexHit) * eyeVertex.rr;
		eyeVertexHit.wo = eyeVertexHit.wi;
		eyeVertexHit.wi = eyeVertex.hit.wo;
		float eyePathPreviousVertexReversePdf_w = getBSDF(eyeVertexHit)->pdf(eyeVertexHit) * eyeVertex.rr;

		float lightPathReversePdf_a = lightPathReversePdf_w * (cosAtLightVertex)* invDistanceSquared;
		float eyePathReversePdf_a = eyePathReversePdf_w * (cosAtEyeVertex)* invDistanceSquared;

		float lightWeight = lightPathReversePdf_a * (lightVertex.vcm + lightPathPreviousVertexReversePdf_w * lightVertex.vc);
		float eyeWeight = eyePathReversePdf_a * (eyeVertex.vcm + eyePathPreviousVertexReversePdf_w * eyeVertex.vc);

		float misWeight = 1.f / (lightWeight + 1.f + eyeWeight);

		return Li * misWeight;
	}

	void splatToImagePlane(const v3f position, int& xPixel, int& yPixel) const {
		float aspectRatio = (float)rgb->width / rgb->height;
		auto worldToCamera = glm::lookAt(scene.config.camera.o, scene.config.camera.at, scene.config.camera.up);
		auto cameraToClip = glm::perspective(deg2rad * scene.config.camera.fov, aspectRatio, 1.f, 1000.f);
		auto NDCToScreen = glm::scale(glm::mat4(1.f), v3f(rgb->width, rgb->height, 1.0f)) * glm::scale(glm::mat4(1.f), v3f(0.5f, -0.5f, 1.f)) * glm::translate(glm::mat4(1.f), v3f(1.f, -1.f, 0.f));
		v4f uv = worldToCamera * v4f(position, 1.f);
		uv = cameraToClip * uv;
		uv /= uv.w;
		uv = NDCToScreen * uv;
		xPixel = static_cast<int>(uv.x);
		yPixel = static_cast<int>(uv.y);
	}

	bool visibilityQuery(const v3f& start, const v3f& end, SurfaceInteraction* hit = nullptr) const {
		v3f dir = end - start;
		float dist = glm::length(dir);
		dir /= dist;
		IntersectionInfo intersectionInfo = IntersectionInfo();
		SurfaceInteraction sInteraction = SurfaceInteraction();
		Ray ray(start, dir, Epsilon, dist - 0.00001f);
		return scene.bvh->bvh->getIntersection(ray, &intersectionInfo, true);

		//USED TO BE://visibility test
		//SurfaceInteraction visibilityQuery = SurfaceInteraction();
		//Ray visibilityRay = Ray(eyeVertex.hit.p, -lightToEyeDir);
		//scene.bvh->intersect(visibilityRay, visibilityQuery);
		////We must also check the primitive ID in case of non-convex objects (concave objects)
		//if (visibilityQuery.shapeID != lightVertex.hit.shapeID || visibilityQuery.primID != lightVertex.hit.primID)
		//	return v3f(0.f);
	}

	int m_rrDepth;      // When to start Russian roulette
	float m_rrProb;     // Russian roulette probability
};

TR_NAMESPACE_END
