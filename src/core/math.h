/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

TR_NAMESPACE_BEGIN

inline float safeSqrt(float v){
    return std::sqrt(std::max(float(0), v));
}

/**
 * Computes barycentric coordinates.
 */
template<class T>
inline T barycentric(const T& a, const T& b, const T& c, const float u, const float v) {
    return a * (1 - u - v) + b * u + c * v;
}

/**
 * Restricts a value to a given interval.
 */
template<class T>
inline T clamp(T v, T min, T max) {
    return std::min(std::max(v, min), max);
}

/**
 * Checks if vector is zero.
 */
inline bool isZero(const v3f v) {
    return glm::dot(v, v) < Epsilon;
}

/**
 * Generates coordinate system.
 */
inline void coordinateSystem(const v3f& a, v3f& b, v3f& c) {
    if (std::abs(a.x) > std::abs(a.y)) {
        float invLen = 1.f / std::sqrt(a.x * a.x + a.z * a.z);
        c = v3f(a.z * invLen, 0.f, -a.x * invLen);
    } else {
        float invLen = 1.f / std::sqrt(a.y * a.y + a.z * a.z);
        c = v3f(0.f, a.z * invLen, -a.y * invLen);
    }
    b = glm::cross(c, a);
}

/**
 * Converts RGB value to luminance.
 */
inline float getLuminance(const v3f& rgb) {
    return glm::dot(rgb, v3f(0.212671f, 0.715160f, 0.072169f));
}

/**
 * Pseudo-random sampler (Mersenne Twister 19937) structure.
 */
struct Sampler {
    std::mt19937 g;
    std::uniform_real_distribution<float> d;
    explicit Sampler(int seed) {
        g = std::mt19937(seed);
        d = std::uniform_real_distribution<float>(0.f, 1.f);
    }
    float next() { return d(g); }
    p2f next2D() { return {d(g), d(g)}; }
    void setSeed(int seed) {
        g.seed(seed);
        d.reset();
    }
};

/**
 * 1D discrete distribution.
 */
struct Distribution1D {
    std::vector<float> cdf{0};
    bool isNormalized = false;

    inline void add(float pdfVal) {
        cdf.push_back(cdf.back() + pdfVal);
    }

    size_t size() {
        return cdf.size() - 1;
    }

    float normalize() {
        float sum = cdf.back();
        for (float& v : cdf) {
            v /= sum;
        }
        isNormalized = true;
        return sum;
    }

    inline float pdf(size_t i) const {
        assert(isNormalized);
        return cdf[i + 1] - cdf[i];
    }

    int sample(float sample) const {
        assert(isNormalized);
        const auto it = std::upper_bound(cdf.begin(), cdf.end(), sample);
        return clamp(int(distance(cdf.begin(), it)) - 1, 0, int(cdf.size()) - 2);
    }
};


/**
 * Warping functions.
 */
namespace Warp {
	inline v3f squareToUniformSphere(const p2f& sample) {
		v3f v(0.f);
		// TODO: Implement this
		float phi = sample.x * M_PI * 2.0f;
		float cosTheta = 1.f - (2.f * sample.y);
		float sinTheta = std::sqrtf(std::fmax(1.f - cosTheta * cosTheta, 0));
		v = v3f(sinTheta * std::cosf(phi), sinTheta * std::sinf(phi), cosTheta);
		return v;
	}

	inline float squareToUniformSpherePdf() {
		float pdf = 0.f;
		// TODO: Implement this
		pdf = INV_FOURPI;
		return pdf;
	}

	inline v3f squareToUniformHemisphere(const p2f& sample) {
		v3f v(0.f);
		// TODO: Implement this
		float phi = sample.x * M_PI * 2.0f;
		float cosTheta = sample.y;
		float sinTheta = std::sqrtf(std::fmax(1.f - (cosTheta * cosTheta), 0)); //fmax is just a precaution, cosTheta*cosTheta should never be more than 1
		v = v3f(sinTheta * std::cosf(phi), sinTheta * std::sinf(phi), cosTheta);
		return v;
	}

	inline float squareToUniformHemispherePdf(const v3f& v) {
		float pdf = 0.f;
		// TODO: Implement this
		pdf = INV_TWOPI;
		return pdf;
	}

	inline v2f squareToUniformDiskConcentric(const p2f& sample) {
		v2f v(0.f);
		// TODO: Implement this (optional)
		float phi, radius;
		float remappedX = (2.f * sample.x) - 1.f;
		float remappedY = (2.f * sample.y) - 1.f;

		//if only one of these is 0, we wont enter the branch where we will need to divide by it, so it is safe
		if (remappedX == 0 && remappedY == 0)
		{
			return v2f(0.f, 0.f);
		}

		if ((remappedX * remappedX) > (remappedY * remappedY)) {
			// Right and left triangles
			radius = remappedX;
			phi = (M_PI * 0.25f) * (remappedY * (1.f / remappedX));
		}
		else {
			// Top and bottom triangles
			radius = remappedY;
			// Minus instead of addition of x/y component to fix a clumping issue with hammersley sequence, need to investigate it though
			phi = (M_PI * 0.5f) - ((M_PI * 0.25f) * (remappedX * (1.f / remappedY)));
		}

		v = v2f(radius * std::cosf(phi), radius * std::sinf(phi));
		return v;
	}

	inline v3f squareToCosineHemisphere(const p2f& sample) {
		v3f v(0.f);
		// TODO: Implement this
		v2f discSample = squareToUniformDiskConcentric(sample);
		float z = 1.0f - glm::length2(discSample);
		z = glm::fmax(z, 0.f);
		z = std::sqrtf(z);

		v = v3f(discSample.x, discSample.y, z);
		return v;
	}

	inline float squareToCosineHemispherePdf(const v3f& v) {
		float pdf = 0.f;
		// TODO: Implement this
		//Remember, cos(theta) is dot product of v with z axis.
		//0 inclusive, as points on the plane can be produced by the distribution, what is the best way to deal with this?
		if (v.z >= 0.f)
		{
			pdf = v.z * INV_PI;
		}
		else
		{
			pdf = 0.f;
		}
		return pdf;
	}

	inline v3f squareToPhongLobe(const p2f& sample, const float& exponent) {
		v3f v(0.f);
		// TODO: Implement this
		float cosTheta = std::powf(sample.x, 1.f / (exponent + 2));
		float sinTheta = std::sqrtf(std::fmax(1.f - (cosTheta * cosTheta), 0));
		float phi = sample.y * 2.f * M_PI;

		v = v3f(sinTheta * std::cosf(phi), sinTheta * std::sinf(phi), cosTheta);
		return v;
	}

	inline float squareToPhongLobePdf(const v3f& v, const float &exponent) {
		float pdf = 0.f;
		// TODO: Implement this
		pdf = v.z >= 0.f ? (exponent + 2) * INV_TWOPI * std::powf(v.z, exponent) : 0.f;

		return pdf;
	}

	inline v2f squareToUniformTriangle(const p2f& sample) {
		v2f v(0.f);
		float u = std::sqrt(1.f - sample.x);
		v = { 1 - u, u * sample.y };
		return v;
	}

	inline v3f squareToUniformCone(const p2f& sample, float cosThetaMax) {
		v3f v(0.f);
		// TODO: Add previous assignment code (if needed)
		float cosTheta = (1.f - sample.x) + sample.x * cosThetaMax;
		float phi = sample.y * M_PI * 2.0f;
		float sinTheta = std::sqrtf(std::fmax(1.f - (cosTheta * cosTheta), 0)); //fmax is just a precaution, cosTheta*cosTheta should never be more than 1
		v = v3f(sinTheta * std::cosf(phi), sinTheta * std::sinf(phi), cosTheta);

		return v;
	}

	inline float squareToUniformConePdf(float cosThetaMax) {
		float pdf = 0.f;
		// TODO: Add previous assignment code (if needed)
		pdf = 1.f / (1.f - cosThetaMax);
		pdf *= INV_TWOPI;

		return pdf;
	}

inline p2f squareToUniformDisk(const p2f& sample) {
    p2f p(0.f);
    // TODO: Add previous assignment code (if needed)
    return p;
}

inline float squareToUniformDiskPdf(const p2f& p) {
    float pdf = 0.f;
    // TODO: Add previous assignment code (if needed)
    return pdf;
}

}

TR_NAMESPACE_END