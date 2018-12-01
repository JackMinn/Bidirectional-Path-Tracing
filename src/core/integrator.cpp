/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#include <core/integrator.h>

#include "tiny_obj_loader.h"

TR_NAMESPACE_BEGIN

Integrator::Integrator(const Scene& scene) : scene(scene) { }

bool Integrator::init() {
    rgb = std::unique_ptr<RenderBuffer>(new RenderBuffer(scene.config.width, scene.config.height));
    rgb->clear();
    return true;
}

void Integrator::cleanUp() {
    save();
}

bool Integrator::save() {
    fs::path p = scene.config.tomlFile;
    saveEXR(rgb->data, p.replace_extension("exr").string(), scene.config.width, scene.config.height);
    return true;
}

const Emitter& Integrator::getEmitterByID(const int emitterID) const {
    return scene.emitters[emitterID];
}
const BSDF* Integrator::getBSDF(const SurfaceInteraction& hit) const {
    assert(hit.shapeID < scene.worldData.shapes.size());
    const tinyobj::shape_t& shape = scene.worldData.shapes[hit.shapeID];
    return scene.bsdfs[shape.mesh.material_ids[hit.primID]].get();
}

v3f Integrator::getEmission(const SurfaceInteraction& hit) const {
    assert(hit.matID < scene.worldData.materials.size());
    return glm::make_vec3(scene.worldData.materials[hit.matID].emission);
}

size_t Integrator::selectEmitter(float sample, float& pdf) const {
    size_t id = size_t(sample * scene.emitters.size());
    id = min(id, scene.emitters.size() - 1); //todo @nico : how can this happen ? (sample ==1)
    pdf = 1.f / scene.emitters.size();
    return id;
}

size_t Integrator::getEmitterIDByShapeID(size_t shapeID) const {
    auto it = std::find_if(scene.emitters.begin(), scene.emitters.end(),
                           [shapeID](const Emitter& s) { return s.shapeID == shapeID; });
    assert(it != scene.emitters.end());
    return (size_t) std::distance(scene.emitters.begin(), it);
}

float Integrator::getEmitterPdf(const Emitter& emitter) const {
    return 1.f / scene.emitters.size();
}

void Integrator::sampleEmitterDirection(Sampler& sampler,
                                        const Emitter& emitter,
                                        const v3f& n,
                                        v3f& d,
                                        float& pdf) const {
    // TODO: Add previous assignment code (if needed)
	
}

void Integrator::sampleEmitterPosition(Sampler& sampler, const Emitter& emitter, v3f& n, v3f& pos, float& pdf, size_t* primitiveID /* = nullptr */) const {
	const tinyobj::shape_t& shape = scene.worldData.shapes[emitter.shapeID];
	const size_t primID = (size_t)emitter.faceAreaDistribution.sample(sampler.next());
	const v2f uv = Warp::squareToUniformTriangle(sampler.next2D());

	if (primitiveID)
		*primitiveID = primID;

	const tinyobj::index_t& idx0 = shape.mesh.indices[3 * primID + 0];
	const tinyobj::index_t& idx1 = shape.mesh.indices[3 * primID + 1];
	const tinyobj::index_t& idx2 = shape.mesh.indices[3 * primID + 2];

	auto& vx = scene.worldData.attrib.vertices;
	const v3f v0{ vx[3 * idx0.vertex_index + 0], vx[3 * idx0.vertex_index + 1], vx[3 * idx0.vertex_index + 2] };
	const v3f v1{ vx[3 * idx1.vertex_index + 0], vx[3 * idx1.vertex_index + 1], vx[3 * idx1.vertex_index + 2] };
	const v3f v2{ vx[3 * idx2.vertex_index + 0], vx[3 * idx2.vertex_index + 1], vx[3 * idx2.vertex_index + 2] };

	pos = barycentric(v0, v1, v2, uv.x, uv.y);

	auto& ns = scene.worldData.attrib.normals;
	const v3f n0{ ns[3 * idx0.normal_index + 0], ns[3 * idx0.normal_index + 1], ns[3 * idx0.normal_index + 2] };
	const v3f n1{ ns[3 * idx1.normal_index + 0], ns[3 * idx1.normal_index + 1], ns[3 * idx1.normal_index + 2] };
	const v3f n2{ ns[3 * idx2.normal_index + 0], ns[3 * idx2.normal_index + 1], ns[3 * idx2.normal_index + 2] };

	n = glm::normalize(barycentric(n0, n1, n2, uv.x, uv.y));

	pdf = 1.f / emitter.area;
}


TR_NAMESPACE_END