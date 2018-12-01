/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include "core.h"
#include "bvh.h"

TR_NAMESPACE_BEGIN

/**
 * Bounding-volume hierarchy (BVH) acceleration structure.
 */
struct AcceleratorBVH {

    struct BVHNode : Object {

        const size_t shapeID, faceID;
        const WorldData& worldData;

        BVHNode(size_t j, size_t i, const WorldData& d) : shapeID(j), faceID(i), worldData(d) { }

        bool getIntersection(const Ray& ray, IntersectionInfo* intersection) const override {
            const tinyobj::attrib_t& a = worldData.attrib;
            const tinyobj::mesh_t& m = worldData.shapes[shapeID].mesh;
            const tinyobj::index_t& idx0 = m.indices[faceID + 0];
            const tinyobj::index_t& idx1 = m.indices[faceID + 1];
            const tinyobj::index_t& idx2 = m.indices[faceID + 2];

            const v3f v0 = {a.vertices[3 * idx0.vertex_index + 0], a.vertices[3 * idx0.vertex_index + 1],
                            a.vertices[3 * idx0.vertex_index + 2]};
            const v3f v1 = {a.vertices[3 * idx1.vertex_index + 0], a.vertices[3 * idx1.vertex_index + 1],
                            a.vertices[3 * idx1.vertex_index + 2]};
            const v3f v2 = {a.vertices[3 * idx2.vertex_index + 0], a.vertices[3 * idx2.vertex_index + 1],
                            a.vertices[3 * idx2.vertex_index + 2]};

            float t, u, v;
            if (rayTriangleIntersect(ray, v0, v1, v2, t, u, v)) {
                if (t > 1e-3) {
                    intersection->t = t;
                    intersection->u = u;
                    intersection->v = v;
                    intersection->object = this;
                    return true;
                }
            }
            return false;
        }

        v3f getNormal(const IntersectionInfo&) const override {
            const tinyobj::mesh_t& m = worldData.shapes[shapeID].mesh;
            const tinyobj::index_t& idx0 = m.indices[faceID + 0];
            const tinyobj::index_t& idx1 = m.indices[faceID + 1];
            const tinyobj::index_t& idx2 = m.indices[faceID + 2];

            const tinyobj::attrib_t& a = worldData.attrib;
            const v3f v0 = {a.normals[3 * idx0.normal_index + 0], a.normals[3 * idx0.normal_index + 1],
                            a.normals[3 * idx0.normal_index + 2]};
            const v3f v1 = {a.normals[3 * idx1.normal_index + 0], a.normals[3 * idx1.normal_index + 1],
                            a.normals[3 * idx1.normal_index + 2]};
            const v3f v2 = {a.normals[3 * idx2.normal_index + 0], a.normals[3 * idx2.normal_index + 1],
                            a.normals[3 * idx2.normal_index + 2]};

            return glm::normalize(glm::cross(v1 - v0, v2 - v0));
        }

        BBox getBBox() const override {
            const tinyobj::mesh_t& m = worldData.shapes[shapeID].mesh;
            const tinyobj::index_t& idx0 = m.indices[faceID + 0];
            const tinyobj::index_t& idx1 = m.indices[faceID + 1];
            const tinyobj::index_t& idx2 = m.indices[faceID + 2];

            const tinyobj::attrib_t& a = worldData.attrib;
            const v3f v0 = {a.vertices[3 * idx0.vertex_index + 0], a.vertices[3 * idx0.vertex_index + 1],
                            a.vertices[3 * idx0.vertex_index + 2]};
            const v3f v1 = {a.vertices[3 * idx1.vertex_index + 0], a.vertices[3 * idx1.vertex_index + 1],
                            a.vertices[3 * idx1.vertex_index + 2]};
            const v3f v2 = {a.vertices[3 * idx2.vertex_index + 0], a.vertices[3 * idx2.vertex_index + 1],
                            a.vertices[3 * idx2.vertex_index + 2]};

            BBox b(v0);
            b.expandToInclude(v1);
            b.expandToInclude(v2);
            return b;
        }

        v3f getCentroid() const override {
            const tinyobj::mesh_t& m = worldData.shapes[shapeID].mesh;
            const tinyobj::index_t& idx0 = m.indices[faceID + 0];
            const tinyobj::index_t& idx1 = m.indices[faceID + 1];
            const tinyobj::index_t& idx2 = m.indices[faceID + 2];

            const tinyobj::attrib_t& a = worldData.attrib;
            const v3f v0 = {a.vertices[3 * idx0.vertex_index + 0], a.vertices[3 * idx0.vertex_index + 1],
                            a.vertices[3 * idx0.vertex_index + 2]};
            const v3f v1 = {a.vertices[3 * idx1.vertex_index + 0], a.vertices[3 * idx1.vertex_index + 1],
                            a.vertices[3 * idx1.vertex_index + 2]};
            const v3f v2 = {a.vertices[3 * idx2.vertex_index + 0], a.vertices[3 * idx2.vertex_index + 1],
                            a.vertices[3 * idx2.vertex_index + 2]};

            return (v0 + v1 + v2) / 3.0f;
        }
    };

    std::unique_ptr<BVH> bvh;
    std::vector<Object*> objects;
    const WorldData& worldData;

    explicit AcceleratorBVH(const WorldData& worldData) : worldData(worldData) { }

    bool build() {
        for (size_t j = 0; j < worldData.shapes.size(); j++) {
            const tinyobj::shape_t& shape = worldData.shapes[j];
            for (size_t i = 0; i < shape.mesh.indices.size(); i += 3)
                objects.emplace_back(new BVHNode(j, i, worldData));
        }
        bvh = std::unique_ptr<BVH>(new BVH(&objects));
        return true;
    }

    bool intersect(const Ray& ray, SurfaceInteraction& info, bool occlusion = false) const {
        IntersectionInfo iInfo{};
        iInfo.object = nullptr;
        const std::vector<tinyobj::shape_t>& ss = worldData.shapes;
        const tinyobj::attrib_t& sa = worldData.attrib;

        if (bvh->getIntersection(ray, &iInfo, false)) {
            info.t = iInfo.t;
            if (iInfo.t <= ray.max_t && iInfo.t >= ray.min_t) {


                const tinyobj::shape_t& s = ss[((BVHNode*) (iInfo.object))->shapeID];
                const size_t i = ((BVHNode*) (iInfo.object))->faceID;
                const tinyobj::index_t& idx0 = s.mesh.indices[i + 0];
                const tinyobj::index_t& idx1 = s.mesh.indices[i + 1];
                const tinyobj::index_t& idx2 = s.mesh.indices[i + 2];

                const v3f v0 = {sa.vertices[3 * idx0.vertex_index + 0], sa.vertices[3 * idx0.vertex_index + 1],
                                sa.vertices[3 * idx0.vertex_index + 2]};
                const v3f v1 = {sa.vertices[3 * idx1.vertex_index + 0], sa.vertices[3 * idx1.vertex_index + 1],
                                sa.vertices[3 * idx1.vertex_index + 2]};
                const v3f v2 = {sa.vertices[3 * idx2.vertex_index + 0], sa.vertices[3 * idx2.vertex_index + 1],
                                sa.vertices[3 * idx2.vertex_index + 2]};

                const v3f n0 = {sa.normals[3 * idx0.normal_index + 0], sa.normals[3 * idx0.normal_index + 1],
                                sa.normals[3 * idx0.normal_index + 2]};
                const v3f n1 = {sa.normals[3 * idx1.normal_index + 0], sa.normals[3 * idx1.normal_index + 1],
                                sa.normals[3 * idx1.normal_index + 2]};
                const v3f n2 = {sa.normals[3 * idx2.normal_index + 0], sa.normals[3 * idx2.normal_index + 1],
                                sa.normals[3 * idx2.normal_index + 2]};

                info.shapeID = ((BVHNode*) (iInfo.object))->shapeID;
                info.primID = ((BVHNode*) (iInfo.object))->faceID / 3;
                info.t = iInfo.t;
                info.u = iInfo.u;
                info.v = iInfo.v;
                info.p = barycentric(v0, v1, v2, iInfo.u, iInfo.v);
                info.frameNg = Frame(glm::normalize(glm::cross(v1 - v0, v2 - v0)));
                info.frameNs = Frame(glm::normalize(barycentric(n0, n1, n2, info.u, info.v)));
                info.wo = info.frameNs.toLocal(-ray.d);
                info.matID = s.mesh.material_ids[info.primID];
                return true;
            }
            return false;
        }
        info.t = std::numeric_limits<float>::max();
        return false;
    }
};

TR_NAMESPACE_END