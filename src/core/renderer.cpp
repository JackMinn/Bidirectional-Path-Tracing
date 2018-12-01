/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#include <core/core.h>
#include <core/accel.h>
#include <core/renderer.h>
#include <GL/glew.h>

#ifdef __APPLE__
#include "SDL.h"
#include <OpenGL/gl.h>
#else
#ifdef _WIN32
#include <GL/gl.h>
#include "SDL.h"
#else
#include <GL/gl.h>
#include "SDL2/SDL.h"
#endif
#endif

#include <bsdfs/diffuse.h>

#include <integrators/normal.h>
#include <renderpasses/normal.h>

#include <integrators/simple.h>
#include <renderpasses/simple.h>
#include <bsdfs/phong.h>

#include <integrators/ao.h>
#include <integrators/ro.h>
#include <renderpasses/ssao.h>

#include <integrators/direct.h>

#include <integrators/path.h>
#include <renderpasses/gi.h>
#include <bsdfs/mixture.h>

#include <integrators/bdpt.h>
#include <bsdfs/perfectmirror.h>
#include <bsdfs/glass.h>

extern std::unique_ptr<std::mutex[]> g_FrameBufferLocks;

TR_NAMESPACE_BEGIN

Renderer::Renderer(const Config& config) : scene(config) { }

bool Renderer::init(const bool isRealTime, bool nogui) {
    realTime = isRealTime;
    this->nogui = nogui;
    realTimeCameraFree = false;

    if (!scene.load(isRealTime)) return false;

    if (realTime) {
        if (scene.config.renderpass == ENormalRenderPass) {
            renderpass = std::unique_ptr<NormalPass>(new NormalPass(scene));
        }
        else if (scene.config.renderpass == EDirectRenderPass) {
            renderpass = std::unique_ptr<SimplePass>(new SimplePass(scene));
        }
        else if (scene.config.renderpass == ESSAORenderPass) {
            renderpass = std::unique_ptr<SSAOPass>(new SSAOPass(scene));
        }
        else if (scene.config.renderpass == EGIRenderPass) {
            renderpass = std::unique_ptr<GIPass>(new GIPass(scene));
        }
        else {
            throw std::runtime_error("Invalid renderpass type");
        }

        bool succ = renderpass.get()->initOpenGL(scene.config.width, scene.config.height);
        if (!succ) return false;

        return renderpass->init(scene.config);
    } else {
        if (scene.config.integrator == ENormalIntegrator) {
            integrator = std::unique_ptr<NormalIntegrator>(new NormalIntegrator(scene));
        }
        else if (scene.config.integrator == EAOIntegrator) {
            integrator = std::unique_ptr<AOIntegrator>(new AOIntegrator(scene));
        } else if (scene.config.integrator == EROIntegrator) {
            integrator = std::unique_ptr<ROIntegrator>(new ROIntegrator(scene));
        }
        else if (scene.config.integrator == ESimpleIntegrator) {
            integrator = std::unique_ptr<SimpleIntegrator>(new SimpleIntegrator(scene));
        }
        else if (scene.config.integrator == EDirectIntegrator) {
            integrator = std::unique_ptr<DirectIntegrator>(new DirectIntegrator(scene));
        }
        else if (scene.config.integrator == EPathTracerIntegrator) {
            integrator = std::unique_ptr<PathTracerIntegrator>(new PathTracerIntegrator(scene));
        }
		else if (scene.config.integrator == EBDPTIntegrator) {
			integrator = std::unique_ptr<BDPTIntegrator>(new BDPTIntegrator(scene));
		}
        else {
            throw std::runtime_error("Invalid integrator type");
        }

        return integrator->init();
    }
}

void Renderer::render() {
	if (realTime) {
		/**
		 * 1) Detect and handle the quit event.
		 * 2) Call the render function using renderpass->render().
		 * 3) Output the rendered image into the GUI window using SDL_GL_SwapWindow(renderpass->window).
		 */
		 // TODO: Implement this

		SDL_Event event = SDL_Event();
		while (event.type != SDL_QUIT) {
			SDL_PollEvent(&event);

			renderpass->updateCamera(event);
			renderpass->render();
			SDL_GL_SwapWindow(renderpass->window);
		}
	}
	else {
		/**
		 * 1) Calculate the camera perspective, the camera-to-world transformation matrix and the aspect ratio.
		 * 2) Clear integral RGB buffer.
		 * 3) Loop over all pixels on the image plane.
		 * 4) Generate a ray through each pixel center.
		 * 5) Splat their contribution onto the image plane.
		 */
		 // TODO: Implement this

		const float near = 1.f;
		const float far = 1000.f;

		auto worldToCamera = glm::lookAt(scene.config.camera.o, scene.config.camera.at, scene.config.camera.up);
		auto cameraToWorld = glm::inverse(worldToCamera);

		float invWidth = 1.f / integrator->rgb->width;
		float invHeight = 1.f / integrator->rgb->height;

		float angle = std::tanf(deg2rad * scene.config.camera.fov * 0.5f);
		float aspectRatio = (float)integrator->rgb->width / integrator->rgb->height;

		auto cameraToClip = glm::perspective(deg2rad * scene.config.camera.fov, aspectRatio, near, far);
		auto NDCToScreen = glm::scale(glm::mat4(1.f), v3f(integrator->rgb->width, integrator->rgb->height, 1.0f)) * glm::scale(glm::mat4(1.f), v3f(0.5f, -0.5f, 1.f)) * glm::translate(glm::mat4(1.f), v3f(1.f, -1.f, 0.f));

		Sampler sampler = Sampler(260450963); //I could converge faster if I used a quasi random sampler perhaps (Sobol sequence?)
		const unsigned int pixelSize = integrator->rgb->height * integrator->rgb->width;

		bool bMultithread = true;
		parallel_for(pixelSize, [&](int start, int end) {
			//We could maybe construct a duplicate of the sampler here, so each batch has its own sampler and no multithreaded race conditions exist...
			for (int pixel = start; pixel < end; ++pixel) {
				int j = pixel % integrator->rgb->width;
				int i = pixel / integrator->rgb->width;

				float y = (1.f - ((float)i + 0.5f) * invHeight) * 2.f - 1.f;
				float x = (((float)j + 0.5f) * invWidth) * 2.f - 1.f;

				v3f averageRadiance = v3f(0.f);
				if (scene.config.spp == 1)
				{
					//The FoV and the aspect ratio scale the image plane. The image plane is a distance of the near plane away, 
					//but while it is originally from [-1,1] in both dimensions, having a more narrow FoV will produce a smaller 
					//angle value (tan of a smaller value), and thus make the image plane smaller in world space. 
					v4f imagePlanePoint = v4f(x * angle * aspectRatio, y * angle, -near, 0);
					v3f rayDir = cameraToWorld * imagePlanePoint;
					rayDir = glm::normalize(rayDir);
					Ray ray = Ray(scene.config.camera.o, rayDir, near, far);

					averageRadiance += integrator->render(ray, sampler);
				}
				else
				{
					for (int s = 0; s < scene.config.spp; s++) {
						p2f randomSample = sampler.next2D();
						randomSample -= 0.5f;
						randomSample.x = randomSample.x * invWidth;
						randomSample.y = randomSample.y * invHeight;

						v4f imagePlanePoint = v4f((x + randomSample.x) * angle * aspectRatio, (y + randomSample.y) * angle, -near, 0);
						v3f rayDir = cameraToWorld * imagePlanePoint;
						rayDir = glm::normalize(rayDir);
						Ray ray = Ray(scene.config.camera.o, rayDir, near, far);

						averageRadiance += integrator->render(ray, sampler);
					}
				}

				if (bMultithread) {
					while (true) {
						std::mutex& lock = g_FrameBufferLocks[pixel];
						if (lock.try_lock()) {
							integrator->rgb->data[pixel] += averageRadiance * (1.f / scene.config.spp);
							lock.unlock();
							break;
						}
					}
				}
				else {
					integrator->rgb->data[pixel] += averageRadiance * (1.f / scene.config.spp);
				}

			}
		}, bMultithread);
	}
}


/**
 * Post-rendering step.
 */
void Renderer::cleanUp() {
    if (realTime) {
        renderpass->cleanUp();
    } else {
        integrator->cleanUp();
    }
}

BSDF::BSDF(const WorldData& d, const Config& c, const size_t matID) : worldData(d), config(c) {
    emission = glm::make_vec3(worldData.materials[matID].emission);
}

Scene::Scene(const Config& config) : config(config) { }

bool Scene::load(bool isRealTime) {
    fs::path file(config.objFile);
    bool ret = false;
    std::string err;

    if (!file.is_absolute())
        file = (config.tomlFile.parent_path() / file).make_preferred();

    tinyobj::attrib_t* attrib_ = &worldData.attrib;
    std::vector<tinyobj::shape_t>* shapes_ = &worldData.shapes;
    std::vector<tinyobj::material_t>* materials_ = &worldData.materials;
    std::string* err_ = &err;
    const string filename_ = file.string();
    const string mtl_basedir_ = file.make_preferred().parent_path().string();
    ret = tinyobj::LoadObj(attrib_, shapes_, materials_, err_, filename_.c_str(), mtl_basedir_.c_str(), true);

    if (!err.empty()) { std::cout << "Error: " << err.c_str() << std::endl; }
    if (!ret) {
        std::cout << "Failed to load scene " << config.objFile << " " << std::endl;
        return false;
    }

    // Build list of BSDFs
    bsdfs = std::vector<std::unique_ptr<BSDF>>(worldData.materials.size());
    for (size_t i = 0; i < worldData.materials.size(); i++) {
        if (worldData.materials[i].illum == 7)
            bsdfs[i] = std::unique_ptr<BSDF>(new DiffuseBSDF(worldData, config, i));
        if (worldData.materials[i].illum == 3)
            bsdfs[i] = std::unique_ptr<BSDF>(new MirrorBSDF(worldData, config, i));
		if (worldData.materials[i].illum == 6)
			bsdfs[i] = std::unique_ptr<BSDF>(new GlassBSDF(worldData, config, i));
        if (worldData.materials[i].illum == 8)
            bsdfs[i] = std::unique_ptr<BSDF>(new MixtureBSDF(worldData, config, i));
		if (worldData.materials[i].illum != 5 && worldData.materials[i].illum != 7 && worldData.materials[i].illum != 8 &&
			worldData.materials[i].illum != 3 && worldData.materials[i].illum != 6)
			bsdfs[i] = std::unique_ptr<BSDF>(new PhongBSDF(worldData, config, i));
    }

    // Build list of emitters (and print what has been loaded)
    std::string nbShapes = worldData.shapes.size() > 1 ? " shapes" : " shape";
    std::cout << "Found " << worldData.shapes.size() << nbShapes << std::endl;
    worldData.shapesCenter.resize(worldData.shapes.size());
    worldData.shapesAABOX.resize(worldData.shapes.size());

    for (size_t i = 0; i < worldData.shapes.size(); i++) {
        const tinyobj::shape_t& shape = worldData.shapes[i];
        const BSDF* bsdf = bsdfs[shape.mesh.material_ids[0]].get();
        std::cout << "Mesh " << i << ": " << shape.name << " ["
                  << shape.mesh.indices.size() / 3 << " primitives | ";

        if (bsdf->isEmissive()) {
            Distribution1D faceAreaDistribution;
            float shapeArea = getShapeArea(i, faceAreaDistribution);
            emitters.emplace_back(Emitter{i, shapeArea, bsdf->emission, faceAreaDistribution});
            std::cout << "Emitter]" << std::endl;
        } else {
            std::cout << bsdf->toString() << "]" << std::endl;
        }

        // Build world AABB and shape centers
        worldData.shapesCenter[i] = v3f(0.0);
        for (auto idx: shape.mesh.indices) {
            v3f p = {worldData.attrib.vertices[3 * idx.vertex_index + 0],
                     worldData.attrib.vertices[3 * idx.vertex_index + 1],
                     worldData.attrib.vertices[3 * idx.vertex_index + 2]};
            worldData.shapesCenter[i] += p;
            worldData.shapesAABOX[i].expandBy(p);
            aabb.expandBy(p);
        }
        worldData.shapesCenter[i] /= float(shape.mesh.indices.size());
    }

    // Build BVH
    bvh = std::unique_ptr<TinyRender::AcceleratorBVH>(new TinyRender::AcceleratorBVH(this->worldData));

    const clock_t beginBVH = clock();
    bvh->build();
    std::cout << "BVH built in " << float(clock() - beginBVH) / CLOCKS_PER_SEC << "s" << std::endl;

    return true;
}

float Scene::getShapeArea(const size_t shapeID, Distribution1D& faceAreaDistribution) {
    const tinyobj::shape_t& s = worldData.shapes[shapeID];

    for (size_t i = 0; i < s.mesh.indices.size(); i += 3) {
        const int i0 = s.mesh.indices[i + 0].vertex_index;
        const int i1 = s.mesh.indices[i + 1].vertex_index;
        const int i2 = s.mesh.indices[i + 2].vertex_index;
        const v3f v0{worldData.attrib.vertices[3 * i0 + 0], worldData.attrib.vertices[3 * i0 + 1],
                     worldData.attrib.vertices[3 * i0 + 2]};
        const v3f v1{worldData.attrib.vertices[3 * i1 + 0], worldData.attrib.vertices[3 * i1 + 1],
                     worldData.attrib.vertices[3 * i1 + 2]};
        const v3f v2{worldData.attrib.vertices[3 * i2 + 0], worldData.attrib.vertices[3 * i2 + 1],
                     worldData.attrib.vertices[3 * i2 + 2]};

        const v3f e1{v1 - v0};
        const v3f e2{v2 - v0};
        const v3f e3{glm::cross(e1, e2)};
        faceAreaDistribution.add(0.5f * std::sqrt(e3.x * e3.x + e3.y * e3.y + e3.z * e3.z));
    }
    const float area = faceAreaDistribution.cdf.back();
    faceAreaDistribution.normalize();
    return area;
}

v3f Scene::getFirstLightPosition() const {
    return worldData.shapesCenter[emitters[0].shapeID];
}

v3f Scene::getFirstLightIntensity() const {
    return emitters[0].getRadiance(); // point lights are defined by intensity not radiance
}

float Scene::getShapeRadius(const size_t shapeID) const {
    assert(shapeID < worldData.shapes.size());
    v3f emitterCenter = worldData.shapesCenter[shapeID];
    return worldData.shapesAABOX[shapeID].max.x - emitterCenter.x;
}

v3f Scene::getShapeCenter(const size_t shapeID) const {
    assert(shapeID < worldData.shapes.size());
    return worldData.shapesCenter[shapeID];
}

size_t Scene::getFirstLight() const {
    if (emitters.size() <= 0) return -1;
    return emitters[0].shapeID;
}

v3f Scene::getObjectVertexPosition(size_t objectIdx, size_t vertexIdx) const {
    const tinyobj::attrib_t& sa = worldData.attrib;
    const tinyobj::shape_t& s = worldData.shapes[objectIdx];

    int idx = s.mesh.indices[vertexIdx].vertex_index;
    float x = sa.vertices[3 * idx + 0];
    float y = sa.vertices[3 * idx + 1];
    float z = sa.vertices[3 * idx + 2];
    return v3f(x,y,z);
}

v3f Scene::getObjectVertexNormal(size_t objectIdx, size_t vertexIdx) const {
    const tinyobj::attrib_t& sa = worldData.attrib;
    const tinyobj::shape_t& s = worldData.shapes[objectIdx];

    int idx_n = s.mesh.indices[vertexIdx].normal_index;
    float nx = sa.normals[3 * idx_n + 0];
    float ny = sa.normals[3 * idx_n + 1];
    float nz = sa.normals[3 * idx_n + 2];
    return glm::normalize(v3f(nx,ny,nz));
}

size_t Scene::getObjectNbVertices(size_t objectIdx) const {
    return worldData.shapes[objectIdx].mesh.indices.size();
}

int Scene::getPrimitiveID(size_t vertexIdx) const {
    return vertexIdx / 3;
}

int Scene::getMaterialID(size_t objectIdx, int primID) const {
    return worldData.shapes[objectIdx].mesh.material_ids[primID];
}

TR_NAMESPACE_END
