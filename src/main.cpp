/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#include <core/core.h>
#include <core/platform.h>
#include <core/renderer.h>
#define TINYEXR_IMPLEMENTATION
#include "tinyexr.h"
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include <chrono>

std::unique_ptr<std::mutex[]> g_FrameBufferLocks;

/**
 * Load TOML scene file and create scene objects.
 */
bool loadTOML(TinyRender::Config& config, const std::string& inputFile) {
    // Scene and Wavefront OBJ files
    const auto data = cpptoml::parse_file(inputFile);
    config.tomlFile = inputFile;
    const auto input = data->get_table("input");
    config.objFile = *input->get_as<std::string>("objfile");

    // Camera settings
    const auto camera = data->get_table("camera");
    config.camera.fov = camera->get_as<double>("fov").value_or(30.);
    auto eye = camera->get_array_of<double>("eye").value_or({1., 1., 0.});
    config.camera.o = v3f(eye[0], eye[1], eye[2]);
    auto at = camera->get_array_of<double>("at").value_or({0., 0., 0.});
    config.camera.at = v3f(at[0], at[1], at[2]);
    auto up = camera->get_array_of<double>("up").value_or({0., 1., 0.});
    config.camera.up = v3f(up[0], up[1], up[2]);

    // Film settings
    const auto film = data->get_table("film");
    config.width = film->get_as<int>("width").value_or(768);
    config.height = film->get_as<int>("height").value_or(576);

    // Renderer settings
    const auto renderer = data->get_table("renderer");
    auto realTime = renderer->get_as<bool>("realtime").value_or(false);
    auto type = renderer->get_as<std::string>("type").value_or("normal");

    // Real-time renderpass
    if (realTime) {
        if (type == "normal") {
            config.renderpass = TinyRender::ENormalRenderPass;
        }
        else if (type == "direct") {
            config.renderpass = TinyRender::EDirectRenderPass;
        }
        else if (type == "ssao") {
            config.renderpass = TinyRender::ESSAORenderPass;
        }
        else if (type == "gi") {
            config.renderpass = TinyRender::EGIRenderPass;
            config.integratorSettings.gi.maxDepth = renderer->get_as<int>("maxDepth").value_or(5);
            config.integratorSettings.gi.rrDepth = renderer->get_as<int>("rrDepth").value_or(5);
            config.integratorSettings.gi.rrProb = renderer->get_as<double>("rrProb").value_or(9.95f);
            config.integratorSettings.gi.samplesByVertex = renderer->get_as<int>("samplesByVertex").value_or(100);
        }
        else {
            throw std::runtime_error("Invalid renderpass type");
        }
    }

    // Offline integrator
    else {
		std::cout << type << std::endl;
        if (type == "normal") {
            config.integrator = TinyRender::ENormalIntegrator;
        }
        else if (type == "simple") {
            config.integrator = TinyRender::ESimpleIntegrator;
        }
        else if (type == "ao") {
            config.integrator = TinyRender::EAOIntegrator;
        }
        else if (type == "ro") {
            config.integrator = TinyRender::EROIntegrator;
            config.integratorSettings.ro.exponent = renderer->get_as<double>("exponent").value_or(30);
        }
        else if (type == "direct") {
            config.integrator = TinyRender::EDirectIntegrator;
            config.integratorSettings.di.emitterSamples = renderer->get_as<size_t>("emitterSamples").value_or(1);
            config.integratorSettings.di.bsdfSamples = renderer->get_as<size_t>("bsdfSamples").value_or(1);
            config.integratorSettings.di.samplingStrategy = renderer->get_as<string>("samplingStrategy").value_or("emitter");
        }
        else if (type == "path") {
            config.integrator = TinyRender::EPathTracerIntegrator;
            config.integratorSettings.pt.isExplicit = renderer->get_as<bool>("isExplicit").value_or(true);
            config.integratorSettings.pt.maxDepth = renderer->get_as<int>("maxDepth").value_or(-1);
            config.integratorSettings.pt.rrDepth = renderer->get_as<int>("rrDepth").value_or(5);
            config.integratorSettings.pt.rrProb = renderer->get_as<double>("rrProb").value_or(0.95f);
			config.integratorSettings.pt.emitterSamples = renderer->get_as<size_t>("emitterSamples").value_or(1);
			config.integratorSettings.pt.bsdfSamples = renderer->get_as<size_t>("bsdfSamples").value_or(0);
        }
		else if (type == "bdpt") {
			config.integrator = TinyRender::EBDPTIntegrator;
			config.integratorSettings.pt.rrDepth = renderer->get_as<int>("rrDepth").value_or(5);
			config.integratorSettings.pt.rrProb = renderer->get_as<double>("rrProb").value_or(0.f);
		}
        else {
            throw std::runtime_error("Invalid integrator type");
        }

        config.spp = renderer->get_as<int>("spp").value_or(1);
    }

    return realTime;
}

/**
 * Launch rendering job.
 */
void run(std::string& inputTOMLFile, bool nogui) {
    TinyRender::Config config;
    bool isRealTime;

    try {
        isRealTime = loadTOML(config, inputTOMLFile);
    } catch (std::exception const& e) {
        std::cerr << "Error while parsing scene file: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }

    TinyRender::Renderer renderer(config);
    renderer.init(isRealTime, nogui);

	//I didnt want to put this here, but const qualification on the integrator functions is a design choice which locked me to this solution

	size_t frameBufferSize = renderer.scene.config.width * renderer.scene.config.height;

	g_FrameBufferLocks.reset(new std::mutex[frameBufferSize]);

	for (int i = 0; i < frameBufferSize; i++) {
		new (&g_FrameBufferLocks[i]) std::mutex();
	}

	//https://www.pluralsight.com/blog/software-development/how-to-measure-execution-time-intervals-in-c--
	auto start = std::chrono::high_resolution_clock::now();

    renderer.render();

	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Render took: " << elapsed.count() << " seconds." << std::endl;

    renderer.cleanUp();
}

/**
 * Main TinyRender program.
 */
int main(int argc, char* argv[]) {
    if (argc != 2 && argc !=3) {
        cerr << "Syntax: " << argv[0] << " <scene.toml>" << endl;
        exit(EXIT_FAILURE);
    }

    bool nogui = false;
    if(argc == 3) {
        if(std::string(argv[2]) == "nogui") {
            nogui = true;
        }
    }

    auto inputTOMLFile = std::string(argv[1]);
    run(inputTOMLFile, nogui);

#ifdef _WIN32
    if(!nogui) system("pause");
#endif

    return EXIT_SUCCESS;
}
