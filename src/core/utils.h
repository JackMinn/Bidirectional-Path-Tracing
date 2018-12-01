/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <GL/glew.h>
#include <core/platform.h>
#include <tinyformat.h>
#include "tinyexr.h"
#include <iterator>
#include <iostream>
#include <iomanip>

TR_NAMESPACE_BEGIN

/**
 * Saves render buffer to .exr image file (OpenGL version).
 */
inline bool saveEXR(const GLfloat* pixels, const std::string& filename, const int width, const int height) {
    EXRHeader header;
    InitEXRHeader(&header);

    EXRImage image;
    InitEXRImage(&image);

    std::vector<float> images[3];
    images[0].resize(width * height);
    images[1].resize(width * height);
    images[2].resize(width * height);

    int cur, pix;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            pix = i * width + j;
            cur = 3 * ((height - i - 1) * width + j);

            v3f rgb = v3f(pixels[cur], pixels[cur + 1], pixels[cur + 2]);
            v3f tmp = convertSRGBToLinear(rgb);

            images[0][pix] = tmp.x;
            images[1][pix] = tmp.y;
            images[2][pix] = tmp.z;
        }
    }

    float* image_ptr[3];
    image_ptr[0] = &(images[2].at(0)); // B
    image_ptr[1] = &(images[1].at(0)); // G
    image_ptr[2] = &(images[0].at(0)); // R

    image.images = (unsigned char**) image_ptr;
    image.width = width;
    image.height = height;

    header.num_channels = 3;
    header.channels = (EXRChannelInfo*) malloc(sizeof(EXRChannelInfo) * header.num_channels);

    // Must be (A)BGR order, since most of EXR viewers expect this channel order
    strncpy(header.channels[0].name, "B", 255);
    header.channels[0].name[strlen("B")] = '\0';
    strncpy(header.channels[1].name, "G", 255);
    header.channels[1].name[strlen("G")] = '\0';
    strncpy(header.channels[2].name, "R", 255);
    header.channels[2].name[strlen("R")] = '\0';

    header.pixel_types = (int*) malloc(sizeof(int) * header.num_channels);
    header.requested_pixel_types = (int*) malloc(sizeof(int) * header.num_channels);
    for (int i = 0; i < header.num_channels; i++) {
        header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT;
        header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF;
    }

    const char* err = nullptr;
    int ret = SaveEXRImageToFile(&image, &header, filename.c_str(), &err);
    if (ret != TINYEXR_SUCCESS) {
        fprintf(stderr, "Save EXR err: %s\n", err);
        FreeEXRErrorMessage(err);
        return false;
    }
    std::cout << "Saved EXR image to " << filename << std::endl;

    free(header.channels);
    free(header.pixel_types);
    free(header.requested_pixel_types);
    return true;
}

/**
 * Saves render buffer to .exr image file.
 */
inline bool saveEXR(const std::unique_ptr<v3f[]>& rgb, const std::string& filename, const int width, const int height) {
    EXRHeader header;
    InitEXRHeader(&header);

    EXRImage image;
    InitEXRImage(&image);

    image.num_channels = 3;

    std::vector<float> images[3];
    images[0].resize(width * height);
    images[1].resize(width * height);
    images[2].resize(width * height);

    // Split RGBRGBRGB... into R, G and B layer
    for (int i = 0; i < width * height; i++) {
        images[0][i] = rgb[i].x;
        images[1][i] = rgb[i].y;
        images[2][i] = rgb[i].z;
    }

    float* image_ptr[3];
    image_ptr[0] = &(images[2].at(0)); // B
    image_ptr[1] = &(images[1].at(0)); // G
    image_ptr[2] = &(images[0].at(0)); // R

    image.images = (unsigned char**) image_ptr;
    image.width = width;
    image.height = height;

    header.num_channels = 3;
    header.channels = (EXRChannelInfo*) malloc(sizeof(EXRChannelInfo) * header.num_channels);

    // Must be (A)BGR order, since most of EXR viewers expect this channel order
    strncpy(header.channels[0].name, "B", 255);
    header.channels[0].name[strlen("B")] = '\0';
    strncpy(header.channels[1].name, "G", 255);
    header.channels[1].name[strlen("G")] = '\0';
    strncpy(header.channels[2].name, "R", 255);
    header.channels[2].name[strlen("R")] = '\0';

    header.pixel_types = (int*) malloc(sizeof(int) * header.num_channels);
    header.requested_pixel_types = (int*) malloc(sizeof(int) * header.num_channels);
    for (int i = 0; i < header.num_channels; i++) {
        header.pixel_types[i] = TINYEXR_PIXELTYPE_FLOAT;
        header.requested_pixel_types[i] = TINYEXR_PIXELTYPE_HALF;
    }

    const char* err = nullptr;
    int ret = SaveEXRImageToFile(&image, &header, filename.c_str(), &err);
    if (ret != TINYEXR_SUCCESS) {
        fprintf(stderr, "Save EXR err: %s\n", err);
        FreeEXRErrorMessage(err);
        return false;
    }
    std::cout << "\nSaved EXR image to " << filename << std::endl;

    free(header.channels);
    free(header.pixel_types);
    free(header.requested_pixel_types);
    return true;
}

/**
 * Variadic template constructor to support printf-style arguments.
 */
class TinyRenderException : public std::runtime_error {
  public:
    template<typename... Args>
    TinyRenderException(const char* fmt, const Args& ... args)
        : std::runtime_error(tfm::format(fmt, args...)) { }
};

TR_NAMESPACE_END