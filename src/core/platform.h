/*
    This file is part of TinyRender, an educative rendering system.

    Designed for ECSE 446/546 Realistic/Advanced Image Synthesis.
    Derek Nowrouzezahrai, McGill University.
*/

#pragma once

#include <iostream>
#include <algorithm>
#include <random>
#include <memory>

#pragma warning(push, 0)
//#define GLM_FORCE_INLINE
//#define GLM_FORCE_AVX2
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/color_space.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/ext.hpp>
#pragma warning(pop)

using namespace std;
#if defined(_WIN32)
#include <experimental/filesystem>
namespace fs = experimental::filesystem;
using I = int;
// Reverses byte order
inline I bswap(I x) { return _byteswap_ulong(x); }
#else
#if defined(__APPLE__)
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#elif defined(__GNUC__)
#include <experimental/filesystem>
namespace fs = experimental::filesystem;
#endif
inline int bswap(int x) { return __builtin_bswap32(x); }
inline string pp(string p) {
    replace(p.begin(), p.end(), '\\', '/');
    return p;
}
#endif

#ifdef M_PI
#undef M_PI
#endif
#define M_PI       3.14159265358979323846f
#define INV_PI     0.31830988618379067154f
#define INV_TWOPI  0.15915494309189533577f
#define INV_FOURPI 0.07957747154594766788f

#define deg2rad M_PI / 180.f
#define Epsilon 1e-8f
typedef glm::fvec2 v2f;
typedef glm::fvec3 v3f;
typedef glm::fvec4 v4f;
typedef glm::fvec2 p2f;
typedef glm::fvec3 p3f;
typedef glm::mat4 mat4f;

#define TR_NAMESPACE_BEGIN namespace TinyRender {
#define TR_NAMESPACE_END }

#ifdef _WIN32
#define sincosf(x,s,c) (*s = sin(x), *c = cos(x))
#else
#if defined(__APPLE__)
#define sincos(x, s, c) __sincos(x, s, c)
#define sincosf(x, s, c) __sincosf(x, s, c)
#endif
#if defined(__FreeBSD__)
#define sincos(x,s,c) (*s = sin(x), *c = cos(x))
#endif
#endif


#define NOTODO