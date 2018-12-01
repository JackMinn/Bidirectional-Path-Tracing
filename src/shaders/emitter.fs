/*
    This file is part of TinyRender, an educative PBR system.
    Designed for ECSE 446/546 Realistic Image Synthesis, McGill University.

    Copyright (c) 2018 by Derek Nowrouzezahrai and others.
*/
#version 330 core


in vec3 vPos;
in vec3 vNormal;

out vec3 color;

void main()
{
    vec3 a = vec3(vPos); // suppress warning
    vec3 b = vec3(vNormal); // suppress warning
    color = vec3(1.0,1.0,1.0);
}
