/*
    This file is part of TinyRender, an educative PBR system.
    Designed for ECSE 446/546 Realistic Image Synthesis, McGill University.

    Copyright (c) 2018 by Derek Nowrouzezahrai and others.
*/
#version 330 core


in vec3 vColor;
out vec4 color;

void main()
{
    color = vec4(vColor,1.0);
}
