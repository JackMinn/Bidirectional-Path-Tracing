/*
    This file is part of TinyRender, an educative PBR system.
    Designed for ECSE 446/546 Realistic Image Synthesis, McGill University.

    Copyright (c) 2018 by Derek Nowrouzezahrai and others.
*/
#version 330 core


in vec3 FragPos;
in vec3 FragNormal;

layout(location = 0) out vec3 gPosition;
layout(location = 1) out vec3 gNormal;

void main()
{
    gPosition = FragPos;
    gNormal = normalize(FragNormal);
}
