/*
    This file is part of TinyRender, an educative PBR system.
    Designed for ECSE 446/546 Realistic Image Synthesis, McGill University.

    Copyright (c) 2018 by Derek Nowrouzezahrai and others.
*/
#version 330 core


layout(location = 0) in vec3 position;
layout(location = 1) in vec3 normal;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec3 FragPos;
out vec3 FragNormal;

void main()
{
    vec4 viewPos = view * model * vec4(position, 1.0);
    FragPos = viewPos.xyz;

    mat4 normalMat = transpose(inverse(view * model));
    FragNormal = (normalMat * vec4(normal, 0.0)).xyz;

    gl_Position = projection * viewPos;
}
