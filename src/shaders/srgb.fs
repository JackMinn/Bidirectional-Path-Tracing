/*
    This file is part of TinyRender, an educative PBR system.
    Designed for ECSE 446/546 Realistic Image Synthesis, McGill University.

    Copyright (c) 2018 by Derek Nowrouzezahrai and others.
*/
#version 330 core

uniform sampler2D textureScreen;

in vec2 texCoords;

out vec3 color;

void main()
{
    // Convert colors from Linear to sRGB (gamma correction)
    // YOU DONT HAVE TO UNDERSTAND ANY OF THIS, IGNORE THIS SHADER

    vec3 colorLinear = texture(textureScreen, texCoords).xyz;
    vec3 S1 = sqrt(colorLinear);
    vec3 S2 = sqrt(S1);
    vec3 S3 = sqrt(S2);
    color = 0.662002687f*S1 + 0.684122060f*S2 - 0.323583601f*S3 - 0.0225411470f*colorLinear;
}
