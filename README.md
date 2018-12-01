# Bidirectional-Path-Tracing

This project extends a rendering framework provided by McGills ECSE 446/546 course to provide Bidirectional Path Tracing while also 
supporting light transport for paths with delta distributions like specular reflection and refraction. The implementation of BDPT
leverages Multiple Importance Sampling as described in Eric Veach's thesis to use multiple path sampling strategies and combine them
to reduce variance. Paths for either the light or eye subpaths are terminated probabalistically based on Russian Roulette where the threshold probability is based on the throughput of the path. 

The reduction in variance allows BDPT to converge to visually acceptable results in a much shorter time than 
traditional path or light tracing. Below are 2 sample images that were rendered that would otherwise contain high frequency noise in
the form of fireflies if done with traditional path tracing.



