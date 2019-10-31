#version 450 core

layout(local_size_x = 1, local_size_y = 1) in;
layout(rgba32f, binding = 0) uniform image2D img_output;

#define pi 3.141

void main()
{
	vec2 dims = vec2(imageSize(img_output));
  vec2 pixel_coord = gl_GlobalInvocationID.xy;

	const int sx = 3;
	const int sy = 3;
	const int sz = sx * sy;
	float r[sz];
	float g[sz];
	float b[sz];
	for (int i = 0; i < sz; ++i) {
		vec3 rgb = imageLoad(img_output, ivec2(pixel_coord + vec2(i % 3 - sx / 2, i / 3 - sy / 2))).rgb;
		r[i] = rgb.r;
		g[i] = rgb.g;
		b[i] = rgb.b;
	}

	for (int i = 0; i < sz; ++i) {
		for (int j = 0; j < sz; ++j) {
			if (r[i] < r[j]) {
				float tmp = r[i];
				r[i] = r[j];
				r[j] = tmp;
			}
			if (g[i] < g[j]) {
				float tmp = g[i];
				g[i] = g[j];
				g[j] = tmp;
			}
			if (b[i] < b[j]) {
				float tmp = b[i];
				b[i] = b[j];
				b[j] = tmp;
			}
		}
	}

	imageStore(img_output, ivec2(pixel_coord), vec4(r[sz / 2], g[sz / 2], b[sz / 2], 1));
}
