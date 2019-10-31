#version 450 core

layout(local_size_x = 1, local_size_y = 1) in;
layout(rgba32f, binding = 0) uniform image2D img_output;

uniform float elapsed_time;
uniform float delta_time;

#define pi 3.141

float n31(vec3 n)
{
	return fract(sin(dot(n, vec3(12.9898,78.233,12.5429))) * 43758.5453);
}

float n11(float n)
{
	return n31(vec3(n + 259.753, n * 87.231, (n - 2.97) * 754.234));
}

vec3 n33(vec3 n)
{
		float r0 = n31(n);
		float r1 = n31(n + r0);
		float r2 = n31(n + r1);
		return vec3(sin(r0), sin(r1), sin(r2));
}

float sphere(vec3 ro, vec3 rd, vec3 c, float r)
{   
    vec3 a = rd * (ro - c);
    float b = dot(rd, rd);
    float p = dot(a, vec3(2)) / b;
    a = ro - c;
    float q = (dot(a, ro - c) - r * r) / b;
    p = p / 2.;
    q = p * p - q;
    
    if (q < 0.)
        return 0.;
    
    q = sqrt(q);
    float l1 = -p + q;
    float l2 = -p - q;
    
    if (l1 > 0. && l2 > 0.)
        return min(l1, l2);
    else
        return max(l1, l2);
}

float plane(vec3 ro, vec3 rd, vec3 o, vec3 n)
{
	/* dot(o + l * v, n) = 0. */
	/* n.x * (o.x + x) +	n.y * (o.y + y) + n.z * (o.z + z) = 0 */
	/* x = ro.x + rd.x * l; */
	/* y = ro.y + rd.y * l; */
	/* z = ro.z + rd.z * l; */
  /*  */
	/* n.x * (ro.x + rd.x * l - o.x) */
	/* n.x * ro.x + n.x * rd.x * l - n.x * o.x */

	/* dot(n, v - o) = 0 */
	/* n.x * (x - o.x) + n.y * (y - o.y) + n.z * (z - o.z) */
	/* n.x * (l * rd.x + ro.x - o.x) */
	return - (n.x * ro.x - n.x * o.x + n.y * ro.y - n.y * o.y + n.z * ro.z - n.z * o.z) / (n.x * rd.x + n.y * rd.y + n.z * rd.z);
	/* return (dot(n, ro) - dot(n, o)) / -dot(n, rd); */

	/* float l = (n.x * ro.x - n.x * o.x + n.y * ro.y - n.y * o.y + n.z * ro.z - n.z * o.z) / (n.x * rd.x + n.y * rd.y + n.z * rd.z); */
	/* return l; */
}

struct material_t {
	bool light;
	bool metal;
	float roughness;
	vec3 color;
};

bool scene(vec3 ro, vec3 rd, out vec3 p, out vec3 n, out material_t m)
{
	int si = -1;
	float lmin;

	const int sz = 9;
	vec3 sc[sz] = vec3[sz] (
			vec3(0, 0, 5),
			vec3(-1, 0, 5),
			vec3(1, 0, 5),
			vec3(1, 1, 5),
			vec3(1, -1, 5),
			vec3(0, -1, 5),
			vec3(0, 1, 5),
			vec3(-1, 1, 5),
			vec3(-1, -1, 5)
	);

	for (int i = 0; i < sz; ++i) {
		float l = sphere(ro, rd, sc[i], .5);
		if (l > 0. && (si == -1 || l < lmin)) {
			lmin = l;
			si = i;
		}
	}

	const int pz = 6;
	vec3 pc[pz] = vec3[pz] (
		vec3(0, 0, 6.5),
		vec3(0, 0, -1), 
		vec3(-1.5, 0, 0),
		vec3(1.5, 0, 0),
		vec3(0, -1.5, 0),
		vec3(0, 1.5, 0)
	);
	vec3 pn[pz] = vec3[pz] (
		vec3(0, 0, -1),
		vec3(0, 0, 1),
		vec3(1, 0, 0),
		vec3(-1, 0, 0),
		vec3(0, 1, 0),
		vec3(0, -1, 0)
	);
	bool b = false;
	for (int i = 0; i < pz; ++i) {
		float l = plane(ro, rd, pc[i], pn[i]);
		if (l > 0. && (si == -1 || l < lmin)) {
			lmin = l;
			si = sz + i;
			b = true;
		}
	}

	if (si == -1)
		return false;

	p = ro + rd * lmin;
	m.color = mix(n33(vec3(si + 1)), vec3(1), .5);
	m.metal = false;
	m.roughness = 0.;
	m.light = n11(floor(elapsed_time * .2) + si + 2) < .1;

	if (si < sz)
		n = normalize(p - sc[si]);
	else
		n = pn[si - sz];

	return true;
}

vec3 rotate(vec3 v, float r, vec3 n)
{
	vec3 v0 = ((v * n) / (v * n)) * n;
	vec3 v1 = v - v0;
	vec3 w = cross(n, v1);
	float x0 = cos(r) / length(v1);
	float x1 = sin(r) / length(w);
	return length(v1) * (x0 * v1 + x1 * w);
}

void main() {
	vec2 dims = vec2(imageSize(img_output));
  vec2 pixel_coord = gl_GlobalInvocationID.xy;

	float t = elapsed_time * pi * 2.;
	vec2 uv = (pixel_coord - dims * .5) / dims.y;
	vec3 c = vec3(0);

	vec3 ro = vec3(0);
	vec3 rd = normalize(vec3(uv, 1));

	vec3 b = vec3(cos(t) * 3., 3, sin(t) * 1.83);

	bool hl = false;
	vec3 rc = vec3(1, 1, 1);
	float rl = 0.;
	vec3 p = ro;
	vec3 n;
	material_t m;

	for (int i = 0; i < 10; ++i) {
		if (scene(ro, rd, p, n, m)) {
			if (m.light) {
				hl = true;
				rc *= m.color;
				rc /= pow((rl * .05 + 1.), 2.);
				break;
			}

			if (!m.metal)
				rc *= m.color;

			rl += length(p - ro);
			ro = p;

			vec3 rr = reflect(rd, n);
			vec3 rn = normalize(n + normalize(n33(p + mod(elapsed_time, 1.))));

			rd = normalize(mix(rr, rn, .3));
			ro += rd * .001;

		} else {
			rc = vec3(0);
			break;
		}
	}

	/* rc += n; */

	if (!hl)
		rc = vec3(0.);
	c += mix(rc, imageLoad(img_output, ivec2(pixel_coord)).rgb, max(0., 1. - delta_time));
	/* c = rc; */
	imageStore(img_output, ivec2(pixel_coord), vec4(c, 1));
}
