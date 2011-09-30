// Packing a [0-1] float value into a 4D vector
// where each component will be a 8-bits integer:

vec4 packFloatToVec4i(const float value)
{
	const vec4 bitSh = vec4((1<<24), (1<<16), (1<<8), 1.0);
	const vec4 bitMsk = vec4(0.0, 1.0/256.0, 1.0/256.0, 1.0/256.0);
	vec4 res = fract(value * bitSh);
	res -= res.xxyz * bitMsk;
	return res;
}

// Unpacking a [0-1] float value from a 4D vector
// where each component was a 8-bits integer:
float unpackFloatFromVec4i(const vec4 value)
{
	const vec4 bitSh = vec4(1.0/(1<<24), 1.0/(1<<16), 1.0/(1<<8), 1.0);
	return(dot(value, bitSh));
}

//=======================================================================
uniform sampler2D texture;
varying vec4 texc;
uniform float width;
uniform float height;

// float z in [0,1]
float get_float32(sampler2D texture, vec2 v)
{
	vec4 z = texture2D(texture, v.xy);
	return z.x/(1<<8) + z.y/(1<<16) + z.z/(1<<24) + z.w/(1<<32);
}

vec4 set_float32(float f)
{
	return vec4(f*(1<<8) + f*(1<<16) + f*(1<<24) + f*(1<<32));
}

// texc in [-1, 1]
// XXX maybe
void main(void)
{
	vec2 dx = vec2(2.0/width, 0.0);
	vec2 dy = vec2(0.0, 2.0/height);
	vec2 z = texc.st;

	float resx = get_float32(texture, texc.st - dx) +
	             get_float32(texture, texc.st)*2 +
	             get_float32(texture, texc.st + dx);
	
	float resy = get_float32(texture, texc.st - dy) +
	             get_float32(texture, texc.st)*2 +
	             get_float32(texture, texc.st + dy);
	
	float res = get_float32(texture, texc.st) +
	            a*a*t/dx/dx * resx +
	            b*b*t/dy/dy * resy;
	
	gl_FragColor = set_float32(res);
}

