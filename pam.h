/**
 ** Copyright (C) 1989, 1991 by Jef Poskanzer.
 ** Copyright (C) 1997, 2000, 2002 by Greg Roelofs; based on an idea by
 **								   Stefan Schneider.
 ** (C) 2011 by Kornel Lesinski.
 **
 ** Permission to use, copy, modify, and distribute this software and its
 ** documentation for any purpose and without fee is hereby granted, provided
 ** that the above copyright notice appear in all copies and that both that
 ** copyright notice and this permission notice appear in supporting
 ** documentation.	This software is provided "as is" without express or
 ** implied warranty.
 */

#include <math.h>
#ifndef MAX
#  define MAX(a,b)	((a) > (b)? (a) : (b))
#  define MIN(a,b)	((a) < (b)? (a) : (b))
#endif

#ifdef __SSE3__
#define USE_SSE
#endif

#ifdef USE_SSE
#include <pmmintrin.h>
#endif


/* from pam.h */

struct rgb_pixel {
	rgb_pixel(unsigned char r, unsigned char g, unsigned char b, unsigned char a)
		:
		r(r),
		g(g),
		b(b),
		a(a)
	{
	}
	
	rgb_pixel() {}

	unsigned char r, g, b, a;
};

struct f_pixel {
	f_pixel(float a, float r, float g, float b)
		:
		a(a),
		r(r),
		g(g),
		b(b)
	{
	}
	
	f_pixel& operator += (const f_pixel& p)
	{
		a += p.a;
		r += p.r;
		g += p.g;
		b += p.b;
		return *this;
	}

	f_pixel& operator - ()
	{
		a = -a;
		r = -r;
		g = -g;
		b = -b;
		return *this;
	}

	f_pixel& operator -= (const f_pixel& p)
	{
		a -= p.a;
		r -= p.r;
		g -= p.g;
		b -= p.b;
		return *this;
	}

	f_pixel& operator *= (float v)
	{
		a *= v;
		r *= v;
		g *= v;
		b *= v;
		return *this;
	}

	f_pixel& operator *= (const f_pixel& p)
	{
		a *= p.a;
		r *= p.r;
		g *= p.g;
		b *= p.b;
		return *this;
	}

	f_pixel& operator /= (float v)
	{
		a /= v;
		r /= v;
		g /= v;
		b /= v;
		return *this;
	}

	f_pixel& square()
	{
		a *= a;
		r *= r;
		g *= g;
		b *= b;
		return *this;
	}

	f_pixel& abs()
	{
		a = fabsf(a);
		r = fabsf(r);
		g = fabsf(g);
		b = fabsf(b);
		return *this;
	}

	f_pixel() {}
	float a, r, g, b;
};

static inline
f_pixel operator + (const f_pixel& l, const f_pixel& r)
{
	f_pixel ret = l;
	ret += r;
	return ret;
}

static inline
f_pixel operator - (const f_pixel& l, const f_pixel& r)
{
	f_pixel ret = l;
	ret -= r;
	return ret;
}

static inline
f_pixel max(const f_pixel& l, const f_pixel& r)
{
	return f_pixel(
		MAX(l.a, r.a),
		MAX(l.r, r.r),
		MAX(l.g, r.g),
		MAX(l.b, r.b)
		);
}

static inline
f_pixel operator * (const f_pixel& l, float v)
{
	f_pixel ret = l;
	ret *= v;
	return ret;
}

static inline
f_pixel operator / (const f_pixel& l, float v)
{
	f_pixel ret = l;
	ret /= v;
	return ret;
}

static inline
bool operator == (const f_pixel& l, const f_pixel& r)
{
	return 1
		&& l.a == r.a
		&& l.r == r.r
		&& l.g == r.g
		&& l.b == r.b
		;
}

static const float internal_gamma = 0.45455;

/**
 Converts scalar color to internal gamma and premultiplied alpha.
 (premultiplied color space is much better for blending of semitransparent colors)
 */
static inline
f_pixel to_f_scalar(float gamma, f_pixel px)
{
	if (gamma != internal_gamma) {
		px.r = powf(px.r, internal_gamma/gamma);
		px.g = powf(px.g, internal_gamma/gamma);
		px.b = powf(px.b, internal_gamma/gamma);
	}

	px.r *= px.a;
	px.g *= px.a;
	px.b *= px.a;

	return px;
}

/**
  Converts 8-bit RGB with given gamma to scalar RGB
 */
static inline
f_pixel to_f(float gamma, rgb_pixel px)
{
	return to_f_scalar(
		gamma,
		f_pixel(
			px.a/255.0f,
			px.r/255.0f,
			px.g/255.0f,
			px.b/255.0f
		)
	);
}

static inline
rgb_pixel to_rgb(float gamma, f_pixel px)
{
	if (px.a < 1.0/256.0) {
		rgb_pixel ret(0,0,0,0);
		return ret;
	}

	float r,g,b,a;

	gamma /= internal_gamma;

	// 256, because numbers are in range 1..255.9999… rounded down
	r = powf(px.r/px.a, gamma)*256.0f;
	g = powf(px.g/px.a, gamma)*256.0f;
	b = powf(px.b/px.a, gamma)*256.0f;
	a = px.a*256.0;

	rgb_pixel ret;
	ret.r = r>=255 ? 255 : (r<=0 ? 0 : r);
	ret.g = g>=255 ? 255 : (g<=0 ? 0 : g);
	ret.b = b>=255 ? 255 : (b<=0 ? 0 : b);
	ret.a = a>=255 ? 255 : a;
	return ret;
}


static inline
float colordifference_stdc(f_pixel px, f_pixel py)
{
	return (px.a - py.a) * (px.a - py.a) * 3.0 +
			(px.r - py.r) * (px.r - py.r) +
			(px.g - py.g) * (px.g - py.g) +
			(px.b - py.b) * (px.b - py.b);
}

static inline
float colordifference(f_pixel px, f_pixel py)
{
#ifdef USE_SSE
	__m128 vpx = _mm_load_ps((const float*)&px);
	__m128 vpy = _mm_load_ps((const float*)&py);

	__m128 tmp = _mm_sub_ps(vpx, vpy); // t = px - py
	tmp = _mm_mul_ps(tmp, tmp); // t = t * t
	tmp = _mm_mul_ss(tmp, _mm_set_ss(3.0)); // alpha * 3.0

	tmp = _mm_hadd_ps(tmp,tmp); // 0+1 2+3 0+1 2+3
	__m128 rev = _mm_shuffle_ps(tmp, tmp, 0x1B); // reverses vector 2+3 0+1 2+3 0+1
	tmp = _mm_add_ss(tmp, rev); // 0+1 + 2+3

	float res = _mm_cvtss_f32(tmp);
	assert(fabs(res - colordifference_stdc(px,py)) < 0.001);
	return res;
#else
	return colordifference_stdc(px,py);
#endif
}

/* from pamcmap.h */

struct hist_item {
	f_pixel acolor;
	float adjusted_weight;
	float perceptual_weight;
};

struct hist {
	hist_item* achv;
	int size;
};

struct colormap_item {
	f_pixel acolor;
	float popularity;
};

struct colormap {
	colormap_item* palette;
	int colors;
};

struct acolorhist_list_item {
	f_pixel acolor;
	acolorhist_list_item* next;
	float perceptual_weight;
};

struct acolorhash_table {
	struct mempool* mempool;
	acolorhist_list_item** buckets;
};


hist* pam_computeacolorhist(const rgb_pixel*const apixels[], int cols, int rows, double gamma, int maxacolors, int ignorebits, int use_contrast);
void pam_freeacolorhist(hist* h);

colormap* pam_colormap(int colors);
void pam_freecolormap(colormap* c);

