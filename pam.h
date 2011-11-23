#pragma once

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

#include <vector>
#include <math.h>
#include <assert.h>

template <typename T>
T min(T a, T b)
{
	if (a < b) {
		return a;
	}else {
		return b;
	}
}

template <typename T>
T min(T a, T b, T c)
{
	return min(min(a,b), c);
}

template <typename T>
T min(T a, T b, T c, T d)
{
	return min(min(a,b), min(c, d));
}

template <typename T>
T min(T a, T b, T c, T d, T e)
{
	return min(min(a,b,c,d), e);
}

template <typename T>
T max(T a, T b)
{
	if (a < b) {
		return b;
	}else {
		return a;
	}
}

template <typename T>
T max(T a, T b, T c)
{
	return max(max(a,b), c);
}

template <typename T>
T max(T a, T b, T c, T d)
{
	return max(max(a, b), max(c, d));
}

template <typename T>
T max(T a, T b, T c, T d, T e)
{
	return max(max(a, b, c, d), e);
}

template <typename T>
T limitValue(T val, T min, T max)
{
	if (val < min) return min;
	if (max < val) return max;
	return val;
}

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

	rgb_pixel& operator *= (int v)
	{
		r *= v;
		g *= v;
		b *= v;
		a *= v;
		return *this;
	}

	rgb_pixel& operator /= (int v)
	{
		r /= v;
		g /= v;
		b /= v;
		a /= v;
		return *this;
	}
};

static inline
rgb_pixel operator * (const rgb_pixel& l, int v)
{
	rgb_pixel ret = l;
	ret *= v;
	return ret;
}

static inline
rgb_pixel operator / (const rgb_pixel& l, int v)
{
	rgb_pixel ret = l;
	ret /= v;
	return ret;
}



struct f_pixel {
	f_pixel(double alpha, double r, double g, double b)
		:
		alpha(alpha),
		r(r),
		g(g),
		b(b)
	{
	}
	
	f_pixel& operator += (const f_pixel& p)
	{
		alpha += p.alpha;
		r += p.r;
		g += p.g;
		b += p.b;
		return *this;
	}

	f_pixel& operator - ()
	{
		alpha = -alpha;
		r = -r;
		g = -g;
		b = -b;
		return *this;
	}

	f_pixel& operator -= (const f_pixel& p)
	{
		alpha -= p.alpha;
		r -= p.r;
		g -= p.g;
		b -= p.b;
		return *this;
	}

	f_pixel& operator *= (double v)
	{
		alpha *= v;
		r *= v;
		g *= v;
		b *= v;
		return *this;
	}

	f_pixel& operator *= (const f_pixel& p)
	{
		alpha *= p.alpha;
		r *= p.r;
		g *= p.g;
		b *= p.b;
		return *this;
	}

	f_pixel& operator /= (double v)
	{
		alpha /= v;
		r /= v;
		g /= v;
		b /= v;
		return *this;
	}

	f_pixel& square()
	{
		alpha *= alpha;
		r *= r;
		g *= g;
		b *= b;
		return *this;
	}

	f_pixel& abs()
	{
		alpha = fabs(alpha);
		r = fabs(r);
		g = fabs(g);
		b = fabs(b);
		return *this;
	}

	f_pixel() {}
	double alpha;
	
	union {
		struct { double r, g, b; };
		struct { double x, y, z; };
		struct { double l, a, b; };
		double arr[3];
	};
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
		max(l.alpha, r.alpha),
		max(l.r, r.r),
		max(l.g, r.g),
		max(l.b, r.b)
		);
}

static inline
f_pixel operator * (const f_pixel& l, double v)
{
	f_pixel ret = l;
	ret *= v;
	return ret;
}

static inline
f_pixel operator / (const f_pixel& l, double v)
{
	f_pixel ret = l;
	ret /= v;
	return ret;
}

static inline
f_pixel operator / (const rgb_pixel& l, double v)
{
	f_pixel ret(l.a, l.r, l.g, l.b);
	ret /= v;
	return ret;
}

static inline
bool operator == (const f_pixel& l, const f_pixel& r)
{
	return 1
		&& l.alpha == r.alpha
		&& l.r == r.r
		&& l.g == r.g
		&& l.b == r.b
		;
}

static const double internal_gamma = 0.45455;

/**
 Converts scalar color to internal gamma and premultiplied alpha.
 (premultiplied color space is much better for blending of semitransparent colors)
 */
static inline
f_pixel to_f_scalar(double gamma, f_pixel px)
{
	if (gamma != internal_gamma) {
		px.r = pow(px.r, internal_gamma/gamma);
		px.g = pow(px.g, internal_gamma/gamma);
		px.b = pow(px.b, internal_gamma/gamma);
	}

	px.r *= px.alpha;
	px.g *= px.alpha;
	px.b *= px.alpha;

	return px;
}

/**
  Converts 8-bit RGB with given gamma to scalar RGB
 */
static inline
f_pixel to_f(double gamma, rgb_pixel px)
{
	return to_f_scalar(
		gamma,
		px / 255.0
	);
}

static inline
rgb_pixel to_rgb(double gamma, f_pixel px)
{
	if (px.alpha < 1.0/256.0) {
		rgb_pixel ret(0,0,0,0);
		return ret;
	}

	double r,g,b,a;

	gamma /= internal_gamma;

	// 256, because numbers are in range 1..255.9999… rounded down
	r = pow(px.r/px.alpha, gamma)*256.0;
	g = pow(px.g/px.alpha, gamma)*256.0;
	b = pow(px.b/px.alpha, gamma)*256.0;
	a = px.alpha*256.0;

	rgb_pixel ret;
	ret.r = r>=255 ? 255 : (r<=0 ? 0 : r);
	ret.g = g>=255 ? 255 : (g<=0 ? 0 : g);
	ret.b = b>=255 ? 255 : (b<=0 ? 0 : b);
	ret.a = a>=255 ? 255 : a;
	return ret;
}

static inline
f_pixel rgb2xyz(f_pixel rgb)
{
	f_pixel xyz;
	xyz.x = 0.4124564 * rgb.r + 0.3575761 * rgb.g + 0.1804375 * rgb.b; // X
	xyz.y = 0.2126729 * rgb.r + 0.7151522 * rgb.g + 0.0721750 * rgb.b; // Y
	xyz.z = 0.0193339 * rgb.r + 0.1191920 * rgb.g + 0.9503041 * rgb.b; // Z
	xyz.alpha = rgb.alpha;
	return xyz;
}

static inline
f_pixel xyz2rgb(f_pixel xyz)
{
	// http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_RGB.html
	f_pixel rgb;
	rgb.r = 3.2404542 * xyz.x -1.5371385 * xyz.y -0.4985314 * xyz.z; // X
	rgb.g = -0.9692660 * xyz.x + 1.8760108 * xyz.y + 0.0415560 * xyz.z; // Y
	rgb.b = 0.0556434 * xyz.x -0.2040259 * xyz.y + 1.0572252 * xyz.z; // Z
	rgb.alpha = xyz.alpha;
	return rgb;
}

static inline
f_pixel lab2xyz(f_pixel lab)
{
	lab.l *= 100.0;
	lab.a *= 354.0;
	lab.b *= 262.0;
	lab.a -= 134.0;
	lab.b -= 140.0;

	// http://www.brucelindbloom.com/index.html?Eqn_Lab_to_XYZ.html
	static const double E = 216.0 / 24389.0;
	static const double K = 24389.0 / 27.0;
	static const double KE = K * E;
	f_pixel f;
	f.y = (lab.l + 16.0) / 116.0;
	f.x = lab.a / 500.0 + f.y;
	f.z = f.y - lab.b / 200.0;
	
	f_pixel f3 = f;
	f3 *= f;
	f3 *= f;
	f_pixel xyz;
	if (f3.x > E) {
		xyz.x = f3.x;
	}else {
		xyz.x = (116.0 * f.x - 16.0) / K;
	}
	if (lab.l > KE) {
		double t = (lab.l + 16.0) / 116.0;
		xyz.y = t * t * t;
	}else {
		xyz.y = lab.l / K;
	}
	if (f3.z > E) {
		xyz.z = f3.z;
	}else {
		xyz.z = (116.0 * f.z - 16.0) / K;
	}
	xyz.alpha = lab.alpha;
	return xyz;
}

static inline
f_pixel xyz2lab(f_pixel xyz)
{
	// http://www.brucelindbloom.com/index.html?Eqn_XYZ_to_Lab.html
	static const double K = 24389.0 / 27.0;
	static const double S = 216.0 / 24389.0;
	f_pixel f;
	f_pixel lab;
	// fx
	if (xyz.x <= S) {
		f.x = (K * xyz.x + 16) / 116;
	}else {
		f.x = pow(xyz.x, 0.3333);
	}
	// fy
	if (xyz.y <= S) {
		f.y = (K * xyz.y + 16) / 116;
	}else {
		f.y = pow(xyz.y, 0.3333);
	}
	// fz
	if (xyz.z <= S) {
		f.z = (K * xyz.z + 16) / 116;
	}else {
		f.z = pow(xyz.z, 0.3333);
	}
	lab.l = 116 * f.y -16;
	lab.a = 500 * (f.x - f.y);
	lab.b = 200 * (f.y - f.z);
	lab.alpha = xyz.alpha;

	// normalize
	lab.l /= 100.0;
	lab.a = (lab.a + 134) / 354.0;
	lab.b = (lab.b + 140) / 262.0;

	return lab;
}

static inline
f_pixel rgb2lab(f_pixel rgb)
{
	return xyz2lab(rgb2xyz(rgb));
}

static inline
f_pixel lab2rgb(f_pixel lab)
{
	return xyz2rgb(lab2xyz(lab));
}

static inline
double colordifference(f_pixel px, f_pixel py)
{
	f_pixel diff = px - py;
	
	diff.square();
	return 
		diff.alpha * 3.0
		+ diff.r
		+ diff.g
		+ diff.b
		;
}

/* from pamcmap.h */

struct hist_item {
	f_pixel acolor;
	double adjusted_weight;
	double perceptual_weight;
};

struct colormap_item {
	f_pixel acolor;
	double popularity;
};

struct acolorhist_list_item {
	f_pixel acolor;
	acolorhist_list_item* next;
	double perceptual_weight;
};

struct acolorhash_table {
	struct mempool* mempool;
	acolorhist_list_item** buckets;
};

int best_color_index(
	const std::vector<colormap_item>& map,
	f_pixel px,
	double* dist_out
	);

std::vector<hist_item> pam_computeacolorhist(
	const f_pixel* input, size_t width, size_t height,
	int maxacolors, int ignorebits,
	const double* importance_map
	);

