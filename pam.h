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
	f_pixel(double a, double r, double g, double b)
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

	f_pixel& operator *= (double v)
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

	f_pixel& operator /= (double v)
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
		a = fabs(a);
		r = fabs(r);
		g = fabs(g);
		b = fabs(b);
		return *this;
	}

	f_pixel() {}
	double a, r, g, b;
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
		max(l.a, r.a),
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
		&& l.a == r.a
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

	px.r *= px.a;
	px.g *= px.a;
	px.b *= px.a;

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
	if (px.a < 1.0/256.0) {
		rgb_pixel ret(0,0,0,0);
		return ret;
	}

	double r,g,b,a;

	gamma /= internal_gamma;

	// 256, because numbers are in range 1..255.9999… rounded down
	r = pow(px.r/px.a, gamma)*256.0;
	g = pow(px.g/px.a, gamma)*256.0;
	b = pow(px.b/px.a, gamma)*256.0;
	a = px.a*256.0;

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
	xyz.r = 0.412453 * rgb.r + 0.35758 * rgb.g + 0.180423 * rgb.b; // X
	xyz.g = 0.212671 * rgb.r + 0.71516 * rgb.g + 0.072169 * rgb.b; // Y
	xyz.b = 0.019334 * rgb.r + 0.119193 * rgb.g + 0.950227 * rgb.b; // Z
	xyz.a = rgb.a;
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
	if (xyz.r <= S) {
		f.r = (K * xyz.r + 16) / 116;
	}else {
		f.r = pow(xyz.r, 0.3333);
	}
	// fy
	if (xyz.g <= S) {
		f.g = (K * xyz.g + 16) / 116;
	}else {
		f.g = pow(xyz.g, 0.3333);
	}
	// fz
	if (xyz.b <= S) {
		f.b = (K * xyz.b + 16) / 116;
	}else {
		f.b = pow(xyz.b, 0.3333);
	}
	lab.r = 116 * f.g -16;
	lab.g = 500 * (f.r - f.g);
	lab.b = 200 * (f.g - f.b);
	lab.a = xyz.a;

	// normalize
	lab.r /= 100.0;
	lab.g = (lab.g + 134) / 354.0;
	lab.b = (lab.b + 140) / 362.0;

	return lab;
}

static inline
f_pixel rgb2lab(f_pixel rgb)
{
	return xyz2lab(rgb2xyz(rgb));
}

static inline
double colordifference(f_pixel px, f_pixel py)
{
//	f_pixel diff = px - py;
//	f_pixel diff = rgb2xyz(px) - rgb2xyz(py);
	f_pixel diff = rgb2lab(px) - rgb2lab(py);
//	diff.square();
//	diff.g *= 0.7; // change this value to adjust balance between chrominance and luminance..
//	diff.b *= 0.7; // change this value to adjust balance between chrominance and luminance..
	
	diff.square();
	return 
		diff.a * 3.0
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

std::vector<hist_item> pam_computeacolorhist(
	const rgb_pixel*const apixels[], int cols, int rows,
	double gamma, int maxacolors, int ignorebits,
	const double* importance_map
	);


