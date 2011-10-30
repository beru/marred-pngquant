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
#ifndef MAX
#  define MAX(a,b)	((a) > (b)? (a) : (b))
#  define MIN(a,b)	((a) < (b)? (a) : (b))
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
		MAX(l.a, r.a),
		MAX(l.r, r.r),
		MAX(l.g, r.g),
		MAX(l.b, r.b)
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
double colordifference(f_pixel px, f_pixel py)
{
	f_pixel diff = px - py;
	diff.square();
	return diff.a * 3.0 + diff.r + diff.g + diff.b;
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


