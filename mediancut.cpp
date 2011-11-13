/*
 **
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

#include <stdlib.h>
#include <assert.h>
#include <stddef.h>

#include "png.h"	/* libpng header; includes zlib.h */
#include "rwpng.h"	/* typedefs, common macros, public prototypes */
#include "pam.h"
#include "mediancut.h"

#include <vector>
#include <algorithm>

#define index_of_channel(ch) (offsetof(f_pixel,ch)/sizeof(double))

static
f_pixel averagepixels(int indx, int clrs, const hist_item achv[], double min_opaque_val)
{
	f_pixel csum(0,0,0,0);
	double sum = 0;
	double maxa = 0;
	int i;

	for (i=0; i<clrs; ++i) {
		double weight = 1.0;
		const hist_item& hist = achv[indx + i];
		f_pixel px = hist.acolor;
		/* give more weight to colors that are further away from average
		 this is intended to prevent desaturation of images and fading of whites
		 */
		f_pixel tmp = f_pixel(0.5f, 0.5f, 0.5f, 0.5f) - px;
		tmp.square();
		weight += tmp.r + tmp.g + tmp.b;
		weight *= hist.adjusted_weight;
		csum += px * weight;
		sum += weight;

		/* find if there are opaque colors, in case we're supposed to preserve opacity exactly (ie_bug) */
		maxa = max(maxa, px.a);
	}

	/* Colors are in premultiplied alpha colorspace, so they'll blend OK
	 even if different opacities were mixed together */
	if (!sum) sum=1;
	csum /= sum;
	
	/** if there was at least one completely opaque color, "round" final color to opaque */
	if (csum.a >= min_opaque_val && maxa >= (255.0/256.0)) csum.a = 1;

	return csum;
}

struct box {
	double variance;
	int sum;
	int ind;
	int colors;
};

struct channelvariance {
	channelvariance(int ch, double var)
		:
		chan(ch),
		variance(var)
	{
	}

	channelvariance() {}
	int chan;
	double variance;
};

static inline
bool operator < (const channelvariance& ch1, const channelvariance& ch2)
{
	return ch1.variance > ch2.variance;
}

static channelvariance channel_sort_order[4];

static inline
bool weightedcompare_other(const f_pixel& p1, const f_pixel& p2)
{
	const double* c1p = (const double*) &p1;
	const double* c2p = (const double*) &p2;
	int chan;
	
	// other channels are sorted backwards
	chan = channel_sort_order[1].chan;
	if (c1p[chan] > c2p[chan]) return true;
	if (c1p[chan] < c2p[chan]) return false;
	
	chan = channel_sort_order[2].chan;
	if (c1p[chan] > c2p[chan]) return true;
	if (c1p[chan] < c2p[chan]) return false;
	
	chan = channel_sort_order[3].chan;
	if (c1p[chan] > c2p[chan]) return true;
	if (c1p[chan] < c2p[chan]) return false;
	
	return true;
}

/** these are specialised functions to make first comparison faster without lookup in channel_sort_order[] */
static
bool weightedcompare_r(const hist_item& lhs, const hist_item& rhs)
{
	const f_pixel& p1 = lhs.acolor;
	const f_pixel& p2 = rhs.acolor;
	double c1 = p1.r;
	double c2 = p2.r;
	if (c1 > c2) return false;
	if (c1 < c2) return true;
	
	return weightedcompare_other(p1, p2);
}

static
bool weightedcompare_g(const hist_item& lhs, const hist_item& rhs)
{
	const f_pixel& p1 = lhs.acolor;
	const f_pixel& p2 = rhs.acolor;
	double c1 = p1.g;
	double c2 = p2.g;
	if (c1 > c2) return false;
	if (c1 < c2) return true;
	
	return weightedcompare_other(p1, p2);
}

static
bool weightedcompare_b(const hist_item& lhs, const hist_item& rhs)
{
	const f_pixel& p1 = lhs.acolor;
	const f_pixel& p2 = rhs.acolor;
	double c1 = p1.b;
	double c2 = p2.b;
	if (c1 > c2) return false;
	if (c1 < c2) return true;
	
	return weightedcompare_other(p1, p2);
}

static
bool weightedcompare_a(const hist_item& lhs, const hist_item& rhs)
{
	const f_pixel& p1 = lhs.acolor;
	const f_pixel& p2 = rhs.acolor;
	double c1 = p1.a;
	double c2 = p2.a;
	if (c1 > c2) return true;
	if (c1 < c2) return false;
	
	return weightedcompare_other(p1, p2);
}

f_pixel channel_variance(const hist_item* achv, int indx, int clrs, double min_opaque_val)
{
	f_pixel mean = averagepixels(indx, clrs, achv, min_opaque_val);
	f_pixel variance(0,0,0,0);

	for (int i=0; i<clrs; ++i) {
		variance += (mean - achv[indx + i].acolor).square();
	}
	return variance;
}

static
void sort_colors_by_variance(f_pixel variance, hist_item* achv, int indx, int clrs)
{
	/*
	 ** Sort dimensions by their variance, and then sort colors first by dimension with highest variance
	 */

	channel_sort_order[0] = channelvariance(index_of_channel(r), variance.r);
	channel_sort_order[1] = channelvariance(index_of_channel(g), variance.g);
	channel_sort_order[2] = channelvariance(index_of_channel(b), variance.b);
	channel_sort_order[3] = channelvariance(index_of_channel(a), variance.a);

	std::sort(channel_sort_order, channel_sort_order+4);
	
	bool (*comp)(const hist_item&, const hist_item&); // comp variable that is a pointer to a function
	const int ch = channel_sort_order[0].chan;
		 if (ch == index_of_channel(r)) comp = weightedcompare_r;
	else if (ch == index_of_channel(g)) comp = weightedcompare_g;
	else if (ch == index_of_channel(b)) comp = weightedcompare_b;
	else comp = weightedcompare_a;

	std::sort(achv+indx, achv+clrs, comp);
}


/*
 ** Find the best splittable box. -1 if no boxes are splittable.
 */
static
int best_splittable_box(const box* bv, int boxes)
{
	int bi = -1;
	double maxsum = 0;
	for (int i=0; i<boxes; i++) {
		const box& b = bv[i];
		if (b.colors < 2) continue;

		double thissum = b.sum * b.variance;
		if (thissum > maxsum) {
			maxsum = thissum;
			bi = i;
		}
	}
	return bi;
}

static inline
double color_weight(f_pixel median, const hist_item& h)
{
	double diff = colordifference(median, h.acolor);
	// if color is "good enough", don't split further
	if (diff < 1.0/256.0) diff /= 2.0;
	return sqrt(diff) * sqrt(h.adjusted_weight);
}

static
std::vector<colormap_item> colormap_from_boxes(const box* bv, int boxes, const hist_item* achv, double min_opaque_val)
{
	/*
	 ** Ok, we've got enough boxes.	 Now choose a representative color for
	 ** each box.  There are a number of possible ways to make this choice.
	 ** One would be to choose the center of the box; this ignores any structure
	 ** within the boxes.  Another method would be to average all the colors in
	 ** the box - this is the method specified in Heckbert's paper.
	 */

	std::vector<colormap_item> map(boxes);

	for (int bi=0; bi<boxes; ++bi) {
		colormap_item& cm = map[bi];
		const box& bx = bv[bi];
		cm.acolor = averagepixels(bx.ind, bx.colors, achv, min_opaque_val);

		/* store total color popularity (perceptual_weight is approximation of it) */
		cm.popularity = 0;
		for (int i=bx.ind; i<bx.ind+bx.colors; i++) {
			cm.popularity += achv[i].perceptual_weight;
		}
	}

	return map;
}

/* increase histogram popularity by difference from the final color (this is used as part of feedback loop) */
static
void adjust_histogram(hist_item* achv, const std::vector<colormap_item>& map, const box* bv, int boxes)
{
	for (int bi=0; bi<boxes; ++bi) {
		const box& bx = bv[bi];
		f_pixel pc = map[bi].acolor;
		for (int i=bx.ind; i<bx.ind+bx.colors; i++) {
			hist_item& hist = achv[i];
			hist.adjusted_weight *= 1.0 + sqrt(colordifference(pc, hist.acolor)) / 2.0;
		}
	}
}

/*
 ** Here is the fun part, the median-cut colormap generator.  This is based
 ** on Paul Heckbert's paper, "Color Image Quantization for Frame Buffer
 ** Display," SIGGRAPH 1982 Proceedings, page 297.
 */
std::vector<colormap_item> mediancut(
	std::vector<hist_item>& hist, double min_opaque_val, int newcolors
	)
{
	std::vector<box> bv(newcolors);

	/*
	 ** Set up the initial box.
	 */
	bv[0].ind = 0;
	bv[0].colors = hist.size();
	bv[0].variance = 1.0;
	for (int i=0; i<bv[0].colors; i++)
		bv[0].sum += hist[i].adjusted_weight;

	int boxes = 1;

	/*
	 ** Main loop: split boxes until we have enough.
	 */
	while (boxes < newcolors) {

		int bi= best_splittable_box(&bv[0], boxes);
		if (bi < 0)
			break;		  /* ran out of colors! */
		
		box& bx = bv[bi];
		int indx = bx.ind;
		int clrs = bx.colors;

		sort_colors_by_variance(channel_variance(&hist[0], indx, clrs, min_opaque_val), &hist[0], indx, clrs);

		/*
		 Classic implementation tries to get even number of colors or pixels in each subdivision.

		 Here, instead of popularity I use (sqrt(popularity)*variance) metric.
		 Each subdivision balances number of pixels (popular colors) and low variance -
		 boxes can be large if they have similar colors. Later boxes with high variance
		 will be more likely to be split.

		 Median used as expected value gives much better results than mean.
		 */

		f_pixel median = averagepixels(indx+(clrs-1)/2, clrs&1 ? 1 : 2, &hist[0], min_opaque_val);

		int lowersum = 0;
		double halfvar = 0, lowervar = 0;
		for (int i=0; i<clrs-1; i++) {
			halfvar += color_weight(median, hist[indx+i]);
		}
		halfvar /= 2.0;

		int break_at;
		for (break_at=0; break_at<clrs-1; ++break_at) {
			if (lowervar >= halfvar)
				break;
			hist_item& hi = hist[indx+break_at];
			lowervar += color_weight(median, hi);
			lowersum += hi.adjusted_weight;
		}

		/*
		 ** Split the box. Sum*variance is then used to find "largest" box to split.
		 */
		int sm = bx.sum;
		bx.colors = break_at;
		bx.sum = lowersum;
		bx.variance = lowervar;
		box& bx2 = bv[boxes];
		bx2.ind = indx + break_at;
		bx2.colors = clrs - break_at;
		bx2.sum = sm - lowersum;
		bx2.variance = halfvar*2.0-lowervar;
		++boxes;
	}

	std::vector<colormap_item>& map = colormap_from_boxes(&bv[0], boxes, &hist[0], min_opaque_val);
	adjust_histogram(&hist[0], map, &bv[0], boxes);

	return map;
}

