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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>

#include "pam.h"

static hist* pam_acolorhashtoacolorhist(acolorhash_table* acht, int maxacolors);
static acolorhash_table* pam_computeacolorhash(const rgb_pixel* const* apixels, int cols, int rows, double gamma, int maxacolors, int ignorebits, int use_contrast, int* acolorsP);
static void pam_freeacolorhash(acolorhash_table* acht);
static acolorhash_table* pam_allocacolorhash(void);

/* libpam3.c - pam (portable alpha map) utility library part 3
 **
 ** Colormap routines.
 **
 ** Copyright (C) 1989, 1991 by Jef Poskanzer.
 ** Copyright (C) 1997 by Greg Roelofs.
 **
 ** Permission to use, copy, modify, and distribute this software and its
 ** documentation for any purpose and without fee is hereby granted, provided
 ** that the above copyright notice appear in all copies and that both that
 ** copyright notice and this permission notice appear in supporting
 ** documentation.	This software is provided "as is" without express or
 ** implied warranty.
 */

#define HASH_SIZE 30029

struct mempool {
	mempool* next;
	size_t used;
};

#define MEMPOOL_RESERVED ((sizeof(mempool)+15) & ~0xF)
#define MEMPOOL_SIZE (1<<18)

static
void* mempool_new(mempool*& rpm, size_t size)
{
	assert(size < MEMPOOL_SIZE-MEMPOOL_RESERVED);
	
	if (rpm && (rpm->used+size) <= MEMPOOL_SIZE) {
		int prevused = rpm->used;
		rpm->used += (size+15) & ~0xF;
		return ((char*)rpm) + prevused;
	}

	mempool* old = rpm ? rpm : NULL;
	char *mem = (char*) calloc(MEMPOOL_SIZE, 1);

	rpm = (mempool*)mem;
	rpm->used = MEMPOOL_RESERVED;
	rpm->next = old;

	return mempool_new(rpm, size);
}

static
void mempool_free(mempool* m)
{
	while (m) {
		mempool* next = m->next;
		free(m);
		m = next;
	}
}

static inline
unsigned long pam_hashapixel(f_pixel px)
{
	unsigned long hash = px.a * 256.0*5.0 + px.r * 256.0*179.0 + px.g * 256.0*17.0 + px.b * 256.0*30047.0;
	return hash % HASH_SIZE;
}

hist* pam_computeacolorhist(const rgb_pixel*const apixels[], int cols, int rows, double gamma, int maxacolors, int ignorebits, int use_contrast)
{
	int hist_size = 0;
	acolorhash_table* acht = pam_computeacolorhash(apixels, cols, rows, gamma, maxacolors, ignorebits, use_contrast, &hist_size);
	if (!acht) return 0;

	hist* achv = pam_acolorhashtoacolorhist(acht, hist_size);
	pam_freeacolorhash(acht);
	return achv;
}

static inline
f_pixel posterize_pixel(rgb_pixel px, int maxval, double gamma)
{
	if (maxval == 255) {
		return to_f(gamma, px);
	}else {
		return to_f_scalar(
			gamma,
			f_pixel(
				(px.a * maxval / 255) / (double)maxval,
				(px.r * maxval / 255) / (double)maxval,
				(px.g * maxval / 255) / (double)maxval,
				(px.b * maxval / 255) / (double)maxval
				)
		);
	}
}

double boost_from_contrast(f_pixel prev, f_pixel fpx, f_pixel next, f_pixel above, f_pixel below, double prev_boost)
{
	f_pixel fpx2 = fpx * 2.0;
	f_pixel c0 = (fpx2 - (prev + next)).abs();
	f_pixel c1 = (fpx2 - (above + below)).abs();
	f_pixel c = max(c0, c1); // maximum of vertical or horizontal contrast. It's better at picking up noise.
	double contrast = c.r + c.g + c.b + 2.0*c.a; // oddly a*2 works better than *3
	// 0.6-contrast avoids boosting flat areas
	contrast = MIN(1.0, fabsf(0.6-contrast/3.0)) * (1.0/0.6);
	// I want only really high contrast (noise, edges) to influence
	contrast *= contrast;

	// it's "smeared" to spread influence of edges to neighboring pixels
	return (contrast < prev_boost) ? contrast : (prev_boost+prev_boost+contrast)/3.0;
}

static
acolorhash_table* pam_computeacolorhash(const rgb_pixel*const* apixels, int cols, int rows, double gamma, int maxacolors, int ignorebits, int use_contrast, int* acolorsP)
{
	acolorhist_list_item* achl;
	const int maxval = 255 >> ignorebits;
	acolorhash_table* acht = pam_allocacolorhash();
	acolorhist_list_item** buckets = acht->buckets;
	int colors = 0;

	/* Go through the entire image, building a hash table of colors. */
	for (int row=0; row<rows; ++row) {
		const rgb_pixel* prevline = apixels[MAX(0,row-1)];
		const rgb_pixel* nextline = apixels[MIN(rows-1,row+1)];
		const rgb_pixel* curline = apixels[row];
		f_pixel curr = posterize_pixel(curline[0], maxval, gamma);
		f_pixel next = posterize_pixel(curline[MIN(cols-1,1)], maxval, gamma);
		f_pixel prev;
		double boost = 0.5;
		for (int col=0; col<cols; ++col) {
			prev = curr;
			curr = next;
			next = posterize_pixel(curline[MIN(cols-1,col+1)], maxval, gamma);

			if (use_contrast) {
				f_pixel above = posterize_pixel(prevline[col], maxval, gamma);
				f_pixel below = posterize_pixel(nextline[col], maxval, gamma);

				boost = boost_from_contrast(prev,curr,next,above,below,boost);
			}

			int hash = pam_hashapixel(curr);

			for (achl=buckets[hash]; achl!=NULL; achl=achl->next) {
				if (achl->acolor == curr) {
					break;
				}
			}
			if (achl != NULL) {
				achl->perceptual_weight += 1.0f + boost;
			}else {
				if (++colors > maxacolors) {
					pam_freeacolorhash(acht);
					return NULL;
				}
				achl = (acolorhist_list_item*) mempool_new(acht->mempool, sizeof(acolorhist_list_item));

				achl->acolor = curr;
				achl->perceptual_weight = 1.0f + boost;
				achl->next = buckets[hash];
				buckets[hash] = achl;
			}
		}

	}
	*acolorsP = colors;
	return acht;
}

static
acolorhash_table* pam_allocacolorhash()
{
	mempool* m = NULL;
	acolorhash_table* t = (acolorhash_table*) mempool_new(m, sizeof(*t));
	t->buckets = (acolorhist_list_item**) mempool_new(m, HASH_SIZE * sizeof(t->buckets[0]));
	t->mempool = m;
	return t;
}

static
hist* pam_acolorhashtoacolorhist(acolorhash_table* acht, int hist_size)
{
	acolorhist_list_item* achl;
	hist* ph = (hist*) malloc(sizeof(hist));
	ph->achv = (hist_item*) malloc(hist_size * sizeof(ph->achv[0]));
	ph->size = hist_size;

	/* Loop through the hash table. */
	int j = 0;
	for (int i=0; i<HASH_SIZE; ++i)
		for (achl=acht->buckets[i]; achl!=NULL; achl=achl->next) {
			hist_item& hist = ph->achv[j];
			/* Add the new entry. */
			hist.acolor = achl->acolor;
			hist.adjusted_weight = hist.perceptual_weight = achl->perceptual_weight;
			++j;
		}

	/* All done. */
	return ph;
}


static
void pam_freeacolorhash(acolorhash_table* acht)
{
	mempool_free(acht->mempool);
}

void pam_freeacolorhist(hist* hist)
{
	free(hist->achv);
	free(hist);
}

colormap* pam_colormap(int colors)
{
	colormap* map = (colormap*) malloc(sizeof(colormap));
	map->palette = (colormap_item*) calloc(colors, sizeof(map->palette[0]));
	map->colors = colors;
	return map;
}

void pam_freecolormap(colormap* c)
{
	free(c->palette);
	free(c);
}


