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

static std::vector<hist_item> pam_acolorhashtoacolorhist(acolorhash_table* acht, int maxacolors);
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
	unsigned long hash = px.alpha * 256.0*5.0 + px.r * 256.0*179.0 + px.g * 256.0*17.0 + px.b * 256.0*30047.0;
	return hash % HASH_SIZE;
}

static inline
f_pixel posterize_pixel(f_pixel px, int maxval)
{
//	if (maxval == 255) {
		return px;
//	}else {
//		return to_f_scalar(
//			gamma,
//			(px * maxval / 255) / (double)maxval
//		);
//	}
}

static
acolorhash_table* pam_computeacolorhash(
	const f_pixel* input, size_t width, size_t height,
	int maxacolors, int ignorebits, const double* importance_map, int* acolorsP
	)
{
	
	acolorhist_list_item* achl;
	const int maxval = 255 >> ignorebits;
	acolorhash_table* acht = pam_allocacolorhash();
	acolorhist_list_item** buckets = acht->buckets;
	int colors = 0;

	/* Go through the entire image, building a hash table of colors. */
	const f_pixel* pLine = input;
	for (size_t y=0; y<height; ++y) {
		for (size_t x=0; x<width; ++x) {
			double boost = 0.5 + *importance_map++;
			f_pixel curr = posterize_pixel(pLine[x], maxval);
			int hash = pam_hashapixel(curr);
			for (achl=buckets[hash]; achl!=NULL; achl=achl->next) {
				if (achl->acolor == curr) {
					break;
				}
			}
			if (achl != NULL) {
				achl->perceptual_weight += boost;
			}else {
				if (++colors > maxacolors) {
					pam_freeacolorhash(acht);
					return NULL;
				}
				achl = (acolorhist_list_item*) mempool_new(acht->mempool, sizeof(acolorhist_list_item));

				achl->acolor = curr;
				achl->perceptual_weight = boost;
				achl->next = buckets[hash];
				buckets[hash] = achl;
			}
		}
		pLine += width;
	}
	*acolorsP = colors;
	return acht;
}

/**
 * Builds color histogram no larger than maxacolors. Ignores (posterizes) ignorebits lower bits in each color.
 * perceptual_weight of each entry is increased by value from importance_map
 */
std::vector<hist_item> pam_computeacolorhist(
	const f_pixel* input, size_t width, size_t height,
	int maxacolors, int ignorebits,
	const double* importance_map
	)
{
	int hist_size = 0;
	acolorhash_table* acht = pam_computeacolorhash(input, width, height, maxacolors, ignorebits, importance_map, &hist_size);
	if (!acht) return std::vector<hist_item>();

	std::vector<hist_item> ret = pam_acolorhashtoacolorhist(acht, hist_size);
	pam_freeacolorhash(acht);
	return ret;
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
std::vector<hist_item> pam_acolorhashtoacolorhist(
	acolorhash_table* acht, int hist_size
	)
{
	std::vector<hist_item> ret(hist_size);
	acolorhist_list_item* achl;
	/* Loop through the hash table. */
	size_t j = 0;
	for (size_t i=0; i<HASH_SIZE; ++i)
		for (achl=acht->buckets[i]; achl!=NULL; achl=achl->next) {
			hist_item& hist = ret[j];
			/* Add the new entry. */
			hist.acolor = achl->acolor;
			hist.adjusted_weight = hist.perceptual_weight = achl->perceptual_weight;
			++j;
		}

	/* All done. */
	return ret;
}


static
void pam_freeacolorhash(acolorhash_table* acht)
{
	mempool_free(acht->mempool);
}

