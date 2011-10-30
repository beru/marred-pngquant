//
//	blur.c
//	pngquant

#include <assert.h>
#include "pam.h"
#include "blur.h"

/*
 Blurs image horizontally (width 2*size+1) and writes it transposed to dst (called twice gives 2d blur)
 */
static
void transposing_1d_blur(const double* src, double* dst, int width, int height, const int size)
{
	const double sizef = size;
	const double sizef2 = 1.0 / (sizef*2.0);
	for (int j=0; j<height; ++j) {
		const double* row = src + j*width;
		double* dstLine = dst + j;

		// accumulate sum for pixels outside line
		double sum;
		sum = row[0] * sizef;
		for (int i=0; i<size; ++i) {
			sum += row[i];
		}
		
		// blur with left side outside line
		for (int i=0; i<size; ++i) {
			sum -= row[0];
			sum += row[i+size];
			dstLine[i*height] = sum * sizef2;
		}
		
		for (int i=size; i<width-size; ++i) {
			sum -= row[i-size];
			sum += row[i+size];
			dstLine[i*height] = sum * sizef2;
		}
		
		// blur with right side outside line
		for (int i=width-size; i<width; ++i) {
			sum -= row[i-size];
			sum += row[width-1];
			dstLine[i*height] = sum * sizef2;
		}
	}
}

void max3(const double* src, double* dst, int width, int height)
{
	for (int j=0; j<height; ++j) {
		const double* row = src + j*width,
		*prevrow = src + max(0,j-1)*width,
		*nextrow = src + min(height-1,j+1)*width;
		
		double prev;
		double curr = row[0];
		double next = row[0];
		
		for (int i=0; i<width-1; ++i) {
			prev = curr;
			curr = next;
			next = row[i+1];
			*dst++ = max(curr, prev, next, nextrow[i], prevrow[i]);
		}
		*dst++ = max(curr, next, nextrow[width-1], prevrow[width-1]);
	}
}

void min3(const double* src, double* dst, int width, int height)
{
	for (int j=0; j<height; ++j) {
		const double* row = src + j*width,
		*prevrow = src + max(0,j-1)*width,
		*nextrow = src + min(height-1,j+1)*width;
		
		double prev;
		double curr = row[0];
		double next = row[0];
		
		for (int i=0; i<width-1; ++i) {
			prev = curr;
			curr = next;
			next = row[i+1];
			
			*dst++ = min(curr, prev, next, nextrow[i], prevrow[i]);
		}
		*dst++ = min(curr, next, nextrow[width-1], prevrow[width-1]);
	}
}

/*
 Filters image with callback and blurs (lousy approximate of gaussian)
 */
void blur(const double* src, double* tmp, double* dst, int width, int height, int size)
{
	transposing_1d_blur(src, tmp, width, height, size);
	transposing_1d_blur(tmp, dst, height, width, size);
}

