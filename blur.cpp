
#include "pam.h"
#include "blur.h"

/*
 Blurs image horizontally (width 2*size+1) and writes it transposed to dst (called twice gives 2d blur)
 */
static
void transposing_1d_blur(
	const double* src,
	double* dst,
	unsigned int width,
	unsigned int height,
	const unsigned int size
	)
{
	const double sizef = size;
	const double sizef2 = 1.0 / (sizef*2.0);
	for (unsigned int j=0; j<height; ++j) {
		const double* row = src + j*width;
		double* dstLine = dst + j;

		// accumulate sum for pixels outside line
		double sum;
		sum = row[0] * sizef;
		for (unsigned int i=0; i<size; ++i) {
			sum += row[i];
		}
		
		// blur with left side outside line
		for (unsigned int i=0; i<size; ++i) {
			sum -= row[0];
			sum += row[i+size];
			dstLine[i*height] = sum * sizef2;
		}
		
		for (unsigned int i=size; i<width-size; ++i) {
			sum -= row[i-size];
			sum += row[i+size];
			dstLine[i*height] = sum * sizef2;
		}
		
		// blur with right side outside line
		for (unsigned int i=width-size; i<width; ++i) {
			sum -= row[i-size];
			sum += row[width-1];
			dstLine[i*height] = sum * sizef2;
		}
	}
}

/**
 * Picks maximum of neighboring pixels (blur + lighten)
 */
void max3(
	const double* src,
	double* dst,
	unsigned int width,
	unsigned int height
	)
{
	for (unsigned int j=0; j<height; ++j) {
		const double* row = src + j*width,
		*prevrow = src + (j>1 ? j-1 : 0)*width,
		*nextrow = src + min(height-1,j+1)*width;
		
		double prev;
		double curr = row[0];
		double next = row[0];
		
		for (unsigned int i=0; i<width-1; ++i) {
			prev = curr;
			curr = next;
			next = row[i+1];
			*dst++ = max(curr, prev, next, nextrow[i], prevrow[i]);
		}
		*dst++ = max(curr, next, nextrow[width-1], prevrow[width-1]);
	}
}

/**
 * Picks minimum of neighboring pixels (blur + darken)
 */
void min3(
	const double* src,
	double* dst,
	unsigned int width,
	unsigned int height
	)
{
	for (unsigned int j=0; j<height; ++j) {
		const double* row = src + j*width,
		*prevrow = src + (j>1 ? j-1 : 0)*width,
		*nextrow = src + min(height-1,j+1)*width;
		
		double prev;
		double curr = row[0];
		double next = row[0];
		
		for (unsigned int i=0; i<width-1; ++i) {
			prev = curr;
			curr = next;
			next = row[i+1];
			
			*dst++ = min(curr, prev, next, nextrow[i], prevrow[i]);
		}
		*dst++ = min(curr, next, nextrow[width-1], prevrow[width-1]);
	}
}

/*
 Filters src image and saves it to dst, overwriting tmp in the process.
 Image must be width*height pixels high. Size controls radius of box blur.
 */
void blur(
	const double* src,
	double* tmp,
	double* dst,
	unsigned int width,
	unsigned int height,
	unsigned int size
	)
{
    if (width<2*size+1 || height<2*size+1) return;
	transposing_1d_blur(src, tmp, width, height, size);
	transposing_1d_blur(tmp, dst, height, width, size);
}

