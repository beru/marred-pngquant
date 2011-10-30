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
void transposing_1d_blur(double* src, double* dst, int width, int height, const int size)
{
	const double sizef = size;

	for (int j=0; j<height; j++) {
		double* row = src + j*width;

		// accumulate sum for pixels outside line
		double sum;
		sum = row[0] * sizef;
		for (int i=0; i<size; i++) {
			sum += row[i];
		}
		
		// blur with left side outside line
		for (int i=0; i<size; i++) {
			sum -= row[0];
			sum += row[i+size];
			
			dst[i*height + j] = sum / (sizef*2.f);
		}
		
		for (int i=size; i<width-size; i++) {
			sum -= row[i-size];
			sum += row[i+size];
			
			dst[i*height + j] = sum / (sizef*2.f);
		}
		
		// blur with right side outside line
		for (int i=width-size; i<width; i++) {
			sum -= row[i-size];
			sum += row[width-1];
			
			dst[i*height + j] = sum/(sizef*2.0);
		}
	}
}

void max3(double* src, double* dst, int width, int height)
{
	for (int j=0; j<height; j++) {
		const double* row = src + j*width,
		*prevrow = src + MAX(0,j-1)*width,
		*nextrow = src + MIN(height-1,j+1)*width;
		
		double prev;
		double curr = row[0];
		double next = row[0];
		
		for (int i=0; i<width-1; i++) {
			prev = curr;
			curr = next;
			next = row[i+1];
			
			double t1 = MAX(prev, next);
			double t2 = MAX(nextrow[i], prevrow[i]);
			*dst++ = MAX(curr,MAX(t1, t2));
		}
		double t1 = MAX(curr, next);
		double t2 = MAX(nextrow[width-1], prevrow[width-1]);
		*dst++ = MAX(t1, t2);
	}
}

void min3(double* src, double* dst, int width, int height)
{
	for (int j=0; j<height; j++) {
		const double* row = src + j*width,
		*prevrow = src + MAX(0,j-1)*width,
		*nextrow = src + MIN(height-1,j+1)*width;
		
		double prev;
		double curr = row[0];
		double next = row[0];
		
		for (int i=0; i<width-1; i++) {
			prev = curr;
			curr = next;
			next = row[i+1];
			
			double t1 = MIN(prev, next);
			double t2 = MIN(nextrow[i], prevrow[i]);
			*dst++ = MIN(curr, MIN(t1,t2));
		}
		double t1 = MIN(curr, next);
		double t2 = MIN(nextrow[width-1], prevrow[width-1]);
		*dst++ = MIN(t1,t2);
	}
}

/*
 Filters image with callback and blurs (lousy approximate of gaussian)
 */
void blur(double* src, double* tmp, double* dst, int width, int height, int size)
{
	transposing_1d_blur(src, tmp, width, height, size);
	transposing_1d_blur(tmp, dst, height, width, size);
}

