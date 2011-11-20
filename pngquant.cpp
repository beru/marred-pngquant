/* pngquant.c - quantize the colors in an alphamap down to a specified number
**
** Copyright (C) 1989, 1991 by Jef Poskanzer.
** Copyright (C) 1997, 2000, 2002 by Greg Roelofs; based on an idea by
**								  Stefan Schneider.
** (C) 2011 by Kornel Lesinski.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
*/

/* GRR TO DO:  "original file size" and "quantized file size" if verbose? */
/* GRR TO DO:  add option to preserve background color (if any) exactly */
/* GRR TO DO:  add mapfile support, but cleanly (build palette in main()) */
/* GRR TO DO:  support 16 bps without down-conversion */
/* GRR TO DO:  if all samples are gray and image is opaque and sample depth
				would be no bigger than palette and user didn't explicitly
				specify a mapfile, switch to grayscale */
/* GRR TO DO:  if all samples are 0 or maxval, eliminate gAMA chunk (rwpng.c) */

#define PNGQUANT_VERSION "0.0.1 (October 2011)"

#define PNGQUANT_USAGE "\
   usage:  pngquant [options] [ncolors] [pngfile [pngfile ...]]\n\n\
   options:\n\
	  -force		overwrite existing output files (synonym: -f)\n\
	  -ext new.png	set custom extension for output filename\n\
	  -nofs			disable dithering (synonyms: -nofloyd, -ordered)\n\
	  -verbose		print status messages (synonyms: -noquiet)\n\
	  -speed N		speed/quality trade-off. 1=slow, 3=default, 10=fast & rough\n\
\n\
   Quantizes one or more 32-bit RGBA PNGs to 8-bit (or smaller) RGBA-palette\n\
   PNGs using Floyd-Steinberg diffusion dithering (unless disabled).\n\
   The output filename is the same as the input name except that\n\
   it ends in \"-fs8.png\", \"-or8.png\" or your custom extension (unless the\n\
   input is stdin, in which case the quantized image will go to stdout).\n\
   The default behavior if the output file exists is to skip the conversion;\n\
   use -force to overwrite.\n"


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdarg.h>
#ifdef WIN32		/* defined in Makefile.w32 (or use _MSC_VER for MSVC) */
#  include <fcntl.h>	/* O_BINARY */
#  include <io.h>	/* setmode() */
#endif

#include <stddef.h>

#include "png.h"	/* libpng header; includes zlib.h */
#include "rwpng.h"	/* typedefs, common macros, public prototypes */
#include "pam.h"
#include "mediancut.h"
#include "remap.h"
#include "blur.h"

#include <vector>
#include <algorithm>

pngquant_error pngquant(read_info* input_image, write_info* output_image, bool floyd, int reqcolors, int speed_tradeoff);
pngquant_error read_image(const char* filename, bool using_stdin, read_info* input_image_p);
pngquant_error write_image(write_info* output_image, const char* filename, const char* newext, bool force, bool using_stdin);

static bool verbose = false;

void verbose_printf(const char* fmt, ...)
{
	va_list va;
	va_start(va, fmt);
	if (verbose) vfprintf(stderr, fmt, va);
	va_end(va);
}

static
void print_full_version(FILE* fd)
{
	fprintf(fd, "pngquant-marred, version %s, by berupon@gmail.com.\n", PNGQUANT_VERSION);
	fprintf(fd, "origina pngquant, by Greg Roelofs, Kornel Lesinski(porneL).\n");
	rwpng_version_info(fd);
	fputs("\n", fd);
}

static
void print_usage(FILE* fd)
{
	fputs(PNGQUANT_USAGE, fd);
}

int main(int argc, char* argv[])
{
	int argn;
	int reqcolors;
	bool floyd = true;
	bool force = false;
	int speed_tradeoff = 3; // 1 max quality, 10 rough & fast. 3 is optimum.
	bool using_stdin = false;
	int latest_error=0, error_count=0, file_count=0;
	const char* filename;
	const char* newext = NULL;

	argn = 1;

	while (argn < argc && argv[argn][0] == '-' && argv[argn][1] != '\0') {
		if (0 == strcmp(argv[argn], "--")) { ++argn;break; }

		if ( 0 == strncmp(argv[argn], "-fs", 3) ||
			 0 == strncmp(argv[argn], "-floyd", 3) )
			floyd = true;
		else if ( 0 == strncmp(argv[argn], "-nofs", 5) ||
				  0 == strncmp(argv[argn], "-nofloyd", 5) ||
				  0 == strncmp(argv[argn], "-ordered", 3) )
			floyd = false;
		else if (0 == strncmp(argv[argn], "-force", 2))
			force = true;
		else if (0 == strncmp(argv[argn], "-noforce", 4))
			force = false;
		else if ( 0 == strcmp(argv[argn], "-verbose") ||
				  0 == strcmp(argv[argn], "-v") ||
				  0 == strncmp(argv[argn], "-noquiet", 4) )
			verbose = true;
		else if ( 0 == strncmp(argv[argn], "-noverbose", 4) ||
				  0 == strncmp(argv[argn], "-quiet", 2) )
			verbose = false;

		else if ( 0 == strcmp(argv[argn], "-version")) {
			puts(PNGQUANT_VERSION);
			return SUCCESS;
		}else if ( 0 == strcmp(argv[argn], "-h") || 0 == strcmp(argv[argn], "--help")) {
			print_full_version(stdout);
			print_usage(stdout);
			return SUCCESS;
		}else if (0 == strcmp(argv[argn], "-ext")) {
			++argn;
			if (argn == argc) {
				print_usage(stderr);
				return MISSING_ARGUMENT;
			}
			newext = argv[argn];
		}else if (0 == strcmp(argv[argn], "-s") ||
				  0 == strcmp(argv[argn], "-speed")) {
			++argn;
			if (argn == argc) {
				print_usage(stderr);
				return MISSING_ARGUMENT;
			}
			speed_tradeoff = atoi(argv[argn]);
		}else {
			print_usage(stderr);
			return MISSING_ARGUMENT;
		}
		++argn;
	}

	if (argn == argc) {
		print_full_version(stderr);
		print_usage(stderr);
		return MISSING_ARGUMENT;
	}
	if (sscanf(argv[argn], "%d", &reqcolors) != 1) {
		reqcolors = 256; argn--;
	}
	if (reqcolors <= 1) {
		fputs("number of colors must be greater than 1\n", stderr);
		return INVALID_ARGUMENT;
	}
	if (reqcolors > 256) {
		fputs("number of colors cannot be more than 256\n", stderr);
		return INVALID_ARGUMENT;
	}
	if (speed_tradeoff < 1 || speed_tradeoff > 10) {
		fputs("speed should be between 1 (slow) and 10 (fast)\n", stderr);
		return INVALID_ARGUMENT;
	}
	++argn;

	if (newext == NULL) {
		newext = floyd ? "-fs8.png" : "-or8.png";
	}

	if (argn == argc || 0==strcmp(argv[argn],"-")) {
		using_stdin = true;
		filename = "stdin";
	}else {
		filename = argv[argn];
		++argn;
	}


	/*=============================	 MAIN LOOP	=============================*/

	while (argn <= argc) {
		int retval;

		verbose_printf("%s:\n", filename);

		read_info input_image = {{0}};
		write_info output_image = {{0}};
		retval = read_image(filename,using_stdin, &input_image);

		if (!retval) {
			retval = pngquant(&input_image, &output_image, floyd, reqcolors, speed_tradeoff);
		}

		/* now we're done with the INPUT data and row_pointers, so free 'em */
		if (input_image.rgba_data) {
			free(input_image.rgba_data);
		}
		if (input_image.row_pointers) {
			free(input_image.row_pointers);
		}

		if (!retval) {
			verbose_printf("  writing %d-color image\n", output_image.num_palette);

			retval = write_image(&output_image,filename,newext,force,using_stdin);
		}

		if (output_image.indexed_data) {
			free(output_image.indexed_data);
		}
		if (output_image.row_pointers) {
			free(output_image.row_pointers);
		}

		if (retval) {
			latest_error = retval;
			++error_count;
		}
		++file_count;

		verbose_printf("\n");

		filename = argv[argn];
		++argn;
	}

	/*=======================================================================*/


	if (error_count)
		verbose_printf("There were errors quantizing %d file%s out of a"
		  " total of %d file%s.\n",
		  error_count, (error_count == 1)? "" : "s",
		  file_count, (file_count == 1)? "" : "s");
	else
		verbose_printf("No errors detected while quantizing %d image%s.\n",
		  file_count, (file_count == 1)? "" : "s");

	return latest_error;
}

static inline
bool compare_popularity(const colormap_item& v1, const colormap_item& v2)
{
	return v1.popularity > v2.popularity;
}

static
void sort_palette(write_info* output_image, std::vector<colormap_item>& map)
{
	assert(output_image);
	
	/*
	** Step 3.5 [GRR]: remap the palette colors so that all entries with
	** the maximal alpha value (i.e., fully opaque) are at the end and can
	** therefore be omitted from the tRNS chunk.
	*/
	
	verbose_printf("  eliminating opaque tRNS-chunk entries...");
	
	/* move transparent colors to the beginning to shrink trns chunk */
	int num_transparent=0;
	for (int i=0; i<map.size(); i++) {
		colormap_item& pal = map[i];
		if (pal.acolor.alpha < 1.0) {
			if (i != num_transparent) {
				std::swap(map[num_transparent], pal);
				--i;
			}
			++num_transparent;
		}
	}
	
	verbose_printf("%d entr%s transparent\n", num_transparent, (num_transparent == 1)? "y" : "ies");
	
	/* colors sorted by popularity make pngs slightly more compressible
	 * opaque and transparent are sorted separately
	 */
	std::sort(map.begin(), map.begin()+num_transparent, compare_popularity);
	if (num_transparent < map.size()) {
		std::sort(map.begin()+num_transparent, map.begin()+map.size()-num_transparent, compare_popularity);
	}
	
	output_image->num_palette = map.size();
	output_image->num_trans = num_transparent;
}

static
void set_palette(write_info* output_image, std::vector<colormap_item>& map)
{
	for (size_t x=0; x<map.size(); ++x) {
		colormap_item& pal = map[x];
		rgb_pixel px = to_rgb(output_image->gamma, lab2rgb(pal.acolor));
		pal.acolor = rgb2lab(to_f(output_image->gamma, px)); /* saves rounding error introduced by to_rgb, which makes remapping & dithering more accurate */
//		rgb_pixel px = to_rgb(output_image->gamma, pal.acolor);
//		pal.acolor = to_f(output_image->gamma, px); /* saves rounding error introduced by to_rgb, which makes remapping & dithering more accurate */
		
		png_color& pc = output_image->palette[x];
		pc.red	 = px.r;
		pc.green = px.g;
		pc.blue	 = px.b;
		output_image->trans[x] = px.a;
	}
}

/* build the output filename from the input name by inserting "-fs8" or
 * "-or8" before the ".png" extension (or by appending that plus ".png" if
 * there isn't any extension), then make sure it doesn't exist already */
char* add_filename_extension(const char* filename, const char* newext)
{
	int x = strlen(filename);

	char* outname = (char*) malloc(x+4+strlen(newext)+1);

	strncpy(outname, filename, x);
	if (strncmp(outname+x-4, ".png", 4) == 0)
		strcpy(outname+x-4, newext);
	else
		strcpy(outname+x, newext);

	return outname;
}

static
void set_binary_mode(FILE* fp)
{
#if defined(MSDOS) || defined(FLEXOS) || defined(OS2) || defined(WIN32)
#if (defined(__HIGHC__) && !defined(FLEXOS))
	setmode(fp, _BINARY);
#else
	_setmode(fp == stdout ? 1 : 0, O_BINARY);
#endif
#endif
}

pngquant_error write_image(
	write_info* output_image, const char* filename,
	const char* newext, bool force, bool using_stdin
	)
{
	FILE* outfile;
	if (using_stdin) {
		set_binary_mode(stdout);
		outfile = stdout;
	}else {
		char* outname = add_filename_extension(filename, newext);
		
		if (!force) {
			if ((outfile = fopen(outname, "rb")) != NULL) {
				fprintf(stderr, "  error:  %s exists; not overwriting\n", outname);
				fclose(outfile);
				free(outname);
				return NOT_OVERWRITING_ERROR;
			}
		}
		if ((outfile = fopen(outname, "wb")) == NULL) {
			fprintf(stderr, "  error:  cannot open %s for writing\n", outname);
			free(outname);
			return CANT_WRITE_ERROR;
		}
		free(outname);
	}
	
	pngquant_error retval = rwpng_write_image_init(outfile, output_image);
	if (retval) {
		fprintf(stderr, "  rwpng_write_image_init() error\n");
		if (!using_stdin)
			fclose(outfile);
		return retval;
	}
	
	/* write entire interlaced palette PNG */
	retval = rwpng_write_image_whole(output_image);
	
	if (!using_stdin)
		fclose(outfile);
	
	/* now we're done with the OUTPUT data and row_pointers, too */
	return retval;
}

static
std::vector<hist_item> histogram(
	const f_pixel* input, size_t width, size_t height,
	int reqcolors, int speed_tradeoff,
	const double* importance_map
	)
{
	int ignorebits = 0;
   /*
	** Step 2: attempt to make a histogram of the colors, unclustered.
	** If at first we don't succeed, increase ignorebits to increase color
	** coherence and try again.
	*/

	if (speed_tradeoff > 7) ignorebits++;
	int maxcolors = (1<<17) + (1<<18)*(10-speed_tradeoff);

	verbose_printf("  making histogram...");
	std::vector<hist_item> hist;
	for (;;) {
		hist = pam_computeacolorhist(input, width, height, maxcolors, ignorebits, importance_map);
		if (hist.size()) {
			break;
		}
		ignorebits++;
		verbose_printf("too many colors!\n	scaling colors to improve clustering...");
	}

	verbose_printf("%d colors found\n", hist.size());
	return hist;
}

double modify_alpha(read_info* input_image)
{
	/* IE6 makes colors with even slightest transparency completely transparent,
	   thus to improve situation in IE, make colors that are less than ~10% transparent
	   completely opaque */

	rgb_pixel** input_pixels = (rgb_pixel**) input_image->row_pointers;
	rgb_pixel* pP;
	const int rows = input_image->height;
	const int cols = input_image->width;
	double gamma = input_image->gamma;
	double min_opaque_val = 1.0;

	for (int row=0; row<rows; ++row) {
		pP = input_pixels[row];
		for (int col=0; col<cols; ++col, ++pP) {
			f_pixel px = to_f(gamma, *pP);

#ifndef NDEBUG
			rgb_pixel rgbcheck = to_rgb(gamma, px);


			if (pP->a && (pP->r != rgbcheck.r || pP->g != rgbcheck.g || pP->b != rgbcheck.b || pP->a != rgbcheck.a)) {
				fprintf(stderr, "Conversion error: expected %d,%d,%d,%d got %d,%d,%d,%d\n",
						pP->r,pP->g,pP->b,pP->a, rgbcheck.r,rgbcheck.g,rgbcheck.b,rgbcheck.a);
				return -1;
			}
#endif
			/* set all completely transparent colors to black */
			if (!pP->a) {
				*pP = rgb_pixel(pP->a,0,0,0);
			}
		}
	}

	return min_opaque_val;
}

pngquant_error read_image(const char* filename, bool using_stdin, read_info* input_image_p)
{
	FILE* infile;

	if (using_stdin) {
		set_binary_mode(stdin);
		infile = stdin;
	}else if ((infile = fopen(filename, "rb")) == NULL) {
		fprintf(stderr, "  error:  cannot open %s for reading\n", filename);
		return READ_ERROR;
	}

	/*
	 ** Step 1: read in the alpha-channel image.
	 */
	/* GRR:	 returns RGBA (4 channels), 8 bps */
	pngquant_error retval = rwpng_read_image(infile, input_image_p);

	if (!using_stdin)
		fclose(infile);

	if (retval) {
		fprintf(stderr, "  rwpng_read_image() error\n");
		return retval;
	}

	return SUCCESS;
}

/**
 Builds two maps:
	noise - approximation of areas with high-frequency noise, except straight edges. 1=flat, 0=noisy.
	edges - noise map including all edges
 */
static
void contrast_maps(const f_pixel* pixels, size_t cols, size_t rows, double* noise, double* edges)
{
	std::vector<double> tmpVec(cols*rows);
	double* tmp = &tmpVec[0];
	
	const f_pixel* pSrc = pixels;
	double* pNoise = noise;
	double* pEdges = edges;
	for (size_t y=0; y<rows; ++y) {
		f_pixel prev, curr = pSrc[0], next=curr;
		const f_pixel* nextLine = (y == 0) ? pSrc : (pSrc - cols);
		const f_pixel* prevLine = (y == rows-1) ? pSrc : (pSrc + cols);
		for (size_t x=0; x<cols; ++x) {
			prev = curr;
			curr = next;
			next = pSrc[min(cols-1,x+1)];

			f_pixel hd = (prev + next - curr * 2.0).abs();
			f_pixel nextl = nextLine[x];
			f_pixel prevl = prevLine[x];
			f_pixel vd = (prevl + nextl - curr * 2.0).abs();
			double horiz = max(hd.alpha, hd.r, hd.g, hd.b);
			double vert = max(vd.alpha, vd.r, vd.g, vd.b);
			double edge = max(horiz, vert);
			double z = edge - fabs(horiz-vert)*.5;
			z = 1.0 - max(z,min(horiz,vert));
			z *= z;
			z *= z;

			pNoise[x] = z;
			pEdges[x] = 1.0 - edge;
		}
		pSrc += cols;
		pNoise += cols;
		pEdges += cols;
	}
	
	max3(noise, tmp, cols, rows);
	max3(tmp, noise, cols, rows);

	blur(noise, tmp, noise, cols, rows, 3);

	max3(noise, tmp, cols, rows);

	min3(tmp, noise, cols, rows);
	min3(noise, tmp, cols, rows);
	min3(tmp, noise, cols, rows);

	min3(edges, tmp, cols, rows);
	max3(tmp, edges, cols, rows);

	for (size_t i=0; i<cols*rows; ++i) {
		edges[i] = min(noise[i], edges[i]);
	}
}

/**
 * Builds map of neighbor pixels mapped to the same palette entry
 *
 * For efficiency/simplicity it mainly looks for same consecutive pixels horizontally
 * and peeks 1 pixel above/below. Full 2d algorithm doesn't improve it significantly.
 * Correct flood fill doesn't have visually good properties.
 */
void update_dither_map(write_info *output_image, double* edges)
{
	const int width = output_image->width;
	const int height = output_image->height;
	const unsigned char* pixels = output_image->indexed_data;

	for (int row=0; row<height; ++row) {
		unsigned char lastpixel = pixels[row*width];
		int lastcol = 0;
		double* pEdge = &edges[row*width];
		const unsigned char* pLine = &pixels[row*width];
		const unsigned char* pPrevLine;
		const unsigned char* pNextLine;
		if (row > 0) {
			pPrevLine = &pixels[(row-1)*width];
		}
		if (row < height-1) {
			pNextLine = &pixels[(row+1)*width];
		}
		for (int col=1; col<width; ++col) {
			unsigned char px = pLine[col];
			if (px != lastpixel || col == width-1) {
				double neighbor_count = 3.0 + col - lastcol;
				for (int i=lastcol; i<col; ++i) {
					if (row > 0) {
						if (pPrevLine[i] == lastpixel) neighbor_count += 1.0;
					}
					if (row < height-1) {
						if (pNextLine[i] == lastpixel) neighbor_count += 1.0;
					}
				}
				while (lastcol < col) {
					pEdge[lastcol++] *= 1.0 - 3.0/neighbor_count;
				}
				lastpixel = px;
			}
		}
	}
}

void convert(const rgb_pixel*const apixels[], size_t cols, size_t rows, double gamma, f_pixel* dest)
{
	f_pixel* pDst = dest;
	for (size_t y=0; y<rows; ++y) {
		const rgb_pixel* pSrc = apixels[y];
		for (size_t x=0; x<cols; ++x) {
			f_pixel lab = rgb2lab(to_f(gamma, pSrc[x]));
			pDst[x] = lab;
//			pDst[x] = to_f(gamma, pSrc[x]);
		}
		pDst += cols;
	}
}

pngquant_error pngquant(
	read_info* input_image, write_info* output_image,
	bool floyd, int reqcolors, int speed_tradeoff
	)
{
	verbose_printf("  reading file corrected for gamma %2.1f\n", 1.0/input_image->gamma);

	double min_opaque_val = modify_alpha(input_image);
	assert(min_opaque_val > 0);
	size_t width = input_image->width;
	size_t height = input_image->height;
	std::vector<f_pixel> input(width * height);
	convert(
		(const rgb_pixel**)input_image->row_pointers,
		width, height,
		input_image->gamma, &input[0]
		);
	
	std::vector<double> noise(width * height);
	std::vector<double> edges(width * height);
	if (speed_tradeoff < 8) {
		contrast_maps(
			&input[0],
			width, height,
			&noise[0], &edges[0]
			);
	}

	// histogram uses noise contrast map for importance. Color accuracy in noisy areas is not very important.
	// noise map does not include edges to avoid ruining anti-aliasing
	std::vector<hist_item> hist = histogram(&input[0], width, height, reqcolors, speed_tradeoff, &noise[0]);
	std::vector<colormap_item> acolormap;
	double least_error = -1;
	int feedback_loop_trials = 56 - 9*speed_tradeoff;
	const double percent = (double)(feedback_loop_trials>0?feedback_loop_trials:1)/100.0;

	do {
		verbose_printf("  selecting colors");

		std::vector<colormap_item> newmap = mediancut(hist, min_opaque_val, reqcolors);
		
		verbose_printf("...");

		double total_error = 0;
		std::vector<f_pixel> average_color(newmap.size());
		std::vector<f_pixel> base_color(newmap.size());
		std::vector<double> average_color_count(newmap.size());
		std::vector<double> base_color_count(newmap.size());

		if (feedback_loop_trials) {

			viter_init(newmap, &average_color[0], &average_color_count[0], &base_color[0], &base_color_count[0]);

			for (size_t i=0; i<hist.size(); i++) {
				double diff;
				hist_item& hi = hist[i];
				int match = best_color_index(newmap, hi.acolor, &diff);
				assert(diff >= 0);
				assert(hi.perceptual_weight > 0);
				total_error += diff * hi.perceptual_weight;
				
				viter_update_color(hi.acolor, hi.perceptual_weight, newmap, match,
								   &average_color[0], &average_color_count[0], &base_color[0], &base_color_count[0]);

				hi.adjusted_weight = (hi.perceptual_weight+hi.adjusted_weight) * (1.0+sqrt(diff));
			}
		}
		
		if (total_error < least_error || acolormap.size() == 0) {
			acolormap = newmap;

			viter_finalize(acolormap, &average_color[0], &average_color_count[0]);

			least_error = total_error;
			feedback_loop_trials -= 1; // asymptotic improvement could make it go on forever
		}else {
			feedback_loop_trials -= 6;
			if (total_error > least_error*4) feedback_loop_trials -= 3;
		}

		verbose_printf("%d%%\n",100-max(0,(int)(feedback_loop_trials/percent)));
	}while (feedback_loop_trials > 0);

	verbose_printf("  moving colormap towards local minimum\n");

	int iterations = max(5-speed_tradeoff, 0);
	iterations *= iterations;
	for (int i=0; i<iterations; i++) {
		viter_do_interation(hist, acolormap, min_opaque_val);
	}

	output_image->width = input_image->width;
	output_image->height = input_image->height;
	output_image->gamma = 0.45455;
	
	/*
	** Step 3.7 [GRR]: allocate memory for the entire indexed image
	*/
	
	output_image->indexed_data = (unsigned char*) malloc(output_image->height * output_image->width);
	output_image->row_pointers = (unsigned char**) malloc(output_image->height * sizeof(output_image->row_pointers[0]));
	
	if (!output_image->indexed_data || !output_image->row_pointers) {
		fprintf(stderr, "  insufficient memory for indexed data and/or row pointers\n");
		return OUT_OF_MEMORY_ERROR;
	}
	
	for (size_t row=0; row<output_image->height; ++row) {
		output_image->row_pointers[row] = output_image->indexed_data + row*output_image->width;
	}
	
	// tRNS, etc.
	sort_palette(output_image, acolormap);

	/*
	 ** Step 4: map the colors in the image to their closest match in the
	 ** new colormap, and write 'em out.
	 */
	verbose_printf("  mapping image to new colors...");

	// If no dithering is required, that's the final remapping.
	// If dithering (with dither map) is required, this image is used to find areas that require dithering
	double remapping_error = remap_to_palette(&input[0], width, height, output_image, acolormap, min_opaque_val);
	update_dither_map(output_image, &edges[0]);
	
	// remapping error from dithered image is absurd, so always non-dithered value is used
	verbose_printf("MSE=%.3f", remapping_error*256.0);
	
	// remapping above was the last chance to do voronoi iteration, hence the final palette is set after remapping
	set_palette(output_image, acolormap);
	
	if (floyd) {
		remap_to_palette_floyd(&input[0], width, height, output_image, acolormap, min_opaque_val, &edges[0]);
	}

	verbose_printf("\n");
	
	return SUCCESS;
}

