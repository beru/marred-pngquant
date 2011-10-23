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
	  -iebug		increase opacity to work around Internet Explorer 6 bug\n\
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

#include <vector>
#include <algorithm>

pngquant_error pngquant(read_info* input_image, write_info* output_image, bool floyd, int reqcolors, bool ie_bug, int speed_tradeoff);
pngquant_error read_image(const char* filename, bool using_stdin, read_info* input_image_p);
pngquant_error write_image(write_info* output_image, const char* filename, const char* newext, bool force, bool using_stdin);

static void viter_init(const colormap* map, f_pixel* average_color, float* average_color_count, f_pixel* base_color, float* base_color_count);
static void viter_update_color(f_pixel acolor, float value, colormap* map, int match,
							   f_pixel* average_color, float* average_color_count,
							   const f_pixel* base_color, const float* base_color_count);
static void viter_finalize(colormap* map, const f_pixel* average_color, const float* average_color_count);

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
	fprintf(fd, "origina pngquant, by porneL, Greg Roelofs, Kornel Lesinski.\n");
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
	bool ie_bug = false;
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
		else if (0 == strcmp(argv[argn], "-iebug"))
			ie_bug = true;
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
		newext = floyd ? "-ie-fs8.png" : "-ie-or8.png";
		if (!ie_bug) newext += 3; /* skip "-ie" */
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
			retval = pngquant(&input_image, &output_image, floyd, reqcolors, ie_bug, speed_tradeoff);
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
	return v1.popularity - v2.popularity;
}

void sort_palette(write_info* output_image, colormap* map)
{
	assert(map); assert(output_image);

	/*
	** Step 3.5 [GRR]: remap the palette colors so that all entries with
	** the maximal alpha value (i.e., fully opaque) are at the end and can
	** therefore be omitted from the tRNS chunk.
	*/

	verbose_printf("  eliminating opaque tRNS-chunk entries...");

	/* move transparent colors to the beginning to shrink trns chunk */
	int num_transparent=0;
	for (int i=0; i<map->colors; i++) {
		colormap_item& pal = map->palette[i];
		rgb_pixel px = to_rgb(output_image->gamma, pal.acolor);
		if (px.a != 255) {
			if (i != num_transparent) {
				std::swap(map->palette[num_transparent], pal);
				i--;
			}
			num_transparent++;
		}
	}

	verbose_printf("%d entr%s transparent\n", num_transparent, (num_transparent == 1)? "y" : "ies");

	/* colors sorted by popularity make pngs slightly more compressible
	 * opaque and transparent are sorted separately
	 */
	std::sort(map->palette, map->palette+num_transparent, compare_popularity);
	std::sort(map->palette+num_transparent, map->palette+map->colors-num_transparent, compare_popularity);

	output_image->num_palette = map->colors;
	output_image->num_trans = num_transparent;
}

void set_palette(write_info* output_image, colormap* map)
{
	for (int x=0; x<map->colors; ++x) {
		colormap_item& pal = map->palette[x];
		rgb_pixel px = to_rgb(output_image->gamma, pal.acolor);
		pal.acolor = to_f(output_image->gamma, px); /* saves rounding error introduced by to_rgb, which makes remapping & dithering more accurate */

		png_color& pc = output_image->palette[x];
		pc.red	 = px.r;
		pc.green = px.g;
		pc.blue	 = px.b;
		output_image->trans[x] = px.a;
	}
}

static
int best_color_index(f_pixel px, const colormap* map, float min_opaque_val, float* dist_out)
{
	const colormap_item* const acolormap = map->palette;
	const int numcolors = map->colors;
	int ind = 0;
	const int iebug = px.a > min_opaque_val;
	float dist = colordifference(px, acolormap[0].acolor);

	for (int i=1; i<numcolors; i++) {
		f_pixel itemCol = acolormap[i].acolor;
		float newdist = colordifference(px, itemCol);
		if (newdist < dist) {
			/* penalty for making holes in IE */
			if (iebug && itemCol.a < 1) {
				if (newdist+1.f/1024.f > dist) continue;
			}
			ind = i;
			dist = newdist;
		}
	}

	if (dist_out) *dist_out = dist;
	return ind;
}

float remap_to_palette(read_info* input_image, write_info* output_image, colormap* map, float min_opaque_val, bool ie_bug)
{
	rgb_pixel** input_pixels = (rgb_pixel**) input_image->row_pointers;
	unsigned char** row_pointers = output_image->row_pointers;
	int rows = input_image->height, cols = input_image->width;
	double gamma = input_image->gamma;

	int remapped_pixels=0;
	float remapping_error=0;
	
	int transparent_ind = best_color_index(f_pixel(0,0,0,0), map, min_opaque_val, NULL);

	std::vector<f_pixel> average_color(map->colors);
	std::vector<float> average_color_count(map->colors);
	viter_init(map, &average_color[0], &average_color_count[0], NULL, NULL);

	for (int row=0; row<rows; ++row) {
		for (int col=0; col<cols; ++col) {

			f_pixel px = to_f(gamma, input_pixels[row][col]);
			int match;

			if (px.a < 1.0/256.0) {
				match = transparent_ind;
			}else {
				float diff;
				match = best_color_index(px, map,min_opaque_val, &diff);

				remapped_pixels++;
				remapping_error += diff;
			}

			row_pointers[row][col] = match;

			viter_update_color(px, 1.0, map, match, &average_color[0], &average_color_count[0], NULL, NULL);
		}
	}

	viter_finalize(map, &average_color[0], &average_color_count[0]);

	return remapping_error / MAX(1,remapped_pixels);
}

float remap_to_palette_floyd(read_info* input_image, write_info* output_image, const colormap* map, float min_opaque_val, bool ie_bug)
{
	rgb_pixel** input_pixels = (rgb_pixel**) input_image->row_pointers;
	unsigned char** row_pointers = output_image->row_pointers;
	int rows = input_image->height, cols = input_image->width;
	double gamma = input_image->gamma;

	int remapped_pixels = 0;
	float remapping_error = 0;

	const colormap_item* acolormap = map->palette;

	int ind = 0;
	int transparent_ind = best_color_index(f_pixel(0,0,0,0), map, min_opaque_val, NULL);

	f_pixel* thiserr = NULL;
	f_pixel* nexterr = NULL;
	float sr=0, sg=0, sb=0, sa=0;
	int fs_direction = 1;

	/* Initialize Floyd-Steinberg error vectors. */
	thiserr = (f_pixel*) malloc((cols + 2) * sizeof(*thiserr));
	nexterr = (f_pixel*) malloc((cols + 2) * sizeof(*thiserr));
	srand(12345); /* deterministic dithering is better for comparing results */

	for (int col=0; col<cols+2; ++col) {
		const double rand_max = RAND_MAX;
		thiserr[col].r = ((double)rand() - rand_max/2.0)/rand_max/255.0;
		thiserr[col].g = ((double)rand() - rand_max/2.0)/rand_max/255.0;
		thiserr[col].b = ((double)rand() - rand_max/2.0)/rand_max/255.0;
		thiserr[col].a = ((double)rand() - rand_max/2.0)/rand_max/255.0;
	}

	for (int row=0; row<rows; ++row) {
		memset(nexterr, 0, (cols + 2) * sizeof(*nexterr));

		int col = (fs_direction) ? 0 : (cols - 1);

		rgb_pixel* input_row = input_pixels[row];
		do {
			f_pixel px = to_f(gamma, input_row[col]);

			/* Use Floyd-Steinberg errors to adjust actual color. */
			sr = px.r + thiserr[col + 1].r;
			sg = px.g + thiserr[col + 1].g;
			sb = px.b + thiserr[col + 1].b;
			sa = px.a + thiserr[col + 1].a;

			if (sr < 0) sr = 0;
			else if (sr > 1) sr = 1;
			if (sg < 0) sg = 0;
			else if (sg > 1) sg = 1;
			if (sb < 0) sb = 0;
			else if (sb > 1) sb = 1;
			if (sa < 0) sa = 0;
			/* when fighting IE bug, dithering must not make opaque areas transparent */
			else if (sa > 1 || (ie_bug && px.a > 255.0/256.0)) sa = 1;

			if (sa < 1.0/256.0) {
				ind = transparent_ind;
			}else {
				float diff;
				ind = best_color_index(f_pixel(sa,sr,sg,sb), map, min_opaque_val, &diff);

				remapped_pixels++;
				remapping_error += diff;
			}

			row_pointers[row][col] = ind;

			float colorimp = (3.0f + acolormap[ind].acolor.a) / 4.0f;
			f_pixel xp = acolormap[ind].acolor;

			f_pixel err;
			err.r = (sr - xp.r) * colorimp;
			err.g = (sg - xp.g) * colorimp;
			err.b = (sb - xp.b) * colorimp;
			err.a = (sa - xp.a);

			/* Propagate Floyd-Steinberg error terms. */
			if (fs_direction) {
				thiserr[col + 2].a += (err.a * 7.0f) / 16.0f;
				thiserr[col + 2].r += (err.r * 7.0f) / 16.0f;
				thiserr[col + 2].g += (err.g * 7.0f) / 16.0f;
				thiserr[col + 2].b += (err.b * 7.0f) / 16.0f;

				nexterr[col	   ].a += (err.a * 3.0f) / 16.0f;
				nexterr[col	   ].r += (err.r * 3.0f) / 16.0f;
				nexterr[col	   ].g += (err.g * 3.0f) / 16.0f;
				nexterr[col	   ].b += (err.b * 3.0f) / 16.0f;

				nexterr[col + 1].a += (err.a * 5.0f) / 16.0f;
				nexterr[col + 1].r += (err.r * 5.0f) / 16.0f;
				nexterr[col + 1].g += (err.g * 5.0f) / 16.0f;
				nexterr[col + 1].b += (err.b * 5.0f) / 16.0f;

				nexterr[col + 2].a += (err.a	   ) / 16.0f;
				nexterr[col + 2].r += (err.r	   ) / 16.0f;
				nexterr[col + 2].g += (err.g	   ) / 16.0f;
				nexterr[col + 2].b += (err.b	   ) / 16.0f;
			}else {
				thiserr[col	   ].a += (err.a * 7.0f) / 16.0f;
				thiserr[col	   ].r += (err.r * 7.0f) / 16.0f;
				thiserr[col	   ].g += (err.g * 7.0f) / 16.0f;
				thiserr[col	   ].b += (err.b * 7.0f) / 16.0f;

				nexterr[col	   ].a += (err.a	   ) / 16.0f;
				nexterr[col	   ].r += (err.r	   ) / 16.0f;
				nexterr[col	   ].g += (err.g	   ) / 16.0f;
				nexterr[col	   ].b += (err.b	   ) / 16.0f;

				nexterr[col + 1].a += (err.a * 5.0f) / 16.0f;
				nexterr[col + 1].r += (err.r * 5.0f) / 16.0f;
				nexterr[col + 1].g += (err.g * 5.0f) / 16.0f;
				nexterr[col + 1].b += (err.b * 5.0f) / 16.0f;

				nexterr[col + 2].a += (err.a * 3.0f) / 16.0f;
				nexterr[col + 2].r += (err.r * 3.0f) / 16.0f;
				nexterr[col + 2].g += (err.g * 3.0f) / 16.0f;
				nexterr[col + 2].b += (err.b * 3.0f) / 16.0f;
			}

			if (fs_direction) {
				++col;
				if (col >= cols) break;
			}else {
				--col;
				if (col < 0) break;
			}
		}while(1);
		std::swap(thiserr, nexterr);
		fs_direction = !fs_direction;
	}

	return remapping_error / MAX(1, remapped_pixels);
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
	setmode(fp == stdout ? 1 : 0, O_BINARY);
#endif
#endif
}

pngquant_error write_image(write_info* output_image, const char* filename, const char* newext, bool force, bool using_stdin)
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

hist* histogram(read_info* input_image, int reqcolors, int speed_tradeoff)
{
	hist* hist;
	int ignorebits = 0;
	const rgb_pixel*const* input_pixels = (const rgb_pixel*const*)input_image->row_pointers;
	int cols = input_image->width, rows = input_image->height;
	double gamma = input_image->gamma;
	assert(gamma > 0);

   /*
	** Step 2: attempt to make a histogram of the colors, unclustered.
	** If at first we don't succeed, increase ignorebits to increase color
	** coherence and try again.
	*/

	if (speed_tradeoff > 7) ignorebits++;
	int maxcolors = (1<<17) + (1<<18)*(10-speed_tradeoff);

	verbose_printf("  making histogram...");
	for (;;) {

		hist = pam_computeacolorhist(input_pixels, cols, rows, gamma, maxcolors, ignorebits, speed_tradeoff < 9);
		if (hist) break;

		ignorebits++;
		verbose_printf("too many colors!\n	scaling colors to improve clustering...");
	}

	verbose_printf("%d colors found\n", hist->size);
	return hist;
}

float modify_alpha(read_info* input_image, bool ie_bug)
{
	/* IE6 makes colors with even slightest transparency completely transparent,
	   thus to improve situation in IE, make colors that are less than ~10% transparent
	   completely opaque */

	rgb_pixel** input_pixels = (rgb_pixel**) input_image->row_pointers;
	rgb_pixel* pP;
	const int rows = input_image->height;
	const int cols = input_image->width;
	double gamma = input_image->gamma;
	float min_opaque_val, almost_opaque_val;

	if (ie_bug) {
		min_opaque_val = 238.0/256.0; /* rest of the code uses min_opaque_val rather than checking for ie_bug */
		almost_opaque_val = min_opaque_val * 169.0/256.0;

		verbose_printf("  Working around IE6 bug by making image less transparent...\n");
	}else {
		min_opaque_val = almost_opaque_val = 1;
	}

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
			/* ie bug: to avoid visible step caused by forced opaqueness, linearily raise opaqueness of almost-opaque colors */
			else if (pP->a < 255 && px.a > almost_opaque_val) {
				assert((min_opaque_val-almost_opaque_val)>0);

				float al = almost_opaque_val + (px.a-almost_opaque_val) * (1-almost_opaque_val) / (min_opaque_val-almost_opaque_val);
				if (al > 1) al = 1;
				px.a = al;
				pP->a = to_rgb(gamma, px).a;
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

/*
 * Voronoi iteration: new palette color is computed from weighted average of colors that map to that palette entry.
 */
static
void viter_init(const colormap* map,
				 f_pixel* average_color, float* average_color_count,
				 f_pixel* base_color, float* base_color_count)
{
	colormap_item* newmap = map->palette;
	const int newcolors = map->colors;
	for (int i=0; i<newcolors; i++) {
		average_color_count[i] = 0;
		average_color[i] = f_pixel(0,0,0,0);
	}

	// Rather than only using separate mapping and averaging steps
	// new palette colors are computed at the same time as mapping is done
	// but to avoid first few matches moving the entry too much
	// some base color and weight is added
	if (base_color) {
		for (int i=0; i<newcolors; i++) {
			const colormap_item& nmi = newmap[i];
			float value = 1.0 + nmi.popularity/2.0;
			base_color_count[i] = value;
			const f_pixel& nc = nmi.acolor;
			base_color[i] = f_pixel(
				nc.a * value,
				nc.r * value,
				nc.g * value,
				nc.b * value
			);
		}
	}
}

static
void viter_update_color(f_pixel acolor, float value, colormap* map, int match,
						 f_pixel* average_color, float* average_color_count,
						 const f_pixel* base_color, const float* base_color_count)
{
	f_pixel& ac = average_color[match];
	float& acc = average_color_count[match];
	const f_pixel& bc = base_color[match];
	const float& bcc = base_color_count[match];
	ac.a += acolor.a * value;
	ac.r += acolor.r * value;
	ac.g += acolor.g * value;
	ac.b += acolor.b * value;
	acc += value;

	if (base_color) {
		float tcc = acc + bcc;
		map->palette[match].acolor = f_pixel(
			(ac.a + bc.a) / tcc,
			(ac.r + bc.r) / tcc,
			(ac.g + bc.g) / tcc,
			(ac.b + bc.b) / tcc
		);
	}
}

static
void viter_finalize(colormap* map, const f_pixel* average_color, const float* average_color_count)
{
	for (int i=0; i<map->colors; i++) {
		colormap_item& pal = map->palette[i];
		const float acc = average_color_count[i];
		if (acc) {
			const f_pixel& ac = average_color[i];
			pal.acolor = f_pixel(
				ac.a / acc,
				ac.r / acc,
				ac.g / acc,
				ac.b / acc
			);
		}
		pal.popularity = acc;
	}
}

void viter_do_interation(const hist* hist, colormap* map, float min_opaque_val)
{
	std::vector<f_pixel> average_color(map->colors);
	std::vector<float> average_color_count(map->colors);

	hist_item* achv = hist->achv;
	viter_init(map, &average_color[0], &average_color_count[0], NULL,NULL);

	for (int j=0; j<hist->size; j++) {
		const hist_item& hist = achv[j];
		int match = best_color_index(hist.acolor, map, min_opaque_val, NULL);
		viter_update_color(hist.acolor, hist.perceptual_weight, map, match, &average_color[0], &average_color_count[0], NULL, NULL);
	}

	viter_finalize(map, &average_color[0], &average_color_count[0]);
}

pngquant_error pngquant(read_info* input_image, write_info* output_image, bool floyd, int reqcolors, bool ie_bug, int speed_tradeoff)
{
	float min_opaque_val;

	verbose_printf("  reading file corrected for gamma %2.1f\n", 1.0/input_image->gamma);

	min_opaque_val = modify_alpha(input_image, ie_bug);
	assert(min_opaque_val>0);

	hist* hist = histogram(input_image, reqcolors, speed_tradeoff);
	hist_item* achv = hist->achv;

	colormap* acolormap = NULL;
	float least_error = -1;
	int feedback_loop_trials = 56 - 9*speed_tradeoff;
	const double percent = (double)(feedback_loop_trials>0?feedback_loop_trials:1)/100.0;

	do {
		verbose_printf("  selecting colors");

		colormap *newmap = mediancut(hist, min_opaque_val, reqcolors);

		verbose_printf("...");

		float total_error = 0;
		std::vector<f_pixel> average_color(newmap->colors);
		std::vector<f_pixel> base_color(newmap->colors);
		std::vector<float> average_color_count(newmap->colors);
		std::vector<float> base_color_count(newmap->colors);

		if (feedback_loop_trials) {

			viter_init(newmap, &average_color[0], &average_color_count[0], &base_color[0], &base_color_count[0]);

			for (int i=0; i<hist->size; i++) {
				float diff;
				int match = best_color_index(achv[i].acolor, newmap, min_opaque_val, &diff);
				assert(diff >= 0);
				assert(achv[i].perceptual_weight > 0);
				total_error += diff * achv[i].perceptual_weight;

				viter_update_color(achv[i].acolor, achv[i].perceptual_weight, newmap, match,
								   &average_color[0], &average_color_count[0], &base_color[0], &base_color_count[0]);

				achv[i].adjusted_weight = (achv[i].perceptual_weight+achv[i].adjusted_weight) * (1.0+sqrtf(diff));
			}
		}

		if (total_error < least_error || !acolormap) {
			if (acolormap) pam_freecolormap(acolormap);

			acolormap = newmap;

			viter_finalize(acolormap, &average_color[0], &average_color_count[0]);

			least_error = total_error;
			feedback_loop_trials -= 1; // asymptotic improvement could make it go on forever
		}else {
			feedback_loop_trials -= 6;
			if (total_error > least_error*4) feedback_loop_trials -= 3;
			pam_freecolormap(newmap);
		}

		verbose_printf("%d%%\n",100-MAX(0,(int)(feedback_loop_trials/percent)));
	}while(feedback_loop_trials > 0);

	verbose_printf("  moving colormap towards local minimum\n");

	int iterations = MAX(5-speed_tradeoff, 0);
	iterations *= iterations;
	for (int i=0; i<iterations; i++) {
		viter_do_interation(hist, acolormap, min_opaque_val);
	}

	pam_freeacolorhist(hist);

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
	
	for (int row=0; row<output_image->height; ++row) {
		output_image->row_pointers[row] = output_image->indexed_data + row*output_image->width;
	}

	// tRNS, etc.
	sort_palette(output_image, acolormap);

	/*
	 ** Step 4: map the colors in the image to their closest match in the
	 ** new colormap, and write 'em out.
	 */
	verbose_printf("  mapping image to new colors...");

	float remapping_error;

	if (floyd) {
		// if dithering, save rounding error and stick to that palette
		// otherwise palette can be improved after remapping
		set_palette(output_image, acolormap);
		remapping_error = remap_to_palette_floyd(input_image, output_image, acolormap, min_opaque_val, ie_bug);
	}else {
		remapping_error = remap_to_palette(input_image, output_image, acolormap, min_opaque_val, ie_bug);
		set_palette(output_image, acolormap);
	}
	
	verbose_printf("MSE=%.3f\n", remapping_error*256.0f);
	
	pam_freecolormap(acolormap);
	
	return SUCCESS;
}


