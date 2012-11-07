#pragma once

//
//  blur.h
//  pngquant

void blur(
	const double* src,
	double* tmp,
	double* dst,
	unsigned int width,
	unsigned int height,
	unsigned int size
);

void max3(
	const double* src,
	double* dst,
	unsigned int width,
	unsigned int height
);

void min3(
	const double* src,
	double* dst,
	unsigned int width,
	unsigned int height
);

