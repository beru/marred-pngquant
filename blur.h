#pragma once

//
//  blur.h
//  pngquant

void blur(const double* src, double* tmp, double* dst, int width, int height, const int size);
void max3(const double* src, double* dst, int width, int height);
void min3(const double* src, double* dst, int width, int height);

