#pragma once

//
//  blur.h
//  pngquant

void blur(double* src, double* tmp, double* dst, int width, int height, int size);
void max3(double* src,double* dst, int width, int height);
void min3(double* src,double* dst, int width, int height);

