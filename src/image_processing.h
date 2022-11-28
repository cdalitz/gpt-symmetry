//
// image_processing.h
//     General image processing routines not directly related to
//     the symmetry transform.
//
// Author:  Christoph Dalitz, Jens Wilberg
// Date:    2019-07-09
// License: see ../LICENSE
//

#ifndef __IMAGE_PROCESSING_H
#define __IMAGE_PROCESSING_H

#include <stddef.h>
#include "image.h"

// scale an image down with pixel averaging
void down_scale_image(const GrayImage& in, GrayImage* out, double scale);

// Sobel filter for computing gradient image
void sobel_x(const FloatImage& in, FloatImage* grad_x);
void sobel_y(const FloatImage& in, FloatImage* grad_y);

// compute absolute values (Euclidean norm) for each point of vector image
// throws invalid_argument when in_x and in_y size differ
void l2_norm(const FloatImage& in_x, const FloatImage& in_y, FloatImage* abs_val);

#endif
