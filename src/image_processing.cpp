//
// image_processing.cpp
//     General image processing routines not directly related to
//     the symmetry transform.
//
// Author:  Christoph Dalitz, Jens Wilberg
// Date:    2019-07-09
// License: see ../LICENSE
//

#include <limits>
#include <math.h>
#include <stdlib.h>
#include <stdexcept>
#include "image.h"

// scale an image down by averaging all values mapped onto the new pixel
// *in* and *out* must be two different objects
void down_scale_image(const GrayImage& in, GrayImage* out, double scale)
{
  size_t height = in.height();
  size_t width = in.width();
  // scale in image
  if (scale < 1) {
    double newval;
    int r = 3;
    size_t ux, vy;
    size_t xx, yy, wcount, factor = 1 / scale;

    // compute new height and width
    size_t newwidth = width * scale;
    size_t newheight = height * scale;

    // compute new radius for the meean filter kernel (default: 3)
    if (factor <= 4) {
      r = 1;
    } else if (factor <= 7) {
      r = 2;
    }
    wcount = (2 * r + 1) * (2 * r + 1);  // number of pixels in the kernel

    out->initialize(newwidth, newheight);

    for (size_t x = 0; x < newwidth; ++x) {
      for (size_t y = 0; y < newheight; ++y) {
        newval = 0.0;
        // compute corresponding pixels in the original image
        yy = y / scale;
        xx = x / scale;
        // iterate over window
        for (int u = -r; u <= r; ++u) {
          for (int v = -r; v <= r; ++v) {
            ux = abs((int)xx + u);
            vy = abs((int)yy + v);
            if (ux >= width) {
              ux = 2 * width - ux - 1;
            }
            if (vy >= height) {
              vy = 2 * height - vy - 1;
            }

            newval += (double)in.get(ux, vy);
          }
        }
        newval /= wcount;
        out->set(x, y, (GrayPixel)newval);
      }
    }
  } else {
    // when the scaling factor is 1, simply copy the original image
    out->initialize(width, height);
    GrayPixel newval;
    for (size_t x = 0; x < in.width(); ++x) {
      for (size_t y = 0; y < in.height(); ++y) {
        newval = in.get(x, y);
        out->set(x, y, newval);
      }
    }
  }
}

// Sobel filter for computing gradient image
void sobel_x(const FloatImage& in, FloatImage* grad_x)
{
  //            -1  0  1
  // kernel =   -2  0  2
  //            -1  0  1
  double sum;
  size_t width = in.width(), height = in.height();
  size_t endx = width-1, endy = height-1;
  grad_x->initialize(width, height);

  // set borders like opencv Sobel
  for (size_t x = 1; x < endx; ++x) {
    sum = -1.0 * in.get(x - 1, 1) - 2.0 * in.get(x - 1, 0) - in.get(x - 1, 1);
    sum += in.get(x + 1, 1) + 2.0 * in.get(x + 1, 0) + in.get(x + 1, 1);
    grad_x->set(x, 0, sum);
    sum = -1.0 * in.get(x - 1, endy - 1) - 2.0 * in.get(x - 1, endy) -
          in.get(x - 1, endy - 1);
    sum += in.get(x + 1, endy - 1) + 2.0 * in.get(x + 1, endy) +
           in.get(x + 1, endy - 1);
    grad_x->set(x, endy, sum);
  }
  for (size_t y = 0; y < height; ++y) {
    grad_x->set(0, y, 0);
    grad_x->set(endx, y, 0);
  }

  for(size_t x = 1; x+1 < in.width(); ++x){
    for(size_t y = 1; y+1 < in.height(); ++y){
      sum = -1.0 * in.get(x-1, y-1) - 2.0*in.get(x-1, y) - in.get(x-1, y+1);
      sum += in.get(x+1, y-1) + 2.0*in.get(x+1, y) + in.get(x+1, y+1);
      grad_x->set(x, y, sum);
    }
  }
}

void sobel_y(const FloatImage& in, FloatImage* grad_y)
{
  //            -1  -2  -1
  // kernel =    0   0   0
  //             1   2   1
  double sum;
  size_t width = in.width(), height = in.height();
  size_t endx = width-1, endy = height-1;
  grad_y->initialize(width, height);

  // set borders like opencv Sobel
  for (size_t x = 0; x < width; ++x) {
    grad_y->set(x, 0, 0);
    grad_y->set(x, endy, 0);
  }
  for (size_t y = 1; y < endy; ++y) {
    sum = -1.0 * in.get(1, y - 1) - 2.0 * in.get(0, y - 1) - in.get(1, y - 1);
    sum += in.get(1, y + 1) + 2.0 * in.get(0, y + 1) + in.get(1, y + 1);
    grad_y->set(0, y, sum);
    sum = -1.0 * in.get(endx - 1, y - 1) - 2.0 * in.get(endx, y - 1) -
          in.get(endx - 1, y - 1);
    sum += in.get(endx - 1, y + 1) + 2.0 * in.get(endx, y + 1) +
           in.get(endx - 1, y + 1);
    grad_y->set(endx, y, sum);
  }

  for(size_t x = 1; x+1 < in.width(); ++x){
    for(size_t y = 1; y+1 < in.height(); ++y){
      sum = -1.0 * in.get(x-1, y-1) - 2.0*in.get(x, y-1) - in.get(x+1, y-1);
      sum += in.get(x-1, y+1) + 2.0*in.get(x, y+1) + in.get(x+1, y+1);
      grad_y->set(x,y, sum);
    }
  }
}

// compute absolute values (Euclidean norm) for each point of vector image
// throws invalid_argument when in_x and in_y size differ
void l2_norm(const FloatImage& in_x, const FloatImage& in_y, FloatImage* abs_val)
{
  size_t height = in_x.height(), width = in_x.width();
  double val;
  if(height != in_y.height() && width != in_y.width()){
    throw std::invalid_argument("in_x and in_y size differ!");
  }
  if(height != abs_val->height() && width != abs_val->width()){
    abs_val->initialize(width, height);
  }

  for(size_t x = 0; x < width; ++x){
    for(size_t y = 0; y < height; ++y){
      val = in_x.get(x, y) * in_x.get(x, y) + in_y.get(x, y) * in_y.get(x, y);
      val = sqrt(val);
      abs_val->set(x, y, val);
    }
  }
}
