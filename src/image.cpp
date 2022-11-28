
//
// image.cpp
//     Basic image utility functions.
//     File I/O is done with png++, but this is wrapped with own
//     image classes in order to minimize dependencies on 3rd party libs.
//
// Author:  Christoph Dalitz
// Date:    2019-07-09
// License: see ../LICENSE
//

#include <iostream>
#include <string>
#include <math.h>
#include "png++/png.hpp"

#include "image.h"


//---------------------------------------------------------------------
// Grayscale image
//---------------------------------------------------------------------


// Constructor (without allocation)
GrayImage::GrayImage()
{
  this->_width = 0;
  this->_height = 0;
  _data = NULL;
}

// Destructor
GrayImage::~GrayImage()
{
  if (_data) delete[] _data;
}

// allocation (deletes old data)
void GrayImage::initialize(size_t width, size_t height)
{
  if (_data) delete[] _data;
  this->_width = width;
  this->_height = height;
  _data = new GrayPixel[width*height];
}

// allocation and initializing all values (deletes old data)
void GrayImage::initialize(size_t width, size_t height, GrayPixel value)
{
  if (_data) delete[] _data;
  this->_width = width;
  this->_height = height;
  _data = new GrayPixel[width*height];
  for (size_t y = 0; y < height; ++y)
    for (size_t x = 0; x < width; ++x)
      this->set(x, y, value);
}


// read PNG file and convert it to grayscale
// Returncode: 0=success, 1=file error, 2=libpng error
// in case of an error, the error message is printed ot stderr
int GrayImage::read_png(const char* filename)
{
  // read PNG file
  png::image<png::rgb_pixel>* pngimg;
  try {
    pngimg = new png::image<png::rgb_pixel>(filename);
  } catch(png::std_error& e) {
    std::cerr <<  "Error while reading '" << std::string(filename) << "': "
              << e.what() << std::endl;
    return 1;
  } catch(png::error& e) {
    std::cerr <<  "Error while reading '" << std::string(filename) << "': "
              << e.what() << std::endl;
    return 2;
  }
  
  // allocate current image
  initialize((size_t)pngimg->get_width(), (size_t)pngimg->get_height());
  
  // copy over PNG image to curent image
  png::rgb_pixel p;
  png::gray_pixel q;
  for (png::uint_32 y = 0; y < pngimg->get_height(); ++y) {
    for (png::uint_32 x = 0; x < pngimg->get_width(); ++x) {
      p = pngimg->get_pixel(x, y);
      q = (p.red + p.green + p.blue) / 3;
      this->set(x, y, q);
    }
  }
  // clean up
  delete pngimg;

  return 0;
}

// write PNG file
// Returncode: 0=success, 1=file error, 2=libpng error
// in case of an error, the error message is printed ot stderr
int GrayImage::write_png(const char* filename)
{
  png::image< png::gray_pixel > pngimg(this->_width, this->_height);
  for (size_t y = 0; y < this->_height; ++y) {
    for (size_t x = 0; x < this->_width; ++x) {
      pngimg.set_pixel(x, y, this->get(x,y));
    }
  }
  try {
    pngimg.write(filename);
  } catch(png::std_error& e) {
    std::cerr <<  "Error writing '" << std::string(filename) << "': "
              << e.what() << std::endl;
    return 1;
  } catch(png::error& e) {
    std::cerr <<  "Error writing '" << std::string(filename) << "': "
              << e.what() << std::endl;
    return 2;
  }
  return 0;
}

// transform to RGB image
void GrayImage::to_rgb(RGBImage* rgb)
{
  rgb->initialize(_width, _height);
  if (_width == 0 || _height == 0) return;

  GrayPixel tmp;
  for (size_t y = 0; y < _height; ++y) {
    for (size_t x = 0; x < _width; ++x) {
      tmp = this->get(x,y);
      rgb->set(x, y, RGBPixel(tmp,tmp,tmp));
    }
  }
}

// transform to float image
void GrayImage::to_float(FloatImage* fimg)
{
  fimg->initialize(_width, _height);
  if (_width == 0 || _height == 0) return;

  for (size_t y = 0; y < _height; ++y) {
    for (size_t x = 0; x < _width; ++x) {
      fimg->set(x, y, (FloatPixel)this->get(x,y));
    }
  }
}

// bilinear interpolation at real values
GrayPixel GrayImage::interpolate(double x, double y) const
{
  double xn, xnn, yn, ynn, val1, val2, val;
  size_t in, inn, jn, jnn;
  xn = floor(x); xnn = xn + 1;
  yn = floor(y); ynn = yn + 1;
  // border extrapolation
  in = (xn > 0.0 ? (size_t)xn : 0);
  inn = (xnn < _width ? (size_t)xnn : (_width - 1));
  jn = (yn > 0.0 ? (size_t)yn : 0);
  jnn = (ynn < _height ? (size_t)ynn : (_height - 1));
  // bilinear interpolation
  val1 = this->get(in,jn)*(ynn-y)/(ynn-yn) + this->get(in,jnn)*(y-yn)/(ynn-yn);
  val2 = this->get(inn,jn)*(ynn-y)/(ynn-yn) + this->get(inn,jnn)*(y-yn)/(ynn-yn);
  val = val1*(xnn-x)/(xnn-xn) + val2*(x-xn)/(xnn-xn);
  return (GrayPixel)val;
}


//---------------------------------------------------------------------
// Float image
//---------------------------------------------------------------------


// Constructor (without allocation)
FloatImage::FloatImage()
{
  this->_width = 0;
  this->_height = 0;
  _data = NULL;
}

// Destructor
FloatImage::~FloatImage()
{
  if (_data) delete[] _data;
}

// allocation (deletes old data)
void FloatImage::initialize(size_t width, size_t height)
{
  if (_data) delete[] _data;
  this->_width = width;
  this->_height = height;
  _data = new FloatPixel[width*height];
}

// allocation and initializing all values (deletes old data)
void FloatImage::initialize(size_t width, size_t height, FloatPixel value)
{
  if (_data) delete[] _data;
  this->_width = width;
  this->_height = height;
  _data = new FloatPixel[width*height];
  for (size_t y = 0; y < height; ++y)
    for (size_t x = 0; x < width; ++x)
      this->set(x, y, value);
}

// transform to grayscale image
// if onlypositive = true, only the positive range is transformed
// if enhance = true, the grayscale image histogram is equalized
void FloatImage::to_gray(GrayImage* gray, bool onlypositive /*=false*/, bool enhance /*=false*/)
{
  gray->initialize(_width, _height);
  if (_width == 0 || _height == 0) return;
  
  // find range of values
  FloatPixel minval, maxval, tmp;
  minval = maxval = this->get(0,0);
  for (size_t y = 0; y < _height; ++y) {
    for (size_t x = 0; x < _width; ++x) {
      tmp = this->get(x,y);
      if (tmp < minval) minval = tmp;
      if (tmp > maxval) maxval = tmp;
    }
  }
  if (onlypositive && minval<0.0)
    minval = 0.0;
  for (size_t y = 0; y < _height; ++y) {
    for (size_t x = 0; x < _width; ++x) {
      tmp = this->get(x,y);
      if (tmp <= minval)
        gray->set(x, y, 0);
      else
        gray->set(x, y, GrayPixel(255*(tmp - minval)/(maxval-minval)));
    }
  }

  // image equalization
  if (enhance) {
    // build histogramm
    std::vector<size_t> hist(256, 0);
    for (size_t y = 0; y < _height; ++y) {
      for (size_t x = 0; x < _width; ++x) {
        hist[gray->get(x,y)]++;
      }
    }
    size_t pixelcount = _width * _height - hist[0];

    hist[0] = 0;
    for(size_t i = 1; i < 256; ++i){
      hist[i] = hist[i-1] + hist[i];
    }
    for (size_t y = 0; y < _height; ++y) {
      for (size_t x = 0; x < _width; ++x) {
        GrayPixel b = (GrayPixel) (hist[gray->get(x,y)] * 255 / pixelcount);
        gray->set(x,y,b);
      }
    }
  }
}

// return mean pixel value
FloatPixel FloatImage::mean_value() const
{
  FloatPixel res = 0.0;
  for (size_t y = 0; y < _height; ++y)
    for (size_t x = 0; x < _width; ++x)
      res += this->get(x, y);
  return res / (_width*_height);
}

// return min pixel value
FloatPixel FloatImage::min_value() const
{
  FloatPixel res = this->get(0,0);
  for (size_t y = 0; y < _height; ++y)
    for (size_t x = 0; x < _width; ++x)
      if (this->get(x,y) < res)
        res = this->get(x,y);
  return res;
}

// bilinear interpolation at real values
FloatPixel FloatImage::interpolate(double x, double y) const
{
  double xn, xnn, yn, ynn;
  FloatPixel val1, val2, val;
  size_t in, inn, jn, jnn;
  xn = floor(x); xnn = xn + 1;
  yn = floor(y); ynn = yn + 1;
  // border extrapolation
  in = (xn > 0.0 ? (size_t)xn : 0);
  inn = (xnn < _width ? (size_t)xnn : (_width - 1));
  jn = (yn > 0.0 ? (size_t)yn : 0);
  jnn = (ynn < _height ? (size_t)ynn : (_height - 1));
  // bilinear interpolation
  val1 = this->get(in,jn)*(ynn-y)/(ynn-yn) + this->get(in,jnn)*(y-yn)/(ynn-yn);
  val2 = this->get(inn,jn)*(ynn-y)/(ynn-yn) + this->get(inn,jnn)*(y-yn)/(ynn-yn);
  val = val1*(xnn-x)/(xnn-xn) + val2*(x-xn)/(xnn-xn);
  return val;
}


//---------------------------------------------------------------------
// RGB image
//---------------------------------------------------------------------


// Constructor (without allocation)
RGBImage::RGBImage()
{
  this->_width = 0;
  this->_height = 0;
  _data = NULL;
}

// Destructor
RGBImage::~RGBImage()
{
  if (_data) delete[] _data;
}

// allocation (deletes old data)
void RGBImage::initialize(size_t width, size_t height)
{
  if (_data) delete[] _data;
  this->_width = width;
  this->_height = height;
  _data = new RGBPixel[width*height];
}

// allocation and initializing all values (deletes old data)
void RGBImage::initialize(size_t width, size_t height, const RGBPixel& value)
{
  if (_data) delete[] _data;
  this->_width = width;
  this->_height = height;
  _data = new RGBPixel[width*height];
  for (size_t y = 0; y < height; ++y)
    for (size_t x = 0; x < width; ++x)
      this->set(x, y, value);
}

// read PNG file
// Returncode: 0=success, 1=file error, 2=libpng error
// in case of an error, the error message is printed ot stderr
int RGBImage::read_png(const char* filename)
{
  // read PNG file
  png::image<png::rgb_pixel>* pngimg;
  try {
    pngimg = new png::image<png::rgb_pixel>(filename);
  } catch(png::std_error& e) {
    std::cerr <<  "Error while reading '" << std::string(filename) << "': "
              << e.what() << std::endl;
    return 1;
  } catch(png::error& e) {
    std::cerr <<  "Error while reading '" << std::string(filename) << "': "
              << e.what() << std::endl;
    return 2;
  }
  
  // allocate current image
  initialize((size_t)pngimg->get_width(), (size_t)pngimg->get_height());
  
  // copy over PNG image to curent image
  png::rgb_pixel p;
  for (png::uint_32 y = 0; y < pngimg->get_height(); ++y) {
    for (png::uint_32 x = 0; x < pngimg->get_width(); ++x) {
      p = pngimg->get_pixel(x, y);
      this->set(x, y, RGBPixel(p.red, p.green, p.blue));
    }
  }
  // clean up
  delete pngimg;

  return 0;
}

// write PNG file
// Returncode: 0=success, 1=file error, 2=libpng error
// in case of an error, the error message is printed ot stderr
int RGBImage::write_png(const char* filename)
{
  png::image< png::rgb_pixel > pngimg(this->_width, this->_height);
  RGBPixel tmp;
  for (size_t y = 0; y < this->_height; ++y) {
    for (size_t x = 0; x < this->_width; ++x) {
      tmp = this->get(x,y);
      pngimg.set_pixel(x, y, png::rgb_pixel(tmp.R, tmp.G, tmp.B));
    }
  }
  try {
    pngimg.write(filename);
  } catch(png::std_error& e) {
    std::cerr <<  "Error writing '" << std::string(filename) << "': "
              << e.what() << std::endl;
    return 1;
  } catch(png::error& e) {
    std::cerr <<  "Error writing '" << std::string(filename) << "': "
              << e.what() << std::endl;
    return 2;
  }
  return 0;
}

// transform to grayscale image
void RGBImage::to_gray(GrayImage* gray)
{
  gray->initialize(_width, _height);
  if (_width == 0 || _height == 0) return;

  RGBPixel tmp;
  for (size_t y = 0; y < _height; ++y) {
    for (size_t x = 0; x < _width; ++x) {
      tmp = this->get(x,y);
      gray->set(x, y, (tmp.R + tmp.G + tmp.B) / 3);
    }
  }
}

// draw rectangle with corners (ul_x, ul_y) and (lr_x, lr_y)
void RGBImage::draw_rect_thickness1(size_t ul_x, size_t ul_y, size_t lr_x, size_t lr_y, const RGBPixel& value) {
  if (ul_x >= _width || lr_x >= _width || ul_y >= _height || lr_y >= _height)
    return;
  size_t x,y;
  for (x = ul_x; x <= lr_x; x++)
    this->set(x, ul_y, value);
  for (x = ul_x; x <= lr_x; x++)
    this->set(x, lr_y, value);
  for (y = ul_y; y <= lr_y; y++)
    this->set(ul_x, y, value);
  for (y = ul_y; y <= lr_y; y++)
    this->set(lr_x, y, value);
}

// draw rectangle with corners (ul_x, ul_y) and (lr_x, lr_y)
// and given border thickness
void RGBImage::draw_rect(size_t ul_x, size_t ul_y, size_t lr_x, size_t lr_y, const RGBPixel& value, size_t thickness /*=1*/) {

  draw_rect_thickness1(ul_x, ul_y, lr_x, lr_y, value);
  size_t x1, y1, x2, y2;
  x1 = ul_x; y1 = ul_y; x2 = lr_x; y2 = lr_y;
  for (size_t i=1; i<=thickness/2; i++) {
    x1++; y1++; x2--; y2--;
    if (x1 >= x2 || y1 >= y2) break;
    draw_rect_thickness1(x1, y1, x2, y2, value);
  }
  x1 = ul_x; y1 = ul_y; x2 = lr_x; y2 = lr_y;
  for (size_t i=1; i<=(thickness-1)/2; i++) {
    if (x1 == 0 || y1 == 0 || x2 == _width || y2 == _height) break;
    x1--; y1--; x2++; y2++;
    draw_rect_thickness1(x1, y1, x2, y2, value);
  }
}

