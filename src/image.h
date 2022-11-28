//
// image.h
//     Basic image utility functions.
//     File I/O is done with png++, but this is wrapped with own
//     image classes in order to minimize dependencies on 3rd party libs.
//
// Author:  Christoph Dalitz
// Date:    2019-07-09
// License: see ../LICENSE
//

#ifndef __IMAGE_H
#define __IMAGE_H

#include <stddef.h>


//---------------------------------------------------------------------
// Pixel types
//---------------------------------------------------------------------

typedef unsigned char GrayPixel;
typedef double FloatPixel;
struct RGBPixel {
  GrayPixel R; GrayPixel G; GrayPixel B;
  inline RGBPixel() {}
  inline RGBPixel(const RGBPixel& p) {R=p.R; G=p.G; B=p.B;}
  inline RGBPixel(GrayPixel r, GrayPixel g, GrayPixel b) {R=r; G=g; B=b;}
};


//---------------------------------------------------------------------
// Grayscale image
//---------------------------------------------------------------------

// forward declaration for conversion methods
class RGBImage;
class FloatImage;

class GrayImage {
private:
  size_t _width, _height;
  GrayPixel* _data;
public:
  GrayImage();
  ~GrayImage();
  void initialize(size_t width, size_t height);
  void initialize(size_t width, size_t height, GrayPixel value);
  int read_png(const char* filename);
  int write_png(const char* filename);
  void to_rgb(RGBImage* rgb);
  void to_float(FloatImage* rgb);
  inline size_t width() const {return _width;}
  inline size_t height() const {return _height;}
  inline GrayPixel get(size_t x, size_t y) const {return _data[y*_width+x];}
  inline void set(size_t x, size_t y, GrayPixel value) {_data[y*_width+x] = value;}
  GrayPixel interpolate(double x, double y) const;
};

//---------------------------------------------------------------------
// Float image
//---------------------------------------------------------------------

class FloatImage {
private:
  size_t _width, _height;
  FloatPixel* _data;
public:
  FloatImage();
  ~FloatImage();
  void initialize(size_t width, size_t height);
  void initialize(size_t width, size_t height, FloatPixel value);
  void to_gray(GrayImage* gray, bool onlypositive=false, bool enhance=false);
  FloatPixel mean_value() const;
  FloatPixel min_value() const;
  inline size_t width() const {return _width;}
  inline size_t height() const {return _height;}
  inline FloatPixel get(size_t x, size_t y) const {return _data[y*_width+x];}
  inline void set(size_t x, size_t y, FloatPixel value) {_data[y*_width+x] = value;}
  FloatPixel interpolate(double x, double y) const;
};

//---------------------------------------------------------------------
// RGB image
//---------------------------------------------------------------------

class RGBImage {
private:
  size_t _width, _height;
  RGBPixel* _data;
  void draw_rect_thickness1(size_t ul_x, size_t ul_y, size_t lr_x, size_t lr_y, const RGBPixel& value);
public:
  RGBImage();
  ~RGBImage();
  void initialize(size_t width, size_t height);
  void initialize(size_t width, size_t height, const RGBPixel& value);
  int read_png(const char* filename);
  int write_png(const char* filename);
  void to_gray(GrayImage* gray);
  void draw_rect(size_t ul_x, size_t ul_y, size_t lr_x, size_t lr_y, const RGBPixel& value, size_t thickness=1);
  inline size_t width() const {return _width;}
  inline size_t height() const {return _height;}
  inline RGBPixel get(size_t x, size_t y) const {return _data[y*_width+x];}
  inline void set(size_t x, size_t y, const RGBPixel& value) {_data[y*_width+x] = value;}
};


#endif
