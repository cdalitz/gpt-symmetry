//
// Computes the symmetry transformation after
// IPOL 2019 paper by Dalitz and Wilberg
// and return the points and regions with score value greter than threshold
//
// Author:  Christoph Dalitz, Jens Wilberg
// Date:    2019-10-29
// License: see ../LICENSE
//

#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "image.h"
#include "image_processing.h"
#include "point.h"
#include "symmetry_transform.h"

#define TRACEPREFIX "trace-"
#define VERSION "1.1"

// command line arguments
char* opt_infile = NULL;
size_t opt_r = 0;
double opt_alpha = 0.5;
double opt_t = -1.0; // impossible value signals: not set by user
size_t opt_k = 5;
int opt_which = 0; // 0 = all, 1 = dark, -1 = light
size_t opt_borderwidth = 1;
bool opt_skipinside = false;
bool opt_trace = true;
bool opt_scale = true;
char* opt_outfile = NULL;
const char* usage = "Usage:\n\tgptsymmetry -r <r> [options] <infile>\nOptions:\n"
  "\t-r <r> maximum radius of symmetry transform\n"
  "\t-alpha <alpha>\n"
  "\t       exponent for normalisation of symmetry score with r^alpha [0.5]\n"
  "\t-t <threshold>\n"
  "\t       threshold t for symmetry score [0.6 if which=all, 0.4 else]\n"
  "\t       when t >= 1, only the first found rotational symmetry is reported\n"
  "\t-k <k>\n"
  "\t       the radius of the window for local maxima search [5]\n"
  "\t-which (all|dark|light)\n"
  "\t       look for symmetries in all/dark/light objects [all]\n"
  "\t-skipinside\n"
  "\t       skip symmetries within box of a previously found symmetry\n"
  "\t-o <outpng>\n"
  "\t       draw detected symmetry as rectangle in file <outpng>\n"
  "\t-bw <borderwidth>\n"
  "\t       width of the borders of symmetry regions drawn in result image [1]\n"
  "\t-notrace\n"
  "\t       do not write images of interim steps with prefix '" TRACEPREFIX "'\n"
  "\t-noscale\n"
  "\t       do not scale large images with large radius down\n"
  "\t-v     print version and exit\n";


int main(int argc, char** argv)
{
  //
  // parse command line
  //
  for (int i=1; i<argc; i++) {
    if (0==strcmp("-r", argv[i])) {
      i++; if (i<argc) opt_r = atoi(argv[i]);
    }
    else if (0==strcmp("-alpha", argv[i])) {
      i++; if (i<argc) opt_alpha = atof(argv[i]);
    }
    else if (0 == strcmp("-t", argv[i])) {
      i++;
      if (i < argc) opt_t = atof(argv[i]);
    }
    else if (0 == strcmp("-k", argv[i])) {
      i++;
      if (i < argc) opt_k = atoi(argv[i]);
    }
    else if (0==strcmp("-which", argv[i])) {
      i++; if (i==argc) {
        fputs(usage, stderr);
        return 1;
      }
      if (0==strcmp("all", argv[i])) {
        opt_which = 0;
      }
      else if (0==strcmp("dark", argv[i])) {
        opt_which = 1;
      }
      else if (0==strcmp("light", argv[i])) {
        opt_which = -1;
      }
      else {
        fputs(usage, stderr);
        return 1;
      }
    }
    else if (0 == strcmp("-skipinside", argv[i])) {
      opt_skipinside = true;
    }
    else if (0==strcmp("-o", argv[i])) {
      i++; if (i<argc) opt_outfile = argv[i];
    }
    else if (0 == strcmp("-bw", argv[i])) {
      i++;
      if (i < argc) opt_borderwidth = atoi(argv[i]);
    }
    else if (0==strcmp("-notrace", argv[i])) {
      opt_trace = false;
    }
    else if (0==strcmp("-noscale", argv[i])) {
      opt_scale = false;
    }
    else if (0==strcmp("-v", argv[i])) {
      printf("%s, compiled with NUM_THREADS=%i\n", VERSION, NUM_THREADS);
      return 0;
    }
    else if (argv[i][0] == '-') {
      fputs(usage, stderr);
      return 1;
    }
    else {
      opt_infile = (char*)argv[i];
    } 
  }
  if (!opt_r || !opt_infile) {
    fputs(usage, stderr);
    return 1;
  }
  // set reasonable default threshold if not set by user
  if (opt_t < 0.0) {
    if (opt_which == 0)
      opt_t = 0.6;
    else
      opt_t = 0.4;
  }

  //
  // load image and convert to grayscale
  //
  GrayImage img;
  if (img.read_png(opt_infile)){
    return 1;
  }

  //
  // if combination of opt_r and image size too large: scale image down
  //
  size_t max_r = std::min(img.width(), img.height()) / 2;
  // check if opt_r is plausible
  if (opt_r > max_r) opt_r = max_r;
  GrayImage gray;
  double scale = 1.0;
  double tmp = sqrt(img.width() * img.height()) * opt_r;
  // check if image needs to be scaled
  if(tmp > 100000){
    // compute the scaling value
    scale = sqrt(100000.0/tmp);
    opt_r *= scale;
  }
  // scale image. If *scale* is >= 1 the original image is copied
  down_scale_image(img, &gray, scale);

  //
  // compute GPT score S and radius (Rx,Ry) (Algorithm 1)
  // ----------------------------------------------------
  //
  FloatImage result, result_rx, result_ry, gradx, grady, tmp_f;

  // compute x,y gradient
  gray.to_float(&tmp_f);
  sobel_x(tmp_f, &gradx);
  sobel_y(tmp_f, &grady);

  // precompute symmetry weights (for speedup)
  std::vector<double> norm_at_r(2*opt_r+1);
  norm_at_r[0] = 1.0; // not needed
  for (size_t r=1; r<=2*opt_r; r++) norm_at_r[r] = pow(r,opt_alpha);

  // call symmetry transform GPT
  symmetry_transform(gray, &result, &result_rx, &result_ry, gradx, grady, opt_r, norm_at_r, opt_which);

  // write result as grayscale image
  GrayImage gray_gpt;
  if (opt_trace) {
    //result.to_gray(&gray_gpt, true, true);
    //gray_gpt.write_png(TRACEPREFIX "symmetry-equalized.png");
    result.to_gray(&gray_gpt, true);
    gray_gpt.write_png(TRACEPREFIX "symmetry.png");
  }

  // precompute some gradients needed by both Algorithm 2 and 3
  FloatImage gpt_gradx, gpt_grady;
  FloatImage grad_abs;
  l2_norm(gradx, grady, &grad_abs);
  sobel_x(result, &gpt_gradx);
  sobel_y(result, &gpt_grady);


  //
  // extract strong symmetry axes with ridge detection (Algorithm 3)
  // ---------------------------------------------------------------
  //

  // we don't need it when no trace image is asked for
  if (opt_trace) {
    // Hesse matrix of GPT score image
    FloatImage gpt_hesse_xx, gpt_hesse_xy, gpt_hesse_yy;
    sobel_x(gpt_gradx, &gpt_hesse_xx);
    sobel_x(gpt_grady, &gpt_hesse_xy);
    sobel_y(gpt_grady, &gpt_hesse_yy);

    // mean gradient value is reference value for threshold Eq. (14)
    double mean_grad2 = grad_abs.mean_value();
    mean_grad2 = mean_grad2*mean_grad2;

    //
    // test for each point whether it is a ridge point
    //
    PointVector ridgepoints;
    double refvalue, rx, ry;
    for (size_t y=1; y+1<result.height(); y++) {
      for (size_t x=1; x+1<result.width(); x++) {
 
        refvalue = result.get(x,y);
        if (refvalue < 0.0) continue;  // shortcut for speedup
        rx = result_rx.get(x,y);
        ry = result_ry.get(x,y);
        // criterion Eq. (14)
        if (refvalue*pow(rx+ry, opt_alpha) < 2*std::max(rx,ry)*mean_grad2)
          continue;
        // criterion Eq. (11)
        if (!staal_criterion(x, y, gpt_gradx, gpt_grady,
                             gpt_hesse_xx, gpt_hesse_xy, gpt_hesse_yy))
          continue;
        // criterion Eq. (12)
        if (!ppa_criterion(x, y, result))
          continue;
        // criterion Eq. (13)
        if (refvalue*pow(rx+ry, opt_alpha) < 0.4 * s_max(grad_abs, x, y, rx, ry))
          continue;
        // all criteria match:
        ridgepoints.push_back(Point(x,y));
      }
    }

    // write ridge image
    RGBImage skel_img;
    gray.to_rgb(&skel_img);
    for (size_t i=0; i<ridgepoints.size(); i++)
      skel_img.set(ridgepoints[i].x, ridgepoints[i].y, RGBPixel(0,220,220));
    skel_img.write_png(TRACEPREFIX "ridges.png");
  }
  
  //
  // detect rotationally symmetric objects (Algorithm 2)
  // ---------------------------------------------------
  //

  // find all local maxima and sort them by symmetry score
  SymmetryPointVector sp;
  if (0 >= points_local_maxima(result, &sp, 2 * opt_k +1)) {
    fputs("no local maxima found\n", stderr);
    return 3;
  }
  std::sort(sp.begin(), sp.end());
  std::reverse(sp.begin(), sp.end());

  // initialize trace images
  RGBImage rgb, rgb_gpt;
  img.to_rgb(&rgb);
  if (opt_trace)
    gray_gpt.to_rgb(&rgb_gpt);

  // print header for output
  printf("rot;x;y;rx;ry;S;s_norm\n");

  // test for each pixel whether it matches the criteria
  // from lines 6-8 of Algorithm 2
  double edge_dir, anti_par, skel_size, cov_ratio, snorm;
  PointVector rot_sym, skeleton;
  bool isaxial;
  for (size_t i=0; i<sp.size(); i++) {
    Point p = sp[i].point;
    int rx = (int)result_rx.get(p.x, p.y);
    int ry = (int)result_ry.get(p.x, p.y);

    // optionally skip symmetries inside of a previously found region
    if (opt_skipinside) {
      int ry2, rx2;
      bool contained = false;
      Point p2;
      for (size_t j = 0; j < rot_sym.size(); j++) {
        p2 = rot_sym[j];
        rx2 = (int)result_rx.get(p2.x, p2.y);
        ry2 = (int)result_ry.get(p2.x, p2.y);
        if (p2.x - rx2 <= p.x && p2.x + rx2 >= p.x &&
            p2.y - rx2 <= p.y && p2.y + ry2 >= p.y) {
          contained = true;
          break;
        }
      }
      if (contained) continue;
    }

    // compute normalized score s_norm
    snorm = s_max(grad_abs, p.x, p.y, rx, ry);
    if (snorm > 1) {
      snorm = sp[i].value * norm_at_r[rx + ry] / snorm;
    } else {
      snorm = 0.0;
    }

    // classification axial versus rotational
    edge_dir = edge_directedness(gpt_gradx, gpt_grady, p, 3, 16);
    get_skeleton(result, p, &skeleton, 0.5, opt_r);
    skel_size = skeleton_size(skeleton, rx, ry);
    antiparallelity(gradx, grady, p, std::min(rx,ry), anti_par);
    cov_ratio = coveigenratio(result, p, 3);
    isaxial = axial_qda(edge_dir, skel_size, anti_par, cov_ratio);

    // axial symmetry: report and move on to next point
    if (isaxial) {
      p.x /= scale;
      p.y /= scale;
      rx /= scale;
      ry /= scale;
      printf("0;%lu;%lu;%i;%i;%f;%f\n", p.x, p.y, rx, ry, sp[i].value, snorm);
      if (opt_trace) {
        rgb_gpt.draw_rect(p.x - rx, p.y - ry, p.x + rx, p.y + ry,
                          RGBPixel(255, 0, 0), opt_borderwidth);
      }
      continue;
    }
 
    // stop if score of point < threshold
    if (snorm < opt_t && opt_t < 1.0) {
      // if (opt_trace) { // might be useful to visualize stopping region
      //   rgb_gpt.draw_rect(p.x - rx, p.y - ry, p.x + rx, p.y + ry,
      //                     RGBPixel(255, 255, 0), opt_borderwidth);
      // }
      break;
    }

    // report rotational symmetry
    rot_sym.push_back(Point(p.x,p.y));
    p.x /= scale;
    p.y /= scale;
    rx /= scale;
    ry /= scale;
    printf("1;%lu;%lu;%i;%i;%f;%f\n", p.x , p.y, rx, ry, sp[i].value, snorm);
    rgb.draw_rect(p.x - rx, p.y - ry, p.x + rx, p.y + ry, 
                  RGBPixel(0, 200, 200), opt_borderwidth);
    if (opt_trace) {
      rgb_gpt.draw_rect(p.x - rx, p.y - ry, p.x + rx, p.y + ry,
                        RGBPixel(0, 200, 200), opt_borderwidth);
    }

    // stop if only asked for first rotational symmetry
    if (opt_t >= 1.0){
      break;
    }
  }

  //
  // write result images
  //
  if (opt_trace) {
    rgb_gpt.write_png(TRACEPREFIX "localmaxpoints.png");
  }
  if (opt_outfile) {
    rgb.write_png(opt_outfile);
  }

  return 0;
}
