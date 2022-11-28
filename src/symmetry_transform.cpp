// 
// symmetry_transform.cpp
//     Functions for computing the symmetry transform and
//     for discriminating symmetry types.
//
// Author:  Christoph Dalitz, Jens Wilberg
// Date:    2019-10-29
// License: see ../LICENSE
//

#define _USE_MATH_DEFINES

#include <algorithm>
#include <math.h>
#include <stdio.h>
#include "point.h"
#include "image.h"
#include "image_processing.h"


//=======================================================================
// Mathematical helper functions for internal use
//=======================================================================

// Heaviside's step function
inline int theta(double x) { return (x > 0 ? 1 : 0); }

// signum function
inline int sign(double x) { return ((x > 0) ? 1 : ((x < 0) ? -1 : 0)); }

// eigenvalues |lambda1| > |lambda2| of symmetric 2x2 matrix
void eigenvalues2x2(double a11, double a12, double a22, double* lambda1, double* lambda2=NULL) {
  double x = (a11 + a22) / 2.0;
  double y = sqrt((a11-a22)*(a11-a22)/4.0 + a12*a12);
  *lambda1 = (x>=0 ? x+y : x-y);
  if (lambda2 != NULL)
    *lambda2 = (x>=0 ? x-y : x+y);
}

// eigenvector v of symmetric 2x2 matrix corresponding to eigenvalue lambda
void eigenvector2x2(double a11, double a12, double a22, double lambda, double* v1, double* v2) {
  if (a12 != 0.0) {
    double length = sqrt((lambda - a22)*(lambda - a22) + a12*a12);
    *v1 = (lambda - a22) / length;
    *v2 = a12 / length;
  } else {
    if (fabs(a11-lambda) < 1E-8) {
      *v1 = 1.0; *v2 = 0.0;
    } else {
      *v1 = 0.0; *v2 = 1.0;
    }
  }
}


//=======================================================================
// Functions for Algorithm 1 (GPT transform)
//=======================================================================

// Compute symmetry score S and radius R at point (x,y) with rectangular
// regions. Arguments:
//  grad(x|y) = gradient image (input)
//  Rx        = x-component of best radius (output)
//  Ry        = y-component of best radius (output)
//  S         = score at best radius (output)
//  radius    = maximum radius examined (input)
//  norm_at_r = precalculated normalisation factors r^alpha (input)
//  x,y       = point position (input)
//  which     = 1: dark objects, -1: light objects, 0: all objects (input)
void s_and_r_at_point(const FloatImage &gradx, const FloatImage &grady, int *Rx, int *Ry, double *S, size_t radius, const std::vector<double> &norm_at_r, size_t x, size_t y, int which = 0)
{
  long rx, ry, max_radius;
  double symmetry;

  // adjust radius so that it does not extend beyond the image border
  max_radius = radius;
  if (y < max_radius) max_radius = y;
  if (y >= gradx.height() - max_radius) max_radius = gradx.height() - y - 1;
  if (x < max_radius) max_radius = x;
  if (x >= gradx.width() - max_radius) max_radius = gradx.width() - x - 1;
  FloatImage weight_at_r;
  weight_at_r.initialize(max_radius + 1, max_radius + 1, 0);

  // compute weight in window with radius max_radius
  if (max_radius < 1) {
    *Rx = 0;
    *Ry = 0;
    *S = 0.0;
  } else {
    if (which == 0) {
      // center row with ry=1
      ry = 1;
      symmetry = -gradx.get(x, y + ry) * gradx.get(x, y - ry) -
                 grady.get(x, y + ry) * grady.get(x, y - ry);
      for (rx = 1; rx <= max_radius; rx++) {
        symmetry += -gradx.get(x + rx, y + ry) * gradx.get(x - rx, y - ry)
                    -grady.get(x + rx, y + ry) * grady.get(x - rx, y - ry)
                    -gradx.get(x + rx, y) * gradx.get(x - rx, y)
                    -grady.get(x + rx, y) * grady.get(x - rx, y)
                    -gradx.get(x + rx, y - ry) * gradx.get(x - rx, y + ry)
                    -grady.get(x + rx, y - ry) * grady.get(x - rx, y + ry);
        weight_at_r.set(rx, ry, symmetry);
      }
      // use recursion formula for rest
      for (ry = 2; ry <= max_radius; ry++) {
        rx = 1;
        symmetry = weight_at_r.get(rx, ry - 1);
        symmetry += -gradx.get(x + rx, y + ry) * gradx.get(x - rx, y - ry)
                    -grady.get(x + rx, y + ry) * grady.get(x - rx, y - ry)
                    -gradx.get(x, y + ry) * gradx.get(x, y - ry)
                    -grady.get(x, y + ry) * grady.get(x, y - ry)
                    -gradx.get(x - rx, y + ry) * gradx.get(x + rx, y - ry)
                    -grady.get(x - rx, y + ry) * grady.get(x + rx, y - ry);
        weight_at_r.set(rx, ry, symmetry);
        for (rx = 2; rx <= max_radius; rx++) {
          symmetry += weight_at_r.get(rx, ry - 1) - weight_at_r.get(rx - 1, ry - 1)
                      -gradx.get(x + rx, y + ry) * gradx.get(x - rx, y - ry)
                      -grady.get(x + rx, y + ry) * grady.get(x - rx, y - ry)
                      -gradx.get(x - rx, y + ry) * gradx.get(x + rx, y - ry)
                      -grady.get(x - rx, y + ry) * grady.get(x + rx, y - ry);
          weight_at_r.set(rx, ry, symmetry);
        }
      }
    } else {  // which != 0
      // center row with ry=1
      ry = 1;
      symmetry = -theta(which * ry * grady.get(x, y + ry)) *
        theta(which * (-ry) * grady.get(x, y - ry)) *
        (gradx.get(x, y + ry) * gradx.get(x, y - ry) +
         grady.get(x, y + ry) * grady.get(x, y - ry));
      for (rx = 1; rx <= max_radius; rx++) {
        symmetry += - theta(which * (rx * gradx.get(x + rx, y + ry)
                     + ry * grady.get(x + rx, y + ry))) *
          theta(which * (-rx * gradx.get(x - rx, y - ry)
                         - ry * grady.get(x - rx, y - ry))) *
          (gradx.get(x + rx, y + ry) * gradx.get(x - rx, y - ry)
           + grady.get(x + rx, y + ry) * grady.get(x - rx, y - ry))
          
          - theta(which * rx * gradx.get(x + rx, y)) *
          theta(which * (-rx) * gradx.get(x - rx, y)) *
          (gradx.get(x + rx, y) * gradx.get(x - rx, y)
           + grady.get(x + rx, y) * grady.get(x - rx, y))
          
          - theta(which * (rx * gradx.get(x + rx, y - ry)
                          -ry * grady.get(x + rx, y - ry))) *
          theta(which * (-rx * gradx.get(x - rx, y + ry)
                          +ry * grady.get(x - rx, y + ry))) *
          (gradx.get(x + rx, y - ry) * gradx.get(x - rx, y + ry)
                     +grady.get(x + rx, y - ry) * grady.get(x - rx, y + ry));
        weight_at_r.set(rx, ry, symmetry);
      }
      // use recursion formula for rest
      for (ry = 2; ry <= max_radius; ry++) {
        rx = 1;
        symmetry = weight_at_r.get(rx, ry - 1);
        symmetry += - theta(which * (rx * gradx.get(x + rx, y + ry)
                                     +ry * grady.get(x + rx, y + ry))) *
          theta(which * (-rx * gradx.get(x - rx, y - ry)
                         -ry * grady.get(x - rx, y - ry))) *
          (gradx.get(x + rx, y + ry) * gradx.get(x - rx, y - ry)
           + grady.get(x + rx, y + ry) * grady.get(x - rx, y - ry))
          
          - theta(which * ry * grady.get(x, y + ry)) *
          theta(which * (-ry) * grady.get(x, y - ry)) *
          (gradx.get(x, y + ry) * gradx.get(x, y - ry)
           + grady.get(x, y + ry) * grady.get(x, y - ry))
          
          - theta(which * (-rx * gradx.get(x - rx, y + ry)
                           +ry * grady.get(x - rx, y + ry))) *
          theta(which * (rx * gradx.get(x + rx, y - ry)
                           -ry * grady.get(x + rx, y - ry))) *
          (gradx.get(x - rx, y + ry) * gradx.get(x + rx, y - ry)
           + grady.get(x - rx, y + ry) * grady.get(x + rx, y - ry));
        weight_at_r.set(rx, ry, symmetry);
        for (rx = 2; rx <= max_radius; rx++) {
          symmetry += weight_at_r.get(rx, ry - 1)
            - weight_at_r.get(rx - 1, ry - 1)
            - theta(which * (rx * gradx.get(x + rx, y + ry)
                            +ry * grady.get(x + rx, y + ry))) *
            theta(which * (-rx * gradx.get(x - rx, y - ry)
                            -ry * grady.get(x - rx, y - ry))) *
            (gradx.get(x + rx, y + ry) * gradx.get(x - rx, y - ry)
             + grady.get(x + rx, y + ry) * grady.get(x - rx, y - ry))
            
            - theta(which * (-rx * gradx.get(x - rx, y + ry)
                             +ry * grady.get(x - rx, y + ry))) *
            theta(which * (rx * gradx.get(x + rx, y - ry)
                           -ry * grady.get(x + rx, y - ry))) *
            (gradx.get(x - rx, y + ry) * gradx.get(x + rx, y - ry)
             + grady.get(x - rx, y + ry) * grady.get(x + rx, y - ry));
          weight_at_r.set(rx, ry, symmetry);
        }
      }
    }

    // find radius of maximum symmetry
    *Rx = 1;
    *Ry = 1;
    *S = weight_at_r.get(1, 1) / norm_at_r[2];
    for (ry = 1; ry <= max_radius; ry++) {
      for (rx = 1; rx <= max_radius; rx++) {
        if (weight_at_r.get(rx, ry) / norm_at_r[rx + ry] > *S) {
          *Rx = rx;
          *Ry = ry;
          *S = weight_at_r.get(rx, ry) / norm_at_r[rx + ry];
        }
      }
    }
  }
}

// Compute symmetry transform from a grayscale image.
// To each point, a symmetry score value and a radius (rx,ry) is assigned.
// Arguments:
//  gray       = grayscale image (input)
//  result     = symmetry scores as FloatImage (output)
//  result_rx  = x-component of symmetry radius as FloatImage (output)
//  result_ry  = y-component of symmetry radius as FloatImage (output)
//  radius     = maximum radius examined (input)
//  norm_at_r  = symmetry weight (input)
//  grad_(x|y) = gradient image (input)
//  which      = 1: dark objects, -1: light objects, 0: all objects (input)
void symmetry_transform(const GrayImage &gray, FloatImage* result, FloatImage* result_rx, FloatImage* result_ry, const FloatImage &grad_x, const FloatImage &grad_y, int radius, const std::vector<double> &norm_at_r, int which=0)
{
  size_t height = gray.height(), width = gray.width();
  // initialize result images
  result->initialize(width, height);
  result_rx->initialize(width, height);
  result_ry->initialize(width, height);

#if (NUM_THREADS > 1)
  #pragma omp parallel num_threads(NUM_THREADS) shared(result, result_rx, result_ry)
  #pragma omp for schedule(dynamic)
#endif
  for (size_t x = 0; x < width; ++x) {
    for (size_t y = 0; y < height; ++y) {
      int Rx, Ry;
      double S;
      s_and_r_at_point(grad_x, grad_y, &Rx, &Ry, &S, radius, norm_at_r, x, y, which);
      result->set(x, y, S);
      result_rx->set(x, y, (double)Rx);
      result_ry->set(x, y, (double)Ry);
    }
  }
}


//=======================================================================
// Functions for Algorithm 2 (symmetric object detection)
//=======================================================================

// Find candidate symmetry centers in symmetry image as local
// maxima of the symmetry score. Returns the number of points found.
// Arguments:
//  symmetry        = symmetry transform image filename (input)
//  candidatepoints = symmetry candidate points (output)
//  windowsize      = window size for local maximum search (input)
size_t points_local_maxima(const FloatImage& symmetry, SymmetryPointVector *candidatepoints, int windowsize /*= 3*/)
{
   if (windowsize % 2 == 0) {
    fprintf(stderr, "odd number for windowsize parameter expected\n");
    return -1;
  }
  candidatepoints->clear();
  double minv = symmetry.min_value();
  size_t s = windowsize / 2;
  int maxindex = -1;
  int maxi = s * windowsize + s;
  for (size_t y = s; y + s + 1 <= symmetry.height(); y++)
    for (size_t x = s; x + s + 1 <= symmetry.width(); x++) {
      double maxivalue = minv;
      for (size_t yy = 0; yy <= 2 * s; yy++)
        for (size_t xx = 0; xx <= 2 * s; xx++) {
          double d = symmetry.get(x + xx - s, y + yy - s);
          if (d > maxivalue) {
            maxindex = yy * windowsize + xx;
            maxivalue = d;
          }
        }
      if ((maxindex == maxi) && (maxivalue > 0.0)) {
        candidatepoints->push_back(
            SymmetryPoint(Point(x, y), symmetry.get(x, y)));
      }
    }

  return candidatepoints->size();
}


// Compute the edge directedness as the frequency of the most frequent
// gradient direction.
// Arguments:
//  gradx    = x-component of gradient image (input)
//  grady    = y-component of gradient image (input)
//  point    = the point to be examined (input)
//  radius   = the radius of the window to be examined (input)
//  nbins    = number of bins in direction histogram (input)
//  cutoff   = ignore an image border of width cutoff (input)
double edge_directedness(const FloatImage &gradx, const FloatImage &grady, const Point& point, int radius = 3, int nbins = 16, int cutoff = 0)
{
  int i, x, y, left, right, top, bot;
  double G, sumG;  // absolute values of gradient
  double gx, gy;   // gradient values
  double angle, maxhist, pi_n;
  double* hist = new double[nbins];
  
  pi_n = M_PI / nbins;
  left  = std::max(cutoff, (int)point.x - radius);
  right = std::min((int)gradx.width() - cutoff - 1, (int)point.x + radius);
  top   = std::max(cutoff, (int)(point.y - radius));
  bot   = std::min((int)gradx.height() - cutoff - 1, (int)point.y + radius);

  sumG = 0.0;
  for (i=0; i<nbins; i++) hist[i] = 0.0;
  for (y = top; y <= bot; y++) {
    for (x = left; x <= right; x++) {
      gx = gradx.get(x, y); gy = grady.get(x, y);
      G = sqrt(gx*gx + gy*gy);
      sumG += G;
      angle = atan2(gy,gx);
      // shift angle bins by pi/n so that zero etc. lie in middle of bin
      if (angle < -pi_n) angle += 2*M_PI;
      i = (int)((angle+pi_n)*nbins/(2*M_PI));
      hist[i] += G;
    }
  }
  maxhist = 0.0;
  for (i=0; i<nbins; i++) 
    if (hist[i] > maxhist) maxhist = hist[i];

  delete[] hist;

  if (sumG == 0.0)
    return 0.0;
  else
    return (maxhist/sumG);
}


// Compute the ratio of the the two eigenvalues of the covariance matrix
// of the window with the given radius around point.
// Arguments:
//  symmetry  = complete image (input)
//  point     = the point to be examined (input)
//  radius    = the radius of the window to be examined (input)
//  cutoff    = ignore an image border of width cutoff (input)
double coveigenratio(const FloatImage& symmetry, const Point& point, int radius = 3, int cutoff = 0)
{
  int left,right,top,bot;
  double cov[2][2];
  double lambda1, lambda2;
  double sumf;  // normalization factor
  double dy,dx;
  double value = 0.0;

  left  = std::max(cutoff, (int)point.x - radius);
  right = std::min((int)symmetry.width() - cutoff - 1, (int)point.x + radius);
  top   = std::max(cutoff, (int)point.y - radius);
  bot   = std::min((int)symmetry.height() - cutoff - 1, (int)point.y + radius);

  // compute covariance matrix
  for (size_t i=0; i<2; i++){
    for (size_t j=0; j<2; j++){
      cov[i][j] = 0.0;
    }
  }
  sumf = 0.0;
  for (int y = top; y <= bot; y++) {
    for (int x = left; x <= right; x++) {
      value = symmetry.get(x, y);
      sumf += value;
      dx = x - (int)point.x; dy = y - (int)point.y;
      cov[0][0] += value*dx*dx ;
      cov[0][1] += value*dx*dy;
      cov[1][0] += value*dy*dx;
      cov[1][1] += value*dy*dy;
    }
  }
  for (size_t i=0; i<2; i++){
    for (size_t j=0; j<2; j++){
      cov[i][j] /= sumf;
    }
  }

  // determine eigenvalue ratio
  eigenvalues2x2(cov[0][0], cov[0][1], cov[1][1], &lambda1, &lambda2);

  // printf("Eigenvalues = (%f,%f)\n", lambda1, lambda2);

  if (lambda1 == 0.0 && lambda2 == 0.0)
    return 0.0;
  else
    return (fabs(lambda2) / fabs(lambda1));
}

// Returns the number of antiparallel gradient bins that have a
// frequency count larger than count_threshold
// inside the window with the given radius around point.
// Arguments:
//  gradx    = x-component of gradient image (input)
//  grady    = y-component of gradient image (input)
//  x,y      = the point to be examined (input)
//  radius   = the radius of the window to be examined (input)
//  max_frequency = frequency count of largest bin (output)
//  weighted = whether edges shall be weighted by gradient strength (input)
//  count_threshold = threshold when a bin is considered "occupied" (input)
//  nbins    = number of bins in direction histogram (input)
double antiparallelity(const FloatImage& gradx, const FloatImage& grady, const Point& point, int radius, double &max_frequency, bool weighted = false, double count_threshold = 0.1, int nbins = 8)
{
  double scalarprod, result;
  int x, y, dx, dy, i, count;
  double v1,v2;  // norms of the opposite gradients
  double angle;  // angle between gradients
  double sumG;   // normalization factor (sum of all gradient norms)
  // threshold when angle between two vectors can be considered antiparallel
  double parallelity_threshold = -0.975;
  // pi divided by number of bins (shift for histogram bins)
  double pi_n = M_PI / nbins;

  // compute antiparallelity
  std::vector<double> bin(nbins);
  for (i=0; i<nbins; i++) bin[i] = 0.0;
  sumG = 0.0;
  x = point.x;
  y = point.y;
  for (int r=1; r<=radius; r++) {
    dy = r;
    for (dx=-r; dx<=r; dx++) {
      scalarprod = gradx.get(x+dx, y+dy)*gradx.get(x-dx, y-dy)
        + grady.get(x+dx, y+dy)*grady.get(x-dx, y-dy);
      if (scalarprod < 0.0) { // only when in opposite directions
        v1 = sqrt(gradx.get(x+dx, y+dy)*gradx.get(x+dx, y+dy) +
                  grady.get(x+dx, y+dy)*grady.get(x+dx, y+dy));
        v2 = sqrt(gradx.get(x-dx, y-dy)*gradx.get(x-dx, y-dy) +
                  grady.get(x-dx, y-dy)*grady.get(x-dx, y-dy));
        if (scalarprod/(v1*v2) < parallelity_threshold) {
          angle = atan2(grady.get(x+dx, y+dy),gradx.get(x+dx, y+dy));
          // shift angle bins by pi/n so that zero etc. lie in middle of bin
          if (angle < -pi_n) angle += 2*M_PI;
          i = (int)((angle+pi_n)*nbins/(2*M_PI));
          if (weighted) {
            bin[i] += v1;
            bin[(i+nbins/2)%nbins] += v2;
            sumG += v1+v2;
          } else {
            bin[i]++;
            bin[(i+nbins/2)%nbins]++;
            sumG += 2.0;
          }
        }
      }
    }
    dx = r;
    for (dy=-r+1; dy<r; dy++) {
      scalarprod = gradx.get(x+dx, y+dy)*gradx.get(x-dx, y-dy)
        + grady.get(x+dx, y+dy)*grady.get(x-dx, y-dy);
      if (scalarprod < 0.0) { // only when in opposite directions
        v1 = sqrt(gradx.get(x+dx, y+dy)*gradx.get(x+dx, y+dy) +
                  grady.get(x+dx, y+dy)*grady.get(x+dx, y+dy));
        v2 = sqrt(gradx.get(x-dx, y-dy)*gradx.get(x-dx, y-dy) +
                  grady.get(x-dx, y-dy)*grady.get(x-dx, y-dy));
        if (scalarprod/(v1*v2) < parallelity_threshold) {
          angle = atan2(grady.get(x+dx, y+dy),gradx.get(x+dx, y+dy));
          // shift angle bins by pi/n so that zero etc. lie in middle of bin
          if (angle < -pi_n) angle += 2*M_PI;
          i = (int)((angle+pi_n)*nbins/(2*M_PI));
          if (weighted) {
            bin[i] += v1;
            bin[(i+nbins/2)%nbins] += v2;
            sumG += v1+v2;
          } else {
            bin[i]++;
            bin[(i+nbins/2)%nbins]++;
            sumG += 2.0;
          }
        }
      }
    }
  }

  for (i=0; i<nbins; i++) bin[i] = bin[i]/sumG;
  count = 0;
  for (i=0; i<nbins; i++)
    if (bin[i] > count_threshold) count++;

  max_frequency = *std::max_element(bin.begin(), bin.end());
  result = count/((double)nbins);
  return result;
}


// Returns the skeleton starting at x until the score value
// falls below the given percentage of S(x).
// Arguments:
//  S          = symmetry transform (input)
//  point      = starting point (input)
//  skeleton   = the points of the skeleton (output)
//  percentage = threshold at which the skeleton is stopped (input)
//  maxlength  = maximum length of the skeleton (input)
void get_skeleton(const FloatImage& symmetry, const Point& point, PointVector* skeleton, double percentage, int maxlength)
{
  int i,nn,r;
  double maxval;
  Point maxpoint, startpoint;
  PointVector neighbors;
  size_t maxx = symmetry.width() - 1;
  size_t maxy = symmetry.height() - 1;
  double threshold = percentage*symmetry.get(point.x, point.y);
  skeleton->clear();
  skeleton->push_back(point);

  // find startpoint in first direction
  nn = get_neighbors(point, &neighbors, maxx, maxy);
  if (!nn) return;
  maxpoint = neighbors[0];
  maxval = symmetry.get(maxpoint.x, maxpoint.y);
  for (i=1; i<nn; i++) {
    if (symmetry.get(neighbors[i].x, neighbors[i].y) > maxval) {
      maxpoint = neighbors[i]; maxval = symmetry.get(maxpoint.x, maxpoint.y);
    }
  }
  startpoint = maxpoint; // remember for opposite direction later
  if (maxval > threshold)
    skeleton->push_back(maxpoint);

  // follow skeleton in first direction
  for (r=2; (r<maxlength) && (maxval>threshold); r++) {
    nn = get_neighbors(maxpoint, &neighbors, maxx, maxy);
    if (!nn) break;
    //maxval = std::numeric_limits<double>::lowest();
    maxval = -1000;
    for (i=0; i<nn; i++) {
      if ((abs((int)neighbors[i].x-(int)point.x)==r || abs((int)neighbors[i].y-(int)point.y)==r)
          && (symmetry.get(neighbors[i].x, neighbors[i].y) > maxval)) {
        maxpoint = neighbors[i]; maxval = symmetry.get(maxpoint.x, maxpoint.y);
      }
    }
    if (maxval > threshold) {
      skeleton->push_back(maxpoint);
    }
  }

  // look for startpoint near mirrored position
  maxpoint = Point(2*point.x-startpoint.x, 2*point.y-startpoint.y);
  maxval = symmetry.get(maxpoint.x, maxpoint.y);
  nn = get_neighbors(maxpoint, &neighbors, maxx, maxy);
  if (!nn) return;
  for (i=0; i<nn; i++) {
    if ((abs((int)neighbors[i].x-(int)point.x)==1 || abs((int)neighbors[i].y-(int)point.y)==1)
        && (symmetry.get(neighbors[i].x, neighbors[i].y) > maxval)) {
      maxpoint = neighbors[i]; maxval = symmetry.get(maxpoint.x, maxpoint.y);
    }
  }
  if (maxval > threshold)
    skeleton->push_back(maxpoint);

  // follow skeleton in opposite direction
  for (r=2; (r<maxlength) && (maxval>threshold); r++) {
    nn = get_neighbors(maxpoint, &neighbors, maxx, maxy);
    if (!nn) break;
    //maxval = std::numeric_limits<double>::lowest();
    maxval = -1000;
    for (i=0; i<nn; i++) {
      if ((abs((int)neighbors[i].x-(int)point.x)==r || abs((int)neighbors[i].y-(int)point.y)==r)
          && (symmetry.get(neighbors[i].x, neighbors[i].y) > maxval)) {
        maxpoint = neighbors[i]; maxval = symmetry.get(maxpoint.x, maxpoint.y);
      }
    }
    if (maxval > threshold) {
      skeleton->push_back(maxpoint);
    }
  }
}

// Returns the ratio skeleton_length/symmetry_region_size
// Arguments:
//  skeleton = the points of the skeleton (input)
//  rx       = x-component of cricumradius
//  ry       = y-component of cricumradius
double skeleton_size(const PointVector &skeleton, double rx, double ry)
{
  double result = sqrt(rx*rx+ry*ry);
  result = double(skeleton.size()) / result;
  if (result > 1.0)
    return 1.0;
  else
    return result;
}

//-----------------------------------------------------------------
// Returns true, when QDA on all four features finds axial symmetry
// Input arguments:
//  edge_dir  = edge directedness
//  skel_size = skeleton size
//  anti_par  = antiparallelity
//  cov_ratio = covariance eigenratio
//-----------------------------------------------------------------
bool axial_qda(double edge_dir, double skel_size, double anti_par, double cov_ratio)
{
  // values measured from test data set
  const double ma[4] = {0.4334325, 0.4974716, 0.4417170, 0.4352612};
  const double mr[4] = {0.18586523, 0.06537124, 0.23017472, 0.71357609};
  const double Da = -15.03338;
  const double Dr = -21.49041;
  const double Sa[4][4] = {{-9.392123, 3.145148, 7.8727492, 0.8677436},
                           {0.000000, -2.436357, 0.3521709, 0.1780096},
                           {0.000000,  0.0000, -18.5275062, 3.1907446},
                           {0.000000,  0.000000, 0.0000000, 4.3364655}};
  const double Sr[4][4] = {{14.31157, -6.191231, 15.545277, -6.7814399},
                           {0.00000,  20.318334,  0.375387, -0.9648399},
                           {0.00000,   0.00000, -20.595436, -3.8729805},
                           {0.00000,   0.00000,   0.000000, -7.7488568}};
  double ga, gr, tmp;
  double x[4];
  size_t j,k;
  // discriminant function for axial
  x[0] = edge_dir  - ma[0];
  x[1] = skel_size - ma[1];
  x[2] = anti_par  - ma[2];
  x[3] = cov_ratio - ma[3];
  ga = -Da;
  for (j=0; j<4; j++) {
    tmp = 0.0;
    for (k=0; k<4; k++) {
      tmp += x[k]*Sa[k][j];
    }
    ga = ga - tmp*tmp;
  }
  // discriminant function for rotational
  x[0] = edge_dir  - mr[0];
  x[1] = skel_size - mr[1];
  x[2] = anti_par  - mr[2];
  x[3] = cov_ratio - mr[3];
  gr = -Dr;
  for (j=0; j<4; j++) {
    tmp = 0.0;
    for (k=0; k<4; k++) {
      tmp += x[k]*Sr[k][j];
    }
    gr = gr - tmp*tmp;
  }
  // decision
  return (ga > gr);
}

// maximum possible GPT value in region (x,y) +/- (rx,ry)
// see Eq. (8) in IPOL paper
// Input arguments:
//  grad_abs  = image with gradient absolute values
//  x,y       = center of region
//  rx,ry     = radius of region
double s_max(const FloatImage& grad_abs, int x, int y, int rx, int ry)
{
  double s = 0.0;
  int dx, dy;
  for (dx=1; dx<=rx; dx++) {
    s += grad_abs.get(x+dx, y)*grad_abs.get(x-dx, y);
  }
  for (dy=1; dy<=ry; dy++) {
    for (dx=-rx; dx<=rx; dx++) {
      s += grad_abs.get(x+dx, y+dy)*grad_abs.get(x-dx, y-dy);
    }
  }
  return s;
}

//=======================================================================
// Functions for Algorithm 3 (medial/ridge axes extraction)
//=======================================================================

// ridge criterion by Staal et al. (2004) (see Eq. (11) in IPOL paper)
// Arguments:
//  x,y              = pixel position (must not be on border!) (input)
//  grad_(x|y)       = gradient (input)
//  hesse_(xx|xy|yy) = Hesse matrix (input)
bool staal_criterion(size_t x, size_t y, const FloatImage& grad_x, const FloatImage& grad_y, const FloatImage& hesse_xx, const FloatImage& hesse_xy, const FloatImage& hesse_yy)
{
  double lambda, ex, ey, staal_crit;
  // compute eigenvector of largest eigenvalue of Hessian
  eigenvalues2x2(hesse_xx.get(x,y), hesse_xy.get(x,y),
                 hesse_yy.get(x,y), &lambda);
  eigenvector2x2(hesse_xx.get(x,y), hesse_xy.get(x,y),
                 hesse_yy.get(x,y), lambda, &ex, &ey);
  
  // Eq. (11)
  staal_crit = -0.5 * sign(lambda) *
    fabs(sign(grad_x.interpolate(double(x)+ex, double(y)+ey)*ex +
              grad_y.interpolate(double(x)+ex, double(y)+ey)*ey)
         - sign(grad_x.interpolate(double(x)-ex, double(y)-ey)*ex +
                grad_y.interpolate(double(x)-ex, double(y)-ey)*ey));

  return (staal_crit > 0.5);
}


// ridge criterion by Chang & Sinha (2007) (see Eq. (12) in IPOL paper)
// Arguments:
//  x,y       = pixel position (must not be on border!) (input)
//  img       = image (input)
bool ppa_criterion(size_t x, size_t y, const FloatImage& img)
{
  double refvalue = img.get(x, y);
  if ((refvalue > img.get(x-1, y-1) &&
       refvalue > img.get(x+1, y+1)) ||
      (refvalue > img.get(x, y-1) &&
       refvalue > img.get(x, y+1)) ||
      (refvalue > img.get(x+1, y-1) &&
       refvalue > img.get(x-1, y+1)) ||
      (refvalue > img.get(x-1, y) &&
       refvalue > img.get(x+1, y)))
    return true;
  else
    return false;
}


// Find ridge points in symmetry transform image with combined methods
// of Staal et al. (2004) and Chang & Sinha (2007).
// Returns the number of points found.
// Arguments:
//  gpt_(s|rx|ry)    = GPT transform images (input)
//  alpha            = power used in symmetry score image (input)
//  grad_abs         = gradient absolute value of original image (input)
//  gpt_grad(x|y)    = gradient of symmetry transform (input)
//  result           = detected ridge points (output)
//  threshold        = minimum value of s_norm required for ridges
size_t ridge_points(const FloatImage& gpt_s, const FloatImage& gpt_rx, const FloatImage& gpt_ry, double alpha, const FloatImage& grad_abs, const FloatImage& gpt_gradx, const FloatImage& gpt_grady, PointVector *result, double threshold /*=0.4*/)
{
  double m, rx, ry, refvalue;

  result->clear();

  // mean gradient value is reference value for threshold Eq. (14)
  m = grad_abs.mean_value();

  // Hesse matrix of GPT score image
  FloatImage gpt_hesse_xx, gpt_hesse_xy, gpt_hesse_yy;
  sobel_x(gpt_gradx, &gpt_hesse_xx);
  sobel_x(gpt_grady, &gpt_hesse_xy);
  sobel_y(gpt_grady, &gpt_hesse_yy);
  
  for (size_t y=1; y+1<gpt_s.height(); y++) {
    for (size_t x=1; x+1<gpt_s.width(); x++) {
 
      refvalue = gpt_s.get(x,y);
      if (refvalue < 0.0)
        continue;
      rx = gpt_rx.get(x,y);
      ry = gpt_ry.get(x,y);
      if (refvalue*pow(rx+ry, alpha) < threshold *
          s_max(grad_abs, x, y, rx, ry))
        continue;

      // test for Staal's criterion
      if (!staal_criterion(x, y, gpt_gradx, gpt_grady,
                           gpt_hesse_xx, gpt_hesse_xy, gpt_hesse_yy))
        continue;

      // test for Chang & Sinha's criterion
      if (!ppa_criterion(x, y, gpt_s))
        continue;
      
      if (refvalue*pow(rx+ry, alpha) < 2*std::max(rx,ry)*m*m) {
        //bgr.at<cv::Vec3b>(y,x) = cv::Vec3b(200,200,0);
      } else {
        result->push_back(Point(x,y));
        //bgr.at<cv::Vec3b>(y,x) = cv::Vec3b(0,255,255);
      }
    }
  }
  
  return result->size();
}

