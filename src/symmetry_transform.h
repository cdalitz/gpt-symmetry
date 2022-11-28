// 
// symmetry_transform.cpp
//     Functions for computing the symmetry transform and
//     for discriminating symmetry types.
//
// Author:  Christoph Dalitz, Jens Wilberg
// Date:    2019-07-09
// License: see ../LICENSE
//

#ifndef _SYMMETRY_TRANSFORM_H
#define _SYMMETRY_TRANSFORM_H

#include "point.h"
#include "image.h"


//=======================================================================
// Functions for Algorithm 1 (GPT transform)
//=======================================================================

// Compute symmetry transform from a greyscale image.
// To each point, a symmetry score value and a radius (rx,ry) is assigned.
// Arguments:
//  grey       = greyscale image (input)
//  result     = symmetry scores as FloatImage (output)
//  result_rx  = x-component of symmetry radius as FloatImage (output)
//  result_ry  = y-component of symmetry radius as FloatImage (output)
//  grad_(x|y) = gradient image (input)
//  radius     = maximum radius examined (input)
//  norm_at_r  = symmetry weight (input)
//  which      = 1: dark objects, -1: light objects, 0: all objects (input)
void symmetry_transform(const GrayImage &grey, FloatImage* result, FloatImage* result_rx, FloatImage* result_ry, const FloatImage &grad_x, const FloatImage &grad_y, int radius, const std::vector<double> &norm_at_r, int which=0);


//=======================================================================
// Functions for Algorithm 2 (symmetric object detection)
//=======================================================================

// Find candidate symmetry centers in symmetry image as local
// maxima of the symmetry score. Returns the number of points found.
// Arguments:
//  symmetry        = symmetry transform image filename (input)
//  candidatepoints = symmetry candidate points (output)
//  windowsize      = window size for local maximum search (input)
size_t points_local_maxima(const FloatImage& symmetry, SymmetryPointVector *candidatepoints, int windowsize = 3);


// Compute the edge directedness as the frequency of the most frequent
// gradient direction.
// Arguments:
//  gradx    = x-component of gradient image (input)
//  grady    = y-component of gradient image (input)
//  point    = the point to be examined (input)
//  radius   = the radius of the window to be examined (input)
//  nbins    = number of bins in direction histogram (input)
//  cutoff   = ignore an image border of width cutoff (input)
double edge_directedness(const FloatImage &gradx, const FloatImage &grady, const Point& point, int radius = 3, int nbins = 16, int cutoff = 0);


// Compute the ratio of the the two eigenvalues of the covariance matrix
// of the window with the given radius around point.
// Arguments:
//  symmetry  = complete image (input)
//  point     = the point to be examined (input)
//  radius    = the radius of the window to be examined (input)
//  cutoff    = ignore an image border of width cutoff (input)
double coveigenratio(const FloatImage& symmetry, const Point& point, int radius = 3, int cutoff = 0);

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
double antiparallelity(const FloatImage& gradx, const FloatImage& grady, const Point& point, int radius, double &max_frequency, bool weighted = false, double count_threshold = 0.1, int nbins = 8);

// Returns the skeleton starting at x until the score value
// falls below the given percentage of S(x).
// Arguments:
//  S          = symmetry transform (input)
//  point      = starting point (input)
//  skeleton   = the points of the skeleton (output)
//  percentage = threshold at which the skeleton is stopped (input)
//  maxlength  = maximum length of the skeleton (input)
void get_skeleton(const FloatImage& symmetry, const Point& point, PointVector* skeleton, double percentage, int maxlength);

// Returns the ratio skeleton_length/symmetry_region_size
// Arguments:
//  skeleton = the points of the skeleton (input)
//  rx       = x-component of cricumradius
//  ry       = y-component of cricumradius
double skeleton_size(const PointVector &skeleton, double rx, double ry);

// Returns true, when QDA on all four features finds axial symmetry
// Input arguments:
//  edge_dir  = edge directedness
//  skel_size = skeleton size
//  anti_par  = antiparallelity
//  cov_ratio = covariance eigenratio
bool axial_qda(double edge_dir, double skel_size, double anti_par, double cov_ratio);

// maximum possible GPT value in region (x,y) +/- (rx,ry)
// Input arguments:
//  grad_abs  = image with gradient absolute values
//  x,y       = center of region
//  rx,ry     = radius of region
double s_max(const FloatImage& grad_abs, int x, int y, int rx, int ry);


//=======================================================================
// Functions for Algorithm 3 (medial/ridge axes extraction)
//=======================================================================

// ridge criterion by Staal et al. (2004) (see Eq. (11) in IPOL paper)
// Arguments:
//  x,y              = pixel position (must not be on border!) (input)
//  grad_(x|y)       = gradient (input)
//  hesse_(xx|xy|yy) = Hesse matrix (input)
bool staal_criterion(size_t x, size_t y, const FloatImage& grad_x, const FloatImage& grad_y, const FloatImage& hesse_xx, const FloatImage& hesse_xy, const FloatImage& hesse_yy);

// ridge criterion by Chang & Sinha (2007) (see Eq. (12) in IPOL paper)
// Arguments:
//  x,y       = pixel position (must not be on border!) (input)
//  img       = image (input)
bool ppa_criterion(size_t x, size_t y, const FloatImage& img);

#endif
