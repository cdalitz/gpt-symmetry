// 
// point.h
//     Convenience data types and functions for processing points.
//
// Author:  Christoph Dalitz, Jens Wilberg
// Date:    2019-07-09
// License: see ../LICENSE
//

#ifndef __POINT_H
#define __POINT_H

#include <stddef.h>
#include <vector>

//-----------------------------------------------------------------
// data types for argument passing
//-----------------------------------------------------------------

// point in image raster
struct Point {
  size_t x; size_t y;
  inline Point() {x=0; y=0;}
  inline Point(size_t xx, size_t yy) {x=xx; y=yy;}
  inline Point(const Point& p) {x=p.x; y=p.y;}
};
typedef std::vector<Point> PointVector;

// point with symmetry score
struct SymmetryPoint {
  Point point;
  double value;
  inline SymmetryPoint(const SymmetryPoint& s) {
    point = s.point; value = s.value;
  }
  inline SymmetryPoint(const Point& p, double v) {
    point = p; value = v;
  }
  inline bool operator<(const SymmetryPoint& s) const {
    return (value < s.value);
  }
};
typedef std::vector<SymmetryPoint> SymmetryPointVector;

// helper function for descending sort
inline bool sortinverse(const SymmetryPoint& p1, const SymmetryPoint& p2)
{
  return (p1.value>p2.value);
}

// computes the neighbors of a given point
// warning: one pixel wide image border is ignored
// return value: number of neighbors
inline size_t get_neighbors(const Point &p, PointVector* neighbors, size_t maxx, size_t maxy)
{
  neighbors->clear();
  if (p.y>0 && p.y<maxy && p.x>0 && p.x<maxx) {
    neighbors->push_back(Point(p.x-1,p.y-1));
    neighbors->push_back(Point(p.x,p.y-1));
    neighbors->push_back(Point(p.x+1,p.y-1));
    neighbors->push_back(Point(p.x+1,p.y));
    neighbors->push_back(Point(p.x+1,p.y+1));
    neighbors->push_back(Point(p.x,p.y+1));
    neighbors->push_back(Point(p.x-1,p.y+1));
    neighbors->push_back(Point(p.x-1,p.y));
  }
  return neighbors->size();
}


#endif
