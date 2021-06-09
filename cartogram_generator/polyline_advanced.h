#ifndef POLYLINE_ADVANCED_H_
#define POLYLINE_ADVANCED_H_

#include <string>
#include <vector>

#include "cgal_typedef.h"

class Polyline_advanced {
  private:
    int pos_;
    Polyline pll_;
    Point v1_, v2_, vl_, v2l_;
    int gd_;
    int pgnwh_;
    bool is_hole_;

  public:
    Polyline_advanced(int pos,
        Polyline pll,
        int gd_v,
        int pgnwh_v,
        bool is_hole);

    int pos();
    Polyline pll();
    Point v1();
    Point v2();
    Point vl();
    Point v2l();
    int gd();
    int pgnwh();
    bool is_hole();
    void pop_front();
    void pop_back();
    void set_pll(Polyline pll);
    void set_is_hole(bool b);
};

#endif
