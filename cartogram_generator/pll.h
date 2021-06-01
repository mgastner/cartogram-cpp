#ifndef PLL_H_
#define PLL_H_

#include <string>
#include <vector>

#include "cgal_typedef.h"

class PLL {
  private:
    int pos_;
    Polyline pll_;
    Point v1_, v2_, vl_, v2l_;
    int gd_;
    int pgnwh_;
    bool is_hole_;

  public:
    PLL(int pos,
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
    void set_pll(Polyline pll);
    void set_is_hole(bool b);
};

#endif
