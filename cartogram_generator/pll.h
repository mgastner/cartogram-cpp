#ifndef PLL_H_
#define PLL_H_

#include <string>
#include <vector>

#include "cgal_typedef.h"

class PLL {
  private:
    int pos;
    Polyline pll;
    Point v1, v2, vl, v2l;
    std::vector<int> gd_v;
    std::vector<int> pgnwh_v;
    bool bool_hole;

  public:
    PLL(int pos_,
        Polyline pll_,
        std::vector<int>
        gd_v_,
        std::vector<int> pgnwh_v_,
        bool bool_hole_);

    int get_pos();
    Polyline get_pll();
    Point get_v1();
    Point get_v2();
    Point get_vl();
    Point get_v2l();
    std::vector<int> get_gd_v();
    std::vector<int> get_pgnwh_v();
    bool get_bool_hole();
    void set_pll(Polyline pll_);
    void set_bool_hole(bool b);
};

#endif
