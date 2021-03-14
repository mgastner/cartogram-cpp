#include <string>
#include <vector>

#include "pll.h"

PLL::PLL(int pos_,
    Polyline pll_,
    int gd_,
    int pgnwh_,
    bool is_hole_) {
  pos = pos_;
  pll = pll_;
  gd = gd_;
  pgnwh = pgnwh_;

  v1 = *pll_.begin();
  auto v2_it = pll_.begin();
  v2_it++;
  v2 = *v2_it;
  vl = *pll_.rbegin();
  auto v2l_it = pll_.rbegin();
  v2l_it++;
  v2l = *v2l_it;

  is_hole = is_hole_;
}

int PLL::get_pos() {
  return pos;
}

Polyline PLL::get_pll() {
  return pll;
}

Point PLL::get_v1() {
  return v1;
}

Point PLL::get_v2() {
  return v2;
}

Point PLL::get_vl() {
  return vl;
}

Point PLL::get_v2l() {
  return v2l;
}

int PLL::get_gd() {
  return gd;
}

int PLL::get_pgnwh() {
  return pgnwh;
}

bool PLL::get_is_hole() {
  return is_hole;
}

void PLL::set_pll(Polyline pll_) {
  pll = pll_;
  v1 = *pll_.begin();
  auto v2_it = pll_.begin();
  v2_it++;
  v2 = *v2_it;
  vl = *pll_.rbegin();
  auto v2l_it = pll_.rbegin();
  v2l_it++;
  v2l = *v2l_it;
}

void PLL::set_is_hole(bool b) {
  is_hole = b;
}
