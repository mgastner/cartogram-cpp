#ifndef DENSIFICATION_POINTS_H_
#define DENSIFICATION_POINTS_H_

std::vector<Point> densification_points(Point a, Point b);

bool point_almost_equal(Point a, Point b);
Point round_point(Point p1);

double half_ceil(double num);
double half_floor(double num);

#endif
