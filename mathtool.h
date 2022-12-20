#ifndef MATHTOOL_H
#define MATHTOOL_H

#include <glm/gtx/transform.hpp>
#include <vector>
#include <iostream>
#include <tuple>
#include <Eigen/Dense>
#include <utility>


namespace MathTool{
double distP2L(glm::f64vec3 la, glm::f64vec3 lb, glm::f64vec3& p, glm::f64vec3& q);//点と線分の距離, s:laからlbへの比率(垂線が内部にあれば0 ~ 1)

double set3pt(glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3);
bool IsIntersect(glm::f64vec3&p1, glm::f64vec3&p2, glm::f64vec3&p3, glm::f64vec3&p4);
glm::f64vec3 getIntersectionPoint(glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3, glm::f64vec3& p4);

void Triangulation(std::vector<glm::f64vec3>&input, std::vector<std::array<glm::f64vec3, 3>>&output);
bool hasPointInTriangle3D(glm::f64vec3 p, std::array<glm::f64vec3, 3>& V);
bool IsAngleLessThan180(glm::f64vec3& o, glm::f64vec3& a, glm::f64vec3& b);

std::vector<double> _bezierclipping(std::vector<glm::f64vec3>&CtrlPts_base, std::vector<glm::f64vec3>&CtrlPts_cur, std::array<glm::f64vec3, 2>& line, int dim);//交点が一つのみの場合
bool is_point_on_line(glm::f64vec3& p, glm::f64vec3& lp1, glm::f64vec3& lp2);
std::pair<std::vector<glm::f64vec3>, std::vector<glm::f64vec3>> BezierSplit(std::vector<glm::f64vec3> CtrlPts, double t, int dim);
std::vector<std::vector<glm::f64vec3>> de_casteljau_algorithm(std::vector<glm::f64vec3> CtrlPts, double t);

std::vector<glm::f64vec3> GrahamScan(std::vector<glm::f64vec3>& Q);//凸包の計算
double SignedArea(glm::f64vec3 a, glm::f64vec3 b, glm::f64vec3 p);

glm::f64vec3 bspline(std::vector<glm::f64vec3>&CtrlPts, double t, int dim, std::vector<double>Knot);
double factorial(int n);
double cmb(int n, int i);
double basis(int i, int p, double u, std::vector<double>& U);

}

#endif // MATHTOOL_H
