#ifndef MATHTOOL_H
#define MATHTOOL_H

#include <cmath>
#include <vector>
#include <iostream>
#include <tuple>
#include <utility>
#include <algorithm>
#include <numbers>
#include <queue>
#include <memory>
#include <Eigen/Dense>

enum class EdgeType{
    none,
    ol,//outline
    r,//ruling(from curve line)
    cl,//curve line
    fl,//fold line
    r_fl,//ruling(fold line)
};

enum class CurveType{
    none,
    bezier3,
    bsp3,
    line,
    arc,
};

struct ColorPoint{
    double color, angle;
    ColorPoint(): color(200), angle(std::numbers::pi/2.0){}
    ColorPoint(double _c, double _a): color(_c), angle(_a){}
};

enum class PaintTool{
    None,
    Reset,
    deform,
    AffinTrans,

    //ruling control curve
    AddCurve,
    Bspline,
    Line,
    Arc,

    SetColor,
    NewGradationMode,

    Rectangle,
    Polygon,
    Polyline,
    EditVertex,
    Move_ol,

    MoveCtrlPt,
    InsertCtrlPt,
    DeleteCtrlPt,
    DeleteCurve,

    Crease,
    Crease_move,

};

namespace DebugMode{
    class Singleton{
    public:
        static Singleton& getInstance();
        bool isdebug() const;
        void switchval();
        int valuecanged(int n);
        void init();
    private:
        Singleton(){}
        Singleton(const Singleton&) = delete;
        void operator=(const Singleton&) = delete;
        bool val;
    };
}

namespace MathTool{
double rad2deg(double a);

bool IsIntersect(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, const Eigen::Vector3d& p4, Eigen::Vector3d& q, bool ConsiderEnd = false);
Eigen::Vector3d calcCrossPoint_2Vector(Eigen::Vector3d p1, Eigen::Vector3d q1, Eigen::Vector3d p2, Eigen::Vector3d q2);
std::vector<double> _bezierclipping(const std::vector<Eigen::Vector3d>&CtrlPts_base, std::vector<Eigen::Vector3d>&CtrlPts_cur, std::array<Eigen::Vector3d, 2>& line, int dim);//交点が一つのみの場合
bool is_point_on_line(Eigen::Vector3d p, Eigen::Vector3d lp1, Eigen::Vector3d lp2);
std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> BezierSplit(const std::vector<Eigen::Vector3d>& CtrlPts, double t);
std::vector<std::vector<Eigen::Vector3d>> de_casteljau_algorithm(const std::vector<Eigen::Vector3d>& CtrlPts, double t);

std::vector<Eigen::Vector3d> GrahamScan(const std::vector<Eigen::Vector3d>& Q);//凸包の計算
double SignedArea(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& p);

Eigen::Vector3d bspline(std::vector<Eigen::Vector3d>&CtrlPts, double t, int dim, std::vector<double>Knot);

double factorial(int n);
double cmb(int n, int i);
double BernsteinBasisFunc(int n, int i, double t);
Eigen::Vector3d bezier(std::vector<Eigen::Vector3d>& CtrlPts, double t, int dim);

double basis(int n, int i, int p, double u, std::vector<double>& U);

Eigen::Vector3d ProjectionVector(const Eigen::Vector3d& v, Eigen::Vector3d n, bool Isnormalize = false);


template<typename T>
void swap(T &a, T& b){ T c = a; a = b; b = c;}

}

#endif // MATHTOOL_H
