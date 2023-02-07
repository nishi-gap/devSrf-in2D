#ifndef FOLDLINE_H
#define FOLDLINE_H
#include <setrulings.h>
#include<glm/gtx/vector_angle.hpp>
#include <cmath>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <numbers>
#include "mathtool.h"
class FoldLine
{
public:
    FoldLine(int crvNum, int rsize, PaintTool _type);
    std::vector<glm::f64vec3> getCtrlPt();
    std::vector<glm::f64vec3> getCtrlPt2d();
    bool addCtrlPt(glm::f64vec3& p, int dim);
    bool delCtrlPt(glm::f64vec3& p, int dim, OUTLINE *outline);
    std::vector<glm::f64vec3> CurvePts;
    std::vector<std::array<glm::f64vec3, 2>> Rulings_3dL, Rulings_3dR, Rulings_2dL, Rulings_2dR;
    bool ChangeColor(OUTLINE *outline, int val, int dim = 3);
    double getColor();
    bool modify2DRulings(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices,const std::vector<HalfEdge*>& edge_outline, int dim, int t_type = 2);
    bool applyCurvedFolding(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, int dim);
    void deform();
    HalfEdge *he, *he2;
    std::vector<glm::f64vec3> point;
    Vertex *vx, *vx2;

    std::vector<glm::f64vec3> CtrlPts_res, Curve_res, CtrlPts_res2d, Curve_res2d;
    std::vector<CrvPt_FL> T_crs;
    std::vector<HalfEdge*> FoldingCurve;
    std::vector<glm::f64vec3>BasisVectors, PointsOnPlane;

    double AngleIn2Edges(HalfEdge *p, HalfEdge *p2, bool Is3d = true);
    void applyAAAMethod(const std::vector<glm::f64vec3>& edge_outline, double a = 2.0*std::numbers::pi/3.0);
private:
    double color;
    bool setCurve(int dim);

    void devide(Vertex *v1, Vertex *v2, std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, EdgeType _type);
    bool BezierCrvOn3dSrf(std::vector<glm::f64vec3>& CtrlPts, double t, int dim, std::vector<Face*>& Faces, glm::f64vec3& v_3d);
    void devide2Faces(std::vector<HalfEdge*>& Inserted, std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces);
    std::vector<glm::f64vec3>CtrlPts;
    bool setPoint(const std::vector<glm::f64vec3>& edge_outline, glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint);
    void cal2VecScale(glm::f64vec3 v1, glm::f64vec3 v2, glm::f64vec3 p, double& s, double& t);
    PaintTool type;
    int maxRsize;

    void ProjectBezierOn3d(int dim);

    void setCoeff(std::vector<double>& a, double t, std::vector<double>& Knot, int j);

    void diff(double t, std::vector<double>& Knot, std::vector<glm::f64vec3>& dP, std::vector<glm::f64vec3>& CtrlPts, int index, const int n_times = 3);
    inline double rad_2d(double k, double tau, double a, double da);

    void LSM_apply_tau(int dim = 1);



    void _FoldingAAAMethod(double & a, double phi02, double phim1, const std::vector<glm::f64vec3>& edge_outline);
    double getCurvetureFromDC(int i);

};

class optimizer{
    double Fvert(std::vector<double>& X);
    double Ffit();
    double Ffair();
    double Fconv();
    void alignedNewRulingDirection();
    std::vector<Face*>Faces;//original
    std::vector<Vertex*>Vertices;
    std::vector<Vertex*> FoldingCurve;
    std::vector<Face*> ptchDevSrf;//after optimization
    HalfEdge* getEdge(Vertex* start, Vertex* end, const std::vector<Face*>& searchedFace);
public:

    optimizer(std::vector<Face*>& _Faces, std::vector<Vertex*>& _Vertices, std::vector<Vertex*>& _FoldingCurve) : Faces{_Faces}, Vertices{_Vertices}, FoldingCurve{_FoldingCurve} {}
    void initialize();
    void apply(double wfit = 1, double wfair = 1e-4);
};
#endif // FOLDLINE_H
