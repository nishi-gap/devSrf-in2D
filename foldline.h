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
    bool moveCtrlPt(glm::f64vec3& p);
    std::vector<glm::f64vec3> CurvePts;
    std::vector<std::array<glm::f64vec3, 2>> Rulings_3dL, Rulings_3dR, Rulings_2dL, Rulings_2dR;
    bool ChangeColor(OUTLINE *outline, int val, int dim = 3);
    double getColor();
    bool modify2DRulings(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices,const std::vector<Vertex*>& Poly_v, int dim, int t_type = 2);

    void deform();
    HalfEdge *he, *he2;
    std::vector<glm::f64vec3> point;
    Vertex *vx, *vx2;

    std::vector<glm::f64vec3> CtrlPts_res, Curve_res, CtrlPts_res2d, Curve_res2d;
    //std::vector<CrvPt_FL> T_crs;
    std::vector<HalfEdge_FL*> FoldingCurve;
    std::vector<glm::f64vec3>BasisVectors, PointsOnPlane;

    double AngleIn2Edges(HalfEdge *p, HalfEdge *p2, bool Is3d = true);
    void applyAAAMethod(std::vector<Vertex*>& Poly_v,  std::vector<Face*>& Faces, std::vector<HalfEdge*>& edges, double a = 2.0*M_PI/3.0);
    void TestFoldingAAAM(double& a, std::vector<Vertex*>& _Vertices, std::vector<HalfEdge*>& _Edges);

    bool SplitFace4DebugAAAMethod(glm::f64vec3& NewPoint, std::vector<Face*> &faces, std::vector<HalfEdge*>& edges, std::vector<Vertex*>& vertices);

    std::vector<std::array<glm::f64vec3, 2>> SingleRuling, AllRulings, NewRuling2d;
    void drawRulingInAllAngles(std::vector<std::array<glm::f64vec3, 2>>& _Rulings);
private:
    double color;
    bool setCurve(int dim);

    void devide(Vertex *v1, Vertex *v2, std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, EdgeType _type);
    bool BezierCrvOn3dSrf(std::vector<glm::f64vec3>& CtrlPts, double t, int dim, std::vector<Face*>& Faces, glm::f64vec3& v_3d);
    void devide2Faces(std::vector<HalfEdge*>& Inserted, std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces);
    std::vector<glm::f64vec3>CtrlPts;
    bool setPoint(const std::vector<Vertex*>& Poly_v, glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint);
    PaintTool type;
    int maxRsize;

    void diff(double t, std::vector<double>& Knot, std::vector<glm::f64vec3>& dP, std::vector<glm::f64vec3>& CtrlPts, int index, const int n_times = 3);
    inline double rad_2d(double k, double tau, double a, double da);


    void _FoldingAAAMethod(double & a, double phi02, double phim1, std::vector<Vertex*>& Poly_v,  std::vector<Face*>& Faces, std::vector<HalfEdge*>& edges);
    void EdgeRecconection(std::vector<Vertex*>& Poly_V,  std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges);

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
