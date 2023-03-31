#ifndef FOLDLINE_H
#define FOLDLINE_H

#include <cmath>
#include <numbers>
#include <setrulings.h>
#include<glm/gtx/vector_angle.hpp>

#include <utility>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include "mathtool.h"
class FoldLine
{
public:
    FoldLine(int crvNum, int rsize, PaintTool _type);
    std::vector<glm::f64vec3> getCtrlPt();
    std::vector<glm::f64vec3> getCtrlPt2d();
    bool addCtrlPt(glm::f64vec3& p, int dim);
    bool delCtrlPt(glm::f64vec3& p, int dim, OUTLINE *outline);
    bool moveCtrlPt(glm::f64vec3& p, int movePtIndex);
    std::vector<glm::f64vec3> CurvePts;
    std::vector<std::array<glm::f64vec3, 2>> Rulings_3dL, Rulings_3dR, Rulings_2dL, Rulings_2dR;
    bool ChangeColor(OUTLINE *outline, int val, int dim = 3);
    double getColor();
    bool modify2DRulings(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices,const std::vector<Vertex*>& Poly_v, int dim, int t_type = 2);

    bool Optimization(std::vector<Vertex*>& Poly_V, std::vector<Face*>&Faces,  std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices);
    void deform();
    HalfEdge *he, *he2;
    std::vector<glm::f64vec3> point;
    Vertex *vx, *vx2;

    std::vector<glm::f64vec3> CtrlPts_res, Curve_res, CtrlPts_res2d, Curve_res2d;
    //std::vector<CrvPt_FL> T_crs;
    std::vector<std::pair<CrvPt_FL*, Vertex*>> FoldingCurve;
    //std::vector<HalfEdge_FL*> FoldingCurve;
    std::vector<glm::f64vec3>BasisVectors, PointsOnPlane;

    double AngleIn2Edges(HalfEdge *p, HalfEdge *p2, bool Is3d = true);
    void applyAAAMethod(std::vector<Vertex*>& Poly_v,  std::vector<Face*>& Faces, std::vector<HalfEdge*>& edges, std::vector<Vertex*>&Vertices, double a = 2.0*M_PI/3.0);
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
    double Fbend(std::vector<Vertex*>& Poly_V, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices);
    double Fruling(std::vector<Vertex*>& Poly_V);
    void _FoldingAAAMethod(double & a, double phi02, double phim1, std::vector<Vertex*>& Poly_v,  std::vector<Vertex*>&Vertices, std::vector<HalfEdge*>& edges);
};


#endif // FOLDLINE_H
