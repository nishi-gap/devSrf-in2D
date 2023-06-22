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

#include <nlopt.hpp>



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
    double getColor();
    bool modify2DRulings(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_v, int dim, int t_type = 2);
    bool RevisionCrosPtsPosition(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, int type, bool TrimMode);
    bool Optimization_FlapAngle(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, int type = 0, bool IsValidSigmoid = false);
    bool Optimization_SmooothSrf(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, int type = 0);
    void Optimization2(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_v, double a);
    void Optimization_Vertices(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, int type, bool IsValidSigmoid);
    void ReassignColor(std::vector<HalfEdge*>& Edges, ColorPoint& CP);
    void SimplifyModel(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, double tol);

    //HalfEdge *he, *he2;
    std::vector<glm::f64vec3> point;
    Vertex *vx, *vx2;

    std::vector<glm::f64vec3> CtrlPts_res, Curve_res, CtrlPts_res2d, Curve_res2d;
    //std::vector<CrvPt_FL> T_crs;
    std::vector<Vertex4d> FoldingCurve;
    //std::vector<HalfEdge_FL*> FoldingCurve;
    void applyAAAMethod(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, double a, int type = 0);
    std::vector<std::array<glm::f64vec3, 2>>  AllRulings, NewRuling2d;
    void drawRulingInAllAngles(std::vector<std::array<glm::f64vec3, 2>>& _Rulings);

private:
    double color;
    bool setCurve(int dim);
    void devide(Vertex *v1, Vertex *v2, std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, EdgeType _type);
    bool BezierCrvOn3dSrf(std::vector<glm::f64vec3>& CtrlPts, double t, int dim, std::vector<Face*>& Faces, glm::f64vec3& v_3d);
    void devide2Faces(std::vector<HalfEdge*>& Inserted, std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces);
    std::vector<glm::f64vec3>CtrlPts;

    void TrimPoints(std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces,std::vector<Vertex*>& Vertices, std::vector<CrvPt_FL*>& T_crs, double tol = 0.2);//kがpiに近い値であるかどうかの検出
    PaintTool type;
    int maxRsize;

    std::vector<CrvPt_FL*> Points_On_Curve;

};


#endif // FOLDLINE_H
