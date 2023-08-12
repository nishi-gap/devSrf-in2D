#ifndef FOLDLINE_H
#define FOLDLINE_H

#include <cmath>
#include <numbers>
#include "mathtool.h"
#include <setrulings.h>
#include<glm/gtx/vector_angle.hpp>

#include <utility>
#include <Eigen/Dense>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>



#include <nlopt.hpp>



class FoldLine
{
public:
    FoldLine(PaintTool _type);
    std::vector<glm::f64vec3>CtrlPts;
    bool addCtrlPt(glm::f64vec3& p, int dim);
    bool delCtrlPt(glm::f64vec3& p, int dim, OUTLINE *outline);
    bool moveCtrlPt(glm::f64vec3& p, int movePtIndex);
    std::vector<glm::f64vec3> CurvePts;
    std::vector<std::array<glm::f64vec3, 2>> Rulings_3dL, Rulings_3dR, Rulings_2dL, Rulings_2dR;
    double getColor();
    bool RevisionCrosPtsPosition();
    bool Optimization_FlapAngle(std::vector<Line*>& Rulings, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, double wb, double wp, bool ConstFunc = true);
    std::vector<std::vector<glm::f64vec3>> Optimization_SmooothSrf(const std::vector<Vertex*>& Poly_v, bool IsConnectEndPoint);
    std::vector<std::vector<glm::f64vec3>> Optimization_PlanaritySrf(const std::vector<Vertex*>& Poly_v);
    void ReassignColor(std::vector<Line*>& Rulings, ColorPoint& CP);
    void SimplifyModel( double tol);
    bool SimpleSmooothSrf(const std::vector<Vertex*>& Poly_v);
    void modifyFoldingCurvePositionOn3d(const std::vector<Line*>& Rulings);

    //HalfEdge *he, *he2;
    std::vector<glm::f64vec3> point;
    Vertex *vx, *vx2;

    std::vector<Vertex4d> FoldingCurve;
    void applyAAAMethod(std::vector<Vertex*>& Poly_V, double a);

    std::vector<std::array<glm::f64vec3, 2>>  AllRulings, NewRuling2d;
    void drawRulingInAllAngles(std::vector<std::array<glm::f64vec3, 2>>& _Rulings);

private:
    double color;
    bool setCurve(int dim);

    PaintTool type;

    std::vector<CrvPt_FL*> Points_On_Curve;

};


#endif // FOLDLINE_H
