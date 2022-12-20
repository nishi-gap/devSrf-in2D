#ifndef FOLDLINE_H
#define FOLDLINE_H
#include <setrulings.h>
#include <cmath>
#include <Eigen/Dense>
#include "mathtool.h"
class FoldLine
{
public:
    FoldLine(int crvNum, int rsize, int _type);
    std::vector<glm::f64vec3> getCtrlPt();
    std::vector<glm::f64vec3> getCtrlPt2d();
    bool addCtrlPt(glm::f64vec3& p, int dim, OUTLINE *outline, std::vector<Face*>& Faces, std::vector<HalfEdge*>&Edges, std::vector<Vertex*>& Vertices);
    bool delCtrlPt(glm::f64vec3& p, int dim, OUTLINE *outline);
    std::vector<glm::f64vec3> CurvePts;
    std::vector<std::array<glm::f64vec3, 2>> Rulings_3dL, Rulings_3dR, Rulings_2dL, Rulings_2dR;
    bool ChangeColor(OUTLINE *outline, int val, int dim = 3);
    double getColor();
    bool modify2DRulings(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices,const std::vector<HalfEdge*>& edge_outline, int dim);
    bool applyCurvedFolding(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, int dim);
    void deform();
    HalfEdge *he, *he2;
    std::vector<glm::f64vec3> point;
    Vertex *vx, *vx2;

    std::vector<glm::f64vec3> CtrlPts_res, Curve_res, CtrlPts_res2d, Curve_res2d;
    std::vector<CrvPt_FL> T_crs;
private:
    double color;
    bool setCurve(int dim);

    void devide(Vertex *v1, Vertex *v2, std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, EdgeType _type);
    bool BezierCrvOn3dSrf(std::vector<glm::f64vec3>& CtrlPts, double t, int dim, std::vector<Face*>& Faces, glm::f64vec3& v_3d);
    void devide2Faces(std::vector<HalfEdge*>& Inserted, std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces);
    std::vector<glm::f64vec3>CtrlPts;
    bool setPoint(const std::vector<glm::f64vec3>& edge_outline, glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint);
    void cal2VecScale(glm::f64vec3 v1, glm::f64vec3 v2, glm::f64vec3 p, double& s, double& t);
    int type;
    int maxRsize;

    void ProjectBezierOn3d(int dim);

    void setCoeff(std::vector<double>& a, double t, std::vector<double>& Knot, int j);

    void diff(double t, std::vector<double>& Knot, std::vector<glm::f64vec3>& dP, std::vector<glm::f64vec3>& CtrlPts, const int n_times = 3);
    inline double rad_2d(double k, double tau, double a, double da);

};

#endif // FOLDLINE_H
