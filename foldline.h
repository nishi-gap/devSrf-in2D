#ifndef FOLDLINE_H
#define FOLDLINE_H
#include <setrulings.h>

#include <cmath>


class FoldLine
{
public:
    FoldLine(int crvNum, int rsize, int _type);
    std::vector<glm::f64vec3> getCtrlPt();
    bool addCtrlPt(glm::f64vec3& p, int dim, OUTLINE *outline, std::vector<Face*>& Faces, std::vector<HalfEdge*>&Edges, std::vector<Vertex*>& Vertices);
    bool delCtrlPt(glm::f64vec3& p, int dim, OUTLINE *outline);
    std::vector<glm::f64vec3> CurvePts;
    std::vector<std::array<glm::f64vec3, 2>> Rulings_3dL, Rulings_3dR, Rulings_2dL, Rulings_2dR;
    bool ChangeColor(OUTLINE *outline, int val, int dim = 3);
    double getColor();
    bool modify2DRulings(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices,int dim);
    bool applyCurvedFolding(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, int dim);
    void deform();
    HalfEdge *he, *he2;
    std::vector<glm::f64vec3> point;
    Vertex *vx, *vx2;

    std::vector<glm::f64vec3> CtrlPts_res, Curve_res;
    std::vector<Vertex*> CrossPts;
private:
    double color;
    bool setCurve(int dim);

    void devide(Vertex *v1, Vertex *v2, std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, EdgeType _type);
    bool BezierCrvOn3dSrf(std::vector<glm::f64vec3>& CtrlPts, double t, int dim, std::vector<Face*>& Faces, glm::f64vec3& v_3d);
    std::array<HalfEdge*, 2> splitOnPoint(Vertex *v, std::vector<HalfEdge*>& Edges);
    void devide2Faces(std::vector<HalfEdge*>& Inserted, std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces);
    std::vector<glm::f64vec3>CtrlPts;
    bool setPoint(std::vector<Face*>& Faces, glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint);
    void cal2VecScale(glm::f64vec3 v1, glm::f64vec3 v2, glm::f64vec3 p, double& s, double& t);
    int type;
    int maxRsize;


};

#endif // FOLDLINE_H
