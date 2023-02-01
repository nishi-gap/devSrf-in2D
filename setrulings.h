#ifndef SETRULINGS_H
#define SETRULINGS_H

#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include "mathtool.h"
#include <glm/gtx/transform.hpp>
#include <tuple>
#include <QDebug>
#include <QPointF>
#include <QString>
#include <Eigen/Dense>
#include <glm/gtx/string_cast.hpp>

constexpr glm::f64vec3 NullVec = glm::f64vec3{-1,-1,-1};

class HalfEdge;
class Model;
class Face;
class ruling;
class OUTLINE;

enum class EdgeType;
enum class CurveType;


class Vertex{
public:
    glm::f64vec3 p;
    glm::f64vec3 p3;
    bool deformed;
    Vertex(glm::f64vec3 _p);
    Vertex(glm::f64vec3 _p2, glm::f64vec3 _p3);
    void addNewEdge(HalfEdge *he);
    std::vector<HalfEdge*> halfedge;
};



class CrvPt_FL : public Vertex{
public:
    double s;
    //void digAlign(double& v, int dig = 11){ int n = pow(10, dig); int tmp = v * n; v = (double)tmp/(double)n; }

    double k2d, k3d, tau, k2d_m, k2d_p, k3d_m, k3d_p;
    double a, da;
    glm::f64vec3 T2d, N2d, B2d, T3d, N3d, B3d;
    glm::f64vec3 Td, Nd, Bd;
    CrvPt_FL(glm::f64vec3 _p2, glm::f64vec3 _p3, double _s) : Vertex(_p2, _p3), s{_s} {}
    bool operator<(const CrvPt_FL& T) const { return s < T.s; }
    double developability();

};

class crvpt{
public:
    double color;
    glm::f64vec3 pt;
    int ind;
    crvpt( int _ind,glm::f64vec3 _pt = NullVec, int _color = 0);
};

class ruling
{
    public:
    std::tuple<Vertex*, Vertex*> r;
    HalfEdge *he[2];
    int IsCrossed; //-1:交差なし, 0:同じruling上で交差, 1:上にあるレイヤー上のrulingと交差
    crvpt *pt;
    double Gradation;
    bool hasGradPt;
    ruling(Vertex *a, Vertex *b, crvpt *_pt);
    ruling();
};

class HalfEdge{
public:
    Vertex *vertex;
    Face *face;
    HalfEdge *prev;
    HalfEdge *next;
    HalfEdge *pair;

    HalfEdge(Vertex *v, EdgeType _type);
    ruling *r;
    EdgeType edgetype;
    std::vector<HalfEdge*> Split(Vertex *v, std::vector<HalfEdge*>& Edges);
private:

};

class Face{
public:
    int edgeNum();
    bool bend;
    //double Gradation;//範囲-255 ~ 255
    bool hasGradPt;
    //std::vector<ruling*> rulings;
    HalfEdge* halfedge;
    Face(HalfEdge *_halfedge);
    bool IsPointInFace(glm::f64vec3 p);
    glm::f64vec3 getNormalVec();
    double sgndist(glm::f64vec3 p);
};

class CRV{
public:
    CRV(int _crvNum, int DivSize);

    void Bspline(int curveDimention,  int crvPtNum);
    void Bezier(int curveDimention, int crvPtNum);
    void BezierRulings(OUTLINE *outline, int& DivSize, int crvPtNum);
    void BsplineRulings(OUTLINE *outline , int& DivSize, int crvPtNum, int curveDimention);
    void Line();
    void LineRulings(OUTLINE *outline, int DivSize);
    void Arc(int crvPtNum);//制御点 0,3,...: 原点. 1,4,...: 始点. 2,5,...: 終点
    void ArcRulings(OUTLINE *outline, int DivSize);

    void addCtrlPt(QPointF p);
    void eraseCtrlPt(int curveDimention, int crvPtNum);
    void movePt(glm::f64vec3 p, int ind);
    int movePtIndex(glm::f64vec3& p, double& dist);
    void ClearPt();
    void InsertControlPoint2(glm::f64vec3& p);
    void SetNewPoint();

    int InsertPointSegment;
    CurveType getCurveType();
    void setCurveType(CurveType n);
    std::vector<glm::f64vec3> ControllPoints;
    std::vector<crvpt> CurvePoints;
    std::vector<ruling*> Rulings;//偶数番目 ruling　奇数番目 グラデーションの多角形に使用
    //std::vector<meshline*> meshLines;
    glm::f64vec3 InsertPoint;
    bool isempty;

    //グラデーション
    void FillColor(int c);
    inline int getColor(int n);
    void setcolor(int ctype, int cval, int n);
    inline void clearColor();
private:

    bool IsInsertNewPoint;
    int OnCurvesORLines(glm::f64vec3& p, int& ind);//-1：どこにものっかっていない　0：曲線上　1：制御点を結んだ線上
    CurveType curveType;
    bool setPoint(std::vector<Vertex*>&outline, glm::f64vec3 N, glm::f64vec3& cp, std::vector<glm::f64vec3>& P);
    inline void swap(glm::f64vec3&a, glm::f64vec3& b);

    double crvStep;
    int curveNum;

};

class OUTLINE{
public:

    OUTLINE();
    QString type;//Rectangle, Polygon, Polyline

    bool IsClosed();
    int VerticesNum;
    void addVertex(Vertex*v, int n);
    void addVertex(glm::f64vec3& p);
    void eraseVertex();
    std::vector<Vertex*> getVertices();
    std::vector<HalfEdge*> getEdges();
    Face* getFace();
    void ConnectEdges(bool IsConnected = true);
    void drawPolygon(glm::f64vec3& p, bool IsClicked);
    void MoveOutline(glm::f64vec3 p);//polygonの移動
    void MoveVertex(glm::f64vec3 p, int ind);
    glm::f64vec2 origin;//polygonの始点
    int hasPtNum; //0 ~ 2 (polygonの点の数)
private:
    std::vector<Vertex*> vertices;
    std::vector<HalfEdge*> edges;
    int movePointIndex(glm::f64vec3 p);
    Face *face;
};

void CrossDetection(Face *f, CRV *crvs);


std::vector<double> BezierClipping(std::vector<glm::f64vec3>&CtrlPts, HalfEdge *line, int dim);
std::vector<glm::f64vec3> ConvertDistBasedBezier(std::vector<glm::f64vec3>& CtrlPts, HalfEdge *line);

std::vector<glm::f64vec3> GlobalSplineInterpolation(std::vector<CrvPt_FL>& Q, std::vector<glm::f64vec3>& CtrlPts_res, std::vector<double>& Knot, double& CurveLen, bool is3d = true, int dim = 3, int t_type = 2);

std::vector<glm::f64vec3> TBCSplineInterpolation(std::vector<CrvPt_FL>& Q, double& CurveLen, bool is3d = true, double ten = 0, double bias = 0, double con = 0);


#endif // SETRULINGS_H
