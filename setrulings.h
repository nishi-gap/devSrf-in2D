#ifndef SETRULINGS_H
#define SETRULINGS_H

#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include "mathtool.h"
#include <glm/gtx/transform.hpp>
#include <tuple>
#include <QDebug>
#include <QPointF>
#include <QString>
#include <Eigen/Dense>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/vector_angle.hpp>

constexpr glm::f64vec3 NullVec = glm::f64vec3{-1,-1,-1};


//class ruling;
class OUTLINE;
class Vertex;
class HalfEdge;
class Face;
enum class EdgeType;
enum class CurveType;


class Vertex{
public:
    glm::f64vec3 p;
    glm::f64vec3 p3;
    glm::f64vec3 p_test;
    std::vector<HalfEdge*> halfedge;
    bool deformed;
    Vertex(glm::f64vec3 _p);
    Vertex(glm::f64vec3 _p2, glm::f64vec3 _p3);
    Vertex(const Vertex* v);
    ~Vertex();

    void addNewEdge(HalfEdge *he);
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
    ruling(Vertex *a, Vertex *b, crvpt *_pt = nullptr);

    ruling();
};

class HalfEdge{
public:
    int IsCrossed;
    Vertex *vertex;
    Face *face;
    HalfEdge *prev;
    HalfEdge *next;
    HalfEdge *pair;

    HalfEdge(Vertex *v, EdgeType _type);
    HalfEdge(const HalfEdge* he);
    ~HalfEdge();
    ruling *r;
    EdgeType edgetype;
    std::vector<HalfEdge*> Split(Vertex *v, std::vector<HalfEdge*>& Edges);
    bool hasCrossPoint2d(glm::f64vec3 p, glm::f64vec3 q, glm::f64vec3& CrossPoint,  bool ConsiderEnd = false);
    double diffEdgeLength();
    HalfEdge* erase(std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces);

protected:
    void edgeSwap(HalfEdge *h);
private:
};

class Face{
public:
    int edgeNum(bool PrintVertex = false);
    bool bend;
    bool hasGradPt;
    HalfEdge* halfedge;

    Face(HalfEdge *_halfedge);
    Face(){}
    Face(const Face& face);
    ~Face(){}

    bool IsPointInFace(glm::f64vec3 p);
    glm::f64vec3 getNormalVec();
    double sgndist(glm::f64vec3 p);
    void ReConnect(HalfEdge *he);
    void TrianglationSplit(std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces);

};


class CrvPt_FL : public Vertex{
public:
    double s, rt;
    Vertex *ve, *vo;
    CrvPt_FL(glm::f64vec3 _p2, glm::f64vec3 _p3,  double _s) : Vertex(_p2, _p3), s{_s} {}
    CrvPt_FL(glm::f64vec3 _p2, double _s) : Vertex(_p2), s{_s} {}
    void set(glm::f64vec3 _p,Vertex *o, Vertex *e);
    double developability();

};

struct PointOnLine{
public:
   double t;
   Vertex *v;
   PointOnLine(double _t, Vertex *_v): t(_t), v(_v){}
   bool operator<(const PointOnLine& P) const { return t < P.t; }
};

class HalfEdge_FL: public HalfEdge{
public:
    CrvPt_FL *v;
    HalfEdge_FL(CrvPt_FL *_v, EdgeType _type): HalfEdge(_v, _type), v(_v){}
    HalfEdge_FL(HalfEdge *h): HalfEdge(h->vertex, h->edgetype){}
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

    //void addCtrlPt(QPointF p);
    void eraseCtrlPt(int curveDimention, int crvPtNum);
    //void movePt(glm::f64vec3 p, int ind);
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
void EdgeRecconection(const std::vector<Vertex*>& Poly_V, std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges);
std::vector<HalfEdge*> EdgeCopy(const std::vector<HalfEdge*> Edges, const std::vector<Vertex*> V);
//std::vector<glm::f64vec3> GlobalSplineInterpolation(std::vector<CrvPt_FL>& Q, std::vector<glm::f64vec3>& CtrlPts_res, std::vector<double>& Knot, double& CurveLen, bool is3d = true, int dim = 3, int t_type = 2);


#endif // SETRULINGS_H
