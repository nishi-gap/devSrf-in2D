#ifndef SETRULINGS_H
#define SETRULINGS_H

#include <iostream>
#include <vector>
#include <glm/glm.hpp>
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
class meshline;
class ruling;
class OUTLINE;

class Vertex{
public:
    glm::f64vec3 p;
    std::vector<HalfEdge*> halfedge;
    bool deformed;
    Vertex(glm::f64vec3 _p);

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
    HalfEdge(Vertex *v);
    ruling *r;
private:

};

class meshline
{
    public:
    Vertex *vertices[2];
    std::vector<std::tuple<Vertex*, Vertex*>> vertices4grad;//rulingは一組だが色の多角形領域を決めるverticesはv複数格納することでrulingが格納される多角形領域を決める際のsafty failを回避.サイズは必ず1として先頭の要素は可能性が最も高い物を残しておく
    ruling *hasRulings[2];
    crvpt *pt;
    meshline(Vertex *a, Vertex *b, crvpt *_pt);
    meshline();
};

class Face{
public:
    bool bend;
    //double Gradation;//範囲-255 ~ 255
    bool hasGradPt;
    //std::vector<ruling*> rulings;
    HalfEdge* halfedge;
    Face(HalfEdge *_halfedge);
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
    void Arc();//制御点 0,3,...: 原点. 1,4,...: 始点. 2,5,...: 終点
    void ArcRulings(OUTLINE *outline, int DivSize);

    void addCtrlPt(QPointF p);
    void eraseCtrlPt(int curveDimention, int crvPtNum);
    void movePt(QPointF p, int ind);
    int movePtIndex(QPointF p, double& dist);
    void ClearPt();
    void InsertControlPoint2(QPointF pt);
    void SetNewPoint();

    int InsertPointSegment;
    int getCurveType();
    void setCurveType(int n);
    std::vector<glm::f64vec3> ControllPoints;
    std::vector<crvpt> CurvePoints;
    std::vector<ruling*> Rulings;//偶数番目 ruling　奇数番目 グラデーションの多角形に使用
    //std::vector<meshline*> meshLines;
    glm::f64vec3 InsertPoint;
    bool isempty;

    //グラデーション
    void FillColor(int c);
    int getColor(int n);
    void setcolor(int ctype, int cval, int n);
    void Interpolation(int method, int g0, int g1);
    void clearColor();

    void CrossDetection();
private:

    bool IsInsertNewPoint;
    int OnCurvesORLines(QPointF& p, int& ind);//-1：どこにものっかっていない　0：曲線上　1：制御点を結んだ線上

    double basis(int j, int k, double t, std::vector<double>& T);
    double factorial(int n);
    double cmb(int n, int i);

    int curveType;
    bool setPoint(std::vector<Vertex*>&outline, glm::f64vec3 N, glm::f64vec3& cp, std::vector<glm::f64vec3>& P);
    glm::f64vec3 getIntersectionPoint(glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3, glm::f64vec3& p4);
    void swap(glm::f64vec3&a, glm::f64vec3& b);

    double crvStep;
    int curveNum;
    double meshLength(int s);
    void SplineInterpolation(std::vector<glm::f64vec2>& cp, std::vector<glm::f64vec2>& CurvePath);
    void LinearInterpolation(int g0, int g1);    

};

class OUTLINE{
public:

    OUTLINE();
    QString type;//Rectangle, Polygon, Polyline

    bool isClosed;
    int VerticesNum;

    void addVertex(QPointF p);
    void eraseVertex();
    std::vector<Vertex*> getVertices();
    void drawPolygon(QPointF p, bool IsClicked);
    void MoveVertex(QPointF p);//polygonの移動
    void EditVertex(QPointF p);
    glm::f64vec2 origin;//polygonの始点
    int hasPtNum; //0 ~ 2 (polygonの点の数)
private:
    std::vector<Vertex*> vertices;
};

double distP2L(glm::f64vec3 la, glm::f64vec3 lb, QPointF pt, glm::f64vec3& q);//点と線分の距離, s:laからlbへの比率(垂線が内部にあれば0 ~ 1)


// -1: 最近傍の点なし
int movePointIndex(QPointF p, std::vector<Vertex*>& V);
double QVdist(QPointF p, glm::f64vec2 v);
void CrossDetection(CRV& crvs);
bool cn(std::vector<glm::f64vec2> &V2, QPointF p);
bool cn(Face *face, glm::f64vec3 p);
double set3pt(glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3);
bool IsIntersect(glm::f64vec3&p1, glm::f64vec3&p2, glm::f64vec3&p3, glm::f64vec3&p4);
void Triangulation(std::vector<glm::f64vec3>&input, std::vector<std::vector<glm::f64vec3>>&output);
bool hasPointInTriangle3D(glm::f64vec3 p, std::vector<glm::f64vec3>& V);
bool hasPointInPolygon(glm::f64vec3 p, std::vector<glm::f64vec3>& V);
bool IsAngleLessThan180(glm::f64vec3& o, glm::f64vec3& a, glm::f64vec3& b);
std::vector<glm::f64vec3> TranslateGLMfromHE(Face *f);
glm::f64vec3 GetCenter(std::vector<glm::f64vec3>& vertices);
#endif // SETRULINGS_H
