#ifndef SETRULINGS_H
#define SETRULINGS_H

#include <iostream>
#include <vector>
//#include <glm/glm.hpp>
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
class meshline;
class ruling;
class OUTLINE;

enum class EdgeType{
    ol,//outline
    r,//ruling(from curve line)
    cl,//curve line
    fl,//fold line
    r_fl,//ruling(fold line)
};

enum class PaintTool{
    None,
    Reset,
    deform,

    AddCurve,
    Bezier_r,
    Bspline_r,
    Line_r,
    Arc_r,

    SetColor,
    NewGradationMode,

    Rectangle_ol,
    Polygon_ol,
    Polyline_ol,
    EditVertex_ol,
    Move_ol,
    Const_ol,
    ConnectVertices_ol,

    MoveCtrlPt,
    InsertCtrlPt,
    DeleteCtrlPt,
    DeleteCurve,

    FoldLine,
    FoldlineColor,
};

class Vertex{
public:
    glm::f64vec3 p;
    glm::f64vec3 p3;
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

    HalfEdge(Vertex *v, EdgeType _type);
    ruling *r;
    EdgeType edgetype;
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
private:

    bool IsInsertNewPoint;
    int OnCurvesORLines(glm::f64vec3& p, int& ind);//-1：どこにものっかっていない　0：曲線上　1：制御点を結んだ線上

    int curveType;
    bool setPoint(std::vector<Vertex*>&outline, glm::f64vec3 N, glm::f64vec3& cp, std::vector<glm::f64vec3>& P);
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

    //bool isClosed;
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
    Face *face;
};

double distP2L(glm::f64vec3 la, glm::f64vec3 lb, glm::f64vec3& p, glm::f64vec3& q);//点と線分の距離, s:laからlbへの比率(垂線が内部にあれば0 ~ 1)


// -1: 最近傍の点なし
int movePointIndex(glm::f64vec3 p, std::vector<Vertex*>& V);

double QVdist(QPointF p, glm::f64vec2 v);
void CrossDetection(Face *f, CRV *crvs);
bool cn(std::vector<glm::f64vec2> &V2, QPointF p);
bool cn(Face *face, glm::f64vec3 p);
double set3pt(glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3);
bool IsIntersect(glm::f64vec3&p1, glm::f64vec3&p2, glm::f64vec3&p3, glm::f64vec3&p4);
glm::f64vec3 getIntersectionPoint(glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3, glm::f64vec3& p4);

void Triangulation(std::vector<glm::f64vec3>&input, std::vector<std::array<glm::f64vec3, 3>>&output);
void Triangulation(Face *f, std::vector<std::array<glm::f64vec3, 3>>& output);
bool hasPointInTriangle3D(glm::f64vec3 p, std::array<glm::f64vec3, 3>& V);
bool hasPointInTriangle3D(glm::f64vec3 p, std::array<HalfEdge*, 3>& V);
bool hasPointInPolygon(glm::f64vec3 p, std::vector<glm::f64vec3>& V);
bool IsAngleLessThan180(glm::f64vec3& o, glm::f64vec3& a, glm::f64vec3& b);

std::vector<double> BezierClipping(std::vector<glm::f64vec3>&CtrlPts, HalfEdge *line, int dim);
std::vector<double> _bezierclipping(std::vector<glm::f64vec3>&CtrlPts_base, std::vector<glm::f64vec3>&CtrlPts_cur, std::array<glm::f64vec3, 2>& line, int dim);//交点が一つのみの場合
bool is_point_on_line(glm::f64vec3& p, glm::f64vec3& lp1, glm::f64vec3& lp2);
std::pair<std::vector<glm::f64vec3>, std::vector<glm::f64vec3>> BezierSplit(std::vector<glm::f64vec3> CtrlPts, double t, int dim);
std::vector<glm::f64vec3> ConvertDistBasedBezier(std::vector<glm::f64vec3>& CtrlPts, HalfEdge *line);
Eigen::MatrixXd create_Ux(int dim);
std::vector<std::vector<glm::f64vec3>> de_casteljau_algorithm(std::vector<glm::f64vec3> CtrlPts, double t);

std::vector<glm::f64vec3> GrahamScan(std::vector<glm::f64vec3>& Q);//凸包の計算
double SignedArea(glm::f64vec3 a, glm::f64vec3 b, glm::f64vec3 p);

std::vector<glm::f64vec3> TranslateGLMfromHE(Face *f);
glm::f64vec3 GetCenter(std::vector<glm::f64vec3>& vertices);
double basis(int j, int k, double t, std::vector<double>& T);
glm::f64vec3 bspline(std::vector<glm::f64vec3>&CtrlPts, double t, int dim, std::vector<double>Knot);
double factorial(int n);
double cmb(int n, int i);

std::vector<glm::f64vec3> GlobalSplineInterpolation(std::vector<glm::f64vec3>& Q, std::vector<glm::f64vec3>& CtrlPts_res);
#endif // SETRULINGS_H
