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
#include <memory>

//class ruling;
class OUTLINE;
class Vertex;
class Line;

class Vertex{
public:
    glm::f64vec3 p;
    glm::f64vec3 p3;
    glm::f64vec3 p3_ori, p2_ori;
    bool deformed;
    Vertex(glm::f64vec3 _p, bool _deformed = false);
    Vertex(glm::f64vec3 _p2, glm::f64vec3 _p3, bool _deformed = false);
    std::shared_ptr<Vertex> deepCopy();
    //~Vertex();
    bool operator != (const Vertex &V)const{return p != V.p || p2_ori != V.p2_ori|| p3 != V.p3 || p3_ori != V.p3_ori || deformed != V.deformed;}
    bool operator == (const Vertex &V)const{return p == V.p && p2_ori == V.p2_ori && p3 == V.p3 && p3_ori == V.p3_ori && deformed == V.deformed;}
};

class CrvPt_FL : public Vertex{
public:
    double s, rt;
    std::shared_ptr<Vertex> ve;
    std::shared_ptr<Vertex> vo;
    bool IsValid;
    CrvPt_FL(glm::f64vec3 _p2, glm::f64vec3 _p3,  double _s) : Vertex(_p2, _p3), s(_s), IsValid(true){}
    CrvPt_FL(glm::f64vec3 _p2, double _s) : Vertex(_p2), s(_s), IsValid(true){}
    bool operator == (const CrvPt_FL &p)const{return s == p.s && rt == p.rt && ve == p.ve && vo == p.vo && IsValid == p.IsValid;}
    bool operator != (const CrvPt_FL &p)const{return !(s == p.s && rt == p.rt && ve == p.ve && vo == p.vo && IsValid == p.IsValid);}
    void set(glm::f64vec3 _p,const std::shared_ptr<Vertex>& o, const std::shared_ptr<Vertex>& e);

};

struct Vertex4d{
    bool IsCalc;
    std::shared_ptr<CrvPt_FL> first;
    std::shared_ptr<Vertex> second;
    std::shared_ptr<Vertex> third;
    Vertex4d(std::shared_ptr<CrvPt_FL>& v, std::shared_ptr<Vertex>& v2, std::shared_ptr<Vertex>& v3);
    Vertex4d(const Vertex4d& V4d);
    Vertex4d();
    void release();
    bool operator == (const Vertex4d &V4d)const{return first == V4d.first && second == V4d.second && third == V4d.third && IsCalc == V4d.IsCalc;}
    bool operator != (const Vertex4d &V4d)const{ return first != V4d.first || second != V4d.second || third != V4d.third || IsCalc != V4d.IsCalc; }
    bool operator == (const std::shared_ptr<Vertex> &V)const{return first == V; }
    bool operator != (const std::shared_ptr<Vertex> &V)const{return first != V; }
    void operator = (const Vertex4d &V){first = V.first; second = V.second; third = V.third; IsCalc = V.IsCalc;}
};


class Line{
public:
    int IsCrossed; //-1:交差なし, 0:同じruling上で交差, 1:上にあるレイヤー上のrulingと交差
    double color;
    Vertex *o, *v;//原則の方向はv[1] - v[0] (v[0]が原点)
    EdgeType et = EdgeType::none;
    Line(Vertex *_o, Vertex *_v, EdgeType _et): o(_o), v(_v), et(_et), IsCrossed(-1), color(0){}
    Line():o(nullptr),v(nullptr), IsCrossed(-1), color(0) {}
    bool operator !=(const Line &l)const{return IsCrossed != l.IsCrossed || color != l.color || (v[0] != l.v[0] && v[0] != l.v[1]);}
    bool operator ==(const Line &l)const{return IsCrossed == l.IsCrossed && color == l.color && ((v[0] == l.v[0] && v[1] == l.v[1]) || (v[1] == l.v[0] && v[0] == l.v[1]));}
    bool is_on_line(glm::f64vec3 p);
};

struct PointOnLine{
public:
   double t;
   Vertex *v;
   PointOnLine(double _t, Vertex *_v): t(_t), v(_v){}
   bool operator<(const PointOnLine& P) const { return t < P.t; }
};


class CRV{
public:
    CRV(int _crvNum, int DivSize);

    bool drawBspline(int curveDimention,  int crvPtNum);
    bool drawBezier(int curveDimention, int crvPtNum);
    void BezierRulings(OUTLINE *outline, int& DivSize, int crvPtNum);
    void BsplineRulings(OUTLINE *outline , int& DivSize, int crvPtNum, int curveDimention);
    bool drawLine();
    void LineRulings(OUTLINE *outline, int DivSize);
    bool drawArc(int crvPtNum);//制御点 0,3,...: 原点. 1,4,...: 始点. 2,5,...: 終点
    void ArcRulings(OUTLINE *outline, int DivSize);

    //void addCtrlPt(QPointF p);
    bool eraseCtrlPt(int curveDimention, int crvPtNum);
    //void movePt(glm::f64vec3 p, int ind);
    int movePtIndex(glm::f64vec3& p, double& dist);
    void ClearPt();
    void InsertControlPoint2(glm::f64vec3& p);
    void SetNewPoint();

    int InsertPointSegment;
    CurveType getCurveType();
    void setCurveType(CurveType n);
    std::vector<glm::f64vec3> ControllPoints;
    std::vector<glm::f64vec3> CurvePoints;
    std::vector<Line*> Rulings;//偶数番目 ruling　奇数番目 グラデーションの多角形に使用
    //std::vector<meshline*> meshLines;
    glm::f64vec3 InsertPoint;
    bool isempty;

    //グラデーション
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
    std::vector<Line*> Lines;
    void ConnectEdges(bool IsConnected = true);
    void drawPolygon(glm::f64vec3& p, bool IsClicked);
    void MoveOutline(glm::f64vec3 p);//polygonの移動
    void MoveVertex(glm::f64vec3 p, int ind);
    glm::f64vec2 origin;//polygonの始点
    int hasPtNum; //0 ~ 2 (polygonの点の数)
    bool IsPointInFace(glm::f64vec3 p);
    glm::f64vec3 getNormalVec();
private:
    std::vector<Vertex*> vertices;
    //std::vector<HalfEdge*> edges;

    int movePointIndex(glm::f64vec3 p);
    //Face *face;
};

void CrossDetection(OUTLINE *outline, CRV *crvs);


std::vector<double> BezierClipping(std::vector<glm::f64vec3>&CtrlPts, const std::shared_ptr<Vertex>& p, const std::shared_ptr<Vertex>& q, int dim);
std::vector<Vertex*> SortPolygon(std::vector<Vertex*>& polygon);
//std::vector<glm::f64vec3> GlobalSplineInterpolation(std::vector<CrvPt_FL>& Q, std::vector<glm::f64vec3>& CtrlPts_res, std::vector<double>& Knot, double& CurveLen, bool is3d = true, int dim = 3, int t_type = 2);


#endif // SETRULINGS_H
