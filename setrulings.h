#ifndef SETRULINGS_H
#define SETRULINGS_H

#include <Eigen/Geometry> //EigenのGeometry関連の関数を使う場合，これが必要
#include "mathtool.h"
#include <QPointF>
#include <QString>

//class ruling;
class OUTLINE;
class Vertex;
class Line;

class Vertex{
public:
    Eigen::Vector3d p;
    Eigen::Vector3d p3;
    Eigen::Vector3d p3_ori, p2_ori;
    bool deformed;
    Vertex(Eigen::Vector3d _p, bool _deformed = false);
    Vertex(Eigen::Vector3d _p2, Eigen::Vector3d _p3, bool _deformed = false);
    std::shared_ptr<Vertex> deepCopy();
    //~Vertex();
    bool operator != (const Vertex &V)const{return p != V.p || p2_ori != V.p2_ori|| p3 != V.p3 || p3_ori != V.p3_ori || deformed != V.deformed;}
    bool operator == (const Vertex &V)const{return p == V.p && p2_ori == V.p2_ori && p3 == V.p3 && p3_ori == V.p3_ori && deformed == V.deformed;}
    Vertex operator-(const Vertex& V)const{return Vertex(p - V.p, p3 - V.p3);}
    Vertex operator+(const Vertex& V)const{return Vertex(p + V.p, p3 + V.p3);}
    friend Vertex operator*(double a, const Vertex& V){return Vertex(a*V.p, a*V.p3);}
};

class CrvPt_FL : public Vertex{
public:
    double s, rt;
    //std::shared_ptr<Vertex> ve;
    //std::shared_ptr<Vertex> vo;
    bool IsValid;
    CrvPt_FL(Eigen::Vector3d _p2, Eigen::Vector3d _p3,  double _s) : Vertex(_p2, _p3), s(_s), IsValid(true){}
    CrvPt_FL(Eigen::Vector3d _p2, double _s) : Vertex(_p2), s(_s), IsValid(true){}
    bool operator == (const CrvPt_FL &p)const{return s == p.s && rt == p.rt && IsValid == p.IsValid;}
    bool operator != (const CrvPt_FL &p)const{return !(s == p.s && rt == p.rt && IsValid == p.IsValid);}
    //void set(Eigen::Vector3d _p,const std::shared_ptr<Vertex>& o, const std::shared_ptr<Vertex>& e);

};

struct Vertex4d{
    bool IsCalc;
    std::shared_ptr<CrvPt_FL> first;
    std::shared_ptr<Vertex> second;
    std::shared_ptr<Vertex> third;
    Vertex4d(const std::shared_ptr<CrvPt_FL>& v, const std::shared_ptr<Vertex>& v2, const std::shared_ptr<Vertex>& v3);

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
    std::shared_ptr<Vertex> o;
    std::shared_ptr<Vertex> v;//原則の方向はv[1] - v[0] (v[0]が原点)
    EdgeType et = EdgeType::none;
    Line(const std::shared_ptr<Vertex>& _o, const std::shared_ptr<Vertex>& _v, EdgeType _et): o(_o), v(_v), et(_et), IsCrossed(-1), color(0){}
    Line():o(nullptr),v(nullptr), IsCrossed(-1), color(0) {}
    //bool operator !=(const Line& l)const{return IsCrossed != l.IsCrossed || color != l.color || (v[0] != l.v[0] && v[0] != l.v[1]);}
    //bool operator ==(const Line &l)const{return IsCrossed == l.IsCrossed && color == l.color && ((v[0] == l.v[0] && v[1] == l.v[1]) || (v[1] == l.v[0] && v[0] == l.v[1]));}
    bool is_on_line(Eigen::Vector3d p);
};

class CRV{
public:
    CRV(int _crvNum, int DivSize);

    bool drawBspline(int curveDimention,  int crvPtNum);
    bool drawBezier(int curveDimention, int crvPtNum);
    void BezierRulings(std::shared_ptr<OUTLINE>& outline, int& DivSize, int crvPtNum);
    void BsplineRulings(std::shared_ptr<OUTLINE>& outline , int& DivSize, int crvPtNum, int curveDimention);
    bool drawLine();
    void LineRulings(std::shared_ptr<OUTLINE>& outline, int DivSize);
    bool drawArc(int crvPtNum);//制御点 0,3,...: 原点. 1,4,...: 始点. 2,5,...: 終点
    void ArcRulings(std::shared_ptr<OUTLINE>& outline, int DivSize);

    bool eraseCtrlPt(int curveDimention, int crvPtNum);
    int movePtIndex(Eigen::Vector3d& p, double& dist);
    void ClearPt();
    void InsertControlPoint2(Eigen::Vector3d& p);
    void SetNewPoint();

    int InsertPointSegment;
    CurveType getCurveType();
    void setCurveType(CurveType n);
    std::vector<Eigen::Vector3d> ControllPoints;
    std::vector<Eigen::Vector3d> CurvePoints;
    std::vector<std::shared_ptr<Line>> Rulings;//偶数番目 ruling　奇数番目 グラデーションの多角形に使用
    //std::vector<meshline*> meshLines;
    Eigen::Vector3d InsertPoint;
    bool isempty;

    //グラデーション
    inline void clearColor();
private:

    bool IsInsertNewPoint;
    int OnCurvesORLines(Eigen::Vector3d& p, int& ind);//-1：どこにものっかっていない　0：曲線上　1：制御点を結んだ線上
    CurveType curveType;
    bool setPoint(const std::vector<std::shared_ptr<Vertex>>&outline, Eigen::Vector3d N, Eigen::Vector3d& cp, std::vector<Eigen::Vector3d>& P);

    double crvStep;
    int curveNum;

};

class OUTLINE{
public:

    OUTLINE();
    QString type;//Rectangle, Polygon, Polyline

    bool IsClosed();
    int VerticesNum;
    void addVertex(const std::shared_ptr<Vertex>&v, int n);
    void addVertex(Eigen::Vector3d& p);
    void eraseVertex();
    std::vector<std::shared_ptr<Vertex>> getVertices();
    std::vector<std::shared_ptr<Line>> Lines;
    void ConnectEdges(bool IsConnected = true);
    void drawPolygon(Eigen::Vector3d& p, bool IsClicked);
    void MoveOutline(Eigen::Vector3d p);//polygonの移動
    void MoveVertex(Eigen::Vector3d p, int ind);
    Eigen::Vector2d origin;//polygonの始点
    int hasPtNum; //0 ~ 2 (polygonの点の数)
    bool IsPointInFace(Eigen::Vector3d p);
    Eigen::Vector3d getNormalVec();
private:
    std::vector<std::shared_ptr<Vertex>> vertices;
    int movePointIndex(Eigen::Vector3d p);
};

void CrossDetection(std::shared_ptr<OUTLINE>& outline, std::shared_ptr<CRV>& crvs);


std::vector<double> BezierClipping(std::vector<Eigen::Vector3d>&CtrlPts, const std::shared_ptr<Vertex>& p, const std::shared_ptr<Vertex>& q, int dim);
std::vector<std::shared_ptr<Vertex>> SortPolygon(std::vector<std::shared_ptr<Vertex>>& polygon);
//std::vector<glm::f64vec3> GlobalSplineInterpolation(std::vector<CrvPt_FL>& Q, std::vector<glm::f64vec3>& CtrlPts_res, std::vector<double>& Knot, double& CurveLen, bool is3d = true, int dim = 3, int t_type = 2);


#endif // SETRULINGS_H
