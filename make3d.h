#ifndef MAKE3D_H
#define MAKE3D_H

#include <QMessageBox>
#include <QString>
#include <iostream>
#include <vector>
#include <tuple>
#include <glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <queue>
#include <list>
#include <numeric>
#include <QPointF>
#include "setrulings.h"

struct LinearRange{
    Face *face;//なし:nullptr
    int CrvInd;//なし:-1
    int RulInd;
    int RulingOnCurve;
    LinearRange();
};

struct FaceGradation{
    HalfEdge *he;
    double *color;
    FaceGradation();
    FaceGradation(HalfEdge *_he, double *_color);
};

class Model{
public:
    std::vector<Vertex*> vertices;
    std::vector<Face*> Faces;
    std::vector<HalfEdge*> Edges;
    OUTLINE *outline;
    //CRV *crv;
    std::vector<CRV*> crvs;

    Model();
    Model(int _crvPtNum);
    bool setRuling(std::vector<ruling*>& Rulings);
    void deform(std::vector<std::vector<glm::f64vec3>>& output, std::vector<ruling*>& Rulings, glm::f64vec3& center);

    void Initialize();
    bool devide(HalfEdge* he1, HalfEdge* he2, std::vector<Face*>& faces);

    void InsertVertex(Vertex *v);
    void ReInsertVertex(HalfEdge *he, std::vector<HalfEdge*>& edges);
    void setGradationValue(int val, int refMeshNum,  int& color, int InterpolationType, std::vector<glm::f64vec2>& CurvePath);
    void SetGradationPoint4linear(std::list<LinearRange>& path, int n = -1);

    void setOutline(std::vector<Vertex*> outline);
    void set2DMesh();
    void addRulings(); //0:move curve point, 1: add(erase, insert) curve point

    void updateRulings();

    void SelectCurve(QPointF pt);
    void AddControlPoint(QPointF pt, int curveDimention, int DivSize);
    int AddNewCurve(int curveType, int DivSize);
    void DeleteControlPoint(QPointF pt, int curveDimention, int DivSize);
    int DeleteCurve();
    void Check4Param(int curveDimention, std::vector<int>&deleteIndex);
    void MoveCurvePoint(QPointF pt, int MoveIndex, int ptInd, int curveDimention, int DivSize);
    bool CrossDection4AllCurve();
    int IsSelectedCurve();
    int getSelectedCurveIndex(QPointF pt);
    int searchPointIndex(QPointF pt, int& ptInd, int type);//type = 0 -> Control Point, 1: Curve Point

private:
    void LinearInterPolation(std::list<LinearRange>& path);
    void SplineInterPolation(std::vector<glm::f64vec2>& CurvePath);
    void setHalfEdgePair(HalfEdge *he);
    
    void LinkRulingAndGradationArea(Face *f);
    void SetGradationArea(ruling *r, meshline *ml);
    void clear();
    void deleteHE(HalfEdge *he);
    void replaceHE(HalfEdge *he);

    std::vector<int> refCrv;//0:未参照　1:参照
    std::vector<FaceGradation> Fgrad;
    LinearRange linearrange[2];
    int crvPtNum;
    int befFaceNum;
    void makePath(Face *start, Face *end, std::list<LinearRange>& path);
};


#endif // MAKE3D_H
