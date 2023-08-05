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
#include "foldline.h"
#include "mathtool.h"


struct FaceGradation{
    HalfEdge *he;
    double *color;
    FaceGradation(): color(0), he(nullptr){}
    FaceGradation(HalfEdge *_he, double *_color): he(_he), color(_color){}
};



class Model{
public:
    std::vector<Vertex*> vertices;
    std::vector<Face*> Faces;
    std::vector<HalfEdge*> Edges;
    OUTLINE *outline;
    std::vector<std::vector<Vertex*>> ol_vertices;
    std::vector<CRV*> crvs;
    std::vector<FoldLine*> FL;
    glm::f64vec3 Axis4Const[2];
    Vertex* Connect2Vertices[2];
    ColorPoint ColorPt;

    Model();
    Model(int _crvPtNum);
    //void deform(std::vector<std::vector<glm::f64vec3>>& output, std::vector<ruling*>& Rulings, glm::f64vec3& center);
    void deform();
    void Initialize();
    bool devide(HalfEdge* he1, HalfEdge* he2, std::vector<Face*>& faces);

    HalfEdge* InsertVertex(Vertex *v);
    void setGradationValue(int val, HalfEdge *refHE,int InterpolationType, std::vector<glm::f64vec2>& CurvePath);
    void SetMaxFold(double val);
    void setOutline();
    void drawOutline(QPointF& cursol, int drawtype, double gridsize, bool IsClicked = true);
    void editOutlineVertex(QPointF& cursol, double gridsize, int event);
    void addConstraint(QPointF& cursol, int type, int gridsize, glm::f64vec3 (&axis)[2]);
    void deleteOutline(QPointF& cursol);
    void ConnectOutline(QPointF& cursol, double gridsize);
    void addRulings(); //0:move curve point, 1: add(erase, insert) curve point

    //FoldLine
    bool AddControlPoint_FL(glm::f64vec3& p, int event, int curveDimention, int FoldCurveIndex);
    void modify2Druling();
    void applyFL();

    std::vector<glm::f64vec3> resPts;

    //Smooth Surface
    void SelectCurve(QPointF pt);
    void AddControlPoint(glm::f64vec3& p, int curveDimention, int DivSize);
    int AddNewCurve(CurveType curveType, int DivSize);
    void DeleteControlPoint(QPointF pt, int curveDimention, int DivSize);
    int DeleteCurve();
    void Check4Param(int curveDimention, std::vector<int>&deleteIndex);
    void MoveCurvePoint(glm::f64vec3& p, int MoveIndex, int ptInd, int curveDimention, int DivSize);
    bool CrossDection4AllCurve();
    int IsSelectedCurve();
    int getSelectedCurveIndex(QPointF pt);
    int searchPointIndex(QPointF pt, int& ptInd, int type);//type = 0 -> Control Point, 1: Curve Point

private:
    void LinearInterPolation(std::vector<HalfEdge*>& path);
    void SplineInterPolation(std::vector<HalfEdge*>& path, std::vector<glm::f64vec2>& CurvePath);
    void setHalfEdgePair(HalfEdge *he);

    inline void clear();
    void ConnectEdge(HalfEdge *he);

    glm::f64vec3 SetOnGrid(QPointF& cursol, double gridsize);

    std::vector<int> refCrv;//0:未参照　1:参照
    std::vector<int> refFL;
    std::vector<FaceGradation> Fgrad;

    int crvPtNum;
    int befFaceNum;

    std::vector<HalfEdge*> GradationPoints;
    std::vector<HalfEdge*> makePath();
};


#endif // MAKE3D_H
