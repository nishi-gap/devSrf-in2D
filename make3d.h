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

class Model{
public:
    std::vector<Vertex*> vertices;
    //std::vector<Face*> Faces;
    //std::vector<HalfEdge*> Edges;
    std::vector<Line*> Rulings;
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

    void setGradationValue(int val, Line *refL, int InterpolationType, std::vector<glm::f64vec2>& CurvePath);
    void SetMaxFold(double val);
    void drawOutline(QPointF& cursol, int drawtype, double gridsize, bool IsClicked = true);
    void editOutlineVertex(QPointF& cursol, double gridsize, int event);
    void addConstraint(QPointF& cursol, int type, int gridsize, glm::f64vec3 (&axis)[2]);
    void deleteOutline(QPointF& cursol);
    void ConnectOutline(QPointF& cursol, double gridsize);
    void addRulings(); //0:move curve point, 1: add(erase, insert) curve point

    //FoldLine
    bool AddControlPoint_FL(glm::f64vec3& p, int event, int curveDimention, int FoldCurveIndex);
    bool SplitRulings(FoldLine *NewFL, int dim);
    bool updateSplitRulings(FoldLine *NewFL, int dim);
    void modify2Druling();
    void applyFL();
    void modifyFoldingCurvePositionOn3d();
    void ChangeFoldLineState();
    std::vector<glm::f64vec3> resPts;

    //Smooth Surface
    void SelectCurve(QPointF pt);
    bool AddControlPoint(glm::f64vec3& p, int curveDimention, int DivSize);
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
    void LinearInterPolation(std::vector<Line*>& path);
    void SplineInterPolation(std::vector<Line*>& path, std::vector<glm::f64vec2>& CurvePath);

    inline void clear();

    glm::f64vec3 SetOnGrid(QPointF& cursol, double gridsize);

    std::vector<int> refCrv;//0:未参照　1:参照
    std::vector<int> refFL;

    int crvPtNum;
    int befFaceNum;
    int FoldCurveIndex;

    std::vector<Line*> GradationPoints;
    //std::vector<Line*> makePath();
};


#endif // MAKE3D_H
