#ifndef MAKE3D_H
#define MAKE3D_H

#include <list>
#include <numeric>
#include <QPointF>
#include "foldline.h"


class Model{
public:
    std::vector<std::shared_ptr<Line>> Rulings;
    std::shared_ptr<OUTLINE> outline;
    std::vector<std::vector<std::shared_ptr<Vertex>>> ol_vertices;
    std::vector<std::shared_ptr<CRV>> crvs;
    std::vector<std::shared_ptr<FoldLine>> FL;
    glm::f64vec3 Axis4Const[2];
    std::shared_ptr<Vertex> Connect2Vertices[2];
    ColorPoint ColorPt;

    Model();
    Model(int _crvPtNum);
    //void deform(std::vector<std::vector<glm::f64vec3>>& output, std::vector<ruling*>& Rulings, glm::f64vec3& center);
    void deform();
    void Initialize();

    void setGradationValue(int val, const std::shared_ptr<Line>& refL, int InterpolationType, std::vector<glm::f64vec2>& CurvePath);
    void SetMaxFold(double val);
    void drawOutline(QPointF& cursol, int drawtype, double gridsize, bool IsClicked = true);
    void editOutlineVertex(QPointF& cursol, double gridsize, int event);
    void addConstraint(QPointF& cursol, int type, int gridsize, glm::f64vec3 (&axis)[2]);
    void deleteOutline(QPointF& cursol);
    void ConnectOutline(QPointF& cursol, double gridsize);
    void addRulings(); //0:move curve point, 1: add(erase, insert) curve point

    //FoldLine
    bool AddControlPoint_FL(glm::f64vec3& p, int event, int curveDimention);

    void UpdateFLOrder(int dim);
    void modify2Druling();
    void applyFL();
    void modifyFoldingCurvePositionOn3d();
    void ChangeFoldLineState();
    void applyAAAMethod(double a, bool begincenter);
    bool RevisionCrosPtsPosition();
    void SortFoldingCurve(int dim);
    bool BendingModel(double wb, double wp, int dim, bool ConstFunc = true);
    bool AssignRuling(int dim);
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

    bool SplitRulings(int dim);
    void LinearInterPolation(const std::vector<std::shared_ptr<Line>>& path);
    void SplineInterPolation(const std::vector<std::shared_ptr<Line>>& path, std::vector<glm::f64vec2>& CurvePath);

    inline void clear();

    glm::f64vec3 SetOnGrid(QPointF& cursol, double gridsize);

    std::vector<int> refCrv;//0:未参照　1:参照
    std::vector<int> refFL;
    NTree<std::shared_ptr<FoldLine>> NTree_fl;
    int crvPtNum;
    int befFaceNum;
    int FoldCurveIndex;

    std::vector<std::shared_ptr<Line>> GradationPoints;
    //std::vector<Line*> makePath();
};


#endif // MAKE3D_H
