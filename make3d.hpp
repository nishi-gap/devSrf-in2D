#ifndef MAKE3D_H
#define MAKE3D_H

#include <list>
#include <numeric>
#include <QPointF>
#include "foldline.hpp"


class Model{
public:
    std::vector<std::shared_ptr<Line>> Rulings;
    std::shared_ptr<OUTLINE> outline;
    std::vector<std::vector<std::shared_ptr<Vertex>>> ol_vertices;
    std::vector<std::shared_ptr<CRV>> crvs;
    std::vector<std::shared_ptr<FoldLine>> FL;
    Eigen::Vector3d Axis4Const[2];
    std::shared_ptr<Vertex> Connect2Vertices[2];
    ColorPoint ColorPt;

    Model();
    Model(int _crvPtNum);
    //void deform(std::vector<std::vector<Eigen::Vector3d>>& output, std::vector<ruling*>& Rulings, Eigen::Vector3d& center);
    void deform();
    void Initialize();

    void setGradationValue(int val, const std::shared_ptr<Line>& refL, int InterpolationType, std::vector<Eigen::Vector2d>& CurvePath);
    void SetMaxFold(double val);
    void drawOutline(QPointF& cursol, int drawtype, double gridsize, bool IsClicked = true);
    void editOutlineVertex(QPointF& cursol, double gridsize, int event);
    void addConstraint(QPointF& cursol, int type, int gridsize, Eigen::Vector3d (&axis)[2]);
    void deleteOutline(QPointF& cursol);
    void ConnectOutline(QPointF& cursol, double gridsize);
    void addRulings(); //0:move curve point, 1: add(erase, insert) curve point

    //FoldLine
    bool AddControlPoint_FL(Eigen::Vector3d& p, int event, int curveDimention);

    void UpdateFLOrder(int dim);
    void modify2Druling();
    void applyFL();
    void modifyFoldingCurvePositionOn3d();
    void ChangeFoldLineState();
    void applyAAAMethod(double a, bool begincenter);
    void SimplifyModel(double tol);
    bool Smoothing();
    bool RevisionCrosPtsPosition();
    void SortFoldingCurve(int dim);
    bool BendingModel(double wb, double wp, int dim, double tol, bool ConstFunc = true);
    bool AssignRuling(int dim, double tol, bool begincenter);
    std::vector<Eigen::Vector3d> resPts;

    //Smooth Surface
    void SelectCurve(QPointF pt);
    bool AddControlPoint(Eigen::Vector3d& p, int curveDimention, int DivSize);
    int AddNewCurve(CurveType curveType, int DivSize);
    void DeleteControlPoint(QPointF pt, int curveDimention, int DivSize);
    int DeleteCurve();
    void Check4Param(int curveDimention, std::vector<int>&deleteIndex);
    bool MoveCurvePoint(Eigen::Vector3d& p, int MoveIndex, int ptInd, int curveDimention, int DivSize);
    bool CrossDection4AllCurve();
    int IsSelectedCurve();
    std::array<int, 2> getSelectedCurveIndex(QPointF pt);
    std::array<int, 2> searchPointIndex(QPointF pt, int& ptInd, int type);//type = 0 -> Control Point, 1: Curve Point

private:

    bool SplitRulings(int dim);
    void LinearInterPolation(const std::vector<std::shared_ptr<Line>>& path);
    void SplineInterPolation(const std::vector<std::shared_ptr<Line>>& path, std::vector<Eigen::Vector2d>& CurvePath);

    inline void clear();

    Eigen::Vector3d SetOnGrid(QPointF& cursol, double gridsize);

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
