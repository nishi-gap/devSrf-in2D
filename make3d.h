#ifndef MAKE3D_H
#define MAKE3D_H

#include <list>
#include <numeric>
#include <QPointF>
#include <QSize>
#include "foldline.h"

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
    void drawOutline(QPointF& cursol, int drawtype, double gridsize, const QSize&S, bool IsClicked = true);
    void editOutlineVertex(QPointF& cursol, double gridsize, const QSize& S, int event);
    void addConstraint(QPointF& cursol, int type, int gridsize, Eigen::Vector3d (&axis)[2], const QSize& S);
    void deleteOutline(QPointF& cursol);
    void ConnectOutline(QPointF& cursol, double gridsize, const QSize& S);
    void addRulings(); //0:move curve point, 1: add(erase, insert) curve point
    void SetOnVertices_outline();

    //FoldLine
    bool AddControlPoint_FL(Eigen::Vector3d& p, int event, int curveDimention);

    void UpdateFLOrder(int dim);
    void modify2Druling();
    void applyFL();
    void modifyFoldingCurvePositionOn3d();
    void ChangeFoldLineState();
    void applyAAAMethod(double a, double tol, bool begincenter);
    void SimplifyModel(double tol);
    void SimplifyModel(int iselim);
    bool Smoothing();
    bool RevisionCrosPtsPosition();
    void SortFoldingCurve(int dim);
    bool BendingModel(double wb, double wp, int dim, double tol, int bendrank, int alg, bool ConstFunc = true);//alg=0:ruling intersection, alg=1:regression curve
    bool AssignRuling(int dim, double tol, bool begincenter);
    int getLayerNum();
    std::vector<Eigen::Vector3d> resPts;

    //Smooth Surface

    bool AddControlPoint(Eigen::Vector3d& p, int curveDimention, int DivSize);
    int AddNewCurve(CurveType curveType, int DivSize);
    void DeleteControlPoint(QPointF pt, int curveDimention, int DivSize);

    void Check4Param(int curveDimention, std::vector<int>&deleteIndex);
    bool CrossDection4AllCurve();
    int IsSelectedCurve();
    std::array<int, 2> getSelectedCurveIndex(QPointF pt);
    std::array<int, 2> searchPointIndex(QPointF pt, int& ptInd, int type);//type = 0 -> Control Point, 1: Curve Point

    bool MoveCurvePoint(Eigen::Vector3d& p, int MoveIndex, int ptInd, int curveDimention, int DivSize);
    void SelectCurve(QPointF pt);
    int DeleteCurve();
private:

    bool SplitRulings(int dim);
    void LinearInterPolation(const std::vector<std::shared_ptr<Line>>& path);
    void SplineInterPolation(const std::vector<std::shared_ptr<Line>>& path, std::vector<Eigen::Vector2d>& CurvePath);

    inline void clear();

    Eigen::Vector3d SetOnGrid(QPointF& cursol, double gridsize, const QSize& S);

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
