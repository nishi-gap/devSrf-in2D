#ifndef MAKE3D_H
#define MAKE3D_H

#include "foldline.h"
#include "transform.h"
#include "fitting_plane.h"

class Model: public std::enable_shared_from_this<Model>{
public:
    std::vector<std::shared_ptr<Line>> Rulings;
    std::shared_ptr<OUTLINE> outline;
    std::vector<std::vector<std::shared_ptr<Vertex>>> ol_vertices;
    std::vector<std::shared_ptr<CRV>> crvs;
    std::vector<std::shared_ptr<FoldLine>> FL;
    Eigen::Vector3d Axis4Const[2];
    std::shared_ptr<Vertex> Connect2Vertices[2];
    ColorPoint ColorPt;

    std::vector<int> refCrv;//0:未参照　1:参照
    std::shared_ptr<FoldLine> refFL;
    NTree NTree_fl;
    int crvPtNum;
    int befFaceNum;
    int FoldCurveIndex;

    Model();
    Model(int _crvPtNum);
    std::shared_ptr<Model> stashcurrentstate();
    std::shared_ptr<Model> deepCopy();
    void detectClickedObj(const QPointF& curPos);//将来的にはfoldlineだけでなくほかのオブジェクトも判定して操作できるようにしたい

    void AffinTranse_Crease(int type, const QPointF& befPos, const QPointF& curPos, const QPointF& basePos);

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
    void SetOnVertices_outline(bool IsupdateEndPt);

    //FoldLine
    void RemoveUnable2GenCurve();
    bool AddControlPoint_FL(Eigen::Vector3d& p, int event, int curveDimention);
    bool AddNewFoldLine(std::shared_ptr<FoldLine>& NewFL);
    bool AssignRuling(int dim, double tol, bool begincenter);
    void applyAAAMethod(double a, bool begincenter);
    void applyFL();
    bool BendingModel(double wb, double wp, double warea, double wsim, int dim, double tol, double bndrange, int bendrank, int alg, bool IsStartEnd, bool OptimizeAngleFor3Rulings);//alg=0:ruling intersection, alg=1:regression curve
    void ChangeFoldLineState();

    int getLayerNum();
    void Interpolation(std::shared_ptr<FoldLine>& FldLine);
    void modify2Druling();
    bool Modify4LastFoldLine(std::shared_ptr<FoldLine>& tar, double warea, double wsim, double bndrange, int alg, bool IsStartEnd);
    void modifyFoldingCurvePositionOn3d();
    void FlattenSpaceCurve(std::shared_ptr<FoldLine>& FldLine, int alg);
    void flatten_lsp(std::shared_ptr<FoldLine>& FldLine);
    void movevertex(std::shared_ptr<FoldLine>& FldLine, double t);
    //void InterpolationTNB();

    void SetEndPoint(std::shared_ptr<Vertex4d>&v4d, const std::vector<std::shared_ptr<Line>>& Surface, const std::vector<std::shared_ptr<Line>>& Rulings, bool IsupdateEndPt);
    void SimplifyModel(double tol);
    void SimplifyModel(int iselim);
    bool Smoothing();
    void SortFoldingCurve(int dim);
    bool RevisionCrosPtsPosition();
    void UpdateFLOrder(int dim);
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

    void MakeTree();
    std::vector<vertexinfo> MappingVertex(bool IsRemoveOverlapping);

    void initializeSurfaceVertices();
    bool SplitRulings(int dim);
    void LinearInterPolation(const std::vector<std::shared_ptr<Line>>& path);
    //void SplineInterPolation(const std::vector<std::shared_ptr<Line>>& path, std::vector<Eigen::Vector2d>& CurvePath);

    inline void clear();

    Eigen::Vector3d SetOnGrid(QPointF& cursol, double gridsize, const QSize& S);
    std::vector<std::shared_ptr<Line>> GradationPoints;
    //std::vector<Line*> makePath();
};


#endif // MAKE3D_H
