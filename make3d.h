#ifndef MAKE3D_H
#define MAKE3D_H

#include "foldline.h"
#include "transform.h"
#include "fitting_plane.h"

class Model: public std::enable_shared_from_this<Model>{
public:
    std::vector<std::shared_ptr<Line>> Rulings;
    std::shared_ptr<OUTLINE> outline;
    //std::vector<std::vector<std::shared_ptr<Vertex>>> ol_vertices;
    std::vector<std::shared_ptr<CRV>> RulingCurve;
    std::vector<std::shared_ptr<FoldLine>> Creases;
    ColorPoint ColorPt;

    std::vector<int> refCrv;//0:未参照　1:参照
    std::shared_ptr<FoldLine> refCreases;
    NTree NTree_Creases;
    int crvPtNum;
    int FoldCurveIndex;

    Model();
    Model(int _crvPtNum);
    std::shared_ptr<Model> stashcurrentstate();
    std::shared_ptr<Model> deepCopy();
    void detectClickedObj(const QPointF& curPos);//将来的にはfoldlineだけでなくほかのオブジェクトも判定して操作できるようにしたい

    void AffinTranse_Crease(int type, const QPointF& befPos, const QPointF& curPos, const QPointF& basePos);

    void deform();
    void Initialize();

    void SetGradationValue(int val, const std::shared_ptr<Line>& refL);
    void SetMaxFold(double val);
    void drawOutline(QPointF& cursol, int drawtype, double gridsize, const QSize&S, bool IsClicked = true);
    void editOutlineVertex(QPointF& cursol, double gridsize, const QSize& S, int event);

    void deleteOutline(QPointF& cursol);
    void addRulings(); //0:move curve point, 1: add(erase, insert) curve point
    void SetOnVertices_outline(bool IsupdateEndPt);

    //FoldLine
    void RemoveUnable2GenCurve();
    bool AddNewFoldLine(std::shared_ptr<FoldLine>& NewFL);
    bool AssignRuling(int dim);
    void applyFL();
    bool BendingModel(double wp, double wsim, int dim, int alg);//alg=0:ruling intersection, alg=1:regression curve
    void ChangeFoldLineState();

    bool Modify4LastFoldLine(std::shared_ptr<FoldLine>& tar, double wp, double wsim);
    void modifyFoldingCurvePositionOn3d();
    void flatten_lsp(std::shared_ptr<FoldLine>& FldLine);

    void SetEndPoint(std::shared_ptr<Vertex4d>&v4d, const std::vector<std::shared_ptr<Line>>& Surface, const std::vector<std::shared_ptr<Line>>& Rulings, bool IsupdateEndPt);
    void SimplifyModel(int iselim);
    bool Smoothing();
    void UpdateTree(int dim);

    //event:0 -> add ,event:1 -> delete
    bool ControlPoint(PaintTool dtype, int event, QPointF p, int curveDimention, int DivSize);
    int AddNewCurve(CurveType curveType, int DivSize);

    void Check4Param(int curveDimention, std::vector<int>&deleteIndex);
    bool CrossDection4AllCurve();
    int IsSelectedCurve();
    std::array<int, 2> getSelectedCurveIndex(QPointF pt);
    std::array<int, 2> searchPointIndex(QPointF _p, int& ptInd, int type);//type = 0 -> Control Point, 1: Curve Point

    bool MoveCurvePoint(QPointF _p, int MoveIndex, int ptInd, int curveDimention, int DivSize);
    void SelectCurve(QPointF pt);
    int DeleteCurve();
private:
    int movingVertexIndex;
    void MakeTree();
    std::vector<vertexinfo> MappingVertex(bool IsRemoveOverlapping);

    void initializeSurfaceVertices();
    bool SplitRulings(int dim);
    void LinearInterPolation(const std::vector<std::shared_ptr<Line>>& path);

    inline void clear();


    std::vector<std::shared_ptr<Line>> GradationPoints;
};


#endif // MAKE3D_H
