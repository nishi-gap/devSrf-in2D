#ifndef GLWIDGET_2D_H
#define GLWIDGET_2D_H

#include <QOpenGLContext>
#include <QOpenGLFunctions_3_0>
#include <QOpenGLWidget>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QWheelEvent>
#include <QPointF>
#include <QMessageBox>
#include "gtoolwnd.h"
#include <QString>
#include <QPushButton>
#include "setrulings.h"
#include <QDebug>
#include <list>

class GLWidget_2D : public QOpenGLWidget, protected QOpenGLFunctions_3_0
{
    Q_OBJECT
public:
    explicit GLWidget_2D(QWidget *parent = 0);
    ~GLWidget_2D();

    //std::vector<QPointF> gradPoints;
    int ctype;//1 red, -1 blue
    int cval;
    int crvPtNum;
    std::shared_ptr<GToolWnd> gw;

    std::list<std::shared_ptr<Model>> model;
    //new gradation mode
    int InterpolationType;//0:直線, 1:spline, 2:B-spline

    bool eraseVec2d, visibleCurve;
    std::vector<std::array<Eigen::Vector3d, 2>> SingleRuling, AllRulings, NewRuling2d;
    std::vector<bool> RulingColor;
    int DivSize;

    //debug for curved folding
    bool IsStopAtFF = false, IsStopAtEq = false, IsStopAtCon = false, IsStop4Debug = false;
    double angle = 0.0;
    void EraseNonFoldEdge(bool state);
    void stashcurrentstate();
    void back2befstate();

    //regression curve
    class drawobj{
    public:
        std::vector<double>c;
        std::vector<std::vector<std::shared_ptr<Vertex>>> V;
        drawobj(const std::vector<double>&_c, const std::vector<std::vector<std::shared_ptr<Vertex>>>& _V){c = _c; V= std::move(_V);}
    };
    std::vector<drawobj> RegressionCurve; 
    void ReceiveRegressionCurve(const std::vector<std::vector<std::vector<std::shared_ptr<Vertex>>>>& RC,const std::vector<std::vector<double>>color);

    //初期状態
    void InitializeDrawMode();
    void switch2AffinMode();
    void switch2VisibleCurve();
    void CopyCurveObj();
    void PasteCurveObj();
    void switch2None();

protected:
    void initializeGL();
    void paintGL();
    void mousePressEvent(QMouseEvent *e);
    void wheelEvent(QWheelEvent *we);
    void mouseMoveEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);

public slots:

    //輪郭関係
    void DrawOutlineRectangle();
    void DrawOutlinePolygon(int state);
    void DrawOutlinePolyline(int state);
    void recieveNewEdgeNum(int num);
    void EditOutlineVertex(int state);
    void MoveOutline(int state);
    void switchGetmetricConstraint(int state);
    void ConnectVertices();

    void Reset();
    void ChangedDivSizeEdit(int n);
    void setColor();
    void receiveColors();

    //new gradation mode
    void setNewGradationMode();
    void ApplyNewGradationMode();
    void getGradationFromSlider(int val);
    void VisualizeMVColor(bool state);

    //add curve
    void AddCurve();
    void MoveCurvePt();

    void InsertNewPoint();
    void DeleteCtrlPt();
    void OpenDebugWindwow();

    //複数の曲線操作
    void cb_ApplyCurveEvent();
    void cb_DeleteCurve();
    void DeleteCurve();
    void changeSelectedCurve(int ind);
    void swapCrvsOnLayer(int n1, int n2);

    //FoldLine
    void changeFoldType(PaintTool state);
    void checkDevelopability(bool state);

    void changeflapgnle(double val,bool begin_center);
    void changeStopAtFF(bool state);
    void changeStopAtCon(bool state);
    void changeStopAtEq(bool state);
    void Start4Debug_CF();

    //Line Width
    void receiveNewLineWidth(double d);
    void switchGrid();

private:

    std::vector<Eigen::Vector3d> tmp_c;
    void draw();

    bool IsEraseNonFoldEdge;
    std::array<int,2> MoveCrvIndex; //-1: 未参照
    double camscale;//0 < camscale < 2

    int KeyEvent; //-1:None  0:Enter  1: Back-Space  2:Other
    //OutlineRectangle, RulingBezier, RulingBspline, RulingLine, RulingArc,  OutlinePolygon, OutlinePolyline, MoveControlPoint, SetColor, NewGradationMode, InsertControlPoint,
    //None(select mode), EditVertex(outline), MoveOutline, DeleteCntrlPt, DeleteCurve, OutlineConst, ConnectVertices, FoldLine, FoldlineColor

    PaintTool drawtype;
    int curveDimention;
    int maxDivSize;
    double standardDist;
    int affinmode;

    //int referencedRuling(QPointF p);
    void addPoints_intplation(QMouseEvent *e, QPointF& p);
    int assignment_refL();
    std::vector<Eigen::Vector2d> CurvePath;

    //std::vector<int> ControllPoints_gradation;//0~510 色の範囲, -1指定なし
    int DiffWheel;
    std::shared_ptr<Line> refL;
    std::shared_ptr<Vertex> refV;

    int movePt;
    CurveType curvetype;
    QList<std::pair<CurveType, PaintTool>> CurveList;

    double gridsize;
    int visibleGrid;//1:enable, -1:disable
    bool IsMVcolor_binary;

    void DrawGrid();
    bool drawpolygon;
    double rulingWidth;

    int constType;
    Eigen::Vector3d SetOnGrid(QPointF& cursol, double gridsize);
    QPointF befCur, basePoint;
    double difCursol_x, difCursol_y;
    bool IsLeftClicked, IsRightClicked, IsCopied, IsAffinMoved;

signals:
    void foldingSignals();
    void ColorChangeFrom(int wd, int color);//0:mouse, 1:LineEdit, 2:Slider
    void signalCurveType(CurveType &ct);
    void getDiviedNumber();
    int getEdgeNum();
    void SendNewActiveCheckBox(PaintTool _drawtype);
    void CurvePathSet(std::vector<Eigen::Vector2d>CurvePath);
    void deleteCrvSignal(std::vector<int> n);
    void getAlphaBeta(double& _alpha, int& _beta, int& _beta2);
    void getFoldParam(double&tol, bool& begincenter);

};

#endif // GLWIDGET_2D_H
