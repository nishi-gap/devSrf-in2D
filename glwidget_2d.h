#ifndef GLWIDGET_2D_H
#define GLWIDGET_2D_H

#include <QOpenGLContext>
#include <QOpenGLFunctions_3_0>
#include <QOpenGLWidget>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QWheelEvent>
#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <tuple>
#include <QPointF>
#include <QMessageBox>
#include <gtoolwnd.h>
#include <QString>

#include <QPushButton>
#include "setrulings.h"
#include "make3d.h"
#include "mathtool.h"

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
    GToolWnd *gw;

    Model *model;
    //new gradation mode
    int InterpolationType;//0:直線, 1:spline, 2:B-spline

    void receiveKeyEvent(QKeyEvent *e);
    bool eraseVec2d, visibleCurve;
    std::vector<std::array<glm::f64vec3, 2>> NewRuling, NewRuling2d;
    std::vector<bool> RulingColor;
    int DivSize;

    //debug for curved folding
    bool IsStopAtFF = false, IsStopAtEq = false, IsStopAtCon = false, IsStop4Debug = false;;
    double angle = 0.0;
protected:
    void initializeGL();
    void paintGL();
    void mousePressEvent(QMouseEvent *e);
    void wheelEvent(QWheelEvent *we);
    void mouseMoveEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);

public slots:
    //初期状態
    void InitializeDrawMode(int state);

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

    void changeBetaValue(int val);
    void changeStopAtFF(bool state);
    void changeStopAtCon(bool state);
    void changeStopAtEq(bool state);
    void Start4Debug_CF();

    //Line Width
    void receiveNewLineWidth(double d);
    void switchGrid();

private:

    std::vector<glm::f64vec3> tmp_c;
    void draw();

    int SelectedCurveIndex; //-1: 未参照
    int KeyEvent; //-1:None  0:Enter  1: Back-Space  2:Other
    //OutlineRectangle, RulingBezier, RulingBspline, RulingLine, RulingArc,  OutlinePolygon, OutlinePolyline, MoveControlPoint, SetColor, NewGradationMode, InsertControlPoint,
    //None(select mode), EditVertex(outline), MoveOutline, DeleteCntrlPt, DeleteCurve, OutlineConst, ConnectVertices, FoldLine, FoldlineColor

    PaintTool drawtype;
    int curveDimention;
    int maxDivSize;
    double standardDist;

    //int referencedRuling(QPointF p);
    void addPoints_intplation(QMouseEvent *e, QPointF& p);
    HalfEdge *assignment_refHE();
    std::vector<glm::f64vec2> CurvePath;

    //std::vector<int> ControllPoints_gradation;//0~510 色の範囲, -1指定なし
    int DiffWheel;
    HalfEdge *refHE;
    int movePt;
    CurveType curvetype;
    QList<std::pair<CurveType, PaintTool>> CurveList;

    double gridsize;
    int visibleGrid;//1:enable, -1:disable

    void DrawGrid();
    bool drawpolygon;
    double rulingWidth;

    int constType;
    glm::f64vec3 SetOnGrid(QPointF& cursol, double gridsize);

signals:
    void foldingSignals();
    void ColorChangeFrom(int wd, int color);//0:mouse, 1:LineEdit, 2:Slider
    void signalCurveType(CurveType &ct);
    void getDiviedNumber();
    int getEdgeNum();
    void SendNewActiveCheckBox(PaintTool _drawtype);
    void CurvePathSet(std::vector<glm::f64vec2>CurvePath);
    void deleteCrvSignal(std::vector<int> n);
    void signalAddRulings_FL();
    void getAlphaBeta(double& _alpha, int& _beta, int& _beta2);
};

#endif // GLWIDGET_2D_H
