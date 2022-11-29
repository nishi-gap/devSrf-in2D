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

protected:
    void initializeGL();
    void paintGL();
    void mousePressEvent(QMouseEvent *e);
    void wheelEvent(QWheelEvent *we);
    void mouseMoveEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);

public slots:
    //初期状態
    inline void InitializeDrawMode(int state);

    //輪郭関係
    inline void DrawOutlineRectangle();
    inline void DrawOutlinePolygon(int state);
    inline void DrawOutlinePolyline(int state);
    inline void recieveNewEdgeNum(int num);
    inline void EditOutlineVertex(int state);
    inline void MoveOutline(int state);
    inline void switchGetmetricConstraint(int state);
    inline void ConnectVertices();

    inline void Reset();
    inline void ChangedDivSizeEdit(int n);
    inline void setColor();
    inline void receiveColors();

    //new gradation mode
    inline void setNewGradationMode();
    inline void ApplyNewGradationMode();
    inline void getGradationFromSlider(int val);

    //add curve
    inline void AddCurve();
    inline void MoveCurvePt();

    inline void InsertNewPoint();
    inline void DeleteCtrlPt();
    inline void OpenDebugWindwow();

    //複数の曲線操作
    inline void cb_ApplyCurveEvent();
    inline void cb_DeleteCurve();
    void DeleteCurve();
    inline void changeSelectedCurve(int ind);
    void swapCrvsOnLayer(int n1, int n2);

    //FoldLine
    inline void changeFoldType(int state);

    //Line Width
    inline void receiveNewLineWidth(double d);
    inline void switchGrid();

private:

    std::vector<glm::f64vec3> tmp_c;
    std::vector<Vertex*> tmp_cp;

    void draw();
    int DivSize;
    int SelectedCurveIndex; //-1: 未参照
    int KeyEvent; //-1:None  0:Enter  1: Back-Space  2:Other
    //OutlineRectangle, RulingBezier, RulingBspline, RulingLine, RulingArc,  OutlinePolygon, OutlinePolyline, MoveControlPoint, SetColor, NewGradationMode, InsertControlPoint,
    //None(select mode), EditVertex(outline), MoveOutline, DeleteCntrlPt, DeleteCurve, OutlineConst, ConnectVertices, FoldLine, FoldlineColor

    PaintTool drawtype;
    int curveDimention;
    int maxDivSize;
    double standardDist;

    int referencedRuling(QPointF p);
    void addPoints_intplation(QMouseEvent *e, QPointF& p);
    HalfEdge *assignment_refHE();
    std::vector<glm::f64vec2> CurvePath;

    //std::vector<int> ControllPoints_gradation;//0~510 色の範囲, -1指定なし
    int DiffWheel;
    HalfEdge *refHE;
    int movePt;
    int curvetype;
    QList<std::pair<int, PaintTool>> CurveList;

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
    int signalCurveType();
    void getDiviedNumber();
    int getEdgeNum();
    void SendNewActiveCheckBox(PaintTool _drawtype);
    void CurvePathSet(std::vector<glm::f64vec2>CurvePath);
    void deleteCrvSignal(std::vector<int> n);
    void signalAddRulings_FL();
};

#endif // GLWIDGET_2D_H
