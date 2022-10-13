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
#include "ui_gtoolwnd.h"
#include "setrulings.h"
#include "make3d.h"
#include "originalbutton.h"

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
    void InitializeDrawMode(int state);

    //輪郭関係
    void DrawOutlineRectangle();
    void DrawOutlinePolygon(int state);
    void DrawOutlinePolyline(int state);
    void recieveNewEdgeNum(int num);
    void EditOutlineVertex(int state);
    void MoveOutline(int state);

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

private:
    void draw();
    int DivSize;
    int SelectedCurveIndex; //-1: 未参照
    int KeyEvent; //-1:None  0:Enter  1: Back-Space  2:Other
    //OutlineRectangle, RulingBezier, RulingBspline, RulingLine,  OutlinePolygon, OutlinePolyline, MoveControlPoint, SetColor, NewGradationMode, InsertControlPoint,
    //None(select mode), EditVertex(outline), MoveOutline, DeleteCntrlPt, DeleteCurve
    QString drawtype;
    int curveDimention;
    int maxDivSize;
    double standardDist;

    int referencedRuling(QPointF p);
    void addPoints_intplation(QMouseEvent *e, QPointF& p);
    void assignment_refHE(QPointF& p);
    std::vector<glm::f64vec2> CurvePath;

    //std::vector<int> ControllPoints_gradation;//0~510 色の範囲, -1指定なし
    int DiffWheel;
    HalfEdge *refHE;
    int movePt;
    int curvetype;
    QList<std::tuple<QString, int, QString >> CurveList;

signals:
    void foldingSignals();
    void ColorChangeFrom(int wd, int color);//0:mouse, 1:LineEdit, 2:Slider
    QString signalCurveType();
    void getDiviedNumber();
    int getEdgeNum();
    void SendNewActiveCheckBox(QString activeDrawType);
    void CurvePathSet(std::vector<glm::f64vec2>CurvePath);
    void deleteCrvSignal(std::vector<int> n);
};

#endif // GLWIDGET_2D_H
