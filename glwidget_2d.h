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
#include "make3d.h"
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

    int crvPtNum;
    std::list<std::shared_ptr<Model>> model;
    bool visibleCurve;
    int DivSize;
    void EraseNonFoldEdge(bool state);

    class drawobj{
    public:
        std::vector<double>c;
        std::vector<std::vector<std::shared_ptr<Vertex>>> V;
        drawobj(const std::vector<double>&_c, const std::vector<std::vector<std::shared_ptr<Vertex>>>& _V){c = _c; V= std::move(_V);}
    };

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

    void Reset();
    void ChangedDivSizeEdit(int n);

    //gradation
    void setGradationMode();
    void DrawGradationMode();
    void GetGradationFromSlider(int val);
    void VisualizeMVColor(bool state);

    //ruling control curve
    void AddCurve();
    void MoveCurvePt();
    void InsertNewPoint();
    void DeleteCurve();
    void changeSelectedCurve(int ind);
    void swapCrvsOnLayer(int n1, int n2);

    void AddNewCrease(); //crease
    void changeflapgnle(double val);

    //Line Width
    void receiveNewLineWidth(double d);
    void switchGrid();

private:

    void draw();
    bool IsEraseNonFoldEdge;
    std::array<int,2> MoveCrvIndex; //-1: 未参照
    double camscale;//0 < camscale < 2
    int KeyEvent; //-1:None  0:Enter  1: Back-Space  2:Other

    PaintTool drawtype;
    int curveDimention;
    int maxDivSize;
    double standardDist;
    int affinmode;

    void AddPoints4Gradation(QMouseEvent *e, QPointF& p);
    std::vector<Eigen::Vector2d> CurvePath;
    CurveType curvetype;
    int DiffWheel;
    std::shared_ptr<Line> refLine;

    int movePt;
    double gridsize;
    int visibleGrid;//1:enable, -1:disable
    bool IsMVcolor_binary;

    void DrawGrid();
    double rulingWidth;
    QPointF befCur, basePoint;
    double difCursol_x, difCursol_y;
    bool IsLeftClicked, IsRightClicked, IsCopied, IsAffinMoved;

signals:
    void foldingSignals();
    void ColorChangeFrom(int wd, int color);//0:mouse, 1:Slider
    void signalCurveType(CurveType &ct);
    void getDiviedNumber();
    int getEdgeNum();
    void SendNewActiveCheckBox(PaintTool _drawtype);
    void deleteCrvSignal(std::vector<int> n);
};

#endif // GLWIDGET_2D_H
