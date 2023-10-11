#ifndef GLWIDGET_3D_H
#define GLWIDGET_3D_H

#include "foldline.h"
#include <QWheelEvent>
#include <QOpenGLFunctions_3_0>
#include <QOpenGLWidget>
#include <QKeyEvent>
#include <QMouseEvent>
#include <QPointF>
#include <GL/glu.h>
#include <QDebug>

typedef std::vector<std::array<Eigen::Vector3d, 2>> Ruling3d;
typedef std::vector<std::shared_ptr<Vertex>> Polygon_V;
typedef std::vector<std::shared_ptr<Line>> Lines;
typedef std::vector<Eigen::Vector3d> Curve3d;
//typedef std::vector<CrvPt_FL> CrvFL3d;
typedef std::vector<std::shared_ptr<FoldLine>> FoldLine3d;

class drawobj{
public:
    std::vector<double>c;
    std::vector<std::vector<Eigen::Vector3d>> V;
    drawobj(const std::vector<double>&_c, const std::vector<std::vector<std::shared_ptr<Vertex>>>& _V){
        c = std::move(_c);
        for(auto&mesh: _V){
            std::vector<Eigen::Vector3d> tmp;
            for(auto&v: mesh)tmp.push_back(v->p3);
            V.push_back(tmp);
        }
    }
};

class GLWidget_3D : public QOpenGLWidget, protected QOpenGLFunctions_3_0
{
    Q_OBJECT
public:
    void reset();
    void setVertices(const Lines Surface = Lines(),  const Lines Rulings = Lines(),  const FoldLine3d FldCrvs = FoldLine3d(), const Ruling3d& _AllRulings = Ruling3d());
    void ReceiveParam(std::vector<std::vector<Eigen::Vector3d>>&_C);
    void ReceiveCurve(std::vector<Eigen::Vector3d>&_C, std::vector<Eigen::Vector3d>& _P);
    void ReceiveRegressionCurve(const std::vector<std::vector<std::vector<std::shared_ptr<Vertex>>>>& _RegCurve, const std::vector<std::vector<double>>&color);
    void receiveKeyEvent(QKeyEvent *e);
    void PlanarityDispay(bool state);
    void EraseNonFoldEdge(bool state);

    std::vector<std::vector<Eigen::Vector3d>>C;
    std::vector<double> PlanarityColor;
    explicit GLWidget_3D(QWidget *parent = 0);
    ~GLWidget_3D();

    std::vector<Eigen::Vector3d> FoldLineVertices;

    //regression curve
    std::vector<drawobj> RegCurve;

    Ruling3d AllRulings;
protected:
    void initializeGL();
    void paintGL(); 
    void mousePressEvent(QMouseEvent *e);
    void wheelEvent(QWheelEvent *we);
    void mouseMoveEvent(QMouseEvent *e);

public slots:


private:
    void DrawGrid();
    void draw();
    void perspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar);

    std::vector<std::vector<Eigen::Vector3d>> Vertices;
    std::vector<std::array<Eigen::Vector3d, 3>> TriMeshs;
    std::vector<Eigen::Vector3d> Curve, Points;


    Eigen::Matrix3d Mirror;
    double Scale;

    double TransX, TransY, TransZ;
    double angleX, angleY;
    bool firstRotate;
    Eigen::Vector3d center;
    void DrawMeshLines();
    void DrawMesh(bool isFront);
    int actionType;//0: None, 1: Left Click, 2: Right click, 3: other

    Eigen::Vector3d getVec(double x, double y);

    QPointF befPos;
    bool eraseMesh, eraseCtrlPt, eraseCrossPt, eraseVec, eraseCurve;
    bool VisiblePlanarity;

    double th_planarity = 1e-3;
    double drawdist;
    int drawEdgePlane;
    bool IsEraseNonFoldEdge;

    inline void dispV(Eigen::Vector3d p);
    void updateRotate();

    //ArcBallCam arccam;

};

#endif // GLWIDGET_3D_H
