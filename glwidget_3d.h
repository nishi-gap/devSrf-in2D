#ifndef GLWIDGET_3D_H
#define GLWIDGET_3D_H

#include <QWheelEvent>
#include <QOpenGLFunctions_3_0>
#include <QOpenGLWidget>
#include <QKeyEvent>
#include <QMouseEvent>
#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <QPointF>
#include<GL/glu.h>
#include <QKeyEvent>
#include <glm/gtx/string_cast.hpp>
//#include "arcball.h"
#include "setrulings.h"
#include "foldline.h"

typedef std::vector<std::array<glm::f64vec3, 2>> Ruling3d;
typedef std::vector<Vertex*> Polygon_V;
typedef std::vector<Line*> Lines;
typedef std::vector<glm::f64vec3> Curve3d;
//typedef std::vector<CrvPt_FL> CrvFL3d;
typedef std::vector<Vertex4d> FoldLine3d;
typedef std::vector<FoldLine3d> FoldLine3ds;

class GLWidget_3D : public QOpenGLWidget, protected QOpenGLFunctions_3_0
{
    Q_OBJECT
public:
    void setVertices(const Lines Surface = Lines(),  const Lines Rulings = Lines(),  const FoldLine3ds FldCrvs = FoldLine3ds(), const Ruling3d& _AllRulings = Ruling3d());
    void ReceiveParam(std::vector<std::vector<glm::f64vec3>>&_C);
    void receiveKeyEvent(QKeyEvent *e);
    void PlanarityDispay(bool state);
    void EraseNonFoldEdge(bool state);

    std::vector<std::vector<glm::f64vec3>>C;
    std::vector<double> PlanarityColor;
    explicit GLWidget_3D(QWidget *parent = 0);
    ~GLWidget_3D();

    std::vector<glm::f64vec3> FoldLineVertices;
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

    std::vector<std::vector<glm::f64vec3>> Vertices;
    std::vector<std::array<glm::f64vec3, 3>> TriMeshs;

    glm::f64mat4x4 Mirror;
    glm::f64mat4x4 Scale;

    double TransX, TransY, TransZ;
    double angleX, angleY;
    bool firstRotate;
    glm::f64vec3 center;
    void DrawMeshLines();
    void DrawMesh(bool isFront);
    int actionType;//0: None, 1: Left Click, 2: Right click, 3: other

    glm::f64vec3 getVec(float x, float y);

    QPointF befPos;
    bool eraseMesh, eraseCtrlPt, eraseCrossPt, eraseVec, eraseCurve;
    bool VisiblePlanarity;

    double th_planarity = 1e-3;
    int switchTNB;
    double drawdist;
    std::vector<std::array<glm::f64vec3, 2>> drawEdges;
    int drawEdgePlane;
    bool IsEraseNonFoldEdge;

    inline void dispV(glm::f64vec3 p);
    void updateRotate();

    //ArcBallCam arccam;

};

#endif // GLWIDGET_3D_H
