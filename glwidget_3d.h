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
#include "setrulings.h"

typedef std::vector<std::array<glm::f64vec3, 2>> Ruling3d;
typedef std::vector<Face*> Faces3d;
typedef std::vector<glm::f64vec3> Curve3d;
typedef std::vector<CrvPt_FL> CrvFL3d;

enum class ArcBallMode{
    none,
    translate,
    rotate,
    scale,
};

class GLWidget_3D : public QOpenGLWidget, protected QOpenGLFunctions_3_0
{
    Q_OBJECT
public:
    void setVertices(const Faces3d& Faces = Faces3d(), const Curve3d& _CtrlPts = Curve3d(),
                     const Curve3d& _Curve = Curve3d(), const CrvFL3d& _CrossPts = CrvFL3d(),
                     const Ruling3d& _Vl = Ruling3d(),const Ruling3d& _Vr = Ruling3d());
    void receive(std::vector<std::vector<glm::f64vec3>>& l, std::vector<std::vector<glm::f64vec3>>& r, glm::f64vec3 center);
    void receiveKeyEvent(QKeyEvent *e);

    std::vector<int>left;
    explicit GLWidget_3D(QWidget *parent = 0);
    ~GLWidget_3D();

    std::vector<glm::f64vec3> CtrlPts, Curve;
    std::vector<CrvPt_FL> CrossPts;
    std::vector<std::array<glm::f64vec3, 2>> Vl, Vr;
protected:
    void initializeGL();
    void paintGL(); 
    void mousePressEvent(QMouseEvent *e);
    void wheelEvent(QWheelEvent *we);
    void mouseMoveEvent(QMouseEvent *e);
    void keyPressEvent(QKeyEvent *e);

public slots:


private:
    void DrawGrid();
    void draw();
    void perspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar);

    std::vector<std::vector<glm::f64vec3>> Vertices;
    std::vector<std::array<glm::f64vec3, 3>> TriMeshs;
    double TransX, TransY, TransZ;
    double angleX, angleY;
    bool firstRotate;
    glm::f64vec3 center;
    void DrawMeshLines();
    void DrawMesh(bool isFront);
    int actionType;//0: None, 1: Left Click, 2: Right click, 3: other

    ArcBallMode arcballMode = ArcBallMode::none;
    glm::f64vec3 getVec(float x, float y);

    QPointF befPos;
    bool eraseMesh, eraseCtrlPt, eraseCrossPt, eraseVec, eraseCurve;
    int switchTNB;
    double drawdist;

    inline void dispV(glm::f64vec3 p);
    void updateRotate();
};

#endif // GLWIDGET_3D_H
