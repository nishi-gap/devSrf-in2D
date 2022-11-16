#ifndef GLWIDGET_3D_H
#define GLWIDGET_3D_H

#include <QWheelEvent>
#include <QOpenGLFunctions_3_0>
#include <QOpenGLWidget>
#include <QMouseEvent>
#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <QPointF>
#include<GL/glu.h>
#include <QKeyEvent>
#include <glm/gtx/string_cast.hpp>
#include "setrulings.h"

class GLWidget_3D : public QOpenGLWidget, protected QOpenGLFunctions_3_0
{
    Q_OBJECT
public:
    void setVertices(std::vector<std::vector<glm::f64vec3>>& _Vertices, glm::f64vec3 center);
    void setVertices(std::vector<Face*>& Faces);
    void receive(std::vector<std::vector<glm::f64vec3>>& l, std::vector<std::vector<glm::f64vec3>>& r, glm::f64vec3 center);
    std::vector<int>left;
    explicit GLWidget_3D(QWidget *parent = 0);
    ~GLWidget_3D();

protected:
    void initializeGL();
    void paintGL(); 
    void mousePressEvent(QMouseEvent *e);
    void wheelEvent(QWheelEvent *we);
    void mouseMoveEvent(QMouseEvent *e);

public slots:
    void RotationX(int _RotX);
    void RotationY(int _RotY);
    void ChangeTranslateX(int _X);
    void ChangeTranslateY(int _Y);
    void ChangeTranslateZ(int _Z);

private:
    void DrawGrid();
    void draw();
    void perspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar);

    std::vector<std::vector<glm::f64vec3>> Vertices;
    std::vector<std::array<glm::f64vec3, 3>> TriMeshs;
    double TransX, TransY, TransZ;
    double RotX, RotY;
    bool firstRotate;
    glm::f64vec3 center;
    void DrawMeshLines();
    void DrawMesh(bool isFront);
    int actionType;//0: None, 1: Left Click, 2: Right click, 3: other

    QPointF befPos;


};

#endif // GLWIDGET_3D_H
