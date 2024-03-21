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

typedef std::vector<std::shared_ptr<Line>> Lines;
typedef std::vector<std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>>> FoldLine3d;

class GLWidget_3D : public QOpenGLWidget, protected QOpenGLFunctions_3_0
{
    Q_OBJECT
public:
    void setVertices(const Lines Surface = Lines(),  const Lines Rulings = Lines(),  const FoldLine3d Creases = FoldLine3d());
    void receiveKeyEvent(QKeyEvent *e);
    void EraseNonFoldEdge(bool state);
    explicit GLWidget_3D(QWidget *parent = 0);
    ~GLWidget_3D();
    std::vector<std::vector<Eigen::Vector3d>> FoldLineVertices;
protected:
    void initializeGL();
    void paintGL(); 
    void mousePressEvent(QMouseEvent *e);
    void wheelEvent(QWheelEvent *we);
    void mouseMoveEvent(QMouseEvent *e);

public slots:

private:
    void perspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar);

    std::vector<std::vector<Eigen::Vector3d>> Vertices;
    double Scale;
    double TransX, TransY, TransZ;
    double angleX, angleY;
    bool firstRotate;
    void DrawMeshLines();
    void DrawMesh(bool isFront);
    int actionType;//0: None, 1: Left Click, 2: Right click, 3: other

    Eigen::Vector3d getVec(double x, double y);

    QPointF befPos;
    bool eraseCtrlPt, eraseCrossPt, eraseVec, eraseCurve;
    double drawdist;
    int drawEdgePlane;
    bool IsEraseNonFoldEdge;

    inline void dispV(Eigen::Vector3d p);
    void updateRotate();

};

#endif // GLWIDGET_3D_H
