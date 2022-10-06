#ifndef GRADATIONWIDGET_H
#define GRADATIONWIDGET_H

#include <QMouseEvent>
#include <iostream>
#include <vector>
#include <tuple>
#include <glm/glm.hpp>
#include <glm/gtx/string_cast.hpp>
#include <QPointF>
#include <QOpenGLFunctions_3_0>
#include <QOpenGLWidget>
#include <Eigen/Dense>
#include "make3d.h"

class GradationWidget : public QOpenGLWidget, protected QOpenGLFunctions_3_0
{
    Q_OBJECT
public:
    explicit GradationWidget(QWidget *parent = 0);
    ~GradationWidget();
    glm::f64vec2 a, b;
    QPointF gradPoints[2];
    void changeDrawMode(QString s);
     QString Cval, Ctype;
    std::vector<glm::f64vec2> CurvePath;
    std::vector<glm::f64vec2> ControllPoints;
    Model *model;
    void GetColorsFromNewGradationMode();

protected:
    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();
    void mousePressEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);


public slots:
    void mouseMoveEvent(QMouseEvent *e);
    void rePaint();

signals:
   void ColorValueChanged();

private:
    int IsClicked;
    bool IsSelected(QPointF p);
    float rad; //クリック時の半径
    QString paintMode;
    void DrawStdLine(int w, int h);
    void DrawCtrlPt();
    void DrawCrvPt();
    int BsplineDim;

    std::vector<double> setKnotVec(int m);;
    double basis(int j, int k, double t, std::vector<double>& T);
    void Bspline(int n, std::vector<glm::f64vec2>& P, std::vector<glm::f64vec2>& curvePt, int ptSize);
    void SplineInterpolation(std::vector<glm::f64vec2>& cp, std::vector<glm::f64vec2>& CurvePath);
};


#endif // GRADATIONWIDGET_H
