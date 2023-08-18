#ifndef GRADATIONWIDGET_H
#define GRADATIONWIDGET_H

#include <QMouseEvent>
#include <iostream>
#include <vector>
#include <tuple>
#include <QPointF>
#include <QOpenGLFunctions_3_0>
#include <QOpenGLWidget>
#include "make3d.h"

class GradationWidget : public QOpenGLWidget, protected QOpenGLFunctions_3_0
{
    Q_OBJECT
public:
    explicit GradationWidget(QWidget *parent = 0);
    ~GradationWidget();
    Eigen::Vector2d a, b;
    QPointF gradPoints[2];
    void changeDrawMode(QString s);
     QString Cval, Ctype;
    std::vector<Eigen::Vector2d> CurvePath;
    std::vector<Eigen::Vector2d> ControllPoints;
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
    void Bspline(int n, std::vector<Eigen::Vector2d>& P, std::vector<Eigen::Vector2d>& curvePt, int ptSize);
    void SplineInterpolation(std::vector<Eigen::Vector2d>& cp, std::vector<Eigen::Vector2d>& CurvePath);
};


#endif // GRADATIONWIDGET_H
