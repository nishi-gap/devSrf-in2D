#include "widget.h"
#include <QApplication>
#include <iomanip>

//#include <QSurfaceFormat>
//https://www.bogotobogo.com/Qt/Qt5_OpenGL_QGLWidget.php

int main(int argc, char *argv[])
{
    int s;
    Eigen::MatrixXd _D(s,3);
    Eigen::MatrixXd N = Eigen::MatrixXd::Zero(s, s);
    Eigen::MatrixXd P = N.inverse() * _D;

    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();

}
