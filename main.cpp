#include "widget.h"
#include "setrulings.h"
#include <glm/glm.hpp>
#include <numbers>
#include <QApplication>
#include <iomanip>
#include <glm/gtx/vector_angle.hpp>
#include <mathtool.h>
//#include <QSurfaceFormat>
//https://www.bogotobogo.com/Qt/Qt5_OpenGL_QGLWidget.php

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();

}
