#include "widget.h"
#include <random>
#include "setrulings.h"
#include <glm/glm.hpp>
#include <fstream>
#include <QApplication>
//#include <QSurfaceFormat>
//https://www.bogotobogo.com/Qt/Qt5_OpenGL_QGLWidget.php

int main(int argc, char *argv[])
{

    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();

}
