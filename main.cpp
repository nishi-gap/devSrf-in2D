#include "widget.h"

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
