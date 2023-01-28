#include "widget.h"
#include "setrulings.h"
#include <glm/glm.hpp>
#include <numbers>
#include <QApplication>
//#include <QSurfaceFormat>
//https://www.bogotobogo.com/Qt/Qt5_OpenGL_QGLWidget.php

int main(int argc, char *argv[])
{

    glm::f64vec3 p{
        6.730273,
                   8.466839,
                   0.022207
                  };
    glm::f64vec3 p2{
        -20.784464,
        10.298959,
        -0.084518
    };

    glm::f64vec3 x{
        -6.521175,
        6.880964,
        0
    };
    glm::f64vec3 x2{
        -4.362441,
        -11.198432,
        0
    };
    glm::f64vec3 axis = glm::normalize(x2 - x);
    glm::f64vec3 e = glm::normalize(p - x);
    glm::f64vec3 e2 = glm::normalize(p2 - x);
    double phi3 = acos(glm::dot(e, axis));
    double phi4 = acos(glm::dot(e2, axis));
    double beta = acos(glm::dot(e,e2));
    double alpha = std::numbers::pi/2.0;
    double k = 2.0 * std::numbers::pi - phi3 - phi4;
    double phi1 = atan2((cos(k) - cos(beta)),(sin(beta)*cos(alpha)- sin(k))), phi2 = k - phi1;
    //std::cout <<"beta " << beta << std::endl;
    //std::cout << cos(k) - cos(beta) <<" , " << sin(beta)*cos(alpha)- sin(k) << std::endl;
    //std::cout << phi1 << ", " << phi2 << ", " << phi3 << ", " << phi4 << std::endl;
    //std::cout << acos(glm::dot(glm::normalize(p - x), glm::normalize(p2 - x))) << std::endl;


    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();

}
