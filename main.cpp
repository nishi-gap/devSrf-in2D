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


#include <iostream>
#include <vector>
#include <nlopt.hpp>

typedef struct {
    double a, b;
} my_constraint_data;

double h = 1e-7;
double myfunc(const std::vector<double> &x, std::vector<double> &grad, void* f_data)
{
    if (!grad.empty()) {
        grad[0] = 0.0;
        grad[1] = (sqrt(x[1] + h) - sqrt(x[1] - h))/(2.0*h);
    }
    //std::cout <<"func " << sqrt(x[1]) << std::endl;
    return sqrt(x[1]);
}

double myconstraint(const std::vector<double> &x, std::vector<double> &grad, void* f_data)
{
    my_constraint_data *d = (my_constraint_data *) f_data;
    double a = d->a, b = d->b;
    if (!grad.empty()) {
        double g0p = std::pow((a*(x[0] + h) + b), 3) - x[1], g0m = std::pow((a*(x[0] - h) + b), 3) - x[1];
        grad[0] = (g0p - g0m)/(2.0*h);
        double g1p = std::pow((a*x[0] + b), 3) - (x[1] + h), g1m = std::pow((a*x[0] + b), 3) - (x[1] - h);
        grad[1] = (g1p - g1m)/(2.0*h);
    }
    //std::cout << "const " << std::pow((a*x[0] + b) * (a*x[0] + b), 3) - x[1] << std::endl;
    return std::pow((a*x[0] + b), 3) - x[1];
 }

int main(int argc, char *argv[])
{

    nlopt::opt opt(nlopt::LD_MMA, 2);
    std::vector<double> lb(2);
    lb[0] = -HUGE_VAL; lb[1] = 0;
    opt.set_lower_bounds(lb);
    opt.set_min_objective(myfunc, NULL);
    my_constraint_data data[2] = { {2,0}, {-1,1} };
    opt.add_inequality_constraint(myconstraint, &data[0]);
    opt.add_inequality_constraint(myconstraint, &data[1]);
    opt.set_xtol_rel(1e-4);
    std::vector<double> x(2);
    x[0] = 1.234; x[1] = 5.678;
    double minf;

    try{
        nlopt::result result = opt.optimize(x, minf);
        std::cout <<"result = " << result << std::endl;
        std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
            << std::setprecision(10) << minf << std::endl;
    }
    catch(std::exception &e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    //return 0;
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();

}
