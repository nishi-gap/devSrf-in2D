#include "widget.h"
#include <QApplication>
#include <iomanip>

#include <iostream>
#include <vector>
#include <gsl/gsl_spline.h>
//#include <QSurfaceFormat>
//https://www.bogotobogo.com/Qt/Qt5_OpenGL_QGLWidget.php


int main(int argc, char *argv[])
{

    // 3次元空間上の点のデータ
    std::vector<double> x_data = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> y_data = {0.0, 1.0, 0.0, -1.0, 0.0};
    std::vector<double> z_data = {1.0, 2.0, 3.0, 2.0, 1.0};
    size_t n = x_data.size(); // データポイントの数

    // x座標、y座標、z座標それぞれに対するスプライン補間のための準備
    gsl_interp_accel *acc_x = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_y = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_z = gsl_interp_accel_alloc();

    gsl_spline *spline_x = gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline *spline_y = gsl_spline_alloc(gsl_interp_cspline, n);
    gsl_spline *spline_z = gsl_spline_alloc(gsl_interp_cspline, n);

    // x座標、y座標、z座標それぞれに対するスプライン補間の計算
    gsl_spline_init(spline_x, x_data.data(), y_data.data(), n);
    gsl_spline_init(spline_y, x_data.data(), y_data.data(), n);
    gsl_spline_init(spline_z, x_data.data(), z_data.data(), n);

    // 3次元空間上の点を評価
    double x_eval = 2.5;
    double y_eval = gsl_spline_eval(spline_x, x_eval, acc_x);
    double z_eval = gsl_spline_eval(spline_z, x_eval, acc_z);

    // 結果の表示
    std::cout << "At x = " << x_eval << ", y = " << y_eval << ", z = " << z_eval << std::endl;

    // メモリの解放
    gsl_spline_free(spline_x);
    gsl_spline_free(spline_y);
    gsl_spline_free(spline_z);
    gsl_interp_accel_free(acc_x);
    gsl_interp_accel_free(acc_y);
    gsl_interp_accel_free(acc_z);
    return 0;


    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();

}
