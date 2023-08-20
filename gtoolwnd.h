#ifndef GTOOLWND_H
#define GTOOLWND_H

#include <QWidget>
#include <QMainWindow>
#include <QPushButton>
#include <QComboBox>
#include <vector>
#include <tuple>
#include "ui_gtoolwnd.h"
#include "setrulings.hpp"
namespace Ui {
class GToolWnd;
}

class GToolWnd : public QMainWindow
{
    Q_OBJECT

public:
    explicit GToolWnd(QWidget *parent = nullptr);
    ~GToolWnd();
    Ui::GToolWnd *gtw;

public slots:
    void set(std::vector<Eigen::Vector2d>CurvePath);

signals:

    void updateCurvePath();
};

#endif // GTOOLWND_H
