#ifndef WIDGET_H
#define WIDGET_H

#include <QString>
#include <QWidget>
#include <QSpinBox>
#include <QPushButton>
#include <QSlider>
#include <QLineEdit>
#include <QOpenGLWidget>
#include <gtoolwnd.h>
#include <QCheckBox>
#include <QAbstractButton>
#include <QList>
#include <tuple>
#include <fstream>
#include <QFileDialog>
#include "make3d.h"
#include "originalbutton.h"
#include "mathtool.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE

class MainWindow : public QWidget
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
    QMenuBar *menubar;

protected:
    void keyPressEvent(QKeyEvent *e);
    void mouseMoveEvent(QMouseEvent *e);

public slots:
    void fold_Sm();

    void ChangedDivSizeEdit();
    void ChangeDivSizeEditFromSlider(int val);
    void ChangeDivSizeEditFromSpinBox(int val);
    void sendCurveType(CurveType &ct);
    void Initialize();

    //new gradation mode
    void changeInterpolationType();
    void ApplyNewColor(int wd, int color);

    //輪郭
    void sendNewEdgeNum();
    void switchActivateCheckBox(PaintTool active);
    void SymmetricConstraint();

    //layer関係
    void addCurveBtn();
    void RemoveBtnFromLayerCrv(std::vector<int> n);
    void SetHandleCrv(Btn4Crv *btn, QMouseEvent *e);

    //fold line
    void addFoldLine_l();
    void addFoldLine_bsp();
    void addFoldLine_bezier();
    void color_FL();
    void fold_FL();

    //Line Width
    void changeLineWidthFromSlider(int n);
    void changeLineWidthFromSpinBox(double d);
private:
    Ui::MainWindow *ui;
    Model *model;
    QList<std::tuple<QCheckBox *, PaintTool >> CBoxlist;
    int crvPtNum;
    std::vector<std::vector<glm::f64vec3>> output;
    void exportobj();

    std::vector<Btn4Crv*> LayerList;

    Btn4Crv *SelectedBtn;
    QPoint dragPos;
    QRect originalPos;
    int CurvesNum[4];

signals:
    void sliderValChanged();
    void makeGradation();
    void signalNewEdgeNum(int num);
    void PressedEnter();
    void PressedBackSpace();
    void signalNewSelectedCrv(int ind);
    void swapIndex(int n1, int n2);
    void constraintType(int state);
    void signalFLtype(int state);
    void signalNewLineWidth(double d);

};
#endif // WIDGET_H
