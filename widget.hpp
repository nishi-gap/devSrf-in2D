#ifndef WIDGET_H
#define WIDGET_H

#include <QDir>
#include <QString>
#include <QWidget>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QSlider>
#include <QLineEdit>
#include <gtoolwnd.hpp>
#include <QCheckBox>
#include <QAbstractButton>
#include <QList>
#include <tuple>
#include <fstream>
#include <QFileDialog>
#include "make3d.hpp"
#include "originalbutton.hpp"


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
    void ModelBack();
    void EraseNonFoldEdge(bool state);

    //色の上限を変更
    void ChangeMaxColor(int val);

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
    void moveCtrlPts_fl();
    void addFoldLine_l();
    void addFoldLine_arc();
    void addFoldLine_bezier();
    void addFoldLine_test();
    void color_FL();
    void changeAngleFromSlider(int val);
    void changeAngleFromSpinBox(double val);

    void changeToleranceValue_Slider(int val);
    void changeToleranceValue_Spin(double val);
    void sendFoldingParam(double &tol, bool &begincenter);

    void ReassinColor();

    //Line Width
    void changeLineWidthFromSlider(int n);
    void changeLineWidthFromSpinBox(double d);

    //optimization or discrete developable surface
    void StartSmoothingSurface();
    void SimpleSmoothing();
    void StartOptimization();
    void StartOptimization_plararity();
private:
    Ui::MainWindow *ui;
    std::shared_ptr<Model> model;
    QList<std::tuple<QCheckBox *, PaintTool >> CBoxlist;
    int crvPtNum;
    std::vector<std::vector<Eigen::Vector3d>> output;
    void exportobj();

    std::vector<std::shared_ptr<Btn4Crv>> LayerList;

    std::shared_ptr<Btn4Crv> SelectedBtn;
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
    void signalFLtype(PaintTool state);
    void signalNewLineWidth(double d);

    void sendAngle(double val, double tol, int keyType);

};
#endif // WIDGET_H
