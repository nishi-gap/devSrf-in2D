#ifndef WIDGET_H
#define WIDGET_H

#include <QDir>
#include <QString>
#include <QWidget>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QPushButton>
#include <QSlider>
#include <QMenuBar>
#include <QLineEdit>
#include <QCheckBox>
#include <QAbstractButton>
#include <QList>
#include <tuple>
#include <fstream>
#include <QFileDialog>
#include <QFileInfo>
#include <utility>
#include "originalbutton.h"


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

    //new gradation mode
    void ApplyNewColor(int wd, int color);
    void FinishColorChange();

    //輪郭
    void sendNewEdgeNum();
    void switchActivateCheckBox(PaintTool active);

    //layer関係
    void addCurveBtn();
    void RemoveBtnFromLayerCrv(std::vector<int> n);
    void SetHandleCrv(Btn4Crv *btn, QMouseEvent *e);

    //fold line
    void changeAngleFromSlider(int val);
    void changeAngleFromSpinBox(double val);
    void changeRulingNum(int val);

    void changeLineWidthFromSlider(int n);//Line Width
    void StartOptimization();//optimization or discrete developable surface

private:
    Ui::MainWindow *ui;
    QList<std::tuple<QCheckBox *, PaintTool >> CBoxlist;
    int crvPtNum;
    std::vector<std::vector<Eigen::Vector3d>> output;
    void exportobj();
    void exportsvg(QString filename);

    std::vector<std::shared_ptr<Btn4Crv>> LayerList;

    std::shared_ptr<Btn4Crv> SelectedBtn;
    QPoint dragPos;
    QRect originalPos;
    int CurvesNum[4];

signals:
    void sliderValChanged();
    void makeGradation();
    void signalNewEdgeNum(int num);
    void signalNewSelectedCrv(int ind);
    void swapIndex(int n1, int n2);
    void signalNewLineWidth(double d);
    void sendAngle(double val);
};
#endif // WIDGET_H
