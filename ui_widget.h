/********************************************************************************
** Form generated from reading UI file 'widget.ui'
**
** Created by: Qt User Interface Compiler version 6.5.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_WIDGET_H
#define UI_WIDGET_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QFrame>
#include <QtWidgets/QGroupBox>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QScrollBar>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include <glwidget_2d.hpp>
#include <glwidget_3d.hpp>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    GLWidget_2D *glWid2dim;
    QFrame *line_4;
    QFrame *line_6;
    QGroupBox *CurveListBox;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout;
    QPushButton *curve_add;
    QPushButton *delete_curve;
    QComboBox *CurveTypeBox;
    QPushButton *move_ctrl_pt;
    QPushButton *insert_ctrl_pt;
    QPushButton *delete_ctrl_pt;
    QGroupBox *GradationBox;
    QPushButton *AddPointsButton;
    QComboBox *CBox_InterpolationType;
    QWidget *horizontalLayoutWidget;
    QHBoxLayout *horizontalLayout_2;
    QLineEdit *CP_colortype;
    QLineEdit *CP_colorval;
    QSlider *CP_colorSlider;
    QCheckBox *BinaryMVColor;
    QGroupBox *OutlineBox;
    QWidget *verticalLayoutWidget_3;
    QVBoxLayout *verticalLayout_3;
    QCheckBox *outline_rectangle;
    QCheckBox *outline_polyline;
    QCheckBox *ConnectVertices;
    QCheckBox *outline_polygon;
    QSpinBox *Polygon_EdgeNum;
    QCheckBox *MoveOutLineButton;
    QCheckBox *EditVertexButton;
    QGroupBox *DivideSizeBox;
    QWidget *verticalLayoutWidget_4;
    QVBoxLayout *verticalLayout_4;
    QSpinBox *DivSizeSpinBox;
    QSlider *DvidedSizeSlider;
    QSpinBox *ColorLimitation;
    QLabel *label;
    QGroupBox *OtherBox;
    QWidget *verticalLayoutWidget_2;
    QVBoxLayout *verticalLayout_2;
    QCheckBox *Reset;
    QCheckBox *PlanarityButton;
    QCheckBox *EraseNonFoldButton;
    QPushButton *DebugWindow;
    QPushButton *SaveButton;
    QPushButton *BackButton;
    QDoubleSpinBox *LineWidthSpinBox;
    QCheckBox *gridCBox;
    QGroupBox *LayerListWidget;
    GLWidget_3D *glWid3dim;
    QFrame *line_5;
    QGroupBox *GeometoryConstraitBox;
    QPushButton *symmetryButton;
    QGroupBox *FoldLineBox;
    QWidget *verticalLayoutWidget_5;
    QVBoxLayout *verticalLayout_5;
    QPushButton *addFL_bezier;
    QLabel *label_3;
    QSlider *ToleranceValue;
    QDoubleSpinBox *TolValue;
    QPushButton *ReassinColorButton;
    QLabel *label_2;
    QPushButton *startButton;
    QPushButton *switchDraw;
    QScrollBar *angleSlider;
    QDoubleSpinBox *angleA;
    QSlider *LineWidthSlider;
    QPushButton *OptBtn;
    QPushButton *SmoothingButton;
    QPushButton *SimpleSmoothButton;
    QDoubleSpinBox *BendWeightButton;
    QDoubleSpinBox *ParalellWeightButton;
    QLabel *label_4;
    QLabel *label_5;
    QCheckBox *ConnectEndPointBox;
    QPushButton *OptPlararity_Button;
    QCheckBox *developablityButton;

    void setupUi(QWidget *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName("MainWindow");
        MainWindow->resize(1402, 847);
        glWid2dim = new GLWidget_2D(MainWindow);
        glWid2dim->setObjectName("glWid2dim");
        glWid2dim->setGeometry(QRect(10, 10, 621, 451));
        line_4 = new QFrame(glWid2dim);
        line_4->setObjectName("line_4");
        line_4->setGeometry(QRect(-10, 0, 20, 451));
        line_4->setFrameShape(QFrame::VLine);
        line_4->setFrameShadow(QFrame::Sunken);
        line_6 = new QFrame(glWid2dim);
        line_6->setObjectName("line_6");
        line_6->setGeometry(QRect(0, -10, 431, 16));
        line_6->setFrameShape(QFrame::HLine);
        line_6->setFrameShadow(QFrame::Sunken);
        CurveListBox = new QGroupBox(MainWindow);
        CurveListBox->setObjectName("CurveListBox");
        CurveListBox->setGeometry(QRect(510, 460, 151, 231));
        verticalLayoutWidget = new QWidget(CurveListBox);
        verticalLayoutWidget->setObjectName("verticalLayoutWidget");
        verticalLayoutWidget->setGeometry(QRect(10, 20, 121, 204));
        verticalLayout = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout->setObjectName("verticalLayout");
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        curve_add = new QPushButton(verticalLayoutWidget);
        curve_add->setObjectName("curve_add");

        verticalLayout->addWidget(curve_add);

        delete_curve = new QPushButton(verticalLayoutWidget);
        delete_curve->setObjectName("delete_curve");

        verticalLayout->addWidget(delete_curve);

        CurveTypeBox = new QComboBox(verticalLayoutWidget);
        CurveTypeBox->addItem(QString());
        CurveTypeBox->addItem(QString());
        CurveTypeBox->addItem(QString());
        CurveTypeBox->setObjectName("CurveTypeBox");

        verticalLayout->addWidget(CurveTypeBox);

        move_ctrl_pt = new QPushButton(verticalLayoutWidget);
        move_ctrl_pt->setObjectName("move_ctrl_pt");

        verticalLayout->addWidget(move_ctrl_pt);

        insert_ctrl_pt = new QPushButton(verticalLayoutWidget);
        insert_ctrl_pt->setObjectName("insert_ctrl_pt");

        verticalLayout->addWidget(insert_ctrl_pt);

        delete_ctrl_pt = new QPushButton(verticalLayoutWidget);
        delete_ctrl_pt->setObjectName("delete_ctrl_pt");

        verticalLayout->addWidget(delete_ctrl_pt);

        GradationBox = new QGroupBox(MainWindow);
        GradationBox->setObjectName("GradationBox");
        GradationBox->setGeometry(QRect(340, 450, 161, 221));
        AddPointsButton = new QPushButton(GradationBox);
        AddPointsButton->setObjectName("AddPointsButton");
        AddPointsButton->setGeometry(QRect(0, 30, 151, 24));
        CBox_InterpolationType = new QComboBox(GradationBox);
        CBox_InterpolationType->addItem(QString());
        CBox_InterpolationType->addItem(QString());
        CBox_InterpolationType->addItem(QString());
        CBox_InterpolationType->setObjectName("CBox_InterpolationType");
        CBox_InterpolationType->setGeometry(QRect(0, 60, 151, 22));
        horizontalLayoutWidget = new QWidget(GradationBox);
        horizontalLayoutWidget->setObjectName("horizontalLayoutWidget");
        horizontalLayoutWidget->setGeometry(QRect(10, 90, 141, 51));
        horizontalLayout_2 = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout_2->setObjectName("horizontalLayout_2");
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        CP_colortype = new QLineEdit(horizontalLayoutWidget);
        CP_colortype->setObjectName("CP_colortype");

        horizontalLayout_2->addWidget(CP_colortype);

        CP_colorval = new QLineEdit(horizontalLayoutWidget);
        CP_colorval->setObjectName("CP_colorval");

        horizontalLayout_2->addWidget(CP_colorval);

        CP_colorSlider = new QSlider(GradationBox);
        CP_colorSlider->setObjectName("CP_colorSlider");
        CP_colorSlider->setGeometry(QRect(10, 170, 141, 21));
        CP_colorSlider->setMinimum(-255);
        CP_colorSlider->setMaximum(255);
        CP_colorSlider->setOrientation(Qt::Horizontal);
        BinaryMVColor = new QCheckBox(GradationBox);
        BinaryMVColor->setObjectName("BinaryMVColor");
        BinaryMVColor->setGeometry(QRect(20, 140, 111, 31));
        OutlineBox = new QGroupBox(MainWindow);
        OutlineBox->setObjectName("OutlineBox");
        OutlineBox->setGeometry(QRect(20, 460, 161, 211));
        verticalLayoutWidget_3 = new QWidget(OutlineBox);
        verticalLayoutWidget_3->setObjectName("verticalLayoutWidget_3");
        verticalLayoutWidget_3->setGeometry(QRect(20, 20, 121, 191));
        verticalLayout_3 = new QVBoxLayout(verticalLayoutWidget_3);
        verticalLayout_3->setObjectName("verticalLayout_3");
        verticalLayout_3->setContentsMargins(0, 0, 0, 0);
        outline_rectangle = new QCheckBox(verticalLayoutWidget_3);
        outline_rectangle->setObjectName("outline_rectangle");

        verticalLayout_3->addWidget(outline_rectangle);

        outline_polyline = new QCheckBox(verticalLayoutWidget_3);
        outline_polyline->setObjectName("outline_polyline");

        verticalLayout_3->addWidget(outline_polyline);

        ConnectVertices = new QCheckBox(verticalLayoutWidget_3);
        ConnectVertices->setObjectName("ConnectVertices");

        verticalLayout_3->addWidget(ConnectVertices);

        outline_polygon = new QCheckBox(verticalLayoutWidget_3);
        outline_polygon->setObjectName("outline_polygon");

        verticalLayout_3->addWidget(outline_polygon);

        Polygon_EdgeNum = new QSpinBox(verticalLayoutWidget_3);
        Polygon_EdgeNum->setObjectName("Polygon_EdgeNum");
        Polygon_EdgeNum->setMinimum(3);
        Polygon_EdgeNum->setMaximum(10);

        verticalLayout_3->addWidget(Polygon_EdgeNum);

        MoveOutLineButton = new QCheckBox(verticalLayoutWidget_3);
        MoveOutLineButton->setObjectName("MoveOutLineButton");

        verticalLayout_3->addWidget(MoveOutLineButton);

        EditVertexButton = new QCheckBox(verticalLayoutWidget_3);
        EditVertexButton->setObjectName("EditVertexButton");

        verticalLayout_3->addWidget(EditVertexButton);

        DivideSizeBox = new QGroupBox(MainWindow);
        DivideSizeBox->setObjectName("DivideSizeBox");
        DivideSizeBox->setGeometry(QRect(190, 530, 121, 131));
        verticalLayoutWidget_4 = new QWidget(DivideSizeBox);
        verticalLayoutWidget_4->setObjectName("verticalLayoutWidget_4");
        verticalLayoutWidget_4->setGeometry(QRect(0, 20, 101, 71));
        verticalLayout_4 = new QVBoxLayout(verticalLayoutWidget_4);
        verticalLayout_4->setObjectName("verticalLayout_4");
        verticalLayout_4->setContentsMargins(0, 0, 0, 0);
        DivSizeSpinBox = new QSpinBox(verticalLayoutWidget_4);
        DivSizeSpinBox->setObjectName("DivSizeSpinBox");
        DivSizeSpinBox->setMinimum(2);
        DivSizeSpinBox->setMaximum(300);
        DivSizeSpinBox->setValue(20);

        verticalLayout_4->addWidget(DivSizeSpinBox);

        DvidedSizeSlider = new QSlider(verticalLayoutWidget_4);
        DvidedSizeSlider->setObjectName("DvidedSizeSlider");
        DvidedSizeSlider->setMinimum(2);
        DvidedSizeSlider->setMaximum(300);
        DvidedSizeSlider->setValue(2);
        DvidedSizeSlider->setOrientation(Qt::Horizontal);

        verticalLayout_4->addWidget(DvidedSizeSlider);

        ColorLimitation = new QSpinBox(DivideSizeBox);
        ColorLimitation->setObjectName("ColorLimitation");
        ColorLimitation->setGeometry(QRect(60, 100, 42, 22));
        ColorLimitation->setMinimum(1);
        ColorLimitation->setMaximum(180);
        ColorLimitation->setValue(90);
        label = new QLabel(DivideSizeBox);
        label->setObjectName("label");
        label->setGeometry(QRect(10, 100, 51, 21));
        OtherBox = new QGroupBox(MainWindow);
        OtherBox->setObjectName("OtherBox");
        OtherBox->setGeometry(QRect(840, 450, 171, 301));
        verticalLayoutWidget_2 = new QWidget(OtherBox);
        verticalLayoutWidget_2->setObjectName("verticalLayoutWidget_2");
        verticalLayoutWidget_2->setGeometry(QRect(10, 20, 151, 267));
        verticalLayout_2 = new QVBoxLayout(verticalLayoutWidget_2);
        verticalLayout_2->setObjectName("verticalLayout_2");
        verticalLayout_2->setContentsMargins(0, 0, 0, 0);
        Reset = new QCheckBox(verticalLayoutWidget_2);
        Reset->setObjectName("Reset");

        verticalLayout_2->addWidget(Reset);

        PlanarityButton = new QCheckBox(verticalLayoutWidget_2);
        PlanarityButton->setObjectName("PlanarityButton");

        verticalLayout_2->addWidget(PlanarityButton);

        EraseNonFoldButton = new QCheckBox(verticalLayoutWidget_2);
        EraseNonFoldButton->setObjectName("EraseNonFoldButton");

        verticalLayout_2->addWidget(EraseNonFoldButton);

        DebugWindow = new QPushButton(verticalLayoutWidget_2);
        DebugWindow->setObjectName("DebugWindow");

        verticalLayout_2->addWidget(DebugWindow);

        SaveButton = new QPushButton(verticalLayoutWidget_2);
        SaveButton->setObjectName("SaveButton");

        verticalLayout_2->addWidget(SaveButton);

        BackButton = new QPushButton(verticalLayoutWidget_2);
        BackButton->setObjectName("BackButton");

        verticalLayout_2->addWidget(BackButton);

        LineWidthSpinBox = new QDoubleSpinBox(verticalLayoutWidget_2);
        LineWidthSpinBox->setObjectName("LineWidthSpinBox");
        LineWidthSpinBox->setDecimals(1);
        LineWidthSpinBox->setMinimum(0.100000000000000);
        LineWidthSpinBox->setMaximum(10.000000000000000);
        LineWidthSpinBox->setSingleStep(0.100000000000000);
        LineWidthSpinBox->setValue(1.000000000000000);

        verticalLayout_2->addWidget(LineWidthSpinBox);

        gridCBox = new QCheckBox(verticalLayoutWidget_2);
        gridCBox->setObjectName("gridCBox");

        verticalLayout_2->addWidget(gridCBox);

        LayerListWidget = new QGroupBox(MainWindow);
        LayerListWidget->setObjectName("LayerListWidget");
        LayerListWidget->setGeometry(QRect(1160, 450, 141, 311));
        glWid3dim = new GLWidget_3D(MainWindow);
        glWid3dim->setObjectName("glWid3dim");
        glWid3dim->setGeometry(QRect(640, 10, 631, 381));
        line_5 = new QFrame(MainWindow);
        line_5->setObjectName("line_5");
        line_5->setGeometry(QRect(-20, 400, 20, 451));
        line_5->setFrameShape(QFrame::VLine);
        line_5->setFrameShadow(QFrame::Sunken);
        GeometoryConstraitBox = new QGroupBox(MainWindow);
        GeometoryConstraitBox->setObjectName("GeometoryConstraitBox");
        GeometoryConstraitBox->setGeometry(QRect(190, 460, 131, 61));
        symmetryButton = new QPushButton(GeometoryConstraitBox);
        symmetryButton->setObjectName("symmetryButton");
        symmetryButton->setGeometry(QRect(10, 30, 111, 21));
        FoldLineBox = new QGroupBox(MainWindow);
        FoldLineBox->setObjectName("FoldLineBox");
        FoldLineBox->setGeometry(QRect(680, 450, 131, 181));
        verticalLayoutWidget_5 = new QWidget(FoldLineBox);
        verticalLayoutWidget_5->setObjectName("verticalLayoutWidget_5");
        verticalLayoutWidget_5->setGeometry(QRect(10, 30, 114, 111));
        verticalLayout_5 = new QVBoxLayout(verticalLayoutWidget_5);
        verticalLayout_5->setObjectName("verticalLayout_5");
        verticalLayout_5->setContentsMargins(0, 0, 0, 0);
        addFL_bezier = new QPushButton(verticalLayoutWidget_5);
        addFL_bezier->setObjectName("addFL_bezier");

        verticalLayout_5->addWidget(addFL_bezier);

        label_3 = new QLabel(verticalLayoutWidget_5);
        label_3->setObjectName("label_3");

        verticalLayout_5->addWidget(label_3);

        ToleranceValue = new QSlider(verticalLayoutWidget_5);
        ToleranceValue->setObjectName("ToleranceValue");
        ToleranceValue->setMaximum(10000);
        ToleranceValue->setSingleStep(1);
        ToleranceValue->setOrientation(Qt::Horizontal);

        verticalLayout_5->addWidget(ToleranceValue);

        TolValue = new QDoubleSpinBox(verticalLayoutWidget_5);
        TolValue->setObjectName("TolValue");
        TolValue->setDecimals(4);
        TolValue->setMaximum(10.000000000000000);
        TolValue->setSingleStep(0.010000000000000);

        verticalLayout_5->addWidget(TolValue);

        ReassinColorButton = new QPushButton(FoldLineBox);
        ReassinColorButton->setObjectName("ReassinColorButton");
        ReassinColorButton->setGeometry(QRect(20, 150, 80, 18));
        label_2 = new QLabel(MainWindow);
        label_2->setObjectName("label_2");
        label_2->setGeometry(QRect(1030, 490, 31, 16));
        startButton = new QPushButton(MainWindow);
        startButton->setObjectName("startButton");
        startButton->setGeometry(QRect(1040, 440, 80, 18));
        switchDraw = new QPushButton(MainWindow);
        switchDraw->setObjectName("switchDraw");
        switchDraw->setGeometry(QRect(1050, 630, 80, 18));
        angleSlider = new QScrollBar(MainWindow);
        angleSlider->setObjectName("angleSlider");
        angleSlider->setGeometry(QRect(1030, 520, 121, 16));
        angleSlider->setMaximum(36000);
        angleSlider->setOrientation(Qt::Horizontal);
        angleA = new QDoubleSpinBox(MainWindow);
        angleA->setObjectName("angleA");
        angleA->setGeometry(QRect(1080, 490, 62, 22));
        angleA->setDecimals(3);
        angleA->setMaximum(360.000000000000000);
        angleA->setSingleStep(0.010000000000000);
        LineWidthSlider = new QSlider(MainWindow);
        LineWidthSlider->setObjectName("LineWidthSlider");
        LineWidthSlider->setGeometry(QRect(1010, 650, 149, 15));
        LineWidthSlider->setMinimum(1);
        LineWidthSlider->setMaximum(100);
        LineWidthSlider->setSingleStep(1);
        LineWidthSlider->setPageStep(0);
        LineWidthSlider->setValue(10);
        LineWidthSlider->setSliderPosition(10);
        LineWidthSlider->setOrientation(Qt::Horizontal);
        OptBtn = new QPushButton(MainWindow);
        OptBtn->setObjectName("OptBtn");
        OptBtn->setGeometry(QRect(1040, 460, 80, 18));
        SmoothingButton = new QPushButton(MainWindow);
        SmoothingButton->setObjectName("SmoothingButton");
        SmoothingButton->setGeometry(QRect(1050, 550, 80, 18));
        SimpleSmoothButton = new QPushButton(MainWindow);
        SimpleSmoothButton->setObjectName("SimpleSmoothButton");
        SimpleSmoothButton->setGeometry(QRect(1050, 570, 91, 21));
        BendWeightButton = new QDoubleSpinBox(MainWindow);
        BendWeightButton->setObjectName("BendWeightButton");
        BendWeightButton->setGeometry(QRect(860, 400, 62, 22));
        BendWeightButton->setSingleStep(0.100000000000000);
        BendWeightButton->setValue(0.000000000000000);
        ParalellWeightButton = new QDoubleSpinBox(MainWindow);
        ParalellWeightButton->setObjectName("ParalellWeightButton");
        ParalellWeightButton->setGeometry(QRect(980, 400, 62, 22));
        ParalellWeightButton->setSingleStep(0.100000000000000);
        ParalellWeightButton->setValue(0.100000000000000);
        label_4 = new QLabel(MainWindow);
        label_4->setObjectName("label_4");
        label_4->setGeometry(QRect(820, 400, 41, 21));
        label_5 = new QLabel(MainWindow);
        label_5->setObjectName("label_5");
        label_5->setGeometry(QRect(940, 400, 41, 21));
        ConnectEndPointBox = new QCheckBox(MainWindow);
        ConnectEndPointBox->setObjectName("ConnectEndPointBox");
        ConnectEndPointBox->setGeometry(QRect(820, 420, 161, 21));
        OptPlararity_Button = new QPushButton(MainWindow);
        OptPlararity_Button->setObjectName("OptPlararity_Button");
        OptPlararity_Button->setGeometry(QRect(1050, 600, 80, 18));
        developablityButton = new QCheckBox(MainWindow);
        developablityButton->setObjectName("developablityButton");
        developablityButton->setGeometry(QRect(1070, 410, 121, 18));

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QWidget *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "MainWindow", nullptr));
        CurveListBox->setTitle(QCoreApplication::translate("MainWindow", "Curve", nullptr));
        curve_add->setText(QCoreApplication::translate("MainWindow", "Add Curve", nullptr));
        delete_curve->setText(QCoreApplication::translate("MainWindow", "Delete Curve", nullptr));
        CurveTypeBox->setItemText(0, QCoreApplication::translate("MainWindow", "Line", nullptr));
        CurveTypeBox->setItemText(1, QCoreApplication::translate("MainWindow", "B-spline", nullptr));
        CurveTypeBox->setItemText(2, QCoreApplication::translate("MainWindow", "Arc", nullptr));

        CurveTypeBox->setCurrentText(QCoreApplication::translate("MainWindow", "Line", nullptr));
        move_ctrl_pt->setText(QCoreApplication::translate("MainWindow", "Move Control Point", nullptr));
        insert_ctrl_pt->setText(QCoreApplication::translate("MainWindow", "Insert Control Point", nullptr));
        delete_ctrl_pt->setText(QCoreApplication::translate("MainWindow", "Delete Control Point", nullptr));
        GradationBox->setTitle(QCoreApplication::translate("MainWindow", "Gradation", nullptr));
        AddPointsButton->setText(QCoreApplication::translate("MainWindow", "add GradationPoints", nullptr));
        CBox_InterpolationType->setItemText(0, QCoreApplication::translate("MainWindow", "linear", nullptr));
        CBox_InterpolationType->setItemText(1, QCoreApplication::translate("MainWindow", "Spline", nullptr));
        CBox_InterpolationType->setItemText(2, QCoreApplication::translate("MainWindow", "B-spline", nullptr));

        BinaryMVColor->setText(QCoreApplication::translate("MainWindow", "binary MV color", nullptr));
        OutlineBox->setTitle(QCoreApplication::translate("MainWindow", "Outline", nullptr));
        outline_rectangle->setText(QCoreApplication::translate("MainWindow", "Rectangle", nullptr));
        outline_polyline->setText(QCoreApplication::translate("MainWindow", "PolyLine", nullptr));
        ConnectVertices->setText(QCoreApplication::translate("MainWindow", "Connect Vertices", nullptr));
        outline_polygon->setText(QCoreApplication::translate("MainWindow", "Polygon", nullptr));
        MoveOutLineButton->setText(QCoreApplication::translate("MainWindow", "Move Outline", nullptr));
        EditVertexButton->setText(QCoreApplication::translate("MainWindow", "Edit Vertex", nullptr));
        DivideSizeBox->setTitle(QCoreApplication::translate("MainWindow", "Divde Size", nullptr));
        label->setText(QCoreApplication::translate("MainWindow", "max fold", nullptr));
        OtherBox->setTitle(QCoreApplication::translate("MainWindow", "Other", nullptr));
        Reset->setText(QCoreApplication::translate("MainWindow", "Reset", nullptr));
        PlanarityButton->setText(QCoreApplication::translate("MainWindow", "planarity", nullptr));
        EraseNonFoldButton->setText(QCoreApplication::translate("MainWindow", "Erase Non Fold", nullptr));
        DebugWindow->setText(QCoreApplication::translate("MainWindow", "Debug", nullptr));
        SaveButton->setText(QCoreApplication::translate("MainWindow", "Save", nullptr));
        BackButton->setText(QCoreApplication::translate("MainWindow", "Back", nullptr));
        gridCBox->setText(QCoreApplication::translate("MainWindow", "erase grid", nullptr));
        LayerListWidget->setTitle(QCoreApplication::translate("MainWindow", "Layer", nullptr));
        GeometoryConstraitBox->setTitle(QCoreApplication::translate("MainWindow", "Geometric Constraint", nullptr));
        symmetryButton->setText(QCoreApplication::translate("MainWindow", "Symmetry", nullptr));
        FoldLineBox->setTitle(QCoreApplication::translate("MainWindow", "Fold Line", nullptr));
        addFL_bezier->setText(QCoreApplication::translate("MainWindow", "Add Bezier", nullptr));
        label_3->setText(QCoreApplication::translate("MainWindow", "tolerance", nullptr));
        ReassinColorButton->setText(QCoreApplication::translate("MainWindow", "Reassin color", nullptr));
        label_2->setText(QCoreApplication::translate("MainWindow", "alpha", nullptr));
        startButton->setText(QCoreApplication::translate("MainWindow", "start", nullptr));
        switchDraw->setText(QCoreApplication::translate("MainWindow", "switch", nullptr));
        OptBtn->setText(QCoreApplication::translate("MainWindow", "optmization", nullptr));
        SmoothingButton->setText(QCoreApplication::translate("MainWindow", "Smoothing", nullptr));
        SimpleSmoothButton->setText(QCoreApplication::translate("MainWindow", "smooth simple", nullptr));
        label_4->setText(QCoreApplication::translate("MainWindow", "bend", nullptr));
        label_5->setText(QCoreApplication::translate("MainWindow", "paralell", nullptr));
        ConnectEndPointBox->setText(QCoreApplication::translate("MainWindow", "ConnectEndPoint", nullptr));
        OptPlararity_Button->setText(QCoreApplication::translate("MainWindow", "planarity", nullptr));
        developablityButton->setText(QCoreApplication::translate("MainWindow", "developability", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_WIDGET_H
