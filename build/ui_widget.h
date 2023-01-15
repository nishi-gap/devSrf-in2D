/********************************************************************************
** Form generated from reading UI file 'widget.ui'
**
** Created by: Qt User Interface Compiler version 6.3.2
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
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QSlider>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include "glwidget_2d.h"
#include "glwidget_3d.h"

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
    QGroupBox *OtherBox;
    QWidget *verticalLayoutWidget_2;
    QVBoxLayout *verticalLayout_2;
    QCheckBox *SelectButton;
    QCheckBox *make_devsrf;
    QCheckBox *Reset;
    QPushButton *DebugWindow;
    QPushButton *SaveButton;
    QDoubleSpinBox *LineWidthSpinBox;
    QSlider *LineWidthSlider;
    QCheckBox *gridCBox;
    QGroupBox *LayerListWidget;
    GLWidget_3D *glWid3dim;
    QFrame *line_5;
    QGroupBox *GeometoryConstraitBox;
    QPushButton *symmetryButton;
    QGroupBox *FoldLineBox;
    QWidget *verticalLayoutWidget_5;
    QVBoxLayout *verticalLayout_5;
    QPushButton *addFL_line;
    QPushButton *addFL_arc;
    QPushButton *addFL_bezier;
    QPushButton *color_FL;

    void setupUi(QWidget *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1402, 847);
        glWid2dim = new GLWidget_2D(MainWindow);
        glWid2dim->setObjectName(QString::fromUtf8("glWid2dim"));
        glWid2dim->setGeometry(QRect(10, 10, 621, 451));
        line_4 = new QFrame(glWid2dim);
        line_4->setObjectName(QString::fromUtf8("line_4"));
        line_4->setGeometry(QRect(-10, 0, 20, 451));
        line_4->setFrameShape(QFrame::VLine);
        line_4->setFrameShadow(QFrame::Sunken);
        line_6 = new QFrame(glWid2dim);
        line_6->setObjectName(QString::fromUtf8("line_6"));
        line_6->setGeometry(QRect(0, -10, 431, 16));
        line_6->setFrameShape(QFrame::HLine);
        line_6->setFrameShadow(QFrame::Sunken);
        CurveListBox = new QGroupBox(MainWindow);
        CurveListBox->setObjectName(QString::fromUtf8("CurveListBox"));
        CurveListBox->setGeometry(QRect(520, 490, 151, 231));
        verticalLayoutWidget = new QWidget(CurveListBox);
        verticalLayoutWidget->setObjectName(QString::fromUtf8("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(10, 20, 121, 204));
        verticalLayout = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        curve_add = new QPushButton(verticalLayoutWidget);
        curve_add->setObjectName(QString::fromUtf8("curve_add"));

        verticalLayout->addWidget(curve_add);

        delete_curve = new QPushButton(verticalLayoutWidget);
        delete_curve->setObjectName(QString::fromUtf8("delete_curve"));

        verticalLayout->addWidget(delete_curve);

        CurveTypeBox = new QComboBox(verticalLayoutWidget);
        CurveTypeBox->addItem(QString());
        CurveTypeBox->addItem(QString());
        CurveTypeBox->addItem(QString());
        CurveTypeBox->setObjectName(QString::fromUtf8("CurveTypeBox"));

        verticalLayout->addWidget(CurveTypeBox);

        move_ctrl_pt = new QPushButton(verticalLayoutWidget);
        move_ctrl_pt->setObjectName(QString::fromUtf8("move_ctrl_pt"));

        verticalLayout->addWidget(move_ctrl_pt);

        insert_ctrl_pt = new QPushButton(verticalLayoutWidget);
        insert_ctrl_pt->setObjectName(QString::fromUtf8("insert_ctrl_pt"));

        verticalLayout->addWidget(insert_ctrl_pt);

        delete_ctrl_pt = new QPushButton(verticalLayoutWidget);
        delete_ctrl_pt->setObjectName(QString::fromUtf8("delete_ctrl_pt"));

        verticalLayout->addWidget(delete_ctrl_pt);

        GradationBox = new QGroupBox(MainWindow);
        GradationBox->setObjectName(QString::fromUtf8("GradationBox"));
        GradationBox->setGeometry(QRect(350, 490, 161, 201));
        AddPointsButton = new QPushButton(GradationBox);
        AddPointsButton->setObjectName(QString::fromUtf8("AddPointsButton"));
        AddPointsButton->setGeometry(QRect(0, 30, 151, 24));
        CBox_InterpolationType = new QComboBox(GradationBox);
        CBox_InterpolationType->addItem(QString());
        CBox_InterpolationType->addItem(QString());
        CBox_InterpolationType->addItem(QString());
        CBox_InterpolationType->setObjectName(QString::fromUtf8("CBox_InterpolationType"));
        CBox_InterpolationType->setGeometry(QRect(0, 60, 151, 22));
        horizontalLayoutWidget = new QWidget(GradationBox);
        horizontalLayoutWidget->setObjectName(QString::fromUtf8("horizontalLayoutWidget"));
        horizontalLayoutWidget->setGeometry(QRect(9, 100, 141, 51));
        horizontalLayout_2 = new QHBoxLayout(horizontalLayoutWidget);
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        horizontalLayout_2->setContentsMargins(0, 0, 0, 0);
        CP_colortype = new QLineEdit(horizontalLayoutWidget);
        CP_colortype->setObjectName(QString::fromUtf8("CP_colortype"));

        horizontalLayout_2->addWidget(CP_colortype);

        CP_colorval = new QLineEdit(horizontalLayoutWidget);
        CP_colorval->setObjectName(QString::fromUtf8("CP_colorval"));

        horizontalLayout_2->addWidget(CP_colorval);

        CP_colorSlider = new QSlider(GradationBox);
        CP_colorSlider->setObjectName(QString::fromUtf8("CP_colorSlider"));
        CP_colorSlider->setGeometry(QRect(0, 170, 141, 21));
        CP_colorSlider->setMinimum(-255);
        CP_colorSlider->setMaximum(255);
        CP_colorSlider->setOrientation(Qt::Horizontal);
        OutlineBox = new QGroupBox(MainWindow);
        OutlineBox->setObjectName(QString::fromUtf8("OutlineBox"));
        OutlineBox->setGeometry(QRect(20, 490, 161, 211));
        verticalLayoutWidget_3 = new QWidget(OutlineBox);
        verticalLayoutWidget_3->setObjectName(QString::fromUtf8("verticalLayoutWidget_3"));
        verticalLayoutWidget_3->setGeometry(QRect(20, 20, 121, 191));
        verticalLayout_3 = new QVBoxLayout(verticalLayoutWidget_3);
        verticalLayout_3->setObjectName(QString::fromUtf8("verticalLayout_3"));
        verticalLayout_3->setContentsMargins(0, 0, 0, 0);
        outline_rectangle = new QCheckBox(verticalLayoutWidget_3);
        outline_rectangle->setObjectName(QString::fromUtf8("outline_rectangle"));

        verticalLayout_3->addWidget(outline_rectangle);

        outline_polyline = new QCheckBox(verticalLayoutWidget_3);
        outline_polyline->setObjectName(QString::fromUtf8("outline_polyline"));

        verticalLayout_3->addWidget(outline_polyline);

        ConnectVertices = new QCheckBox(verticalLayoutWidget_3);
        ConnectVertices->setObjectName(QString::fromUtf8("ConnectVertices"));

        verticalLayout_3->addWidget(ConnectVertices);

        outline_polygon = new QCheckBox(verticalLayoutWidget_3);
        outline_polygon->setObjectName(QString::fromUtf8("outline_polygon"));

        verticalLayout_3->addWidget(outline_polygon);

        Polygon_EdgeNum = new QSpinBox(verticalLayoutWidget_3);
        Polygon_EdgeNum->setObjectName(QString::fromUtf8("Polygon_EdgeNum"));
        Polygon_EdgeNum->setMinimum(3);
        Polygon_EdgeNum->setMaximum(10);

        verticalLayout_3->addWidget(Polygon_EdgeNum);

        MoveOutLineButton = new QCheckBox(verticalLayoutWidget_3);
        MoveOutLineButton->setObjectName(QString::fromUtf8("MoveOutLineButton"));

        verticalLayout_3->addWidget(MoveOutLineButton);

        EditVertexButton = new QCheckBox(verticalLayoutWidget_3);
        EditVertexButton->setObjectName(QString::fromUtf8("EditVertexButton"));

        verticalLayout_3->addWidget(EditVertexButton);

        DivideSizeBox = new QGroupBox(MainWindow);
        DivideSizeBox->setObjectName(QString::fromUtf8("DivideSizeBox"));
        DivideSizeBox->setGeometry(QRect(210, 560, 121, 131));
        verticalLayoutWidget_4 = new QWidget(DivideSizeBox);
        verticalLayoutWidget_4->setObjectName(QString::fromUtf8("verticalLayoutWidget_4"));
        verticalLayoutWidget_4->setGeometry(QRect(0, 20, 101, 71));
        verticalLayout_4 = new QVBoxLayout(verticalLayoutWidget_4);
        verticalLayout_4->setObjectName(QString::fromUtf8("verticalLayout_4"));
        verticalLayout_4->setContentsMargins(0, 0, 0, 0);
        DivSizeSpinBox = new QSpinBox(verticalLayoutWidget_4);
        DivSizeSpinBox->setObjectName(QString::fromUtf8("DivSizeSpinBox"));
        DivSizeSpinBox->setMinimum(10);
        DivSizeSpinBox->setMaximum(1000);
        DivSizeSpinBox->setValue(30);

        verticalLayout_4->addWidget(DivSizeSpinBox);

        DvidedSizeSlider = new QSlider(verticalLayoutWidget_4);
        DvidedSizeSlider->setObjectName(QString::fromUtf8("DvidedSizeSlider"));
        DvidedSizeSlider->setMinimum(10);
        DvidedSizeSlider->setMaximum(1000);
        DvidedSizeSlider->setValue(30);
        DvidedSizeSlider->setOrientation(Qt::Horizontal);

        verticalLayout_4->addWidget(DvidedSizeSlider);

        OtherBox = new QGroupBox(MainWindow);
        OtherBox->setObjectName(QString::fromUtf8("OtherBox"));
        OtherBox->setGeometry(QRect(840, 490, 171, 241));
        verticalLayoutWidget_2 = new QWidget(OtherBox);
        verticalLayoutWidget_2->setObjectName(QString::fromUtf8("verticalLayoutWidget_2"));
        verticalLayoutWidget_2->setGeometry(QRect(10, 20, 151, 216));
        verticalLayout_2 = new QVBoxLayout(verticalLayoutWidget_2);
        verticalLayout_2->setObjectName(QString::fromUtf8("verticalLayout_2"));
        verticalLayout_2->setContentsMargins(0, 0, 0, 0);
        SelectButton = new QCheckBox(verticalLayoutWidget_2);
        SelectButton->setObjectName(QString::fromUtf8("SelectButton"));

        verticalLayout_2->addWidget(SelectButton);

        make_devsrf = new QCheckBox(verticalLayoutWidget_2);
        make_devsrf->setObjectName(QString::fromUtf8("make_devsrf"));

        verticalLayout_2->addWidget(make_devsrf);

        Reset = new QCheckBox(verticalLayoutWidget_2);
        Reset->setObjectName(QString::fromUtf8("Reset"));

        verticalLayout_2->addWidget(Reset);

        DebugWindow = new QPushButton(verticalLayoutWidget_2);
        DebugWindow->setObjectName(QString::fromUtf8("DebugWindow"));

        verticalLayout_2->addWidget(DebugWindow);

        SaveButton = new QPushButton(verticalLayoutWidget_2);
        SaveButton->setObjectName(QString::fromUtf8("SaveButton"));

        verticalLayout_2->addWidget(SaveButton);

        LineWidthSpinBox = new QDoubleSpinBox(verticalLayoutWidget_2);
        LineWidthSpinBox->setObjectName(QString::fromUtf8("LineWidthSpinBox"));
        LineWidthSpinBox->setDecimals(1);
        LineWidthSpinBox->setMinimum(0.100000000000000);
        LineWidthSpinBox->setMaximum(10.000000000000000);
        LineWidthSpinBox->setSingleStep(0.100000000000000);
        LineWidthSpinBox->setValue(1.000000000000000);

        verticalLayout_2->addWidget(LineWidthSpinBox);

        LineWidthSlider = new QSlider(verticalLayoutWidget_2);
        LineWidthSlider->setObjectName(QString::fromUtf8("LineWidthSlider"));
        LineWidthSlider->setMinimum(1);
        LineWidthSlider->setMaximum(100);
        LineWidthSlider->setSingleStep(1);
        LineWidthSlider->setPageStep(0);
        LineWidthSlider->setValue(10);
        LineWidthSlider->setSliderPosition(10);
        LineWidthSlider->setOrientation(Qt::Horizontal);

        verticalLayout_2->addWidget(LineWidthSlider);

        gridCBox = new QCheckBox(verticalLayoutWidget_2);
        gridCBox->setObjectName(QString::fromUtf8("gridCBox"));

        verticalLayout_2->addWidget(gridCBox);

        LayerListWidget = new QGroupBox(MainWindow);
        LayerListWidget->setObjectName(QString::fromUtf8("LayerListWidget"));
        LayerListWidget->setGeometry(QRect(1140, 460, 141, 311));
        glWid3dim = new GLWidget_3D(MainWindow);
        glWid3dim->setObjectName(QString::fromUtf8("glWid3dim"));
        glWid3dim->setGeometry(QRect(640, 10, 631, 441));
        line_5 = new QFrame(MainWindow);
        line_5->setObjectName(QString::fromUtf8("line_5"));
        line_5->setGeometry(QRect(-20, 400, 20, 451));
        line_5->setFrameShape(QFrame::VLine);
        line_5->setFrameShadow(QFrame::Sunken);
        GeometoryConstraitBox = new QGroupBox(MainWindow);
        GeometoryConstraitBox->setObjectName(QString::fromUtf8("GeometoryConstraitBox"));
        GeometoryConstraitBox->setGeometry(QRect(200, 490, 131, 61));
        symmetryButton = new QPushButton(GeometoryConstraitBox);
        symmetryButton->setObjectName(QString::fromUtf8("symmetryButton"));
        symmetryButton->setGeometry(QRect(10, 30, 111, 21));
        FoldLineBox = new QGroupBox(MainWindow);
        FoldLineBox->setObjectName(QString::fromUtf8("FoldLineBox"));
        FoldLineBox->setGeometry(QRect(680, 490, 131, 151));
        verticalLayoutWidget_5 = new QWidget(FoldLineBox);
        verticalLayoutWidget_5->setObjectName(QString::fromUtf8("verticalLayoutWidget_5"));
        verticalLayoutWidget_5->setGeometry(QRect(10, 30, 111, 121));
        verticalLayout_5 = new QVBoxLayout(verticalLayoutWidget_5);
        verticalLayout_5->setObjectName(QString::fromUtf8("verticalLayout_5"));
        verticalLayout_5->setContentsMargins(0, 0, 0, 0);
        addFL_line = new QPushButton(verticalLayoutWidget_5);
        addFL_line->setObjectName(QString::fromUtf8("addFL_line"));

        verticalLayout_5->addWidget(addFL_line);

        addFL_arc = new QPushButton(verticalLayoutWidget_5);
        addFL_arc->setObjectName(QString::fromUtf8("addFL_arc"));

        verticalLayout_5->addWidget(addFL_arc);

        addFL_bezier = new QPushButton(verticalLayoutWidget_5);
        addFL_bezier->setObjectName(QString::fromUtf8("addFL_bezier"));

        verticalLayout_5->addWidget(addFL_bezier);

        color_FL = new QPushButton(verticalLayoutWidget_5);
        color_FL->setObjectName(QString::fromUtf8("color_FL"));

        verticalLayout_5->addWidget(color_FL);


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

        OutlineBox->setTitle(QCoreApplication::translate("MainWindow", "Outline", nullptr));
        outline_rectangle->setText(QCoreApplication::translate("MainWindow", "Rectangle", nullptr));
        outline_polyline->setText(QCoreApplication::translate("MainWindow", "PolyLine", nullptr));
        ConnectVertices->setText(QCoreApplication::translate("MainWindow", "Connect Vertices", nullptr));
        outline_polygon->setText(QCoreApplication::translate("MainWindow", "Polygon", nullptr));
        MoveOutLineButton->setText(QCoreApplication::translate("MainWindow", "Move Outline", nullptr));
        EditVertexButton->setText(QCoreApplication::translate("MainWindow", "Edit Vertex", nullptr));
        DivideSizeBox->setTitle(QCoreApplication::translate("MainWindow", "Divde Size", nullptr));
        OtherBox->setTitle(QCoreApplication::translate("MainWindow", "Other", nullptr));
        SelectButton->setText(QCoreApplication::translate("MainWindow", "Select", nullptr));
        make_devsrf->setText(QCoreApplication::translate("MainWindow", "Apply", nullptr));
        Reset->setText(QCoreApplication::translate("MainWindow", "Reset", nullptr));
        DebugWindow->setText(QCoreApplication::translate("MainWindow", "Debug", nullptr));
        SaveButton->setText(QCoreApplication::translate("MainWindow", "Save", nullptr));
        gridCBox->setText(QCoreApplication::translate("MainWindow", "erase grid", nullptr));
        LayerListWidget->setTitle(QCoreApplication::translate("MainWindow", "Layer", nullptr));
        GeometoryConstraitBox->setTitle(QCoreApplication::translate("MainWindow", "Geometric Constraint", nullptr));
        symmetryButton->setText(QCoreApplication::translate("MainWindow", "Symmetry", nullptr));
        FoldLineBox->setTitle(QCoreApplication::translate("MainWindow", "Fold Line", nullptr));
        addFL_line->setText(QCoreApplication::translate("MainWindow", "Add Line", nullptr));
        addFL_arc->setText(QCoreApplication::translate("MainWindow", "Arc", nullptr));
        addFL_bezier->setText(QCoreApplication::translate("MainWindow", "Add Bezier", nullptr));
        color_FL->setText(QCoreApplication::translate("MainWindow", "Color", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_WIDGET_H
