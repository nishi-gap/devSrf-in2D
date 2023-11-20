#include "widget.h"
#include <QDebug>
#include "ui_widget.h"
MainWindow::MainWindow(QWidget *parent)
    : QWidget(parent), ui(new Ui::MainWindow)
{

    crvPtNum = 3000;
    output.clear();

    ui->setupUi(this);
    ui->glWid2dim->DivSize = ui->DivSizeSpinBox->value();
    //model = new Model(crvPtNum);
    //ui->glWid2dim->model = model;
    CBoxlist = {{ui->outline_rectangle, PaintTool::Rectangle_ol}, {ui->outline_polygon,  PaintTool::Polygon_ol},{ui->outline_polyline, PaintTool::Polyline_ol},
                {ui->EditVertexButton, PaintTool::EditVertex_ol}, {ui->MoveOutLineButton, PaintTool::Move_ol},{ui->Reset, PaintTool::Reset}};

    setGeometry(0,0,1250, 620);
    //色の最大値
    connect(ui->ColorLimitation, &QSpinBox::valueChanged, this, &MainWindow::ChangeMaxColor);
    connect(ui->BinaryMVColor, &QCheckBox::clicked, ui->glWid2dim, &GLWidget_2D::VisualizeMVColor);//山谷の色を二値化するかグラデーションを描画するか切り替える
    connect(ui->Reset, &QCheckBox::stateChanged,ui->glWid2dim,&GLWidget_2D::Reset);
    connect(ui->Reset, &QCheckBox::stateChanged, this, &MainWindow::Initialize);
    connect(ui->DvidedSizeSlider, &QSlider::valueChanged, this, &MainWindow::ChangeDivSizeEditFromSlider);
    connect(ui->DivSizeSpinBox, &QSpinBox::valueChanged, this, &MainWindow::ChangeDivSizeEditFromSpinBox);
    connect(ui->glWid2dim, &GLWidget_2D::foldingSignals, this , &MainWindow::fold_Sm);

    //line width
    connect(ui->LineWidthSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::changeLineWidthFromSpinBox);
    connect(ui->LineWidthSlider, &QSlider::valueChanged, this, &MainWindow::changeLineWidthFromSlider);
    connect(this, &MainWindow::signalNewLineWidth, ui->glWid2dim, &GLWidget_2D::receiveNewLineWidth);

    //visualize Grid
    connect(ui->gridCBox, &QCheckBox::clicked, ui->glWid2dim, &GLWidget_2D::switchGrid);

    connect(ui->AddPointsButton,&QPushButton::clicked,ui->glWid2dim,&GLWidget_2D::setNewGradationMode);
    connect(ui->CBox_InterpolationType,&QComboBox::currentIndexChanged, this, &MainWindow::changeInterpolationType);
    connect(ui->glWid2dim, &GLWidget_2D::ColorChangeFrom, this, &MainWindow::ApplyNewColor);
    connect(this, &MainWindow::makeGradation, ui->glWid2dim, &GLWidget_2D::ApplyNewGradationMode);
    connect(ui->CP_colorSlider, &QSlider::valueChanged, ui->glWid2dim, &GLWidget_2D::getGradationFromSlider);
    connect(ui->btn_ColorFinish, &QPushButton::clicked, this, &MainWindow::FinishColorChange);

    //曲線
    connect(ui->curve_add, &QPushButton::clicked, ui->glWid2dim, &GLWidget_2D::AddCurve);
    connect(ui->curve_add, &QPushButton::clicked, this, &MainWindow::addCurveBtn);
    connect(ui->move_ctrl_pt, &QPushButton::clicked, ui->glWid2dim, &GLWidget_2D::MoveCurvePt);
    //connect(ui->CurveTypeBox, &QComboBox::currentTextChanged, ui->glWid2dim, &GLWidget_2D::ChangeCurveType);
    connect(ui->glWid2dim, &GLWidget_2D::signalCurveType, this, &MainWindow::sendCurveType);
    connect(ui->insert_ctrl_pt, &QPushButton::clicked, ui->glWid2dim, &GLWidget_2D::InsertNewPoint);
    connect(ui->delete_ctrl_pt, &QPushButton::clicked, ui->glWid2dim, &GLWidget_2D::DeleteCtrlPt);
    connect(this, &MainWindow::PressedEnter, ui->glWid2dim, &GLWidget_2D::cb_ApplyCurveEvent);
    connect(this, &MainWindow::PressedBackSpace, ui->glWid2dim, &GLWidget_2D::cb_DeleteCurve);
    connect(ui->delete_curve, &QPushButton::clicked, ui->glWid2dim, &GLWidget_2D::DeleteCurve);

    //ruling
    connect(ui->glWid2dim, &GLWidget_2D::getDiviedNumber, this, &MainWindow::ChangedDivSizeEdit);

    //輪郭関係
    connect(ui->outline_rectangle, &QCheckBox::clicked, ui->glWid2dim, &GLWidget_2D::DrawOutlineRectangle);
    connect(ui->outline_polygon, &QCheckBox::clicked, ui->glWid2dim, &GLWidget_2D::DrawOutlinePolygon);
    connect(ui->outline_polyline, &QCheckBox::clicked, ui->glWid2dim, &GLWidget_2D::DrawOutlinePolyline);
    connect(ui->EditVertexButton, &QCheckBox::clicked, ui->glWid2dim, &GLWidget_2D::EditOutlineVertex);
    connect(ui->MoveOutLineButton, &QCheckBox::clicked, ui->glWid2dim, &GLWidget_2D::MoveOutline);
    connect(ui->Polygon_EdgeNum, &QSpinBox::valueChanged, this, &MainWindow::sendNewEdgeNum);
    connect(this, &MainWindow::signalNewEdgeNum, ui->glWid2dim, &GLWidget_2D::recieveNewEdgeNum);
    connect(ui->glWid2dim, &GLWidget_2D::getEdgeNum, this, &MainWindow::sendNewEdgeNum);
    connect(ui->glWid2dim, &GLWidget_2D::SendNewActiveCheckBox, this, &MainWindow::switchActivateCheckBox);
    connect(ui->symmetryButton, &QPushButton::clicked, this, &MainWindow::SymmetricConstraint);
    connect(this, &MainWindow::constraintType, ui->glWid2dim, &GLWidget_2D::switchGetmetricConstraint);
    connect(ui->ConnectVertices, &QPushButton::clicked, ui->glWid2dim, &GLWidget_2D::ConnectVertices);

    //FoldLine
    connect(ui->addFL_bezier, &QPushButton::clicked, this, &MainWindow::addFoldLine_bezier);
    connect(this, &MainWindow::signalFLtype, ui->glWid2dim, &GLWidget_2D::changeFoldType);
    connect(ui->ToleranceValue, &QSlider::valueChanged, this, &MainWindow::changeToleranceValue_Slider);
    connect(ui->elimRuling, &QSpinBox::valueChanged, this, &MainWindow::changeRulingNum);
    connect(ui->spinBox_bendcurve, &QSpinBox::valueChanged, this, &MainWindow::BendCurve);

    connect(ui->TolValue, &QDoubleSpinBox::valueChanged, this, &MainWindow::changeToleranceValue_Spin);
    connect(ui->ReassinColorButton, &QPushButton::clicked, this, &MainWindow::ReassinColor);
    connect(ui->glWid2dim, &GLWidget_2D::getFoldParam, this, &MainWindow::sendFoldingParam);

    //fold line debug
    connect(ui->startButton, &QPushButton::clicked, ui->glWid2dim, &GLWidget_2D::Start4Debug_CF);
    connect(ui->developablityButton, &QCheckBox::clicked, ui->glWid2dim, &GLWidget_2D::checkDevelopability);

    connect(ui->angleSlider, &QSlider::sliderMoved, this, &MainWindow::changeAngleFromSlider);
    connect(ui->angleA, &QDoubleSpinBox::valueChanged, this, &MainWindow::changeAngleFromSpinBox);
    connect(this, &MainWindow::sendAngle, ui->glWid2dim, &GLWidget_2D::changeflapgnle);

    connect(ui->DebugWindow, &QPushButton::clicked,ui->glWid2dim,&GLWidget_2D::OpenDebugWindwow);
    connect(ui->SaveButton, &QPushButton::clicked,this, &MainWindow::exportobj);

    LayerList.clear();
    connect(ui->glWid2dim, &GLWidget_2D::deleteCrvSignal, this, &MainWindow::RemoveBtnFromLayerCrv);
    connect(this, &MainWindow::signalNewSelectedCrv, ui->glWid2dim, &GLWidget_2D::changeSelectedCurve);
    connect(this, &MainWindow::swapIndex, ui->glWid2dim, &GLWidget_2D::swapCrvsOnLayer);

    //planarity
    connect(ui->SimpleSmoothButton, &QPushButton::clicked, this, &MainWindow::SimpleSmoothing);
    connect(ui->PlanarityButton, &QCheckBox::clicked, ui->glWid3dim, &GLWidget_3D::PlanarityDispay);

    //optimization or discrete developable surface
    connect(ui->OptBtn, &QPushButton::clicked, this, &MainWindow::StartOptimization);

    //Erase Non Fold Edge
    connect(ui->EraseNonFoldButton, &QCheckBox::clicked, this, &MainWindow::EraseNonFoldEdge);

    //regression curve
    connect(ui->VisualizeRegCrv, &QCheckBox::clicked, this, &MainWindow::SwitchingVisualization_RegCurve);

    CurvesNum[0] = MainWindow::CurvesNum[1] = MainWindow::CurvesNum[2] = MainWindow::CurvesNum[3] = 0;
    SelectedBtn = nullptr;

}

//regression curve
static bool IsvisibleRegCrv = false;
static bool IsStartEnd = false;
static int befNum = 0;
static bool OptimizeAngleFor3Rulings = false;

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::SwitchingVisualization_RegCurve(){
    IsvisibleRegCrv  = !IsvisibleRegCrv;
    qDebug() << "IsvisibleRegCrv = " << IsvisibleRegCrv;
    std::vector<std::vector<std::shared_ptr<Vertex>>> empty;
    if(!IsvisibleRegCrv){
        std::vector<double> color_reg{0.8, 0, 0.}, color_fixreg{0,0,0.8};
        ui->glWid3dim->ReceiveRegressionCurve({}, {});
        ui->glWid2dim->ReceiveRegressionCurve({}, {});
    }
}

void MainWindow::SymmetricConstraint(){ emit constraintType(0);}

void MainWindow::Initialize(){
    LayerList.clear();
    ui->spinBox_bendcurve->setValue(0);
    ui->glWid3dim->reset();
    update();
}

void MainWindow::FinishColorChange(){
    qDebug()<< "gradation mode finish";
    ui->glWid2dim->InitializeDrawMode();
}
void MainWindow::ChangeMaxColor(int val){ui->glWid2dim->model.back()->SetMaxFold((double)val);}

void MainWindow::addFoldLine_l(){emit signalFLtype(PaintTool::FoldLine_test);}
void MainWindow::addFoldLine_arc(){emit signalFLtype(PaintTool::FoldLine_arc);}
void MainWindow::addFoldLine_bezier(){emit signalFLtype(PaintTool::FoldLine_bezier);}
void MainWindow::addFoldLine_test(){emit signalFLtype(PaintTool::FoldLine_test);}
void MainWindow::moveCtrlPts_fl(){emit signalFLtype(PaintTool::FoldLine_move);}
void MainWindow::color_FL(){emit signalFLtype(PaintTool::FoldLine_move);}

void MainWindow::ModelBack(){
    ui->glWid2dim->update();
    ui->glWid3dim->setVertices(ui->glWid2dim->model.back()->outline->Lines, ui->glWid2dim->model.back()->Rulings, ui->glWid2dim->model.back()->FL, ui->glWid2dim->AllRulings);
}

void MainWindow::changeRulingNum(int n){
    int iselim = (n < befNum)? -1: (n == befNum)? 0: 1;
    if(iselim == -1 || iselim == 1)ui->glWid2dim->model.back()->SimplifyModel(iselim);
    befNum = n;
    ui->glWid2dim->update();
    ui->glWid3dim->setVertices(ui->glWid2dim->model.back()->outline->Lines, ui->glWid2dim->model.back()->Rulings, ui->glWid2dim->model.back()->FL, ui->glWid2dim->AllRulings);
}

void MainWindow::EraseNonFoldEdge(bool state){
    ui->glWid2dim->EraseNonFoldEdge(state);
    ui->glWid3dim->EraseNonFoldEdge(state);
    if(ui->glWid2dim->model.back()->FL.empty() || ui->glWid2dim->model.back()->FL[0]->FoldingCurve.empty())return;
    ui->glWid3dim->setVertices(ui->glWid2dim->model.back()->outline->Lines, ui->glWid2dim->model.back()->Rulings, ui->glWid2dim->model.back()->FL, ui->glWid2dim->AllRulings);
}

void MainWindow::sendFoldingParam(double &tol, bool &begincenter){
    tol = ui->TolValue->value();
    begincenter = false;
}

void MainWindow::changeToleranceValue_Slider(int val){
    double maxSpin = ui->TolValue->maximum(), maxSlider = ui->ToleranceValue->maximum();
    if(ui->glWid2dim->model.back()->FL.empty() || ui->glWid2dim->model.back()->FL[0]->FoldingCurve.empty())return;
    if(!ui->BinaryMVColor->isChecked())ui->BinaryMVColor->setChecked(true);
    ui->glWid2dim->setNewGradationMode();
    ui->glWid2dim->VisualizeMVColor(true);
    double tol = maxSpin * (double)val/(double)maxSlider;
    ui->TolValue->setValue(tol);
    ui->glWid2dim->model.back()->SimplifyModel(tol);
    auto fl = ui->glWid2dim->model.back()->FL[0];
    ui->glWid2dim->update();
     ui->glWid3dim->setVertices(ui->glWid2dim->model.back()->outline->Lines, ui->glWid2dim->model.back()->Rulings, ui->glWid2dim->model.back()->FL, ui->glWid2dim->AllRulings);
}

void MainWindow::changeToleranceValue_Spin(double val){
    double maxSpin = ui->TolValue->maximum(), maxSlider = ui->ToleranceValue->maximum();
    if(ui->glWid2dim->model.back()->FL.empty() || ui->glWid2dim->model.back()->FL[0]->FoldingCurve.empty())return;
    if(!ui->BinaryMVColor->isChecked())ui->BinaryMVColor->setChecked(true);
    ui->glWid2dim->setNewGradationMode();
    ui->glWid2dim->VisualizeMVColor(true);

    ui->ToleranceValue->setValue(int(val/maxSpin * maxSlider));
    //ui->glWid2dim->model.back()->SimplifyModel(val);
    fold_Sm();
}



void MainWindow::StartOptimization(){
    if(ui->glWid2dim->model.back()->FL.empty() || ui->glWid2dim->model.back()->FL[0]->FoldingCurve.empty())return;
    double tol = ui->TolValue->value();
    double wb = ui->BendWeightButton->value(), wp = ui->ParalellWeightButton->value();
    double warea = ui->TriAreaWeight->value(), wsim = ui->NormErrorWeight->value();
    double bndrange = ui->BoundaryRange->value();
    int layerNum = ui->glWid2dim->model.back()->getLayerNum();
    ui->glWid2dim->model.back()->BendingModel(wb, wp, warea, wsim, 3, tol, bndrange, layerNum, 1, IsStartEnd, OptimizeAngleFor3Rulings);
    fold_Sm();

}

void MainWindow::BendCurve(int num){
    if(ui->glWid2dim->model.back()->FL.empty() || ui->glWid2dim->model.back()->FL[0]->FoldingCurve.empty())return;
    double tol = ui->TolValue->value();
    double wb = ui->BendWeightButton->value(), wp = ui->ParalellWeightButton->value();
    double warea = ui->TriAreaWeight->value(), wsim = ui->NormErrorWeight->value();
    double bndrange = ui->BoundaryRange->value();
    qDebug() << "the number of bend curve is " << num;
    ui->glWid2dim->model.back()->BendingModel(wb, wp,warea, wsim, 3, tol, bndrange, num, 1, IsStartEnd ,OptimizeAngleFor3Rulings);
    fold_Sm();
}


void MainWindow::SimpleSmoothing(){
    auto Poly_V = ui->glWid2dim->model.back()->outline->getVertices();
    if(ui->glWid2dim->model.back()->FL.empty() || ui->glWid2dim->model.back()->FL[0]->FoldingCurve.empty())return;
    bool res = ui->glWid2dim->model.back()->Smoothing();
    if(res){
        ui->glWid2dim->update();
        ui->glWid2dim->model.back()->SetOnVertices_outline(false);
         ui->glWid3dim->setVertices(ui->glWid2dim->model.back()->outline->Lines, ui->glWid2dim->model.back()->Rulings, ui->glWid2dim->model.back()->FL, ui->glWid2dim->AllRulings);
    }
}

void MainWindow::changeAngleFromSlider(int val){
    ui->angleA->setValue((double)val/100);
    double a = (double)val/18000.0 * std::numbers::pi;
    emit sendAngle(a, IsStartEnd);

    if(IsvisibleRegCrv){
         std::vector<double> color_reg{0.8, 0, 0.}, color_fixreg{0,0,0.8};
         std::vector<std::shared_ptr<Vertex>> Poly_V = ui->glWid2dim->model.back()->outline->getVertices();
         std::vector<std::vector<std::shared_ptr<Vertex>>> Tri_fixside;
         auto Triangles = ui->glWid2dim->model.back()->FL.back()->CalclateRegressionCurve(static_cast<double>(val)/100.0, Poly_V, false, IsStartEnd, Tri_fixside);
         ui->glWid3dim->ReceiveRegressionCurve({Triangles, {}}, {color_reg,color_fixreg});
         ui->glWid2dim->ReceiveRegressionCurve({Triangles, {}}, {color_reg,color_fixreg});
    }
}

void MainWindow::changeAngleFromSpinBox(double val){
    ui->angleSlider->setValue(val*100);
    double a = (double)val*std::numbers::pi/180.0;
    emit sendAngle(a, IsStartEnd);
    if(IsvisibleRegCrv){
        std::vector<std::shared_ptr<Vertex>> Poly_V = ui->glWid2dim->model.back()->outline->getVertices();
        std::vector<std::vector<std::shared_ptr<Vertex>>> Tri_fixside;
        std::vector<double> color_reg{0.8, 0, 0.}, color_fixreg{0,0,0.8};
        std::vector<std::vector<std::shared_ptr<Vertex>>> Triangles = ui->glWid2dim->model.back()->FL.back()->CalclateRegressionCurve(val, Poly_V, false, IsStartEnd, Tri_fixside);
        ui->glWid3dim->ReceiveRegressionCurve({Triangles}, {color_reg});
        ui->glWid2dim->ReceiveRegressionCurve({Triangles}, {color_reg});
    }
}

void MainWindow::changeLineWidthFromSlider(int n){
    double d = (double)n/10.;
    ui->LineWidthSpinBox->setValue(d);
    emit signalNewLineWidth(d);
}
void MainWindow::changeLineWidthFromSpinBox(double d){
    ui->LineWidthSlider->setValue((int)(d * 10));
    emit signalNewLineWidth(d);
}

void MainWindow::fold_Sm(){
    ui->glWid2dim->update();
    if(!ui->glWid2dim->model.back()->outline->IsClosed())ui->glWid3dim->setVertices();

    else{
        ui->glWid2dim->model.back()->SetOnVertices_outline(false);
        auto Surface = ui->glWid2dim->model.back()->outline->Lines;
        auto Rulings = ui->glWid2dim->model.back()->Rulings;
        auto FL = ui->glWid2dim->model.back()->FL;
        auto AllRulings = ui->glWid2dim->AllRulings;
        ui->glWid3dim->setVertices(Surface, Rulings , FL , AllRulings);
    }
}

void MainWindow::ChangedDivSizeEdit(){
    int val = ui->DivSizeSpinBox->value();
    this->ui->glWid2dim->ChangedDivSizeEdit(val);
}

void MainWindow::ChangeDivSizeEditFromSlider(int val){
    ui->DivSizeSpinBox->setValue(val);
    this->ui->glWid2dim->ChangedDivSizeEdit(val);
}

void MainWindow::ChangeDivSizeEditFromSpinBox(int val){
    ui->DvidedSizeSlider->setValue(val);
    this->ui->glWid2dim->ChangedDivSizeEdit(val);
}

void MainWindow::changeInterpolationType(){
    ui->glWid2dim->InterpolationType = this->ui->CBox_InterpolationType->currentIndex();
}

void MainWindow::ApplyNewColor(int wd, int color){
    QString ctype = (color > 0) ? "Red": "Blue" , cval = QString::number(abs(color));
    if(wd == 0){//mouse
        this->ui->CP_colortype->setText(ctype);
        this->ui->CP_colorval->setText(cval);
        this->ui->CP_colorSlider->setValue(color);
    }else if(wd == 1){//line edit
        this->ui->CP_colorSlider->setValue(color);

    }else{//slider
        this->ui->CP_colortype->setText(ctype);
        this->ui->CP_colorval->setText(cval);
    }
    emit makeGradation();
}

void MainWindow::sendNewEdgeNum(){
    emit signalNewEdgeNum(ui->Polygon_EdgeNum->value());
}

void MainWindow::switchActivateCheckBox(PaintTool active){
    for(auto& T: CBoxlist){
        if(active != std::get<1>(T))std::get<0>(T)->setCheckState(Qt::Unchecked);
    }

}
void MainWindow::ReassinColor(){
    if(!ui->glWid2dim->model.back()->outline->IsClosed())ui->glWid3dim->setVertices();
    else if(!ui->glWid2dim->model.back()->FL.empty() && !ui->glWid2dim->model.back()->FL.back()->FoldingCurve.empty()){
        ui->glWid2dim->model.back()->FL.back()->ReassignColor();
        ui->glWid2dim->model.back()->deform();
        ui->glWid2dim->update();
        ui->glWid2dim->model.back()->modifyFoldingCurvePositionOn3d();
         ui->glWid3dim->setVertices(ui->glWid2dim->model.back()->outline->Lines, ui->glWid2dim->model.back()->Rulings, ui->glWid2dim->model.back()->FL, ui->glWid2dim->AllRulings);
    }
    else{
         ui->glWid3dim->setVertices(ui->glWid2dim->model.back()->outline->Lines, ui->glWid2dim->model.back()->Rulings, ui->glWid2dim->model.back()->FL, ui->glWid2dim->AllRulings);
    }
}



void MainWindow::keyPressEvent(QKeyEvent *e){
    ui->glWid3dim->receiveKeyEvent(e);
    double bndrange = ui->BoundaryRange->value();
    double warea = ui->TriAreaWeight->value(), wsim = ui->NormErrorWeight->value();
    double wb = ui->BendWeightButton->value(), wp = ui->ParalellWeightButton->value();
    int layerNum = ui->glWid2dim->model.back()->getLayerNum();
    int dim = 3;

    if(e->key() == Qt::Key_0){
         qDebug()<<"optimization for last adding folding curve";
         ui->glWid2dim->model.back()->FL.back()->Optimization_FlapAngle(ui->glWid2dim->model.back()->outline->vertices,  wb,  wp,  0,  1,  IsStartEnd, 0, OptimizeAngleFor3Rulings);
    }
    if(e->key() == Qt::Key_1){
        qDebug() <<"use ruling intersection";
        if(ui->glWid2dim->model.back()->FL.empty() || ui->glWid2dim->model.back()->FL[0]->FoldingCurve.empty())return;
        double tol = ui->TolValue->value();
        ui->glWid2dim->model.back()->BendingModel(wb, wp, warea, wsim, dim, tol, bndrange, layerNum, 0, IsStartEnd, OptimizeAngleFor3Rulings);
        fold_Sm();
        qDebug()<<"/////////////////////////";
    }
    if(e->key() == Qt::Key_2){
        qDebug()<<"use regression curve and triangle area";
        if(ui->glWid2dim->model.back()->FL.empty() || ui->glWid2dim->model.back()->FL[0]->FoldingCurve.empty())return;
        double tol = ui->TolValue->value();
        ui->glWid2dim->model.back()->BendingModel(wb, wp, warea, wsim, dim, tol, bndrange, layerNum, 1, IsStartEnd, OptimizeAngleFor3Rulings);
        fold_Sm();
        qDebug()<<"/////////////////////////";
    }
    if(e->key() == Qt::Key_3){
        qDebug()<<"vertices moving";
        if(e->modifiers().testFlag(Qt::ControlModifier)){
            ui->glWid2dim->model.back()->BendingModel(0, 0, warea, wsim, dim, 0, bndrange, layerNum, 12, IsStartEnd, OptimizeAngleFor3Rulings);
        }
        else  ui->glWid2dim->model.back()->BendingModel(0, 0, warea, wsim, dim, 0, bndrange, layerNum, 2, IsStartEnd, OptimizeAngleFor3Rulings) ;
        qDebug()<<"/////////////////////////";
    }
    if(e->key() == Qt::Key_4){
        qDebug()<<"propagate optimization vertex points from center to end point ";
        if(e->modifiers().testFlag(Qt::ControlModifier)){
            ui->glWid2dim->model.back()->BendingModel(0, 0, warea, wsim, dim, 0, bndrange, layerNum, 13, IsStartEnd, OptimizeAngleFor3Rulings);
        }
        else  ui->glWid2dim->model.back()->BendingModel(0, 0, warea, wsim, dim, 0, bndrange, layerNum, 3, IsStartEnd, OptimizeAngleFor3Rulings);
        qDebug()<<"/////////////////////////";
    }

    if(e->key() == Qt::Key_5){
        ui->glWid2dim->model.back()->BendingModel(0, 0, warea, wsim, dim, 0, bndrange, layerNum, 4, IsStartEnd, OptimizeAngleFor3Rulings);
    }
    if(e->key() == Qt::Key_6){
        ui->glWid2dim->model.back()->BendingModel(0, 0, warea, wsim, dim, 0, bndrange, layerNum, 5, IsStartEnd, OptimizeAngleFor3Rulings);
    }
    if(e->key() == Qt::Key_7){
        ui->glWid2dim->model.back()->BendingModel(0, 0, warea, wsim, dim, 0, bndrange, layerNum, 6, IsStartEnd, OptimizeAngleFor3Rulings);
    }
    if(e->key() == Qt::Key_8){
        if(e->modifiers().testFlag(Qt::ControlModifier)){
            ui->glWid2dim->model.back()->BendingModel(0, 0, warea, wsim, dim, 0, bndrange, layerNum, 7, IsStartEnd, OptimizeAngleFor3Rulings);
        }
        else  ui->glWid2dim->model.back()->BendingModel(0, 0, warea, wsim, dim, 0, bndrange, layerNum, 17, IsStartEnd, OptimizeAngleFor3Rulings);
    }
    if(e->key() == Qt::Key_Return){//
        ui->glWid2dim->stashcurrentstate();
    }

    if(e->key() == Qt::Key_A){
        if(!e->modifiers().testFlag(Qt::ControlModifier))ui->glWid2dim->switch2AffinMode();
        else ui->glWid2dim->switch2VisibleCurve();
    }
    if(e->key() == Qt::Key_B){
        qDebug()<<"if optimization of flap angle failed, optimize only 3 rulings and flap angle";
        OptimizeAngleFor3Rulings = !OptimizeAngleFor3Rulings;
    }
    if(e->key() == Qt::Key_C){
        if(e->modifiers().testFlag(Qt::ControlModifier))ui->glWid2dim->CopyCurveObj();
        else{
            if(!ui->glWid2dim->model.back()->outline->IsClosed())ui->glWid3dim->setVertices();
            else if(!ui->glWid2dim->model.back()->FL.empty()){
                ui->glWid2dim->model.back()->FL[0]->ReassignColor();
                ui->glWid2dim->model.back()->deform();
                ui->glWid2dim->model.back()->modifyFoldingCurvePositionOn3d();
                ui->glWid2dim->update();
                ui->glWid3dim->setVertices(ui->glWid2dim->model.back()->outline->Lines, ui->glWid2dim->model.back()->Rulings, ui->glWid2dim->model.back()->FL, ui->glWid2dim->AllRulings);
            }
            else ui->glWid3dim->setVertices(ui->glWid2dim->model.back()->outline->Lines, ui->glWid2dim->model.back()->Rulings, ui->glWid2dim->model.back()->FL, ui->glWid2dim->AllRulings);
        }

    }

    if(e->key() == Qt::Key_D){
        DebugMode::Singleton::getInstance().switchval();
        if(DebugMode::Singleton::getInstance().isdebug())qDebug() << "DebugMode On";
        else qDebug() << "DebugMode off";
    }

    if(e->key() == Qt::Key_E){
        if(e->modifiers().testFlag(Qt::ControlModifier)){
            qDebug()<<"optimization base";
            ui->glWid2dim->model.back()->Modify4LastFoldLine(ui->glWid2dim->model.back()->FL.back(), warea, wsim, bndrange, 0, IsStartEnd);
        }
        else{
            qDebug()<<"fullsearch";
            ui->glWid2dim->model.back()->Modify4LastFoldLine(ui->glWid2dim->model.back()->FL.back(), warea, wsim, bndrange, 1, IsStartEnd);
        }
    }

    if(e->key() == Qt::Key_N){
        ui->glWid2dim->InitializeDrawMode();
    }

    if(e->key() == Qt::Key_M){
        if(ui->glWid2dim->model.back()->FL.empty())return;
        double tol = ui->TolValue->value();
        if(e->modifiers().testFlag(Qt::ControlModifier)){
            auto NewFL = ui->glWid2dim->model.back()->FL.back();
            ui->glWid2dim->model.back()->AddNewFoldLine(NewFL);
        }
        else ui->glWid2dim->model.back()->AssignRuling(3, tol, false);
        auto NewFL = ui->glWid2dim->model.back()->FL.back();
        ui->glWid2dim->model.back()->AddNewFoldLine(NewFL);
        ui->glWid2dim->update();
        fold_Sm();
    }

    if( e->key() == Qt::Key_P){
        if(!ui->glWid2dim->model.back()->outline->IsClosed())ui->glWid3dim->setVertices();
        else if(!ui->glWid2dim->model.back()->FL.empty()){
             ui->glWid3dim->setVertices(ui->glWid2dim->model.back()->outline->Lines, ui->glWid2dim->model.back()->Rulings, ui->glWid2dim->model.back()->FL, ui->glWid2dim->AllRulings);
        }
        else{
             ui->glWid3dim->setVertices(ui->glWid2dim->model.back()->outline->Lines, ui->glWid2dim->model.back()->Rulings, ui->glWid2dim->model.back()->FL, ui->glWid2dim->AllRulings);
        }
    }

    if(e->key() == Qt::Key_Q)IsStartEnd = !IsStartEnd;

    if(e->key() == Qt::Key_R){
        if(e->modifiers().testFlag(Qt::ControlModifier))ui->glWid2dim->model.back()->BendingModel(0, 0, warea, wsim, 3, 0, bndrange, layerNum, -1, IsStartEnd, OptimizeAngleFor3Rulings);//initialization
        else{
             double a = static_cast<double>(ui->angleSlider->value())/100.0;
             std::vector<std::shared_ptr<Vertex>> Poly_V = ui->glWid2dim->model.back()->outline->getVertices();
             std::vector<std::vector<std::shared_ptr<Vertex>>> Tri_fixside;
             auto Triangles = ui->glWid2dim->model.back()->FL.back()->CalclateRegressionCurve(a,Poly_V,true, IsStartEnd, Tri_fixside);
             qDebug() <<"csv file export";
        }
    }

    if(e->key() == Qt::Key_S){
        if(e->modifiers().testFlag(Qt::ControlModifier)){
             exportobj();
        }
    }

    if(e->key() == Qt::Key_T){
        static int AlgType = 0;
        std::vector<std::shared_ptr<Vertex>> Poly_V = ui->glWid2dim->model.back()->outline->getVertices();
        ui->glWid2dim->model.back()->FL.back()->ReviseCenterFlapAngle(Poly_V, IsStartEnd, AlgType);
        AlgType = (AlgType + 1) % 2;
    }

    if(e->key() == Qt::Key_V){
        if(e->modifiers().testFlag(Qt::ControlModifier)) ui->glWid2dim->PasteCurveObj();
    }

    if(e->key() == Qt::Key_W){
        auto NewFL = ui->glWid2dim->model.back()->FL.back();
        ui->glWid2dim->model.back()->Interpolation(NewFL);
    }

    if(e->key() == Qt::Key_Z){
        if(e->modifiers().testFlag(Qt::ControlModifier)) ui->glWid2dim->back2befstate();
    }

    fold_Sm();
    ui->glWid2dim->update();
    ui->glWid3dim->update();
}

void MainWindow::mouseMoveEvent(QMouseEvent *e){
    if(SelectedBtn == nullptr)return;
    std::vector<std::shared_ptr<Btn4Crv>>::iterator itr = std::find(LayerList.begin(), LayerList.end(), SelectedBtn);
    if(itr == LayerList.end())return;
    int n = std::distance(LayerList.begin(),itr);
    QPoint p = ui->LayerListWidget->mapFromGlobal(QCursor::pos());
    int pad = 5, btn_h = SelectedBtn->height(), top = (p.y() < dragPos.y()) ? 0 : p.y() - dragPos.y();

    SelectedBtn->move(pad, top);
    int bottom = SelectedBtn->y() + SelectedBtn->height();
    for(int i = 0; i < (int)LayerList.size(); i++){
        if(LayerList[i] == SelectedBtn)continue;
        if((top < LayerList[i]->geometry().y() && originalPos.y() > LayerList[i]->geometry().y())
                || (LayerList[i]->geometry().y() + LayerList[i]->geometry().height() < bottom && originalPos.y() < LayerList[i]->geometry().y())){
            std::iter_swap(LayerList.begin() + i, LayerList.begin() + n);
            emit swapIndex(n, i);
            n = i;
            SelectedBtn = LayerList[n];
            for(int j = 0; j < (int)LayerList.size(); j++){
                if(j == n)continue;
                LayerList[j]->move(pad, btn_h + j * (btn_h + pad));
            }

        }
    }
}

void MainWindow::exportsvg(QString filename){
    auto getcolor = [](double c, double a, double y)->double{
        if(y < a)return c/a * y/255.0;
        return ((255.0 - c)*(y - a)/(std::numbers::pi - a) + c)/255.0;
    };

    auto filesvg = filename.replace(".obj", ".svg");
    QStringList WriteList;
    WriteList.append("<svg xmlns=\"http://www.w3.org/2000/svg\">\n");
    WriteList.append("<g>");

    std::string strokewidth = "\"0.5\"", str, strokecolor;
    //折曲線と交点を持たないruling
    for(auto itr_r = ui->glWid2dim->model.back()->Rulings.begin(); itr_r != ui->glWid2dim->model.back()->Rulings.end(); itr_r++){
        if((*itr_r)->hasCrossPoint)continue;

        if((*itr_r)->et == EdgeType::r){
            double color = getcolor(ui->glWid2dim->model.back()->ColorPt.color, ui->glWid2dim->model.back()->ColorPt.angle, (*itr_r)->color/255.0);
            if(color > 1e-3){//mount
                strokecolor = "\"red\"";
            }else if(color < -1e-3){//valley
                strokecolor = "\"blue\"";
            }else{
                strokecolor ="\"black\"";
            }
            double x1 = (*itr_r)->o->p.x(), y1 = (*itr_r)->o->p.y(), x2 = (*itr_r)->v->p.x(), y2 = (*itr_r)->v->p.y();
            str = "<line x1 = \"" + std::to_string(x1) + "\" y1 = \"" + std::to_string(y1) + "\" x2= \"" + std::to_string(x2) + "\" y2= \"" + std::to_string(y2)+
                  "\" fill = \"none\" stroke=" + strokecolor + " stroke-width=" + strokewidth + "/>\n";
            WriteList.append(QString::fromStdString(str));
        }
    }

    //可展面の輪郭描画
    {
        strokecolor ="\"black\"";
        auto  Vertices = ui->glWid2dim->model.back()->outline->getVertices();
        int vsize = Vertices.size();
        for(int i = 0; i < vsize; i++){
            double x1 = Vertices[i]->p.x(), y1 = Vertices[i]->p.y(), x2 = Vertices[(i + 1) % vsize]->p.x(), y2 = Vertices[(i + 1) % vsize]->p.y();
            str = "<line x1 = \"" + std::to_string(x1) + "\" y1 = \"" + std::to_string(y1) + "\" x2= \"" + std::to_string(x2) + "\" y2= \"" + std::to_string(y2)+
                  "\" fill = \"none\" stroke=" + strokecolor + " stroke-width=" + strokewidth + "/>\n";
            WriteList.append(QString::fromStdString(str));
        }
    }

    //折曲線上の4価頂点の描画
    for(auto& FL: ui->glWid2dim->model.back()->FL){
        if(FL->FoldingCurve.empty())continue;
        std::vector<std::shared_ptr<Vertex4d>> ValidFC;
        for(auto&fc: FL->FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}

        //p0:描画するエッジの始点, p1: 描画するエッジの端点, p2: 片側の平面上の点, p3: もう片方の平面上の点
        auto DrawEdge = [&](const std::shared_ptr<Vertex>& p0, const std::shared_ptr<Vertex>& p1, const std::shared_ptr<Vertex>& p2, const std::shared_ptr<Vertex>& p3, bool IsReverse){
            Eigen::Vector3d SpinAxis = (p1->p3 - p0->p3).normalized();
            Eigen::Vector3d f_nv = (SpinAxis.cross(p2->p3 - p0->p3)).normalized(),fp_nv = ((p3->p3 - p0->p3).cross(SpinAxis)).normalized();
            double f = (!IsReverse)?SpinAxis.dot(f_nv.cross(fp_nv)) :-SpinAxis.dot(f_nv.cross(fp_nv));
            double th = -1e-5;
            if(f <th)strokecolor = "\"red\"";//mount
            else if(f > -th)strokecolor = "\"blue\"";//valley
            else strokecolor ="\"black\"";
            double x1 = p1->p.x(), y1 = p1->p.y(), x2 = p0->p.x(), y2 = p0->p.y();
            str = "<line x1 = \"" + std::to_string(x1) + "\" y1 = \"" + std::to_string(y1) + "\" x2= \"" + std::to_string(x2) + "\" y2= \"" + std::to_string(y2)+
                  "\" fill = \"none\" stroke=" + strokecolor + " stroke-width=" + strokewidth + "/>\n";
            WriteList.append(QString::fromStdString(str));
        };

        for(int i = 1; i < (int)ValidFC.size()-1; i++)
            DrawEdge(ValidFC[i]->first, ValidFC[i-1]->first, ValidFC[i]->second, ValidFC[i]->third, false);
        DrawEdge(ValidFC.end()[-2]->first, ValidFC.back()->first, ValidFC.end()[-2]->second, ValidFC.end()[-2]->third, true);
        for(int i = 1; i < (int)ValidFC.size() - 1; i++){
            DrawEdge(ValidFC[i]->first, ValidFC[i]->second, ValidFC[i+1]->first, ValidFC[i-1]->first, false);
            DrawEdge(ValidFC[i]->first, ValidFC[i]->third, ValidFC[i-1]->first, ValidFC[i+1]->first, false);
        }
    }
    WriteList.append("</g>\n");
    WriteList.append("</svg>\n");

    QFile file(filesvg);
    if (!file.open(QIODevice::WriteOnly)) {
        QMessageBox::information(this, tr("Unable to open file"), file.errorString());
        return;
    }

    QTextStream out(&file);
    foreach (QString item, WriteList) {
        out << item;
    }
    return;
}

//https://nprogram.hatenablog.com/entry/2017/08/10/082236
void MainWindow::exportobj(){
    if(!ui->glWid2dim->model.back()->outline->IsClosed()){
        QMessageBox msgBox;
        msgBox.setText("There is no developable surface.");
        msgBox.exec();
        return;
    }
    QStringList WriteList, WriteList_tri;
    std::vector<Eigen::Vector3d> Normals;
    Eigen::Vector3d befN = {0,0,0}, N;
    Eigen::Transform<double, 3, Eigen::Affine> Mirror;
    Mirror = Eigen::Matrix3d::Identity(); Mirror(1, 1) = -1;
    Eigen::Matrix3d m = Eigen::Matrix3d::Identity(); m(1,1) = -1;  Mirror = m;
    std::vector<std::vector<Eigen::Vector3d>> Vertices;
    std::vector<std::vector<std::shared_ptr<Vertex>>> Polygons;
    std::vector<std::shared_ptr<Vertex>> polygon;
    std::vector<double> planerity_value;

    auto Planerity  = [](const std::vector<std::shared_ptr<Vertex>>& vertices, const std::vector<std::shared_ptr<Line>>& lines)->double{
        if((int)vertices.size() <= 3)return 0.0;
       else{
           std::vector<Eigen::Vector3d> QuadPlane;
           for(auto&v: vertices){
               bool IsOutlineVertices = false;
               for(auto&l: lines){
                   if(l->v->p3 == v->p3)IsOutlineVertices = true;
               }
               if(!IsOutlineVertices)QuadPlane.push_back(v->p3);
           }
           if((int)QuadPlane.size() <= 3)return 0.0;
           double l_avg = ((QuadPlane[0] - QuadPlane[2]).norm() + (QuadPlane[1] - QuadPlane[3]).norm())/2.0;
           double d;
           Eigen::Vector3d u1 = (QuadPlane[0] - QuadPlane[2]).normalized(), u2 = (QuadPlane[1]-  QuadPlane[3]).normalized();
           if((u1.cross(u2)).norm() < 1e-9){
               Eigen::Vector3d H = QuadPlane[3] + u2.dot(QuadPlane[1] - QuadPlane[3])*u2;
               d = (H - QuadPlane[1]).norm();
           }else{
               auto u = u1.cross(u2);
               d = (u.dot(QuadPlane[2] - QuadPlane[3]))/(u1.cross(u2)).norm();
           }
           return d/l_avg;
       }
   };

    std::vector<std::vector<std::shared_ptr<Vertex4d>>> FoldingCurves;
    for(auto&FldCrv: ui->glWid2dim->model.back()->FL)FoldingCurves.push_back(FldCrv->FoldingCurve);
    Polygons = MakeModel(ui->glWid2dim->model.back()->outline->Lines, ui->glWid2dim->model.back()->Rulings, FoldingCurves);

    for(auto&poly: Polygons){
        std::vector<Eigen::Vector3d> V;
        auto poly_sort = SortPolygon(poly);
        for(auto&p: poly_sort)V.push_back(Mirror * p->p3);
        Vertices.push_back(V);
        planerity_value.push_back(Planerity(poly, ui->glWid2dim->model.back()->outline->Lines));
    }

    for(const auto& mesh: Vertices){
        for(const auto&v: mesh)WriteList.append("v " + QString::number(v.x()) + " " + QString::number(v.y()) + " " + QString::number(v.z()) + "\n");
        N = ((mesh[1] - mesh[0]).cross(mesh[2] - mesh[0])).normalized();
        Normals.push_back(N);
        befN = N;
    }
    for(const auto&n : Normals) WriteList.append("vn " + QString::number(n.x()) + " " + QString::number(n.y()) + " " + QString::number(n.z()) + "\n");

    int cnt = 1;
    for(int i = 0; i < (int)Vertices.size(); i++){
        QString s = "f ";
        for(const auto& v: Vertices[i]){
            s += QString::number(cnt) + "//" + QString::number(i+1) + " ";
            cnt++;
        }
        s += "\n";
        WriteList.append(s);
    }
    /*
    Eigen::Vector3d Zaxis(0,0,1);
    std::vector<std::array<Eigen::Vector3d, 3>> TriMeshs;
    auto Triangulation = [](std::vector<std::shared_ptr<Vertex>>&input, std::vector<std::array<Eigen::Vector3d, 3>>&output){
        output.clear();
        std::vector<std::shared_ptr<Vertex>> Edges;
        std::copy(input.begin(), input.end(), back_inserter(Edges) );
        int n = Edges.size();
        while(Edges.size() >= 3){
            n = Edges.size();
            for(int i = 0; i < n; i++){
               int prev = (n + i - 1) % n;
               int next = (n + i + 1) % n;
               std::array<std::shared_ptr<Vertex>, 3> tri = {Edges[prev], Edges[i], Edges[next]};
               bool elimTriMesh = true;
               for(int j = 0; j < n - 3; j++){
                   Eigen::Vector3d p = Edges[(next + 1 + j) % n]->p;
                   auto _tri = std::array{tri[0]->p, tri[1]->p, tri[2]->p};
                   bool check1 = MathTool::hasPointInTriangle3D(p, _tri);
                   bool check2 = MathTool::IsAngleLessThan180(tri[1]->p, tri[0]->p, tri[2]->p);
                   if(check1 || !check2){
                       elimTriMesh = false;
                       break;
                   }
               }
               if(elimTriMesh){
                   Edges.erase(Edges.begin() + i);
                   Eigen::Vector3d N = ((tri[1]->p3 - tri[0]->p3).cross(tri[2]->p3 -tri[0]->p3)).normalized();
                   if(N.dot(Eigen::Vector3d::UnitZ()) < 0){ std::swap(tri[0], tri[2]);}
                   output.push_back({tri[0]->p3, tri[1]->p3, tri[2]->p3});
                   break;
               }
            }

        }
    };

    for(auto&poly: Polygons){
        std::vector<std::array<Eigen::Vector3d, 3>> trimesh;
        Triangulation(poly, trimesh);
        TriMeshs.insert(TriMeshs.end(), trimesh.begin(), trimesh.end());
    }
    for(auto& mesh: TriMeshs){
        for(const auto&v: mesh)WriteList_tri.append("v " + QString::number(v.x()) + " " + QString::number(v.y()) + " " + QString::number(v.z()) + "\n");
        N = ((mesh[1] - mesh[0]).cross(mesh[2] - mesh[0])).normalized();
        Normals.push_back(N);
    }
    for(const auto&n : Normals) WriteList_tri.append("vn " + QString::number(n.x()) + " " + QString::number(n.y()) + " " + QString::number(n.z()) + "\n");
    cnt = 1;
    for(int i = 0; i < (int)TriMeshs.size(); i++){
        QString s = "f ";
        for(int j = 0; j < 3; j++)s += QString::number(cnt++) + "//" + QString::number(i+1) + " ";
        s += "\n";
        WriteList_tri.append(s);
    }
    */
    const QString DirName = "./OBJ";
    const QDir dir; dir.mkdir(DirName);
    QDir CurDir = QDir::current(); CurDir.cd(DirName);
    const QString CSVDir = "./ResultCSV";
    const QDir ChildDir; ChildDir.mkdir(CSVDir);
    QString fileName = QFileDialog::getSaveFileName(this, tr("save as obj"), "", tr("テキストファイル(*.obj)") );
    if(fileName.isEmpty())return;

    QFile file(fileName);
    if (!file.open(QIODevice::WriteOnly)) {
            QMessageBox::information(this, tr("Unable to open file"), file.errorString());
            return;
    }

    QTextStream out(&file);
    foreach (QString item, WriteList) {
        out << item;
    }
    exportsvg(fileName);

    CurDir.cd(CSVDir);
    //平面性の結果
    QFile QuantitativeResult_file(fileName.split(u'.')[0] + "QuantitativeResult.csv");
    if (!QuantitativeResult_file.open(QIODevice::WriteOnly)) {
            QMessageBox::information(this, tr("Unable to open file"), file.errorString());
            return;
    }

    QTextStream QuantitativeResult(&QuantitativeResult_file);
    QuantitativeResult << "Planarity\n" ;
    for(auto&c: planerity_value)QuantitativeResult << c << ", ";
    QuantitativeResult << "\nDevelopability\n" ;
    if(!(ui->glWid2dim->model.back()->FL.empty() || ui->glWid2dim->model.back()->FL[0]->FoldingCurve.empty())){
        for(auto&FL: ui->glWid2dim->model.back()->FL){
            if(FL->FoldingCurve.empty())continue;
            for(auto itr = FL->FoldingCurve.begin() + 1; itr != FL->FoldingCurve.end() - 1; itr++){
                Eigen::Vector3d et = (*itr)->second->p3 - (*itr)->first->p3, er = (*(itr - 1))->first->p3 - (*itr)->first->p3, eb = (*itr)->third->p3 - (*itr)->first->p3, el = (*(itr + 1))->first->p3 - (*itr)->first->p3;
                et = et.normalized(); er = er.normalized(); eb = eb.normalized(); el = el.normalized();
                double phi1 = std::acos(et.dot(er)), phi2 = std::acos(et.dot(el)), phi3 = std::acos(eb.dot(el)), phi4 = std::acos(eb.dot(er));
                double a = abs(2.0*std::numbers::pi - (phi1 + phi2 + phi3 + phi4));
                if(a != -1)QuantitativeResult << a << ", ";
            }
        }
    }
    QuantitativeResult << "\nPlanarity(adjacent ruling)\n" ;
    for(auto&FL: ui->glWid2dim->model.back()->FL){
        if(FL->FoldingCurve.empty())continue;
        for(auto itr = FL->FoldingCurve.begin() + 1; itr != FL->FoldingCurve.end() - 2; itr++){
            double l_avg = (((*itr)->first->p3 - (*(itr+1))->second->p3).norm() + ((*itr)->second->p3 - (*(itr+1))->first->p3).norm())/2.0;
            double d;
            Eigen::Vector3d u1 = ((*itr)->first->p3 - (*(itr+1))->second->p3).normalized(), u2 = ((*itr)->second->p3-  (*(itr+1))->first->p3).normalized();
            if((u1.cross(u2)).norm() < 1e-9){
                Eigen::Vector3d H = (*(itr+1))->first->p3 + u2.dot((*itr)->second->p3 - (*(itr+1))->first->p3)*u2;
                d = (H - (*itr)->second->p3).norm();
            }else{
                auto u = u1.cross(u2);
                d = (u.dot((*(itr+1))->second->p3 - (*(itr+1))->first->p3))/(u1.cross(u2)).norm();
            }
            QuantitativeResult << d/l_avg << ", ";
        }
    }

}

void MainWindow::addCurveBtn(){
    QRect geo = ui->LayerListWidget->geometry();
    int pad = 5;
    int btn_w = geo.width() - 2 * pad, btn_h = 25;
    QString text;
    if(ui->glWid2dim->model.back()->crvs[0]->getCurveType() == CurveType::bezier3)text = "Bezier" + QString::number(++CurvesNum[0]);
    else if(ui->glWid2dim->model.back()->crvs[0]->getCurveType()  == CurveType::bsp3)text = "Bspline" + QString::number(++CurvesNum[1]);
    else if(ui->glWid2dim->model.back()->crvs[0]->getCurveType() == CurveType::line)text = "Line" + QString::number(++CurvesNum[2]);
    else if(ui->glWid2dim->model.back()->crvs[0]->getCurveType() == CurveType::arc)text = "Arc"+QString::number(++CurvesNum[3]);

    std::shared_ptr<Btn4Crv> newbtn = std::make_shared<Btn4Crv>(ui->glWid2dim->model.back()->crvs[0],text, ui->LayerListWidget);
    connect(newbtn.get(), &Btn4Crv::clicked, this, &MainWindow::SetHandleCrv);
    for(int i = 0; i < (int)LayerList.size(); i++) LayerList[i]->setGeometry(pad, btn_h + (i + 1) * (btn_h + pad), btn_w, btn_h);
    LayerList.insert(LayerList.begin(), newbtn);

    newbtn->setGeometry(pad, btn_h, btn_w, btn_h);
    newbtn->show();

}

void MainWindow::RemoveBtnFromLayerCrv(std::vector<int> deleteIndex){
    if(LayerList.empty())return;
    int pad = 5, btn_h = LayerList[0]->size().height();
    for(auto& n : deleteIndex){
        if(n == -1)continue;
        for(int i = n + 1; i < (int)LayerList.size(); i++) LayerList[i]->move(pad, btn_h + (i - 1) * (btn_h + pad));
        LayerList[n]->hide();
        LayerList.erase(LayerList.begin() + n);
        LayerList.shrink_to_fit();
    }
}

void MainWindow::SetHandleCrv(Btn4Crv *btn, QMouseEvent *e){
    if(LayerList.empty())return;

    int ind;
    if(btn == nullptr){
        ind = -1;
        int pad = 5, btn_h = LayerList[0]->height();
        for(int i = 0; i < (int)LayerList.size(); i++) LayerList[i]->move(pad, btn_h + i * (btn_h + pad));
        SelectedBtn = nullptr;
    }
    else{
        if(e->button() == Qt::LeftButton){
            dragPos = btn->mapFromGlobal(QCursor::pos());
            originalPos = btn->geometry();
            std::vector<std::shared_ptr<Btn4Crv>>::iterator itr = std::find_if(LayerList.begin(), LayerList.end(), [&btn](std::shared_ptr<Btn4Crv>& b){return b.get() == btn;});
            if(itr == LayerList.end())ind = -1;
            else{
                ind = std::distance(LayerList.begin(), itr);
                SelectedBtn = *itr;
            }
        }
        else if(e->button() == Qt::RightButton){
            ind = -1;
        }else ind = -1;
    }
    emit signalNewSelectedCrv(ind);
}

void MainWindow::sendCurveType(CurveType &ct){
    int n = ui->CurveTypeBox->currentIndex();
    switch(n){
    case 0:
        ct = CurveType::line;
        break;
    case 1:
        ct = CurveType::bsp3;
        break;
    case 2:
        ct = CurveType::arc;
        break;
    default:
        ct = CurveType::none;
        break;
    }
}
