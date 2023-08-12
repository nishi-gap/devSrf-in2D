#include "widget.h"
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
    connect(ui->TolValue, &QDoubleSpinBox::valueChanged, this, &MainWindow::changeToleranceValue_Spin);
    connect(ui->ReassinColorButton, &QPushButton::clicked, this, &MainWindow::ReassinColor);

    //fold line debug
    connect(ui->startButton, &QPushButton::clicked, ui->glWid2dim, &GLWidget_2D::Start4Debug_CF);
    connect(ui->developablityButton, &QCheckBox::clicked, ui->glWid2dim, &GLWidget_2D::checkDevelopability);

    connect(ui->angleSlider, &QSlider::sliderMoved, this, &MainWindow::changeAngleFromSlider);
    connect(ui->angleA, &QDoubleSpinBox::valueChanged, this, &MainWindow::changeAngleFromSpinBox);
    connect(this, &MainWindow::sendAngle, ui->glWid2dim, &GLWidget_2D::changeBetaValue);

    connect(ui->DebugWindow, &QPushButton::clicked,ui->glWid2dim,&GLWidget_2D::OpenDebugWindwow);
    connect(ui->SaveButton, &QPushButton::clicked,this, &MainWindow::exportobj);

    LayerList.clear();
    connect(ui->glWid2dim, &GLWidget_2D::deleteCrvSignal, this, &MainWindow::RemoveBtnFromLayerCrv);
    connect(this, &MainWindow::signalNewSelectedCrv, ui->glWid2dim, &GLWidget_2D::changeSelectedCurve);
    connect(this, &MainWindow::swapIndex, ui->glWid2dim, &GLWidget_2D::swapCrvsOnLayer);


    //planarity
    connect(ui->SmoothingButton, &QPushButton::clicked, this, &MainWindow::StartSmoothingSurface);
    connect(ui->SimpleSmoothButton, &QPushButton::clicked, this, &MainWindow::SimpleSmoothing);
    connect(ui->PlanarityButton, &QCheckBox::clicked, ui->glWid3dim, &GLWidget_3D::PlanarityDispay);
    connect(ui->OptPlararity_Button, &QPushButton::clicked, this, &MainWindow::StartOptimization_plararity);

    //optimization or discrete developable surface
    connect(ui->OptBtn, &QPushButton::clicked, this, &MainWindow::StartOptimization);

    //Erase Non Fold Edge
    connect(ui->EraseNonFoldButton, &QCheckBox::clicked, this, &MainWindow::EraseNonFoldEdge);

    CurvesNum[0] = MainWindow::CurvesNum[1] = MainWindow::CurvesNum[2] = MainWindow::CurvesNum[3] = 0;
    SelectedBtn = nullptr;
}

static int keyType = 0;

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::SymmetricConstraint(){ emit constraintType(0);}

void MainWindow::Initialize(){
    for(auto&v: LayerList)delete v;
    LayerList.clear();
    update();
}
void MainWindow::ChangeMaxColor(int val){model->SetMaxFold((double)val);}

void MainWindow::addFoldLine_l(){emit signalFLtype(PaintTool::FoldLine_test);}
void MainWindow::addFoldLine_arc(){emit signalFLtype(PaintTool::FoldLine_arc);}
void MainWindow::addFoldLine_bezier(){emit signalFLtype(PaintTool::FoldLine_bezier);}
void MainWindow::addFoldLine_test(){emit signalFLtype(PaintTool::FoldLine_test);}
void MainWindow::moveCtrlPts_fl(){emit signalFLtype(PaintTool::FoldLine_move);}
void MainWindow::color_FL(){emit signalFLtype(PaintTool::FoldLine_move);}

void MainWindow::ModelBack(){
    for(auto&v: ui->glWid2dim->model->vertices){
        v->p = v->p2_ori;
        v->p3 = v->p3_ori;
    }
    ui->glWid2dim->update();
    ui->glWid3dim->setVertices(ui->glWid2dim->model->outline->Lines, ui->glWid2dim->model->Rulings, ui->glWid2dim->model->FL, ui->glWid2dim->AllRulings);
}

void MainWindow::EraseNonFoldEdge(bool state){
    ui->glWid2dim->EraseNonFoldEdge(state);
    ui->glWid3dim->EraseNonFoldEdge(state);
    if(ui->glWid2dim->model->FL.empty() || ui->glWid2dim->model->FL[0]->FoldingCurve.empty())return;
    ui->glWid3dim->setVertices(ui->glWid2dim->model->outline->Lines, ui->glWid2dim->model->Rulings, ui->glWid2dim->model->FL, ui->glWid2dim->AllRulings);
}

void MainWindow::changeToleranceValue_Slider(int val){
    double maxSpin = ui->TolValue->maximum(), maxSlider = ui->ToleranceValue->maximum();
    if(ui->glWid2dim->model->FL.empty() || ui->glWid2dim->model->FL[0]->FoldingCurve.empty())return;
    if(!ui->BinaryMVColor->isChecked())ui->BinaryMVColor->setChecked(true);
    ui->glWid2dim->setNewGradationMode();
    ui->glWid2dim->VisualizeMVColor(true);
    double tol = maxSpin * (double)val/(double)maxSlider;
    ui->TolValue->setValue(tol);
    ui->glWid2dim->model->FL[0]->SimplifyModel(tol);
    ui->glWid2dim->update();
     ui->glWid3dim->setVertices(ui->glWid2dim->model->outline->Lines, ui->glWid2dim->model->Rulings, ui->glWid2dim->model->FL, ui->glWid2dim->AllRulings);
}

void MainWindow::changeToleranceValue_Spin(double val){
    double maxSpin = ui->TolValue->maximum(), maxSlider = ui->ToleranceValue->maximum();
    if(ui->glWid2dim->model->FL.empty() || ui->glWid2dim->model->FL[0]->FoldingCurve.empty())return;
    if(!ui->BinaryMVColor->isChecked())ui->BinaryMVColor->setChecked(true);
    ui->glWid2dim->setNewGradationMode();
    ui->glWid2dim->VisualizeMVColor(true);

    ui->ToleranceValue->setValue(int(val/maxSpin * maxSlider));
    ui->glWid2dim->model->FL[0]->SimplifyModel(val);
    fold_Sm();
}


void MainWindow::StartOptimization(){
    if(ui->glWid2dim->model->FL.empty() || ui->glWid2dim->model->FL[0]->FoldingCurve.empty())return;
    auto Poly_V = ui->glWid2dim->model->outline->getVertices();
    double wb = ui->BendWeightButton->value(), wp = ui->ParalellWeightButton->value();
    bool res = ui->glWid2dim->model->FL[0]->Optimization_FlapAngle(ui->glWid2dim->model->Rulings, ui->glWid2dim->model->vertices, Poly_V, wb, wp);
    if(res)fold_Sm();

}

void MainWindow::StartOptimization_plararity(){
    if(ui->glWid2dim->model->FL.empty() || ui->glWid2dim->model->FL[0]->FoldingCurve.empty())return;
    auto Poly_V = ui->glWid2dim->model->outline->getVertices();
    auto res = ui->glWid2dim->model->FL[0]->Optimization_PlanaritySrf(Poly_V);
    fold_Sm();
    ui->glWid3dim->ReceiveParam(res);
}

void MainWindow::StartSmoothingSurface(){
    auto Poly_V = ui->glWid2dim->model->outline->getVertices();
    if(ui->glWid2dim->model->FL.empty() || ui->glWid2dim->model->FL[0]->FoldingCurve.empty())return;
    bool ICE = ui->ConnectEndPointBox->isChecked();
    auto res = ui->glWid2dim->model->FL[0]->Optimization_SmooothSrf(Poly_V, ICE);
    if(!res.empty()){
        //ui->glWid2dim->update();
        //ui->glWid3dim->setVertices(ui->glWid2dim->model->Faces, ui->glWid2dim->model->outline->getVertices(), ui->glWid2dim->model->Edges, ui->glWid2dim->model->vertices, ui->glWid2dim->model->FL[0]->FoldingCurve, ui->glWid2dim->AllRulings);
        fold_Sm();
        ui->glWid3dim->ReceiveParam(res);
    }
}

void MainWindow::SimpleSmoothing(){
    auto Poly_V = ui->glWid2dim->model->outline->getVertices();
    if(ui->glWid2dim->model->FL.empty() || ui->glWid2dim->model->FL[0]->FoldingCurve.empty())return;
    bool res = ui->glWid2dim->model->FL[0]->SimpleSmooothSrf(Poly_V);
    if(res){
        ui->glWid2dim->update();
         ui->glWid3dim->setVertices(ui->glWid2dim->model->outline->Lines, ui->glWid2dim->model->Rulings, ui->glWid2dim->model->FL, ui->glWid2dim->AllRulings);
    }
}

void MainWindow::changeAngleFromSlider(int val){
    ui->angleA->setValue((double)val/100);
    double a = (double)val/18000.0 * std::numbers::pi;
    emit sendAngle(a, keyType);
}

void MainWindow::changeAngleFromSpinBox(double val){
    ui->angleSlider->setValue(val*100);
    double a = (double)val*std::numbers::pi/180.0;
    emit sendAngle(a, keyType);

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
    if(!ui->glWid2dim->model->outline->IsClosed())ui->glWid3dim->setVertices();
    else if(!ui->glWid2dim->model->FL.empty()){
         ui->glWid3dim->setVertices(ui->glWid2dim->model->outline->Lines, ui->glWid2dim->model->Rulings, ui->glWid2dim->model->FL, ui->glWid2dim->AllRulings);

    }
    else{
         ui->glWid3dim->setVertices(ui->glWid2dim->model->outline->Lines, ui->glWid2dim->model->Rulings, ui->glWid2dim->model->FL, ui->glWid2dim->AllRulings);
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
    if(!ui->glWid2dim->model->outline->IsClosed())ui->glWid3dim->setVertices();
    else if(!ui->glWid2dim->model->FL.empty() && !ui->glWid2dim->model->FL[0]->FoldingCurve.empty()){
        ui->glWid2dim->model->FL[0]->ReassignColor(ui->glWid2dim->model->Rulings, ui->glWid2dim->model->ColorPt);
        ui->glWid2dim->model->deform();
        ui->glWid2dim->update();
        ui->glWid2dim->model->modifyFoldingCurvePositionOn3d();
         ui->glWid3dim->setVertices(ui->glWid2dim->model->outline->Lines, ui->glWid2dim->model->Rulings, ui->glWid2dim->model->FL, ui->glWid2dim->AllRulings);
    }
    else{
         ui->glWid3dim->setVertices(ui->glWid2dim->model->outline->Lines, ui->glWid2dim->model->Rulings, ui->glWid2dim->model->FL, ui->glWid2dim->AllRulings);
    }
}

void MainWindow::keyPressEvent(QKeyEvent *e){
    static bool switchDraw = false;
    ui->glWid3dim->receiveKeyEvent(e);
    ui->glWid2dim->receiveKeyEvent(e);
    if(e->key() == Qt::Key_W){
        if(ui->glWid2dim->model->FL.empty() || ui->glWid2dim->model->FL[0]->FoldingCurve.empty())return;
        auto Poly_V = ui->glWid2dim->model->outline->getVertices();
        double wb = ui->BendWeightButton->value(), wp = ui->ParalellWeightButton->value();
        bool res = ui->glWid2dim->model->FL[0]->Optimization_FlapAngle(ui->glWid2dim->model->Rulings, ui->glWid2dim->model->vertices, Poly_V, wb, wp, false);
        fold_Sm();
    }
    else if(e->key() == Qt::Key_Return){emit PressedEnter();}
    else if(e->key() == Qt::Key_C){

        if(!ui->glWid2dim->model->outline->IsClosed())ui->glWid3dim->setVertices();
        else if(!ui->glWid2dim->model->FL.empty()){
            ui->glWid2dim->model->FL[0]->ReassignColor(ui->glWid2dim->model->Rulings, ui->glWid2dim->model->ColorPt);
            ui->glWid2dim->model->deform();
            ui->glWid2dim->model->modifyFoldingCurvePositionOn3d();
            ui->glWid2dim->update();
             ui->glWid3dim->setVertices(ui->glWid2dim->model->outline->Lines, ui->glWid2dim->model->Rulings, ui->glWid2dim->model->FL, ui->glWid2dim->AllRulings);
        }
        else{
             ui->glWid3dim->setVertices(ui->glWid2dim->model->outline->Lines, ui->glWid2dim->model->Rulings, ui->glWid2dim->model->FL, ui->glWid2dim->AllRulings);
        }
    }
    else if( e->key() == Qt::Key_P){

        if(!ui->glWid2dim->model->outline->IsClosed())ui->glWid3dim->setVertices();
        else if(!ui->glWid2dim->model->FL.empty()){
             ui->glWid3dim->setVertices(ui->glWid2dim->model->outline->Lines, ui->glWid2dim->model->Rulings, ui->glWid2dim->model->FL, ui->glWid2dim->AllRulings);
        }
        else{
             ui->glWid3dim->setVertices(ui->glWid2dim->model->outline->Lines, ui->glWid2dim->model->Rulings, ui->glWid2dim->model->FL, ui->glWid2dim->AllRulings);
        }
    }
    else if(e->key() == Qt::Key_3 ||e->key() == Qt::Key_4){
    }else if(e->modifiers().testFlag(Qt::ControlModifier)){
        if(e->key() == Qt::Key_S)exportobj();
    }
    else{

    }
    update();
}

void MainWindow::mouseMoveEvent(QMouseEvent *e){
    if(SelectedBtn == nullptr)return;
    std::vector<Btn4Crv*>::iterator itr = std::find(LayerList.begin(), LayerList.end(), SelectedBtn);
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

//https://nprogram.hatenablog.com/entry/2017/08/10/082236
void MainWindow::exportobj(){
    if(!ui->glWid2dim->model->outline->IsClosed()){
        QMessageBox msgBox;
        msgBox.setText("There is no developable surface.");
        msgBox.exec();
        return;
    }
    QStringList WriteList;
    std::vector<glm::f64vec3> Normals;
    glm::f64vec3 befN = {0,0,0}, N;
    glm::f64mat4x4 Mirror = glm::mat4(1.0f); Mirror[1][1] = -1;
   std::vector<std::vector<glm::f64vec3>> Vertices;
   std::vector<std::vector<Vertex*>> Polygons;
   std::vector<Vertex*> polygon;
   std::vector<double> planerity_value;

   auto Planerity  = [](const std::vector<glm::f64vec3>& vertices, const std::vector<Line*> lines)->double{
       if(vertices.size() == 3)return 0.0;
       else{
           std::vector<glm::f64vec3> QuadPlane;
           for(auto&v: vertices){
               bool IsOutlineVertices = false;
               for(auto&l: lines){
                   if(l->v->p3 == v)IsOutlineVertices = true;
               }
               if(!IsOutlineVertices)QuadPlane.push_back(v);
           }
           double l_avg = (glm::distance(QuadPlane[0], QuadPlane[2]) + glm::distance(QuadPlane[1], QuadPlane[3]))/2.0;
           double d;
           glm::f64vec3 u1 = glm::normalize(QuadPlane[0] - QuadPlane[2]), u2 = glm::normalize(QuadPlane[1]-  QuadPlane[3]);
           if(glm::length(glm::cross(u1, u2)) < DBL_EPSILON){
               glm::f64vec3 H = QuadPlane[3] + glm::dot(QuadPlane[1] - QuadPlane[3], u2)*u2;
               d = glm::distance(H, QuadPlane[1]);
           }else{
               d = glm::length(glm::dot(glm::cross(u1,u2),  QuadPlane[2] - QuadPlane[3]))/glm::length(glm::cross(u1, u2));
           }
           return d/l_avg;
       }
   };

   for(auto& l: ui->glWid2dim->model->outline->Lines) polygon.push_back(l->v);
   Polygons.push_back(polygon);

   if(!ui->glWid2dim->model->FL.empty() ||!ui->glWid2dim->model->FL[0]->FoldingCurve.empty()){
       for(auto&FC: ui->glWid2dim->model->FL){
           std::vector<Vertex4d> FldCrv = FC->FoldingCurve;
           glm::f64vec3 CrvDir = glm::normalize(FldCrv.back().first->p - FldCrv.front().first->p);
           for(auto&P: Polygons){
               int ind_fr = -1, ind_bc = -1;
               for(int i = 0; i < (int)P.size(); i++){
                   if(MathTool::is_point_on_line(FldCrv.front().first->p, P[i]->p, P[(i + 1) % (int)P.size()]->p))ind_fr = i;
                   if(MathTool::is_point_on_line(FldCrv.back().first->p, P[i]->p, P[(i + 1)  % (int)P.size()]->p))ind_bc = i;
               }
               if(ind_fr != -1 && ind_bc != -1){
                   std::vector<Vertex*> InsertedV, InsertedV_inv;
                   for(auto&v: FldCrv){InsertedV.push_back(v.first);InsertedV_inv.insert(InsertedV_inv.begin(), v.first);}

                   int i_min = std::min(ind_fr, ind_bc) + 1, i_max = std::max(ind_fr, ind_bc) + 1;
                   glm::f64vec3 Dir_prev = glm::normalize(P[i_min]->p - P[(i_min - 1) % (int)P.size()]->p);
                   std::vector<Vertex*> poly2 = {P.begin() + i_min, P.begin() + i_max};
                   P.erase(P.begin() + i_min, P.begin() + i_max);
                   if(glm::dot(glm::cross(CrvDir, Dir_prev), glm::f64vec3{0,0,1}) < 0){
                       P.insert(P.begin() + i_min, InsertedV.begin(), InsertedV.end());
                       poly2.insert(poly2.end(), InsertedV_inv.begin(), InsertedV_inv.end());
                   }else{
                       P.insert(P.begin() + i_min, InsertedV_inv.begin(), InsertedV_inv.end());
                       poly2.insert(poly2.end(), InsertedV.begin(), InsertedV.end());
                   }
                   Polygons.push_back(poly2);
                   break;
               }
           }
       }

       auto SplitPolygon = [](std::vector<std::vector<Vertex*>>& Polygons, Vertex *o, Vertex *v){//v:新たに挿入したいvertex, o:基本的にfirstを与える
           for(auto& Poly :Polygons){
               for(int i = 0; i < (int)Poly.size(); i++){
                   if(std::find(Poly.begin(), Poly.end(), v) != Poly.end())break;
                   if(!MathTool::is_point_on_line(v->p, Poly[i]->p, Poly[(i + 1) % (int)Poly.size()]->p))continue;
                   Poly.insert(Poly.begin() + i + 1, v);
                   break;
               }
               int f_ind = -1, s_ind = -1;
               for(int i = 0; i < (int)Poly.size(); i++){
                   if(Poly[i] == o)f_ind = i;
                   if(Poly[i] == v)s_ind = i;
               }
               if(f_ind == -1 || s_ind == -1)continue;
               int i_min = std::min(f_ind, s_ind), i_max = std::max(f_ind, s_ind);
               std::vector<Vertex*> poly2(Poly.begin() + i_min, Poly.begin() + i_max +1);
               Poly.erase(Poly.begin() + i_min + 1, Poly.begin() + i_max);
               Polygons.push_back(poly2);
               return;
           }
       };

       for(auto&FC: ui->glWid2dim->model->FL){
           for(auto itr = FC->FoldingCurve.begin() + 1; itr != FC->FoldingCurve.end() - 1; itr++){
               SplitPolygon(Polygons, itr->first, itr->second);
               SplitPolygon(Polygons, itr->first, itr->third);
           }
      }
   }else{
       for(auto itr_r = ui->glWid2dim->model->Rulings.begin(); itr_r != ui->glWid2dim->model->Rulings.end(); itr_r++){
           for(auto&P: Polygons){
               int vind = -1, oind = -1;
               for(int i = 0; i < (int)P.size(); i++){
                   //std::cout << glm::to_string((*itr_r)->o->p) << "  ,  "<< glm::to_string((*itr_r)->v->p) << std::endl;
                   if(MathTool::is_point_on_line((*itr_r)->o->p, P[i]->p, P[(i + 1) % (int)P.size()]->p))oind = i;
                   if(MathTool::is_point_on_line((*itr_r)->v->p, P[i]->p, P[(i + 1)  % (int)P.size()]->p))vind = i;
               }
               if(vind == -1 || oind == -1)continue;
               int i_min = std::min(vind, oind) + 1, i_max = std::max(vind, oind) + 1;
               std::vector<Vertex*> poly2 = {P.begin() + i_min, P.begin() + i_max};
               P.erase(P.begin() + i_min, P.begin() + i_max);
               P.push_back((*itr_r)->o); P.push_back((*itr_r)->v); P = SortPolygon(P);
               poly2.push_back((*itr_r)->o); poly2.push_back((*itr_r)->v); poly2 = SortPolygon(poly2);
               Polygons.push_back(poly2);
           }
       }
   }

    for(auto&poly: Polygons){
        std::vector<glm::f64vec3> V;
        poly = SortPolygon(poly);

        for(auto&p: poly)V.push_back(Mirror * glm::f64vec4{p->p3, 1});
        Vertices.push_back(V);
        planerity_value.push_back(Planerity(V,ui->glWid2dim->model->outline->Lines));
    }

    for(const auto& mesh: Vertices){
        for(const auto&v: mesh)WriteList.append("v " + QString::number(v.x) + " " + QString::number(v.y) + " " + QString::number(v.z) + "\n");
        glm::f64vec3 v =  glm::f64vec4{mesh[0], 1};
        glm::f64vec3 v2 =  glm::f64vec4{mesh[1], 1};
        glm::f64vec3 v3 =   glm::f64vec4{mesh[2], 1};
        N = -glm::normalize(glm::cross(v3 - v, v2 - v));
        Normals.push_back(N);
        befN = N;
    }
    for(const auto&n : Normals) WriteList.append("vn " + QString::number(n.x) + " " + QString::number(n.y) + " " + QString::number(n.z) + "\n");

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
    if(!(ui->glWid2dim->model->FL.empty() || ui->glWid2dim->model->FL[0]->FoldingCurve.empty())){
        for(auto itr = ui->glWid2dim->model->FL[0]->FoldingCurve.begin() + 1; itr != ui->glWid2dim->model->FL[0]->FoldingCurve.end() - 1; itr++){
            glm::f64vec3 et = itr->second->p3 - itr->first->p3, er = (itr - 1)->first->p3 - itr->first->p3, eb = itr->third->p3 - itr->first->p3, el = (itr + 1)->first->p3 - itr->first->p3;
            et /= glm::length(et); er /= glm::length(er); eb /= glm::length(eb); el /= glm::length(el);
            double phi1 = std::acos(glm::dot(et, er)), phi2 = std::acos(glm::dot(et, el)), phi3 = std::acos(glm::dot(eb, el)), phi4 = std::acos(glm::dot(eb, er));
            double a = abs(2.0*std::numbers::pi - (phi1 + phi2 + phi3 + phi4));
            if(a != -1)QuantitativeResult << a << ", ";
        }
    }


}

void MainWindow::addCurveBtn(){
    QRect geo = ui->LayerListWidget->geometry();
    int pad = 5;
    int btn_w = geo.width() - 2 * pad, btn_h = 25;
    QString text;
    if(ui->glWid2dim->model->crvs[0]->getCurveType() == CurveType::bezier3)text = "Bezier" + QString::number(++CurvesNum[0]);
    else if(ui->glWid2dim->model->crvs[0]->getCurveType()  == CurveType::bsp3)text = "Bspline" + QString::number(++CurvesNum[1]);
    else if(ui->glWid2dim->model->crvs[0]->getCurveType() == CurveType::line)text = "Line" + QString::number(++CurvesNum[2]);
    else if(ui->glWid2dim->model->crvs[0]->getCurveType() == CurveType::arc)text = "Arc"+QString::number(++CurvesNum[3]);

    Btn4Crv *newbtn = new Btn4Crv(ui->glWid2dim->model->crvs[0],text, ui->LayerListWidget);
    connect(newbtn, &Btn4Crv::clicked, this, &MainWindow::SetHandleCrv);
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
    SelectedBtn = btn;
    int ind;
    if(btn == nullptr){
        ind = -1;
        int pad = 5, btn_h = LayerList[0]->height();
        for(int i = 0; i < (int)LayerList.size(); i++) LayerList[i]->move(pad, btn_h + i * (btn_h + pad));
    }
    else{
        if(e->button() == Qt::LeftButton){
            dragPos = btn->mapFromGlobal(QCursor::pos());
            originalPos = btn->geometry();
            std::vector<Btn4Crv*>::iterator itr = std::find(LayerList.begin(), LayerList.end(), btn);
            if(itr == LayerList.end())ind = -1;
            else ind = std::distance(LayerList.begin(), itr);
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
