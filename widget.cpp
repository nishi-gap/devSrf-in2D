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
    CBoxlist = {{ui->SelectButton, PaintTool::None}, {ui->outline_rectangle, PaintTool::Rectangle_ol}, {ui->outline_polygon,  PaintTool::Polygon_ol},
                {ui->outline_polyline, PaintTool::Polyline_ol}, {ui->EditVertexButton, PaintTool::EditVertex_ol}, {ui->MoveOutLineButton, PaintTool::Move_ol},
                 {ui->Reset, PaintTool::Reset}
               };

    connect(ui->SelectButton, &QCheckBox::stateChanged, ui->glWid2dim, &GLWidget_2D::InitializeDrawMode);

    //色の最大値
    connect(ui->ColorLimitation, &QSpinBox::valueChanged, this, &MainWindow::ChangeMaxColor);

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
    connect(ui->addFL_line, &QPushButton::clicked, this, &MainWindow::addFoldLine_l);
    connect(ui->addFL_arc, &QPushButton::clicked, this, &MainWindow::addFoldLine_arc);
    connect(ui->addFL_bezier, &QPushButton::clicked, this, &MainWindow::addFoldLine_bezier);
    connect(this, &MainWindow::signalFLtype, ui->glWid2dim, &GLWidget_2D::changeFoldType);
    connect(ui->color_FL, &QPushButton::clicked, this, &MainWindow::color_FL);  
    connect(ui->move_ctrl_pt_fl, &QPushButton::clicked, this, &MainWindow::moveCtrlPts_fl);


    //fold line debug
    connect(ui->addFL_test, &QPushButton::clicked, this, &MainWindow::addFoldLine_test);    
    connect(ui->startButton, &QPushButton::clicked, ui->glWid2dim, &GLWidget_2D::Start4Debug_CF);
    connect(ui->OptBtn, &QPushButton::clicked, this, &MainWindow::StartOptimization);
    connect(ui->angleSlider, &QSlider::sliderMoved, this, &MainWindow::changeAngleFromSlider);
    connect(ui->angleA, &QDoubleSpinBox::valueChanged, this, &MainWindow::changeAngleFromSpinBox);
    connect(this, &MainWindow::sendAngle, ui->glWid2dim, &GLWidget_2D::changeBetaValue);

    connect(ui->DebugWindow, &QPushButton::clicked,ui->glWid2dim,&GLWidget_2D::OpenDebugWindwow);
    connect(ui->SaveButton, &QPushButton::clicked,this, &MainWindow::exportobj);

    LayerList.clear();
    connect(ui->glWid2dim, &GLWidget_2D::deleteCrvSignal, this, &MainWindow::RemoveBtnFromLayerCrv);
    connect(this, &MainWindow::signalNewSelectedCrv, ui->glWid2dim, &GLWidget_2D::changeSelectedCurve);
    connect(this, &MainWindow::swapIndex, ui->glWid2dim, &GLWidget_2D::swapCrvsOnLayer);

    //developability
    connect(ui->DevelopabilityButton, &QCheckBox::clicked, ui->glWid2dim, &GLWidget_2D::checkDevelopability);
    //planarity
    connect(ui->PlanarityButton, &QCheckBox::clicked, ui->glWid3dim, &GLWidget_3D::PlanarityDispay);


    CurvesNum[0] = MainWindow::CurvesNum[1] = MainWindow::CurvesNum[2] = MainWindow::CurvesNum[3] = 0;
    SelectedBtn = nullptr;
}

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

static int keyType = 0;

void MainWindow::StartOptimization(){
    if(ui->glWid2dim->model->FL.empty() || ui->glWid2dim->model->FL[0]->FoldingCurve.empty())return;
    auto Poly_V = ui->glWid2dim->model->outline->getVertices();
    bool res = ui->glWid2dim->model->FL[0]->Optimization(ui->glWid2dim->model->Edges, ui->glWid2dim->model->vertices, Poly_V, keyType);
    if(res){
        ui->glWid3dim->setVertices(ui->glWid2dim->model->Faces, ui->glWid2dim->model->outline->getVertices(), ui->glWid2dim->model->Edges, ui->glWid2dim->model->vertices, ui->glWid2dim->model->FL[0]->SingleRuling, ui->glWid2dim->AllRulings);
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

    //switchActivateCheckBox("MakeDevSrf");
    if(!ui->glWid2dim->model->outline->IsClosed())ui->glWid3dim->setVertices();
    else if(!ui->glWid2dim->model->FL.empty()){
        ui->glWid3dim->setVertices(ui->glWid2dim->model->Faces, ui->glWid2dim->model->outline->getVertices(), ui->glWid2dim->model->Edges, ui->glWid2dim->model->vertices, ui->glWid2dim->model->FL[0]->SingleRuling, ui->glWid2dim->AllRulings);
    }
    else{
        ui->glWid3dim->setVertices(ui->glWid2dim->model->Faces, ui->glWid2dim->model->outline->getVertices(), ui->glWid2dim->model->Edges, ui->glWid2dim->model->vertices);
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


void MainWindow::keyPressEvent(QKeyEvent *e){
    static bool switchDraw = false;
    switchDraw = (!switchDraw)? true: false;
    auto Poly_V = ui->glWid2dim->model->outline->getVertices();
    ui->glWid3dim->receiveKeyEvent(e);
    ui->glWid2dim->receiveKeyEvent(e);
    if(e->key() == Qt::Key_Backspace){
        emit PressedBackSpace();
    }
    else if(e->key() == Qt::Key_Return){emit PressedEnter();}
    else if(e->key() == Qt::Key_Q || e->key() == Qt::Key_P){
        if(!ui->glWid2dim->model->outline->IsClosed())ui->glWid3dim->setVertices();
        else if(!ui->glWid2dim->model->FL.empty()){
            ui->glWid3dim->setVertices(ui->glWid2dim->model->Faces, ui->glWid2dim->model->outline->getVertices(), ui->glWid2dim->model->Edges,
                                       ui->glWid2dim->model->vertices, ui->glWid2dim->model->FL[0]->SingleRuling, ui->glWid2dim->AllRulings, switchDraw);
        }
        else{
            ui->glWid3dim->setVertices(ui->glWid2dim->model->Faces, ui->glWid2dim->model->outline->getVertices(), ui->glWid2dim->model->Edges, ui->glWid2dim->model->vertices);
        }
    }
    else if(e->key() == Qt::Key_3 ||e->key() == Qt::Key_4){
        keyType = (e->key() == Qt::Key_3)? 0: 1;
        if(ui->glWid2dim->model->FL.empty() || ui->glWid2dim->model->FL[0]->FoldingCurve.empty())return;

        bool res = ui->glWid2dim->model->FL[0]->Optimization(ui->glWid2dim->model->Edges, ui->glWid2dim->model->vertices, Poly_V, keyType);
        if(res){
            ui->glWid3dim->setVertices(ui->glWid2dim->model->Faces, ui->glWid2dim->model->outline->getVertices(), ui->glWid2dim->model->Edges, ui->glWid2dim->model->vertices, ui->glWid2dim->model->FL[0]->SingleRuling, ui->glWid2dim->AllRulings);
        }
    }
    if(e->key() == Qt::Key_W){
        double a = (double)ui->angleSlider->value()/100.0;
        auto Edges = ui->glWid2dim->model->Edges;
        auto vertices = ui->glWid2dim->model->vertices;
        ui->glWid2dim->model->FL[0]->Optimization2(Edges, vertices, Poly_V, a);
        ui->glWid3dim->setVertices(ui->glWid2dim->model->Faces, ui->glWid2dim->model->outline->getVertices(), ui->glWid2dim->model->Edges,
                                   ui->glWid2dim->model->vertices, ui->glWid2dim->model->FL[0]->SingleRuling, ui->glWid2dim->model->FL[0]->AllRulings);
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
    if(ui->glWid2dim->model->Faces.empty()){
        QMessageBox msgBox;
        msgBox.setText("There is no developable surface.");
        msgBox.exec();
        return;
    }

    QStringList WriteList;
    std::vector<glm::f64vec3> Normals;
    glm::f64vec3 befN = {0,0,0}, N;
    glm::f64mat4x4 Mirror = glm::mat4(1.0f); Mirror[1][1] = -1;
   std::vector<Face*> outputFace;
   auto Poly_V = ui->glWid2dim->model->outline->getVertices();
   auto _edges = EdgeCopy(ui->glWid2dim->model->Edges, ui->glWid2dim->model->vertices);
   EdgeRecconection(Poly_V, outputFace, _edges);

    for(auto& f: outputFace){
        HalfEdge *h = f->halfedge;
        do{
            glm::f64vec3 v = Mirror * glm::f64vec4{h->vertex->p3, 1};
            WriteList.append("v " + QString::number(v.x) + " " + QString::number(v.y) + " " + QString::number(v.z) + "\n");
            h = h->next;
        }while(h != f->halfedge);
        glm::f64vec3 v = Mirror * glm::f64vec4{h->vertex->p3, 1};
        glm::f64vec3 v2 = Mirror * glm::f64vec4{h->next->vertex->p3, 1};
        glm::f64vec3 v3 = Mirror * glm::f64vec4{h->next->next->vertex->p3, 1};
        N = -glm::normalize(glm::cross(v3 - v, v2 - v));
        //if(glm::dot(befN, N) < 0) N *= -1;
        Normals.push_back(N);
        befN = N;
    }
    for(auto&n : Normals) WriteList.append("vn " + QString::number(n.x) + " " + QString::number(n.y) + " " + QString::number(n.z) + "\n");

    int cnt = 1;
    for(int i = 0; i < (int)outputFace.size(); i++){
        QString s = "f ";
        HalfEdge *h = outputFace[i]->halfedge;
        do{
            s += QString::number(cnt) + "//" + QString::number(i+1) + " ";
            cnt++;
            h = h->next;
        }while(h != outputFace[i]->halfedge);
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
    for(auto&c: ui->glWid3dim->PlanarityColor)QuantitativeResult << c << ", ";
    QuantitativeResult << "\nDevelopability\n" ;
    for(auto&v: ui->glWid2dim->model->vertices){
        double a = v->developability();
        if(a != -1)QuantitativeResult << a << ", ";
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
