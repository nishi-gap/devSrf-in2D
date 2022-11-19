#include "widget.h"
#include "ui_widget.h"

MainWindow::MainWindow(QWidget *parent)
    : QWidget(parent), ui(new Ui::MainWindow)
{

    crvPtNum = 300;
    output.clear();

    ui->setupUi(this);

    model = new Model(crvPtNum);
    ui->glWid2dim->model = model;
    CBoxlist = {{ui->SelectButton, PaintTool::None}, {ui->outline_rectangle, PaintTool::Rectangle_ol}, {ui->outline_polygon,  PaintTool::Polygon_ol},
                {ui->outline_polyline, PaintTool::Polyline_ol}, {ui->EditVertexButton, PaintTool::EditVertex_ol}, {ui->MoveOutLineButton, PaintTool::Move_ol},
                {ui->make_devsrf, PaintTool::deform}, {ui->Reset, PaintTool::Reset}
               };

    connect(ui->SelectButton, &QCheckBox::stateChanged, ui->glWid2dim, &GLWidget_2D::InitializeDrawMode);

    connect(ui->make_devsrf, &QCheckBox::stateChanged, this , &MainWindow::fold_Sm);
    connect(ui->Reset, &QCheckBox::stateChanged,ui->glWid2dim,&GLWidget_2D::Reset);

    connect(ui->DvidedSizeSlider, &QSlider::valueChanged, this, &MainWindow::ChangeDivSizeEditFromSlider);
    connect(ui->DivSizeSpinBox, &QSpinBox::valueChanged, this, &MainWindow::ChangeDivSizeEditFromSpinBox);
    connect(ui->glWid2dim, &GLWidget_2D::foldingSignals, this , &MainWindow::fold_Sm);

    //line width
    connect(ui->LineWidthSpinBox, &QDoubleSpinBox::valueChanged, this, &MainWindow::changeLineWidthFromSpinBox);
    connect(ui->LineWidthSlider, &QSlider::valueChanged, this, &MainWindow::changeLineWidthFromSlider);
    connect(this, &MainWindow::signalNewLineWidth, ui->glWid2dim, &GLWidget_2D::receiveNewLineWidth);

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
    connect(ui->glWid2dim, &GLWidget_2D::signalCurveType, ui->CurveTypeBox, &QComboBox::currentIndex);
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
    connect(ui->addFL_bsp, &QPushButton::clicked, this, &MainWindow::addFoldLine_bsp);
    connect(ui->addFL_bezier, &QPushButton::clicked, this, &MainWindow::addFoldLine_bezier);
    connect(this, &MainWindow::signalFLtype, ui->glWid2dim, &GLWidget_2D::changeFoldType);
    connect(ui->glWid2dim, &GLWidget_2D::signalAddRulings_FL, this, &MainWindow::fold_FL);
    connect(ui->color_FL, &QPushButton::clicked, this, &MainWindow::color_FL);

    connect(ui->DebugWindow, &QPushButton::clicked,ui->glWid2dim,&GLWidget_2D::OpenDebugWindwow);
    connect(ui->SaveButton, &QPushButton::clicked,this, &MainWindow::exportobj);

    LayerList.clear();
    connect(ui->glWid2dim, &GLWidget_2D::deleteCrvSignal, this, &MainWindow::RemoveBtnFromLayerCrv);
    connect(this, &MainWindow::signalNewSelectedCrv, ui->glWid2dim, &GLWidget_2D::changeSelectedCurve);
    connect(this, &MainWindow::swapIndex, ui->glWid2dim, &GLWidget_2D::swapCrvsOnLayer);
    CurvesNum[0] = MainWindow::CurvesNum[1] = MainWindow::CurvesNum[2] = MainWindow::CurvesNum[3] = 0;
    SelectedBtn = nullptr;
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::SymmetricConstraint(){ emit constraintType(0);}

void MainWindow::addFoldLine_l(){
    emit signalFLtype(0);
}
void MainWindow::addFoldLine_bsp(){
    emit signalFLtype(1);
}
void MainWindow::addFoldLine_bezier(){
    emit signalFLtype(2);
}
void MainWindow::color_FL(){emit signalFLtype(3);}

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
    if(!model->outline->IsClosed())return;

    ui->glWid3dim->setVertices(model->Faces);
}

void MainWindow::fold_FL(){
    glm::f64mat4x4 Scale = glm::scale(glm::f64vec3{0.05, 0.05, 0.05});
    glm::f64mat4x4 Mirror = glm::mat4(1.0f); Mirror[1][1] = -1;
    std::vector<std::array<glm::f64vec3, 2>> Ruling_l = model->FL[0]->Rulings_3dL;
    std::vector<std::array<glm::f64vec3, 2>> Ruling_r = model->FL[0]->Rulings_3dR;
    output.clear();
    for(auto& m : Ruling_l){
        for(auto&v: m)
            v = Mirror * Scale * glm::f64vec4{v,1};
    }
    for(auto& m : Ruling_r){
        for(auto&v: m)
            v = Mirror * Scale * glm::f64vec4{v,1};
    }
    std::vector<std::vector<glm::f64vec3>> Left, Right;
    for(int i = 1; i < (int)Ruling_l.size(); i++){
        std::vector<glm::f64vec3> tmp{Ruling_l[i - 1][0], Ruling_l[i][0], Ruling_l[i][1], Ruling_l[i - 1][1]};
        Left.push_back(tmp);
        output.push_back(tmp);
        tmp = {Ruling_r[i - 1][0], Ruling_r[i][0], Ruling_r[i][1], Ruling_r[i - 1][1]};
        output.push_back(tmp);
        Right.push_back(tmp);
    }

    glm::f64vec3 center{0,0,0};
    for(auto& m : Left){
        double s = m.size();
        glm::f64vec3 c{0,0,0};
        for(auto&v: m){
            c += v;
        }
        center += c /s;
    }
    for(auto& m : Right){
        double s = m.size();
        glm::f64vec3 c{0,0,0};
        for(auto&v: m){
            c += v;
        }
        center += c /s;
    }
    center /= (Left.size() + Right.size());
    ui->glWid3dim->receive(Left, Right, center);
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

    if(e->key() == Qt::Key_Backspace){
        emit PressedBackSpace();
    }
    else if(e->key() == Qt::Key_Return){emit PressedEnter();}
    else{

    }
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
    if(model->Faces.empty()){
        QMessageBox msgBox;
        msgBox.setText("There is no developable surface.");
        msgBox.exec();
        return;
    }

    QStringList WriteList;
    std::vector<std::array<glm::f64vec3, 3>> TriMeshs;
    std::vector<std::array<glm::f64vec3, 3>> trimesh;
    std::vector<glm::f64vec3> Normals;
    glm::f64vec3 befN = {0,0,0};

    for(auto& f: model->Faces){
        std::vector<glm::f64vec3> polygon;
        HalfEdge *h = f->halfedge;
        do{
            polygon.push_back(h->vertex->p3);
            h = h->next;
        }while(h != f->halfedge);
        Triangulation(polygon, trimesh);
        glm::f64vec3 N = glm::normalize(glm::cross(trimesh[0][2] - trimesh[0][0], trimesh[0][1] - trimesh[0][0]));
        if(befN == glm::f64vec3{0,0,0}) Normals.push_back(N);
        else{
            if(glm::dot(befN, N) < 0) N *= -1;
        }
        befN = N;
        std::vector<glm::f64vec3> _Normals((int)trimesh.size(), N);
        Normals.insert(Normals.end(), _Normals.begin(), _Normals.end());
        TriMeshs.insert(TriMeshs.end(), trimesh.begin(), trimesh.end());
    }

    for(auto& mesh: TriMeshs){
        for(auto&v: mesh){
            WriteList.append("v " + QString::number(v.x) + " " + QString::number(v.y) + " " + QString::number(v.z) + "\n");
        }
    }
    for(auto&n : Normals){
        WriteList.append("vn " + QString::number(n.x) + " " + QString::number(n.y) + " " + QString::number(n.z) + "\n");
    }
    int cnt = 1;
    for(int i = 1; i <= (int)TriMeshs.size(); i++){
        QString s = "f ";
        for(int j = 0; j < (int)TriMeshs[i - 1].size(); j++){
            s += QString::number(cnt) + "//" + QString::number(i) + " ";
            cnt++;
        }
        s += "\n";
        WriteList.append(s);
    }
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


}


void MainWindow::addCurveBtn(){

    QRect geo = ui->LayerListWidget->geometry();
    int pad = 5;
    int btn_w = geo.width() - 2 * pad, btn_h = 25;
    QString text;
    if(model->crvs[0]->getCurveType() == 0)text = "Bezier" + QString::number(++CurvesNum[0]);
    else if(model->crvs[0]->getCurveType()  == 1)text = "Bspline" + QString::number(++CurvesNum[1]);
    else if(model->crvs[0]->getCurveType() == 2)text = "Line" + QString::number(++CurvesNum[2]);
    else if(model->crvs[0]->getCurveType() == 3)text = "Arc"+QString::number(++CurvesNum[3]);

    Btn4Crv *newbtn = new Btn4Crv(model->crvs[0],text, ui->LayerListWidget);
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

