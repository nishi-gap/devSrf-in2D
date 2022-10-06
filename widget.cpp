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
    CBoxlist = {{ui->SelectButton, "None"}, {ui->outline_rectangle, "Rectangle"}, {ui->outline_polygon,  "Polygon"},
                {ui->outline_polyline, "Polyline"}, {ui->EditVertexButton, "EditVertex"}, {ui->MoveOutLineButton, "MoveOutline"},
                {ui->make_devsrf, "MakeDevSrf"}, {ui->Reset, "Reset"}
               };

    connect(ui->SelectButton, &QCheckBox::stateChanged, ui->glWid2dim, &GLWidget_2D::InitializeDrawMode);

    connect(ui->make_devsrf, &QCheckBox::stateChanged, this , &MainWindow::FoldingModel);
    connect(ui->Reset, &QCheckBox::stateChanged,ui->glWid2dim,&GLWidget_2D::Reset);

    connect(ui->DvidedSizeSlider, &QSlider::valueChanged, this, &MainWindow::ChangeDivSizeEditFromSlider);
    connect(ui->DivSizeSpinBox, &QSpinBox::valueChanged, this, &MainWindow::ChangeDivSizeEditFromSpinBox);
    connect(ui->glWid2dim, &GLWidget_2D::foldingSignals, this , &MainWindow::FoldingModel);

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
    connect(ui->glWid2dim, &GLWidget_2D::signalCurveType, ui->CurveTypeBox, &QComboBox::currentText);
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

    connect(ui->DebugWindow, &QPushButton::clicked,ui->glWid2dim,&GLWidget_2D::OpenDebugWindwow);
    connect(ui->SaveButton, &QPushButton::clicked,this, &MainWindow::exportobj);

    LayerList.clear();
    connect(ui->glWid2dim, &GLWidget_2D::deleteCrvSignal, this, &MainWindow::RemoveBtnFromLayerCrv);
    connect(this, &MainWindow::signalNewSelectedCrv, ui->glWid2dim, &GLWidget_2D::changeSelectedCurve);
    connect(this, &MainWindow::swapIndex, ui->glWid2dim, &GLWidget_2D::swapCrvsOnLayer);
    CurvesNum[0] = MainWindow::CurvesNum[1] = MainWindow::CurvesNum[2] = 0;
    SelectedBtn = nullptr;
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::FoldingModel(){

    switchActivateCheckBox("MakeDevSrf");
    std::vector<ruling*> Rulings;
    for(auto&f: model->Faces){
        for(auto& ruling: f->rulings){
         if(ruling->IsCrossed == -1) Rulings.push_back(ruling);
        }
    }

    std::vector<Vertex*>outline = model->outline->getVertices();
    if(outline.empty() || Rulings.empty())return;
    glm::f64vec3 center;
    Model m = Model();
    m.setOutline(outline);
    m.deform(output, Rulings, center);
    ui->glWid3dim->setVertices(output, center);
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

void MainWindow::switchActivateCheckBox(QString active){
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
    if(output.size() == 0){
        QMessageBox msgBox;
        msgBox.setText("There is no developable surface.");
        msgBox.exec();
        return;
    }

    QStringList WriteList;
    std::vector<std::vector<glm::f64vec3>> TriMeshs;
    std::vector<std::vector<glm::f64vec3>> trimesh;
    std::vector<glm::f64vec3> Normals;
    glm::f64vec3 befN = {0,0,0};

    for(auto& mesh: output){
        Triangulation(mesh, trimesh);
        glm::f64vec3 N = glm::normalize(glm::cross(mesh[2] - mesh[0], mesh[1] - mesh[0]));
        if(befN == glm::f64vec3{0,0,0}) Normals.push_back(N);
        else{
            if(glm::dot(befN, N) < 0) N *= -1;
        }
        befN = N;
        for(auto& m: trimesh){
            TriMeshs.push_back(m);
            Normals.push_back(N);
        }
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

