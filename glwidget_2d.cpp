#define PI 3.14159265359
#include "glwidget_2d.h"

GLWidget_2D::GLWidget_2D(QWidget *parent):QOpenGLWidget(parent)
{

    CurveList = {{"Bezier", 0, "RulingBezier"}, {"B-spline", 1, "RulingBspline"}, {"Line", 2, "RulingLine"}};

    curveDimention = 3;
    drawtype = "OutlineRectangle";
    DivSize = 30;
    crvPtNum = 300;
    maxDivSize = 100;
    standardDist = 2.5;
    gw = new GToolWnd(this);
    gw->hide();

    cval = 128;
    ctype = 1;
    connect(this , &GLWidget_2D::CurvePathSet, gw, &GToolWnd::set);

    InterpolationType = 0;
    //ControllPoints_gradation.clear();
    DiffWheel = 0;
    refMeshNum = -1;
    movePt = -1;
    SelectedCurveIndex = -1;
    KeyEvent = -1;
    curvetype = 0;

}
GLWidget_2D::~GLWidget_2D(){}

void GLWidget_2D::InitializeDrawMode(int state){ if(state == 0)return; drawtype = "None"; SelectedCurveIndex = -1; emit SendNewActiveCheckBox("None");}

void GLWidget_2D::AddCurve(){
    std::vector<int>deleteIndex;
    QString type = emit signalCurveType();
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    for(auto& c: CurveList){
        if(std::get<0>(c) == type){
            drawtype = std::get<2>(c);
            curvetype = std::get<1>(c);
        }
    }

    SelectedCurveIndex = model->AddNewCurve(curvetype, DivSize);
    setMouseTracking(true);
    emit SendNewActiveCheckBox("AddCurve");
}

void GLWidget_2D::InsertNewPoint(){
    drawtype = "InsertControlPoint";
    SelectedCurveIndex = -1;
    setMouseTracking(true);
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    emit SendNewActiveCheckBox("InsertCtrlPt");
}

void GLWidget_2D::MoveCurvePt(){
    //curveNum = 1;
    drawtype = "MoveControlPoint";
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    SelectedCurveIndex = -1;
    setMouseTracking(true);
    emit SendNewActiveCheckBox("MoveControlPoint");
}

void GLWidget_2D::DrawOutlineRectangle(){
    //drawtype = -1;
    drawtype = "OutlineRectangle";
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    model->outline = new OUTLINE();
    model->outline->type = "Rectangle";
    this->setMouseTracking(false);
    emit SendNewActiveCheckBox("Rectangle");
}

void GLWidget_2D::DrawOutlinePolygon(int state){
    if(state == 0)return;
    drawtype = "OutlinePolygon";
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    model->outline = new OUTLINE();
    model->outline->type = "Polygon";
    for(auto&c: model->crvs)c->isempty = true;
    emit getEdgeNum();
    emit SendNewActiveCheckBox("Polygon");
    this->setMouseTracking(true);
}
void GLWidget_2D::DrawOutlinePolyline(int state){
    if(state == 0)return;
    drawtype = "OutlinePolyline";
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    model->outline = new OUTLINE();
    model->outline->type = "Polyline";
    for(auto&c: model->crvs)c->isempty = true;
    emit SendNewActiveCheckBox("Polyline");
}

void GLWidget_2D::MoveOutline(int state){
    if(state == 0)return; drawtype = "MoveOutline";
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    emit SendNewActiveCheckBox("MoveOutline");
}

void GLWidget_2D::EditOutlineVertex(int state){
    if(state == 0)return; drawtype = "EditVertex";
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    emit SendNewActiveCheckBox("EditVertex");
}

void GLWidget_2D::DeleteCtrlPt(){
    drawtype = "DeleteCntrlPt";
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    emit deleteCrvSignal(deleteIndex);
    SelectedCurveIndex = -1;
    setMouseTracking(true);
    emit SendNewActiveCheckBox("DeleteCntrlPt");
}

void GLWidget_2D::DeleteCurve(){
    drawtype = "DeleteCurve";
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    setMouseTracking(true);
    emit SendNewActiveCheckBox("DeleteCurve");

}

void GLWidget_2D::setColor(){
    drawtype = "SetColor";
}

void GLWidget_2D::setNewGradationMode(){
    std::vector<int> deleteIndex;
    drawtype = "NewGradationMode";
    model->Check4Param(curveDimention, deleteIndex); SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
}
void GLWidget_2D::recieveNewEdgeNum(int num){model->outline->VerticesNum = num; update();}

void GLWidget_2D::ChangedDivSizeEdit(int n){
    this->DivSize = (n < 0)? DivSize: (maxDivSize < n)? maxDivSize: n;
    if(SelectedCurveIndex == -1) return;
    if(curvetype == 0){
        model->crvs[SelectedCurveIndex]->Bezier(curveDimention, crvPtNum);
        //if(model->outline->isClosed  && model->crv->CurvePoints[0].pt != glm::f64vec2{-1,-1})model->crv->BezierRulings(model->outline->vertices,DivSize,crvPtNum);
    }
    if(curvetype == 1){
        model->crvs[SelectedCurveIndex]->Bspline(curveDimention, crvPtNum);
        if(model->outline->isClosed && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1, -1})model->crvs[SelectedCurveIndex]->BsplineRulings(model->outline,DivSize,crvPtNum, curveDimention);
    }
    if(curvetype == 2){
        model->crvs[SelectedCurveIndex]->Line();
        if(model->outline->isClosed && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1})model->crvs[SelectedCurveIndex]->LineRulings(model->outline,DivSize);
    }

    if(model->outline->isClosed && model->CrossDection4AllCurve()){
        model->updateRulings();
        emit foldingSignals();
    }

    update();
}

void GLWidget_2D::changeSelectedCurve(int ind){SelectedCurveIndex = ind; drawtype = "None"; update();}

void GLWidget_2D::swapCrvsOnLayer(int n1, int n2){
    if(n1 == n2)return;
    std::iter_swap(model->crvs.begin() + n1, model->crvs.begin() + n2);
    if(SelectedCurveIndex == n1)SelectedCurveIndex = n2;
    else SelectedCurveIndex = n1;
    model->updateRulings();
    emit foldingSignals();
    update();
}

void GLWidget_2D::initializeGL(){
    makeCurrent();
    initializeOpenGLFunctions();
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    QSize s = this->size();
    glViewport(0,0,s.width(),s.height());
}

void GLWidget_2D::paintGL(){
    makeCurrent();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    QSize s = this->size();
    glLoadIdentity();
    glOrtho(-0.5, (float)s.width() -0.5, (float)s.height() -0.5, -0.5, -1, 1);

    //グラデーションの描画
    float r, g, b = 0;
    std::vector<glm::f64vec3> Mesh;
    std::vector<std::vector<glm::f64vec3>> TriMeshs;

    for(auto& f : model->Faces){
        Mesh.clear();
        HalfEdge *h = f->halfedge;
        do{
            Mesh.push_back(h->vertex->p);
            h = h->next;
        }while(h != f->halfedge);
        Triangulation(Mesh, TriMeshs);
        if(f->Gradation == 0){
            r = g = b = 1;
        } else if(f->Gradation > 0){
            r = 1; g = b = 1 - (float)f->Gradation/255.f;
        }else{
            b = 1;
            g = r = 1 + (float)f->Gradation/255.f;
        }
        glColor3d(r,g,b);
        for(auto& tri: TriMeshs){
            glBegin(GL_POLYGON);
            for(auto&v: tri){
                glVertex2d(v.x, v.y);
            }
            glEnd();
        }
    }

    //rulingの描画
    glm::f64vec2 p, p2;
    glPolygonOffset(0.0,1.0);
    for(auto& curve: model->crvs){
        if(curve->isempty)continue;
        for (int i = 0; i < (int)curve->Rulings.size(); i++) {
            if(curve->Rulings[i]->IsCrossed == -1)glColor3d(0,0,0);
            else if(curve->Rulings[i]->IsCrossed == 0) glColor3d(0,1,0);
            else continue;
            glBegin(GL_LINES);
            glLineWidth(5.0f);
            p = std::get<0>(curve->Rulings[i]->r)->p, p2 = std::get<1>(curve->Rulings[i]->r)->p;
            glVertex2d(p.x, p.y);
            glVertex2d(p2.x, p2.y);
            glEnd();
        }
    }


    //曲線の制御点
    for(int i = 0; i < (int)model->crvs.size(); i++){
        if((i == model->getSelectedCurveIndex(mapFromGlobal(QCursor::pos()))) || ((drawtype == "RulingBezier" || drawtype == "RulingBspline" || drawtype == "RulingLine")
                 && SelectedCurveIndex == i) || drawtype == "MoveControlPoint" || drawtype == "InsertControlPoint" ||  drawtype == "DeleteCntrlPt" || drawtype == "None"){
            glPolygonOffset(0.f,1.f);
            for (auto& c: model->crvs[i]->ControllPoints) {
                glColor3d(1, 0, 0);
                glPointSize(5.0f);
                glBegin(GL_POINTS);
                glVertex2d( c.x, c.y);
                glEnd();
            }
            glEnable(GL_LINE_STIPPLE);
            glLineStipple(1 , 0xF0F0);
            glBegin(GL_LINE_STRIP);
            glColor3d(0.4, 0.4, 0.4);
            for (auto& c: model->crvs[i]->ControllPoints)glVertex2d( c.x, c.y);
            glEnd();
            glDisable(GL_LINE_STIPPLE);
            glPolygonOffset(0.f,0.f);
        }
    }

    //可展面の輪郭描画
    if(model->outline->type == "Rectangle" || model->outline->type == "Polyline"){
        glColor3d(0, 0, 0);
        glPointSize(5.0f);
        for(auto& v: model->outline->getVertices()){
            glBegin(GL_POINTS);
            glVertex2d(v->p.x, v->p.y);
            glEnd();
        }

        if(model->outline->isClosed){
            glBegin(GL_LINE_LOOP);
            for(auto& r: model->outline->getVertices())glVertex2d(r->p.x, r->p.y);
            glEnd();
        }else{
            glBegin(GL_LINE_STRIP);
            for(auto& v: model->outline->getVertices()) glVertex2d(v->p.x, v->p.y);
            glEnd();
        }
    }
    if(model->outline->type == "Polygon"){
        if(model->outline->hasPtNum > 0){
            glColor3d(0, 0, 0);
            glPointSize(5.0f);
            glBegin(GL_POINTS);
            glVertex2d(model->outline->origin.x, model->outline->origin.y);
            glEnd();

            glColor3d(0, 0, 0);
            glPointSize(5.0f);
            for(auto& v: model->outline->getVertices()){
                glBegin(GL_POINTS);
                glVertex2d(v->p.x, v->p.y);
                glEnd();
            }

            glBegin(GL_LINE_LOOP);
            for(auto& r: model->outline->getVertices())glVertex2d(r->p.x, r->p.y);
            glEnd();
        }

    }

    //曲線の描画
    for(int i = 0; i < (int)model->crvs.size(); i++){
        if(i == SelectedCurveIndex)glColor3d(1, 0, 0);
        else glColor3d(0.4, 0.4, 0.4);
        glBegin(GL_LINE_STRIP);
        for(auto &c: model->crvs[i]->CurvePoints){
            glVertex2s(c.pt.x, c.pt.y);
        }
        glEnd();
    }

    //挿入する制御点
    if(drawtype == "InsertControlPoint" && SelectedCurveIndex != -1 && model->crvs[SelectedCurveIndex]->InsertPointSegment != -1){
        glColor3d(0,0,1);
        glPointSize(5.0f);
        glBegin(GL_POINTS);
        glVertex2d(model->crvs[SelectedCurveIndex]->InsertPoint.x, model->crvs[SelectedCurveIndex]->InsertPoint.y);
        glEnd();
    }
}

void GLWidget_2D::mousePressEvent(QMouseEvent *e){
    QPointF p = this->mapFromGlobal(QCursor::pos());
    if(drawtype == "None") {return;}
    if(e->button() ==Qt::LeftButton){
        if(drawtype == "OutlineRectangle"  || drawtype == "OutlinePolyline"){
            std::vector<Vertex*> vertices = model->outline->getVertices();
            if(model->outline->hasPtNum != 2)model->outline->addVertex(p);
            for(auto&c: model->crvs){
                if(c->getCurveType() == 0 )c->Bezier(curveDimention, crvPtNum);
                else if(c->getCurveType() == 1)c->Bspline(curveDimention,crvPtNum);
                else if(c->getCurveType() == 2)c->Line();
            }
        }
        else if(drawtype == "OutlinePolygon"){
            model->outline->drawPolygon(p, true);
            for(auto&c: model->crvs){
                if(c->getCurveType() == 0 )c->Bezier(curveDimention, crvPtNum);
                else if(c->getCurveType() == 1)c->Bspline(curveDimention,crvPtNum);
                else if(c->getCurveType() == 2)c->Line();
            }
        }
        else if(drawtype == "EditVertex") model->outline->EditVertex(p);
        else if(drawtype == "MoveOutline") model->outline->MoveVertex(p);
        else if(drawtype == "NewGradationMode")addPoints_intplation(e, p);
        else if(drawtype == "DeleteCurve"){
            model->SelectCurve(p);
            std::vector<int> n(1);
            n[0] = model->DeleteCurve();
            deleteCrvSignal(n);
            SelectedCurveIndex = -1;
        }
        if(drawtype == "RulingBezier" || drawtype == "RulingBspline" || drawtype == "RulingLine"){model->AddControlPoint(p, curveDimention, DivSize); }
        else if(drawtype == "InsertControlPoint"){
            if(SelectedCurveIndex == -1){
                model->SelectCurve(p);
                SelectedCurveIndex = model->IsSelectedCurve();
            }else{
                if(model->crvs[SelectedCurveIndex]->getCurveType() == 1 && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){
                    model->crvs[SelectedCurveIndex]->InsertControlPoint2(p);
                    model->crvs[SelectedCurveIndex]->SetNewPoint();
                    model->crvs[SelectedCurveIndex]->Bspline(curveDimention,crvPtNum);
                    if(model->outline->isClosed){
                        model->crvs[SelectedCurveIndex]->BsplineRulings(model->outline,DivSize,crvPtNum, curveDimention);
                        if(!model->crvs[SelectedCurveIndex]->isempty){
                            model->updateRulings();
                            emit foldingSignals();
                        }
                    }
                }
            }
        }
        else if(drawtype == "DeleteCntrlPt"){
            model->SelectCurve(p);
            model->DeleteControlPoint(p, curveDimention, DivSize);
        }else if(drawtype == "MoveControlPoint"){
            SelectedCurveIndex = model->searchPointIndex(p, movePt, 0);

        }

    }else if(e->button() == Qt::RightButton){
        if(drawtype == "OutlineRectangle" || drawtype == "OutlinePolyline") model->outline->eraseVertex();
        if(drawtype == "OutlinePolygon") model->outline->eraseVertex();
        if(drawtype == "RulingBezier" || drawtype == "RulingBspline" || drawtype == "RulingLine"){
            if(SelectedCurveIndex != -1) model->crvs[SelectedCurveIndex]->eraseCtrlPt(curveDimention, crvPtNum);
        }

        if(model->outline->isClosed)model->addRulings();
    }
    if(model->outline->isClosed) emit foldingSignals();
    update();
}

void GLWidget_2D::mouseMoveEvent(QMouseEvent *e){
    if(drawtype == "None") {return;}
    QPointF p = this->mapFromGlobal(QCursor::pos());
    if(drawtype == "OutlinePolygon"){
        if(model->outline->hasPtNum == 1){
            model->outline->drawPolygon(p, false);
            model->outline->isClosed = false;
            model->outline->hasPtNum = 1;
        }
    }
    else if(drawtype == "OutlineRectangle" && movePt != -1)model->outline->getVertices()[movePt]->p = glm::f64vec3{p.x(), p.y(),0};
    else if(drawtype == "EditVertex") model->outline->EditVertex(p);
    else if(drawtype == "MoveOutline") model->outline->MoveVertex(p);
    //else if(model->getSelectedCurveIndex(p) == -1 && (drawtype == "DeleteCurve" || drawtype == "MoveControlPoint" || drawtype == "InsertControlPoint" || drawtype == "DeleteCntrlPt"))
        //SelectedCurveIndex = model->getSelectedCurveIndex(p);

    if(drawtype == "MoveControlPoint")model->MoveCurvePoint(p,SelectedCurveIndex, movePt, curveDimention, DivSize);

    if(SelectedCurveIndex != -1){
        if(drawtype == "InsertControlPoint" &&  model->crvs[SelectedCurveIndex]->getCurveType() == 1 && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){//制御点の挿入(B-spline)
            model->crvs[SelectedCurveIndex]->InsertControlPoint2(p);
        }
        
        if(model->outline->isClosed && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){
            //if(model->crvs[SelectedCurveIndex]->getCurveType() == 0)model->crvs[SelectedCurveIndex]->BezierRulings(model->outline->vertices,DivSize,crvPtNum);
            if(model->crvs[SelectedCurveIndex]->getCurveType() == 1) model->crvs[SelectedCurveIndex]->BsplineRulings(model->outline,DivSize,crvPtNum, curveDimention);
            if(model->crvs[SelectedCurveIndex]->getCurveType() == 2){
                model->crvs[SelectedCurveIndex]->LineRulings(model->outline,DivSize);
            }
            model->updateRulings();
        }
    }

    if(model->CrossDection4AllCurve()){
        //model->addRulings();
        emit foldingSignals();
    }
    update();

}

void GLWidget_2D::mouseReleaseEvent(QMouseEvent * e){
    movePt = -1;
    if(SelectedCurveIndex != -1 && !model->crvs.empty()) model->crvs[SelectedCurveIndex]->InsertPointSegment = -1;
}

void GLWidget_2D::cb_ApplyCurveEvent(){
    if(KeyEvent == 0){
        std::cout << "finished " << std::endl; KeyEvent = -1;
    }else if(KeyEvent == -1){
        if(SelectedCurveIndex == -1) std::cout << "no curve is selected"<<std::endl;
        else {
            std::cout << "drawtype " << drawtype.toStdString() << std::endl;
        }
    }
    update();
}

void GLWidget_2D::cb_DeleteCurve(){
    KeyEvent = 1;
    if(SelectedCurveIndex != -1){
        model->DeleteCurve();
    }
    update();
}

void GLWidget_2D::Reset(){
    model->Initialize();
    emit SendNewActiveCheckBox("Reset");
    SelectedCurveIndex = -1;
    update();
}

void GLWidget_2D::receiveColors(){
    emit foldingSignals();
    update();
}

int GLWidget_2D::referencedRuling(QPointF p){
    //if(refcurve == -1){qDebug()<<"there is no curve"; return -1;}
    glm::f64vec2 a,b;
    int n = -1;
    glm::f64vec2 v = {p.x(), p.y()};
    double minDist = standardDist;
    //auto Rulings = CRVS[refcurve].getRulings();
    for(auto&c: model->crvs){
        for(int i = 0; i < (int)c->Rulings.size(); i++){
            a = std::get<0>(c->Rulings[i]->r)->p, b = std::get<1>(c->Rulings[i]->r)->p;
            double d = abs((v.x - a.x) * (a.y - b.y) - (v.y - a.y) * (a.x - b.x)) / sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
            if(d < minDist){
                n = i;
                minDist = d;
            }
        }
    }

    return n;
}

void GLWidget_2D::wheelEvent(QWheelEvent *we){
    if(drawtype != "NewGradationMode")return;

    DiffWheel = (we->angleDelta().y() > 0) ? 1 : -1;
    int color;
    model->setGradationValue(DiffWheel, refMeshNum, color, InterpolationType, CurvePath);
    emit ColorChangeFrom(0, color);
    emit foldingSignals();

    if(isVisibleTo(gw)) emit CurvePathSet(CurvePath);
    update();
}

void GLWidget_2D::addPoints_intplation(QMouseEvent *e, QPointF& p){
    int i = 0;
    std::vector<glm::f64vec2> mesh;
    HalfEdge *h;
    for(auto& f: model->Faces){
        h = f->halfedge;
        mesh.clear();
        do{
            mesh.push_back(h->vertex->p);
            h = h->next;
        }while(h != f->halfedge);
        if(cn(mesh, p)){
            refMeshNum = i;

            if(e->button() == Qt::RightButton)model->Faces[i]->hasGradPt = false;
            else{

                f->hasGradPt = true;
                f->Gradation = f->rulings[0]->pt->color;
                f->Gradation = (f->Gradation < -255)? -255 : (255 < f->Gradation)? 255 : f->Gradation;
                f->rulings[0]->pt->color = f->Gradation;
                int color;
                model->setGradationValue(0, refMeshNum, color, InterpolationType, CurvePath);
                emit ColorChangeFrom(0, f->Gradation);
                emit foldingSignals();
            }
        }
        i++;
    }


    update();
}

void GLWidget_2D::ApplyNewGradationMode(){
    //if(refcurve == -1){qDebug() << "there is no curve"; return;}
    //auto CurvePoints = CRVS[refcurve].getCurvePoints();
    if(refMeshNum == -1){std::cout << "you neeed to add gradation point"<<std::endl; return;}
    //QPointF p{-1,-1};
    //int color;
    //model->setGradationValue(p, 0, refMeshNum, color, InterpolationType, CurvePath);
    emit foldingSignals();
    if(isVisibleTo(gw)) emit CurvePathSet(CurvePath);
    update();
}

void GLWidget_2D::getGradationFromSlider(int val){
    if(refMeshNum == -1) return;
    model->Faces[refMeshNum]->Gradation = val;
    for(auto& r: model->Faces[refMeshNum]->rulings)r->pt->color = val;
    //model->crv->setcolor(1, val, refMeshNum);
    emit ColorChangeFrom(2, val);
}

void GLWidget_2D::OpenDebugWindwow(){
    if(!isVisibleTo(gw)) gw = new GToolWnd(this);
    emit CurvePathSet(CurvePath);

    gw->show();
}
