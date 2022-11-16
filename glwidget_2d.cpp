#define PI 3.14159265359
#include "glwidget_2d.h"

GLWidget_2D::GLWidget_2D(QWidget *parent):QOpenGLWidget(parent)
{

    CurveList = {{"Bezier", 0, "RulingBezier"}, {"B-spline", 1, "RulingBspline"}, {"Line", 2, "RulingLine"}, {"Arc", 3, "RulingArc"}};

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
    movePt = -1;
    SelectedCurveIndex = -1;
    refHE = nullptr;
    KeyEvent = -1;
    curvetype = 0;

    gridsize = 10;

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
    drawtype = "EditVertex";
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

void GLWidget_2D::changeFoldType(int state){
    if(state < 3){
        drawtype = "FoldLine";
        int rulingNum = 30;
        int crvNum = 1000;
        FoldLine *fl = new FoldLine(crvNum, rulingNum, state);
        model->FL.push_back(fl);
    }else drawtype = "FoldlineColor";
    setMouseTracking(false);
}

void GLWidget_2D::setColor(){
    drawtype = "SetColor";
}

void GLWidget_2D::ConnectVertices(){drawtype = "ConnectVertices";}

void GLWidget_2D::switchGetmetricConstraint(int state){
    drawtype = "OutlineConst";
    constType = state;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
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
        if(model->outline->IsClosed() && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1, -1})model->crvs[SelectedCurveIndex]->BsplineRulings(model->outline,DivSize,crvPtNum, curveDimention);
    }
    if(curvetype == 2){
        model->crvs[SelectedCurveIndex]->Line();
        if(model->outline->IsClosed() && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1})model->crvs[SelectedCurveIndex]->LineRulings(model->outline,DivSize);
    }

    if(model->outline->IsClosed() && model->CrossDection4AllCurve()){
        model->addRulings();
        model->deform();
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
    model->addRulings();
    model->deform();
    emit foldingSignals();
    update();
}

void GLWidget_2D::initializeGL(){
    makeCurrent();
    initializeOpenGLFunctions();
    glClearColor(0.99f, 0.99f, 0.99f, 1.f);
    QSize s = this->size();
    glViewport(0,0,s.width(),s.height());
}

void GLWidget_2D::paintGL(){
    makeCurrent();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    QSize s = this->size();
    glLoadIdentity();
    glOrtho(-0.5, (float)s.width() -0.5, (float)s.height() -0.5, -0.5, -1, 1);

    DrawGrid();
    //rulingの描画
    {
        glm::f64vec2 p, p2;
        float r, g, b = 0;
        glPolygonOffset(0.0,1.0);
        HalfEdge *_refHE = assignment_refHE();
        for(auto& curve: model->crvs){
            if(curve->isempty)continue;
            for(auto&rl: curve->Rulings){
                if(_refHE != nullptr && _refHE->r == rl)glColor3d(1,1,0);
                else{
                    if(rl->IsCrossed == 1)continue;
                    if(rl->IsCrossed == 0)glColor3d(0,1,0);
                    else{
                        if(rl->Gradation == 0) r = g = b = 0.4;
                        else if(rl->Gradation > 0){
                            r = 1; g = b = 1 - (float)rl->Gradation/255.f;
                        }else{
                            b = 1;
                            g = r = 1 + (float)rl->Gradation/255.f;
                        }
                        glColor3d(r,g,b);
                    }
                }

                glBegin(GL_LINES);
                glLineWidth(5.0f);
                p = std::get<0>(rl->r)->p, p2 = std::get<1>(rl->r)->p;
                glVertex2d(p.x, p.y);
                glVertex2d(p2.x, p2.y);
                glEnd();
            }
        }
    }

    //曲線の制御点
    for(int i = 0; i < (int)model->crvs.size(); i++){
        if((i == model->getSelectedCurveIndex(mapFromGlobal(QCursor::pos()))) || ((drawtype == "RulingBezier" || drawtype == "RulingBspline" || drawtype == "RulingLine" || drawtype == "RulingArc")
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
    {
        if(model->outline->type == "Rectangle" || model->outline->type == "Polyline"){
            glColor3d(0, 0, 0);
            glPointSize(3.0f);
            for(auto& v: model->outline->getVertices()){
                glBegin(GL_POINTS);
                glVertex2d(v->p.x, v->p.y);
                glEnd();
            }

            if(model->outline->IsClosed()){
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
                glPointSize(4.0f);
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
        glPointSize(4);
        for(auto& Vertices: model->ol_vertices){
            for(auto&v: Vertices){
                glBegin(GL_POINTS);
                glVertex2d(v->p.x, v->p.y);
                glEnd();
            }
        }

        glColor3d(0, 0, 0);
        for(auto& Vertices: model->ol_vertices){
            glBegin(GL_LINE_STRIP);
            for(auto&v: Vertices){
                glVertex2d(v->p.x, v->p.y);
            }
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

    if(drawtype == "OutlineConst"){
        glPolygonOffset(0.f,0.5f);
        if(model->Axis4Const[0] != glm::f64vec3{-1,-1,0} && model->Axis4Const[1] != glm::f64vec3{-1,-1,0}){
            glColor3d(1,0,0);
            glLineWidth(2);
            glBegin(GL_LINES);
            glVertex2d(model->Axis4Const[0].x, model->Axis4Const[0].y);
            glVertex2d(model->Axis4Const[1].x, model->Axis4Const[1].y);
            glEnd();
        }else if(model->Axis4Const[0] != glm::f64vec3{-1,-1,0} && model->Axis4Const[1] == glm::f64vec3{-1,-1,0}){
            glColor3d(1,0,0);
            glPointSize(3.0f);
            glBegin(GL_POINTS);
            glVertex2d(model->Axis4Const[0].x, model->Axis4Const[0].y);
            glEnd();
        }
    }

    //折り線の描画

    for(auto&fl: model->FL){
        glColor3d(0,0,0);
        glBegin(GL_LINE_STRIP);
        for(auto&v: fl->CurvePts)glVertex2d(v.x, v.y);
        glEnd();
        std::vector<glm::f64vec3> Pts = fl->getCtrlPt();
        glColor3d(1,0,0);
        glPointSize(5);
        for(auto&v: Pts){
            glBegin(GL_POINTS);
            glVertex2d(v.x,v.y);
            glEnd();
        }
        glColor3d(0,0,0);
        for(auto&l2: fl->Rulings_2dL){
            glBegin(GL_LINES);
            glVertex2d(l2[0].x, l2[0].y);
            glVertex2d(l2[1].x, l2[1].y);
            glEnd();
        }
        for(auto&r2: fl->Rulings_2dR){
            glBegin(GL_LINES);
            glVertex2d(r2[0].x, r2[0].y);
            glVertex2d(r2[1].x, r2[1].y);
            glEnd();
        }

        if(fl->he2 != nullptr){
            glColor3d(0,0,0);
            glBegin(GL_LINES);
            glVertex2d(fl->he->vertex->p.x, fl->he->vertex->p.y);
            glVertex2d(fl->he2->vertex->p.x, fl->he2->vertex->p.y);
            glEnd();
        }
        glPointSize(6);
        glColor3d(0,0,1);
        for(auto&p: model->resPts){
            glBegin(GL_POINTS);
            glVertex2d(p.x, p.y);
            glEnd();
        }
    }

}

void GLWidget_2D::DrawGrid(){
    QRect geo = this->geometry();
    int x = gridsize/2, y = gridsize/2;
    //int x = 0, y = 0;
    glColor3d(0.85, 0.85, 0.85);
    while(x < geo.width() - gridsize/2){
        glBegin(GL_LINES);
        glVertex2d(x, 0);
        glVertex2d(x,geo.height());
        glEnd();
        x += gridsize;
    }

    while(y < geo.height() - gridsize/2){
        glBegin(GL_LINES);
        glVertex2d(0, y);
        glVertex2d(geo.width(), y);
        glEnd();
        y += gridsize;
    }

}



void GLWidget_2D::mousePressEvent(QMouseEvent *e){
    QPointF p = this->mapFromGlobal(QCursor::pos());
    glm::f64vec3 p_ongrid = SetOnGrid(p, gridsize);
    if(drawtype == "None") {return;}
    if(e->button() ==Qt::LeftButton){
        if(drawtype == "OutlineRectangle")model->drawOutline(p, 0, gridsize);
        else if(drawtype == "OutlinePolyline")model->drawOutline(p, 1, gridsize);
        else if(drawtype == "OutlinePolygon")model->drawOutline(p, 2, gridsize);
        else if(drawtype == "EditVertex")model->editOutlineVertex(p, gridsize, 0);

        else if(drawtype == "MoveOutline") model->outline->MoveOutline(p_ongrid);
        else if(drawtype == "OutlineConst") model->addConstraint(p, 0, gridsize, model->Axis4Const);
        else if(drawtype =="ConnectVertices")model->ConnectOutline(p, gridsize);
        else if(drawtype == "NewGradationMode" || drawtype =="FoldlineColor")addPoints_intplation(e, p);
        else if(drawtype == "FoldLine"){
            bool hasRulings = model->AddControlPoint_FL(p_ongrid, 0, curveDimention);
            update();
            if(model->outline->IsClosed() && hasRulings){
                bool res;
                //res = model->FL[0]->applyCurvedFolding(model->Faces, model->Edges, model->vertices, curveDimention);
                res = model->FL[0]->modify2DRulings(model->Faces, model->Edges, model->vertices, curveDimention);
                if(res) emit foldingSignals();
            }
        }
        else if(drawtype == "DeleteCurve"){
            model->SelectCurve(p);
            std::vector<int> n(1);
            n[0] = model->DeleteCurve();
            deleteCrvSignal(n);
            SelectedCurveIndex = -1;
        }
        if(drawtype == "RulingBezier" || drawtype == "RulingBspline" || drawtype == "RulingLine" || drawtype == "RulingArc"){model->AddControlPoint(p_ongrid, curveDimention, DivSize); }
        else if(drawtype == "InsertControlPoint"){
            if(SelectedCurveIndex == -1){
                model->SelectCurve(p);
                SelectedCurveIndex = model->IsSelectedCurve();
            }else{
                if(model->crvs[SelectedCurveIndex]->getCurveType() == 1 && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){
                    model->crvs[SelectedCurveIndex]->InsertControlPoint2(p_ongrid);
                    model->crvs[SelectedCurveIndex]->SetNewPoint();
                    model->crvs[SelectedCurveIndex]->Bspline(curveDimention,crvPtNum);
                    if(model->outline->IsClosed()){
                        model->crvs[SelectedCurveIndex]->BsplineRulings(model->outline,DivSize,crvPtNum, curveDimention);
                        if(!model->crvs[SelectedCurveIndex]->isempty){
                            model->addRulings();
                            model->deform();
                            //emit foldingSignals();
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
        if(drawtype == "RulingBezier" || drawtype == "RulingBspline" || drawtype == "RulingLine" || drawtype == "RulingArc"){
            if(SelectedCurveIndex != -1) model->crvs[SelectedCurveIndex]->eraseCtrlPt(curveDimention, crvPtNum);
        }
        else if(drawtype == "FoldLine"){
            bool hasRulings = model->AddControlPoint_FL(p_ongrid, 1, curveDimention);
        }

    }
    if(model->outline->IsClosed()){
        model->addRulings();
        model->deform();
        emit foldingSignals();
    }
    update();
}

void GLWidget_2D::mouseMoveEvent(QMouseEvent *e){
    if(drawtype == "None") {return;}
    QPointF p = this->mapFromGlobal(QCursor::pos());
    glm::f64vec3 p_ongrid = SetOnGrid(p, gridsize);
    if(drawtype == "OutlinePolygon"){
        if(model->outline->hasPtNum == 1){
            model->drawOutline(p, 2, gridsize, false);
            model->outline->hasPtNum = 1;
        }
    }
    //else if(drawtype == "OutlineRectangle" && movePt != -1)model->outline->getVertices()[movePt]->p = SetOnGrid(p);
    if(drawtype == "MoveOutline") model->outline->MoveOutline(p_ongrid);
    else if(drawtype == "EditVertex") model->editOutlineVertex(p, gridsize, 1);
    else if(drawtype == "NewGradationMode"){}
    //else if(model->getSelectedCurveIndex(p) == -1 && (drawtype == "DeleteCurve" || drawtype == "MoveControlPoint" || drawtype == "InsertControlPoint" || drawtype == "DeleteCntrlPt"))
    //SelectedCurveIndex = model->getSelectedCurveIndex(p);

    if(drawtype == "MoveControlPoint")model->MoveCurvePoint(p_ongrid,SelectedCurveIndex, movePt, curveDimention, DivSize);

    if(SelectedCurveIndex != -1){
        if(drawtype == "InsertControlPoint" &&  model->crvs[SelectedCurveIndex]->getCurveType() == 1 && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){//制御点の挿入(B-spline)
            model->crvs[SelectedCurveIndex]->InsertControlPoint2(p_ongrid);
        }
        
        if(model->outline->IsClosed() && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){
            //if(model->crvs[SelectedCurveIndex]->getCurveType() == 0)model->crvs[SelectedCurveIndex]->BezierRulings(model->outline->vertices,DivSize,crvPtNum);
            if(model->crvs[SelectedCurveIndex]->getCurveType() == 1) model->crvs[SelectedCurveIndex]->BsplineRulings(model->outline,DivSize,crvPtNum, curveDimention);
            if(model->crvs[SelectedCurveIndex]->getCurveType() == 2) model->crvs[SelectedCurveIndex]->LineRulings(model->outline,DivSize);
            if(model->crvs[SelectedCurveIndex]->getCurveType() == 3) model->crvs[SelectedCurveIndex]->ArcRulings(model->outline,DivSize);
            model->addRulings();
            model->deform();
        }
    }

    if(model->CrossDection4AllCurve()) emit foldingSignals();
    update();

}

void GLWidget_2D::mouseReleaseEvent(QMouseEvent * e){
    movePt = -1;
    QPointF p = this->mapFromGlobal(QCursor::pos());
    if(SelectedCurveIndex != -1 && !model->crvs.empty()) model->crvs[SelectedCurveIndex]->InsertPointSegment = -1;
    if(drawtype == "EditVertex") model->editOutlineVertex(p, gridsize, 2);
    if(model->outline->IsClosed())model->deform();
    update();
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
    model->deform();
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
    DiffWheel = (we->angleDelta().y() > 0) ? 1 : -1;
    if(drawtype == "FoldlineColor"){
        bool hasRulings = model->FL[0]->ChangeColor(model->outline, DiffWheel, curveDimention);
        emit ColorChangeFrom(0, model->FL[0]->getColor());
        if(hasRulings){
            bool res;
            model->FL[0]->applyCurvedFolding(model->Faces, model->Edges, model->vertices, curveDimention);
            emit foldingSignals();

        }
    }else if(drawtype == "NewGradationMode"){
        if(refHE == nullptr)return;

        model->setGradationValue(DiffWheel, refHE, InterpolationType, CurvePath);
        emit ColorChangeFrom(0, refHE->r->Gradation);
        model->deform();
        emit foldingSignals();
        if(isVisibleTo(gw)) emit CurvePathSet(CurvePath);
    }

    update();
}

void GLWidget_2D::addPoints_intplation(QMouseEvent *e, QPointF& p){
    glm::f64vec3 curPos{p.x(), p.y(), 0};
    double dist = 10;
    if(drawtype == "NewGradationMode"){
        refHE = nullptr;
        for(auto&crv: model->crvs){
            if(crv->isempty)continue;
            for(auto&rl: crv->Rulings){
                if(rl->IsCrossed != -1)continue;
                double d = glm::length(glm::cross((curPos - std::get<0>(rl->r)->p), (std::get<0>(rl->r)->p) - std::get<1>(rl->r)->p))/glm::length(std::get<0>(rl->r)->p - std::get<1>(rl->r)->p);
                if(d < dist){
                    dist = d; refHE = rl->he[0];
                }
            }
        }
        if(refHE == nullptr)return;

        model->setGradationValue(0, refHE, InterpolationType, CurvePath);
        emit ColorChangeFrom(0, refHE->r->Gradation);
        model->deform();
        emit foldingSignals();
    }else if(drawtype == "FoldlineColor"){
    }


    update();
}

void GLWidget_2D::ApplyNewGradationMode(){
    if(refHE == nullptr){std::cout << "you neeed to add gradation point"<<std::endl; return;}
    //QPointF p{-1,-1};
    //int color;
    //model->setGradationValue(p, 0, refMeshNum, color, InterpolationType, CurvePath);
    emit foldingSignals();
    if(isVisibleTo(gw)) emit CurvePathSet(CurvePath);
    update();
}

void GLWidget_2D::getGradationFromSlider(int val){
    if(refHE == nullptr) return;
    refHE->r->Gradation = val;
    emit ColorChangeFrom(2, val);
    model->setGradationValue(0, refHE, InterpolationType, CurvePath);
    emit foldingSignals();
    update();
}

void GLWidget_2D::OpenDebugWindwow(){
    if(!isVisibleTo(gw)) gw = new GToolWnd(this);
    emit CurvePathSet(CurvePath);

    gw->show();
}

HalfEdge *GLWidget_2D::assignment_refHE(){

    QPointF p = mapFromGlobal(QCursor::pos());
    glm::f64vec3 curPos{p.x(), p.y(), 0};
    HalfEdge *he = nullptr;
    double dist = 10;
    for(auto&crv: model->crvs){
        if(crv->isempty)continue;
        for(auto&rl: crv->Rulings){
            if(rl->IsCrossed != -1)continue;
            double d = glm::length(glm::cross((curPos - std::get<0>(rl->r)->p), (std::get<0>(rl->r)->p) - std::get<1>(rl->r)->p))/glm::length(std::get<0>(rl->r)->p - std::get<1>(rl->r)->p);
            if(d < dist){
                dist = d; he = rl->he[0];
            }
        }
    }
    return he;
}

glm::f64vec3 GLWidget_2D::SetOnGrid(QPointF& cursol, double gridsize){
    int x = (int)cursol.x() % (int)gridsize, y = (int)cursol.y() % (int)gridsize;
    x = (cursol.x() - x + gridsize/2);
    y = (cursol.y() - y + gridsize/2);
    return glm::f64vec3{x,y,0};
}
