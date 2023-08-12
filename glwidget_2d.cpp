#define PI 3.14159265359
#include "glwidget_2d.h"

GLWidget_2D::GLWidget_2D(QWidget *parent):QOpenGLWidget(parent)
{

    CurveList = {{CurveType::bezier3, PaintTool::Bezier_r}, {CurveType::bsp3, PaintTool::Bspline_r}, {CurveType::line, PaintTool::Line_r}, {CurveType::arc, PaintTool::Arc_r}};

    curveDimention = 3;
    drawtype = PaintTool::None;
    crvPtNum = 3000;
    maxDivSize = 3000;
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
    SmoothCurveIndex = -1;
    refV = nullptr;
    refL = nullptr;
    KeyEvent = -1;
    curvetype = CurveType::none;

    gridsize = 10;
    visibleGrid = 1;
    IsMVcolor_binary = false;
    eraseVec2d = false;
    IsEraseNonFoldEdge = false;
    visibleCurve = true;
    model = new Model(crvPtNum);

}
GLWidget_2D::~GLWidget_2D(){}

void GLWidget_2D::InitializeDrawMode(int state){
    if(state == 0)return; drawtype = PaintTool::None; SmoothCurveIndex = -1; emit SendNewActiveCheckBox(PaintTool::None);
}
void GLWidget_2D::VisualizeMVColor(bool state){IsMVcolor_binary = state;update();}

void GLWidget_2D::AddCurve(){
    std::vector<int>deleteIndex;
    CurveType _type;
    emit signalCurveType(_type);
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SmoothCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    for(auto& c: CurveList){
        if(std::get<0>(c) == _type){
            drawtype = std::get<1>(c);
            curvetype = _type;
        }
    }

    SmoothCurveIndex = model->AddNewCurve(curvetype, DivSize);
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::AddCurve);
}

void GLWidget_2D::InsertNewPoint(){
    drawtype = PaintTool::InsertCtrlPt;
    SmoothCurveIndex = -1;
    setMouseTracking(true);
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SmoothCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    emit SendNewActiveCheckBox(PaintTool::InsertCtrlPt);
}

void GLWidget_2D::EraseNonFoldEdge(bool state){
    IsEraseNonFoldEdge = state;
    update();
}

void GLWidget_2D::MoveCurvePt(){
    //curveNum = 1;
    drawtype = PaintTool::MoveCtrlPt;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SmoothCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    SmoothCurveIndex = -1;
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::MoveCtrlPt);
}

void GLWidget_2D::DrawOutlineRectangle(){
    //drawtype = -1;
    drawtype = PaintTool::Rectangle_ol;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SmoothCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    model->outline = new OUTLINE();
    model->outline->type = "Rectangle";
    this->setMouseTracking(false);
    emit SendNewActiveCheckBox(PaintTool::Rectangle_ol);
}

void GLWidget_2D::DrawOutlinePolygon(int state){
    if(state == 0)return;
    drawtype = PaintTool::Polygon_ol;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SmoothCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    model->outline = new OUTLINE();
    model->outline->type = "Polygon";
    for(auto&c: model->crvs)c->isempty = true;
    emit getEdgeNum();
    emit SendNewActiveCheckBox(PaintTool::Polygon_ol);
    this->setMouseTracking(true);
}

void GLWidget_2D::DrawOutlinePolyline(int state){
    if(state == 0)return;
    drawtype = PaintTool::Polyline_ol;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SmoothCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    model->outline = new OUTLINE();
    model->outline->type = "Polyline";
    for(auto&c: model->crvs)c->isempty = true;
    emit SendNewActiveCheckBox(PaintTool::Polyline_ol);
}

void GLWidget_2D::MoveOutline(int state){
    if(state == 0)return; drawtype = PaintTool::Move_ol;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SmoothCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    emit SendNewActiveCheckBox(PaintTool::Move_ol);
}

void GLWidget_2D::EditOutlineVertex(int state){
    drawtype = PaintTool::EditVertex_ol;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SmoothCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    emit SendNewActiveCheckBox(PaintTool::EditVertex_ol);
    update();
}

void GLWidget_2D::DeleteCtrlPt(){
    drawtype = PaintTool::DeleteCtrlPt;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    emit deleteCrvSignal(deleteIndex);
    SmoothCurveIndex = -1;
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::DeleteCtrlPt);
}

void GLWidget_2D::DeleteCurve(){
    drawtype = PaintTool::DeleteCurve;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SmoothCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::DeleteCurve);

}

void GLWidget_2D::changeFoldType(PaintTool state){
    drawtype = state;
    if(state != PaintTool::FoldlineColor){
        FoldLine *fl = new FoldLine(state);
        model->FL.push_back(fl);
        model->ChangeFoldLineState();
    }
    setMouseTracking(false);
    update();
}

void GLWidget_2D::setColor(){
    drawtype = PaintTool::SetColor;
}

void GLWidget_2D::switchGrid(){
    visibleGrid *= -1;
    update();
}

void GLWidget_2D::ConnectVertices(){drawtype = PaintTool::ConnectVertices_ol;}

void GLWidget_2D::switchGetmetricConstraint(int state){
    drawtype = PaintTool::Const_ol;
    constType = state;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
}

void GLWidget_2D::setNewGradationMode(){
    std::vector<int> deleteIndex;
    drawtype = PaintTool::NewGradationMode;
    model->Check4Param(curveDimention, deleteIndex); SmoothCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    update();
}
void GLWidget_2D::recieveNewEdgeNum(int num){model->outline->VerticesNum = num; update();}

void GLWidget_2D::ChangedDivSizeEdit(int n){
    this->DivSize = (n < 0)? DivSize: (maxDivSize < n)? maxDivSize: n;
    if(SmoothCurveIndex == -1 || model->crvs.empty()) return;
    bool res = false;
    if(curvetype == CurveType::bezier3){
       res = model->crvs[SmoothCurveIndex]->drawBezier(curveDimention, crvPtNum);
        //if(model->outline->isClosed  && model->crv->CurvePoints[0].pt != glm::f64vec2{-1,-1})model->crv->BezierRulings(model->outline->vertices,DivSize,crvPtNum);
    }
    if(curvetype == CurveType::bsp3){
        res = model->crvs[SmoothCurveIndex]->drawBspline(curveDimention, crvPtNum);
        if(model->outline->IsClosed() && res)model->crvs[SmoothCurveIndex]->BsplineRulings(model->outline,DivSize,crvPtNum, curveDimention);
    }
    if(curvetype == CurveType::line){
        res = model->crvs[SmoothCurveIndex]->drawLine();
        if(model->outline->IsClosed() && res)model->crvs[SmoothCurveIndex]->LineRulings(model->outline,DivSize);
    }

    if(model->outline->IsClosed() && model->CrossDection4AllCurve()){
        model->addRulings();
        model->deform();
        emit foldingSignals();
    }

    update();
}

void GLWidget_2D::changeSelectedCurve(int ind){SmoothCurveIndex = ind; drawtype = PaintTool::None; update();}

void GLWidget_2D::swapCrvsOnLayer(int n1, int n2){
    if(n1 == n2)return;
    std::iter_swap(model->crvs.begin() + n1, model->crvs.begin() + n2);
    if(SmoothCurveIndex == n1)SmoothCurveIndex = n2;
    else SmoothCurveIndex = n1;
    model->addRulings();
    model->deform();
    emit foldingSignals();
    update();
}

void GLWidget_2D::receiveNewLineWidth(double d){
    rulingWidth = d;
    update();
}

void GLWidget_2D::changeStopAtFF(bool state){IsStopAtFF = state;}
void GLWidget_2D::changeStopAtCon(bool state){IsStopAtCon = state;}
void GLWidget_2D::changeStopAtEq(bool state){IsStopAtEq = state;}
void GLWidget_2D::checkDevelopability(bool state){
    if(state)drawtype = PaintTool::CheckDevelopability;
    else{
        drawtype = PaintTool::None;
        refV = nullptr;
    }

}
void GLWidget_2D::Start4Debug_CF(){
    //if(model->Faces.size() < 2 )return;
    if(model->FL.empty())return;
    model->FL[0]->drawRulingInAllAngles(AllRulings);
    update();
    emit foldingSignals();
}


void GLWidget_2D::changeBetaValue(double val, int keyType){
    //if(model->Faces.size() < 2 || IsStop4Debug || model->FL.empty())return;
    if(model->FL.empty())return;
    if(model->FL.empty())model->FL.push_back(new FoldLine(PaintTool::FoldLine_test) );
    model->applyAAAMethod(val);
    emit foldingSignals();
    update();
}

void GLWidget_2D::initializeGL(){
    makeCurrent();
    initializeOpenGLFunctions();
    glClearColor(1.f, 1.f, 1.f, 1.f);
    QSize s = this->size();
    glViewport(0,0,s.width(),s.height());
}

void GLWidget_2D::paintGL(){
    makeCurrent();
    glClearColor(1.0,1.0,1.0,1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    QSize s = this->size();
    glLoadIdentity();
    glOrtho(-0.5, (float)s.width() -0.5, (float)s.height() -0.5, -0.5, -1, 1);

    if(visibleGrid == 1)DrawGrid();

    //折曲線の描画
    {
        if(visibleCurve){
            for(auto&fl: model->FL){
                glBegin(GL_LINE_STRIP);
                glColor3d(0,0,0); for(auto&v: fl->CurvePts)glVertex2d(v.x, v.y);
                glEnd();
                glColor3d(0,0.3,0.3);
                glPointSize(5);
                for(auto&v: fl->CtrlPts){
                    glBegin(GL_POINTS);
                    glVertex2d(v.x,v.y);
                    glEnd();
                }

                glColor3d(0,1,0);
                glPointSize(5);
                for(auto&h: fl->FoldingCurve){
                    if(IsEraseNonFoldEdge && !h.IsCalc)continue;
                    glBegin(GL_POINTS);
                    glVertex2d(h.first->p.x, h.first->p.y);
                    glEnd();
                }

            }
        }

    }

    glPolygonOffset(0.0,1.0);
    auto getcolor = [](double c, double a, double y)->double{
        if(y < a)return c/a * y/255.0;
        return ((255.0 - c)*(y - a)/(std::numbers::pi - a) + c)/255.0;
    };
    if(!model->FL.empty() && !model->FL[0]->FoldingCurve.empty()){
        for(auto& FL: model->FL){
            for(const auto& fc: FL->FoldingCurve){
                glColor3d(0, 0, 0);
                glPointSize(3.0f);
                glBegin(GL_POINTS);
                glVertex2d(fc.first->p.x, fc.first->p.y);
                glEnd();
            }
            //p0:描画するエッジの始点, p1: 描画するエッジの端点, p2: 片側の平面上の点, p3: もう片方の平面上の点
            auto DrawEdge = [&](Vertex *p0, Vertex *p1, Vertex *p2, Vertex *p3, double LineWidth, bool IsGradation, bool IsRuling){
                glLineWidth(LineWidth);
                glm::f64vec3 f_nv = glm::normalize(glm::cross(p1->p3 - p0->p3, p2->p3 - p0->p3)),fp_nv = glm::normalize(glm::cross(p3->p3 - p0->p3, p1->p3 - p0->p3));
                glm::f64vec3 SpinAxis = glm::normalize(p1->p3 - p0->p3);
                if(IsGradation){
                    double color = getcolor(model->ColorPt.color, model->ColorPt.angle, std::acos(glm::dot(f_nv, fp_nv)));
                    if(glm::dot(SpinAxis, glm::cross(f_nv, fp_nv)) <-1e-5){//mount
                       if(!IsMVcolor_binary)glColor3d(1,1.0 - color,1.0 - color);
                       else glColor3d(1,0,0);
                    }else if(glm::dot(SpinAxis, glm::cross(f_nv, fp_nv)) > 1e-5){//valley
                        if(!IsMVcolor_binary)glColor3d(1.0 - color,1.0 - color,1);
                        else glColor3d(0,0,1);
                    }else{
                        if(IsEraseNonFoldEdge &&  IsRuling)return;
                        glColor3d(0,0,0);
                    }
                }else{ glLineWidth(1.f); glColor3d(0,0,0); }
                glBegin(GL_LINES);
                glVertex2d(p1->p.x, p1->p.y);  glVertex2d(p0->p.x, p0->p.y);
                glEnd();

                glColor3d(0, 0, 0);
                glPointSize(3.0f); glBegin(GL_POINTS); glVertex2d(p0->p.x, p0->p.y);
                glEnd();

            };
            bool IsGradation = (drawtype == PaintTool::NewGradationMode)?true: false;
            for(int i = 1; i < (int)FL->FoldingCurve.size(); i++)
                DrawEdge(FL->FoldingCurve[i].first, FL->FoldingCurve[i-1].first, FL->FoldingCurve[i].second, FL->FoldingCurve[i].third, rulingWidth, IsGradation, false);
            for(int i = 1; i < (int)FL->FoldingCurve.size() - 1; i++){
                DrawEdge(FL->FoldingCurve[i].first, FL->FoldingCurve[i].second, FL->FoldingCurve[i+1].first, FL->FoldingCurve[i-1].first, rulingWidth, IsGradation, true);
                DrawEdge(FL->FoldingCurve[i].first, FL->FoldingCurve[i].third, FL->FoldingCurve[i-1].first, FL->FoldingCurve[i+1].first, rulingWidth, IsGradation, true);
            }


            glColor3d(0,0,0);
        }
        if(refV != nullptr){
            glBegin(GL_POINT);
            glPointSize(6.f);
            glColor3d(1,0,0);
            glVertex2d(refV->p.x, refV->p.y);
            glEnd();
        }
    }else{
        for(auto itr_r = model->Rulings.begin(); itr_r != model->Rulings.end(); itr_r++){
            glColor3d(0,0,0);
            glLineWidth(1);
            if((*itr_r)->et == EdgeType::r){
                if(drawtype == PaintTool::NewGradationMode){
                    glLineWidth(rulingWidth);
                    double color = getcolor(model->ColorPt.color, model->ColorPt.angle, (*itr_r)->color/255.0);
                    if(color > 1e-3){//mount
                       if(!IsMVcolor_binary)glColor3d(1,1.0 - color,1.0 - color);
                       else glColor3d(1,0,0);
                    }else if(color < -1e-3){//valley
                        if(!IsMVcolor_binary)glColor3d(1.0 - color,1.0 - color,1);
                        else glColor3d(0,0,1);
                    }else{
                        if(IsEraseNonFoldEdge && (*itr_r)->et == EdgeType::r)continue;
                        glColor3d(0,0,0);
                    }
                }else{ glLineWidth(1.f); glColor3d(0,0,0); }
                glBegin(GL_LINES);
                glVertex2d((*itr_r)->o->p.x, (*itr_r)->o->p.y);
                glVertex2d((*itr_r)->v->p.x, (*itr_r)->v->p.y);
                glEnd();

                glColor3d(0, 0, 0);
                glPointSize(3.0f);
                glBegin(GL_POINTS);
                glVertex2d((*itr_r)->o->p.x, (*itr_r)->o->p.y);
                glEnd();
                glBegin(GL_POINTS);
                glVertex2d((*itr_r)->v->p.x, (*itr_r)->v->p.y);
                glEnd();
                glColor3d(0, 0, 0);
            }
        }
    }

    //可展面の輪郭描画
    {
        glLineWidth(1);
        if(model->outline->type == "Rectangle" || model->outline->type == "Polyline"){
            glColor3d(0, 0, 0);
            glPointSize(3.0f);
            auto  Vertices = model->outline->getVertices();
            for(auto& v: Vertices){
                glBegin(GL_POINTS);
                glVertex2d(v->p.x, v->p.y);
                glEnd();
            }

            if(model->outline->IsClosed()){
                glBegin(GL_LINE_LOOP);
                for(auto& r: Vertices)glVertex2d(r->p.x, r->p.y);
                glEnd();
            }else{
                glBegin(GL_LINE_STRIP);
                for(auto& v: Vertices) glVertex2d(v->p.x, v->p.y);
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

    glLineWidth(1.0);
    if(!visibleCurve || drawtype == PaintTool::FoldLine_bezier)return;
    //ruling制御曲線の制御点
    if(drawtype != PaintTool::NewGradationMode){
        for(int i = 0; i < (int)model->crvs.size(); i++){
            if((i == model->getSelectedCurveIndex(mapFromGlobal(QCursor::pos()))) ||
                    ((drawtype == PaintTool::Bezier_r || drawtype == PaintTool::Bspline_r || drawtype == PaintTool::Line_r || drawtype == PaintTool::Arc_r)
                     && SmoothCurveIndex == i)
                    || drawtype == PaintTool::MoveCtrlPt || drawtype == PaintTool::InsertCtrlPt ||  drawtype == PaintTool::DeleteCtrlPt || drawtype == PaintTool::None){
                //glPolygonOffset(0.f,1.f);
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

        //ruling制御曲線の描画
        for(int i = 0; i < (int)model->crvs.size(); i++){
            if(i == SmoothCurveIndex)glColor3d(1, 0, 0);
            else glColor3d(0.4, 0.4, 0.4);
            glBegin(GL_LINE_STRIP);
            for(auto &c: model->crvs[i]->CurvePoints){
                glVertex2s(c.x, c.y);
            }
            glEnd();
        }
    }


    if(model->FL.empty())return;

    for(const auto&fl : model->FL){
        for(int i = 0; i < fl->NewRuling2d.size(); i++){
            auto r = fl->NewRuling2d[i];
            if(i % 3 == 0)glColor3d(1,0,0);
            else if(i % 3 == 1)glColor3d(0,1,0);
            else glColor3d(0,0,1);
            glBegin(GL_LINES);
            glVertex2d(r[0].x, r[0].y);
            glVertex2d(r[1].x, r[1].y);
            glEnd();
        }
    }
    glColor3d(0,0,0);

}

void GLWidget_2D::DrawGrid(){
    if(drawtype == PaintTool::NewGradationMode)return;
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

void GLWidget_2D::receiveKeyEvent(QKeyEvent *e){
    auto Poly_v = model->outline->getVertices();
    bool res;
    if(e->key() == Qt::Key_V)eraseVec2d = !eraseVec2d;
    if(e->key() == Qt::Key_A) visibleCurve = !visibleCurve;
    if(e->key() == Qt::Key_2){
        if(model->FL.empty())return;
        res = model->SplitRulings(curveDimention);
        if(res) emit foldingSignals();
    }
    if(e->key() == Qt::Key_Q){
        int type = 1;
        res = model->RevisionCrosPtsPosition();
        if(res) emit foldingSignals();
    }
    update();
}


void GLWidget_2D::mousePressEvent(QMouseEvent *e){
    QPointF p = this->mapFromGlobal(QCursor::pos());
    glm::f64vec3 p_ongrid = SetOnGrid(p, gridsize);
    if(drawtype == PaintTool::None) {return;}
    if(e->button() ==Qt::LeftButton){
        if(drawtype == PaintTool::Rectangle_ol)model->drawOutline(p, 0, gridsize);
        else if(drawtype == PaintTool::Polyline_ol){
            model->drawOutline(p, 1, gridsize);
        }
        else if(drawtype == PaintTool::Polygon_ol)model->drawOutline(p, 2, gridsize);
        else if(drawtype == PaintTool::EditVertex_ol)model->editOutlineVertex(p, gridsize, 0);

        else if(drawtype == PaintTool::Move_ol) model->outline->MoveOutline(p_ongrid);
        else if(drawtype == PaintTool::Const_ol) model->addConstraint(p, 0, gridsize, model->Axis4Const);
        else if(drawtype ==PaintTool::ConnectVertices_ol)model->ConnectOutline(p, gridsize);
        else if(drawtype == PaintTool::NewGradationMode || drawtype ==PaintTool::FoldlineColor)addPoints_intplation(e, p);
        else if(drawtype == PaintTool::FoldLine_bezier || drawtype == PaintTool::FoldLine_arc || drawtype == PaintTool::FoldLine_line ){
            bool hasRulings = model->AddControlPoint_FL(p_ongrid, 0, curveDimention);
            update();
        }
        else if(drawtype == PaintTool::FoldLine_test){
            std::cout << "can't use test now" << std::endl;
        }
        else if(drawtype == PaintTool::DeleteCurve){
            model->SelectCurve(p);
            std::vector<int> n(1);
            n[0] = model->DeleteCurve();
            emit deleteCrvSignal(n);
            SmoothCurveIndex = -1;
        }
        if(drawtype == PaintTool::Bezier_r || drawtype == PaintTool::Bspline_r || drawtype == PaintTool::Line_r || drawtype == PaintTool::Arc_r){
            bool res = model->AddControlPoint(p_ongrid, curveDimention, DivSize);
            if(res)model->addRulings();
        }
        else if(drawtype == PaintTool::InsertCtrlPt){
            if(SmoothCurveIndex == -1){
                model->SelectCurve(p);
                SmoothCurveIndex = model->IsSelectedCurve();
            }else{
                if(model->crvs[SmoothCurveIndex]->getCurveType() == CurveType::bsp3 && !model->crvs[SmoothCurveIndex]->isempty){
                    model->crvs[SmoothCurveIndex]->InsertControlPoint2(p_ongrid);
                    model->crvs[SmoothCurveIndex]->SetNewPoint();
                    model->crvs[SmoothCurveIndex]->drawBspline(curveDimention,crvPtNum);
                    if(model->outline->IsClosed()){
                        model->crvs[SmoothCurveIndex]->BsplineRulings(model->outline,DivSize,crvPtNum, curveDimention);
                        if(model->crvs[SmoothCurveIndex]->Rulings.front()->v != nullptr){
                            model->addRulings();
                            model->deform();
                            //emit foldingSignals();
                        }
                    }
                }
            }
        }
        else if(drawtype == PaintTool::DeleteCtrlPt){
            model->SelectCurve(p);
            model->DeleteControlPoint(p, curveDimention, DivSize);
        }else if(drawtype == PaintTool::MoveCtrlPt){
            SmoothCurveIndex = model->searchPointIndex(p, movePt, 0);

        }else if(drawtype == PaintTool::CheckDevelopability){
            int ind = -1; FoldLine *_fl = nullptr; refV = nullptr;
            double dist = 3.0;
            if(model->FL.empty())return;
            for(const auto&fl: model->FL){
                for(auto itr = fl->FoldingCurve.begin() + 1; itr != fl->FoldingCurve.end() - 1; itr++){
                    if(glm::length(glm::f64vec3{p.x(), p.y(), 0} - (*itr).first->p) < dist){
                        ind = std::distance(fl->FoldingCurve.begin(), itr); dist = glm::length(glm::f64vec3{p.x(), p.y(), 0} - (*itr).first->p);
                        _fl = fl;
                    }
                }
            }
            if(ind != -1){
                auto v4d = _fl->FoldingCurve[ind];
                glm::f64vec3 et = v4d.second->p3 -v4d.first->p3, er = _fl->FoldingCurve[ind-1].first->p3 -v4d.first->p3, eb = v4d.third->p3 -v4d.first->p3, el = _fl->FoldingCurve[ind+1].first->p3 -v4d.first->p3;
                et /= glm::length(et); er /= glm::length(er); eb /= glm::length(eb); el /= glm::length(el);
                double phi1 = std::acos(glm::dot(et, er)), phi2 = std::acos(glm::dot(et, el)), phi3 = std::acos(glm::dot(eb, el)), phi4 = std::acos(glm::dot(eb, er));
                refV = _fl->FoldingCurve[ind].first;
                std::cout << "developability  = " <<  abs(2.0*std::numbers::pi - phi1 - phi2 - phi3 - phi4) << ", phi1 = " << glm::degrees(phi1) << " , phi2 = " << glm::degrees(phi2) << ", phi3 = " << glm::degrees(phi3) << ", phi4 = " << glm::degrees(phi4) << std::endl;
            }

        }

    }else if(e->button() == Qt::RightButton){
        if(drawtype == PaintTool::Rectangle_ol || drawtype == PaintTool::Polyline_ol || drawtype == PaintTool::Polygon_ol) model->outline->eraseVertex();
        if(drawtype == PaintTool::Bezier_r || drawtype == PaintTool::Bspline_r || drawtype == PaintTool::Line_r || drawtype == PaintTool::Arc_r){
            if(SmoothCurveIndex != -1) model->crvs[SmoothCurveIndex]->eraseCtrlPt(curveDimention, crvPtNum);
        }
        else if(drawtype == PaintTool::FoldLine_line || drawtype == PaintTool::FoldLine_arc || drawtype == PaintTool::FoldLine_bezier){
            bool hasRulings = model->AddControlPoint_FL(p_ongrid, 1, curveDimention);
        }

    }
    if(model->outline->IsClosed()){
        emit foldingSignals();
    }
    update();
}

void GLWidget_2D::mouseMoveEvent(QMouseEvent *e){
    if(drawtype == PaintTool::None) {return;}
    QPointF p = this->mapFromGlobal(QCursor::pos());
    glm::f64vec3 p_ongrid = SetOnGrid(p, gridsize);
    if(drawtype == PaintTool::Polygon_ol){
        if(model->outline->hasPtNum == 1){
            model->drawOutline(p, 2, gridsize, false);
            model->outline->hasPtNum = 1;
        }
    }
    if(drawtype == PaintTool::Move_ol) model->outline->MoveOutline(p_ongrid);
    else if(drawtype == PaintTool::EditVertex_ol){
        model->editOutlineVertex(p, gridsize, 1);
        model->addRulings();
        model->deform();
        emit foldingSignals();
    }
    else if(drawtype == PaintTool::NewGradationMode){}

    if(drawtype == PaintTool::MoveCtrlPt)model->MoveCurvePoint(p_ongrid,SmoothCurveIndex, movePt, curveDimention, DivSize);
    else if(drawtype == PaintTool::FoldLine_bezier){
        //bool res = model->updateSplitRulings(model->FL[FoldCurveIndex], curveDimention);
        //if(res) emit foldingSignals();
        update();
    }
    if(SmoothCurveIndex != -1){
        if(drawtype == PaintTool::InsertCtrlPt &&  model->crvs[SmoothCurveIndex]->getCurveType() == CurveType::bsp3 && !model->crvs[SmoothCurveIndex]->isempty){//制御点の挿入(B-spline)
            model->crvs[SmoothCurveIndex]->InsertControlPoint2(p_ongrid);
        }
        
        if(drawtype != PaintTool::NewGradationMode && model->outline->IsClosed() && !model->crvs[SmoothCurveIndex]->isempty){

            if(model->crvs[SmoothCurveIndex]->getCurveType() == CurveType::bsp3) model->crvs[SmoothCurveIndex]->BsplineRulings(model->outline,DivSize,crvPtNum, curveDimention);
            if(model->crvs[SmoothCurveIndex]->getCurveType() == CurveType::line) model->crvs[SmoothCurveIndex]->LineRulings(model->outline,DivSize);
            if(model->crvs[SmoothCurveIndex]->getCurveType() == CurveType::arc) model->crvs[SmoothCurveIndex]->ArcRulings(model->outline,DivSize);
            model->addRulings();
            model->deform();
        }
    }

    if(model->CrossDection4AllCurve()) emit foldingSignals();
    update();

}

void GLWidget_2D::mouseReleaseEvent(QMouseEvent * e){
    movePt = -1;
    QPointF p = mapFromGlobal(QCursor::pos());
    if(SmoothCurveIndex != -1 && !model->crvs.empty()) model->crvs[SmoothCurveIndex]->InsertPointSegment = -1;
    if(drawtype == PaintTool::EditVertex_ol){
        model->editOutlineVertex(p, gridsize, 2);
        if(model->outline->IsClosed())model->deform();
    }
    update();
}

void GLWidget_2D::cb_ApplyCurveEvent(){
    if(KeyEvent == 0){
        std::cout << "finished " << std::endl; KeyEvent = -1;
    }else if(KeyEvent == -1){
        if(SmoothCurveIndex == -1) std::cout << "no curve is selected"<<std::endl;
    }
    update();
}

void GLWidget_2D::cb_DeleteCurve(){
    KeyEvent = 1;
    if(SmoothCurveIndex != -1){
        model->DeleteCurve();
    }
    update();
}

void GLWidget_2D::Reset(){
    delete model;
    tmp_c.clear();
    model = new Model(crvPtNum);
    emit SendNewActiveCheckBox(PaintTool::Reset);
    SmoothCurveIndex = -1;
    emit foldingSignals();
    update();
}

void GLWidget_2D::receiveColors(){
    model->deform();
    emit foldingSignals();
    update();
}

/*
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
}*/
void GLWidget_2D::wheelEvent(QWheelEvent *we){
    DiffWheel = (we->angleDelta().y() > 0) ? 1 : -1;
    if(drawtype == PaintTool::FoldlineColor){
        //emit ColorChangeFrom(0, model->FL[FoldCurveIndex]->getColor());
    }else if(drawtype == PaintTool::NewGradationMode){
        if(refL == nullptr || std::find(model->Rulings.begin(), model->Rulings.end(), refL) == model->Rulings.end()){
            std::cout<<"there is no referenced Ruling"<<std::endl;
            return;
        }
        model->setGradationValue(DiffWheel, refL, InterpolationType, CurvePath);
        emit ColorChangeFrom(0, refL->color);
        model->deform();
        if(isVisibleTo(gw)) emit CurvePathSet(CurvePath);
    }
    emit foldingSignals();
    update();
}

void GLWidget_2D::addPoints_intplation(QMouseEvent *e, QPointF& p){
    glm::f64vec3 curPos{p.x(), p.y(), 0};
    double dist = 10;
    if(drawtype == PaintTool::NewGradationMode){
        refL = nullptr;
        for(auto* r: model->Rulings){
            if(r->et == EdgeType::ol || r->et == EdgeType::cl)continue;
            double d = glm::length(glm::cross((curPos - r->o->p), r->o->p - r->v->p))/glm::length(r->o->p - r->v->p);
            if(d < dist){
                dist = d; refL = r;
            }
        }
        if(refL == nullptr || std::find(model->Rulings.begin(), model->Rulings.end(), refL) == model->Rulings.end()){std::cout<<"not found" << std::endl; return;}
        model->setGradationValue(0, refL, InterpolationType, CurvePath);
        emit ColorChangeFrom(0, refL->color);
        model->deform();
        if(model->outline->IsClosed() && !model->FL.empty()){
            //auto oriedge = model->outline->getEdges();
            //model->FL[0]->modify2DRulings(model->Faces, model->Edges, model->vertices, oriedge, curveDimention);
        }
    }else if(drawtype == PaintTool::FoldlineColor){}
    emit foldingSignals();
    update();
}

void GLWidget_2D::ApplyNewGradationMode(){
    if(refL == nullptr){std::cout << "you neeed to add gradation point"<<std::endl; return;}
    emit foldingSignals();
    if(isVisibleTo(gw)) emit CurvePathSet(CurvePath);
    update();
}

void GLWidget_2D::getGradationFromSlider(int val){
    if(refL == nullptr) return;
    refL->color = val;
    emit ColorChangeFrom(2, val);
    model->setGradationValue(0, refL, InterpolationType, CurvePath);
    emit foldingSignals();
    update();
}

void GLWidget_2D::OpenDebugWindwow(){
    if(!isVisibleTo(gw)) gw = new GToolWnd(this);
    emit CurvePathSet(CurvePath);

    gw->show();
}

int GLWidget_2D::assignment_refL(){

    QPointF p = mapFromGlobal(QCursor::pos());
    glm::f64vec3 curPos{p.x(), p.y(), 0};
    int ind = -1;
    double dist = 10;
    for(int i = 0; i < (int)model->Rulings.size(); i++){
        Line *r = model->Rulings[i];
        if(r->et != EdgeType::r)continue;
        double d = glm::length(glm::cross((curPos - r->o->p), r->o->p - r->v->p))/glm::length(r->o->p - r->v->p);
        if(d < dist){
            dist = d; ind = i;
        }
    }
    return ind;
}

glm::f64vec3 GLWidget_2D::SetOnGrid(QPointF& cursol, double gridsize){
    int x = (int)cursol.x() % (int)gridsize, y = (int)cursol.y() % (int)gridsize;
    x = (cursol.x() - x + gridsize/2);
    y = (cursol.y() - y + gridsize/2);
    return glm::f64vec3{x,y,0};
}
