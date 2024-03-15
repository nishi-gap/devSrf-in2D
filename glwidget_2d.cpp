#include "glwidget_2d.h"
#include <GL/glu.h>

GLWidget_2D::GLWidget_2D(QWidget *parent):QOpenGLWidget(parent)
{
    curveDimention = 3;
    drawtype = PaintTool::None;
    crvPtNum = 3000;
    maxDivSize = 3000;
    standardDist = 2.5;

    DiffWheel = 0;
    movePt = -1;
    MoveCrvIndex = {-1, -1};
    refLine = std::shared_ptr<Line>(nullptr);
    KeyEvent = -1;
    curvetype = CurveType::none;
    basePoint = QPointF(-1,-1);
    gridsize = 10;
    visibleGrid = 1;
    IsMVcolor_binary = false;
    IsEraseNonFoldEdge = false;
    IsAffinMoved = false;
    visibleCurve = true;
    model.push_back(std::make_shared<Model>(crvPtNum));
    difCursol_x = difCursol_y = 0.0; camscale = -1;
}

GLWidget_2D::~GLWidget_2D(){}

void GLWidget_2D::InitializeDrawMode(){
    drawtype = PaintTool::None; MoveCrvIndex = {-1, -1};
    emit SendNewActiveCheckBox(PaintTool::None);
}
void GLWidget_2D::VisualizeMVColor(bool state){IsMVcolor_binary = state;update();}

void GLWidget_2D::switch2AffinMode(){ drawtype = PaintTool::AffinTrans;}

void GLWidget_2D::switch2VisibleCurve(){
    visibleCurve = !visibleCurve;
    update();}

void GLWidget_2D::CopyCurveObj(){
    if(model.back()->refCreases == nullptr || std::find(model.back()->Creases.begin(), model.back()->Creases.end(), model.back()->refCreases) == model.back()->Creases.end())IsCopied = false;
    IsCopied = true;
}

void GLWidget_2D::PasteCurveObj(){
    if(!IsCopied || std::find(model.back()->Creases.begin(), model.back()->Creases.end(), model.back()->refCreases) == model.back()->Creases.end())return;
    auto newCrease = model.back()->refCreases->deepCopy();
    //上方向に10移動させる
    double y = -30;
    for(auto&p: newCrease->CtrlPts)p.y() += y;
    newCrease->setCurve(3);
    newCrease->FoldingCurve.clear(); newCrease->a_flap = -1;
    model.back()->Creases.push_back(newCrease);
}

void GLWidget_2D::AddCurve(){
    std::vector<int>deleteIndex;

    emit signalCurveType(curvetype);
    model.back()->Check4Param(curveDimention, deleteIndex);
    if(model.back()->RulingCurve.empty()) MoveCrvIndex = {-1,-1};
    emit deleteCrvSignal(deleteIndex);

    if(curvetype == CurveType::bsp3)drawtype = PaintTool::Bspline;
    else if(curvetype == CurveType::line)drawtype = PaintTool::Line;
    else if(curvetype == CurveType::arc)drawtype = PaintTool::Arc;

    MoveCrvIndex[0] = model.back()->AddNewCurve(curvetype, DivSize);
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::AddCurve);
}

void GLWidget_2D::InsertNewPoint(){
    drawtype = PaintTool::InsertCtrlPt;
    MoveCrvIndex[0] = -1;
    setMouseTracking(true);
    std::vector<int>deleteIndex;
    model.back()->Check4Param(curveDimention, deleteIndex);
    if(model.back()->RulingCurve.empty()) MoveCrvIndex[0] = -1;
    emit deleteCrvSignal(deleteIndex);
    emit SendNewActiveCheckBox(PaintTool::InsertCtrlPt);
}

void GLWidget_2D::EraseNonFoldEdge(bool state){
    IsEraseNonFoldEdge = state;
    update();
}

void GLWidget_2D::MoveCurvePt(){
    drawtype = PaintTool::MoveCtrlPt;
    std::vector<int>deleteIndex;
    model.back()->Check4Param(curveDimention, deleteIndex);
    emit deleteCrvSignal(deleteIndex);
    MoveCrvIndex = {-1, -1};
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::MoveCtrlPt);
}

void GLWidget_2D::DrawOutlineRectangle(){
    //drawtype = -1;
    drawtype = PaintTool::Rectangle;
    std::vector<int>deleteIndex;
    model.back()->Check4Param(curveDimention, deleteIndex);
    if(model.back()->RulingCurve.empty()) MoveCrvIndex[0] = -1;
    if(model.back()->Creases.empty())MoveCrvIndex[1] = -1;
    emit deleteCrvSignal(deleteIndex);
    model.back()->outline = std::make_shared<OUTLINE>();
    model.back()->outline->type = "Rectangle";
    this->setMouseTracking(false);
    emit SendNewActiveCheckBox(PaintTool::Rectangle);
}

void GLWidget_2D::DrawOutlinePolygon(int state){
    if(state == 0)return;
    drawtype = PaintTool::Polygon;
    std::vector<int>deleteIndex;
    model.back()->Check4Param(curveDimention, deleteIndex);
    if(model.back()->RulingCurve.empty()) MoveCrvIndex[0] = -1;
    if(model.back()->Creases.empty()) MoveCrvIndex[1] = -1;
    emit deleteCrvSignal(deleteIndex);
    model.back()->outline = std::make_shared<OUTLINE>();
    model.back()->outline->type = "Polygon";
    for(auto&c: model.back()->RulingCurve)c->isempty = true;
    emit getEdgeNum();
    emit SendNewActiveCheckBox(PaintTool::Polygon);
    this->setMouseTracking(true);
}

void GLWidget_2D::DrawOutlinePolyline(int state){
    if(state == 0)return;
    drawtype = PaintTool::Polyline;
    std::vector<int>deleteIndex;
    model.back()->Check4Param(curveDimention, deleteIndex);
    if(model.back()->RulingCurve.empty()) MoveCrvIndex[0] = -1;
    if(model.back()->Creases.empty()) MoveCrvIndex[1] = -1;
    emit deleteCrvSignal(deleteIndex);
    model.back()->outline = std::make_shared<OUTLINE>();
    model.back()->outline->type = "Polyline";
    for(auto&c: model.back()->RulingCurve)c->isempty = true;
    emit SendNewActiveCheckBox(PaintTool::Polyline);
}

void GLWidget_2D::DeleteCtrlPt(){
    drawtype = PaintTool::DeleteCtrlPt;
    std::vector<int>deleteIndex;
    model.back()->Check4Param(curveDimention, deleteIndex);
    emit deleteCrvSignal(deleteIndex);
    MoveCrvIndex = {-1, -1};
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::DeleteCtrlPt);
}

void GLWidget_2D::DeleteCurve(){
    drawtype = PaintTool::DeleteCurve;
    std::vector<int>deleteIndex;
    model.back()->Check4Param(curveDimention, deleteIndex);
    if(model.back()->RulingCurve.empty()) MoveCrvIndex[0] = -1;
    if(model.back()->Creases.empty()) MoveCrvIndex[1] = -1;
    emit deleteCrvSignal(deleteIndex);
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::DeleteCurve);

}

void GLWidget_2D::AddNewCrease(){
    drawtype = PaintTool::Crease;
    IsMVcolor_binary = true;
    model.back()->RemoveUnable2GenCurve();

    std::shared_ptr<FoldLine> crease = std::make_shared<FoldLine>(PaintTool::Crease);
    model.back()->Creases.push_back(crease);
    model.back()->ChangeFoldLineState();
    setMouseTracking(false);
    update();
}

void GLWidget_2D::setColor(){ drawtype = PaintTool::SetColor; }

void GLWidget_2D::switchGrid(){
    visibleGrid *= -1;
    update();
}


void GLWidget_2D::setNewGradationMode(){
    std::vector<int> deleteIndex;
    drawtype = PaintTool::NewGradationMode;
    model.back()->Check4Param(curveDimention, deleteIndex); MoveCrvIndex = {-1, -1};
    emit deleteCrvSignal(deleteIndex);
    update();
}
void GLWidget_2D::recieveNewEdgeNum(int num){model.back()->outline->VerticesNum = num; update();}

void GLWidget_2D::ChangedDivSizeEdit(int n){
    this->DivSize = (n < 0)? DivSize: (maxDivSize < n)? maxDivSize: n;
    if(MoveCrvIndex[0] == -1 || model.back()->RulingCurve.empty()) return;
    bool res = false;
    if(curvetype == CurveType::bezier3){
       res = model.back()->RulingCurve[MoveCrvIndex[0]]->drawBezier(curveDimention, crvPtNum);
        //if(model.back()->outline->isClosed  && model.back()->crv->CurvePoints[0].pt != glm::f64vec2{-1,-1})model.back()->crv->BezierRulings(model.back()->outline->vertices,DivSize,crvPtNum);
    }
    if(curvetype == CurveType::bsp3){
        res = model.back()->RulingCurve[MoveCrvIndex[0]]->drawBspline(curveDimention, crvPtNum);
        if(model.back()->outline->IsClosed() && res)model.back()->RulingCurve[MoveCrvIndex[0]]->BsplineRulings(model.back()->outline,DivSize,crvPtNum, curveDimention);
    }
    if(curvetype == CurveType::line){
        res = model.back()->RulingCurve[MoveCrvIndex[0]]->drawLine();
        if(model.back()->outline->IsClosed() && res)model.back()->RulingCurve[MoveCrvIndex[0]]->LineRulings(model.back()->outline,DivSize);
    }

    if(model.back()->outline->IsClosed() && model.back()->CrossDection4AllCurve()){
        model.back()->addRulings();
        model.back()->deform();
        emit foldingSignals();
    }

    update();
}

void GLWidget_2D::stashcurrentstate(){
    auto m = model.back()->stashcurrentstate();
    model.push_back(m);
    update();
    emit foldingSignals();
    qDebug()<<"stashed size is " << model.size();
}

void GLWidget_2D::back2befstate(){
    if(model.empty())return;
    model.pop_back();
    if(model.empty())model.push_back(std::make_shared<Model>(crvPtNum));
    qDebug()<<"back and stashed size is " << model.size();
    update();
    emit foldingSignals();
}

void GLWidget_2D::changeSelectedCurve(int ind){MoveCrvIndex[0] = ind; drawtype = PaintTool::None; update();}

void GLWidget_2D::swapCrvsOnLayer(int n1, int n2){
    if(n1 == n2)return;
    std::iter_swap(model.back()->RulingCurve.begin() + n1, model.back()->RulingCurve.begin() + n2);
    if(MoveCrvIndex[0] == n1)MoveCrvIndex[0] = n2;
    else MoveCrvIndex[0] = n1;
    model.back()->addRulings();
    model.back()->deform();
    emit foldingSignals();
    update();
}

void GLWidget_2D::receiveNewLineWidth(double d){
    rulingWidth = d;
    update();
}

void GLWidget_2D::changeflapgnle(double val){
    if(model.back()->Creases.empty())return;
    model.back()->Creases.back()->applyAAAMethod(model.back()->outline->vertices, val);
    emit foldingSignals();
    update();
}


void GLWidget_2D::initializeGL(){
    makeCurrent();
    initializeOpenGLFunctions();
    glClearColor(1.f, 1.f, 1.f, 1.f);
    QSize s = this->size();
    glViewport(0,0,s.width(),s.height());

    //glEnable(GL_DEPTH_TEST);
    // アルファブレンディングを有効にする
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void GLWidget_2D::paintGL(){
    makeCurrent();
    glClearColor(1.0,1.0,1.0,1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-0.5, (float)size().width() -0.5, (float)size().height() -0.5, -0.5, -1, 100);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslated(difCursol_x, difCursol_y, camscale);

    if(visibleGrid == 1)DrawGrid();

    //折曲線の描画
    {
        if(visibleCurve || drawtype == PaintTool::DeleteCtrlPt || drawtype == PaintTool::MoveCtrlPt || drawtype == PaintTool::None){
            for(auto & crease: model.back()->Creases){
                glColor3d(0,0.3,0.3);
                glPointSize(5);
                for(auto&v: crease->CtrlPts){
                    glBegin(GL_POINTS);
                    glVertex3d(v.x(),v.y(), 0);
                    glEnd();
                }

                glEnable(GL_LINE_STIPPLE);
                glLineStipple(1 , 0xF0F0);
                glBegin(GL_LINE_STRIP);
                glColor3d(0,0,0);
                for(auto&v: crease->CtrlPts)glVertex3d(v.x(), v.y(), 0);
                glEnd();
                glDisable(GL_LINE_STIPPLE);

                if(crease->CurvePts.empty())continue;
                glBegin(GL_LINE_STRIP);
                for(auto&v: crease->CurvePts) glVertex3d(v.x(),v.y(), 0);
                glEnd();

                glColor3d(0,1,0);
                glPointSize(5);
                for(auto&h: crease->FoldingCurve){
                    if(IsEraseNonFoldEdge && !h->IsCalc)continue;
                    glBegin(GL_POINTS);
                    glVertex3d(h->first->p.x(), h->first->p.y(),0);
                    glEnd();
                }
            }
        }

    }

    if(drawtype == PaintTool::AffinTrans){
        if(affinmode != 0){
            glPointSize(5);
            glBegin(GL_POINTS);
            glVertex2d(basePoint.x(), basePoint.y());
            glEnd();
        }
    }

    glPolygonOffset(0.0,1.0);
    auto getcolor = [](double c, double a, double y)->double{
        if(y < a)return c/a * y/255.0;
        return ((255.0 - c)*(y - a)/(std::numbers::pi - a) + c)/255.0;
    };
    //折曲線上の4価頂点の描画
    for(auto& crease: model.back()->Creases){
        if(crease->FoldingCurve.empty())continue;
        std::vector<int> Vertices_Ind;
        for(int i = 0; i < (int)crease->FoldingCurve.size(); i++){if(crease->FoldingCurve[i]->IsCalc)Vertices_Ind.push_back(i);}
        for(const auto& fc: crease->FoldingCurve){
            glColor3d(0, 0, 0);
            glPointSize(3.0f);
            glBegin(GL_POINTS);
            glVertex3d(fc->first->p.x(), fc->first->p.y(),0);
            glEnd();
        }
        //p0:描画するエッジの始点, p1: 描画するエッジの端点, p2: 片側の平面上の点, p3: もう片方の平面上の点
        auto DrawEdge = [&](const std::shared_ptr<Vertex>& p0, const std::shared_ptr<Vertex>& p1, const std::shared_ptr<Vertex>& p2, const std::shared_ptr<Vertex>& p3, double LineWidth, bool IsGradation, bool IsRuling, bool IsReverse, bool DashedLine){
            glLineWidth(LineWidth);
            Eigen::Vector3d SpinAxis = (p1->p3 - p0->p3).normalized();
            Eigen::Vector3d f_nv = (SpinAxis.cross(p2->p3 - p0->p3)).normalized(),fp_nv = ((p3->p3 - p0->p3).cross(SpinAxis)).normalized();
            if(IsGradation){
                double color = getcolor(model.back()->ColorPt.color, model.back()->ColorPt.angle, std::acos(f_nv.dot(fp_nv)));
                double th = -1e-5;
                double f = SpinAxis.dot(f_nv.cross(fp_nv));
                if(IsReverse) f *= -1;
                if(f <th){//mount
                    if(!IsMVcolor_binary)glColor3d(1, 1.0 - color, 1.0 - color);
                    else glColor3d(1,0,0);
                }else if(f > -th){//valley
                    if(!IsMVcolor_binary)glColor3d(1.0 - color,1.0 - color,1);
                    else glColor3d(0,0,1);
                }else{
                    if(IsEraseNonFoldEdge &&  IsRuling)return;
                    glColor3d(0,0,0);
                }
            }else{ glLineWidth(1.f); glColor3d(0,0,0); }
            if(DashedLine){
                glEnable(GL_LINE_STIPPLE);
                glLineStipple(1 , 0xF0F0);
            }
            glBegin(GL_LINES);
            glVertex3d(p1->p.x(), p1->p.y(),0);  glVertex3d(p0->p.x(), p0->p.y(),0);
            glEnd();
            if(DashedLine)glDisable(GL_LINE_STIPPLE);

            glColor3d(0, 0, 0);
            glPointSize(3.0f); glBegin(GL_POINTS); glVertex3d(p0->p.x(), p0->p.y(),0);
            glEnd();

        };
        bool IsGradation = (drawtype == PaintTool::NewGradationMode || drawtype == PaintTool::Crease)?true: false;

        for(int i = 1; i < (int)crease->FoldingCurve.size()-1; i++)
            DrawEdge(crease->FoldingCurve[i]->first, crease->FoldingCurve[i-1]->first, crease->FoldingCurve[i]->second, crease->FoldingCurve[i]->third, rulingWidth, IsGradation, false, false, false);
        DrawEdge(crease->FoldingCurve.end()[-2]->first, crease->FoldingCurve.back()->first, crease->FoldingCurve.end()[-2]->second, crease->FoldingCurve.end()[-2]->third, rulingWidth, IsGradation, true, true, false);
        for(int i = 1; i < (int)crease->FoldingCurve.size() - 1; i++){
            bool DashedLine = (crease->FoldingCurve[i]->IsCalc)? false: true;
            DrawEdge(crease->FoldingCurve[i]->first, crease->FoldingCurve[i]->second, crease->FoldingCurve[i+1]->first, crease->FoldingCurve[i-1]->first, rulingWidth, IsGradation, true, false, DashedLine);
            DrawEdge(crease->FoldingCurve[i]->first, crease->FoldingCurve[i]->third, crease->FoldingCurve[i-1]->first, crease->FoldingCurve[i+1]->first, rulingWidth, IsGradation, true, false, DashedLine);
        }



        glColor3d(0,0,0);
    }

    //折曲線との交差がないrulingの描画
    for(auto itr_r = model.back()->Rulings.begin(); itr_r != model.back()->Rulings.end(); itr_r++){
        if((*itr_r)->hasCrossPoint)continue;
        glColor3d(0,0,0);
        glLineWidth(1);
        if((*itr_r)->et == EdgeType::r){
            if(drawtype == PaintTool::NewGradationMode|| drawtype == PaintTool::Crease){
                glLineWidth(rulingWidth);
                double color = getcolor(model.back()->ColorPt.color, model.back()->ColorPt.angle, (*itr_r)->color/255.0);
                if(color > 1e-3){//mount
                    if(!IsMVcolor_binary)glColor4d(1, 0, 0, color);
                    else glColor3d(1,0,0);
                }else if(color < -1e-3){//valley
                    if(!IsMVcolor_binary)glColor4d(0, 0, 1, -color);
                    else glColor3d(0,0,1);
                }else{
                    if(IsEraseNonFoldEdge && (*itr_r)->et == EdgeType::r)continue;
                    glColor3d(0,0,0);
                }
            }else{ glLineWidth(1.f); glColor3d(0,0,0); }
            if((*itr_r)->IsCrossed == 0)glColor3d(0,1,0);
            glBegin(GL_LINES);
            glVertex3d((*itr_r)->o->p.x(), (*itr_r)->o->p.y(),0);
            glVertex3d((*itr_r)->v->p.x(), (*itr_r)->v->p.y(),0);
            glEnd();

            glColor3d(0, 0, 0);
            glPointSize(3.0f);
            glBegin(GL_POINTS);
            glVertex3d((*itr_r)->o->p.x(), (*itr_r)->o->p.y(),0);
            glEnd();
            glBegin(GL_POINTS);
            glVertex3d((*itr_r)->v->p.x(), (*itr_r)->v->p.y(),0);
            glEnd();
            glColor3d(0, 0, 0);
        }
    }


    //可展面の輪郭描画
    {
        glLineWidth(1);
        if(model.back()->outline->type == "Rectangle" || model.back()->outline->type == "Polyline"){
            glColor3d(0, 0, 0);
            glPointSize(3.0f);
            auto  Vertices = model.back()->outline->getVertices();
            for(auto& v: Vertices){
                glBegin(GL_POINTS);
                glVertex3d(v->p.x(), v->p.y(),0);
                glEnd();
            }

            if(model.back()->outline->IsClosed()){
                glBegin(GL_LINE_LOOP);
                for(auto& r: Vertices)glVertex3d(r->p.x(), r->p.y(),0);
                glEnd();
            }else{
                glBegin(GL_LINE_STRIP);
                for(auto& v: Vertices) glVertex3d(v->p.x(), v->p.y(),0);
                glEnd();
            }
        }
        if(model.back()->outline->type == "Polygon"){
            if(model.back()->outline->hasPtNum > 0){
                glColor3d(0, 0, 0);
                glPointSize(5.0f);
                glBegin(GL_POINTS);
                glVertex3d(model.back()->outline->origin.x(), model.back()->outline->origin.y(),0);
                glEnd();

                glColor3d(0, 0, 0);
                glPointSize(4.0f);
                for(auto& v: model.back()->outline->getVertices()){
                    glBegin(GL_POINTS);
                    glVertex3d(v->p.x(), v->p.y(),0);
                    glEnd();
                }

                glBegin(GL_LINE_LOOP);
                for(auto& r: model.back()->outline->getVertices())glVertex3d(r->p.x(), r->p.y(),0);
                glEnd();
            }

        }

    }

    glLineWidth(1.0);
    if(!visibleCurve || drawtype == PaintTool::Crease)return;
    //ruling制御曲線の制御点
    if(drawtype != PaintTool::NewGradationMode){
        for(int i = 0; i < (int)model.back()->RulingCurve.size(); i++){
            if((i == model.back()->getSelectedCurveIndex(mapFromGlobal(QCursor::pos()))[0]) ||
                    ((drawtype == PaintTool::Bspline || drawtype == PaintTool::Line || drawtype == PaintTool::Arc)
                     && MoveCrvIndex[0] == i)
                    || drawtype == PaintTool::MoveCtrlPt || drawtype == PaintTool::InsertCtrlPt ||  drawtype == PaintTool::DeleteCtrlPt || drawtype == PaintTool::None){
                //glPolygonOffset(0.f,1.f);
                for (auto& c: model.back()->RulingCurve[i]->ControllPoints) {
                    glColor3d(1, 0, 0);
                    glPointSize(5.0f);
                    glBegin(GL_POINTS);
                    glVertex3d( c.x(), c.y(),0);

                    glEnd();
                }
                glEnable(GL_LINE_STIPPLE);
                glLineStipple(1 , 0xF0F0);
                glBegin(GL_LINE_STRIP);
                glColor3d(0.4, 0.4, 0.4);
                for (auto& c: model.back()->RulingCurve[i]->ControllPoints)glVertex3d( c.x(), c.y(),0);
                glEnd();
                glDisable(GL_LINE_STIPPLE);
                glPolygonOffset(0.f,0.f);
            }
        }

        //ruling制御曲線の描画
        for(int i = 0; i < (int)model.back()->RulingCurve.size(); i++){
            if(i == MoveCrvIndex[0])glColor3d(1, 0, 0);
            else glColor3d(0.4, 0.4, 0.4);
            glBegin(GL_LINE_STRIP);
            for(auto &c: model.back()->RulingCurve[i]->CurvePoints){
                glVertex2s(c.x(), c.y());
            }
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

void GLWidget_2D::mousePressEvent(QMouseEvent *e){
    QPointF p = mapFromGlobal(QCursor::pos());
    QPointF p_ongrid = SetOnGrid(p, gridsize, this->size());
    if(drawtype == PaintTool::None) {
        if(e->button() == Qt::LeftButton){
            IsLeftClicked = true;
            model.back()->detectClickedObj(p);
        }
        if(e->button() == Qt::RightButton){IsRightClicked = true; difCursol_x = difCursol_y = 0.0;}
        befCur = p;
        update();
        return;
    }if(drawtype == PaintTool::AffinTrans){
        if(e->button() == Qt::LeftButton)affinmode = 0;//カーソルを動かした方向への平行移動
        if(e->button() == Qt::RightButton) affinmode = 1;//比率を維持した拡縮
        if(e->button() == Qt::MiddleButton)affinmode = 2;//回転
        if(e->button() == Qt::RightButton || e->button() == Qt::MiddleButton){
            if(basePoint == QPointF(-1,-1) && !IsAffinMoved){
                qDebug() << "select second point to determine base axis,  base point is (" << p.x() << ", " << p.y() << ")";
                basePoint = p;
            }
        }
        befCur = p;
    }
    if(e->button() ==Qt::LeftButton){
        if(drawtype == PaintTool::Rectangle)model.back()->drawOutline(p, 0, gridsize, size());
        else if(drawtype == PaintTool::Polyline){
            model.back()->drawOutline(p, 1, gridsize, size());
        }
        else if(drawtype == PaintTool::Polygon)model.back()->drawOutline(p, 2, gridsize, this->size());
        else if(drawtype == PaintTool::EditVertex)model.back()->editOutlineVertex(p, gridsize, this->size(), 0);
        else if(drawtype == PaintTool::NewGradationMode)AddPoints4Gradation(e, p);
        else if(drawtype == PaintTool::Crease) model.back()->ControlPoint(drawtype, 0, p_ongrid, curveDimention, 0 ,model.back()->refCreases);
        else if(drawtype == PaintTool::DeleteCurve){
            model.back()->SelectCurve(p);
            std::vector<int> n(1);
            n[0] = model.back()->DeleteCurve();
            emit deleteCrvSignal(n);
            MoveCrvIndex = {-1, -1};
        }
        if(drawtype == PaintTool::Bspline || drawtype == PaintTool::Line || drawtype == PaintTool::Arc){
            bool res = model.back()->ControlPoint(drawtype, 0, p_ongrid, curveDimention, DivSize, model.back()->refCreases);
            if(res)model.back()->addRulings();
        }
        else if(drawtype == PaintTool::InsertCtrlPt){
            if(MoveCrvIndex[0] == -1){
                model.back()->SelectCurve(p);
                MoveCrvIndex[0] = model.back()->IsSelectedCurve();
            }else{
                if(model.back()->RulingCurve[MoveCrvIndex[0]]->getCurveType() == CurveType::bsp3 && !model.back()->RulingCurve[MoveCrvIndex[0]]->isempty){
                    model.back()->RulingCurve[MoveCrvIndex[0]]->InsertControlPoint(p_ongrid);
                    model.back()->RulingCurve[MoveCrvIndex[0]]->SetNewPoint();
                    model.back()->RulingCurve[MoveCrvIndex[0]]->drawBspline(curveDimention,crvPtNum);
                    if(model.back()->outline->IsClosed()){
                        model.back()->RulingCurve[MoveCrvIndex[0]]->BsplineRulings(model.back()->outline,DivSize,crvPtNum, curveDimention);
                        if(model.back()->RulingCurve[MoveCrvIndex[0]]->Rulings.front()->v != nullptr){
                            model.back()->addRulings();
                            model.back()->deform();
                            //emit foldingSignals();
                        }
                    }
                }
            }
        }
        else if(drawtype == PaintTool::DeleteCtrlPt){
            model.back()->SelectCurve(p);
            model.back()->ControlPoint(drawtype, 1, p, curveDimention, DivSize, model.back()->refCreases);
            model.back()->AssignRuling(3);

        }else if(drawtype == PaintTool::MoveCtrlPt){
            MoveCrvIndex = model.back()->searchPointIndex(p, movePt, 0);

        }

    }else if(e->button() == Qt::RightButton){
        if(drawtype == PaintTool::Rectangle || drawtype == PaintTool::Polyline || drawtype == PaintTool::Polygon) model.back()->outline->eraseVertex();
        if(drawtype == PaintTool::Bspline || drawtype == PaintTool::Line || drawtype == PaintTool::Arc){
            if(MoveCrvIndex[0] != -1) model.back()->RulingCurve[MoveCrvIndex[0]]->eraseCtrlPt(curveDimention, crvPtNum);
        }
        else if(drawtype == PaintTool::Crease) model.back()->ControlPoint(drawtype, 1 , p_ongrid, curveDimention, DivSize, model.back()->refCreases);

    }
    if(model.back()->outline->IsClosed()){
        emit foldingSignals();
    }
    update();
}

void GLWidget_2D::mouseMoveEvent(QMouseEvent *e){
    QPointF p = this->mapFromGlobal(QCursor::pos());

    if(drawtype == PaintTool::None) {
        if(IsLeftClicked){
            difCursol_x += 0.5*(p.x() - befCur.x()); difCursol_y += 0.5*(p.y() - befCur.y());
        }
        update();
        befCur = p; return;
    }

    if(drawtype == PaintTool::AffinTrans){
        IsAffinMoved = true;
        model.back()->AffinTranse_Crease(affinmode, befCur, p, basePoint);

        befCur = p;
    }

    QPointF p_ongrid = SetOnGrid(p, gridsize, this->size());
    if(drawtype == PaintTool::Polygon){
        if(model.back()->outline->hasPtNum == 1){
            model.back()->drawOutline(p, 2, gridsize, size(), false);
            model.back()->outline->hasPtNum = 1;
        }
    }

    if(drawtype == PaintTool::EditVertex){
        model.back()->editOutlineVertex(p, gridsize, size(), 1);
        model.back()->addRulings();
        model.back()->deform();
        emit foldingSignals();
    }
    else if(drawtype == PaintTool::NewGradationMode){}

    if(drawtype == PaintTool::MoveCtrlPt){
        auto res = model.back()->MoveCurvePoint(p_ongrid,MoveCrvIndex[0], movePt, curveDimention, DivSize);
        if(res){
            model.back()->AssignRuling(3);
            update();
            emit foldingSignals();
        }
    }

    if(MoveCrvIndex[0] != -1){
        if(drawtype == PaintTool::InsertCtrlPt &&  model.back()->RulingCurve[MoveCrvIndex[0]]->getCurveType() == CurveType::bsp3 && !model.back()->RulingCurve[MoveCrvIndex[0]]->isempty){//制御点の挿入(B-spline)
            model.back()->RulingCurve[MoveCrvIndex[0]]->InsertControlPoint(p_ongrid);
        }
        
        if(drawtype != PaintTool::NewGradationMode && model.back()->outline->IsClosed() && !model.back()->RulingCurve[MoveCrvIndex[0]]->isempty){

            if(model.back()->RulingCurve[MoveCrvIndex[0]]->getCurveType() == CurveType::bsp3) model.back()->RulingCurve[MoveCrvIndex[0]]->BsplineRulings(model.back()->outline,DivSize,crvPtNum, curveDimention);
            if(model.back()->RulingCurve[MoveCrvIndex[0]]->getCurveType() == CurveType::line) model.back()->RulingCurve[MoveCrvIndex[0]]->LineRulings(model.back()->outline,DivSize);
            if(model.back()->RulingCurve[MoveCrvIndex[0]]->getCurveType() == CurveType::arc) model.back()->RulingCurve[MoveCrvIndex[0]]->ArcRulings(model.back()->outline,DivSize);
            model.back()->addRulings();
            model.back()->deform();
        }
    }

    if(model.back()->CrossDection4AllCurve()) emit foldingSignals();
    update();

}

void GLWidget_2D::mouseReleaseEvent(QMouseEvent * e){
    movePt = -1; affinmode = -1;
    if(IsAffinMoved){basePoint = QPointF(-1,-1); IsAffinMoved = false;}
    QPointF p = mapFromGlobal(QCursor::pos());
    if(drawtype == PaintTool::None){
        IsLeftClicked = IsRightClicked = false;
    }


    if(MoveCrvIndex[0] != -1 && !model.back()->RulingCurve.empty()) model.back()->RulingCurve[MoveCrvIndex[0]]->InsertPointSegment = -1;
    if(drawtype == PaintTool::EditVertex){
        model.back()->editOutlineVertex(p, gridsize, size(), 2);
        if(model.back()->outline->IsClosed())model.back()->deform();
    }
    update();
}

void GLWidget_2D::Reset(){
    model.clear();
    model.push_back(std::make_shared<Model>(crvPtNum));
    emit SendNewActiveCheckBox(PaintTool::Reset);
    MoveCrvIndex = {-1, -1};
    emit foldingSignals();
    update();
}

void GLWidget_2D::wheelEvent(QWheelEvent *we){
    DiffWheel = (we->angleDelta().y() > 0) ? 1 : -1;
    if(drawtype == PaintTool::NewGradationMode){
        if(refLine == nullptr || std::find(model.back()->Rulings.begin(), model.back()->Rulings.end(), refLine) == model.back()->Rulings.end()){
            return;
        }
        model.back()->SetGradationValue(DiffWheel, refLine);
        emit ColorChangeFrom(0, refLine->color);
        model.back()->deform();

    }
    if(drawtype == PaintTool::None){
        camscale += DiffWheel;
        camscale = (camscale < -100)?-100.0: (camscale > 200)? 200: camscale;
        qDebug() << "camscale  " << camscale;
    }
    emit foldingSignals();
    update();
}

void GLWidget_2D::AddPoints4Gradation(QMouseEvent *e, QPointF& p){
    Eigen::Vector3d curPos{p.x(), p.y(), 0};
    double dist = 10;
    if(drawtype == PaintTool::NewGradationMode){
        refLine = nullptr;
        for(auto& r: model.back()->Rulings){

            if(r->et == EdgeType::ol || r->et == EdgeType::cl)continue;
            double d = ((curPos - r->o->p).cross(r->o->p - r->v->p)).norm()/(r->o->p - r->v->p).norm();
            if(d < dist){
                dist = d; refLine = r;
            }
        }
        if(refLine == nullptr || std::find(model.back()->Rulings.begin(), model.back()->Rulings.end(), refLine) == model.back()->Rulings.end()){qDebug()<<"not found" ; return;}
        model.back()->SetGradationValue(0, refLine);
        emit ColorChangeFrom(0, refLine->color);
        model.back()->deform();

    }
    emit foldingSignals();
    update();
}

void GLWidget_2D::DrawGradationMode(){
    if(refLine == nullptr)return;
    emit foldingSignals();
    update();
}

void GLWidget_2D::GetGradationFromSlider(int val){
    if(refLine == nullptr) return;
    refLine->color = val;
    emit ColorChangeFrom(1, val);
    model.back()->SetGradationValue(0, refLine);
    emit foldingSignals();
    update();
}

