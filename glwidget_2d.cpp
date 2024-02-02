#include "glwidget_2d.h"
#include <GL/glu.h>

GLWidget_2D::GLWidget_2D(QWidget *parent):QOpenGLWidget(parent)
{

    CurveList = {{CurveType::bezier3, PaintTool::Bezier_r}, {CurveType::bsp3, PaintTool::Bspline_r}, {CurveType::line, PaintTool::Line_r}, {CurveType::arc, PaintTool::Arc_r}};

    curveDimention = 3;
    drawtype = PaintTool::None;
    crvPtNum = 3000;
    maxDivSize = 3000;
    standardDist = 2.5;
    gw = std::make_shared<GToolWnd>(this);
    gw->hide();

    cval = 128;
    ctype = 1;
    connect(this , &GLWidget_2D::CurvePathSet, gw.get(), &GToolWnd::set);

    InterpolationType = 0;
    //ControllPoints_gradation.clear();
    DiffWheel = 0;
    movePt = -1;
    MoveCrvIndex = {-1, -1};
    refV = std::shared_ptr<Vertex>(nullptr);
    refL = std::shared_ptr<Line>(nullptr);
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
    RegressionCurve.clear();
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
    if(model.back()->refFL == nullptr || std::find(model.back()->FL.begin(), model.back()->FL.end(), model.back()->refFL) == model.back()->FL.end())IsCopied = false;
    IsCopied = true;
}

void GLWidget_2D::PasteCurveObj(){
    if(!IsCopied || std::find(model.back()->FL.begin(), model.back()->FL.end(), model.back()->refFL) == model.back()->FL.end())return;
    auto newfl = model.back()->refFL->deepCopy();
    //上方向に10移動させる
    double y = -30;
    for(auto&p: newfl->CtrlPts)p.y() += y;
    newfl->setCurve(3);
    newfl->FoldingCurve.clear(); newfl->a_flap = -1;
    model.back()->FL.push_back(newfl);
}

void GLWidget_2D::AddCurve(){
    std::vector<int>deleteIndex;
    CurveType _type;
    emit signalCurveType(_type);
    model.back()->Check4Param(curveDimention, deleteIndex);
    if(model.back()->crvs.empty()) MoveCrvIndex = {-1,-1};
    emit deleteCrvSignal(deleteIndex);
    for(auto& c: CurveList){
        if(std::get<0>(c) == _type){
            drawtype = std::get<1>(c);
            curvetype = _type;
        }
    }

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
    if(model.back()->crvs.empty()) MoveCrvIndex[0] = -1;
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
    model.back()->Check4Param(curveDimention, deleteIndex);
    emit deleteCrvSignal(deleteIndex);
    MoveCrvIndex = {-1, -1};
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::MoveCtrlPt);
}

void GLWidget_2D::DrawOutlineRectangle(){
    //drawtype = -1;
    drawtype = PaintTool::Rectangle_ol;
    std::vector<int>deleteIndex;
    model.back()->Check4Param(curveDimention, deleteIndex);
    if(model.back()->crvs.empty()) MoveCrvIndex[0] = -1;
    if(model.back()->FL.empty())MoveCrvIndex[1] = -1;
    emit deleteCrvSignal(deleteIndex);
    model.back()->outline = std::make_shared<OUTLINE>();
    model.back()->outline->type = "Rectangle";
    this->setMouseTracking(false);
    emit SendNewActiveCheckBox(PaintTool::Rectangle_ol);
}

void GLWidget_2D::DrawOutlinePolygon(int state){
    if(state == 0)return;
    drawtype = PaintTool::Polygon_ol;
    std::vector<int>deleteIndex;
    model.back()->Check4Param(curveDimention, deleteIndex);
    if(model.back()->crvs.empty()) MoveCrvIndex[0] = -1;
    if(model.back()->FL.empty()) MoveCrvIndex[1] = -1;
    emit deleteCrvSignal(deleteIndex);
    model.back()->outline = std::make_shared<OUTLINE>();
    model.back()->outline->type = "Polygon";
    for(auto&c: model.back()->crvs)c->isempty = true;
    emit getEdgeNum();
    emit SendNewActiveCheckBox(PaintTool::Polygon_ol);
    this->setMouseTracking(true);
}

void GLWidget_2D::DrawOutlinePolyline(int state){
    if(state == 0)return;
    drawtype = PaintTool::Polyline_ol;
    std::vector<int>deleteIndex;
    model.back()->Check4Param(curveDimention, deleteIndex);
    if(model.back()->crvs.empty()) MoveCrvIndex[0] = -1;
    if(model.back()->FL.empty()) MoveCrvIndex[1] = -1;
    emit deleteCrvSignal(deleteIndex);
    model.back()->outline = std::make_shared<OUTLINE>();
    model.back()->outline->type = "Polyline";
    for(auto&c: model.back()->crvs)c->isempty = true;
    emit SendNewActiveCheckBox(PaintTool::Polyline_ol);
}

void GLWidget_2D::MoveOutline(int state){
    if(state == 0)return; drawtype = PaintTool::Move_ol;
    std::vector<int>deleteIndex;
    model.back()->Check4Param(curveDimention, deleteIndex);
    if(model.back()->crvs.empty()) MoveCrvIndex[0] = -1;
    if(model.back()->FL.empty()) MoveCrvIndex[1] = -1;
    emit deleteCrvSignal(deleteIndex);
    emit SendNewActiveCheckBox(PaintTool::Move_ol);
}

void GLWidget_2D::EditOutlineVertex(int state){
    drawtype = PaintTool::EditVertex_ol;
    std::vector<int>deleteIndex;
    model.back()->Check4Param(curveDimention, deleteIndex);
    if(model.back()->crvs.empty()) MoveCrvIndex[0] = -1;
    if(model.back()->FL.empty()) MoveCrvIndex[1] = -1;
    emit deleteCrvSignal(deleteIndex);
    emit SendNewActiveCheckBox(PaintTool::EditVertex_ol);
    update();
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
    if(model.back()->crvs.empty()) MoveCrvIndex[0] = -1;
    if(model.back()->FL.empty()) MoveCrvIndex[1] = -1;
    emit deleteCrvSignal(deleteIndex);
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::DeleteCurve);

}

void GLWidget_2D::changeFoldType(PaintTool state){
    drawtype = state;
    IsMVcolor_binary = true;
    model.back()->RemoveUnable2GenCurve();

    if(state != PaintTool::FoldlineColor){
        std::shared_ptr<FoldLine> fl = std::make_shared<FoldLine>(state);
        model.back()->FL.push_back(fl);
        model.back()->ChangeFoldLineState();
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
    model.back()->Check4Param(curveDimention, deleteIndex);
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
    if(MoveCrvIndex[0] == -1 || model.back()->crvs.empty()) return;
    bool res = false;
    if(curvetype == CurveType::bezier3){
       res = model.back()->crvs[MoveCrvIndex[0]]->drawBezier(curveDimention, crvPtNum);
        //if(model.back()->outline->isClosed  && model.back()->crv->CurvePoints[0].pt != glm::f64vec2{-1,-1})model.back()->crv->BezierRulings(model.back()->outline->vertices,DivSize,crvPtNum);
    }
    if(curvetype == CurveType::bsp3){
        res = model.back()->crvs[MoveCrvIndex[0]]->drawBspline(curveDimention, crvPtNum);
        if(model.back()->outline->IsClosed() && res)model.back()->crvs[MoveCrvIndex[0]]->BsplineRulings(model.back()->outline,DivSize,crvPtNum, curveDimention);
    }
    if(curvetype == CurveType::line){
        res = model.back()->crvs[MoveCrvIndex[0]]->drawLine();
        if(model.back()->outline->IsClosed() && res)model.back()->crvs[MoveCrvIndex[0]]->LineRulings(model.back()->outline,DivSize);
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
    std::iter_swap(model.back()->crvs.begin() + n1, model.back()->crvs.begin() + n2);
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
    //if(model.back()->Faces.size() < 2 )return;
    if(model.back()->FL.empty())return;
    model.back()->FL[0]->drawRulingInAllAngles(AllRulings);
    update();
    emit foldingSignals();
}


void GLWidget_2D::changeflapgnle(double val, bool begin_center){
    //if(model.back()->Faces.size() < 2 || IsStop4Debug || model.back()->FL.empty())return;
    if(model.back()->FL.empty())return;
    if(model.back()->FL.empty())model.back()->FL.push_back(std::make_shared<FoldLine>(PaintTool::FoldLine_test) );
    model.back()->applyAAAMethod(val, begin_center);
    //model.back()->FL[0]->test_rotate(val);
    emit foldingSignals();
    update();
}

void GLWidget_2D::ReceiveRegressionCurve(const std::vector<std::vector<std::vector<std::shared_ptr<Vertex>>>>& RC, const std::vector<std::vector<double>>color){
    RegressionCurve.clear();
    for(int i = 0; i < (int)RC.size(); i++)RegressionCurve.push_back(drawobj(color[i], RC[i]));
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
            for(auto&fl: model.back()->FL){
                glColor3d(0,0.3,0.3);
                glPointSize(5);
                for(auto&v: fl->CtrlPts){
                    glBegin(GL_POINTS);
                    glVertex3d(v.x(),v.y(), 0);
                    glEnd();
                }

                glEnable(GL_LINE_STIPPLE);
                glLineStipple(1 , 0xF0F0);
                glBegin(GL_LINE_STRIP);
                glColor3d(0,0,0);
                for(auto&v: fl->CtrlPts)glVertex3d(v.x(), v.y(), 0);
                glEnd();
                glDisable(GL_LINE_STIPPLE);

                if(fl->CurvePts.empty())continue;
                glBegin(GL_LINE_STRIP);
                for(auto&v: fl->CurvePts) glVertex3d(v.x(),v.y(), 0);
                glEnd();

                glColor3d(0,1,0);
                glPointSize(5);
                for(auto&h: fl->FoldingCurve){
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

    //regression curve
    for(const auto&RC: RegressionCurve){
        glColor3d(RC.c[0], RC.c[1], RC.c[2]);
        glPointSize(4);
        for(auto&Vertices: RC.V){
            for(auto&v: Vertices){
                glBegin(GL_POINTS);
                glVertex3d(v->p.x(), v->p.y(),0);
                glEnd();
            }
        }
        glColor3d(0,0,0);
        glEnable(GL_LINE_STIPPLE);
        glLineStipple(1 , 0x00F0);
        for(const auto&Vertices: RC.V){
            glBegin(GL_LINE_LOOP);
            for(auto&v: Vertices)glVertex3d(v->p.x(), v->p.y(),0);
            glEnd();
        }
        glDisable(GL_LINE_STIPPLE);

        glColor4d(1,0,0, 0.4);
        glBegin(GL_LINE_STRIP);
        for(const auto&V: RC.V) glVertex3d(V.back()->p.x(), V.back()->p.y(),0);
        glEnd();
    }

    glPolygonOffset(0.0,1.0);
    auto getcolor = [](double c, double a, double y)->double{
        if(y < a)return c/a * y/255.0;
        return ((255.0 - c)*(y - a)/(std::numbers::pi - a) + c)/255.0;
    };
    //折曲線上の4価頂点の描画
    for(auto& FL: model.back()->FL){
        if(FL->FoldingCurve.empty())continue;
        std::vector<int> Vertices_Ind;
        for(int i = 0; i < (int)FL->FoldingCurve.size(); i++){if(FL->FoldingCurve[i]->IsCalc)Vertices_Ind.push_back(i);}
        for(const auto& fc: FL->FoldingCurve){
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
        bool IsGradation = (drawtype == PaintTool::NewGradationMode || drawtype == PaintTool::FoldLine_bezier)?true: false;

        for(int i = 1; i < (int)FL->FoldingCurve.size()-1; i++)
            DrawEdge(FL->FoldingCurve[i]->first, FL->FoldingCurve[i-1]->first, FL->FoldingCurve[i]->second, FL->FoldingCurve[i]->third, rulingWidth, IsGradation, false, false, false);
        DrawEdge(FL->FoldingCurve.end()[-2]->first, FL->FoldingCurve.back()->first, FL->FoldingCurve.end()[-2]->second, FL->FoldingCurve.end()[-2]->third, rulingWidth, IsGradation, true, true, false);
        for(int i = 1; i < (int)FL->FoldingCurve.size() - 1; i++){
            bool DashedLine = (FL->FoldingCurve[i]->IsCalc)? false: true;
            DrawEdge(FL->FoldingCurve[i]->first, FL->FoldingCurve[i]->second, FL->FoldingCurve[i+1]->first, FL->FoldingCurve[i-1]->first, rulingWidth, IsGradation, true, false, DashedLine);
            DrawEdge(FL->FoldingCurve[i]->first, FL->FoldingCurve[i]->third, FL->FoldingCurve[i-1]->first, FL->FoldingCurve[i+1]->first, rulingWidth, IsGradation, true, false, DashedLine);
        }



        glColor3d(0,0,0);
    }
    if(refV != nullptr){
        glBegin(GL_POINT);
        glPointSize(6.f);
        glColor3d(1,0,0);
        glVertex3d(refV->p.x(), refV->p.y(),0);
        glEnd();
    }

    //折曲線との交差がないrulingの描画
    for(auto itr_r = model.back()->Rulings.begin(); itr_r != model.back()->Rulings.end(); itr_r++){
        if((*itr_r)->hasCrossPoint)continue;
        glColor3d(0,0,0);
        glLineWidth(1);
        if((*itr_r)->et == EdgeType::r){
            if(drawtype == PaintTool::NewGradationMode|| drawtype == PaintTool::FoldLine_bezier){
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
    if(!visibleCurve || drawtype == PaintTool::FoldLine_bezier)return;
    //ruling制御曲線の制御点
    if(drawtype != PaintTool::NewGradationMode){
        for(int i = 0; i < (int)model.back()->crvs.size(); i++){
            if((i == model.back()->getSelectedCurveIndex(mapFromGlobal(QCursor::pos()))[0]) ||
                    ((drawtype == PaintTool::Bezier_r || drawtype == PaintTool::Bspline_r || drawtype == PaintTool::Line_r || drawtype == PaintTool::Arc_r)
                     && MoveCrvIndex[0] == i)
                    || drawtype == PaintTool::MoveCtrlPt || drawtype == PaintTool::InsertCtrlPt ||  drawtype == PaintTool::DeleteCtrlPt || drawtype == PaintTool::None){
                //glPolygonOffset(0.f,1.f);
                for (auto& c: model.back()->crvs[i]->ControllPoints) {
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
                for (auto& c: model.back()->crvs[i]->ControllPoints)glVertex3d( c.x(), c.y(),0);
                glEnd();
                glDisable(GL_LINE_STIPPLE);
                glPolygonOffset(0.f,0.f);
            }
        }

        //ruling制御曲線の描画
        for(int i = 0; i < (int)model.back()->crvs.size(); i++){
            if(i == MoveCrvIndex[0])glColor3d(1, 0, 0);
            else glColor3d(0.4, 0.4, 0.4);
            glBegin(GL_LINE_STRIP);
            for(auto &c: model.back()->crvs[i]->CurvePoints){
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
    Eigen::Vector3d p_ongrid = SetOnGrid(p, gridsize);
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
        if(drawtype == PaintTool::Rectangle_ol)model.back()->drawOutline(p, 0, gridsize, size());
        else if(drawtype == PaintTool::Polyline_ol){
            model.back()->drawOutline(p, 1, gridsize, size());
        }
        else if(drawtype == PaintTool::Polygon_ol)model.back()->drawOutline(p, 2, gridsize, size());
        else if(drawtype == PaintTool::EditVertex_ol)model.back()->editOutlineVertex(p, gridsize, size(), 0);

        else if(drawtype == PaintTool::Move_ol) model.back()->outline->MoveOutline(p_ongrid);
        else if(drawtype == PaintTool::Const_ol) model.back()->addConstraint(p, 0, gridsize, model.back()->Axis4Const, size());
        else if(drawtype ==PaintTool::ConnectVertices_ol)model.back()->ConnectOutline(p, gridsize, size());
        else if(drawtype == PaintTool::NewGradationMode || drawtype ==PaintTool::FoldlineColor)addPoints_intplation(e, p);
        else if(drawtype == PaintTool::FoldLine_bezier || drawtype == PaintTool::FoldLine_arc || drawtype == PaintTool::FoldLine_line ){
            model.back()->AddControlPoint_FL(p_ongrid, 0, curveDimention);
            update();
        }
        else if(drawtype == PaintTool::FoldLine_test){
            qDebug() << "can't use test now";
        }
        else if(drawtype == PaintTool::DeleteCurve){
            model.back()->SelectCurve(p);
            std::vector<int> n(1);
            n[0] = model.back()->DeleteCurve();
            emit deleteCrvSignal(n);
            MoveCrvIndex = {-1, -1};
        }
        if(drawtype == PaintTool::Bezier_r || drawtype == PaintTool::Bspline_r || drawtype == PaintTool::Line_r || drawtype == PaintTool::Arc_r){
            bool res = model.back()->AddControlPoint(p_ongrid, curveDimention, DivSize);
            if(res)model.back()->addRulings();
        }
        else if(drawtype == PaintTool::InsertCtrlPt){
            if(MoveCrvIndex[0] == -1){
                model.back()->SelectCurve(p);
                MoveCrvIndex[0] = model.back()->IsSelectedCurve();
            }else{
                if(model.back()->crvs[MoveCrvIndex[0]]->getCurveType() == CurveType::bsp3 && !model.back()->crvs[MoveCrvIndex[0]]->isempty){
                    model.back()->crvs[MoveCrvIndex[0]]->InsertControlPoint2(p_ongrid);
                    model.back()->crvs[MoveCrvIndex[0]]->SetNewPoint();
                    model.back()->crvs[MoveCrvIndex[0]]->drawBspline(curveDimention,crvPtNum);
                    if(model.back()->outline->IsClosed()){
                        model.back()->crvs[MoveCrvIndex[0]]->BsplineRulings(model.back()->outline,DivSize,crvPtNum, curveDimention);
                        if(model.back()->crvs[MoveCrvIndex[0]]->Rulings.front()->v != nullptr){
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
            model.back()->DeleteControlPoint(p, curveDimention, DivSize);
            double tol; bool begincenter;
            emit getFoldParam(tol, begincenter);
            model.back()->AssignRuling(3,tol, begincenter);

        }else if(drawtype == PaintTool::MoveCtrlPt){
            MoveCrvIndex = model.back()->searchPointIndex(p, movePt, 0);

        }else if(drawtype == PaintTool::CheckDevelopability){
            int ind = -1;
            std::shared_ptr<FoldLine> _fl = std::shared_ptr<FoldLine>(nullptr);
            refV = nullptr;
            double dist = 3.0;
            if(model.back()->FL.empty())return;
            for(auto&fl: model.back()->FL){
                for(auto itr = fl->FoldingCurve.begin() + 1; itr != fl->FoldingCurve.end() - 1; itr++){
                    if((Eigen::Vector3d(p.x(), p.y(), 0) - (*itr)->first->p).norm() < dist){
                        ind = std::distance(fl->FoldingCurve.begin(), itr); dist = (Eigen::Vector3d(p.x(), p.y(), 0) - (*itr)->first->p).norm();
                        _fl = fl;
                    }
                }
            }
            if(ind != -1){
                auto v4d = _fl->FoldingCurve[ind];
                Eigen::Vector3d et = (v4d->second->p3 -v4d->first->p3).normalized(), er = (_fl->FoldingCurve[ind-1]->first->p3 -v4d->first->p3).normalized(),
                        eb = (v4d->third->p3 -v4d->first->p3).normalized(), el = (_fl->FoldingCurve[ind+1]->first->p3 -v4d->first->p3).normalized();
                double phi1 = std::acos(et.dot(er)), phi2 = std::acos(et.dot(el)), phi3 = std::acos(eb.dot(el)), phi4 = std::acos(eb.dot(er));
                refV = _fl->FoldingCurve[ind]->first;
                if(DebugMode::Singleton::getInstance().isdebug())
                qDebug() << "developability  = " <<  abs(2.0*std::numbers::pi - phi1 - phi2 - phi3 - phi4) << ", phi1 = " << MathTool::rad2deg(phi1) << " , phi2 = " << MathTool::rad2deg(phi2) << ", phi3 = " << MathTool::rad2deg(phi3) << ", phi4 = " << MathTool::rad2deg(phi4);
            }

        }

    }else if(e->button() == Qt::RightButton){
        if(drawtype == PaintTool::Rectangle_ol || drawtype == PaintTool::Polyline_ol || drawtype == PaintTool::Polygon_ol) model.back()->outline->eraseVertex();
        if(drawtype == PaintTool::Bezier_r || drawtype == PaintTool::Bspline_r || drawtype == PaintTool::Line_r || drawtype == PaintTool::Arc_r){
            if(MoveCrvIndex[0] != -1) model.back()->crvs[MoveCrvIndex[0]]->eraseCtrlPt(curveDimention, crvPtNum);
        }
        else if(drawtype == PaintTool::FoldLine_line || drawtype == PaintTool::FoldLine_arc || drawtype == PaintTool::FoldLine_bezier){
            model.back()->AddControlPoint_FL(p_ongrid, 1, curveDimention);
        }

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

    Eigen::Vector3d p_ongrid = SetOnGrid(p, gridsize);
    if(drawtype == PaintTool::Polygon_ol){
        if(model.back()->outline->hasPtNum == 1){
            model.back()->drawOutline(p, 2, gridsize, size(), false);
            model.back()->outline->hasPtNum = 1;
        }
    }
    if(drawtype == PaintTool::Move_ol) model.back()->outline->MoveOutline(p_ongrid);
    else if(drawtype == PaintTool::EditVertex_ol){
        model.back()->editOutlineVertex(p, gridsize, size(), 1);
        model.back()->addRulings();
        model.back()->deform();
        emit foldingSignals();
    }
    else if(drawtype == PaintTool::NewGradationMode){}

    if(drawtype == PaintTool::MoveCtrlPt){
        auto res = model.back()->MoveCurvePoint(p_ongrid,MoveCrvIndex[0], movePt, curveDimention, DivSize);
        if(res){
            double tol; bool begincenter;
            emit getFoldParam(tol, begincenter);
            model.back()->AssignRuling(3,tol, begincenter);
            update();
            emit foldingSignals();
        }
    }
    else if(drawtype == PaintTool::FoldLine_bezier){
        //bool res = model.back()->updateSplitRulings(model.back()->FL[FoldCurveIndex], curveDimention);
        //if(res) emit foldingSignals();
        update();
    }
    if(MoveCrvIndex[0] != -1){
        if(drawtype == PaintTool::InsertCtrlPt &&  model.back()->crvs[MoveCrvIndex[0]]->getCurveType() == CurveType::bsp3 && !model.back()->crvs[MoveCrvIndex[0]]->isempty){//制御点の挿入(B-spline)
            model.back()->crvs[MoveCrvIndex[0]]->InsertControlPoint2(p_ongrid);
        }
        
        if(drawtype != PaintTool::NewGradationMode && model.back()->outline->IsClosed() && !model.back()->crvs[MoveCrvIndex[0]]->isempty){

            if(model.back()->crvs[MoveCrvIndex[0]]->getCurveType() == CurveType::bsp3) model.back()->crvs[MoveCrvIndex[0]]->BsplineRulings(model.back()->outline,DivSize,crvPtNum, curveDimention);
            if(model.back()->crvs[MoveCrvIndex[0]]->getCurveType() == CurveType::line) model.back()->crvs[MoveCrvIndex[0]]->LineRulings(model.back()->outline,DivSize);
            if(model.back()->crvs[MoveCrvIndex[0]]->getCurveType() == CurveType::arc) model.back()->crvs[MoveCrvIndex[0]]->ArcRulings(model.back()->outline,DivSize);
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


    if(MoveCrvIndex[0] != -1 && !model.back()->crvs.empty()) model.back()->crvs[MoveCrvIndex[0]]->InsertPointSegment = -1;
    if(drawtype == PaintTool::EditVertex_ol){
        model.back()->editOutlineVertex(p, gridsize, size(), 2);
        if(model.back()->outline->IsClosed())model.back()->deform();
    }
    update();
}

void GLWidget_2D::cb_ApplyCurveEvent(){
    if(KeyEvent == 0){
        qDebug() << "finished "; KeyEvent = -1;
    }else if(KeyEvent == -1){
        if(MoveCrvIndex[0] == -1) qDebug() << "no curve is selected";
    }
    update();
}

void GLWidget_2D::cb_DeleteCurve(){
    KeyEvent = 1;
    if(MoveCrvIndex[0] != -1){
        model.back()->DeleteCurve();
    }
    update();
}

void GLWidget_2D::Reset(){
    model.clear();
    model.push_back(std::make_shared<Model>(crvPtNum));
    tmp_c.clear();
    emit SendNewActiveCheckBox(PaintTool::Reset);
    MoveCrvIndex = {-1, -1};
    emit foldingSignals();
    update();
}

void GLWidget_2D::receiveColors(){
    model.back()->deform();
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
    for(auto&c: model.back()->crvs){
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
        //emit ColorChangeFrom(0, model.back()->FL[FoldCurveIndex]->getColor());
    }else if(drawtype == PaintTool::NewGradationMode){
        if(refL == nullptr || std::find(model.back()->Rulings.begin(), model.back()->Rulings.end(), refL) == model.back()->Rulings.end()){
            //qDebug()<<"there is no referenced Ruling";
            return;
        }
        model.back()->setGradationValue(DiffWheel, refL, InterpolationType, CurvePath);
        emit ColorChangeFrom(0, refL->color);
        model.back()->deform();
        if(isVisibleTo(gw.get())) emit CurvePathSet(CurvePath);
    }
    if(drawtype == PaintTool::None){
        camscale += DiffWheel;
        camscale = (camscale < -100)?-100.0: (camscale > 200)? 200: camscale;
        qDebug() << "camscale  " << camscale;
    }
    emit foldingSignals();
    update();
}

void GLWidget_2D::addPoints_intplation(QMouseEvent *e, QPointF& p){
    Eigen::Vector3d curPos{p.x(), p.y(), 0};
    double dist = 10;
    if(drawtype == PaintTool::NewGradationMode){
        refL = nullptr;
        for(auto& r: model.back()->Rulings){

            if(r->et == EdgeType::ol || r->et == EdgeType::cl)continue;
            double d = ((curPos - r->o->p).cross(r->o->p - r->v->p)).norm()/(r->o->p - r->v->p).norm();
            if(d < dist){
                dist = d; refL = r;
            }
        }
        if(refL == nullptr || std::find(model.back()->Rulings.begin(), model.back()->Rulings.end(), refL) == model.back()->Rulings.end()){qDebug()<<"not found" ; return;}
        model.back()->setGradationValue(0, refL, InterpolationType, CurvePath);
        emit ColorChangeFrom(0, refL->color);
        model.back()->deform();
        if(model.back()->outline->IsClosed() && !model.back()->FL.empty()){
            //auto oriedge = model.back()->outline->getEdges();
            //model.back()->FL[0]->modify2DRulings(model.back()->Faces, model.back()->Edges, model.back()->vertices, oriedge, curveDimention);
        }
    }else if(drawtype == PaintTool::FoldlineColor){}
    emit foldingSignals();
    update();
}

void GLWidget_2D::ApplyNewGradationMode(){
    if(refL == nullptr){qDebug() << "you neeed to add gradation point"; return;}
    emit foldingSignals();
    if(isVisibleTo(gw.get())) emit CurvePathSet(CurvePath);
    update();
}

void GLWidget_2D::getGradationFromSlider(int val){
    if(refL == nullptr) return;
    refL->color = val;
    emit ColorChangeFrom(2, val);
    model.back()->setGradationValue(0, refL, InterpolationType, CurvePath);
    emit foldingSignals();
    update();
}

void GLWidget_2D::OpenDebugWindwow(){
    if(!isVisibleTo(gw.get())) gw = std::make_shared<GToolWnd>(this);
    emit CurvePathSet(CurvePath);

    gw->show();
}

int GLWidget_2D::assignment_refL(){

    QPointF p = mapFromGlobal(QCursor::pos());
    Eigen::Vector3d curPos{p.x(), p.y(), 0};
    int ind = -1;
    double dist = 10;
    for(int i = 0; i < (int)model.back()->Rulings.size(); i++){
        std::shared_ptr<Line> r = model.back()->Rulings[i];
        if(r->et != EdgeType::r)continue;
        double d = ((curPos - r->o->p).cross(r->o->p - r->v->p)).norm()/(r->o->p - r->v->p).norm();
        if(d < dist){
            dist = d; ind = i;
        }
    }
    return ind;
}

Eigen::Vector3d GLWidget_2D::SetOnGrid(QPointF& cursol, double gridsize){
    int x = (int)cursol.x() % (int)gridsize, y = (int)cursol.y() % (int)gridsize;
    x = (cursol.x() - x + gridsize/2);
    y = (cursol.y() - y + gridsize/2);
    return Eigen::Vector3d{static_cast<double>(x),static_cast<double>(y),0};
}
