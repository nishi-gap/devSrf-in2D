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
    SelectedCurveIndex = -1;
    refHE = nullptr;
    KeyEvent = -1;
    curvetype = CurveType::none;

    gridsize = 10;
    visibleGrid = 1;

    eraseVec2d = false;
    visibleCurve = true;
    model = new Model(crvPtNum);

}
GLWidget_2D::~GLWidget_2D(){}

void GLWidget_2D::InitializeDrawMode(int state){
    if(state == 0)return; drawtype = PaintTool::None; SelectedCurveIndex = -1; emit SendNewActiveCheckBox(PaintTool::None);
}

void GLWidget_2D::AddCurve(){
    std::vector<int>deleteIndex;
    CurveType _type;
    emit signalCurveType(_type);
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    for(auto& c: CurveList){
        if(std::get<0>(c) == _type){
            drawtype = std::get<1>(c);
            curvetype = _type;
        }
    }

    SelectedCurveIndex = model->AddNewCurve(curvetype, DivSize);
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::AddCurve);
}

void GLWidget_2D::InsertNewPoint(){
    drawtype = PaintTool::InsertCtrlPt;
    SelectedCurveIndex = -1;
    setMouseTracking(true);
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    emit SendNewActiveCheckBox(PaintTool::InsertCtrlPt);
}

void GLWidget_2D::MoveCurvePt(){
    //curveNum = 1;
    drawtype = PaintTool::MoveCtrlPt;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    SelectedCurveIndex = -1;
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::MoveCtrlPt);
}

void GLWidget_2D::DrawOutlineRectangle(){
    //drawtype = -1;
    drawtype = PaintTool::Rectangle_ol;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
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
    if(model->crvs.empty()) SelectedCurveIndex = -1;
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
    if(model->crvs.empty()) SelectedCurveIndex = -1;
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
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    emit SendNewActiveCheckBox(PaintTool::Move_ol);
}

void GLWidget_2D::EditOutlineVertex(int state){
    drawtype = PaintTool::EditVertex_ol;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    emit SendNewActiveCheckBox(PaintTool::EditVertex_ol);
    update();
}

void GLWidget_2D::DeleteCtrlPt(){
    drawtype = PaintTool::DeleteCtrlPt;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    emit deleteCrvSignal(deleteIndex);
    SelectedCurveIndex = -1;
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::DeleteCtrlPt);
}

void GLWidget_2D::DeleteCurve(){
    drawtype = PaintTool::DeleteCurve;
    std::vector<int>deleteIndex;
    model->Check4Param(curveDimention, deleteIndex);
    if(model->crvs.empty()) SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    setMouseTracking(true);
    emit SendNewActiveCheckBox(PaintTool::DeleteCurve);

}

void GLWidget_2D::changeFoldType(PaintTool state){
    drawtype = state;
    if(state != PaintTool::FoldlineColor){
        int rulingNum = 30;
        int crvNum = 1000;
        FoldLine *fl = new FoldLine(crvNum, rulingNum, state);
        model->FL.push_back(fl);
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
    model->Check4Param(curveDimention, deleteIndex); SelectedCurveIndex = -1;
    emit deleteCrvSignal(deleteIndex);
    update();
}
void GLWidget_2D::recieveNewEdgeNum(int num){model->outline->VerticesNum = num; update();}

void GLWidget_2D::ChangedDivSizeEdit(int n){
    this->DivSize = (n < 0)? DivSize: (maxDivSize < n)? maxDivSize: n;
    if(SelectedCurveIndex == -1 || model->crvs.empty()) return;
    if(curvetype == CurveType::bezier3){
        model->crvs[SelectedCurveIndex]->Bezier(curveDimention, crvPtNum);
        //if(model->outline->isClosed  && model->crv->CurvePoints[0].pt != glm::f64vec2{-1,-1})model->crv->BezierRulings(model->outline->vertices,DivSize,crvPtNum);
    }
    if(curvetype == CurveType::bsp3){
        model->crvs[SelectedCurveIndex]->Bspline(curveDimention, crvPtNum);
        if(model->outline->IsClosed() && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1, -1})model->crvs[SelectedCurveIndex]->BsplineRulings(model->outline,DivSize,crvPtNum, curveDimention);
    }
    if(curvetype == CurveType::line){
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

void GLWidget_2D::changeSelectedCurve(int ind){SelectedCurveIndex = ind; drawtype = PaintTool::None; update();}

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

void GLWidget_2D::receiveNewLineWidth(double d){
    rulingWidth = d;
    update();
}

void GLWidget_2D::changeStopAtFF(bool state){IsStopAtFF = state;}
void GLWidget_2D::changeStopAtCon(bool state){IsStopAtCon = state;}
void GLWidget_2D::changeStopAtEq(bool state){IsStopAtEq = state;}
void GLWidget_2D::Start4Debug_CF(){
    if(model->Faces.size() < 2)return;

    auto IsSame = [](double c, double b){return (abs(c - b) <= 1e-6)? true: false;};
    std::string file = "result_AAADebug.csv"; std::ofstream res(file);
    bool hasFoldLine = false;
    std::vector<HalfEdge*> _axis;
    for(const auto& h: model->Edges){
        if(h->edgetype == EdgeType::r){
            hasFoldLine = true;
            if(std::find(_axis.begin(), _axis.end(), h) == _axis.end() && std::find(_axis.begin(), _axis.end(), h->pair) == _axis.end()){
                _axis.push_back(h);
            }
        }
    }
    if(!hasFoldLine)return;
    model->deform();

    angle = 0.0;
    NewRuling.clear();NewRuling2d.clear(); RulingColor.clear();
    while(angle <= 2.0*std::numbers::pi){
        double a2 = 0;
        double a = angle;
        glm::f64vec3 e, e2, x;
        glm::f64vec3 N, r, n, r2;
        glm::f64vec3 e_bef;
        double vecLen = 50;

        for(int i = 0; i < (int)model->Faces.size() - 1; i++){
            x = _axis[i]->vertex->p3;
            e = glm::normalize(_axis[i]->prev->vertex->p3 - x);
            e2 = glm::normalize(_axis[i]->pair->next->next->vertex->p3 - x);
            double beta = (1.0 - abs(_axis[i]->color/255.0)) * std::numbers::pi;
            if(i != 0){
                glm::f64vec3 e_next = glm::normalize(MathTool::ProjectionVector(e2,-e));
                double tau = glm::angle(e_bef, e_next);
                if(glm::dot(glm::cross(e_bef, e_next), e) > 0) tau = 2.0*std::numbers::pi - tau;
                a = (a2 - tau < 0) ? a2 - tau + 2.0 * std::numbers::pi: a2 - tau;
                //std::cout << i << " : tau = " << tau << " , a = " <<  glm::degrees(a) << " , a' = " << glm::degrees(a2) << std::endl;

            }

            double phi3 = glm::angle(e2, glm::normalize(_axis[i]->next->vertex->p3 - x));
            double phi4 = glm::angle(e, glm::normalize(_axis[i]->next->vertex->p3 - x));
            double k = 2.0 * std::numbers::pi - phi3 - phi4;
            if(!(beta <= phi3 + phi4 && phi3 + phi4 <= 2.0 * std::numbers::pi - beta)){std::cout <<"exceed boundary condition   " << beta << " , " << phi3+phi4<< std::endl; exit(0);}
            double phi1 = IsSame(sin(beta)*cos(a), sin(k))? std::numbers::pi/2.0: atan((cos(k) - cos(beta))/(sin(beta)*cos(a)- sin(k)));
            if(phi1 < 0)phi1 += std::numbers::pi;
            double phi2 = k - phi1;

            if((abs(phi1 + phi3 - std::numbers::pi) < 1e-4 && abs(phi2 + phi4 - std::numbers::pi) < 1e-4)){IsStopAtFF = true; std::cout << a << "  flat foldable state " << std::endl;}
            if((abs(phi1 - k/2.0) < 1e-4 && abs(phi2 -  k/2.0) < 1e-4)){IsStopAtEq = true; std::cout << a << "  equal degin angles " << std::endl;}
            if((abs((phi1 + phi4) - std::numbers::pi) < 1e-4 && abs((phi2 + phi3) - std::numbers::pi) < 1e-4)){IsStopAtCon = true; std::cout << a << "  continuation" << std::endl;}
            if(phi1 < std::numbers::pi/2.0)RulingColor.push_back(true); else RulingColor.push_back(false);
            glm::f64vec3 axis_alpha = glm::normalize(glm::cross(e, e2));
            N = glm::normalize(glm::rotate (a, e)  * glm::f64vec4{axis_alpha,1});//phi1の回転軸

            r = glm::rotate(phi1, N) * glm::f64vec4{e,1};r = glm::normalize(r);//新しいruling方向
            n = glm::rotate(phi1, glm::f64vec3{0,0,1})*glm::f64vec4{glm::normalize(_axis[i]->prev->vertex->p - _axis[i]->vertex->p),1};n = glm::normalize(n);//展開図のruling方向


            std::array<glm::f64vec3, 2> tmp = {vecLen*r + x, x}; NewRuling.push_back(tmp);
            tmp = {vecLen*n + _axis[i]->vertex->p, _axis[i]->vertex->p}; //NewRuling2d.push_back(tmp);

            double sin_a = (sin(phi1)*sin(a)/sin(phi2) > 1) ? 1: (sin(phi1)*sin(a)/sin(phi2) < -1)? -1: (sin(phi1)*sin(a)/sin(phi2));
            r2 = MathTool::ProjectionVector(r,e2); r2 = glm::normalize(r2);
            glm::f64vec3 ax_beta = e;
            ax_beta = MathTool::ProjectionVector(ax_beta, e2); ax_beta = glm::normalize(ax_beta);
            double cos_a = (glm::dot(r2, ax_beta) > 1) ? 1: (glm::dot(r2, ax_beta) < -1)? -1: glm::dot(r2, ax_beta);//値が期待しているのと違う
            cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
            a2 = (sin_a >= 0 && cos_a >= 0)? asin(sin_a): (sin_a >= 0 && cos_a < 0)?std::numbers::pi - asin(sin_a): (sin_a < 0 && cos_a < 0)? std::numbers::pi + asin(sin_a): 2.0*std::numbers::pi + asin(sin_a);
            e_bef = (MathTool::ProjectionVector(e,e2)); e_bef = glm::normalize(e_bef);
        res << beta << ", " << a <<  ", " << phi1 << ", " << phi2 <<", " << phi3 << ", " << phi4 << ", " << k <<
               ", " << (phi1 + phi3- std::numbers::pi) << ", " << (phi2 + phi4 - std::numbers::pi) << ", " << (phi1- k/2.0) <<", "<< (phi2- k/2.0) <<", " << (phi1 + phi4- std::numbers::pi)<<", " << (phi2 + phi3- std::numbers::pi)<<  std::endl;
        }

        angle += 1e-4;
    }
    update();
    emit foldingSignals();
}


void GLWidget_2D::changeBetaValue(int val){
    if(model->Faces.size() < 2 || IsStop4Debug)return;

    bool hasFoldLine = false;

    double a = (double)val * std::numbers::pi/180000.0;
    std::vector<HalfEdge*> _axis;
    for(const auto& h: model->Edges){
        if(h->edgetype == EdgeType::r){
            hasFoldLine = true;
            if(std::find(_axis.begin(), _axis.end(), h) == _axis.end() && std::find(_axis.begin(), _axis.end(), h->pair) == _axis.end()){
                _axis.push_back(h);
            }
        }
    }
    if(!hasFoldLine)return;
    model->deform();

    auto IsSame = [](double c, double b){
        return (abs(c - b) <= 1e-8)? true: false;
    };
    double a2 = 0;
    glm::f64vec3 e, e2, x;
    glm::f64vec3 N, r, n, r2;
    glm::f64vec3 e_bef;
    double vecLen = 50;
    NewRuling.clear();NewRuling2d.clear(); RulingColor.clear();
    for(int i = 0; i < (int)model->Faces.size() - 1; i++){
        x = _axis[i]->vertex->p3;
        bool IsVallyFoldLine = false;
        e = glm::normalize(_axis[i]->prev->vertex->p3 - x);
        e2 = glm::normalize(_axis[i]->pair->next->next->vertex->p3 - x);
        double beta = (1.0 - abs(_axis[i]->color/255.0)) * std::numbers::pi;
        if(_axis[i]->color < 0)IsVallyFoldLine = true;
        if(i != 0){
            glm::f64vec3 e_next = glm::normalize(MathTool::ProjectionVector(e2,-e));
            double tau = (glm::dot(e_bef, e_next) > 1) ? 0: (glm::dot(e_bef, e_next) < -1) ? std::numbers::pi: glm::angle(e_bef,e_next);
            tau = glm::angle(e_bef, e_next);
            if(glm::dot(glm::cross(e_bef, e_next), e) > 0) tau = 2.0*std::numbers::pi - tau;
            a = (a2 - tau < 0) ? a2 - tau + 2.0 * std::numbers::pi: a2 - tau;
            //std::cout << "---------------------------------" << std::endl;
            if(!std::isfinite(a))std::cout << "a is not finite " << i << std::endl;
            if(!std::isfinite(a2))std::cout << "a' is not finite " << i << std::endl;
            std::cout << i << " : tau = " << tau << " , a = " <<  glm::degrees(a) << " , a' = " << glm::degrees(a2) << std::endl;
        }
        double phi3 = glm::angle(e2, glm::normalize(_axis[i]->next->vertex->p3 - x));
        double phi4 = glm::angle(e, glm::normalize(_axis[i]->next->vertex->p3 - x));
        double k = 2.0 * std::numbers::pi - phi3 - phi4;
        if(!(beta <= phi3 + phi4 && phi3 + phi4 <= 2.0 * std::numbers::pi - beta)){std::cout <<"exceed boundary condition   " << beta << " , " << phi3+phi4<< std::endl; exit(0);}
        double phi1 = IsSame(glm::sin(beta)*glm::cos(a), glm::sin(k))? std::numbers::pi/2.0: atan2((glm::cos(k) - glm::cos(beta)),(glm::sin(beta)*glm::cos(a)- glm::sin(k)));
        if(phi1 < 0){phi1 = std::numbers::pi + phi1;}
        double phi2 = k - phi1;
        std::cout <<"angle "<< i << " : " << a << ", " <<  std::setprecision(10) <<  glm::degrees(phi1) << ", " << glm::degrees(phi2) << ", " << glm::degrees(phi3) << ", " << glm::degrees(phi4)<< " : " <<  glm::degrees(a) <<  std::endl;
        if(IsSame(phi1 + phi3, std::numbers::pi) && IsSame(phi2 + phi4, std::numbers::pi)){IsStopAtFF = true; std::cout << i << "  flat foldable state " << std::endl;}
        if(IsSame(phi1, std::numbers::pi/2.0) && IsSame(phi2, std::numbers::pi/2.0)){IsStopAtEq = true; std::cout << i << "  equal degin angles " << std::endl;}
        if(IsSame(phi1 + phi4, std::numbers::pi) && IsSame(phi2 + phi3, std::numbers::pi)){IsStopAtCon = true; std::cout << i << "  continuation" << std::endl;}
        if(phi1 < std::numbers::pi/2.0)RulingColor.push_back(true); else RulingColor.push_back(false);
        glm::f64vec3 axis_alpha = glm::normalize(glm::cross(e, e2));
        if(IsVallyFoldLine){N = glm::normalize(glm::rotate (a, e2)  * glm::f64vec4{-axis_alpha,1});std::cout <<"Vally" << std::endl;}//phi1の回転軸
        else {N = glm::normalize(glm::rotate (a, e)  * glm::f64vec4{axis_alpha,1});std::cout << "Mount" << std::endl;}//phi1の回転軸

        r = glm::rotate(phi1, N) * glm::f64vec4{e,1};r = glm::normalize(r);//新しいruling方向
        n = glm::rotate(phi1, glm::f64vec3{0,0,1})*glm::f64vec4{glm::normalize(_axis[i]->prev->vertex->p - _axis[i]->vertex->p),1};n = glm::normalize(n);//展開図のruling方向

        std::array<glm::f64vec3, 2> tmp = {vecLen*r + x, x}; NewRuling.push_back(tmp);
        //tmp = {vecLen*N + x, x}; NewRuling.push_back(tmp);

        tmp = {vecLen*n + _axis[i]->vertex->p, _axis[i]->vertex->p}; NewRuling2d.push_back(tmp);

        double sin_a = (sin(phi1)*sin(a)/sin(phi2) > 1) ? 1: (sin(phi1)*sin(a)/sin(phi2) < -1)? -1: (sin(phi1)*sin(a)/sin(phi2));
        r2 = MathTool::ProjectionVector(r,e2); r2 = glm::normalize(r2);
        glm::f64vec3 ax_beta = e;
        ax_beta = MathTool::ProjectionVector(ax_beta, e2); ax_beta = glm::normalize(ax_beta);
        double cos_a = (glm::dot(r2, ax_beta) > 1) ? 1: (glm::dot(r2, ax_beta) < -1)? -1: glm::dot(r2, ax_beta);//値が期待しているのと違う
        cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
        a2 = (sin_a >= 0 && cos_a >= 0)? asin(sin_a): (sin_a >= 0 && cos_a < 0)?std::numbers::pi - asin(sin_a): (sin_a < 0 && cos_a < 0)? std::numbers::pi + abs(asin(sin_a)): 2.0*std::numbers::pi + asin(sin_a);
        //std::cout <<"sin " << glm::degrees(sin_a) << " , cos " << glm::degrees(cos_a) << " , a' =  " << glm::degrees(a2) << std::endl;
        e_bef = (MathTool::ProjectionVector(e,e2)); e_bef = glm::normalize(e_bef);
       }
        std::cout <<"||||||||||||||||||||||||||||||" << std::endl;

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

    //折り線の描画
    {
        for(auto&fl: model->FL){

            glBegin(GL_LINE_STRIP);
            if(visibleCurve){glColor3d(0,0,0); for(auto&v: fl->CurvePts)glVertex2d(v.x, v.y);}
            else{glColor3d(0,0,1); for(auto&v: fl->Curve_res2d)glVertex2d(v.x, v.y);}
            glEnd();
            std::vector<glm::f64vec3> Pts;
            Pts = fl->getCtrlPt();
            glColor3d(0,0.3,0.3);
            glPointSize(5);
            for(auto&v: Pts){
                glBegin(GL_POINTS);
                glVertex2d(v.x,v.y);
                glEnd();
            }

            glColor3d(0,1,0);
            glPointSize(5);
            for(auto&h: fl->FoldingCurve){
                glBegin(GL_POINTS);
                glVertex2d(h->vertex->p.x, h->vertex->p.y);
                glEnd();
            }

        }
    }

    if(model->outline->IsClosed()){
        float r, g, b = 0;
        glPolygonOffset(0.0,1.0);
        HalfEdge *_refHE = assignment_refHE();
        for(auto&edge: model->Edges){
            glColor3d(0,0,0);
            glLineWidth(1);

            if(edge->edgetype == EdgeType::r){
                if(drawtype == PaintTool::NewGradationMode){
                    glLineWidth(rulingWidth);
                    if(edge->IsCrossed != -1)glColor3d(0,1,0);
                    else {
                        if(edge->color == 0) r = g = b = 0.4;
                        else if(edge->color > 0){
                            r = 1; g = b =1 - edge->color/255.0;
                        }else{
                            b = 1;
                            g = r = 1 + edge->color/255.0;
                        }
                        glColor3d(r,g,b);
                    }
                }else{
                    if(edge->IsCrossed  != -1)glColor3d(0,1,0);
                    else glColor3d(0.4,0.4,0.4);
                    glLineWidth(1.f);
                }
                if(edge == _refHE || edge->pair == _refHE)glColor3d(1,1,0);
            }else if(edge->edgetype == EdgeType::ol || edge->edgetype == EdgeType::fl){
                glColor3d(0, 0, 0);
                glPointSize(3.0f);

                glBegin(GL_POINTS);
                glVertex2d(edge->vertex->p.x, edge->vertex->p.y);
                glEnd();
            }
            glBegin(GL_LINES);
            glVertex2d(edge->vertex->p.x, edge->vertex->p.y);
            glVertex2d(edge->next->vertex->p.x, edge->next->vertex->p.y);
            glEnd();

        }
    }else{
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

    }

    glLineWidth(1.0);
    //曲線の制御点
    if(drawtype != PaintTool::NewGradationMode){
        for(int i = 0; i < (int)model->crvs.size(); i++){
            if((i == model->getSelectedCurveIndex(mapFromGlobal(QCursor::pos()))) ||
                    ((drawtype == PaintTool::Bezier_r || drawtype == PaintTool::Bspline_r || drawtype == PaintTool::Line_r || drawtype == PaintTool::Arc_r)
                     && SelectedCurveIndex == i)
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
    }

    for(int i = 0; i < (int)NewRuling2d.size(); i++){
        auto r = NewRuling2d[i];
        if(RulingColor[i])glColor3d(1,0,0);
        else glColor3d(0,0,1);
        glBegin(GL_LINES);
        glVertex2d(r[0].x, r[0].y);
        glVertex2d(r[1].x, r[1].y);
        glEnd();
    }

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
    auto oriedge = model->outline->getEdges();
    if(e->key() == Qt::Key_V)eraseVec2d = !eraseVec2d;
    if(e->key() == Qt::Key_A) visibleCurve = !visibleCurve;
    //if(e->key() == Qt::Key_0){model->FL[0]->modify2DRulings(model->Faces, model->Edges, model->vertices, oriedge, curveDimention, 0); emit foldingSignals();}
    //if(e->key() == Qt::Key_1){model->FL[0]->modify2DRulings(model->Faces, model->Edges, model->vertices, oriedge, curveDimention, 1); emit foldingSignals();}
    if(e->key() == Qt::Key_2){model->FL[0]->modify2DRulings(model->Faces, model->Edges, model->vertices, oriedge, curveDimention, 2); emit foldingSignals();}
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
            bool res = model->FL[0]->SplitFace4DebugAAAMethod(p_ongrid, model->Faces, model->Edges);
        }
        else if(drawtype == PaintTool::DeleteCurve){
            model->SelectCurve(p);
            std::vector<int> n(1);
            n[0] = model->DeleteCurve();
            emit deleteCrvSignal(n);
            SelectedCurveIndex = -1;
        }
        if(drawtype == PaintTool::Bezier_r || drawtype == PaintTool::Bspline_r || drawtype == PaintTool::Line_r || drawtype == PaintTool::Arc_r){
            model->AddControlPoint(p_ongrid, curveDimention, DivSize);
            model->addRulings();
        }
        else if(drawtype == PaintTool::InsertCtrlPt){
            if(SelectedCurveIndex == -1){
                model->SelectCurve(p);
                SelectedCurveIndex = model->IsSelectedCurve();
            }else{
                if(model->crvs[SelectedCurveIndex]->getCurveType() == CurveType::bsp3 && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){
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
        else if(drawtype == PaintTool::DeleteCtrlPt){
            model->SelectCurve(p);
            model->DeleteControlPoint(p, curveDimention, DivSize);
        }else if(drawtype == PaintTool::MoveCtrlPt){
            SelectedCurveIndex = model->searchPointIndex(p, movePt, 0);

        }

    }else if(e->button() == Qt::RightButton){
        if(drawtype == PaintTool::Rectangle_ol || drawtype == PaintTool::Polyline_ol || drawtype == PaintTool::Polygon_ol) model->outline->eraseVertex();
        if(drawtype == PaintTool::Bezier_r || drawtype == PaintTool::Bspline_r || drawtype == PaintTool::Line_r || drawtype == PaintTool::Arc_r){
            if(SelectedCurveIndex != -1) model->crvs[SelectedCurveIndex]->eraseCtrlPt(curveDimention, crvPtNum);
        }
        else if(drawtype == PaintTool::FoldLine_line || drawtype == PaintTool::FoldLine_arc || drawtype == PaintTool::FoldLine_bezier){
            bool hasRulings = model->AddControlPoint_FL(p_ongrid, 1, curveDimention);
        }

    }
    if(model->outline->IsClosed()){
        //model->deform();
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

    if(drawtype == PaintTool::MoveCtrlPt)model->MoveCurvePoint(p_ongrid,SelectedCurveIndex, movePt, curveDimention, DivSize);

    if(SelectedCurveIndex != -1){
        if(drawtype == PaintTool::InsertCtrlPt &&  model->crvs[SelectedCurveIndex]->getCurveType() == CurveType::bsp3 && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){//制御点の挿入(B-spline)
            model->crvs[SelectedCurveIndex]->InsertControlPoint2(p_ongrid);
        }
        
        if(drawtype != PaintTool::NewGradationMode && model->outline->IsClosed() && model->crvs[SelectedCurveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){
            //if(model->crvs[SelectedCurveIndex]->getCurveType() == 0)model->crvs[SelectedCurveIndex]->BezierRulings(model->outline->vertices,DivSize,crvPtNum);
            if(model->crvs[SelectedCurveIndex]->getCurveType() == CurveType::bsp3) model->crvs[SelectedCurveIndex]->BsplineRulings(model->outline,DivSize,crvPtNum, curveDimention);
            if(model->crvs[SelectedCurveIndex]->getCurveType() == CurveType::line) model->crvs[SelectedCurveIndex]->LineRulings(model->outline,DivSize);
            if(model->crvs[SelectedCurveIndex]->getCurveType() == CurveType::arc) model->crvs[SelectedCurveIndex]->ArcRulings(model->outline,DivSize);
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
        if(SelectedCurveIndex == -1) std::cout << "no curve is selected"<<std::endl;
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
    delete model;
    tmp_c.clear();
    model = new Model(crvPtNum);
    emit SendNewActiveCheckBox(PaintTool::Reset);
    SelectedCurveIndex = -1;
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
        bool hasRulings = model->FL[0]->ChangeColor(model->outline, DiffWheel, curveDimention);
        emit ColorChangeFrom(0, model->FL[0]->getColor());
        if(hasRulings){
            bool res;
            model->FL[0]->applyCurvedFolding(model->Faces, model->Edges, model->vertices, curveDimention);
        }
    }else if(drawtype == PaintTool::NewGradationMode){
        if(refHE == nullptr || std::find(model->Edges.begin(), model->Edges.end(), refHE) == model->Edges.end()){
            std::cout<<"there is no refHE"<<std::endl;
            return;
        }
        model->setGradationValue(DiffWheel, refHE, InterpolationType, CurvePath);
        emit ColorChangeFrom(0, refHE->color);
        model->deform();
        if(model->outline->IsClosed() && !model->FL.empty()){
            auto oriedge = model->outline->getEdges();
            //model->FL[0]->modify2DRulings(model->Faces, model->Edges, model->vertices, oriedge, curveDimention);
        }
        if(isVisibleTo(gw)) emit CurvePathSet(CurvePath);
    }
    emit foldingSignals();
    update();
}

void GLWidget_2D::addPoints_intplation(QMouseEvent *e, QPointF& p){
    glm::f64vec3 curPos{p.x(), p.y(), 0};
    double dist = 10;
    if(drawtype == PaintTool::NewGradationMode){
        refHE = nullptr;
        for(auto* _he: model->Edges){
            if(_he->edgetype == EdgeType::ol || _he->edgetype == EdgeType::cl)continue;
            double d = glm::length(glm::cross((curPos - _he->vertex->p), _he->vertex->p - _he->next->vertex->p))/glm::length(_he->vertex->p - _he->next->vertex->p);
            if(d < dist){
                dist = d; refHE = _he;
            }
        }
        if(refHE == nullptr || std::find(model->Edges.begin(), model->Edges.end(), refHE) == model->Edges.end()){std::cout<<"not found" << std::endl; return;}
        model->setGradationValue(0, refHE, InterpolationType, CurvePath);
        emit ColorChangeFrom(0, refHE->color);
        model->deform();
        if(model->outline->IsClosed() && !model->FL.empty()){
            auto oriedge = model->outline->getEdges();
            model->FL[0]->modify2DRulings(model->Faces, model->Edges, model->vertices, oriedge, curveDimention);
        }
    }else if(drawtype == PaintTool::FoldlineColor){}
    emit foldingSignals();
    update();
}

void GLWidget_2D::ApplyNewGradationMode(){
    if(refHE == nullptr){std::cout << "you neeed to add gradation point"<<std::endl; return;}
    emit foldingSignals();
    if(isVisibleTo(gw)) emit CurvePathSet(CurvePath);
    update();
}

void GLWidget_2D::getGradationFromSlider(int val){
    if(refHE == nullptr) return;
    refHE->color = val;
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
    for(auto& _he: model->Edges){
        if(_he->edgetype != EdgeType::r)continue;
        //if(_he->next == nullptr)continue;
        double d = glm::length(glm::cross((curPos - _he->vertex->p), _he->vertex->p - _he->next->vertex->p))/glm::length(_he->vertex->p - _he->next->vertex->p);
        if(d < dist){
            dist = d; he = _he;
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
