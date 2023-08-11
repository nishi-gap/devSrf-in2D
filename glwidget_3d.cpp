#include "glwidget_3d.h"
using namespace MathTool;

GLWidget_3D::GLWidget_3D(QWidget *parent):QOpenGLWidget(parent)
{

    this->Vertices.clear();
    TriMeshs.clear();
    this->firstRotate = true;
    actionType = 0;
    center = glm::f64vec3{0,0,0};
    eraseMesh = eraseCtrlPt = eraseCrossPt = eraseVec = eraseCurve = false;
    VisiblePlanarity = false;
    switchTNB = 0;
    glm::f64vec3 up{0,1,0};
    //arccam = ArcBallCam(glm::f64vec3{0,0,-100}, center, up);
    drawdist = 0.0;
    drawEdgePlane = -1;
    IsEraseNonFoldEdge = false;

    Mirror = glm::mat4(1.0f); Mirror[1][1] = -1;
    Scale = glm::scale(glm::f64vec3{0.1, 0.1, 0.1});
}
GLWidget_3D::~GLWidget_3D(){

}

void GLWidget_3D::initializeGL(){
    initializeOpenGLFunctions();
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    QSize s = this->size();
    TransX = 50.f, TransY = 0.f, TransZ = -200.f;
    angleY = 90;
    glm::f64vec3 up{0,1,0};
    //arccam = ArcBallCam(glm::f64vec3{0,0,-100}, center, up);
    glViewport(s.width(),0,s.width(),s.height());
    glEnable(GL_POLYGON_OFFSET_FILL);
    glEnable(GL_DEPTH_TEST);
}

void GLWidget_3D::EraseNonFoldEdge(bool state){
    IsEraseNonFoldEdge = state;
    update();
}

void GLWidget_3D::setVertices(const Lines Surface,  const Lines Rulings,  const FoldLine3d FldCrvs, const Ruling3d& _AllRulings){

    auto Planerity  = [](const std::vector<Vertex*>& vertices, const Lines Poly_V)->double{
        if(vertices.size() == 3)return 0.0;
        else{
            std::vector<glm::f64vec3> QuadPlane;
            for(auto&v: vertices){
                bool IsOutlineVertices = false;
                for(auto&p: Poly_V){
                    if(p->v->p3 == v->p3)IsOutlineVertices = true;
                }
                if(!IsOutlineVertices)QuadPlane.push_back(v->p3);
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


    drawdist = 0.0;
    Vertices.clear();
    TriMeshs.clear();
    PlanarityColor.clear();
    C.clear();
    glm::f64vec3 _center;

    std::vector<std::array<glm::f64vec3, 3>> trimesh;
    FoldLineVertices.clear();

    AllRulings.clear();
    for(auto&v: _AllRulings){
        std::array<glm::f64vec3, 2> tmpV{glm::f64vec3(Scale * Mirror * glm::f64vec4(v[0],1)), glm::f64vec3(Scale * Mirror * glm::f64vec4(v[1],1))};
        AllRulings.push_back(tmpV);
    }

    TriMeshs.clear();
    drawEdges.clear();
    std::vector<std::vector<Vertex*>> Polygons;
    std::vector<Vertex*> polygon;
    for(auto& l: Surface) polygon.push_back(l->v);
    Polygons.push_back(polygon);
    if(!FldCrvs.empty() &&!FldCrvs[0]->FoldingCurve.empty()){
        for(auto&_FC: FldCrvs){
            std::vector<Vertex4d> FldCrv = _FC->FoldingCurve;
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
                    if(!MathTool::is_point_on_line(v->p, Poly[i]->p, Poly[(i + 1) % (int)Poly.size()]->p))continue;
                    Poly.insert(Poly.begin() + i, v);
                    break;
                }
                int f_ind, s_ind;
                for(int i = 0; i < (int)Poly.size(); i++){
                    if(Poly[i] == o)f_ind = i;
                    if(Poly[i] == v)s_ind = i;
                }
                int i_min = std::min(f_ind, s_ind), i_max = std::max(f_ind, s_ind);
                std::vector<Vertex*> poly2(Poly.begin() + i_min, Poly.begin() + i_max + 1);
                Poly.erase(Poly.begin() + i_min, Poly.begin() + i_max -1);
                Polygons.push_back(poly2);
            }
        };

        for(auto&FC: FldCrvs){
            for(auto itr = FC->FoldingCurve.begin() + 1; itr != FC->FoldingCurve.end() - 1; itr++){
                SplitPolygon(Polygons, itr->first, itr->second);
                SplitPolygon(Polygons, itr->first, itr->third);
            }
       }
    }else{
        for(auto itr_r = Rulings.begin(); itr_r != Rulings.end(); itr_r++){
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
                drawEdges.push_back({(*itr_r)->v->p3, (*itr_r)->o->p3});
            }
        }
    }

    for(auto& polygon: Polygons){
        //std::vector<Vertex*> p_sort = SortPolygon(polygon);
        if(!Rulings.empty())PlanarityColor.push_back(Planerity(polygon, Surface));
        else PlanarityColor.push_back(0);
        std::vector<glm::f64vec3> vertices;
        for(auto& p: polygon)vertices.push_back(glm::f64vec3(Scale * Mirror * glm::f64vec4(p->p3,1)));
        Vertices.push_back(vertices);
    }

    for(auto&V: Vertices){
        Triangulation(V, trimesh);
        TriMeshs.insert(TriMeshs.end(), trimesh.begin(), trimesh.end());
    }
    for(auto&tri : TriMeshs){
        _center = (tri[2] + tri[1] + tri[0])/3.0;
        double area = glm::length(glm::cross(tri[2] - tri[0], tri[1] - tri[0]))/2.0;
        center += _center * area;
    }
    for(auto&e: drawEdges){
        e[0] = glm::f64vec3(Scale * Mirror * glm::f64vec4(e[0],1)); e[1] = glm::f64vec3(Scale * Mirror * glm::f64vec4(e[1],1));
    }

    center = glm::f64vec3{0,0,0};

    update();
}

inline void GLWidget_3D::dispV(glm::f64vec3 p){
    glVertex3d(p.x, p.y, p.z);
}
void GLWidget_3D::ReceiveParam(std::vector<std::vector<glm::f64vec3>>&_C){
    C = _C;
    for(auto&c: C){
        for(auto&p: c) p = glm::f64vec3(Scale * Mirror * glm::f64vec4(p,1));
    }
}

void GLWidget_3D::paintGL(){
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    QSize s = this->size();

    glLoadIdentity();
    perspective(30.0f, (float)s.width() / (float)s.height(), 1.f, 100.f);
    //glm::f64mat4 RotY = glm::rotate(angleY, glm::f64vec3{0,1,0}), RotX = glm::rotate(angleX, glm::f64vec3(1,0,0));
    //glm::f64mat4 Trans = glm::translate(glm::f64vec3{- TransX, - TransY, - TransZ});
    //glm::f64vec3 camPos = RotY * RotX * Trans * glm::f64vec4{center, 1};


    glScaled(0.1, 0.1, 0.1);
    glTranslated(-center.x - TransX, -center.y - TransY, -center.z -drawdist + TransZ);
    glRotated(0.2 * angleX, 0.0, 1.0, 0.0);
    glRotated(0.2 * angleY, 1.0, 0.0, 0.0);

    if(!eraseMesh){
        DrawMeshLines();
        DrawMesh(true);
        DrawMesh(false);
    }
    //

    glPolygonOffset(0.5f,1.f);
    for(int i = 0; i < (int)AllRulings.size(); i++){
        glColor3d(0, (double)i/AllRulings.size(), (double)i/AllRulings.size());
        //if(i % 2 == 0)glColor3d(1,0,0);
        //else glColor3d(0,1,0);
        glBegin(GL_LINES);
        glVertex3d(AllRulings[i][0].x, AllRulings[i][0].y, AllRulings[i][0].z);
        glVertex3d(AllRulings[i][1].x, AllRulings[i][1].y, AllRulings[i][1].z);
        glEnd();

        glPointSize(5);
        glBegin(GL_POINTS);
        glVertex3d(AllRulings[i][0].x, AllRulings[i][0].y, AllRulings[i][0].z);
        glEnd();
    }

    glPolygonOffset(1.f,2.f);



    for(const auto&c: C){
        glPointSize(6);
        glColor3d(1,0,0);
        for(const auto&v: c){
            glBegin(GL_POINTS);
            glVertex3d(v.x, v.y, v.z);
            glEnd();
        }
        glEnable(GL_LINE_STIPPLE);
        glLineStipple(1 , 0xF0F0);
        glBegin(GL_LINE_STRIP);
        glColor3d(0.4, 0.4, 0.4);
        for (const auto& v: c)glVertex3d( v.x, v.y, v.z);
        glEnd();
        glDisable(GL_LINE_STIPPLE);

    }

    glPolygonOffset(0.f,0.f);

}

void GLWidget_3D::DrawMesh(bool isFront){
    glPolygonOffset(1.f,0.5f);

    if(VisiblePlanarity){
        for(int i = 0; i < (int)Vertices.size(); i++){
            if(drawEdgePlane == i)continue;
            if(PlanarityColor[i] >= th_planarity){
                glColor3d(1,0,0);
            }else{
                glColor3d(PlanarityColor[i]/th_planarity, 1 - PlanarityColor[i]/th_planarity, 0);
            }

            glBegin(GL_POLYGON);
            for (auto& v : Vertices[i]) {
                glVertex3d(v.x, v.y, v.z);
            }
            glEnd();
        }
    }else{
        glEnable(GL_CULL_FACE);
        if(isFront){
            glCullFace(GL_FRONT);
            glColor3d(0.6,0.6,0.6);
        }
        else {
            glCullFace(GL_BACK);
            glColor3d(0.9,0.9,0.9);

        }
        for(int i = 0; i < (int)Vertices.size(); i++){
            if(drawEdgePlane == i)continue;
            glBegin(GL_POLYGON);
            for (auto& v : Vertices[i]) {
                glVertex3d(v.x, v.y, v.z);
            }
            glEnd();
        }
        glDisable(GL_LIGHTING);
        glDisable(GL_LIGHT0);
        glDisable(GL_CULL_FACE);
    }

    //for(int i = 0; i < (int)TriMeshs.size(); i++){
      //  glBegin(GL_POLYGON);
        //for (auto& v : TriMeshs[i]) {
          //  glVertex3d(v.x, v.y, v.z);
        //}
        //glEnd();
    //}


    glPolygonOffset(0.f,0.f);
}

void GLWidget_3D::DrawMeshLines(){
    //std::vector<std::vector<glm::f64vec3>> TriMeshs;
    glColor3d(0.f, 0.f, 0.f);
    glLineWidth(1.0f);
    for(const auto& e: drawEdges){
        glBegin(GL_LINES);
        glVertex3d(e[0].x, e[0].y, e[0].z); glVertex3d(e[1].x, e[1].y, e[1].z);
        glEnd();
    }
}

void GLWidget_3D::PlanarityDispay(bool state){
    VisiblePlanarity = !VisiblePlanarity;
    //for(auto&c: PlanarityColor)
    //if(c > th_planarity)std::cout << "planarity : " << c << std::endl;
    update();
}

void GLWidget_3D::DrawGrid(){
    double ground_max_x = 300.0;
    double ground_max_y = 300.0;
    glColor3d(0.8, 0.8, 0.8);
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    for (double ly = -ground_max_y; ly <= ground_max_y; ly += 10.0) {
        glVertex3d(-ground_max_x, 0, ly);
        glVertex3d(ground_max_x, 0, ly);
    }
    for (double lx = -ground_max_x; lx <= ground_max_x; lx += 10.0) {
        glVertex3d(lx, 0, ground_max_y);
        glVertex3d(lx, 0, -ground_max_y);
    }
    glEnd();
}

void GLWidget_3D::perspective(GLdouble fovy, GLdouble aspect, GLdouble zNear, GLdouble zFar){
    GLdouble xmin, xmax, ymin, ymax;

    ymax = zNear * tan( fovy * M_PI / 360.0 );
    ymin = -ymax;
    xmin = ymin * aspect;
    xmax = ymax * aspect;

    glFrustum( xmin, xmax, ymin, ymax, zNear, zFar );
}


void GLWidget_3D::mousePressEvent(QMouseEvent *e){
    actionType = (e->button() == Qt::LeftButton)? 1 : (e->button() == Qt::RightButton)? 2: 0;

    befPos = this->mapFromGlobal(QCursor::pos());
    update();
}

void GLWidget_3D::wheelEvent(QWheelEvent *we){
    double z = (we->angleDelta().y() > 0) ? 5 : -5;
    //arccam.zoom(z);
    TransZ += z;
    update();
}

void GLWidget_3D::receiveKeyEvent(QKeyEvent *e){
    if(e->key() == Qt::Key_M)eraseMesh = !eraseMesh;
    if(e->key() == Qt::Key_C)eraseCtrlPt = !eraseCtrlPt;
    if(e->key() == Qt::Key_X)eraseCrossPt = !eraseCrossPt;
    if(e->key() == Qt::Key_D)eraseCurve = !eraseCurve;
    if(e->key() == Qt::Key_S){
        switchTNB = (switchTNB + 1) % 3;
       if(switchTNB == 1) std::cout<<"Now analysis "<<std::endl;
       if(switchTNB == 2) std::cout<<"Now diff "<<std::endl;
    }
    if(e->key() == Qt::Key_O)drawEdgePlane = -1;


    update();
}

void GLWidget_3D::mouseMoveEvent(QMouseEvent *e){

    QPointF curPos = this->mapFromGlobal(QCursor::pos());
    QPointF diff = curPos - befPos;

    if(actionType == 1){//平行移動
        TransX -=  0.3 * diff.x();
        TransY +=  0.3 * diff.y();
        //arccam.pan(glm::f64vec2{diff.x(), diff.y()});
    }else if(actionType == 2){//回転
        angleX += diff.x();
        angleY += diff.y();
        //arccam.rotate(glm::f64vec2{befPos.x(), befPos.y()}, glm::f64vec2{curPos.x(), curPos.y()});
    }
    befPos = curPos;
    update();
}

glm::f64vec3 GLWidget_3D::getVec(float x, float y){
    QSize s = this->size();
    double shortSide = std::min(s.width(), s.height());
    glm::f64vec3 pt(2. * x / shortSide - 1.0, -2.0 * y / shortSide + 1.0, 0.0);
    // z座標の計算
    const double xySquared = pt.x * pt.x + pt.y * pt.y;
    if (xySquared <= 1.0) pt.z = std::sqrt(1.0 - xySquared);// 単位円の内側ならz座標を計算
    else pt = glm::normalize(pt); // 外側なら球の外枠上にあると考える
    return pt;
}

void GLWidget_3D::updateRotate() {
    QPointF curPos = this->mapFromGlobal(QCursor::pos());
    // マウスクリック位置をアークボール球上の座標に変換
    const glm::f64vec3 u = getVec(befPos.x(), befPos.y());
    const glm::f64vec3 v = getVec(curPos.x(), curPos.y());

    // カメラ座標における回転量 (=世界座標における回転量)
    const double angle = std::acos(std::max(-1.0, std::min(glm::dot(u, v), 1.0)));

    // カメラ空間における回転軸
    const glm::f64vec3 rotAxis = glm::cross(u, v);

    // カメラ座標の情報を世界座標に変換する行列
    //const glm::f64mat4 c2wMat = glm::inverse(viewMat);

    // 世界座標における回転軸
    //const glm::f64vec3 rotAxisWorldSpace = glm::vec3(c2wMat * glm::vec4(rotAxis, 0.0f));

    // 回転行列の更新
    //acRotMat = glm::rotate(((4.0 * angle), rotAxisWorldSpace) * acRotMat;
}

