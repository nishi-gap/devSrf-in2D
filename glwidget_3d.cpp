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
    switchTNB = 0;
    glm::f64vec3 up{0,1,0};
    arccam = ArcBallCam(glm::f64vec3{0,0,-100}, center, up);
    drawdist = 0.0;
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
    arccam = ArcBallCam(glm::f64vec3{0,0,-100}, center, up);
    glViewport(s.width(),0,s.width(),s.height());
    glEnable(GL_POLYGON_OFFSET_FILL);
    glEnable(GL_DEPTH_TEST);
}

void EdgeRecconection(const std::vector<Vertex*>& Poly_V, std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges){
    int n = Poly_V.size();
    std::vector<std::vector<Vertex>> surface;
    for(auto&f: Faces){
        std::vector<Vertex> mesh;
        HalfEdge *he = f->halfedge;
        do{
            Vertex v(he->vertex->p, he->vertex->p3);
            mesh.push_back(v);
            he = he->next;
        }while(he != f->halfedge);
        surface.push_back(mesh);
    }
    for(int i = 0; i < n; i++){
        Vertex *v1 = Poly_V[i], *v2 = Poly_V[(i + 1) % n];
        HalfEdge *h1 = nullptr, *h2 = nullptr;
        std::vector<PointOnLine> PoL;
        for(auto&e: Edges){
            if((e->vertex == v1) && e->edgetype == EdgeType::ol)h1 = e;
            if((e->vertex == v2) && e->edgetype == EdgeType::ol)h2 = e;

            if((e->vertex == v1 || e->vertex == v2) || e->edgetype != EdgeType::ol)continue;
            if(abs(glm::distance(e->vertex->p, v1->p) + glm::distance(e->vertex->p, v2->p) - glm::distance(v1->p, v2->p)) < 1e-5){
                double t = glm::distance(e->vertex->p, v2->p)/glm::distance(v1->p, v2->p);
                auto p = PointOnLine(t, e->vertex);
                PoL.push_back(p);
            }
        }
        std::sort(PoL.begin(), PoL.end());

        if(PoL.empty()){//polygon同士の頂点で一度つなぎなおす
            h1->prev = h2; h2->next = h1;
        }else{
            for(int j = 1; j < (int)PoL.size(); j++){
                HalfEdge *edge1 = nullptr, *edge2= nullptr;
                for(auto&h: PoL[j-1].v->halfedge){
                    if(h->edgetype == EdgeType::ol)edge1 = h;
                }
                for(auto&h: PoL[j].v->halfedge){
                    if(h->edgetype != EdgeType::ol)edge2 = h;
                }
                edge1->next = edge2; edge2->prev = edge1;
            }
            for(auto&e: PoL.front().v->halfedge){
                if(e->edgetype == EdgeType::ol)continue;
                e->prev = h2;
                h2->next = e;
            }
            for(auto&e: PoL.back().v->halfedge){
                if(e->edgetype != EdgeType::ol)continue;
                e->next = h1;
                h1->prev = e;
            }
        }
    }
    std::vector<bool> IsReconnected(Edges.size(), false);
    int faceIndex = 0;

    for(int i = 0; i < (int)IsReconnected.size(); i++){
        if(IsReconnected[i])continue;
        HalfEdge *h = Edges[i];
        do{
            int j = std::distance(Edges.begin(), std::find(Edges.begin(), Edges.end(), h));
            IsReconnected[j] = true;
            Edges[j]->face = Faces[faceIndex];
            h = h->next;
        }while(h != Edges[i]);
        faceIndex++;
    }
}

void GLWidget_3D::setVertices(const Faces3d Faces, const Polygon_V Poly_V, const HalfEdges Edges, const Ruling3d& _SingleRuling, const Ruling3d& _AllRulings){
    std::vector<glm::f64vec3> vertices, centers;
    drawdist = 0.0;
    Vertices.clear();
    TriMeshs.clear();
    glm::f64vec3 _center;
    glm::f64mat4x4 Mirror = glm::mat4(1.0f); Mirror[1][1] = -1;
    glm::f64mat4x4 Scale = glm::scale(glm::f64vec3{0.1, 0.1, 0.1});
    std::vector<std::array<glm::f64vec3, 3>> trimesh;
    FoldLineVertices.clear();

    SingleRuling.clear();
    for(auto&v: _SingleRuling){
        std::array<glm::f64vec3, 2> tmpV{glm::f64vec3(Scale * Mirror * glm::f64vec4(v[0],1)), glm::f64vec3(Scale * Mirror * glm::f64vec4(v[1],1))};
        SingleRuling.push_back(tmpV);
    }

    AllRulings.clear();
    for(auto&v: _AllRulings){
        std::array<glm::f64vec3, 2> tmpV{glm::f64vec3(Scale * Mirror * glm::f64vec4(v[0],1)), glm::f64vec3(Scale * Mirror * glm::f64vec4(v[1],1))};
        AllRulings.push_back(tmpV);
    }

    TriMeshs.clear();
    Faces3d _faces;
    Face *f;
    for(auto _f: Faces){
        f = _f;
        _faces.push_back(f);
    }
    HalfEdges _edges;
    HalfEdge *h;
    for(auto _h: Edges){
        h = _h;
        _edges.push_back(h);
    }
    //EdgeRecconection(Poly_V, _faces, _edges);
    for(auto&f: _faces){
        vertices.clear();
        HalfEdge *he = f->halfedge;       
        do{
            glm::f64vec3 v = glm::f64vec3(Scale * Mirror * glm::f64vec4(he->vertex->p3,1));
            vertices.push_back(v);
            he = he->next;           
        }while(he != f->halfedge);
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

    center = glm::f64vec3{0,0,0};
    update();
}

void GLWidget_3D::receive(std::vector<std::vector<glm::f64vec3>>& l, std::vector<std::vector<glm::f64vec3>>& r, glm::f64vec3 center){
    std::vector<std::array<glm::f64vec3, 3>> TriMesh;
    TriMeshs.clear();
    Vertices = l;
    for(auto& V : l){
        Triangulation(V, TriMesh);
        for(auto&tri: TriMesh){
            TriMeshs.push_back(tri);
            left.push_back(1);
        }
    }
    for(auto&V: r){
        Triangulation(V, TriMesh);
        for(auto&tri: TriMesh){
            TriMeshs.push_back(tri);
            left.push_back(0);
        }
        Vertices.push_back(V);
    }

    TransX = 0.1 * center.x;
    TransY = 0.05 * center.y;
    this->center = 0.1 * center;


    update();
}

inline void GLWidget_3D::dispV(glm::f64vec3 p){
    glVertex3d(p.x, p.y, p.z);
}
void GLWidget_3D::paintGL(){
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    QSize s = this->size();

    glLoadIdentity();
    perspective(30.0f, (float)s.width() / (float)s.height(), 1.f, 100.f);
    glm::f64mat4 RotY = glm::rotate(angleY, glm::f64vec3{0,1,0}), RotX = glm::rotate(angleX, glm::f64vec3(1,0,0));
    glm::f64mat4 Trans = glm::translate(glm::f64vec3{- TransX, - TransY, - TransZ});
    glm::f64vec3 camPos = RotY * RotX * Trans * glm::f64vec4{center, 1};


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
        glColor3d(0,0, (double)i/AllRulings.size());
        glBegin(GL_LINES);
        glVertex3d(AllRulings[i][0].x, AllRulings[i][0].y, AllRulings[i][0].z);
        glVertex3d(AllRulings[i][1].x, AllRulings[i][1].y, AllRulings[i][1].z);
        glEnd();
    }

    glPolygonOffset(1.f,2.f);

    for(int i = 0; i < (int)SingleRuling.size(); i++){
        if(SingleRuling.size()% 3 == 0){
            if(i% 3 == 0)glColor3d(1,0, 0);//e
            else if(i% 3 == 1)glColor3d(0,1, 0);//N
            else if(i% 3 == 2)glColor3d(0,0, 1);//r
        }else glColor3d(0,0, (double)i/SingleRuling.size());
        glLineWidth(3);
        glBegin(GL_LINES);

        glVertex3d(SingleRuling[i][0].x, SingleRuling[i][0].y, SingleRuling[i][0].z);
        glVertex3d(SingleRuling[i][1].x, SingleRuling[i][1].y, SingleRuling[i][1].z);
        glEnd();
    }
    glPolygonOffset(0.f,0.f);

}

void GLWidget_3D::DrawMesh(bool isFront){
    glPolygonOffset(1.f,0.5f);
    glEnable(GL_CULL_FACE);
    if(isFront){
        glCullFace(GL_FRONT);
        glColor3d(0.6,0.6,0.6);
    }
    else {
        glCullFace(GL_BACK);
        glColor3d(0.9,0.9,0.9);
    }
    for(int i = 0; i < (int)TriMeshs.size(); i++){
        //if(left[i] == 1)glColor3d(0.8, 0.2, 0.2);
        //else glColor3d(0.2, 0.2, 0.8);
        glBegin(GL_POLYGON);
        for (auto& v : TriMeshs[i]) {
            glVertex3d(v.x, v.y, v.z);
        }
        glEnd();
    }
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glDisable(GL_CULL_FACE);
    glPolygonOffset(0.f,0.f);
}

void GLWidget_3D::DrawMeshLines(){
    //std::vector<std::vector<glm::f64vec3>> TriMeshs;
    glColor3d(0.f, 0.f, 0.f);
    glLineWidth(1.0f);

    for(auto& f: Vertices){
        glBegin(GL_LINE_LOOP);
        for(auto& v: f){ glVertex3d(v.x, v.y, v.z);}
        glEnd();
    }
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
    arccam.zoom(z);
    TransZ += z;
    update();
}

void GLWidget_3D::keyPressEvent(QKeyEvent *e){
    if(e->key() == Qt::Key_M)eraseMesh = !eraseMesh;
    if(e->key() == Qt::Key_C) eraseCtrlPt = !eraseCtrlPt;
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
    update();
}

void GLWidget_3D::mouseMoveEvent(QMouseEvent *e){

    QPointF curPos = this->mapFromGlobal(QCursor::pos());
    QPointF diff = curPos - befPos;

    if(actionType == 1){//平行移動
        TransX -=  0.3 * diff.x();
        TransY +=  0.3 * diff.y();
        arccam.pan(glm::f64vec2{diff.x(), diff.y()});
    }else if(actionType == 2){//回転
        angleX += diff.x();
        angleY += diff.y();
        arccam.rotate(glm::f64vec2{befPos.x(), befPos.y()}, glm::f64vec2{curPos.x(), curPos.y()});
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
