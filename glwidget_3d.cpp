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
    drawdist = 0.0;
}
GLWidget_3D::~GLWidget_3D(){

}

void GLWidget_3D::initializeGL(){
    initializeOpenGLFunctions();
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    QSize s = this->size();
    TransX = 0.f, TransY = 0.f, TransZ = 0.f;
    angleY = 90;
    glViewport(s.width(),0,s.width(),s.height());
    glEnable(GL_POLYGON_OFFSET_FILL);
    glEnable(GL_DEPTH_TEST);
}

void GLWidget_3D::setVertices(const Faces3d& Faces, const Curve3d& _CtrlPts, const Curve3d& _Curve, const CrvFL3d& _CrossPts, const Ruling3d& _Vl,const Ruling3d& _Vr){
    std::vector<glm::f64vec3> vertices, centers;
    drawdist = 0.0;
    Vertices.clear();
    TriMeshs.clear();
    glm::f64vec3 _center;
    glm::f64mat4x4 Mirror = glm::mat4(1.0f); Mirror[1][1] = -1;
    glm::f64mat4x4 Scale = glm::scale(glm::f64vec3{0.1, 0.1, 0.1});
    std::vector<std::array<glm::f64vec3, 3>> trimesh;
    CtrlPts.clear(); Curve.clear(); CrossPts.clear();
    //for(auto&v: _CtrlPts)CtrlPts.push_back(glm::f64vec3(Scale * Mirror * glm::f64vec4(v,1)));
    //for(auto&v: _Curve)Curve.push_back(glm::f64vec3(Scale * Mirror * glm::f64vec4(v,1)));
    copy(_CrossPts.begin(), _CrossPts.end(), back_inserter(CrossPts) );
    for(auto&v: CrossPts){
        v.p3 = Scale * Mirror * glm::f64vec4(v.p3,1);
        v.T3d = Mirror * glm::f64vec4{v.T3d, 1};  v.N3d = Mirror * glm::f64vec4{v.N3d, 1}; v.B3d = Mirror * glm::f64vec4{v.B3d, 1};
        v.Td = Mirror * glm::f64vec4{v.Td, 1};  v.Nd = Mirror * glm::f64vec4{v.Nd, 1}; v.Bd = Mirror * glm::f64vec4{v.Bd, 1};
    }
    TriMeshs.clear();
    double Area = 0.0; center = glm::f64vec3{0,0,0};
    for(auto&f: Faces){
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
        Area += area;
    }
    center /= Area;

    for(auto&f: Faces){
        HalfEdge *he = f->halfedge;
        do{
            if(glm::distance(he->vertex->p3, center) > drawdist)drawdist = glm::distance(he->vertex->p3, center);
        }while(he != f->halfedge);
    }
    //std::cout<<glm::to_string(center)<< drawdist <<std::endl;
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
        DrawMesh(true);
        DrawMesh(false);
        glPolygonOffset(0.2f,0.2f);
        DrawMeshLines();
    }
    //glPolygonOffset(0.f,0.5f);

    glPointSize(5.f);
    glColor3d(1,0,0);
    glBegin(GL_POINT);
    glVertex3d(center.x, center.y, center.z);
    glEnd();

    if(!eraseCurve){
        glColor3d(0,1,0);
        glLineWidth(2);
        glBegin(GL_LINE_STRIP);
        for(auto&v: Curve)dispV(v);
        glEnd();
    }

    if(!eraseCtrlPt){
        glColor3d(1,0,0);
        glPointSize(5);
        for(auto&v: CtrlPts){
            glBegin(GL_POINTS);
            dispV(v);
            glEnd();
        }

        glColor3d(0.4,0.4,0.4);
        glLineWidth(1);
        glEnable(GL_LINE_STIPPLE);
        glLineStipple(1 , 0xF0F0);
        glBegin(GL_LINE_STRIP);
        for(auto&v: CtrlPts) dispV(v);
        glEnd();
        glDisable(GL_LINE_STIPPLE);
    }

    if(!eraseCrossPt){
        glColor3d(0,0,1);
        glPointSize(5);
        for(auto&v: CrossPts){
            glBegin(GL_POINTS);
            dispV(v.p3);
            glEnd();
        }
    }

    for(auto&c: CrossPts){
        glPolygonOffset(1.f, 01.f);
        if(switchTNB == 1){
            glColor3d(1,0,0);//T
            glBegin(GL_LINES);
            dispV(c.p3);
            dispV(c.T3d + c.p3);
            glEnd();

            glColor3d(0,1,0);//N
            glBegin(GL_LINES);
            dispV(c.p3);
            dispV(c.N3d + c.p3);
            glEnd();

            glColor3d(0,0,1);//B
            glBegin(GL_LINES);
            dispV(c.p3);
            dispV(c.B3d + c.p3);
            glEnd();
        }else if(switchTNB == 2){
            glColor3d(1,0,0);//T
            glBegin(GL_LINES);
            dispV(c.p3);
            dispV(c.Td + c.p3);
            glEnd();

            glColor3d(0,1,0);//N
            glBegin(GL_LINES);
            dispV(c.p3);
            dispV(c.Nd + c.p3);
            glEnd();
            glColor3d(0,0,1);//B
            glBegin(GL_LINES);
            dispV(c.p3);
            dispV(c.Bd + c.p3);
            glEnd();
        }

    }
}

void GLWidget_3D::DrawMesh(bool isFront){
    glEnable(GL_CULL_FACE);
    if(isFront){
        glCullFace(GL_FRONT);
        glColor3d(0.9,0.9,0.9);
    }
    else {
        glCullFace(GL_BACK);
        glColor3d(0.6,0.6,0.6);
    }
    for(int i = 0; i < (int)TriMeshs.size(); i++){
        //if(left[i] == 1)glColor3d(0.8, 0.2, 0.2);
        //else glColor3d(0.2, 0.2, 0.8);
        glBegin(GL_POLYGON);
        for (auto& v : TriMeshs[i]) glVertex3d(v.x, v.y, v.z);
        glEnd();
    }
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glDisable(GL_CULL_FACE);
}

void GLWidget_3D::DrawMeshLines(){
    //std::vector<std::vector<glm::f64vec3>> TriMeshs;

    glLineWidth(1.0f);
    glPolygonOffset(1.0f, 1.f);
    glColor3d(0.f, 0.f, 0.f);

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
    arcballMode = (e->button() == Qt::LeftButton) ? ArcBallMode::translate : (e->button() == Qt::RightButton)? ArcBallMode::rotate: ArcBallMode::none;
    befPos = this->mapFromGlobal(QCursor::pos());
    update();
}

void GLWidget_3D::wheelEvent(QWheelEvent *we){
    double z = (we->angleDelta().y() > 0) ? 5 : -5;
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

    if(actionType == 1){//????????????
        TransX -=  0.3 * diff.x();
        TransY +=  0.3 * diff.y();
    }else if(actionType == 2){//??????
        angleX += diff.x();
        angleY += diff.y();
    }
    befPos = curPos;
    update();
}

glm::f64vec3 GLWidget_3D::getVec(float x, float y){
    QSize s = this->size();
    double shortSide = std::min(s.width(), s.height());
    glm::f64vec3 pt(2. * x / shortSide - 1.0, -2.0 * y / shortSide + 1.0, 0.0);
    // z???????????????
    const double xySquared = pt.x * pt.x + pt.y * pt.y;
    if (xySquared <= 1.0) pt.z = std::sqrt(1.0 - xySquared);// ????????????????????????z???????????????
    else pt = glm::normalize(pt); // ????????????????????????????????????????????????
    return pt;
}

void GLWidget_3D::updateRotate() {
    QPointF curPos = this->mapFromGlobal(QCursor::pos());
    // ????????????????????????????????????????????????????????????????????????
    const glm::f64vec3 u = getVec(befPos.x(), befPos.y());
    const glm::f64vec3 v = getVec(curPos.x(), curPos.y());

    // ???????????????????????????????????? (=?????????????????????????????????)
    const double angle = std::acos(std::max(-1.0, std::min(glm::dot(u, v), 1.0)));

    // ????????????????????????????????????
    const glm::f64vec3 rotAxis = glm::cross(u, v);

    // ????????????????????????????????????????????????????????????
    //const glm::f64mat4 c2wMat = glm::inverse(viewMat);

    // ?????????????????????????????????
    //const glm::f64vec3 rotAxisWorldSpace = glm::vec3(c2wMat * glm::vec4(rotAxis, 0.0f));

    // ?????????????????????
    //acRotMat = glm::rotate(((4.0 * angle), rotAxisWorldSpace) * acRotMat;
}
