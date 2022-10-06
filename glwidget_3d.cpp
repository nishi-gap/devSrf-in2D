#include "glwidget_3d.h"

GLWidget_3D::GLWidget_3D(QWidget *parent):QOpenGLWidget(parent)
{
    this->Vertices.clear();
    TriMeshs.clear();
    this->firstRotate = true;
    actionType = 0;
    center = glm::f64vec3{0,0,0};
}
GLWidget_3D::~GLWidget_3D(){

}

void GLWidget_3D::initializeGL(){
    initializeOpenGLFunctions();
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    QSize s = this->size();
    TransX = 0.f, TransY = 0.f, TransZ = -6.f;
    RotY = 90;
    glViewport(s.width(),0,s.width(),s.height());
    glEnable(GL_POLYGON_OFFSET_FILL);
    glEnable(GL_DEPTH_TEST);
}

void GLWidget_3D::setVertices(std::vector<std::vector<glm::f64vec3>>& _Vertices, glm::f64vec3 center){
    this->Vertices = _Vertices;
    std::vector<std::vector<glm::f64vec3>> TriMesh;
    TriMeshs.clear();
    for(auto& V : Vertices){
       // std::cout << "vertices " << V.size() << std::endl;
        Triangulation(V, TriMesh);
        for(auto&tri: TriMesh)TriMeshs.push_back(tri);
    }
    //std::cout <<"\n";
    TransX = 0.1 * center.x;
    TransY = 0.05 * center.y;
    this->center = 0.1 * center;
    update();
}

void GLWidget_3D::paintGL(){
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    QSize s = this->size();

    glLoadIdentity();
    perspective(30.0f, (float)s.width() / (float)s.height(), 1.f, 100.f);

    glTranslated(-TransX, -TransY, -7.f + TransZ);
    glRotated(0.2 * RotX, 0.0, 1.0, 0.0);
    glRotated(0.2 * RotY, 1.0, 0.0, 0.0);
    glScaled(0.1, 0.1, 0.1);

    DrawGrid();
    DrawMesh(true);
    DrawMesh(false);
    DrawMeshLines();

}

void GLWidget_3D::DrawMesh(bool isFront){
    //std::vector<std::vector<glm::f64vec3>> TriMeshs;
    glEnable(GL_CULL_FACE);
    if(isFront){
        glCullFace(GL_FRONT);
        glColor3d(0.8,0.8,0.8);
    }else{
        glCullFace(GL_BACK);
        glColor3d(0.8,0.8,0.8);
    }
    for(auto& tri: TriMeshs){
        glBegin(GL_TRIANGLES);
        for (auto& v : tri) glVertex3d(v.x, v.y, v.z);
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

void GLWidget_3D::RotationX(int _RotX){
    RotX = _RotX;
    update();
}

void GLWidget_3D::RotationY(int _RotY){
    RotY = _RotY;
    update();
}

void GLWidget_3D::ChangeTranslateX(int _X){
    TransX = 0.1 * (double)_X;
    update();
}

void GLWidget_3D::ChangeTranslateY(int _Y){
    TransY = 0.1 * (double)_Y;
    update();
}

void GLWidget_3D::ChangeTranslateZ(int _Z){
    TransZ = 0.1 * (double)_Z;
    update();
}

void GLWidget_3D::mousePressEvent(QMouseEvent *e){
    actionType = (e->button() == Qt::LeftButton)? 1 : (e->button() == Qt::RightButton)? 2: 0;
    befPos = this->mapFromGlobal(QCursor::pos());
}

void GLWidget_3D::wheelEvent(QWheelEvent *we){
    double z = (we->angleDelta().y() > 0) ? 0.1 : -0.1;
    TransZ += z;
    update();
}

void GLWidget_3D::mouseMoveEvent(QMouseEvent *e){

    QPointF curPos = this->mapFromGlobal(QCursor::pos());
    QPointF diff = curPos - befPos;

    if(actionType == 1){//平行移動
        TransX -= 0.01 * diff.x();
        TransY += 0.01 * diff.y();
    }else if(actionType == 2){//回転
        RotX += diff.x();
        RotY += diff.y();
    }
    befPos = curPos;
    update();
}
