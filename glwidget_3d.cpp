#include "glwidget_3d.h"
using namespace MathTool;

GLWidget_3D::GLWidget_3D(QWidget *parent):QOpenGLWidget(parent)
{
    this->Vertices.clear();
    this->firstRotate = true;
    actionType = 0;
    eraseCtrlPt = eraseCrossPt = eraseVec = eraseCurve = false;
    drawEdgePlane = -1;
    IsEraseNonFoldEdge = false;
    Scale = 0.1;
}

GLWidget_3D::~GLWidget_3D(){}

void GLWidget_3D::initializeGL(){
    initializeOpenGLFunctions();
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    QSize s = this->size();
    TransX = 25.f, TransY = -16.f, TransZ = -100.f;
    angleY = 90;
    glViewport(s.width(),0,s.width(),s.height());
    glMatrixMode(GL_PROJECTION);
    perspective(30.0f, (float)s.width() / (float)s.height(), 0.1, 100.f);
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glEnable(GL_DEPTH_TEST);

    // アルファブレンディングを有効にする
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

}

void GLWidget_3D::EraseNonFoldEdge(bool state){
    IsEraseNonFoldEdge = state;
    update();
}

void GLWidget_3D::setVertices(const Lines Surface,  const Lines Rulings,  const FoldLine3d Creases){

    auto Planerity  = [](const std::vector<std::shared_ptr<Vertex>>& vertices, const Lines Poly_V)->double{
        if((int)vertices.size() == 3)return 0.0;
        else{
            std::vector<Eigen::Vector3d> QuadPlane;
            for(auto&v: vertices){
                if(std::find_if(Poly_V.begin(), Poly_V.end(),[&](const std::shared_ptr<Line>& L){return (L->v->p - v->p).norm() < 1e-9; }) == Poly_V.end())
                    QuadPlane.push_back(v->p3);
            }
            if((int)QuadPlane.size() <= 3)return 0;
            double l_avg = ((QuadPlane[0] - QuadPlane[2]).norm() + (QuadPlane[1] - QuadPlane[3]).norm())/2.0;
            double d;
            Eigen::Vector3d u1 = (QuadPlane[0] - QuadPlane[2]).normalized(), u2 = (QuadPlane[1]-  QuadPlane[3]).normalized();
            if((u1.cross(u2)).norm() < 1e-9){
                Eigen::Vector3d H = QuadPlane[3] + u2.dot(QuadPlane[1] - QuadPlane[3]) * u2;
                d = (H - QuadPlane[1]).norm();
            }else{
                d = ((u1.cross(u2)).dot(QuadPlane[2] - QuadPlane[3]))/(u1.cross(u2)).norm();
            }
            return d/l_avg;
        }
    };
    Eigen::Matrix3d Mirror = Eigen::Matrix3d::Identity();
    Mirror(1,1) = -1;

    Vertices.clear();
    FoldLineVertices.clear();

    std::vector<std::vector<std::shared_ptr<Vertex4d>>> FoldingCurves;
    for(auto&crease: Creases){
        if(crease->data != nullptr)FoldingCurves.push_back(crease->data->FoldingCurve);
    }
    std::vector<std::vector<std::shared_ptr<Vertex>>> Polygons = MakeModel(Surface, Rulings, FoldingCurves);
    for(auto& polygon: Polygons){
        auto p_sort = SortPolygon(polygon);
        std::vector<Eigen::Vector3d> vertices;
        for(auto& p: p_sort)vertices.push_back(Scale * (Mirror * p->p3));
        Vertices.push_back(vertices);
    }
    update();
}

inline void GLWidget_3D::dispV(Eigen::Vector3d p){
    glVertex3d(p.x(), p.y(), p.z());
}

void GLWidget_3D::paintGL(){
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    glLoadIdentity();
    glScaled(0.1, 0.1, 0.1);
    glTranslated(-TransX, -TransY, TransZ);
    glRotated(0.2 * angleX, 0.0, 1.0, 0.0);
    glRotated(0.2 * angleY, 1.0, 0.0, 0.0);


    DrawMeshLines();
    DrawMesh(true);
    DrawMesh(false);

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

    for(int i = 0; i < (int)Vertices.size(); i++){
        if(drawEdgePlane == i)continue;
        glBegin(GL_POLYGON);
        for (auto& v : Vertices[i]) {
            glVertex3d(v.x(), v.y(), v.z());
        }
        glEnd();
    }
    glDisable(GL_LIGHTING);
    glDisable(GL_LIGHT0);
    glDisable(GL_CULL_FACE);

    glPolygonOffset(0.f,0.f);
}

void GLWidget_3D::DrawMeshLines(){
    glColor3d(0.f, 0.f, 0.f);
    glLineWidth(1.0f);
    for(auto&V: Vertices){
        glBegin(GL_LINE_LOOP);
        for(auto&v: V)glVertex3d(v.x(), v.y(), v.z());
        glEnd();
    }
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
    if(e->key() == Qt::Key_C)eraseCtrlPt = !eraseCtrlPt;
    if(e->key() == Qt::Key_X)eraseCrossPt = !eraseCrossPt;
    if(e->key() == Qt::Key_D)eraseCurve = !eraseCurve;
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

Eigen::Vector3d GLWidget_3D::getVec(double x, double y){
    QSize s = this->size();
    double shortSide = std::min(s.width(), s.height());
    Eigen::Vector3d pt(2. * x / shortSide - 1.0, -2.0 * y / shortSide + 1.0, 0.0);
    // z座標の計算
    const double xySquared = pt.x() * pt.x() + pt.y() * pt.y();
    if (xySquared <= 1.0) pt.z() = std::sqrt(1.0 - xySquared);// 単位円の内側ならz座標を計算
    else pt = pt.normalized(); // 外側なら球の外枠上にあると考える
    return pt;
}

void GLWidget_3D::updateRotate() {
    QPointF curPos = this->mapFromGlobal(QCursor::pos());
    // マウスクリック位置をアークボール球上の座標に変換
    const Eigen::Vector3d u = getVec(befPos.x(), befPos.y());
    const Eigen::Vector3d v = getVec(curPos.x(), curPos.y());

    // カメラ座標における回転量 (=世界座標における回転量)
    const double angle = std::acos(std::max(-1.0, std::min(u.dot(v), 1.0)));

    // カメラ空間における回転軸
    const Eigen::Vector3d rotAxis = u.cross(v);

    // カメラ座標の情報を世界座標に変換する行列
    //const glm::f64mat4 c2wMat = glm::inverse(viewMat);

    // 世界座標における回転軸
    //const Eigen::Vector3d rotAxisWorldSpace = glm::vec3(c2wMat * glm::vec4(rotAxis, 0.0f));

    // 回転行列の更新
    //acRotMat = glm::rotate(((4.0 * angle), rotAxisWorldSpace) * acRotMat;
}

