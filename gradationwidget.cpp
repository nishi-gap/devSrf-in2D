#include "gradationwidget.h"

GradationWidget::GradationWidget(QWidget *parent):QOpenGLWidget(parent)
{
    IsClicked = -1;
    rad = 5.0;
    glm::f64vec2 v = {0,0};
    a = v;
    b = v;
    Ctype = Cval = "";
    paintMode = "linear";
    CurvePath.clear();
    ControllPoints.clear();
    BsplineDim = 3;
}

GradationWidget::~GradationWidget(){

}

void GradationWidget::initializeGL(){
    initializeOpenGLFunctions();
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    QSize s = this->size();
    glViewport(0,0,s.width(),s.height());
}

void GradationWidget::resizeGL(int width, int height){
    glViewport(0,0,width,height);
    glLoadIdentity();
    glOrtho(-0.5, (float)width -0.5, (float)height -0.5, -0.5, -1, 1);
}

void GradationWidget::paintGL(){
   glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
   QSize s = this->size();
   glLoadIdentity();
   glOrtho(-0.5, (float)s.width() -0.5, (float)s.height() -0.5, -0.5, -1, 1);
   DrawStdLine(s.width(),s.height()/2);
   /*
    if(paintMode == "linear"){
        glPointSize(5.0);
        glColor3d(0,0,0);
        glBegin(GL_POINTS);

        glVertex2d(a.x,a.y);
        glEnd();

        glBegin(GL_POINTS);
        glVertex2d(b.x, b.y);
        glEnd();

        glColor3d(0.3, 0.3,0.3);
        glBegin(GL_LINES);
        glVertex2d(a.x,a.y);
        glVertex2d(b.x, b.y);
        glEnd();

    }else if(paintMode == "SplineInterpolation"){
        DrawCtrlPt();
        DrawCrvPt();
    }else if(paintMode == "FreeCurve"){
        DrawCtrlPt();

    }else if(paintMode == "B-spline"){
        DrawCtrlPt();
        DrawCrvPt();

    }*/
   DrawCrvPt();
   update();
}
void GradationWidget::rePaint(){
    update();
}

void GradationWidget::DrawStdLine(int w, int h){
    glColor3d(0.f, 0.f,0.f);
    glBegin(GL_LINES);
    glVertex2d(0,h);
    glVertex2d(w,h);
    glEnd();
}
void GradationWidget::DrawCtrlPt(){
    glColor3d(0.f,0.f,0.f);
    glBegin(GL_POINTS);
    for(auto& cp: ControllPoints) glVertex2d(cp.x(), cp.y());
    glEnd();

    glColor3d(1.f, 0.f, 0.f);
    glBegin(GL_LINE_STRIP);
    for(auto&cp: ControllPoints) glVertex2d(cp.x(), cp.y());
    glEnd();
}
void GradationWidget::DrawCrvPt(){
    glColor3d(0,0,1);
    glBegin(GL_LINE_STRIP);
    for(auto& p: CurvePath) glVertex2d(p.x(), p.y());
    glEnd();
}

void GradationWidget::mousePressEvent(QMouseEvent *e){
    QPointF p = this->mapFromGlobal(QCursor::pos());

    if(paintMode == "linear"){
    }
    if(e->button() ==Qt::LeftButton){
        if((p.x() - a.x) * (p.x() - a.x) + (p.y() - a.y) * (p.y() - a.y) < rad * rad) IsClicked = 0;
        if((p.x() - b.x) * (p.x() - b.x) + (p.y() - b.y) * (p.y() - b.y) < rad * rad) IsClicked = 1;
    }

    else if(paintMode == "SplineInterpolation"){
        if(e->button()==Qt::RightButton){
            ControllPoints.clear();
            CurvePath.clear();
        }
    }
    else if(paintMode == "FreeCurve"){
        if(e->button() == Qt::RightButton){
            ControllPoints.clear();
            CurvePath.clear();
        }
    }

    update();
}

void GradationWidget::mouseMoveEvent(QMouseEvent *e){
    QPointF p = this->mapFromGlobal(QCursor::pos());
    QSize s = this->size();
    if(paintMode == "linear"){
        p.setY((p.y() < 0)? 0: (p.y() > s.height()) ? s.height() : p.y());
        if(IsClicked == 0)a.y() = p.y();
        if(IsClicked == 1)b.y() = p.y();

        //cm.colorupdate((double)s.height()/2 - a.y, (double)s.height()/2 - b.y, gradPoints);
        Cval = QString::number(abs((double)s.height()/2 - p.y()));
        Ctype = (p.y() <= (double)s.height()/2)? "Red": "Blue";

    }else if(paintMode == "SplineInterpolation"){

    }else if(paintMode == "FreeCurve"){
        ControllPoints.push_back(Eigen::Vector2d(p.x(), p.y()));
    }else if(paintMode == "B-spline"){

    }

    emit ColorValueChanged();
    update();
}

void GradationWidget::mouseReleaseEvent(QMouseEvent *e){
    IsClicked = -1;
    if(e->button() ==Qt::RightButton){
        if(ControllPoints.size() > 0)ControllPoints.erase(ControllPoints.end() - 1);
        if(paintMode == "linear"){}
        else if(paintMode == "SplineInterpolation"){
            if(ControllPoints.size() > 1) {
                SplineInterpolation(ControllPoints, CurvePath);
                emit ColorValueChanged();
            }
        }else if(paintMode == "FreeCurve"){
            if(ControllPoints.size() > 0) emit ColorValueChanged();
        }else if(paintMode == "B-spline"){
            if((int)ControllPoints.size() > BsplineDim){
                Bspline(BsplineDim, ControllPoints, CurvePath, 100);             
                if(CurvePath.size() > 1)emit ColorValueChanged();
            }
        }

      }else{
        QPointF p = this->mapFromGlobal(QCursor::pos());
        if(paintMode == "linear"){

        }else if(paintMode == "SplineInterpolation"){
            ControllPoints.push_back(Eigen::Vector2d(p.x(), p.y()));
            if(ControllPoints.size() > 1){
                SplineInterpolation(ControllPoints, CurvePath);
                emit ColorValueChanged();
            }
        }else if(paintMode == "FreeCurve"){
            if(ControllPoints.size() != 0) emit ColorValueChanged();

        }else if(paintMode == "B-spline"){
            ControllPoints.push_back(Eigen::Vector2d(p.x(), p.y()));
            if((int)ControllPoints.size() > BsplineDim){          
                Bspline(BsplineDim, ControllPoints, CurvePath, 100);

                emit ColorValueChanged();
            }
        }
    }

    update();
}

void GradationWidget::changeDrawMode(QString s){
    paintMode = s;
    CurvePath.clear();
    ControllPoints.clear();
    update();
}

double GradationWidget::basis(int j, int k, double t, std::vector<double>& T){
    if(k == 0) return (T[j] <= t && t < T[j + 1]) ? 1.0: 0.0;
    double a = T[j + k] - T[j], a2 = T[j + k + 1] - T[j + 1];
    double b = 0, b2 = 0;
    if(a != 0) b = (t - T[j])/a;
    if(a2 != 0)b2 = (T[j + k + 1] - t)/a2;
    return b * basis(j,k-1,t,T) + b2 * basis(j+1,k-1,t,T);
}

void GradationWidget::Bspline(int n, std::vector<Eigen::Vector2d>& P, std::vector<Eigen::Vector2d>& curvePt, int ptSize){
    if((int)P.size() < n) return;
    curvePt.clear();
    int knotSize = (int)P.size() + n + 1;
    std::vector<double>T(knotSize);
    for(int j = 0; j < knotSize; j++)T[j] = (double)j/(double)knotSize;

    //端点をくっつける
    for(int j = 0;j<n+1;j++){ T[j] = 0; T[knotSize - 1 - j] = 1;}
    double t, b;

    for(t = T[n]; t <= T[(int)P.size()]; t += 1.0/(double)(ptSize)){
        Eigen::Vector2d vec(0,0);

        for(int j = 0; j < (int)P.size();j++){
            b = basis(j,n,t,T); if(std::isnan(b)) {qDebug()<<"stop"; break;}
                vec += P[j] *  b;
        }
        curvePt.push_back(vec);
        //qDebug() <<QString::fromStdString(glm::to_string(vec)) << b;
    }

}

void GradationWidget::SplineInterpolation(std::vector<Eigen::Vector2d>& cp, std::vector<Eigen::Vector2d>& CurvePath){

    if(cp.size() < 2)return;
    int curveNum = 300;
    curveNum = ((int)cp.size() < curveNum) ? curveNum: (int)cp.size();
    CurvePath.clear();

    int N =(int)cp.size() - 1;
    Eigen::VectorXd v(N - 1);
    std::vector<double>h(N);
    for(int i = 1; i < (int)cp.size(); i++)h[i - 1] = cp[i].x() - cp[i - 1].x();
    for(int i = 1; i < N; i++){
        double a = (h[i] != 0) ? (cp[i + 1].y() - cp[i].y())/h[i] : 0, b = (h[i - 1] != 0) ? (cp[i].y() - cp[i - 1].y())/h[i - 1]: 0;
        v(i - 1) = 6 * (a - b);
    }
    Eigen::VectorXd u;
    if(N == 1){
        u = Eigen::VectorXd::Zero(2);
    }else if(N == 2){
        u = Eigen::VectorXd::Zero(3);
        double dx1, dx2, dx3;
        dx1 = cp[2].x() - cp[0].x(); dx2 = cp[2].x() - cp[1].x(); dx3 = cp[1].x() - cp[0].x();
        if(abs(dx1) < 1e-7) u(1) = 0;
        else{
            double a = (dx2 == 0) ? 0: (cp[2].y() - cp[1].y())/dx2;
            double b = (dx3 == 0) ? 0: (cp[1].y() - cp[0].y())/dx3;
            u(1) = 3 * (a - b)/dx1;
        }
    }
    else if(N > 2){
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N-1,N-1);
        for(int i = 0; i < N-1;i++){
            if(i == 0) A(0,1) = h[1];
            else if(i == N-2) A(N-2,N-3) = h[N-2];
            else{
                A(i,i-1) = h[i];  A(i,i+1) = h[i+1];
            }
            A(i,i) = 2 * (h[i] + h[i + 1]);
        }
        Eigen::VectorXd t = A.colPivHouseholderQr().solve(v);
        u = Eigen::VectorXd::Zero(N+1);
        for(int i = 0; i < (int)t.size(); i++)u(i+1) = t(i);
    }

    double x, a, b, c, d;
    for(double i = 0; i <= N; i += 1/(double)curveNum){

        int j = floor(i);
        double den = cp[j + 1].x - cp[j].x;
        a = (den != 0)? (u(j + 1) - u(j))/(6 * den) : 0;
        b = u(j)/2;
        c = - den * (2 * u(j) + u(j + 1))/6.0;
        c += (den != 0)? (cp[j + 1].y - cp[j].y)/den : 0;
        d = cp[j].y;
        x = (i - j) * den;
        double y = a * std::pow(x, 3) + b * std::pow(x, 2) + c * x + d;
        x += cp[j].x();
        //qDebug() << i << " " << j << " " << x << " " << y << " " <<  QString::fromStdString(glm::to_string(cp[j]));
        CurvePath.push_back(Eigen::Vector2d(x,y));
    }

}

void GradationWidget::GetColorsFromNewGradationMode(){

}
