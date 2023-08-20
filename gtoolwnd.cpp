#include "gtoolwnd.hpp"

GToolWnd::GToolWnd(QWidget *parent) :
    QMainWindow(parent),
    gtw(new Ui::GToolWnd)
{
   gtw->setupUi(this); 
   //this->mode = "linear";
   //connect(gtw->Gview, &GradationWidget::ColorValueChanged,  this , &GToolWnd::ColorChanged);
   connect(this, &GToolWnd::updateCurvePath, gtw->Gview, &GradationWidget::rePaint);
   //curvePt = nullptr;
}

GToolWnd::~GToolWnd()
{
    delete gtw;
}

void GToolWnd::set(std::vector<Eigen::Vector2d>CurvePath){
    if(CurvePath.size() == 0)return;
    double left = CurvePath[0].x();
    double right = CurvePath[(int)CurvePath.size() - 1].x();
    double width = 500;
    double r = width/(right - left);
    for(auto& c: CurvePath){
        c.y() -= 255.0;
        c.y() *= -1;
        c.x() = (c.x() - left + 5) * r;
    }
    gtw->Gview->CurvePath = CurvePath;
    emit updateCurvePath();
    return;
}
