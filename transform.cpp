#include "transform.h"

void AffinTrans(std::shared_ptr<FoldLine>& Crease, const QPointF& befPos, const QPointF& curPos){//回転と拡縮においてbefPosは軸として扱う
    if(Crease == nullptr)return;
    Eigen::Vector3d move{(curPos - befPos).x(), (curPos - befPos).y(), 0};
    for(auto&p: Crease->CtrlPts)p += move;
    Crease->a_flap = -1;
    Crease->setCurve(3);
}

//最も遠い位置にある制御点の二点の距離が閾値よりも小さい場合は縮小処理しない
void AffinScale(std::shared_ptr<FoldLine>& Crease, const QPointF& basePos, const QPointF& befPos, const QPointF& curPos){
    double fardist = 0,  mindist = 100;//mindistが閾値
    for(int i = 0; i < (int)Crease->CtrlPts.size(); i++){
        for(int j = 0; j < (int)Crease->CtrlPts.size(); j++)fardist = std::max(fardist, (Crease->CtrlPts[i] - Crease->CtrlPts[j]).norm());
    }

    Eigen::Vector3d a{(befPos - basePos).x(), (befPos - basePos).y(), 0}, b{(curPos - basePos).x(), (curPos - basePos).y(), 0}, basePt{basePos.x(), basePos.y(),0};
    double scale = b.norm()/a.norm();
    if(fardist < mindist && scale < 1.0)return;
    for(auto&p: Crease->CtrlPts)p = scale * (p - basePt) + basePt;
    Crease->setCurve(3);
    Crease->a_flap = -1;
}

void AffinRotate(std::shared_ptr<FoldLine>& Crease, const QPointF& basePos, const QPointF& befPos, const QPointF& curPos){
    Eigen::Vector3d basePt{befPos.x(), befPos.y(), 0}, b{befPos.x() - basePos.x(), befPos.y() - basePos.y(), 0}, c{curPos.x() - basePos.x(), curPos.y() - basePos.y(), 0};
    b = b.normalized(); c = c.normalized();
    double r = b.dot(c); r = (r < -1)? std::numbers::pi: (r > 1)? 0: std::acos(r);
    Eigen::Vector3d axis = (b.cross(c)).normalized();
    Eigen::AngleAxisd R = Eigen::AngleAxisd(r, axis);
    for(auto&p: Crease->CtrlPts)   p = R * (p - basePt) + basePt;

    Crease->setCurve(3);
    Crease->a_flap = -1;
}
