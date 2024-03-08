#include "setrulings.h"
#include <QDebug>
using namespace MathTool;

Vertex::Vertex(Eigen::Vector3d _p, bool _deformed){
    p2_ori = p = _p;
    p3_ori = p3 = _p;
    deformed = _deformed;
}
Vertex::Vertex(Eigen::Vector3d _p2, Eigen::Vector3d _p3, bool _deformed){
    p2_ori = p = _p2;
    p3_ori = p3 = _p3;
    deformed = _deformed;
}

void Vertex::update(){
    p2_ori = p; p3_ori = p3;
}

std::shared_ptr<Vertex> Vertex::deepCopy(){ return std::make_shared<Vertex>(p, p3, deformed);}

std::shared_ptr<CrvPt_FL> CrvPt_FL::deepCopy(){return std::make_shared<CrvPt_FL>(p,p3,s);}

std::shared_ptr<Line> Line::deepCopy(){
    std::shared_ptr<Line> L = std::make_shared<Line>(o,v,et);
    L->color = color; L->IsCrossed = IsCrossed; L->hasCrossPoint = hasCrossPoint;
    return L;
}

std::shared_ptr<Vertex4d> Vertex4d::deepCopy(){ return std::make_shared<Vertex4d>(first->deepCopy(),second->deepCopy(),third->deepCopy());}

bool Line::is_on_line(Eigen::Vector3d p){ return MathTool::is_point_on_line(p, v->p, o->p);}

CRV::CRV(int _crvNum, int DivSize){
    curveNum = _crvNum;
    CurvePoints.resize(curveNum);
    for(int i = 0; i < DivSize-1;i++)Rulings.push_back(std::make_shared<Line>());
    isempty = true;
    curveType = CurveType::none;
    IsInsertNewPoint = false;
    InsertPointSegment = -1;
}

std::shared_ptr<CRV> CRV::deepCopy(){
    std::shared_ptr<CRV> c = std::make_shared<CRV>(curveNum, 1);
    c->isempty = isempty; c->curveType = curveType;
    c->IsInsertNewPoint = IsInsertNewPoint; c->InsertPointSegment = InsertPointSegment;
    c->CurvePoints = CurvePoints; c->Rulings.resize(Rulings.size());
    c->ControllPoints = ControllPoints;
    for(int i = 0; i < (int)Rulings.size(); i++){
        c->Rulings[i] = Rulings[i]->deepCopy();
    }
    return c;
}

bool CRV::eraseCtrlPt(int curveDimention, int crvPtNum){
    if((int)ControllPoints.size() > 0)ControllPoints.erase(ControllPoints.end() - 1);
    ControllPoints.shrink_to_fit();
    if(curveType == CurveType::bezier3){
        //return drawBezier(curveDimention, crvPtNum);
    }else if(curveType == CurveType::bsp3){
        return drawBspline(curveDimention, crvPtNum);
    }else if(curveType == CurveType::line){
        return drawLine();
    }
    return false;
}

//void CRV::movePt(glm::f64vec3 p, int ind){ControllPoints[ind] = p;}

int CRV::movePtIndex(Eigen::Vector3d& p, double& dist){
    int n = ControllPoints.size();
    dist = 5.0;
    int ind = -1;
    for(int i = 0; i < n; i++){
        double d = (p - ControllPoints[i]).norm();
        if(d < dist){
            dist = d;
            ind = i;
        }
    }
    return ind;

}

CurveType CRV::getCurveType(){return curveType;}
void CRV::setCurveType(CurveType n){curveType = n;}

bool CRV::drawBspline(int curveDimention,  int crvPtNum){
    //Rulings.clear();
    if((int)ControllPoints.size() <= curveDimention) return false;

    int knotSize = (int)ControllPoints.size() + curveDimention + 1;
    std::vector<double>T(knotSize);
    for(int j = 0; j < knotSize; j++)T[j] = (double)j/(double)knotSize;

    //端点をくっつける
    for(int j = 0;j< curveDimention + 1;j++){ T[j] = 0; T[knotSize - 1 - j] = 1;}
    double t;
    t = T[curveDimention];
    for(int i = 0; i < crvPtNum; i++){
        Eigen::Vector3d vec =  bspline(ControllPoints,t  + 1e-9, curveDimention, T);
        CurvePoints[i] = vec;
        t += (T[(int)ControllPoints.size()] - T[curveDimention])/(double)(crvPtNum);
    }
    isempty = false;
    return true;
}

bool CRV::drawBezier(int curveDimention, int crvPtNum){
    if((int)ControllPoints.size() < curveDimention){return false;}
    double t = 0.0;
    for(int n = 0; n < crvPtNum; n++){
        Eigen::Vector3d v = Eigen::Vector3d(0,0, 0);
        for (int i = 0; i < int(ControllPoints.size()); i++)v += MathTool::BernsteinBasisFunc(curveDimention, i, t) * ControllPoints[i];
        CurvePoints[n] = v;
        t += 1/(double)crvPtNum;
    }
    isempty = false;
    return true;
}

bool CRV::setPoint(const std::vector<std::shared_ptr<Vertex>>&outline, Eigen::Vector3d N, Eigen::Vector3d& cp, std::vector<Eigen::Vector3d>& P){
    Eigen::Vector3d N0 = N + cp, N1 = -N + cp;
    Eigen::Vector3d v, v2;
    std::vector<Eigen::Vector3d> crossPoint;
    P.clear();
    bool IsIntersected = false;
    int onum = outline.size();
    for(int i = 0; i < onum; i++){
        v = outline[i]->p;
        v2 = outline[(i + 1) % onum]->p;
        if (IsIntersect(v, v2, N0, N1, true)) {
            auto tmp = getIntersectionPoint(v, v2, N0, N1);
            if(std::find_if(crossPoint.begin(), crossPoint.end(), [&](Eigen::Vector3d c){return (c - tmp).norm() < 1e-9;}) == crossPoint.end()) crossPoint.push_back(tmp);
            IsIntersected = true;
        }
    }

    std::vector<double>dist;
    for(auto& p: crossPoint)dist.push_back((p - N1).norm());
    for(int i = 0; i < (int)crossPoint.size(); i++){
        for(int j = i + 1; j < (int)crossPoint.size(); j++){
            if(dist[j] < dist[i]){
                MathTool::swap(crossPoint[i], crossPoint[j]); MathTool::swap(dist[i], dist[j]);
            }
        }
    }
    Eigen::Vector3d minPt = N0;
    for(auto&p: crossPoint){
        minPt = ((p - cp).norm() < (minPt - cp).norm()) ? p : minPt;
    }
    auto itr = std::find(crossPoint.begin(), crossPoint.end(), minPt);
    if(itr == crossPoint.end())return IsIntersected;
    //int minInd = std::distance(crossPoint.begin(), itr);
    //minInd /= 2;
    //P.push_back(crossPoint[2 * minInd]); P.push_back(crossPoint[2 * minInd + 1]);//原因ここ?
    for(int i = 0; i < (int)crossPoint.size(); i++){
        P.push_back(crossPoint[i]);
    }
    std::sort(P.begin(), P.end(), [](auto& pl, auto& pr){return pl.y() < pr.y();});
    return IsIntersected;
}

void CRV::BezierRulings(std::shared_ptr<OUTLINE>& outline, int& DivSize, int crvPtNum){
    qDebug() << "you can't use now";
    return;
    double l = 1000; //適当に大きな値
    Eigen::Vector3d N, T;
    std::vector<Eigen::Vector3d> crossPoint;
    double t = 0.0;
    //DivSize = (DivSize > crvPtNum)? crvPtNum: DivSize;

    //int m = std::min(2 * DivSize - 3, crvPtNum);//2 * DivSize - 1 ... DivSize: the number of mesh , DivSize - 1: rulings(folding line)
    int m = DivSize - 1;
    int i = 0, rind = 0;
    std::vector<std::shared_ptr<Vertex>> vertices = outline->getVertices();

    while(i < m){
        T = 3.0 * (-ControllPoints[0] + 3.0 * ControllPoints[1] - 3.0 * ControllPoints[2] + ControllPoints[3]) * t * t
                + 6.0 * (ControllPoints[0] - 2.0 * ControllPoints[1] + ControllPoints[2]) * t + 3.0 * (-ControllPoints[0] + ControllPoints[1]);//一階微分(tについて)
        N = l * Eigen::Vector3d(-T.y(), T.x(), 0);
        //if (glm::dot(glm::normalized(N), glm::f64vec3{0,-1, 0}) < 0) N *= -1;
        int n = std::min(int(((double)i/(double)(m - 1)) * crvPtNum), crvPtNum - 1);
        bool IsIntersected = setPoint(vertices, N, CurvePoints[n], crossPoint);

        if(IsIntersected){
            for(int j = 0; j < (int)crossPoint.size(); j += 2){
                //Rulings[rind] = std::make_tuple(new Vertex(crossPoint[j]), new Vertex(crossPoint[j+1]));
                i++;
            }
        }
        t += 1.0 / (double)(std::min(DivSize, crvPtNum) - 1);
    }
    return;
}

void CRV::BsplineRulings(std::shared_ptr<OUTLINE>& outline, int& DivSize, int crvPtNum, int curveDimention){
    if((int)ControllPoints.size() <= curveDimention) return;
    double l = 1000;//適当に大きな値
    Eigen::Vector3d N, T;
    std::vector<Eigen::Vector3d> crossPoint;
    std::vector<std::vector<Eigen::Vector3d>> CrossPoints;
    double t = 0.0;
    int knotSize = (int)ControllPoints.size() + curveDimention + 1;
    std::vector<double>Knot(knotSize);
    for(int j = 0; j < knotSize; j++)Knot[j] = (double)j/(double)knotSize;
    //端点をくっつける
    for(int j = 0;j< curveDimention + 1;j++){ Knot[j] = 0; Knot[knotSize - 1 - j] = 1;}
    int i = 0;
    int sind = -1, eind = crvPtNum;
    t = Knot[curveDimention] + (Knot[(int)ControllPoints.size()] - Knot[curveDimention]) * (1.0 / (double)crvPtNum);
    while(i < crvPtNum){
        Eigen::Vector3d vec, vec2;
        vec = bspline(ControllPoints,t  + 1e-9, curveDimention, Knot);
        vec2 = bspline(ControllPoints,t  - 1e-9, curveDimention, Knot);
        T = (vec - vec2);
        T /= 2 * 1e-9;
        T = T.normalized();
        N = l * Eigen::Vector3d(-T.y(), T.x(), 0);
        if (N.dot(Eigen::Vector3d(0,1,0)) < 0) N *= -1;
        bool IsIntersected = setPoint(outline->vertices, N, CurvePoints[i], crossPoint);
        CrossPoints.push_back(crossPoint);
        if(sind == -1 && IsIntersected)sind = i;
        if(!IsIntersected && sind != -1)eind = std::min(eind, i);
        i++;
        t += (Knot[(int)ControllPoints.size()] - Knot[curveDimention]) * (1.0/ (double)crvPtNum);
    }
    double n = double(sind);
    int rind = 0;
    //ruling *r;
    while(n < eind){
        Rulings[rind]->et = EdgeType::r;
        Rulings[rind]->v = std::make_shared<Vertex>(Vertex(CrossPoints[floor(n)][0]));
        Rulings[rind]->o = std::make_shared<Vertex>(Vertex(CrossPoints[floor(n)][1]));
        rind++;
        n += (double)(eind - sind + 1)/(double)(DivSize - 1);
    }
    return;
}

//制御点の個数は2で固定
bool CRV::drawLine(){
    if((int)ControllPoints.size() < 2)return false;
    while(ControllPoints.size() != 2){
        ControllPoints.erase(ControllPoints.end() - 1);
        ControllPoints.shrink_to_fit();
    }
    Eigen::Vector3d V = ControllPoints[1] - ControllPoints[0];
    for(int i = 0; i < curveNum; i++){
        CurvePoints[i] = (double)i/(double)(curveNum - 1) * V + ControllPoints[0];
    }
    isempty = false;
    return true;
}

void CRV::LineRulings(std::shared_ptr<OUTLINE>& outline, int DivSize){
    if(ControllPoints.size() != 2)return;
    double l = 1000;//適当に大きな値

    std::vector<Eigen::Vector3d> crossPoint;
    std::vector<std::vector<Eigen::Vector3d>> CrossPoints;
    Eigen::Vector3d V = (ControllPoints[1] - ControllPoints[0]).normalized();
    Eigen::Vector3d N = l * Eigen::Vector3d(-V.y(), V.x(), 0);
    int i = 0;
    int sind = -1, eind = curveNum - 1;

    for(i = 0; i < curveNum; i++){
        bool IsIntersected = setPoint(outline->vertices, N, CurvePoints[i], crossPoint);
        CrossPoints.push_back(crossPoint);
        if(sind == -1 && IsIntersected)sind = i;
        if(!IsIntersected && sind != -1)eind = std::min(eind, i);
    }
    //bool RulingSet = false;
    double n = double(sind);
    int rind = 0;
    while(n < eind){
        Rulings[rind]->et = EdgeType::r;
        Rulings[rind]->v = std::make_shared<Vertex>(Vertex(CrossPoints[floor(n)][0]));
        Rulings[rind]->o = std::make_shared<Vertex>(Vertex(CrossPoints[floor(n)][1]));
        rind++;
        n += (double)(eind - sind + 1)/(double)(DivSize - 1);
    }
    return;
}

//制御点 0: 原点. 1,2 終点
bool CRV::drawArc(int crvPtNum){
    if((int)ControllPoints.size() < 3)return false;
    auto v1 = (ControllPoints[1] - ControllPoints[0]).normalized(), v2 = (ControllPoints[2] - ControllPoints[0]).normalized();
    double phi = std::acos(v1.dot(v2));
    Eigen::Vector3d axis = v1.cross(v2);
    Eigen::Translation3d T(ControllPoints[0]), invT(-ControllPoints[0]);// 平行移動ベクトルを作成
    for(int i = 0; i < crvPtNum; i++){
        Eigen::AngleAxisd R = Eigen::AngleAxisd((phi * (double)i/(double)crvPtNum), axis); // 回転行列を作成
        Eigen::Transform<double, 3, Eigen::Affine> transform = T * R * invT;// Transform行列を作成
        CurvePoints[i] = transform * ControllPoints[1];
    }
    isempty = false;
    return true;
}

void CRV::ArcRulings(std::shared_ptr<OUTLINE>& outline, int DivSize){
    if(ControllPoints.size() != 3)return;
    double l = 1000;//適当に大きな値
    std::vector<Eigen::Vector3d> crossPoint;
    std::vector<std::vector<Eigen::Vector3d>> CrossPoints;
    int sind = -1, eind = curveNum - 1;
    Eigen::Vector3d V, N;
    for(int i = 0; i < curveNum; i++){

        //if(i == 0)V = (CurvePoints[1] - CurvePoints[0]).normalized();
        //else V = (CurvePoints[i] - CurvePoints[i - 1]).normalized();
        N = l * Eigen::Vector3d(-V.y(), V.x(), 0);
        N = l * (ControllPoints[0] - CurvePoints[i]).normalized();
        bool IsIntersected = setPoint(outline->vertices, N, CurvePoints[i], crossPoint);
        CrossPoints.push_back(crossPoint);
        if(sind == -1 && IsIntersected)sind = i;
        if(!IsIntersected && sind != -1)eind = std::min(eind, i);
    }
    double n = double(sind);
    int rind = 0;
    while(n < eind){
        Rulings[rind]->et = EdgeType::r;
        Rulings[rind]->v = std::make_shared<Vertex>(Vertex(CrossPoints[floor(n)][0]));
        Rulings[rind]->o = std::make_shared<Vertex>(Vertex(CrossPoints[floor(n)][1]));
        rind++;
        n += (double)(eind - sind + 1)/(double)(DivSize - 1);
    }
    return;
}

void CRV::InsertControlPoint(QPointF _p){
    int ind, d;
    Eigen::Vector3d p(_p.x(), _p.y(), 0);
    //d = -1：どこにものっかっていない　d = 0：曲線上　d = 1：制御点を結んだ線上
    if(isempty)d = -1;
    double lc = 5, lp = 5;
    for(int i = 0; i < curveNum; i++){
        if((p - CurvePoints[i]).norm() < lc){
            lc = (p - CurvePoints[i]).norm(); ind = i;
        }
    }
    Eigen::Vector3d q;
    for(int i = 0; i < (int)ControllPoints.size() - 1; i++) {
        double l = distP2L(ControllPoints[i], ControllPoints[i + 1], p, q);
        if(l != -1 && l < lp){
            lp = l; ind = i;
        }
    }
    if(lc == 5 && lp == 5)d = -1;
    else if(lc < lp) d = 0;
    else d = 1;

    double dist;
    std::vector<int>CurveIndexs((int)ControllPoints.size(), -1);
    for(int j = 0; j < (int)ControllPoints.size(); j++){
        dist = 100;
        for(int i = 0; i < curveNum; i++){
            if(dist > (CurvePoints[i] - ControllPoints[j]).norm()){
                dist = (CurvePoints[i] - ControllPoints[j]).norm(); CurveIndexs[j] = i;
            }
        }
    }
    if(d == -1){
        double minDist = 200.0;
        double l;
        Eigen::Vector3d q;
        for(int i = 0; i < (int)ControllPoints.size() - 1; i++){
            l = distP2L(ControllPoints[i], ControllPoints[i + 1], p, q);
            if(l != -1 && l < minDist){ minDist = l; InsertPoint = p; InsertPointSegment = i;}
        }
        if(minDist == 100.0){IsInsertNewPoint = false;InsertPointSegment = -1;}
    }
    else if(d == 0){//曲線上
        for(int i = 0; i < (int)CurveIndexs.size() - 1; i++){
            if(ind >=CurveIndexs[i] && ind < CurveIndexs[i + 1]){
                double s = (double)(ind - CurveIndexs[i])/(double)(CurveIndexs[i + 1] - CurveIndexs[i]);
                InsertPoint = s * (ControllPoints[i + 1] - ControllPoints[i]) + ControllPoints[i];
                InsertPointSegment = i;
                return;
            }
        }
    }else if(d == 1){//直線上
        Eigen::Vector3d q;
        double l = distP2L(ControllPoints[ind], ControllPoints[ind + 1], p, q);
        if(l != -1){
            InsertPoint = q;
            InsertPointSegment = ind;
        }else{
            InsertPointSegment = -1;
        }
    }
    return;
}

void CRV::SetNewPoint(){
    if(InsertPointSegment == -1)return;

    ControllPoints.insert(ControllPoints.begin() + InsertPointSegment + 1, InsertPoint);
}

OUTLINE::OUTLINE(){
    type = "Rectangle";
    vertices.clear();
    Lines.clear();
    VerticesNum = 3;
    origin = Eigen::Vector2d(-1,-1);
    hasPtNum = 0;
}

std::shared_ptr<OUTLINE> OUTLINE::deepCopy(){
    std::shared_ptr<OUTLINE> o = std::make_shared<OUTLINE>();
    o->VerticesNum = VerticesNum; o->hasPtNum = hasPtNum; o->origin = origin; o->type = type;
    o->vertices.resize(vertices.size()); o->Lines.resize(Lines.size());
    for(int i = 0; i < (int)vertices.size(); i++){
        o->vertices[i] = vertices[i]->deepCopy();
    }
    for(int i = 0; i < (int)Lines.size(); i++){
        o->Lines[i] = Lines[i]->deepCopy();
    }
    return o;
}

void OUTLINE::addVertex(const std::shared_ptr<Vertex>& v, int n){
    if(n > (int)vertices.size())vertices.push_back(v);
    else vertices.insert(vertices.begin() + n, v);
}

void OUTLINE::addVertex(Eigen::Vector3d p){
    if(IsClosed()) return;

    if(type == "Rectangle"){
        vertices.push_back(std::make_shared<Vertex>(p));
        if((int)vertices.size() == 2){
            vertices.insert(vertices.begin() + 1, std::make_shared<Vertex>(Eigen::Vector3d(p.x(), vertices[0]->p.y(), 0)));
            vertices.push_back(std::make_shared<Vertex>(Eigen::Vector3d(vertices[0]->p.x(), p.y(), 0)));
            ConnectEdges();
            if(Eigen::Vector3d::UnitX().dot(getNormalVec()) > 0){
                //vertices[1]->p = vertices[1]->p3 = glm::f64vec3{p.x, vertices[0]->p.y, 0}; vertices[3]->p = vertices[3]->p3 = glm::f64vec3{vertices[0]->p.x, p.y, 0};
            }
        }
    }else if(type == "Polyline"){
        double d = 5;
        int ind = -1;
        for(int i = 0; i < (int)vertices.size(); i++){
            double dist = (p - vertices[i]->p).norm();
            if(dist < d){
                ind = i; d = dist;
            }
        }
        if(ind != -1){
            ConnectEdges();
            return;
        }
        else vertices.push_back(std::make_shared<Vertex>(p));
    }
}

void OUTLINE::eraseVertex(){
    if(type == "Polyline"){
        if(vertices.size() == 0)return;
        vertices.erase(vertices.end() - 1);
    }
    if(type == "Rectangle"){
        if(vertices.size() == 0)return;
        if(vertices.size() != 1){
            while(vertices.size() != 1){
                vertices.erase(vertices.end() - 1);
            }
        }else vertices.erase(vertices.begin());
    }
    if(type == "Polygon"){
        if(hasPtNum == 2){
            vertices.clear();
            hasPtNum = 1;
        }
        else if(hasPtNum == 1){
            origin = Eigen::Vector2d(-1,-1);
            hasPtNum = 0;
        }
    }
    vertices.shrink_to_fit();
    ConnectEdges(false);
}

void OUTLINE::drawPolygon(Eigen::Vector3d p, bool IsClicked){
    if(hasPtNum == 0 && IsClicked){
        origin = Eigen::Vector2d(p.x(), p.y());
        hasPtNum = 1;
        return;
    }
    if(hasPtNum == 1){
        if((p.x() - origin.x()) * (p.x() - origin.x()) + (p.y() - origin.y()) * (p.y() - origin.y()) < 16) return;
        vertices.clear();
        Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d invT = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();

        vertices.push_back(std::make_shared<Vertex>(p));
        double a = 2 * std::numbers::pi /VerticesNum;

        R(0,0) = R(1,1) = cos(a); R(0,1) = -sin(a); R(1,0) = sin(a);
        T(0,2) = origin.x(); T(1,2) = origin.y();
        invT(0,2) = -origin.x(); invT(1,2) = -origin.y();

        Eigen::Vector3d x;
        for(int n = 0; n < VerticesNum - 1; n++){
            x = Eigen::Vector3d(vertices[n]->p.x(), vertices[n]->p.y(), 1);
            x = T * R * invT * x; x.z() = 0.0;
            vertices.push_back(std::make_shared<Vertex>(x));
        }
        ConnectEdges();
        hasPtNum = 2;
    }

}

void OUTLINE::MoveOutline(Eigen::Vector3d p){
    double dist = 5;
    int ind = -1;
    Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
    Eigen::Vector3d x;
    if(type == "Polygon"){
        double d = (p - Eigen::Vector3d(origin.x(), origin.y(), 0)).norm();
        ind = movePointIndex(p);
        if(ind != -1){
            T(0,2) = p.x() - vertices[ind]->p.x(); T(1,2) = p.y() - vertices[ind]->p.y();
            for(auto& v: vertices){
                x = T * Eigen::Vector3d(v->p.x(), v->p.y(), 1);
                 v->p2_ori.x() = v->p.x() = x(0);  v->p2_ori.y() = v->p.y() = x(1);
            }
            x = T * Eigen::Vector3d(origin.x(), origin.y(), 1);
            origin(0) = x(0); origin(1) = x(1);
        }else if(d < dist){
            T(0,2) = p.x() - origin.x(); T(1,2) = p.y() - origin.y();
            for(auto& v: vertices){
                x = T * Eigen::Vector3d(v->p.x(), v->p.y(), 1);
                v->p2_ori.x() = v->p.x() = x(0);  v->p2_ori.y() = v->p.y() = x(1);
            }
            x = T * Eigen::Vector3d(origin.x(), origin.y(), 1);
            origin.x() = x(0); origin.y() = x(1);
        }
    }else {
        ind = movePointIndex(p);
        if(ind == -1)return;
        T(0,2) = p.x() - vertices[ind]->p.x(); T(1,2) = p.y() - vertices[ind]->p.y();
        for(auto& v: vertices){
            x = T * Eigen::Vector3d(v->p.x(), v->p.y(), 1);
            v->p2_ori.x() = v->p.x() = x(0);  v->p2_ori.y() = v->p.y() = x(1);
        }
    }
}

void OUTLINE::MoveVertex(Eigen::Vector3d p, int ind){
    if(ind < 0 || (int)vertices.size() < ind)return;
    vertices[ind]->p2_ori = vertices[ind]->p = p;
}

std::vector<std::shared_ptr<Vertex>> OUTLINE::getVertices(){return vertices;}

void OUTLINE::ConnectEdges(bool IsConnected){
    //if(!isClosed)return;
    if(vertices.size() < 2)return;
    Lines.clear();
    if(IsConnected){
        for(int i = 0; i < (int)vertices.size(); i++)Lines.push_back(std::make_shared<Line>(Line(vertices[i], vertices[(i+1) % (int)vertices.size()], EdgeType::ol)));
        IsConnected = true;
    }else{

    }
}

bool OUTLINE::IsPointInFace(Eigen::Vector3d p){
    if(!IsClosed())return false;
    int cnt = 0;
    double vt;
    for(auto&l: Lines){
        if(l->o->p.y() <= p.y() && l->v->p.y() > p.y()){
            vt = (p.y() - l->o->p.y()) / (l->v->p.y() - l->o->p.y());
            if(p.x() < (l->o->p.x() + vt * (l->v->p.x() - l->o->p.x())))cnt++;
        }else if(l->o->p.y() > p.y() && l->v->p.y() <= p.y()){
            vt = (p.y() - l->o->p.y())/(l->v->p.y() - l->o->p.y());
            if(p.x() < (l->o->p.x() + vt * (l->v->p.x() - l->o->p.x())))cnt--;
        }
    }
    if (cnt == 0) return false;
    else return true;
}


bool OUTLINE::IsClosed(){
    if(Lines.empty())return false;
    for(auto& l: Lines)if(l == nullptr)return false;
    return true;
}


void CrossDetection(std::shared_ptr<OUTLINE>& outline, std::shared_ptr<CRV>& crvs){
    if(crvs->Rulings.front()->v == nullptr)return;
    for(auto&r: crvs->Rulings)r->IsCrossed = -1;
    if(crvs->getCurveType() == CurveType::arc || crvs->getCurveType() == CurveType::line)return;
    for(int in = 0; in < (int)crvs->Rulings.size(); in ++){
        for(int inn = in+1; inn < (int)crvs->Rulings.size(); inn++){
            bool rs = IsIntersect(crvs->Rulings[in]->v->p, crvs->Rulings[in]->o->p, crvs->Rulings[inn]->v->p,crvs->Rulings[inn]->o->p, false);
            if(rs){
                auto p = getIntersectionPoint(crvs->Rulings[in]->v->p, crvs->Rulings[in]->o->p, crvs->Rulings[inn]->v->p,crvs->Rulings[inn]->o->p);
                bool PointOnLines = false;
                bool PointInFace = outline->IsPointInFace(p);
                for(auto&l: outline->Lines){if(l->is_on_line(p))PointOnLines = true;}
                if(PointInFace)crvs->Rulings[in]->IsCrossed = crvs->Rulings[inn]->IsCrossed = 0;
            }
        }
    }
}

Eigen::Vector3d OUTLINE::getNormalVec(){
    Eigen::Vector3d N(0,0,0);
    auto prev = Lines.end() - 1;
    for(auto cur = Lines.begin(); cur != Lines.end(); cur++){
        auto v = ((*cur)->v->p - (*cur)->o->p).normalized(), v2 = ((*prev)->v->p - (*prev)->o->p).normalized();
        if(abs(v.dot(v2)) < 1.0 - 1e-5){
            return (v.cross(v2)).normalized();
        }
        prev = cur;
    }
    return N;
}


int OUTLINE::movePointIndex(Eigen::Vector3d p){
    int n = vertices.size();
    double dist = 5.0;
    int ind = -1;
    for(int i = 0; i < n; i++){
        double d = (p - vertices[i]->p).norm();
        if(d < dist){
            dist = d;
            ind = i;
        }
    }
    return ind;
}


std::vector<double> BezierClipping(std::vector<Eigen::Vector3d>&CtrlPts, const std::shared_ptr<Vertex>& p, const std::shared_ptr<Vertex>& q, int dim){
    double a, b, c;
    if(p->p.x() <= q->p.x()){
        a = q->p.y() - p->p.y(), b = p->p.x() - q->p.x(), c = q->p.x() * p->p.y() - p->p.x() * q->p.y();
    }
    else{a = p->p.y() - q->p.y(), b = q->p.x() - p->p.x(), c = p->p.x() * q->p.y() - q->p.x() * p->p.y();}

    int n = CtrlPts.size();
    std::vector<Eigen::Vector3d> base(n);
    for(int i = 0; i < n; i++){
        double d = -(a * CtrlPts[i].x() + b * CtrlPts[i].y() + c)/sqrt(a*a + b*b);
        base[i] = Eigen::Vector3d((double)i/(double)(n - 1), d, 0);
    }
    std::vector<Eigen::Vector3d> current;
    std::copy(base.begin(), base.end(), std::back_inserter(current));
    std::array<Eigen::Vector3d, 2> _line{Eigen::Vector3d::Zero(), Eigen::Vector3d::UnitX()};
    auto res = _bezierclipping(base, current, _line, dim);

    return res;
}

std::vector<std::shared_ptr<Vertex>> ConvexHull_polygon(const std::vector<std::shared_ptr<Vertex>>& Q){
    std::vector<std::shared_ptr<Vertex>> P;
    if((int)Q.size() < 3)return Q;
    std::shared_ptr<Vertex> p_ml = Q[0];
    for(auto&p: Q){
        if(p_ml->p.y() > p->p.y())p_ml = p;
        else if(p_ml->p.y() == p->p.y() && p_ml->p.x() > p->p.x()) p_ml = p;
    }
    std::vector<std::pair<double, std::shared_ptr<Vertex>>>Args;
    for(auto&p: Q){
        if(p->p == p_ml->p)continue;
        double phi = atan2(p->p.y() - p_ml->p.y(), p->p.x() - p_ml->p.x());
        bool hasSameAngle = false;
        for(auto& X: Args){
            if(X.first == phi){
                if((X.second->p - p_ml->p).norm() < (p->p - p_ml->p).norm())X.second = p;
                hasSameAngle = true;
            }
        }
        if(!hasSameAngle)Args.push_back(std::make_pair(phi, p));
    }
    // compare only the first value
    std::sort(Args.begin(), Args.end(),[](auto const& x, auto const& y) {return x.first < y.first; });
    if((int)Args.size() >= 1)P.push_back(Args[Args.size() - 1].second);
    P.push_back(p_ml);
    std::shared_ptr<Vertex> top, next;
    for(int i = 0; i < (int)Args.size(); i++){
        do{
            top = P.back();
            P.pop_back();
            next = P.back();
            if(SignedArea(next->p, Args[i].second->p,top->p) <= 0){
                P.push_back(top);
                P.push_back(Args[i].second);
                break;
            }
        }while(1);
    }
    P.pop_back();
    return P;
}
std::vector<std::shared_ptr<Vertex>> SortPolygon(std::vector<std::shared_ptr<Vertex>>& polygon){
    std::vector<std::shared_ptr<Vertex>> poly_sort;
    if(polygon.size() > 3) poly_sort = ConvexHull_polygon(polygon);
    else poly_sort = polygon;
    if((int)poly_sort.size()< 3){return{};}
    auto N = (poly_sort[0]->p - poly_sort[1]->p).cross(poly_sort[2]->p - poly_sort[1]->p);
    if(N.dot(Eigen::Vector3d::UnitZ()) < 0) std::reverse(poly_sort.begin(), poly_sort.end());
    return poly_sort;
}

std::shared_ptr<Vertex> getClosestVertex(const std::shared_ptr<Vertex>& v, const std::shared_ptr<Vertex>& o,  const std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, bool SkipTrimedPoint){
    std::shared_ptr<Vertex> V_max = nullptr;
    double t_max = -1;
    for(auto&fc: FoldingCurve){
        if(SkipTrimedPoint && !fc->IsCalc)continue;
        if(((fc->second->p - v->p).norm() < 1e-7|| (fc->second->p - o->p).norm() < 1e-7))continue;
        if(MathTool::is_point_on_line(fc->second->p, v->p, o->p)){
            double t = (fc->second->p - o->p).norm()/(v->p - o->p).norm();
            if(t > t_max){
                t_max = t; V_max = fc->second;
            }
        }
    }
    return V_max;
}

std::vector<std::vector<std::shared_ptr<Vertex>>> MakeModel(const std::vector<std::shared_ptr<Line>>& Surface,
                                                            const std::vector<std::shared_ptr<Line>>& Rulings,  const std::vector<std::vector<std::shared_ptr<Vertex4d>>>& Creases){
    //firstの頂点(折曲線上の点)と輪郭の頂点が重なる場合輪郭の頂点を取り除く
    auto RemoveOverlappedVertex = [&](std::vector<std::shared_ptr<Vertex>>& mesh){
        int i = 0;
        while(i < (int)mesh.size()){
            bool isOverlapped = false;
            for(int j = i+1; j < (int)mesh.size() && !isOverlapped; j++){
                if((mesh[i]->p - mesh[j]->p).norm() < 1e-6){mesh.erase(mesh.begin() + j); isOverlapped = true;}
            }
            i = (isOverlapped)? 0: i+1;
        }
    };

    auto SplitPolygon = [&](std::vector<std::vector<std::shared_ptr<Vertex>>>& Polygons, const std::shared_ptr<Vertex>& o, const std::shared_ptr<Vertex>& v){//v:新たに挿入したいvertex, o:基本的にfirstを与える
        for(auto& Poly :Polygons){
            int IsOnLine_v = -1, ind_o = -1, ind_v = -1;
            for(int j = 0; j < (int)Poly.size(); j++){
                if(MathTool::is_point_on_line(v->p, Poly[j]->p, Poly[(j + 1) % (int)Poly.size()]->p))IsOnLine_v = j;
                if((o->p - Poly[j]->p).norm() < 1e-6)ind_o = j;
                if((v->p - Poly[j]->p).norm() < 1e-6)ind_v = j;
            }
            if(ind_v != -1 && ind_o != -1){
                int i_min = std::min(ind_v, ind_o), i_max = std::max(ind_v, ind_o);
                std::vector<std::shared_ptr<Vertex>> poly2(Poly.begin() + i_min, Poly.begin() + i_max +1);
                if((i_max + 1 - i_min) < 3 || ((int)Poly.size() - (i_max - i_min - 1)) < 3)return;
                Poly.erase(Poly.begin() + i_min + 1, Poly.begin() + i_max);
                Polygons.push_back(poly2);
                return;
            }
            if(IsOnLine_v == -1)continue;
            Poly.insert(Poly.begin() + IsOnLine_v + 1, v);
            int f_ind = -1, s_ind = -1;
            for(int i = 0; i < (int)Poly.size(); i++){
                if((Poly[i]->p - o->p).norm() < 1e-6)f_ind = i;
                if((Poly[i]->p - v->p).norm() < 1e-6)s_ind = i;
            }
            if(f_ind == -1 || s_ind == -1 || f_ind == s_ind)continue;
            int i_min = std::min(f_ind, s_ind), i_max = std::max(f_ind, s_ind);
            std::vector<std::shared_ptr<Vertex>> poly2(Poly.begin() + i_min, Poly.begin() + i_max +1);
            Poly.erase(Poly.begin() + i_min + 1, Poly.begin() + i_max);
            Polygons.push_back(poly2);
            return;
        }
    };
    
    std::vector<std::vector<std::shared_ptr<Vertex>>> Polygons;
    std::vector<std::shared_ptr<Vertex>> polygon;
    for(auto& l: Surface) polygon.push_back(l->v);
    Polygons.push_back(polygon);
    std::vector<std::vector<std::shared_ptr<Vertex4d>>> DrawCrvs;

    for(auto&FC: Creases){
        if((int)FC.size() <= 2)continue;
        std::vector<std::shared_ptr<Vertex4d>> DrawCrv;
        for(auto& v: FC){if(v->IsCalc)DrawCrv.push_back(v);}
        DrawCrvs.push_back(DrawCrv);
        Eigen::Vector3d CrvDir = (DrawCrv.back()->first->p - DrawCrv.front()->first->p).normalized();
        for(auto&P: Polygons){
            int ind_fr = -1, ind_bc = -1;
            for(int i = 0; i < (int)P.size(); i++){
                if(MathTool::is_point_on_line(DrawCrv.front()->first->p, P[i]->p, P[(i + 1) % (int)P.size()]->p))ind_fr = i;
                if(MathTool::is_point_on_line(DrawCrv.back()->first->p, P[i]->p, P[(i + 1)  % (int)P.size()]->p))ind_bc = i;
            }
            if(ind_fr != -1 && ind_bc != -1){
                std::vector<std::shared_ptr<Vertex>> InsertedV, InsertedV_inv;
                for(auto&v: DrawCrv){InsertedV.push_back(v->first);InsertedV_inv.insert(InsertedV_inv.begin(), v->first);}
                int i_min = std::min(ind_fr, ind_bc) + 1, i_max = std::max(ind_fr, ind_bc) + 1;
                Eigen::Vector3d Dir_prev;
                std::vector<std::shared_ptr<Vertex>> poly2;
                if(i_min == i_max){//折曲線の端点が同じ辺上にある
                    Dir_prev = (P[(ind_bc + 1) % (int)P.size()]->p - P[ind_bc]->p).normalized();
                    if(CrvDir.dot(Dir_prev) > 0){
                        P.insert(P.begin() + i_min, InsertedV.begin(), InsertedV.end());
                        poly2.insert(poly2.end(), InsertedV_inv.begin(), InsertedV_inv.end());
                    }else{
                        P.insert(P.begin() + i_min, InsertedV_inv.begin(), InsertedV_inv.end());
                        poly2.insert(poly2.end(), InsertedV.begin(), InsertedV.end());
                    }
                }else{
                    Dir_prev = (P[i_min]->p - P[(i_min - 1) % (int)P.size()]->p).normalized();
                    poly2 = {P.begin() + i_min, P.begin() + i_max};
                    P.erase(P.begin() + i_min, P.begin() + i_max);
                    if(Eigen::Vector3d(0,0,1).dot(CrvDir.cross(Dir_prev)) < 0){
                        P.insert(P.begin() + i_min, InsertedV.begin(), InsertedV.end());
                        poly2.insert(poly2.end(), InsertedV_inv.begin(), InsertedV_inv.end());
                    }else{
                        P.insert(P.begin() + i_min, InsertedV_inv.begin(), InsertedV_inv.end());
                        poly2.insert(poly2.end(), InsertedV.begin(), InsertedV.end());
                    }
                }
                RemoveOverlappedVertex(P); RemoveOverlappedVertex(poly2);
                Polygons.push_back(poly2);
                break;
            }
        }
    }
    for(auto&DC: DrawCrvs){
        for(auto itr = DC.begin() + 1; itr != DC.end() - 1; itr++){
            SplitPolygon(Polygons, (*itr)->first, (*itr)->second);
            SplitPolygon(Polygons, (*itr)->first, (*itr)->third);
        }
    }

    for(auto itr_r = Rulings.begin(); itr_r != Rulings.end(); itr_r++){
        if((*itr_r)->hasCrossPoint)continue;
        for(auto&P: Polygons){
            int vind = -1, oind = -1;
            for(int i = 0; i < (int)P.size(); i++){
                double dif = (P[i]->p - P[(i + 1) % (int)P.size()]->p).norm() - (P[i]->p - (*itr_r)->o->p).norm() + (P[(i + 1) % (int)P.size()]->p - (*itr_r)->o->p).norm();
                if(MathTool::is_point_on_line((*itr_r)->o->p, P[i]->p, P[(i + 1) % (int)P.size()]->p))oind = i;
                if(MathTool::is_point_on_line((*itr_r)->v->p, P[i]->p, P[(i + 1)  % (int)P.size()]->p))vind = i;
            }
            if(vind == -1 || oind == -1)continue;
            int i_min = std::min(vind, oind) + 1, i_max = std::max(vind, oind) + 1;
            std::vector<std::shared_ptr<Vertex>> poly2 = {P.begin() + i_min, P.begin() + i_max};
            P.erase(P.begin() + i_min, P.begin() + i_max);
            P.push_back((*itr_r)->o); P.push_back((*itr_r)->v); P = SortPolygon(P);
            poly2.push_back((*itr_r)->o); poly2.push_back((*itr_r)->v); poly2 = SortPolygon(poly2);
            Polygons.push_back(poly2);
        }
    }
    return Polygons;
}

QPointF SetOnGrid(QPointF& cursol, double gridsize, const QSize& S){
    int x = (int)cursol.x() % (int)gridsize, y = (int)cursol.y() % (int)gridsize;
    x = (cursol.x() - x + gridsize/2);
    y = (cursol.y() - y + gridsize/2);
    QPointF p; p.setX(x); p.setY(y);
    return p;
}

