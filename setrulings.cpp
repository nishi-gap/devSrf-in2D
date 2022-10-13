#include "setrulings.h"

Vertex::Vertex(glm::f64vec3 _p){
    p = _p;
    halfedge.clear();
    deformed = false;
}

HalfEdge::HalfEdge(Vertex *v){
    vertex = v;
    pair = nullptr;
    r = nullptr;
    v->halfedge.push_back(this);
}

crvpt::crvpt(int _ind, glm::f64vec3 _pt, int _color){
    pt = _pt;
    color = _color;
    ind = _ind;
}

ruling::ruling(Vertex *a, Vertex *b, crvpt *_pt) {
    IsCrossed = -1;
    r = std::make_tuple(a, b);
    Gradation = 0;
    hasGradPt = false;
    he[0] = he[1] = nullptr;
    pt = _pt;
}

ruling::ruling(){
    IsCrossed = -1;
    he[0] = he[1] = nullptr;
    Gradation = 0;
    hasGradPt = false;
}

meshline::meshline(){
    vertices[0] = new Vertex(NullVec);
    vertices[1] = new Vertex(NullVec);
    vertices4grad.clear();
    hasRulings[0] = nullptr;
    hasRulings[1] = nullptr;
    pt = nullptr;
}

meshline::meshline(Vertex *a, Vertex *b, crvpt *_pt){
    vertices4grad.clear();
    vertices[0] = a;
    vertices[1] = b;
    hasRulings[0] = nullptr;
    hasRulings[1] = nullptr;
    pt = _pt;
}

Face::Face(HalfEdge *_halfedge){
    //rulings.clear();
    halfedge = _halfedge;
    bend = false;
    hasGradPt = false;
    //Gradation = 0;
    //GradPt = false;
}

CRV::CRV(int _crvNum, int DivSize){

    curveNum = _crvNum;
    ControllPoints.clear();
    CurvePoints.clear();
    for(int i = 0; i < curveNum; i++)CurvePoints.push_back(crvpt(i));
    for(int i = 0; i < DivSize-1;i++){
        ruling *r = new ruling();
        Rulings.push_back(r);
    }
    //for(int i = 0; i < DivSize-2;i++){
        //meshline *l = new meshline();
        //meshLines.push_back(l);
    //}

    isempty = true;
    curveType = -1;
    IsInsertNewPoint = false;
    InsertPointSegment = -1;
}

void CRV::clearColor(){
    for(int i = 0; i < curveNum; i++)CurvePoints[i].color = 0;
}

void CRV::ClearPt(){
    for(auto& cp: CurvePoints)cp.pt = glm::f64vec3{-1,-1, -1};
    isempty = true;
}


void CRV::addCtrlPt(QPointF p){ControllPoints.push_back(glm::f64vec3{p.x(), p.y(), 0}); qDebug()<<"addCtrlPt called and size is " << ControllPoints.size(); }

void CRV::eraseCtrlPt(int curveDimention, int crvPtNum){
    if((int)ControllPoints.size() > 0)ControllPoints.erase(ControllPoints.end() - 1);
    ControllPoints.shrink_to_fit();
    if(curveType == 0){
        Bezier(curveDimention, crvPtNum);
    }else if(curveType == 1){
        Bspline(curveDimention, crvPtNum);
    }else if(curveType == 2){
        Line();
    }
}

void CRV::movePt(QPointF p, int ind){ControllPoints[ind] = glm::f64vec3{p.x(), p.y(), 0};}

int CRV::movePtIndex(QPointF p, double& dist){
    int n = ControllPoints.size();
    dist = 5.0;
    int ind = -1;
    for(int i = 0; i < n; i++){
        double d = QVdist(p,ControllPoints[i]);
        if(d < dist){
            dist = d;
            ind = i;
        }
    }
    return ind;

}

int CRV::getCurveType(){return curveType;}
void CRV::setCurveType(int n){curveType = n;}

double CRV::basis(int j, int k, double t, std::vector<double>& T){
    if(k == 0) return (T[j] <= t && t < T[j + 1]) ? 1.0: 0.0;
    double a = T[j + k] - T[j], a2 = T[j + k + 1] - T[j + 1];
    double b = 0, b2 = 0;
    if(a != 0) b = (t - T[j])/a;
    if(a2 != 0)b2 = (T[j + k + 1] - t)/a2;
    return b * basis(j,k-1,t,T) + b2 * basis(j+1,k-1,t,T);
}

void CRV::Bspline(int curveDimention,  int crvPtNum){
    //Rulings.clear();
    ClearPt();
    if((int)ControllPoints.size() <= curveDimention) return;

    int knotSize = (int)ControllPoints.size() + curveDimention + 1;
    std::vector<double>T(knotSize);
    for(int j = 0; j < knotSize; j++)T[j] = (double)j/(double)knotSize;

    //端点をくっつける
    for(int j = 0;j< curveDimention + 1;j++){ T[j] = 0; T[knotSize - 1 - j] = 1;}
    double t, b;
    t = T[curveDimention];
    for(int i = 0; i < crvPtNum; i++){
        glm::f64vec3 vec = glm::f64vec3{0,0,0};

        for(int j = 0; j < (int)ControllPoints.size();j++){
            b = basis(j,curveDimention,t,T); if(std::isnan(b)) {qDebug()<<"stop"; break;}
            vec += ControllPoints[j] *  b;
        }
        CurvePoints[i].pt = vec;
        t += (T[(int)ControllPoints.size()] - T[curveDimention])/(double)(crvPtNum);
    }
    return;
}

double CRV::factorial(int n){
    return (n > 1)? factorial(n - 1) * n : 1;
}

double CRV::cmb(int n , int i){
    return factorial(n) / (factorial(i) * factorial(n - i));
}

void CRV::Bezier(int curveDimention, int crvPtNum){

    if((int)ControllPoints.size() < curveDimention){return;}

    double t = 0.0;
    glm::f64vec3 v;
    for(int n = 0; n < crvPtNum; n++){
        v = {0.f,0.f, 0.f};
        for (int i = 0; i < int(ControllPoints.size()); i++) {
            v += cmb(curveDimention, i) * std::pow(t, i) * std::pow(1 - t, curveDimention - i) * ControllPoints[i];
        }
        CurvePoints[n].pt = v;
        t += 1/(double)crvPtNum;
    }
    return;
}

void CRV::swap(glm::f64vec3&a, glm::f64vec3& b){
    auto c = a;
    a = b; b = c;
    return;
}

glm::f64vec3 CRV::getIntersectionPoint(glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3, glm::f64vec3& p4){
    double det = (p1.x - p2.x) * (p4.y - p3.y) - (p4.x - p3.x) * (p1.y - p2.y);
    double t = ((p4.y - p3.y) * (p4.x - p2.x) + (p3.x - p4.x) * (p4.y - p2.y)) / det;

    return glm::f64vec3{ t * p1.x + (1.0 - t) * p2.x, t * p1.y + (1.0 - t) * p2.y, 0};
}

bool CRV::setPoint(std::vector<Vertex*>&outline, glm::f64vec3 N, glm::f64vec3& cp, std::vector<glm::f64vec3>& P){
    glm::f64vec3 N0 = N + cp, N1 = -N + cp;
    glm::f64vec3 v, v2;
    std::vector<glm::f64vec3> crossPoint;
    P.clear();
    bool IsIntersected = false;
    int onum = outline.size();
    for(int i = 0; i < onum; i++){
        v = outline[i]->p;
        v2 = outline[(i + 1) % onum]->p;
        if (IsIntersect(v, v2, N0, N1)) {
            crossPoint.push_back(getIntersectionPoint(v, v2, N0, N1));
            IsIntersected = true;
        }
    }

    std::vector<double>dist;
    for(auto& p: crossPoint)dist.push_back(glm::distance(p, N1));
    for(int i = 0; i < (int)crossPoint.size(); i++){
        for(int j = i + 1; j < (int)crossPoint.size(); j++){
            if(dist[j] < dist[i]){
                std::swap(crossPoint[i], crossPoint[j]); std::swap(dist[i], dist[j]);
            }
        }
    }
    glm::f64vec3 minPt = N0;
    for(auto&p: crossPoint){
        minPt = (glm::distance(p, cp) < glm::distance(minPt, cp)) ? p : minPt;
    }
    std::vector<glm::f64vec3>::iterator itr = std::find(crossPoint.begin(), crossPoint.end(), minPt);
    if(itr == crossPoint.end())return IsIntersected;
    int minInd = std::distance(crossPoint.begin(), itr);
    minInd /= 2;
    //P.push_back(crossPoint[2 * minInd]); P.push_back(crossPoint[2 * minInd + 1]);//原因ここ?
    for(int i = 0; i < (int)crossPoint.size(); i++){
        //if(i != minInd || i != minInd + 1)
            P.push_back(crossPoint[i]);
    }
    return IsIntersected;
}

void CRV::BezierRulings(OUTLINE *outline, int& DivSize, int crvPtNum){
    std::cout << "you can't use now"<<std::endl;
    return;
    double l = 1000; //適当に大きな値
    glm::f64vec3 N, T;
    std::vector<glm::f64vec3> crossPoint;
    double t = 0.0;
    //DivSize = (DivSize > crvPtNum)? crvPtNum: DivSize;

    //int m = std::min(2 * DivSize - 3, crvPtNum);//2 * DivSize - 1 ... DivSize: the number of mesh , DivSize - 1: rulings(folding line)
    int m = DivSize - 1;
    int i = 0, rind = 0, mind = 0;
    bool RulingSet = false;
    std::vector<Vertex*> vertices = outline->getVertices();

    while(i < m){
        T = 3.0 * (-ControllPoints[0] + 3.0 * ControllPoints[1] - 3.0 * ControllPoints[2] + ControllPoints[3]) * t * t
                + 6.0 * (ControllPoints[0] - 2.0 * ControllPoints[1] + ControllPoints[2]) * t + 3.0 * (-ControllPoints[0] + ControllPoints[1]);//一階微分(tについて)
        N = l * glm::f64vec3{ -T.y, T.x, 0 };
        //if (glm::dot(glm::normalize(N), glm::f64vec3{0,-1, 0}) < 0) N *= -1;
        int n = std::min(int(((double)i/(double)(m - 1)) * crvPtNum), crvPtNum - 1);
        bool IsIntersected = setPoint(vertices, N, CurvePoints[n].pt, crossPoint);

        if(IsIntersected){
            /*
            if(RulingSet){
                for(int j = 0; j < (int)crossPoint.size(); j += 2){
                    meshLines[mind]->vertices[0]->p = crossPoint[j];
                    meshLines[mind]->vertices[1]->p = crossPoint[j + 1];
                }
                RulingSet = false;
                mind++;
            }else{
                for(int j = 0; j < (int)crossPoint.size(); j += 2){
                    Rulings[rind]->r = std::make_tuple(new Vertex(crossPoint[j]), new Vertex(crossPoint[j+1]));
                    Rulings[rind]->pt = &CurvePoints[n];
                    //Rulings.push_back(new ruling(crossPoint[j],crossPoint[j + 1], &CurvePoints[n]));
                    i++;
                }
                RulingSet = true;
                rind++;
            }*/
            for(int j = 0; j < (int)crossPoint.size(); j += 2){
                Rulings[rind]->r = std::make_tuple(new Vertex(crossPoint[j]), new Vertex(crossPoint[j+1]));
                Rulings[rind]->pt = &CurvePoints[n];
                //Rulings.push_back(new ruling(crossPoint[j],crossPoint[j + 1], &CurvePoints[n]));
                i++;
            }
        }
        t += 1.0 / (double)(std::min(DivSize, crvPtNum) - 1);
    }
    return;
}

void CRV::BsplineRulings(OUTLINE *outline, int& DivSize, int crvPtNum, int curveDimention){
    if((int)ControllPoints.size() <= curveDimention) return;
    double l = 1000;//適当に大きな値
    glm::f64vec3 N, T;
    std::vector<glm::f64vec3> crossPoint;
    std::vector<std::vector<glm::f64vec3>> CrossPoints;
    double t = 0.0;
    std::vector<Vertex*> vertices = outline->getVertices();
    int knotSize = (int)ControllPoints.size() + curveDimention + 1;
    std::vector<double>Knot(knotSize);
    for(int j = 0; j < knotSize; j++)Knot[j] = (double)j/(double)knotSize;
    //端点をくっつける
    for(int j = 0;j< curveDimention + 1;j++){ Knot[j] = 0; Knot[knotSize - 1 - j] = 1;}
    double b, b2;
    int i = 0;
    int sind = -1, eind = crvPtNum;
    t = Knot[curveDimention] + (Knot[(int)ControllPoints.size()] - Knot[curveDimention]) * (1.0 / (double)crvPtNum);
    while(i < crvPtNum){
        glm::f64vec3 vec = glm::f64vec3{0,0,0}, vec2 = glm::f64vec3{0,0,0};
        for(int j = 0; j < (int)ControllPoints.size();j++){
            b = basis(j,curveDimention, t + FLT_EPSILON,Knot); if(std::isnan(b)) {qDebug()<<"stop"; return;}
            b2 = basis(j,curveDimention, t - FLT_EPSILON,Knot); if(std::isnan(b)) {qDebug()<<"stop"; return;}
            vec += ControllPoints[j] *  b;
            vec2 += ControllPoints[j] * b2;
        }
        T = (vec - vec2);
        T /= 2 * FLT_EPSILON;
        T = glm::normalize(T);
        N = l * glm::f64vec3{ -T.y, T.x, 0};
        if (glm::dot(glm::normalize(N), glm::f64vec3{0,1,0}) < 0) N *= -1;
        bool IsIntersected = setPoint(vertices, N, CurvePoints[i].pt, crossPoint);
        CrossPoints.push_back(crossPoint);
        if(sind == -1 && IsIntersected)sind = i;
        if(!IsIntersected && sind != -1)eind = std::min(eind, i);
        i++;
        t += (Knot[(int)ControllPoints.size()] - Knot[curveDimention]) * (1.0/ (double)crvPtNum);
    }

    bool RulingSet = false;
    double n = double(sind);
    int rind = 0, mind = 0;
    //ruling *r;
    std::vector<ruling*> bef;
    while(n < eind){
        /*
        if(RulingSet){
            for(int j = 0; j < (int)CrossPoints[floor(n)].size(); j += 2){
                meshLines[mind]->vertices[0]->p = CrossPoints[floor(n)][j];
                meshLines[mind]->vertices[1]->p = CrossPoints[floor(n)][j + 1];
                //std::vector<Vertex*> tmp {new Vertex(crossPoint[j]), new Vertex(crossPoint[j + 1])};
                //meshLines[mind]->vertices4grad.push_back(tmp);
            }
            meshLines[mind]->pt = &CurvePoints[floor(n)];
            meshLines[mind]->hasRulings[0] = Rulings[rind - 1];
            meshLines[mind]->hasRulings[1] = Rulings[rind];
            RulingSet = false;
            mind++;
        }else{
            Rulings[rind]->r = std::make_tuple(new Vertex(CrossPoints[floor(n)][0]), new Vertex(CrossPoints[floor(n)][1]));
            Rulings[rind]->pt = &CurvePoints[floor(n)];
            for(int j = 0; j < (int)CrossPoints[floor(n)].size(); j += 2){

            }
            rind++;
            RulingSet = true;
        }*/
        Rulings[rind]->r = std::make_tuple(new Vertex(CrossPoints[floor(n)][0]), new Vertex(CrossPoints[floor(n)][1]));
        Rulings[rind]->pt = &CurvePoints[floor(n)];
        for(int j = 0; j < (int)CrossPoints[floor(n)].size(); j += 2){

        }
        rind++;
        n += (double)(eind - sind + 1)/(double)(DivSize - 1);
    }
    isempty = false;
    return;
}

//制御点の個数は2で固定
void CRV::Line(){
    if(ControllPoints.size() < 2)return;
    while(ControllPoints.size() != 2){
        ControllPoints.erase(ControllPoints.end() - 1);
        ControllPoints.shrink_to_fit();
    }
    glm::f64vec3 V = ControllPoints[1] - ControllPoints[0];
    for(int i = 0; i < curveNum; i++){
        CurvePoints[i].pt = (double)i/(double)(curveNum - 1) * V + ControllPoints[0];
        CurvePoints[i].ind = i;
    }
    
}

void CRV::LineRulings(OUTLINE *outline, int DivSize){
    if(ControllPoints.size() != 2)return;
    double l = 1000;//適当に大きな値

    std::vector<glm::f64vec3> crossPoint;
    std::vector<std::vector<glm::f64vec3>> CrossPoints;
    glm::f64vec3 V = glm::normalize(ControllPoints[1] - ControllPoints[0]);
    glm::f64vec3 N = l * glm::f64vec3{-V.y, V.x, 0};
    int i = 0;
    int sind = -1, eind = curveNum - 1;
    std::vector<Vertex*> vertices = outline->getVertices();

    for(i = 0; i < curveNum; i++){
        bool IsIntersected = setPoint(vertices, N, CurvePoints[i].pt, crossPoint);
        CrossPoints.push_back(crossPoint);
        if(sind == -1 && IsIntersected)sind = i;
        if(!IsIntersected && sind != -1)eind = std::min(eind, i);
    }
    bool RulingSet = false;
    double n = double(sind);
    int mind = 0, rind = 0;
    while(n < eind){
        /*
        if(RulingSet){
            meshLines[mind]->vertices4grad.clear();
            for(int j = 0; j < (int)CrossPoints[floor(n)].size(); j += 2){
                meshLines[mind]->vertices[0]->p = CrossPoints[floor(n)][j];
                meshLines[mind]->vertices[1]->p = CrossPoints[floor(n)][j + 1];
                std::tuple<Vertex*, Vertex*> tmp = std::make_tuple(new Vertex(CrossPoints[floor(n)][j]), new Vertex(CrossPoints[floor(n)][j + 1]));
                meshLines[mind]->vertices4grad.push_back(tmp);
            }
            meshLines[mind]->pt = &CurvePoints[floor(n)];
            meshLines[mind]->hasRulings[0] = Rulings[rind - 1];
            meshLines[mind]->hasRulings[1] = Rulings[rind];
            RulingSet = false;
            mind++;
        }else{
            Rulings[rind]->r = std::make_tuple(new Vertex(CrossPoints[floor(n)][0]), new Vertex(CrossPoints[floor(n)][1]));
            Rulings[rind]->pt = &CurvePoints[floor(n)];

            rind++;
            RulingSet = true;
        }*/
        Rulings[rind]->r = std::make_tuple(new Vertex(CrossPoints[floor(n)][0]), new Vertex(CrossPoints[floor(n)][1]));
        Rulings[rind]->pt = &CurvePoints[floor(n)];

        rind++;
        n += (double)(eind - sind + 1)/(double)(DivSize - 1);
    }
    isempty = false;
    return;
}

double CRV::meshLength(int s){
    double l = 0.0;
    int ns = int(s * crvStep), ne = (int((s + 1) * crvStep) > curveNum) ? curveNum: int((s + 1) * crvStep);

    glm::f64vec3 bef = CurvePoints[ns].pt;
    for(int n = ns; n < ne; n++){
        l += glm::distance(CurvePoints[n].pt, bef);
        bef = CurvePoints[n].pt;
    }
    return l;
}

void CRV::InsertControlPoint2(QPointF pt){
    int ind;
    int d = OnCurvesORLines(pt, ind);
    double dist;
    std::vector<int>CurveIndexs((int)ControllPoints.size(), -1);
    for(int j = 0; j < (int)ControllPoints.size(); j++){
        dist = 100;
        for(int i = 0; i < curveNum; i++){
            if(dist > glm::distance(CurvePoints[i].pt, ControllPoints[j])){
                dist = glm::distance(CurvePoints[i].pt, ControllPoints[j]); CurveIndexs[j] = i;
            }
        }
    }
    if(d == -1){
        double minDist = 200.0;
        double l;
        glm::f64vec3 q;
        for(int i = 0; i < (int)ControllPoints.size() - 1; i++){
            l = distP2L(ControllPoints[i], ControllPoints[i + 1], pt, q);
            if(l != -1 && l < minDist){ minDist = l; InsertPoint = glm::f64vec3{pt.x(), pt.y(), 0}; InsertPointSegment = i;}
        }
        if(minDist == 100.0){IsInsertNewPoint = false;InsertPointSegment = -1;}
    }
    else if(d == 0){//曲線上
        for(int i = 0; i < (int)CurveIndexs.size() - 1; i++){
            if(CurveIndexs[i] <= ind && ind < CurveIndexs[i + 1]){
                double s = (double)(ind - CurveIndexs[i])/(double)(CurveIndexs[i + 1] - CurveIndexs[i]);
                InsertPoint = s * (ControllPoints[i + 1] - ControllPoints[i]) + ControllPoints[i];
                InsertPointSegment = i;
                return;
            }
        }
    }else if(d == 1){//直線上
        glm::f64vec3 q;
        double l = distP2L(ControllPoints[ind], ControllPoints[ind + 1], pt, q);
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

int CRV::OnCurvesORLines(QPointF& p, int& ind){
    if(CurvePoints[0].pt == NullVec)return -1;
    double lc = 5, lp = 5;
    for(int i = 0; i < curveNum; i++){
        if(QVdist(p,CurvePoints[i].pt) < lc){
            lc = QVdist(p,CurvePoints[i].pt); ind = i;
        }
    }
    glm::f64vec3 q;
    for(int i = 0; i < (int)ControllPoints.size() - 1; i++) {
        double l = distP2L(ControllPoints[i], ControllPoints[i + 1], p, q);
        if(l != -1 && l < lp){
            lp = l; ind = i;
        }
    }
    if(lc == 5 && lp == 5)return -1;
    else if(lc < lp) return 0;
    else return 1;
}

int CRV::getColor(int n){
    return CurvePoints[floor(n * crvStep)].color;
}

void CRV::setcolor(int ctype, int cval, int n){
    CurvePoints[floor(n * crvStep)].color = ctype * cval;
}
void CRV::Interpolation(int method, int g0, int g1){
    if(method == 0){ //linear
        LinearInterpolation(g0, g1);
    }
}

void CRV::LinearInterpolation(int g0, int g1){
    if(g0 == -1){qDebug() << "There is no specification regarding the range"; return;}
    double alpha, beta;
    double xdist = 0.0;

    beta = CurvePoints[floor(g0 * crvStep)].color;
    double l = 0.0;
    for(int i = g0 ; i <= g1 ; i++)l += meshLength(i);
    alpha = (double)(CurvePoints[floor(g1 * crvStep)].color - CurvePoints[floor(g0 * crvStep)].color)/l;
    glm::f64vec3 bef = CurvePoints[floor(g0 * crvStep)].pt;
    for(int n = floor(2 * g0 * crvStep); n < floor(2 * g1 * crvStep); n++){
        CurvePoints[n].color = alpha * xdist + beta;
        xdist += glm::distance(CurvePoints[n].pt, bef);
        bef = CurvePoints[n].pt;
    }
}


void CRV::CrossDetection(){
    if(isempty)return;
    for(auto&r:Rulings)r->IsCrossed = -1;
    for(int in = 0; in < (int)Rulings.size(); in ++){
        for(int inn = in+1; inn < (int)Rulings.size(); inn++){
            bool rs = IsIntersect(std::get<0>(Rulings[in]->r)->p, std::get<1>(Rulings[in]->r)->p, std::get<0>(Rulings[inn]->r)->p,std::get<1>(Rulings[inn]->r)->p);
            if(rs)Rulings[in]->IsCrossed = Rulings[inn]->IsCrossed = 0;
        }
    }
}

void CRV::FillColor(int c){
    for(auto& cp: CurvePoints) cp.color = c;
}

OUTLINE::OUTLINE(){
    type = "Rectangle";
    vertices.clear();
    isClosed = false;
    VerticesNum = 3;
    origin = glm::f64vec2{-1,-1};
    hasPtNum = 0;
}

void OUTLINE::addVertex(QPointF p){
    if(isClosed) return;
    if(type == "Rectangle"){
        vertices.push_back(new Vertex(glm::f64vec3{p.x(), p.y(), 0}));
        if((int)vertices.size() == 2){
            vertices.insert(vertices.begin() + 1, new Vertex(glm::f64vec3{vertices[0]->p.x, p.y(), 0}));
            vertices.push_back(new Vertex(glm::f64vec3{p.x(), vertices[0]->p.y, 0}));
            isClosed = true;
        }
    }else if(type == "Polyline"){
        double d = 5;
        int ind = -1;
        for(int i = 0; i < (int)vertices.size(); i++){
            double dist = QVdist(p, vertices[i]->p);
            if(dist < d){
                ind = i; d = dist;
            }
        }
        if(ind != -1) isClosed = true;
        else vertices.push_back(new Vertex(glm::f64vec3{p.x(), p.y(), 0}));
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
            origin = glm::f64vec2{-1,-1};
            hasPtNum = 0;
        }

    }
    vertices.shrink_to_fit();
    isClosed = false;
}

void OUTLINE::drawPolygon(QPointF p, bool IsClicked){
    if(hasPtNum == 0 && IsClicked){
        origin = glm::f64vec2{p.x(), p.y()};
        hasPtNum = 1;
        return;
    }
    if(hasPtNum == 1){
        if((p.x() - origin.x) * (p.x() - origin.x) + (p.y() - origin.y) * (p.y() - origin.y) < 16) return;
        vertices.clear();
        Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d invT = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();

        glm::f64vec2 v = glm::f64vec2{p.x(), p.y()};
        vertices.push_back(new Vertex(glm::f64vec3{v, 0}));
        double a = 2 * M_PI /VerticesNum;

        R(0,0) = R(1,1) = cos(a); R(0,1) = -sin(a); R(1,0) = sin(a);
        T(0,2) = origin.x; T(1,2) = origin.y;
        invT(0,2) = -origin.x; invT(1,2) = -origin.y;

        Eigen::Vector3d x;
        for(int n = 0; n < VerticesNum - 1; n++){
            x = Eigen::Vector3d(vertices[n]->p.x, vertices[n]->p.y, 1);
            x = T * R * invT * x;
            vertices.push_back(new Vertex(glm::f64vec3{x(0), x(1), 0}));
        }
        isClosed = true;
        hasPtNum = 2;
    }
}

void OUTLINE::MoveVertex(QPointF p){
    double dist = 5;
    int ind = -1;
    Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
    Eigen::Vector3d x;
    if(type == "Polygon"){
        double d = QVdist(p, origin);
        ind = movePointIndex(p,vertices);
        if(ind != -1){
            T(0,2) = p.x() - vertices[ind]->p.x; T(1,2) = p.y() - vertices[ind]->p.y;
            for(auto& v: vertices){
                x = T * Eigen::Vector3d(v->p.x, v->p.y, 1);
                v->p.x = x(0); v->p.y = x(1);
            }
            x = T * Eigen::Vector3d(origin.x, origin.y, 1);
            origin.x = x(0); origin.y = x(1);
        }else if(d < dist){
            T(0,2) = p.x() - origin.x; T(1,2) = p.y() - origin.y;
            for(auto& v: vertices){
                x = T * Eigen::Vector3d(v->p.x, v->p.y, 1);
                v->p.x = x(0); v->p.y = x(1);
            }
            x = T * Eigen::Vector3d(origin.x, origin.y, 1);
            origin.x = x(0); origin.y = x(1);
        }
    }else {
        ind = movePointIndex(p,vertices);
        if(ind == -1)return;
        T(0,2) = p.x() - vertices[ind]->p.x; T(1,2) = p.y() - vertices[ind]->p.y;
        for(auto& v: vertices){
            x = T * Eigen::Vector3d(v->p.x, v->p.y, 1);
            v->p.x = x(0); v->p.y = x(1);
        }
    }
}

void OUTLINE::EditVertex(QPointF p){

}

std::vector<Vertex*> OUTLINE::getVertices(){
    return vertices;
}

double distP2L(glm::f64vec3 la, glm::f64vec3 lb, QPointF pt, glm::f64vec3& q){
    glm::f64vec3 p{pt.x(), pt.y(), 0};

    double s = glm::dot((p - la), (lb - la))/glm::dot(lb - la, lb - la);
    if(0 <= s && s <= 1){
        q = la + s * (lb - la);
        return glm::distance(p,q);
    }
    return -1;
}

double QVdist(QPointF p, glm::f64vec2 v){
    return std::sqrt((p.x() - v.x) * (p.x() - v.x) + (p.y() - v.y) * (p.y() - v.y));
}

int movePointIndex(QPointF p, std::vector<Vertex*>& V){
    int n = V.size();
    double dist = 5.0;
    int ind = -1;
    for(int i = 0; i < n; i++){
        double d = QVdist(p,V[i]->p);
        if(d < dist){
            dist = d;
            ind = i;
        }
    }
    return ind;
}

/*
void CrossDetection(CRV& crvs){
    for(int in = 0; in < (int)crvs.Rulings.size(); in ++){
        for(int inn = in+1; inn < (int)crvs.Rulings.size(); inn++){
            bool rs = IsIntersect(std::get<0>(crvs.Rulings[in]->r)->p, std::get<1>(crvs.Rulings[in]->r)->p, std::get<0>(crvs.Rulings[inn]->r)->p,std::get<1>(crvs.Rulings[inn]->r)->p);
            if(rs){
                crvs.Rulings[in]->IsCrossed = crvs.Rulings[inn]->IsCrossed = 0;
            }
        }
    }
}
*/
bool cn(std::vector<glm::f64vec2> &V2, QPointF p) {
    int cnt = 0;
    auto V = V2;
    V.push_back(V[0]); std::reverse(V.begin(), V.end());
    double vt;
    for (int i = 0; i < (int)V.size() - 1; i++) {
        if (V[i].y <= p.y() && V[i + 1].y > p.y()) {
            vt = ((double)p.y() - V[i].y) / (V[i + 1].y - V[i].y);
            if (p.x() < (V[i].x + (vt * (V[i + 1].x - V[i].x)))) cnt++;
        }
        else if (V[i].y > p.y() && V[i + 1].y <= p.y()) {
            vt = (p.y() - V[i].y) / (V[i + 1].y - V[i].y);
            if (p.x() < (V[i].x + (vt * (V[i + 1].x - V[i].x)))) cnt--;
        }
    }
    if (cnt == 0) return false;
    else return true;
}

bool cn(Face *face, glm::f64vec3 p){
    int cnt = 0;

    double vt;
    HalfEdge *he = face->halfedge;
    do{
        if(he->vertex->p.y <= p.y && he->next->vertex->p.y > p.y){
            vt = (p.y - he->vertex->p.y) / (he->next->vertex->p.y - he->vertex->p.y);
            if(p.x < (he->vertex->p.x + vt * (he->next->vertex->p.x - he->vertex->p.x)))cnt++;
        }else if(he->vertex->p.y > p.y && he->next->vertex->p.y <= p.y){
            vt = (p.y - he->vertex->p.y)/(he->next->vertex->p.y - he->vertex->p.y);
            if(p.x < (he->vertex->p.x + vt * (he->next->vertex->p.x - he->vertex->p.x)))cnt--;
        }

    }while(he != face->halfedge);
    if (cnt == 0) return false;
    else return true;
}

double set3pt(glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3) {
    return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

bool IsIntersect(glm::f64vec3&p1, glm::f64vec3&p2, glm::f64vec3&p3, glm::f64vec3&p4){
    double t1 = set3pt(p1, p2, p3);
    double t2 = set3pt(p1, p2, p4);
    double t3 = set3pt(p3, p4, p1);
    double t4 = set3pt(p3, p4, p2);

    if (t1 * t2 < 0 && t3 * t4 < 0) return true;//交点を持つ
    return false;
}

void Triangulation(std::vector<glm::f64vec3>&input, std::vector<std::vector<glm::f64vec3>>&output){
    output.clear();
    std::vector<glm::f64vec3> Edges = input;
    //std::reverse(Edges.begin(), Edges.end());
    int n = Edges.size();
    while(Edges.size() >= 3){
        n = Edges.size();
        for(int i = 0; i < n; i++){
            int prev = (n + i - 1) % n;
            int next = (n + i + 1) % n;
            std::vector<glm::f64vec3> tri = {Edges[prev], Edges[i], Edges[next]};
            bool elimTriMesh = true;
            for(int j = 0; j < n - 3; j++){
                glm::f64vec3 p = Edges[(next + 1 + j) % n];
                bool check1 = hasPointInTriangle3D(p, tri);
                bool check2 = IsAngleLessThan180(tri[1], tri[0], tri[2]);
                if(check1 || !check2){
                    elimTriMesh = false;
                    break;
                }
            }
            if(elimTriMesh){
                Edges.erase(Edges.begin() + i);
                output.push_back(tri);
                break;
            }
        }

    }
}

bool IsAngleLessThan180(glm::f64vec3& o, glm::f64vec3& a, glm::f64vec3& b){
    glm::f64vec3 ao = glm::normalize(a - o);
    glm::f64vec3 bo = glm::normalize(b - o);

    glm::f64vec3 Normal = glm::cross(ao, bo);
    glm::f64vec3 BiNormal = glm::cross(ao, Normal);
    if(glm::dot(BiNormal, bo) < 0){return true;}
    return false;
}

bool hasPointInTriangle3D(glm::f64vec3 p, std::vector<glm::f64vec3>& V){
    double angle = 0.0;
    glm::f64vec3 v1, v2;
    int n = V.size();
    for(int i = 0; i < n; i++){
        v1 = glm::normalize(V[i] - p);
        v2 = glm::normalize(V[(i + 1) % n] - p);
        angle += std::acos(glm::dot(v1, v2));
    }
    if(angle >= 2 * M_PI - FLT_EPSILON)return true;
    return false;

}

bool hasPointInPolygon(glm::f64vec3 p, std::vector<glm::f64vec3>& V){
    //std::cout<<"ca't use hasPointInPolygon Function now"<<std::endl; exit(0);
    int n = V.size();
    if(n < 3){std::cout << "not enough polygon size" << std::endl; exit(0);}
    glm::f64vec3 normal = glm::normalize(glm::cross(V[0] - V[1], V[2] - V[1]));
    glm::f64vec3 v1, v2;
    double res = 0.0;
    double angle;
    for(int i = 0; i < n; i++){
        int prev = (n + i - 1) % n;
        v1 = glm::normalize(V[n] - p);
        v2 = glm::normalize(V[prev] - p);
        angle = glm::dot(v1, v2);
        glm::f64vec3 c = glm::cross(v1, v2);
        if(glm::dot(c, normal) < 0) angle *= -1;
        res += angle;
    }
    //res /= 360.f;

    return (abs(res) >= FLT_EPSILON) ? true : false;
}

std::vector<glm::f64vec3> TranslateGLMfromHE(Face *f){
    std::vector<glm::f64vec3> vertices;
    HalfEdge *he = f->halfedge;
    do{
        vertices.push_back(he->vertex->p);
        he = he->next;
    }while(he != f->halfedge);
    return vertices;
}

glm::f64vec3 GetCenter(std::vector<glm::f64vec3>& vertices){
    glm::f64vec3 center = glm::f64vec3{0,0,0};
    for(auto& v: vertices)center += v;
    center/= vertices.size();
    return center;
}
