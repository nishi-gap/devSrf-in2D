#include "setrulings.h"

Vertex::Vertex(glm::f64vec3 _p){
    p = _p;
    p3 = _p;
    halfedge.clear();
    deformed = false;
}

HalfEdge::HalfEdge(Vertex *v, EdgeType _type){
    vertex = v;
    pair = nullptr;
    prev = nullptr;
    next = nullptr;
    r = nullptr;
    edgetype = _type;
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

void CRV::movePt(glm::f64vec3 p, int ind){ControllPoints[ind] = p;}

int CRV::movePtIndex(glm::f64vec3& p, double& dist){
    int n = ControllPoints.size();
    dist = 5.0;
    int ind = -1;
    for(int i = 0; i < n; i++){
        double d = glm::distance(p,ControllPoints[i]);
        if(d < dist){
            dist = d;
            ind = i;
        }
    }
    return ind;

}

int CRV::getCurveType(){return curveType;}
void CRV::setCurveType(int n){curveType = n;}

double basis(int j, int k, double t, std::vector<double>& T){
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
    double t;
    t = T[curveDimention];
    for(int i = 0; i < crvPtNum; i++){
        glm::f64vec3 vec =  bspline(ControllPoints,t  + FLT_EPSILON,curveDimention, T);
        CurvePoints[i].pt = vec;
        t += (T[(int)ControllPoints.size()] - T[curveDimention])/(double)(crvPtNum);
    }
    return;
}

double factorial(int n){
    return (n > 1)? factorial(n - 1) * n : 1;
}

double cmb(int n , int i){
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
                //Rulings[rind]->pt = &CurvePoints[n];
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
    int i = 0;
    int sind = -1, eind = crvPtNum;
    t = Knot[curveDimention] + (Knot[(int)ControllPoints.size()] - Knot[curveDimention]) * (1.0 / (double)crvPtNum);
    while(i < crvPtNum){
        glm::f64vec3 vec, vec2;
        vec = bspline(ControllPoints,t  + FLT_EPSILON,curveDimention, Knot);
        vec2 = bspline(ControllPoints,t  - FLT_EPSILON,curveDimention, Knot);
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
        //Rulings[rind]->pt = &CurvePoints[floor(n)];
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
        //Rulings[rind]->pt = &CurvePoints[floor(n)];

        rind++;
        n += (double)(eind - sind + 1)/(double)(DivSize - 1);
    }
    isempty = false;
    return;
}

//制御点 0: 原点. 1,2 終点
void CRV::Arc(int crvPtNum){
    if(ControllPoints.size() < 3)return;
    double phi = acos(glm::dot(glm::normalize(ControllPoints[1] - ControllPoints[0]), glm::normalize(ControllPoints[2] - ControllPoints[0])));
    glm::f64vec3 axis = glm::cross(glm::normalize(ControllPoints[1] - ControllPoints[0]), glm::normalize(ControllPoints[2] - ControllPoints[0]));
    glm::f64mat4x4  T = glm::translate(ControllPoints[0]);
    glm::f64mat4x4 invT = glm::translate(-ControllPoints[0]);
    glm::f64mat4x4  R;
    for(int i = 0; i < crvPtNum; i++){
        R = glm::rotate(phi * (double)i/(double)crvPtNum, axis);
        CurvePoints[i].pt = T * R * invT * glm::f64vec4{ControllPoints[1],1};
    }

}

void CRV::ArcRulings(OUTLINE *outline, int DivSize){
    if(ControllPoints.size() != 3)return;
    double l = 1000;//適当に大きな値
    std::vector<glm::f64vec3> crossPoint;
    std::vector<std::vector<glm::f64vec3>> CrossPoints;
    int sind = -1, eind = curveNum - 1;
    std::vector<Vertex*> vertices = outline->getVertices();
    glm::f64vec3 V, N;
    for(int i = 0; i < curveNum; i++){
        if(i == 0)V = glm::normalize(CurvePoints[1].pt - CurvePoints[0].pt);
        else V = glm::normalize(CurvePoints[i].pt - CurvePoints[i - 1].pt);
        N = l * glm::f64vec3{-V.y, V.x, 0};
        bool IsIntersected = setPoint(vertices, N, CurvePoints[i].pt, crossPoint);
        CrossPoints.push_back(crossPoint);
        if(sind == -1 && IsIntersected)sind = i;
        if(!IsIntersected && sind != -1)eind = std::min(eind, i);
    }
    double n = double(sind);
    int rind = 0;
    while(n < eind){
        Rulings[rind]->r = std::make_tuple(new Vertex(CrossPoints[floor(n)][0]), new Vertex(CrossPoints[floor(n)][1]));
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

void CRV::InsertControlPoint2(glm::f64vec3& p){
    int ind;
    int d = OnCurvesORLines(p, ind);
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
            l = distP2L(ControllPoints[i], ControllPoints[i + 1], p, q);
            if(l != -1 && l < minDist){ minDist = l; InsertPoint = p; InsertPointSegment = i;}
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

int CRV::OnCurvesORLines(glm::f64vec3& p, int& ind){
    if(CurvePoints[0].pt == NullVec)return -1;
    double lc = 5, lp = 5;
    for(int i = 0; i < curveNum; i++){
        if(glm::distance(p,CurvePoints[i].pt) < lc){
            lc = glm::distance(p,CurvePoints[i].pt); ind = i;
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

void CRV::FillColor(int c){
    for(auto& cp: CurvePoints) cp.color = c;
}

OUTLINE::OUTLINE(){
    type = "Rectangle";
    vertices.clear();
    edges.clear();
    face = nullptr;
    VerticesNum = 3;
    origin = glm::f64vec2{-1,-1};
    hasPtNum = 0;
}

void OUTLINE::addVertex(Vertex*v, int n){
    if(n > (int)vertices.size())vertices.push_back(v);
    else vertices.insert(vertices.begin() + n, v);
}

void OUTLINE::addVertex(glm::f64vec3& p){
    if(IsClosed()) return;
    if(type == "Rectangle"){
        vertices.push_back(new Vertex(p));
        if((int)vertices.size() == 2){
            vertices.insert(vertices.begin() + 1, new Vertex(glm::f64vec3{vertices[0]->p.x, p.y, 0}));
            vertices.push_back(new Vertex(glm::f64vec3{p.x, vertices[0]->p.y, 0}));
            ConnectEdges();
        }
    }else if(type == "Polyline"){
        double d = 5;
        int ind = -1;
        for(int i = 0; i < (int)vertices.size(); i++){
            double dist = glm::distance(p, vertices[i]->p);
            if(dist < d){
                ind = i; d = dist;
            }
        }
        if(ind != -1){
            ConnectEdges();
            return;
        }
        else vertices.push_back(new Vertex(p));
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
    ConnectEdges(false);
}

void OUTLINE::drawPolygon(glm::f64vec3& p, bool IsClicked){
    if(hasPtNum == 0 && IsClicked){
        origin = glm::f64vec2{p.x, p.y};
        hasPtNum = 1;
        return;
    }
    if(hasPtNum == 1){
        if((p.x - origin.x) * (p.x - origin.x) + (p.y - origin.y) * (p.y - origin.y) < 16) return;
        vertices.clear();
        Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d invT = Eigen::Matrix3d::Identity();
        Eigen::Matrix3d R = Eigen::Matrix3d::Identity();

        vertices.push_back(new Vertex(p));
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
        ConnectEdges();
        hasPtNum = 2;
    }

}

void OUTLINE::MoveOutline(glm::f64vec3 p){
    double dist = 5;
    int ind = -1;
    Eigen::Matrix3d T = Eigen::Matrix3d::Identity();
    Eigen::Vector3d x;
    if(type == "Polygon"){
        double d = glm::distance(p, glm::f64vec3{origin, 0});
        ind = movePointIndex(p,vertices);
        if(ind != -1){
            T(0,2) = p.x - vertices[ind]->p.x; T(1,2) = p.y - vertices[ind]->p.y;
            for(auto& v: vertices){
                x = T * Eigen::Vector3d(v->p.x, v->p.y, 1);
                v->p.x = x(0); v->p.y = x(1);
            }
            x = T * Eigen::Vector3d(origin.x, origin.y, 1);
            origin.x = x(0); origin.y = x(1);
        }else if(d < dist){
            T(0,2) = p.x - origin.x; T(1,2) = p.y - origin.y;
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
        T(0,2) = p.x - vertices[ind]->p.x; T(1,2) = p.y - vertices[ind]->p.y;
        for(auto& v: vertices){
            x = T * Eigen::Vector3d(v->p.x, v->p.y, 1);
            v->p.x = x(0); v->p.y = x(1);
        }
    }
}

void OUTLINE::MoveVertex(glm::f64vec3 p, int ind){
    if(ind < 0 || vertices.size() < ind)return;
    vertices[ind]->p = p;
}

std::vector<Vertex*> OUTLINE::getVertices(){return vertices;}
std::vector<HalfEdge*> OUTLINE::getEdges(){return edges;}
Face* OUTLINE::getFace(){return face;}

void OUTLINE::ConnectEdges(bool IsConnected){
    //if(!isClosed)return;
    if(vertices.size() < 2)return;
    edges.clear();
    for(auto&v: vertices){
        HalfEdge *he = new HalfEdge(v, EdgeType::ol);
        edges.push_back(he);
    }
    if(IsConnected){
        face = new Face(edges[0]);
        for(int i = 0; i < (int)edges.size(); i++){
            edges[i]->next = edges[(i + 1) % edges.size()];
            edges[i]->prev = edges[(i - 1) % edges.size()];
            edges[i]->face = face;
        }
        //isClosed = true;
    }else{
        edges[0]->next = edges[1];
        for(int i = 1; i < (int)edges.size() - 1; i++){
            edges[i]->next = edges[(i + 1) % edges.size()];
            edges[i]->prev = edges[(i - 1) % edges.size()];
        }
        edges[edges.size() - 1]->prev = edges[edges.size() - 2];
    }

}


bool OUTLINE::IsClosed(){
    if(edges.empty())return false;
    HalfEdge *h = edges[0];
    do{
        if(h->next == nullptr)return false;
        h = h->next;
    }while(h != edges[0]);
    return true;
}

glm::f64vec3 bspline(std::vector<glm::f64vec3>&CtrlPts, double t, int dim, std::vector<double>Knot){
    glm::f64vec3 vec{0,0,0};
    for(int j = 0; j < (int)CtrlPts.size();j++){
        double b = basis(j,dim, t + FLT_EPSILON,Knot);
        vec += CtrlPts[j] *  b;
    }
    return vec;
}

double distP2L(glm::f64vec3 la, glm::f64vec3 lb, glm::f64vec3& p, glm::f64vec3& q){
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

int movePointIndex(glm::f64vec3 p, std::vector<Vertex*>& V){
    int n = V.size();
    double dist = 5.0;
    int ind = -1;
    for(int i = 0; i < n; i++){
        double d = glm::distance(p,V[i]->p);
        if(d < dist){
            dist = d;
            ind = i;
        }
    }
    return ind;
}

glm::f64vec3 getIntersectionPoint(glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3, glm::f64vec3& p4){
    double det = (p1.x - p2.x) * (p4.y - p3.y) - (p4.x - p3.x) * (p1.y - p2.y);
    double t = ((p4.y - p3.y) * (p4.x - p2.x) + (p3.x - p4.x) * (p4.y - p2.y)) / det;

    return glm::f64vec3{ t * p1.x + (1.0 - t) * p2.x, t * p1.y + (1.0 - t) * p2.y, 0};
}

void CrossDetection(Face *f, CRV *crvs){
    if(crvs->isempty)return;
    for(auto&r: crvs->Rulings)r->IsCrossed = -1;
    for(int in = 0; in < (int)crvs->Rulings.size(); in ++){
        for(int inn = in+1; inn < (int)crvs->Rulings.size(); inn++){
            bool rs = IsIntersect(std::get<0>(crvs->Rulings[in]->r)->p, std::get<1>(crvs->Rulings[in]->r)->p, std::get<0>(crvs->Rulings[inn]->r)->p,std::get<1>(crvs->Rulings[inn]->r)->p);
            if(rs){
                auto p = getIntersectionPoint(std::get<0>(crvs->Rulings[in]->r)->p, std::get<1>(crvs->Rulings[in]->r)->p, std::get<0>(crvs->Rulings[inn]->r)->p,std::get<1>(crvs->Rulings[inn]->r)->p);
                if(cn(f,p)) crvs->Rulings[in]->IsCrossed = crvs->Rulings[inn]->IsCrossed = 0;
            }
        }
    }
}

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
        he = he->next;
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

void Triangulation(std::vector<glm::f64vec3>&input, std::vector<std::array<glm::f64vec3, 3>>&output){
    output.clear();
    std::vector<glm::f64vec3> Edges;
    std::copy(input.begin(), input.end(), back_inserter(Edges) );
    int n = Edges.size();
    while(Edges.size() >= 3){
        n = Edges.size();
        for(int i = 0; i < n; i++){
            int prev = (n + i - 1) % n;
            int next = (n + i + 1) % n;
            std::array<glm::f64vec3, 3> tri = {Edges[prev], Edges[i], Edges[next]};
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

void Triangulation(Face *f, std::vector<std::array<glm::f64vec3, 3>>& output){
    std::vector<HalfEdge*> Edges;
    HalfEdge *h = f->halfedge;
    do{
        Edges.push_back(h);
        h = h->next;
    }while(h != f->halfedge);

    int n = Edges.size();
    while(Edges.size() >= 3){
        n = Edges.size();
        for(int i = 0; i < n; i++){
            std::array<HalfEdge*, 3> tri = {Edges[i]->prev, Edges[i], Edges[i]->next};
            bool elimTriMesh = true;
            for(int j = 0; j < n - 3; j++){
                glm::f64vec3 p = Edges[(i + 2 + j) % n]->vertex->p3;
                bool check1 = hasPointInTriangle3D(p, tri);
                bool check2 = IsAngleLessThan180(tri[1]->vertex->p3, tri[0]->vertex->p3, tri[2]->vertex->p3);
                if(check1 || !check2){
                    elimTriMesh = false;
                    break;
                }
            }
            if(elimTriMesh){
                Edges.erase(Edges.begin() + i);
                output.push_back(std::array{tri[0]->vertex->p3, tri[1]->vertex->p3, tri[2]->vertex->p3});
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

bool hasPointInTriangle3D(glm::f64vec3 p, std::array<glm::f64vec3, 3>& V){
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

bool hasPointInTriangle3D(glm::f64vec3 p, std::array<HalfEdge*, 3>& V){
    double angle = 0.0;
    glm::f64vec3 v1, v2;
    int n = V.size();
    for(int i = 0; i < n; i++){
        v1 = glm::normalize(V[i]->vertex->p - p);
        v2 = glm::normalize(V[i]->next->vertex->p - p);
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

std::vector<glm::f64vec3> ConvertDistBasedBezier(std::vector<glm::f64vec3>& CtrlPts, HalfEdge *line){
    glm::f64vec3 p, q;
    if(line->vertex->p.x <= line->next->vertex->p.x){p = line->next->vertex->p; q = line->vertex->p;}
    else{q = line->next->vertex->p; p = line->vertex->p;}
    double a = p.y - q.y, b = q.x - p.x, c = p.x * q.y - q.x * p.y;
    int n = CtrlPts.size();
    std::vector<glm::f64vec3> D(n);
    for(int i = 0; i < n; i++){
        double d = -(a * CtrlPts[i].x + b * CtrlPts[i].y + c)/sqrt(a*a + b*b);
        D[i] = glm::f64vec3{(double)i/(double)(n - 1), d, 0};
    }
    return D;
}

//交点が一つのみの場合
std::vector<double> BezierClipping(std::vector<glm::f64vec3>&CtrlPts, HalfEdge *line, int dim){
    std::vector<glm::f64vec3> base = ConvertDistBasedBezier(CtrlPts, line);
    std::vector<glm::f64vec3> current;
    std::copy(base.begin(), base.end(), std::back_inserter(current));
    std::array<glm::f64vec3, 2> _line{glm::f64vec3{0.,0,0}, glm::f64vec3{1.,0,0}};
    auto res = _bezierclipping(base, current, _line, dim);
    //std::cout<<"res  ";
    //for(auto&t: res)std::cout<< t << " , ";
    //std::cout<<std::endl;
    return res;
}

std::vector<double> _bezierclipping(std::vector<glm::f64vec3>&CtrlPts_base, std::vector<glm::f64vec3>&CtrlPts_cur, std::array<glm::f64vec3, 2>& line, int dim){
    double t_min = 1, t_max = 0;
    glm::f64vec3 p = line[0], q = line[1];   
    std::vector<glm::f64vec3> D = GrahamScan(CtrlPts_cur);
    //std::cout<< "line " << glm::to_string(p) << " , " << glm::to_string(q) << std::endl;
    int n = D.size();
    std::vector<double> T;
    for(int i = 0; i < n - 1; i++){
        glm::f64vec3 v = D[i], v2 = D[i + 1];

        double xp = ((p.y * q.x - p.x * q.y)*(v2.x - v.x) - (v.y * v2.x - v.x * v2.y)*(q.x - p.x)) / ((v2.y - v.y)*(q.x - p.x) - (v2.x - v.x)*(q.y - p.y));
        double yp = ((p.y * q.x - p.x * q.y)*(v2.y - v.y) - (v.y * v2.x - v.x * v2.y)*(q.y - p.y)) / ((v2.y - v.y)*(q.x - p.x) - (v2.x - v.x)*(q.y - p.y));
        glm::f64vec3 pt_new{xp, yp, 0};
        //std::cout<<glm::to_string(pt_new) << " , " << glm::to_string(v) << " , " << glm::to_string(v2) << std::endl;
        if(is_point_on_line(pt_new, p, q) && is_point_on_line(pt_new, v, v2)){ T.push_back(xp); }
    }
    static int cnt = 0;
    if(T.size() == 0){
        //std::cout<<"not found"<<std::endl;
        return {};
    }
    if(T.size() == 1){
        //std::cout<<"midpoint " << (p.x + q.x)/2 << std::endl;
        return{(p.x + q.x)/2};
    }
    t_min = *std::min_element(T.begin(), T.end()); t_min = (t_min < 0) ? 0: (t_min > 1)? 1: t_min;
    t_max = *std::max_element(T.begin(), T.end()); t_max = (t_max < 0) ? 0: (t_max > 1) ? 1:  t_max;
    std::array<glm::f64vec3, 2> next_line = std::array{glm::f64vec3{t_min, 0,0}, glm::f64vec3{t_max, 0,0}};

    if(abs(t_max - t_min) < 1e-7){
        //std::cout<<"success " << t_min << std::endl;
        return {(t_max + t_min)/2};
    }
    //std::cout<<"cnt "<<cnt++ << "  , " << T.size() << " min " << t_min << " , t_max " << t_max <<  std::endl;
    std::pair<std::vector<glm::f64vec3>, std::vector<glm::f64vec3>> _bez = BezierSplit(CtrlPts_base, t_max, dim);
    double bez_t = t_min / (t_max);
    //for(auto&b: _bez.first)std::cout<<"bez " << glm::to_string(b)<<std::endl;
    _bez = BezierSplit(_bez.first, bez_t, dim);
    //for(auto&b: _bez.second)std::cout<<"bez_next " << glm::to_string(b)<<std::endl;
    //std::cout<<"\n";
    if(abs(glm::distance(p,q) - abs(t_max - t_min)) < 1e-7){

        std::pair<std::vector<glm::f64vec3>, std::vector<glm::f64vec3>> bez_spl = BezierSplit(_bez.second, 0.5,  dim);
        std::vector<glm::f64vec3> b1 = bez_spl.first, b2 = bez_spl.second;
        std::array<glm::f64vec3, 2> next_line2; std::copy(next_line.begin(), next_line.end(), next_line2.begin());
        //for(auto&b: b1)std::cout<<"b1 " << glm::to_string(b)<<std::endl;
        //for(auto&b: b2)std::cout<<"b2 " << glm::to_string(b)<<std::endl;
        std::vector<double>t1 = _bezierclipping(CtrlPts_base, b1, next_line, dim), t2 = _bezierclipping(CtrlPts_base, b2, next_line2, dim);
        for(auto&t: t2)t1.push_back(t);
        return t1;
    }

    return _bezierclipping(CtrlPts_base, _bez.second, next_line, dim);
}

//http://www-ikn.ist.hokudai.ac.jp/~k-sekine/slides/convexhull.pdf
//https://kajindowsxp.com/graham-algo/
std::vector<glm::f64vec3> GrahamScan(std::vector<glm::f64vec3>& Q){
    std::vector<glm::f64vec3> S;
    if(Q.size() < 3)return S;
    glm::f64vec3 p_ml = Q[0];
    for(auto&p: Q){
        if(p_ml.y > p.y)p_ml = p;
        else if(p_ml.y == p.y && p_ml.x > p.x) p_ml = p;
    }
    std::vector<std::pair<double, glm::f64vec3>>Args;
    for(auto&p: Q){
        if(p == p_ml)continue;
        double phi = atan2(p.y - p_ml.y, p.x - p_ml.x);
        //if(phi < 0) phi = 2 * M_PI + phi;
        bool hasSameAngle = false;
        for(auto& X: Args){
            if(X.first == phi){
                if(glm::distance(X.second, p_ml) < glm::distance(p, p_ml))X.second = p;
                hasSameAngle = true;
            }
        }
        if(!hasSameAngle)Args.push_back(std::make_pair(phi, p));
    }
    // compare only the first value
    std::sort(Args.begin(), Args.end(),[](auto const& x, auto const& y) {return x.first < y.first; });
    auto cross3 = [](glm::f64vec3 a, glm::f64vec3 b, glm::f64vec3 c){return (b.x- a.x)*(c.y-a.y) - (b.y-a.y)*(c.x-a.x);};

    if(Args.size() >= 1)S.push_back(Args[Args.size() - 1].second);
    S.push_back(p_ml);
    glm::f64vec3 top, next;

    for(int i = 0; i < (int)Args.size(); i++){
        do{
            top = S.back();
            S.pop_back();
            next = S.back();
            if(SignedArea(next, Args[i].second,top) <= 0){
                S.push_back(top);
                S.push_back(Args[i].second);
                break;
            }
        }while(1);

    }
    //S.pop_back();
    return S;
}
//xy平面上に乗っていると仮定
double SignedArea(glm::f64vec3 a, glm::f64vec3 b, glm::f64vec3 p){
    glm::f64vec3 v = a-p, v2 = b-p;
    return v.x * v2.y - v.y*v2.x;
}

bool is_point_on_line(glm::f64vec3& p, glm::f64vec3& lp1, glm::f64vec3& lp2){
    double ac = glm::distance(p, lp1), bc = glm::distance(p, lp2), lp = glm::distance(lp1, lp2);
    if(abs(lp - ac - bc) < 1e-5) return true;
    return false;
}

//split at t for bezier curve
#include <iomanip>
std::pair<std::vector<glm::f64vec3>, std::vector<glm::f64vec3>> BezierSplit(std::vector<glm::f64vec3> CtrlPts, double t, int dim){
    std::vector<glm::f64vec3>lp, rp;
    std::vector<std::vector<glm::f64vec3>> Tree = de_casteljau_algorithm(CtrlPts, t);
    Tree.insert(Tree.begin(), CtrlPts);

    for(auto& d: Tree){
        lp.push_back(d[0]);
        rp.push_back(d[d.size() - 1]);
    }
    return {lp, rp};

    {
    using namespace Eigen;
    MatrixXd Ux = create_Ux(dim);
    std::cout<<Ux<<std::endl;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Ux,ComputeFullU | ComputeFullV);
    VectorXd s = svd.singularValues();
    s = s.array().inverse();
    Eigen::MatrixXd Uinv = svd.matrixV() * s.asDiagonal() * svd.matrixU().transpose();
    MatrixXd Z = Eigen::MatrixXd::Identity(dim+1,dim+1);
    for (int i = 0; i < dim + 1; i++)Z(i,i) = Z(i,i) * std::pow(t,i);

    MatrixXd Q = Uinv * Z * Ux;
    MatrixXd X(CtrlPts.size(), 2);
    for(int i = 0; i < (int)CtrlPts.size(); i++){
        X(i,0) = CtrlPts[i].x; X(i,1) = CtrlPts[i].y;
    }

    MatrixXd _Q = MatrixXd::Zero(dim + 1, dim + 1);
    for(int i = dim; i >= 0; i--)_Q.col(i) = Q.col(dim - i);
    for(int i = 0; i < dim+1; i++)_Q.row(i).swap(_Q.row(dim - i));
    MatrixXd left_bezier = Q * X;
    MatrixXd right_bezier = _Q * X;

    auto cast_Eig2GLM = [](VectorXd v){return glm::f64vec3{v(0), v(1), 0};};
    std::vector<glm::f64vec3> bezier_l(dim+1), bezier_r(dim+1);
    for(int i = 0; i < dim + 1; i++){
        bezier_l[i] = cast_Eig2GLM(left_bezier.row(i));
        bezier_r[i] = cast_Eig2GLM(right_bezier.row(i));
    }
    return {bezier_l, bezier_r};
    }
}

Eigen::MatrixXd create_Ux(int dim){
    Eigen::MatrixXd U = Eigen::MatrixXd::Zero(dim, dim);
    auto lmbd = [](int n,int i,int j) { return cmb(n, j) * cmb(n-j, i-j) * std::pow(-1, i - j); };
    for(int i = 0; i < dim; i++){
        for(int j = 0; j < i+1; j++)U(i, j) = lmbd(dim, i, j);
    }
    return U;
}

std::vector<std::vector<glm::f64vec3>> de_casteljau_algorithm(std::vector<glm::f64vec3> CtrlPts, double t){
    std::vector<glm::f64vec3> Q;
    std::vector<std::vector<glm::f64vec3>> _Q;
    glm::f64vec3 prev = CtrlPts[0];
    for(int i = 1; i < (int)CtrlPts.size(); i++){
        glm::f64vec3 new_p = (1. - t) * prev + t * CtrlPts[i];
        Q.push_back(new_p);
        prev = CtrlPts[i];
    }
    if(Q.size() == 0)return _Q;

    _Q.push_back(Q);
    if(Q.size() == 1)return _Q;

    std::vector<std::vector<glm::f64vec3>> tmp = de_casteljau_algorithm(Q,t);
    _Q.insert(_Q.end(), tmp.begin(), tmp.end());
    return _Q;
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

//https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/INT-APP/CURVE-INT-global.html
//http://www.cad.zju.edu.cn/home/zhx/GM/009/00-bsia.pdf
std::vector<glm::f64vec3> GlobalSplineInterpolation(std::vector<glm::f64vec3>& Q, std::vector<glm::f64vec3>& CtrlPts_res){
    using namespace Eigen;
    double L= 0.0;
    int n = Q.size();
    MatrixXd D(n,3);
    for(int i = 0; i < n; i++){
        D(i,0) = Q[i].x; D(i,1) = Q[i].y; D(i,2) = Q[i].z;
    }

    //The Centripetal Method
    for(int i = 1; i < n; i++)L += std::sqrt(glm::distance(Q[i], Q[i - 1]));
    std::vector<double> T(n); T[0] = 0; T[n - 1] = 1;
    for(int i = 1; i < n - 1; i++)T[i] += T[i - 1] + std::sqrt(glm::distance(Q[i], Q[i-1]))/L;

    //average method
    int dim = 3;
    int m = n + dim + 1;
    std::vector<double> U(m, 0);
    for(int i = 0; i <= dim; i++) U[m - i - 1] = 1;
    for(int j = 1; j < n - dim; j++){
        for(int i = j; i < j + dim; i++)U[j + dim] += T[i]/(double)dim;
    }
    //Knot Matrix
    MatrixXd N(n, n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++)N(i,j) = basis(j, dim, T[i], U);
    }

    MatrixXd P = N.fullPivLu().solve(D);

    //cast
    CtrlPts_res.resize(n);
    for(int i = 0; i < n; i++){
        CtrlPts_res[i].x = P(i,0); CtrlPts_res[i].y = P(i,1); CtrlPts_res[i].z = P(i,2);
    }

    int num = 1000;
    double t = T[dim];
    std::vector<glm::f64vec3> BCurve(num);
    for(int i = 0; i < num; i++){
        BCurve[i] = bspline(CtrlPts_res, t, dim, T);
        t += (T[n - dim] - T[dim])/(double)(num);
    }
    return BCurve;
}
