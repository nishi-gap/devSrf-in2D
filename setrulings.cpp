#include "setrulings.h"
using namespace MathTool;

Vertex::Vertex(glm::f64vec3 _p){
    p2_ori = p = _p;
    p3_ori = p3 = _p;
    halfedge.clear();
    deformed = false;
}
Vertex::Vertex(glm::f64vec3 _p2, glm::f64vec3 _p3){
    p2_ori = p = _p2;
    p3_ori = p3 = _p3;
    halfedge.clear();
    deformed = false;
}
void Vertex::addNewEdge(HalfEdge *he){
    if(std::find(halfedge.begin(), halfedge.end(), he) == halfedge.end())halfedge.push_back(he);
}

Vertex::Vertex(const Vertex* v){
    p2_ori = p = v->p; p3_ori = p3 = v->p3;
    halfedge = v->halfedge;
    deformed = v->deformed;
}

Vertex::~Vertex(){ for(auto* h: halfedge){if(h == nullptr)delete h;}}

HalfEdge::HalfEdge(Vertex *v, EdgeType _type){
    vertex = v;
    pair = prev = next =  nullptr;
    IsCrossed = -1;
    angle = 0.0;
    //r = nullptr;
    edgetype = _type;
    v->addNewEdge(this);
}

HalfEdge::HalfEdge(const HalfEdge* he){
    vertex = new Vertex(he->vertex);
    pair = prev = next =  nullptr;
    IsCrossed = he->IsCrossed;
    angle = he->angle;
    edgetype = he->edgetype;
    vertex->addNewEdge(this);
}

HalfEdge::~HalfEdge(){
}
std::vector<HalfEdge*> HalfEdge::Split(Vertex *v, std::vector<HalfEdge*>& Edges){
    std::vector<HalfEdge*> res;
    if(!MathTool::is_point_on_line(v->p, this->vertex->p, this->next->vertex->p))return res;
    double t = glm::length(v->p - vertex->p)/glm::length(next->vertex->p - vertex->p);
    v->p3 = vertex->p3 + t * (next->vertex->p3 - vertex->p3);
    HalfEdge *h_new = new HalfEdge(v, edgetype); h_new->r = r;
    h_new->face = face; h_new->next = next; next->prev = h_new;
    h_new->prev = this; next = h_new;
    res.push_back(h_new);
    Edges.push_back(h_new);
    if(pair != nullptr){
        HalfEdge *h2_new = new HalfEdge(v, edgetype); h2_new->r = r;
        h2_new->face = pair->face;
        h2_new->next = pair->next; h2_new->prev = pair;
        pair->next->prev = h2_new; pair->next = h2_new;
        h_new->pair = pair; h2_new->pair = this;
        pair->pair = h_new; pair = h2_new;
        Edges.push_back(h2_new);
        res.push_back(h2_new);
    }
    return res;
}

bool HalfEdge::hasCrossPoint2d(glm::f64vec3 p, glm::f64vec3 q, glm::f64vec3& CrossPoint, bool ConsiderEnd){
    bool res = MathTool::IsIntersect(p,q, vertex->p, next->vertex->p, ConsiderEnd);
    if(res)CrossPoint = MathTool::getIntersectionPoint(p,q, vertex->p, next->vertex->p);
    return res;
}

double HalfEdge::diffEdgeLength(){
    return abs(glm::length(next->vertex->p - vertex->p) - glm::length(next->vertex->p3 - vertex->p3));
}

HalfEdge* HalfEdge::erase(std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces){
    std::vector<HalfEdge*>::iterator itr = std::find(Edges.begin(), Edges.end(), this);
    if(itr == Edges.end())return nullptr;
    HalfEdge *next = this->next;
    if(this->pair != nullptr){
        pair->prev->next = pair->next; pair->next->prev = pair->prev;
        pair->face->ReConnect(pair->prev);
        std::vector<HalfEdge*>::iterator itr_p = std::find(Edges.begin(), Edges.end(), pair);
        std::vector<Face*>::iterator itr_f = std::find(Faces.begin(), Faces.end(), this->face);
        Faces.erase(itr_f);
        Edges.erase(itr_p);
        delete r;
        delete pair;
    }else{
        prev->next = this->next; this->next->prev = prev;
        face->ReConnect(prev);
    }

    Edges.erase(itr);
    delete this;
    return next;
}


void CrvPt_FL::set(glm::f64vec3 _p, Vertex *o, Vertex *e){
    double sa = glm::distance(_p, o->p), sc = glm::distance(o->p, e->p);
    ve = e; vo = o;
    rt = sa/sc;
    p3 = rt * (e->p3 - o->p3) + o->p3;
    IsValid = true;
    this->p = _p;
}

Vertex4d::Vertex4d(CrvPt_FL *v, Vertex *v2, Vertex *v3){
    first = v; second = v2; third = v3; IsCalc = true;
}
Vertex4d::Vertex4d(const Vertex4d& V4d){
    first = V4d.first; second = V4d.second; third = V4d.third; IsCalc = V4d.IsCalc;
}
Vertex4d::Vertex4d(){first = nullptr; second = nullptr; third = nullptr; IsCalc = false;}

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


Face::Face(HalfEdge *_halfedge){
    //rulings.clear();
    halfedge = _halfedge;
    _halfedge->face = this;
    bend = false;
    hasGradPt = false;
}

Face::Face(const Face& face){
  auto _h = face.halfedge;
  halfedge->edgetype = _h->edgetype;
}

int Face::edgeNum(bool PrintVertex){
    int cnt = 0;
    HalfEdge *h = halfedge;
    do{
        cnt++;
        if(PrintVertex)std::cout << glm::to_string(h->vertex->p) << " , " << h << " , " << h->vertex << std::endl;
        h = h->next;

    }while(h != halfedge);
    if(PrintVertex)std::cout << std::endl;
    return cnt;
}

void Face::ReConnect(HalfEdge *he){    
    halfedge = he;
    HalfEdge *h = he;
   do{
        h->face = this;
        h = h->next;
    }while(he != h);
}

void Face::TrianglationSplit(std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces){

    int n = edgeNum(false);
    HalfEdge *h = halfedge;
    do{
        if(glm::length(h->vertex->p - h->next->vertex->p) < 1e-6)h = h->erase(Edges, Faces);
        else h = h->next;
    }while(h != halfedge);
    while(n > 3){
        HalfEdge *next = halfedge->next, *cur = halfedge, *prev = halfedge->prev;

        n = edgeNum(false);
        do{
            std::array<glm::f64vec3, 3> tri = {prev->vertex->p, cur->vertex->p, next->vertex->p};
            bool elimTriMesh = true;
            HalfEdge *p = next;
            do{
                bool check1 = MathTool::hasPointInTriangle3D(p->vertex->p, tri);
                bool check2 = MathTool::IsAngleLessThan180(tri[1], tri[0], tri[2]);
                if(check1 || !check2){
                    elimTriMesh = false;
                    break;
                }
                p = p->next;
            }while(p != next);
            if(elimTriMesh){
                HalfEdge *h = new HalfEdge(cur->vertex, EdgeType::ol), *h2 = new HalfEdge(next->next->vertex, EdgeType::ol);
                h->pair = h2; h2->pair = h;
                h2->next = cur; h2->prev = next;
                h->next = next->next; h->prev = prev;
                prev->next  = next->next->prev = h;
                cur->prev = next->next = h2;
                Edges.push_back(h); Edges.push_back(h2);
                Face *f = new Face(h); f->ReConnect(h); cur->face->ReConnect(cur);
                Faces.push_back(f);
                cur = h;
                break;
            }
            cur = cur->next;
        }while(cur != halfedge);
        n = edgeNum(false);
    }
}

CRV::CRV(int _crvNum, int DivSize){

    curveNum = _crvNum;
    ControllPoints.clear();
    CurvePoints.clear();
    for(int i = 0; i < curveNum; i++)CurvePoints.push_back(crvpt(i));


    for(int i = 0; i < DivSize-1;i++){
        Rulings.push_back(new ruling());
    }

    isempty = true;
    curveType = CurveType::none;
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


//void CRV::addCtrlPt(QPointF p){ControllPoints.push_back(glm::f64vec3{p.x(), p.y(), 0}); qDebug()<<"addCtrlPt called and size is " << ControllPoints.size(); }


void CRV::eraseCtrlPt(int curveDimention, int crvPtNum){
    if((int)ControllPoints.size() > 0)ControllPoints.erase(ControllPoints.end() - 1);
    ControllPoints.shrink_to_fit();
    if(curveType == CurveType::bezier3){
        Bezier(curveDimention, crvPtNum);
    }else if(curveType == CurveType::bsp3){
        Bspline(curveDimention, crvPtNum);
    }else if(curveType == CurveType::line){
        Line();
    }
}

//void CRV::movePt(glm::f64vec3 p, int ind){ControllPoints[ind] = p;}

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

CurveType CRV::getCurveType(){return curveType;}
void CRV::setCurveType(CurveType n){curveType = n;}

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

void CRV::Bezier(int curveDimention, int crvPtNum){

    if((int)ControllPoints.size() < curveDimention){return;}

    double t = 0.0;
    glm::f64vec3 v;
    for(int n = 0; n < crvPtNum; n++){
        v = {0.f,0.f, 0.f};
        for (int i = 0; i < int(ControllPoints.size()); i++)v += MathTool::BernsteinBasisFunc(curveDimention, i, t) * ControllPoints[i];

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
        if (IsIntersect(v, v2, N0, N1, true)) {
            auto tmp = getIntersectionPoint(v, v2, N0, N1);
            bool hasSamePoint = false;
            for(auto&c: crossPoint){
                if(glm::distance(c, tmp) < 1e-5)hasSamePoint = true;
            }
            if(!hasSamePoint)crossPoint.push_back(tmp);
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
    int i = 0, rind = 0;
    std::vector<Vertex*> vertices = outline->getVertices();

    while(i < m){
        T = 3.0 * (-ControllPoints[0] + 3.0 * ControllPoints[1] - 3.0 * ControllPoints[2] + ControllPoints[3]) * t * t
                + 6.0 * (ControllPoints[0] - 2.0 * ControllPoints[1] + ControllPoints[2]) * t + 3.0 * (-ControllPoints[0] + ControllPoints[1]);//一階微分(tについて)
        N = l * glm::f64vec3{ -T.y, T.x, 0 };
        //if (glm::dot(glm::normalize(N), glm::f64vec3{0,-1, 0}) < 0) N *= -1;
        int n = std::min(int(((double)i/(double)(m - 1)) * crvPtNum), crvPtNum - 1);
        bool IsIntersected = setPoint(vertices, N, CurvePoints[n].pt, crossPoint);

        if(IsIntersected){
            for(int j = 0; j < (int)crossPoint.size(); j += 2){
                Rulings[rind]->r = std::make_tuple(new Vertex(crossPoint[j]), new Vertex(crossPoint[j+1]));
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
    double n = double(sind);
    int rind = 0;
    //ruling *r;
    while(n < eind){
        Rulings[rind]->r = std::make_tuple(new Vertex(CrossPoints[floor(n)][0]), new Vertex(CrossPoints[floor(n)][1]));
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
    //bool RulingSet = false;
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
            if(ind >=CurveIndexs[i] && ind < CurveIndexs[i + 1]){
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

void CRV::FillColor(int c){
    for(auto& cp: CurvePoints) cp.color = c;
}

OUTLINE::OUTLINE(){
    type = "Rectangle";
    vertices.clear();
    //edges.clear();
    Lines.clear();
    //face = nullptr;
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
            vertices.insert(vertices.begin() + 1, new Vertex(glm::f64vec3{p.x, vertices[0]->p.y, 0}));
            vertices.push_back(new Vertex(glm::f64vec3{vertices[0]->p.x, p.y, 0}));
            ConnectEdges();
            if(glm::dot(getNormalVec(), glm::f64vec3{0,0,1}) > 0){
                //vertices[1]->p = vertices[1]->p3 = glm::f64vec3{p.x, vertices[0]->p.y, 0}; vertices[3]->p = vertices[3]->p3 = glm::f64vec3{vertices[0]->p.x, p.y, 0};
            }
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
        double a = 2 * std::numbers::pi /VerticesNum;

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
        ind = movePointIndex(p);
        if(ind != -1){
            T(0,2) = p.x - vertices[ind]->p.x; T(1,2) = p.y - vertices[ind]->p.y;
            for(auto& v: vertices){
                x = T * Eigen::Vector3d(v->p.x, v->p.y, 1);
                 v->p2_ori.x = v->p.x = x(0);  v->p2_ori.y = v->p.y = x(1);
            }
            x = T * Eigen::Vector3d(origin.x, origin.y, 1);
            origin.x = x(0); origin.y = x(1);
        }else if(d < dist){
            T(0,2) = p.x - origin.x; T(1,2) = p.y - origin.y;
            for(auto& v: vertices){
                x = T * Eigen::Vector3d(v->p.x, v->p.y, 1);
                v->p2_ori.x = v->p.x = x(0);  v->p2_ori.y = v->p.y = x(1);
            }
            x = T * Eigen::Vector3d(origin.x, origin.y, 1);
            origin.x = x(0); origin.y = x(1);
        }
    }else {
        ind = movePointIndex(p);
        if(ind == -1)return;
        T(0,2) = p.x - vertices[ind]->p.x; T(1,2) = p.y - vertices[ind]->p.y;
        for(auto& v: vertices){
            x = T * Eigen::Vector3d(v->p.x, v->p.y, 1);
            v->p2_ori.x = v->p.x = x(0);  v->p2_ori.y = v->p.y = x(1);
        }
    }
}

void OUTLINE::MoveVertex(glm::f64vec3 p, int ind){
    if(ind < 0 || (int)vertices.size() < ind)return;
    vertices[ind]->p2_ori = vertices[ind]->p = p;
}

std::vector<Vertex*> OUTLINE::getVertices(){return vertices;}
//Face* OUTLINE::getFace(){return face;}

void OUTLINE::ConnectEdges(bool IsConnected){
    //if(!isClosed)return;
    if(vertices.size() < 2)return;
    Lines.clear();
    if(IsConnected){
        for(int i = 0; i < (int)vertices.size(); i++)Lines.push_back(new Line(vertices[i], vertices[(i+1) % (int)vertices.size()], EdgeType::ol));
        IsConnected = true;
    }else{

    }
    /*
    for(auto&v: vertices){
        HalfEdge *he = new HalfEdge(v, EdgeType::ol);
        edges.push_back(he);
    }
    if(IsConnected){
        face = new Face(edges[0]);
        for(int i = 0; i < (int)edges.size(); i++){
            edges[i]->prev = edges[(i + 1) % edges.size()];
            edges[i]->next = edges[(i - 1) % edges.size()];
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
    }*/

}

bool OUTLINE::IsPointInFace(glm::f64vec3 p){
    if(!IsClosed())return false;
    int cnt = 0;
    double vt;
    for(auto&l: Lines){
        if(l->o->p.y <= p.y && l->v->p.y > p.y){
            vt = (p.y - l->o->p.y) / (l->v->p.y - l->o->p.y);
            if(p.x < (l->o->p.x + vt * (l->v->p.x - l->o->p.x)))cnt++;
        }else if(l->o->p.y > p.y && l->v->p.y <= p.y){
            vt = (p.y - l->o->p.y)/(l->v->p.y - l->o->p.y);
            if(p.x < (l->o->p.x + vt * (l->v->p.x - l->o->p.x)))cnt--;
        }
    }
    if (cnt == 0) return false;
    else return true;
}


bool OUTLINE::IsClosed(){
    if(Lines.empty())return false;
    /*
    HalfEdge *h = edges[0];
    do{
        if(h->next == nullptr)return false;
        h = h->next;
    }while(h != edges[0]);
    */
    for(auto& l: Lines)if(l == nullptr)return false;
    return true;
}


void CrossDetection(OUTLINE *outline, CRV *crvs){
    if(crvs->isempty)return;
    for(auto&r: crvs->Rulings)r->IsCrossed = -1;
    if(crvs->getCurveType() == CurveType::arc || crvs->getCurveType() == CurveType::line)return;
    for(int in = 0; in < (int)crvs->Rulings.size(); in ++){
        for(int inn = in+1; inn < (int)crvs->Rulings.size(); inn++){
            bool rs = IsIntersect(std::get<0>(crvs->Rulings[in]->r)->p, std::get<1>(crvs->Rulings[in]->r)->p, std::get<0>(crvs->Rulings[inn]->r)->p,std::get<1>(crvs->Rulings[inn]->r)->p, false);
            if(rs){
                glm::f64vec3 p = getIntersectionPoint(std::get<0>(crvs->Rulings[in]->r)->p, std::get<1>(crvs->Rulings[in]->r)->p, std::get<0>(crvs->Rulings[inn]->r)->p,std::get<1>(crvs->Rulings[inn]->r)->p);
                bool PointOnLines = false;
                bool PointInFace = outline->IsPointInFace(p);
                for(auto&l: outline->Lines){if(is_point_on_line(p, l->o->p, l->v->p))PointOnLines = true;}
                if(PointInFace)crvs->Rulings[in]->IsCrossed = crvs->Rulings[inn]->IsCrossed = 0;
            }
        }
    }
}

glm::f64vec3 OUTLINE::getNormalVec(){
    glm::f64vec3 N{0,0,0};
    auto prev = Lines.end() - 1;
    for(auto cur = Lines.begin(); cur != Lines.end(); cur++){
        if(abs(glm::dot(glm::normalize((*cur)->v->p - (*cur)->o->p), glm::normalize((*prev)->v->p - (*prev)->o->p))) < 0.9){
            return glm::normalize(glm::cross((*cur)->v->p - (*cur)->o->p, (*prev)->o->p - (*prev)->v->p));
        }
        prev = cur;

    }
    return N;
}


int OUTLINE::movePointIndex(glm::f64vec3 p){
    int n = vertices.size();
    double dist = 5.0;
    int ind = -1;
    for(int i = 0; i < n; i++){
        double d = glm::distance(p,vertices[i]->p);
        if(d < dist){
            dist = d;
            ind = i;
        }
    }
    return ind;
}

std::vector<glm::f64vec3> ConvertDistBasedBezier(std::vector<glm::f64vec3>& CtrlPts, Line *l){
    glm::f64vec3 p, q;
    if(l->o->p.x <= l->v->p.x){p = l->v->p; q = l->v->p;}
    else{q = l->v->p; p = l->v->p;}
    double a = p.y - q.y, b = q.x - p.x, c = p.x * q.y - q.x * p.y;
    int n = CtrlPts.size();
    std::vector<glm::f64vec3> D(n);
    for(int i = 0; i < n; i++){
        double d = -(a * CtrlPts[i].x + b * CtrlPts[i].y + c)/sqrt(a*a + b*b);
        D[i] = glm::f64vec3{(double)i/(double)(n - 1), d, 0};       
    }
    return D;
}

std::vector<double> BezierClipping(std::vector<glm::f64vec3>&CtrlPts, Line *l, int dim){
    std::vector<glm::f64vec3> base = ConvertDistBasedBezier(CtrlPts, l);
    std::vector<glm::f64vec3> current;
    std::copy(base.begin(), base.end(), std::back_inserter(current));
    std::array<glm::f64vec3, 2> _line{glm::f64vec3{0.,0.0,0.0}, glm::f64vec3{1.,0.0,0.0}};
    auto res = _bezierclipping(base, current, _line, dim);

    return res;
}

std::vector<Vertex*> ConvexHull_polygon(std::vector<Vertex*>& Q){
    std::vector<Vertex*> P;
    if((int)Q.size() < 3)return Q;
    Vertex* p_ml = Q[0];
    for(auto&p: Q){
        if(p_ml->p.y > p->p.y)p_ml = p;
        else if(p_ml->p.y == p->p.y && p_ml->p.x > p->p.x) p_ml = p;
    }
    std::vector<std::pair<double, Vertex*>>Args;
    for(auto&p: Q){
        if(p->p == p_ml->p)continue;
        double phi = atan2(p->p.y - p_ml->p.y, p->p.x - p_ml->p.x);
        bool hasSameAngle = false;
        for(auto& X: Args){
            if(X.first == phi){
                if(glm::distance(X.second->p, p_ml->p) < glm::distance(p->p, p_ml->p))X.second = p;
                hasSameAngle = true;
            }
        }
        if(!hasSameAngle)Args.push_back(std::make_pair(phi, p));
    }
    // compare only the first value
    std::sort(Args.begin(), Args.end(),[](auto const& x, auto const& y) {return x.first < y.first; });
    if(Args.size() >= 1)P.push_back(Args[Args.size() - 1].second);
    P.push_back(p_ml);
    Vertex *top, *next;
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
    return P;
}
std::vector<Vertex*> SortPolygon(std::vector<Vertex*>& polygon){
    glm::f64vec3 Zaxis{0,0,1};//面の表裏の基準
    auto poly_sort = ConvexHull_polygon(polygon);
    if((int)poly_sort.size()< 3){std::cout<<"not enought vertices size"<<std::endl; return{};}
    if(glm::dot(glm::cross(poly_sort[0]->p - poly_sort[1]->p, poly_sort[2]->p - poly_sort[1]->p), Zaxis) < 0)
        std::reverse(poly_sort.begin(), poly_sort.end());
    return poly_sort;
}

