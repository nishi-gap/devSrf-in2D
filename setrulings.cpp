#include "setrulings.h"
using namespace MathTool;

Vertex::Vertex(glm::f64vec3 _p){
    p = _p;
    p3 = _p;
    halfedge.clear();
    deformed = false;
}
Vertex::Vertex(glm::f64vec3 _p2, glm::f64vec3 _p3){
    p = _p2;
    p3 = _p3;
    halfedge.clear();
    deformed = false;
}
void Vertex::addNewEdge(HalfEdge *he){
    halfedge.push_back(he);
}

HalfEdge::HalfEdge(Vertex *v, EdgeType _type){
    vertex = v;
    pair = nullptr;
    prev = nullptr;
    next = nullptr;
    r = nullptr;
    edgetype = _type;
    v->addNewEdge(this);
}

std::vector<HalfEdge*> HalfEdge::Split(Vertex *v, std::vector<HalfEdge*>& Edges){
    std::vector<HalfEdge*> res;
    if(!MathTool::is_point_on_line(v->p, this->vertex->p, this->next->vertex->p))return res;
    double t = glm::length(v->p - vertex->p)/glm::length(next->vertex->p - vertex->p);
    v->p3 = vertex->p3 + t * (next->vertex->p3 - vertex->p3);
    HalfEdge *h_new = new HalfEdge(v, edgetype); h_new->r = r;
    h_new->face = face; h_new->next = next; h_new->prev = this; next = h_new;
    HalfEdge *h2_new = new HalfEdge(v, edgetype); h2_new->r = r;
    res.push_back(h_new);
    Edges.push_back(h_new);
    if(pair != nullptr){

        h2_new->face = pair->face;
        h2_new->next = pair->next; h2_new->prev = pair; pair->next = h2_new;
        h2_new->pair = this; h_new->pair = pair; pair->pair = h_new; pair = h2_new;
        Edges.push_back(h2_new);
        res.push_back(h2_new);
    }
    return res;
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


Face::Face(HalfEdge *_halfedge){
    //rulings.clear();
    halfedge = _halfedge;
    bend = false;
    hasGradPt = false;
}

bool Face::IsPointInFace(glm::f64vec3 p){
    int cnt = 0;
    double vt;
    HalfEdge *he = halfedge;
    do{
        if(he->vertex->p.y <= p.y && he->next->vertex->p.y > p.y){
            vt = (p.y - he->vertex->p.y) / (he->next->vertex->p.y - he->vertex->p.y);
            if(p.x < (he->vertex->p.x + vt * (he->next->vertex->p.x - he->vertex->p.x)))cnt++;
        }else if(he->vertex->p.y > p.y && he->next->vertex->p.y <= p.y){
            vt = (p.y - he->vertex->p.y)/(he->next->vertex->p.y - he->vertex->p.y);
            if(p.x < (he->vertex->p.x + vt * (he->next->vertex->p.x - he->vertex->p.x)))cnt--;
        }
        he = he->next;
    }while(he != halfedge);
    if (cnt == 0) return false;
    else return true;
}

glm::f64vec3 Face::getNormalVec(){
    glm::f64vec3 N{0,0,0};
    HalfEdge *h = this->halfedge;
    do{
        if(!is_point_on_line(h->vertex->p3, h->prev->vertex->p3, h->next->vertex->p3))
            return glm::normalize(glm::cross(h->prev->vertex->p3 - h->vertex->p3, h->next->vertex->p3 - h->vertex->p3));
        h = h->next;
    }while(h != halfedge);
    return N;
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


void CRV::addCtrlPt(QPointF p){ControllPoints.push_back(glm::f64vec3{p.x(), p.y(), 0}); qDebug()<<"addCtrlPt called and size is " << ControllPoints.size(); }

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

CurveType CRV::getCurveType(){return curveType;}
void CRV::setCurveType(CurveType n){curveType = n;}

void CRV::Bspline(int curveDimention,  int crvPtNum){
    //Rulings.clear();
    ClearPt();
    if((int)ControllPoints.size() <= curveDimention) return;

    int knotSize = (int)ControllPoints.size() + curveDimention + 1;
    std::vector<double>T(knotSize);
    for(int j = 0; j < knotSize; j++)T[j] = (double)j/(double)knotSize;

    //????????????????????????
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
    //P.push_back(crossPoint[2 * minInd]); P.push_back(crossPoint[2 * minInd + 1]);//?????????????
    for(int i = 0; i < (int)crossPoint.size(); i++){
        //if(i != minInd || i != minInd + 1)
        P.push_back(crossPoint[i]);
    }
    return IsIntersected;
}

void CRV::BezierRulings(OUTLINE *outline, int& DivSize, int crvPtNum){
    std::cout << "you can't use now"<<std::endl;
    return;
    double l = 1000; //?????????????????????
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
                + 6.0 * (ControllPoints[0] - 2.0 * ControllPoints[1] + ControllPoints[2]) * t + 3.0 * (-ControllPoints[0] + ControllPoints[1]);//????????????(t????????????)
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
    double l = 1000;//?????????????????????
    glm::f64vec3 N, T;
    std::vector<glm::f64vec3> crossPoint;
    std::vector<std::vector<glm::f64vec3>> CrossPoints;
    double t = 0.0;
    std::vector<Vertex*> vertices = outline->getVertices();
    int knotSize = (int)ControllPoints.size() + curveDimention + 1;
    std::vector<double>Knot(knotSize);
    for(int j = 0; j < knotSize; j++)Knot[j] = (double)j/(double)knotSize;
    //????????????????????????
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

//?????????????????????2?????????
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
    double l = 1000;//?????????????????????

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

//????????? 0: ??????. 1,2 ??????
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
    double l = 1000;//?????????????????????
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
    else if(d == 0){//?????????
        for(int i = 0; i < (int)CurveIndexs.size() - 1; i++){
            if(ind >=CurveIndexs[i] && ind < CurveIndexs[i + 1]){
                double s = (double)(ind - CurveIndexs[i])/(double)(CurveIndexs[i + 1] - CurveIndexs[i]);
                InsertPoint = s * (ControllPoints[i + 1] - ControllPoints[i]) + ControllPoints[i];
                InsertPointSegment = i;
                return;
            }
        }
    }else if(d == 1){//?????????
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
        ind = movePointIndex(p);
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
        ind = movePointIndex(p);
        if(ind == -1)return;
        T(0,2) = p.x - vertices[ind]->p.x; T(1,2) = p.y - vertices[ind]->p.y;
        for(auto& v: vertices){
            x = T * Eigen::Vector3d(v->p.x, v->p.y, 1);
            v->p.x = x(0); v->p.y = x(1);
        }
    }
}

void OUTLINE::MoveVertex(glm::f64vec3 p, int ind){
    if(ind < 0 || (int)vertices.size() < ind)return;
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

void CrossDetection(Face *f, CRV *crvs){
    if(crvs->isempty)return;
    for(auto&r: crvs->Rulings)r->IsCrossed = -1;
    if(crvs->getCurveType() == CurveType::arc || crvs->getCurveType() == CurveType::line)return;
    for(int in = 0; in < (int)crvs->Rulings.size(); in ++){
        for(int inn = in+1; inn < (int)crvs->Rulings.size(); inn++){
            bool rs = IsIntersect(std::get<0>(crvs->Rulings[in]->r)->p, std::get<1>(crvs->Rulings[in]->r)->p, std::get<0>(crvs->Rulings[inn]->r)->p,std::get<1>(crvs->Rulings[inn]->r)->p);
            if(rs){
                glm::f64vec3 p = getIntersectionPoint(std::get<0>(crvs->Rulings[in]->r)->p, std::get<1>(crvs->Rulings[in]->r)->p, std::get<0>(crvs->Rulings[inn]->r)->p,std::get<1>(crvs->Rulings[inn]->r)->p);
                bool PointOnLines = false;
                bool PointInFace = f->IsPointInFace(p);
                HalfEdge *h = f->halfedge;
                do{
                    if(is_point_on_line(p, h->vertex->p, h->next->vertex->p))PointOnLines = true;
                    h = h->next;
                }while(h != f->halfedge);
                if(PointInFace)crvs->Rulings[in]->IsCrossed = crvs->Rulings[inn]->IsCrossed = 0;
            }
        }
    }
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

std::vector<double> BezierClipping(std::vector<glm::f64vec3>&CtrlPts, HalfEdge *line, int dim){
    std::vector<glm::f64vec3> base = ConvertDistBasedBezier(CtrlPts, line);
    std::vector<glm::f64vec3> current;
    std::copy(base.begin(), base.end(), std::back_inserter(current));
    std::array<glm::f64vec3, 2> _line{glm::f64vec3{0.,0,0}, glm::f64vec3{1.,0,0}};
    auto res = _bezierclipping(base, current, _line, dim);
    return res;
}


//https://pages.mtu.edu/~shene/COURSES/cs3621/NOTES/INT-APP/CURVE-INT-global.html
//http://www.cad.zju.edu.cn/home/zhx/GM/009/00-bsia.pdf
std::vector<glm::f64vec3> GlobalSplineInterpolation(std::vector<CrvPt_FL>& Q, std::vector<glm::f64vec3>& CtrlPts_res, std::vector<double>& Knot, bool is3d){
    using namespace Eigen;
    double L= 0.0;
    int n = Q.size() - 1;
    if(n < 0)return std::vector<glm::f64vec3>{};
    MatrixXd D(n + 1,3);
    for(int i = 0; i <= n; i++){
        if(is3d){D(i,0) = Q[i].p3.x; D(i,1) = Q[i].p3.y; D(i,2) = Q[i].p3.z; }
        else{D(i,0) = Q[i].p.x; D(i,1) = Q[i].p.y; D(i,2) = Q[i].p.z; }
    }

    //The Centripetal Method
    std::vector<double> T(n+1, 0); T[n] = 1;
    if(is3d){
        for(int i = 1; i <= n; i++)L += sqrt(glm::distance(Q[i].p3, Q[i - 1].p3));
        for(int i = 1; i < n; i++)T[i] = T[i - 1] + sqrt(glm::distance(Q[i].p3, Q[i-1].p3))/L;
    }else{
        for(int i = 1; i <= n; i++)L += sqrt(glm::distance(Q[i].p, Q[i - 1].p));
        for(int i = 1; i < n; i++)T[i] = T[i - 1] + sqrt(glm::distance(Q[i].p, Q[i-1].p))/L;
    }

    //average method
    int dim = 3;
    int m = n + dim + 1;
    Knot.assign(m+1,0);
    for(int j = 1; j <= n - dim; j++){
        for(int i = j; i < j + dim; i++)Knot[j + dim] += T[i]/(double)dim;
    }
    for(int i = 0; i <= dim; i++) Knot[m - i] = 1;
    //Knot Matrix
    MatrixXd N = MatrixXd::Zero(n+1, n+1);
    for(int i = 0; i <= n; i++){
        if(T[i] == Knot[0]){ N(i,0) = 1;continue;}
        if(T[i] == Knot[m]){N(i,n) = 1; continue;}
        for(int j = 0; j <= n; j++)N(i,j) = basis(j,dim,T[i],Knot);
    }
    MatrixXd P = N.fullPivLu().solve(D);
    //cast
    CtrlPts_res.resize(n+1);
    for(int i = 0; i <= n; i++){
        CtrlPts_res[i].x = P(i,0); CtrlPts_res[i].y = P(i,1); CtrlPts_res[i].z = P(i,2);
    }

    int num = 5000;
    double t = Knot[dim];
    std::vector<glm::f64vec3> BCurve(num);
    for(int i = 0; i < num; i++){
        BCurve[i] = bspline(CtrlPts_res, t, dim, Knot);
        t += (Knot[CtrlPts_res.size()] - Knot[dim])/(double)(num);
    }
    return BCurve;
}
