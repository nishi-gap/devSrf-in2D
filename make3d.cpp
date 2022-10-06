#define PI 3.14159265359
#include "make3d.h"

FaceGradation::FaceGradation(){
    color = 0;
    he = nullptr;
}

FaceGradation::FaceGradation(HalfEdge *_he, double *_color){
    he = _he;
    color = _color;
}

Model::Model(){
    clear();
    outline = new OUTLINE();
}

Model::Model(int _crvPtNum){
    crvPtNum = _crvPtNum;
    clear();
    outline = new OUTLINE();
    refCrv.clear();
    Fgrad.clear();
    befFaceNum = 0;
}

void Model::clear(){
    vertices.clear();
    Edges.clear();
    Faces.clear();
}

void Model::Initialize(){
    clear();
    outline = new OUTLINE();
    refCrv.clear();
}

LinearRange::LinearRange(){
    face = nullptr;
    CrvInd = -1;
    RulingOnCurve = -1;
    RulInd = -1;
}

void Model::InsertVertex(Vertex *v){
    for(auto& f: Faces){
        HalfEdge *he = f->halfedge;
        do{
            if ((glm::dot(glm::normalize(he->vertex->p - v->p), glm::normalize(he->next->vertex->p - v->p))) <= -1 + FLT_EPSILON){
                HalfEdge *NewEdge = new HalfEdge(v);
                HalfEdge *next = he->next;
                he->next = NewEdge;
                next->prev = NewEdge;
                NewEdge->prev = he;
                NewEdge->next = next;
                NewEdge->face = f;
                Edges.push_back(NewEdge);
                return;
            }
            he = he->next;
        }while(he != f->halfedge);
    }

    return;
}

void Model::ReInsertVertex(HalfEdge *he, std::vector<HalfEdge*>& edges){
    for(auto& f: Faces){
        HalfEdge *h = f->halfedge;
        do{
            if ((glm::dot(glm::normalize(h->vertex->p - he->vertex->p), glm::normalize(h->next->vertex->p - he->vertex->p))) <= -1 + FLT_EPSILON){
                HalfEdge *NewEdge = new HalfEdge(he->vertex);
                HalfEdge *next = h->next;
                h->next = NewEdge;
                next->prev = NewEdge;
                NewEdge->prev = h;
                NewEdge->next = next;
                NewEdge->face = f;
                edges.push_back(NewEdge);
                return;
            }
            h = h->next;
        }while(h != f->halfedge);
    }
    std::cout <<"not found" << std::endl;
    return;
}

bool Model::devide(HalfEdge *he1, HalfEdge *he2, std::vector<Face*> &faces){
    HalfEdge *h1;
    HalfEdge *h2;

    for(auto& f: faces){
        h1 = nullptr; h2 = nullptr;
        HalfEdge *he = f->halfedge;
        int cnt = 0;
        do{
            if(he->vertex->p == he1->vertex->p){ h1 = he;}
            if(he->vertex->p == he2->vertex->p){ h2 = he;}
            he = he->next;
            cnt++;
        }while(he != f->halfedge);
        if(h1 != nullptr && h2 != nullptr){
            he2->pair = he1;
            he1->pair = he2;

            HalfEdge *h1_prev = h1->prev;
            HalfEdge *h2_prev = h2->prev;
            he1->prev = h1_prev;
            he1->next = h2;
            he2->prev = h2_prev;
            he2->next = h1;
            h1->prev = he2;
            h2->prev = he1;
            h1_prev->next = he1;
            h2_prev->next = he2;

            h1 = he1;
            h2 = he2;
            f->halfedge = h1;
            do{
                h1->face = f;
                h1 = h1->next;
            }while(h1 != he1);
            Face* face2 = new Face(he2);
            faces.push_back(face2);
            face2->halfedge = he2;
            do{
                h2->face = face2;
                h2 = h2->next;
            }while(h2 != he2);
            //Edges.push_back(NewEdge1); Edges.push_back(NewEdge2);
            return true;
        }
    }
    return false;
}

bool Model::setRuling(std::vector<ruling*>& Rulings){
    if(Rulings.empty()) return false;
    for(auto& r: Rulings)if(r->IsCrossed != -1)return false;
    Vertex *v1;
    Vertex *v2;
    for(auto* r: Rulings){
        v1 = std::get<0>(r->r);
        v2 = std::get<1>(r->r);

        InsertVertex(v1);
        InsertVertex(v2);
        HalfEdge *he1 = new HalfEdge(v1, (double)(r->Gradation) * M_PI / 255.f);
        HalfEdge *he2 = new HalfEdge(v2, (double)(r->Gradation) * M_PI / 255.f);
        Edges.push_back(he1); Edges.push_back(he2);
        devide(he1, he2, Faces);
        vertices.push_back(v1); vertices.push_back(v2);
    }

    return true;
}

void Model::setHalfEdgePair(HalfEdge*he){
    for(auto& f: Faces){
        HalfEdge*he_if = f->halfedge;
        do{
            if(he->vertex == he_if->next->vertex && he->next->vertex == he_if->vertex){
                he->pair = he_if; he_if->pair = he;
                return;
            }
            he_if = he_if->next;
        }while(he_if != f->halfedge);
    }
}

void Model::deform(std::vector<std::vector<glm::f64vec3>>& output, std::vector<ruling*>& Rulings, glm::f64vec3& center){

    glm::f64mat4x4 T, R, A;

    std::queue<glm::f64mat4x4> MatQue;
    std::queue<HalfEdge*> FacesQue;
    std::vector<glm::f64vec3> mesh;
    HalfEdge*Pos;
    glm::f64vec3 v1, v2;
    glm::f64mat4x4 Mirror = glm::mat4(1.0f); Mirror[1][1] = -1;
    glm::f64mat4x4 Scale = glm::scale(glm::f64vec3{0.1f, 0.1f, 0.1f});

    HalfEdge *he;

    bool res = setRuling(Rulings);
    if(!res){clear(); return;}
    output.clear();
    for(auto& h: Edges)h->vertex->deformed = false;
    he = (*Faces.begin())->halfedge;
    do{
        T = glm::translate(he->vertex->p  - (*Faces.begin())->halfedge->vertex->p);
        mesh.push_back(T * glm::f64vec4{0,0,0, 1});
        if(he->pair != nullptr){
            FacesQue.push(he->pair);
            MatQue.push(T);
        }
        he = he->next;
        he->vertex->deformed = true;
    }while(he != (*Faces.begin())->halfedge);
    he->face->bend = true;
    output.push_back(mesh);
    while(!FacesQue.empty()){
        
        Pos = FacesQue.front(); A = MatQue.front();
        FacesQue.pop(); MatQue.pop();
        if(Pos->face->bend)continue;
        he = Pos;
        mesh.clear();
        v1 = Pos->pair->vertex->p;
        v2 = Pos->vertex->p;
        glm::f64vec3 axis = glm::normalize(v2 - v1);
        //if(glm::dot(axis, glm::f64vec3{0,1,0}) < 0)axis *= -1;
        R = glm::rotate(he->angle, axis);
        do{
            T = glm::translate(he->vertex->p - v1);
            mesh.push_back(A * R * T * glm::f64vec4{0,0,0, 1});
            he->vertex->deformed = true;
            if(he->pair != nullptr && !he->pair->face->bend){
                FacesQue.push(he->pair);
                MatQue.push(A * R * T);
            }
            he = he->next;
        }while(he != Pos);
        Pos->face->bend = true;
        output.push_back(mesh);
    }
    std::vector<glm::f64vec3>c;
    for(auto& m:output){
        center = {0,0,0};
        for(auto&v:m){
            v =   Scale * Mirror *  glm::f64vec4(v,1);
            center += v;
        }
        center /= double(m.size());
        c.push_back(center);
    }
    center = {0,0,0};
    for(auto& v: c)center += v;
    center /= double(c.size());
    return;
}

void Model::setOutline(std::vector<Vertex*> outline){

    int n = outline.size();
    if(n < 3)return;
    clear();
    HalfEdge*he;
    for(auto& p: outline){
        vertices.push_back(p);
        he = new HalfEdge(p);
        Edges.push_back(he);
    }
    for(int i = 0; i < n; i++){
        Edges[i]->prev = Edges[(i + 1) % n];
        Edges[(i + 1) % n]->next = Edges[i];
    }

    Face *face = new Face(Edges[0]);
    Faces.push_back(face);
    for(auto& he: Edges)he->face = face;
}

void Model::LinearInterPolation(std::list<LinearRange>& path){
    if(path.empty())return;

    glm::f64vec3 befcenter = glm::f64vec3{-1,-1,-1}, center;
    std::vector<glm::f64vec3> vertices;
    double len = 0.0;
    for(auto& l: path){
        if(befcenter == glm::f64vec3{-1,-1,-1}){
            vertices = TranslateGLMfromHE(l.face);
            befcenter = GetCenter(vertices);
            continue;
        }
        vertices = TranslateGLMfromHE(l.face);
        center = GetCenter(vertices);
        len += glm::distance(center, befcenter);
        befcenter = center;
    }
    double r = (linearrange[1].face->Gradation - linearrange[0].face->Gradation)/len;
    Face *bef = nullptr;
    for(auto&l : path){
        if(bef == nullptr){
            bef = l.face;
            vertices = TranslateGLMfromHE(l.face);
            befcenter = GetCenter(vertices);
            continue;
        }
        vertices = TranslateGLMfromHE(l.face);
        center = GetCenter(vertices);
        l.face->Gradation = r * glm::distance(center, befcenter) + bef->Gradation;

        bef = l.face;
        befcenter = center;
    }
}

void Model::SplineInterPolation(std::vector<glm::f64vec2>& CurvePath){//凹型に関してはとりあえず虫
    std::cout <<"can't use now" <<std::endl; return;
    std::vector<crvpt*> Curve;
    std::vector<double>color;
    for(int i = 0; i < (int)Faces.size(); i++){
        if(Faces[i]->hasGradPt){
            if(Curve.size() == 0){
                //Curve.push_back(Faces[i]->rhe->pt);
                color.push_back(Faces[i]->Gradation);
            }else{
                bool hasPoint = false;
                for(auto& p: Curve){
                    for(auto& r: Faces[i]->rulings)
                        if(p->pt == r->pt->pt){
                            hasPoint = true;
                        }
                }
                if(!hasPoint){
                    //Curve.push_back(Faces[i]->rhe->pt);
                    color.push_back(Faces[i]->Gradation);
                }
            }

        }
    }
    if(Curve.size() < 2)return;

    for(int i = 0; i < (int)Curve.size(); i++){
        for(int j = i + 1; j < (int)Curve.size(); j++){
            if(Curve[j]->pt.x < Curve[i]->pt.x){
                std::swap(Curve[i], Curve[j]);
                std::swap(color[i], color[j]);
            }
        }
    }

    CurvePath.clear();
    int N =(int)Curve.size() - 1;

    Eigen::VectorXd v(N - 1);
    std::vector<double>h(N);
    for(int i = 1; i < (int)Curve.size(); i++)h[i - 1] = Curve[i]->pt.x - Curve[i - 1]->pt.x;
    for(int i = 1; i < N; i++){
        double a = (h[i] != 0) ? (double(color[i + 1] - color[i]))/h[i] : 0, b = (h[i - 1] != 0) ? (double(color[i] - color[i - 1]))/h[i - 1]: 0;
        v(i - 1) = 6 * (a - b);
    }
    Eigen::VectorXd u;
    if(N == 1){
        u = Eigen::VectorXd::Zero(2);
    }
    if(N == 2){
        u = Eigen::VectorXd::Zero(3);
        double dx1, dx2, dx3;
        dx1 = Curve[2]->pt.x - Curve[0]->pt.x; dx2 = Curve[2]->pt.x - Curve[1]->pt.x; dx3 = Curve[1]->pt.x - Curve[0]->pt.x;
        if(abs(dx1) < FLT_EPSILON) u(1) = 0;
        else{
            double a = (dx2 == 0) ? 0: (color[2] - color[1])/dx2;
            double b = (dx3 == 0) ? 0: (color[1] - color[0])/dx3;
            u(1) = 3 * (a - b)/dx1;
        }
    }
    if(N > 2){
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N-1,N-1);
        for(int i = 0; i < N-1;i++){
            A(i,i) = 2 * (h[i] + h[i + 1]);
            if(i == 0) A(0,1) = h[1];
            else if(i == N-2) A(N-2,N-3) = h[N-2];
            else{
                A(i,i-1) = h[i];  A(i,i+1) = h[i+1];
            }

        }
        Eigen::VectorXd t = A.colPivHouseholderQr().solve(v);
        u = Eigen::VectorXd::Zero(N+1);
        for(int i = 0; i < (int)t.size(); i++)u(i+1) = t(i);
    }

    double x, a, b, c, d, y;

    int ns = Curve[0]->ind, ne = Curve[N]->ind;
    for(int i = 0; i < N; i++){
        ns = std::min(Curve[i]->ind, Curve[i + 1]->ind);
        ne = std::max(Curve[i]->ind, Curve[i + 1]->ind);
        for(int k = ns + 1; k < ne; k++){
            //std::cout <<k << "  " << std::min(Curve[i]->ind, Curve[i + 1]->ind) + 1 << "  " << std::max(Curve[i]->ind, Curve[i + 1]->ind) << std::endl;
            double den = Curve[i + 1]->pt.x - Curve[i]->pt.x;
            a = (den != 0)? (u(i + 1) - u(i))/(6 * den) : 0;
            b = u(i)/2;
            c = - den * (2 * u(i) + u(i + 1))/6.0;
            c += (den != 0)? (color[i + 1] - color[i])/den : 0;
            d = color[i];
            //x =  crv->CurvePoints[k].pt.x - Curve[i]->pt.x;
            y = a * std::pow(x, 3) + b * std::pow(x, 2) + c * x + d;
            //crv->CurvePoints[k].color = y;
            //x = crv->CurvePoints[k].pt.x;
            CurvePath.push_back(glm::f64vec2{x,y});
        }
        //CurvePath.push_back(glm::f64vec2{crv->CurvePoints[ne].pt.x, color[i + 1]});
    }
    for(auto&f: Faces){
        for (auto& r : f->rulings) {
            f->Gradation = r->Gradation = r->pt->color;
        }

    }
}

void Model::setGradationValue(int val, int refMeshNum, int& color, int InterpolationType, std::vector<glm::f64vec2>& CurvePath){
    if(Faces.size() == 0 || refMeshNum == -1){std::cout<<"no selected" << std::endl; return;}
    color = -256;
    std::vector<glm::f64vec2> mesh;
    Face *f = Faces[refMeshNum];
    f->Gradation += val;
    f->Gradation = (f->Gradation < -255)? -255 : (255 < f->Gradation)? 255 : f->Gradation;
    for (auto& r : f->rulings) {
        r->Gradation = r->pt->color = f->Gradation;
    }
    color = f->Gradation;
    std::list<LinearRange> path;
    if(InterpolationType == 0){
        SetGradationPoint4linear(path, refMeshNum);
        LinearInterPolation(path);
    }
    if(InterpolationType == 1)SplineInterPolation(CurvePath);

    for(auto&f: Faces){
        for(auto&r: f->rulings) r->Gradation = f->Gradation;
    }

    return;
}

void Model::SetGradationPoint4linear(std::list<LinearRange>& path, int n){

    if(n == -1) return;
    if(linearrange[0].face == nullptr && linearrange[1].face == nullptr){ linearrange[0].face = Faces[n];}
    else if(linearrange[0].face != nullptr && linearrange[0].face != Faces[n] && linearrange[1].face == nullptr){linearrange[1].face = Faces[n];}
    if(linearrange[1].face == nullptr) return;
    makePath(linearrange[0].face, linearrange[1].face, path);
}

void Model::LinkRulingAndGradationArea(Face *f){        
    std::vector<glm::f64vec2> mesh;
    HalfEdge *h = f->halfedge;
    do{
        mesh.push_back(h->vertex->p);
        h = h->next;
    }while(h != f->halfedge);
    for(auto&c: crvs){
        if(c->isempty)continue;
        for(auto& r: c->Rulings){
            glm::f64vec2 p = (std::get<0>(r->r)->p + std::get<1>(r->r)->p);
            p /= 2;
            QPointF v{p.x, p.y};
            if(cn(mesh, v)){
                f->rulings.push_back(r);
                f->Gradation = r->Gradation = r->pt->color;
            }
        }
    }
    ruling *tmp;
    if(f->rulings.size() == 0){//safty的な感じ。面の中心とrulingの中点の距離が最小となるrulingを選択
        std::cout<<"safty fail in LinkRulingAndGradationArea"<<std::endl;
        std::vector<glm::f64vec2> mesh;
        HalfEdge *h = f->halfedge;
        do{
            mesh.push_back(h->vertex->p);
            h = h->next;
        }while(h != f->halfedge);
        glm::f64vec2 center{0,0};
        for(auto& v: mesh)center += v;
        center /= (double)mesh.size();
        double mindist = 100;
        for(auto& c: crvs){
            for(auto& r: c->Rulings){
                glm::f64vec3 p = std::get<0>(r->r)->p + std::get<1>(r->r)->p;
                p /= 2;
                if(glm::distance(p, glm::f64vec3{center,0}) < mindist){
                    mindist = glm::distance(p, glm::f64vec3{center,0});
                    tmp = r;
                }
            }
        }
        f->rulings.push_back(tmp);
    }
}

void Model::SetGradationArea(ruling *r, meshline *ml){
    Face *f;
}

void Model::deleteHE(HalfEdge *he){
    HalfEdge *prev = he->prev, *next = he->next;
    prev->next = next;
    next->prev = prev;
    Vertex *v = he->vertex;
    std::vector<HalfEdge*>::iterator itr_v = std::find(v->halfedge.begin(), v->halfedge.end(), he);
    if(itr_v != v->halfedge.end())v->halfedge.erase(itr_v);
    std::vector<HalfEdge*>::iterator itr = std::find(Edges.begin(), Edges.end(), he);
    if(itr != Edges.end())Edges.erase(itr);
}

void Model::replaceHE(HalfEdge *he){

}

void Model::addRulings(){
    if(!outline->isClosed)return;
    for(auto&c: crvs)c->CrossDetection();
    CrossDection4AllCurve();
    setOutline(outline->getVertices());
    for(auto& c: crvs){
        if(c->isempty)continue;
        for(auto l: c->meshLines){
            if(l->hasRulings[0]->IsCrossed != -1 || l->hasRulings[1]->IsCrossed != -1)continue;
            InsertVertex(l->vertices[0]);
            InsertVertex(l->vertices[1]);
            HalfEdge *he1 = new HalfEdge(l->vertices[0]);
            HalfEdge *he2 = new HalfEdge(l->vertices[1]);
            devide(he1, he2, Faces);
            Edges.push_back(he1); Edges.push_back(he2);
            vertices.push_back(l->vertices[0]);
            vertices.push_back(l->vertices[1]);
            std::vector<HalfEdge*>::iterator itr0 = Edges.end(); itr0 -= 2;
            std::vector<HalfEdge*>::iterator itr1 = Edges.end(); itr1 -= 1;
            l->vertices[0]->halfedge.push_back(*itr0);
            l->vertices[1]->halfedge.push_back(*itr1);
        }
    }

    for(auto&f : Faces)LinkRulingAndGradationArea(f);
    Fgrad.clear();
    for(auto&f: Faces){
        Fgrad.push_back(FaceGradation(f->halfedge, &f->Gradation));
    }

    /*
    for(auto& c: crvs){
        if(c->isempty)continue;
        for(auto ruling: c->Rulings){
            if(ruling->IsCrossed != -1 || ruling->IsCrossed != -1)continue;
            InsertVertex(std::get<0>(ruling->r));
            InsertVertex(std::get<1>(ruling->r));
            HalfEdge *he1 = new HalfEdge(std::get<0>(ruling->r));
            HalfEdge *he2 = new HalfEdge(std::get<1>(ruling->r));
            devide(he1, he2, Faces);
            Edges.push_back(he1); Edges.push_back(he2);
            vertices.push_back(std::get<0>(ruling->r));
            vertices.push_back(std::get<1>(ruling->r));
        }
    }

    for(auto&f : Faces)LinkRulingAndGradationArea(f);
    Fgrad.clear();
    for(auto&f: Faces){
        Fgrad.push_back(FaceGradation(f->halfedge, &f->Gradation));
    }
    */
}

void Model::updateRulings(){
    if(outline->getVertices().empty() || !outline->isClosed || crvs.empty())return;
    for(auto&c: crvs)c->CrossDetection();
    CrossDection4AllCurve();
    std::vector<Vertex*> V = outline->getVertices();
    int n =  V.size();
    std::vector<HalfEdge*>edges;
    for(int i = 0; i < (int)Fgrad.size(); i++){
        Fgrad[i].color = &Faces[i]->Gradation;
    }
    Faces.clear();
    vertices.clear();
    HalfEdge*he;
    for(auto& p: V){
        vertices.push_back(p);
        he = new HalfEdge(p);
        edges.push_back(he);
    }
    for(int i = 0; i < n; i++){
        edges[i]->prev = edges[(i + 1) % n];
        edges[(i + 1) % n]->next = edges[i];
    }

    Face *face = new Face(edges[0]);
    Faces.push_back(face);
    for(auto& he: edges)he->face = face;

    for(auto& crv: crvs){
        for(auto& l: crv->meshLines){
            if(l->hasRulings[0]->IsCrossed != -1 || l->hasRulings[1]->IsCrossed != -1)continue;
            l->vertices[0]->halfedge[0]->vertex = l->vertices[0];
            l->vertices[1]->halfedge[0]->vertex = l->vertices[1];
            vertices.push_back(l->vertices[0]->halfedge[0]->vertex);
            vertices.push_back(l->vertices[0]->halfedge[1]->vertex);
            ReInsertVertex(l->vertices[0]->halfedge[0], edges);
            ReInsertVertex(l->vertices[1]->halfedge[0], edges);
            edges.push_back(l->vertices[0]->halfedge[0]);
            edges.push_back(l->vertices[1]->halfedge[0]);
            devide(l->vertices[0]->halfedge[0], l->vertices[1]->halfedge[0], Faces);
        }
    }

    Edges = edges;
    for(auto&f : Faces) LinkRulingAndGradationArea(f);

    for(int i = 0; i < (int)Fgrad.size(); i++) Faces[i]->Gradation = *(Fgrad[i].color);
    for(auto&f: Faces){
        for(auto&r: f->rulings)r->Gradation = f->Gradation;
    }
}

void Model::SelectCurve(QPointF pt){
    int ptInd;
    int crvInd = searchPointIndex(pt, ptInd, 1);
    std::fill(refCrv.begin(), refCrv.end(), 0);
    if(crvInd != -1)refCrv[crvInd] = 1;
}

int Model::IsSelectedCurve(){
    auto itr = std::find(refCrv.begin(), refCrv.end(), 1);
    if(itr == refCrv.end())return -1;
    return std::distance(refCrv.begin(), itr);
}

int Model::getSelectedCurveIndex(QPointF pt){
    int ptInd;
    return searchPointIndex(pt, ptInd, 1);
}

//type = 0 -> Control Point, 1: Curve Point
int Model::searchPointIndex(QPointF pt, int& ptInd, int type){
    double dist = 5.0;
    ptInd = -1;
    int crvInd = -1;
    if(type == 0){
        for(int j = 0; j < (int)crvs.size(); j++){
            for(int i = 0; i < (int)crvs[j]->ControllPoints.size(); i++){
                glm::f64vec2 cp = crvs[j]->ControllPoints[i];
                if(dist * dist > (cp.x - pt.x()) * (cp.x - pt.x()) + (cp.y - pt.y()) * (cp.y - pt.y())){
                    dist = sqrt((cp.x - pt.x()) * (cp.x - pt.x()) + (cp.y - pt.y()) * (cp.y - pt.y()));
                    ptInd = i;
                    crvInd = j;
                }
            }
        }

    }else{
        for(int j = 0; j < (int)crvs.size(); j++){
            for(int i = 0; i < (int)crvs[j]->CurvePoints.size(); i++){
                glm::f64vec2 cp = crvs[j]->CurvePoints[i].pt;
                if(dist * dist > (cp.x - pt.x()) * (cp.x - pt.x()) + (cp.y - pt.y()) * (cp.y - pt.y())){
                    dist = sqrt((cp.x - pt.x()) * (cp.x - pt.x()) + (cp.y - pt.y()) * (cp.y - pt.y()));
                    ptInd = i;
                    crvInd = j;
                }
            }
        }
    }
    return crvInd;
}

int Model::DeleteCurve(){
    int DelIndex = IsSelectedCurve();
    if(DelIndex == -1)return -1;
    crvs.erase(crvs.begin() + DelIndex);
    refCrv.erase(refCrv.begin() + DelIndex);
    crvs.shrink_to_fit();
    refCrv.shrink_to_fit();
    return DelIndex;
}

void Model::DeleteControlPoint(QPointF pt, int curveDimention, int DivSize){
    int ptInd;
    int DelIndex = searchPointIndex(pt, ptInd, 0);
    if(DelIndex == -1)return;
    crvs[DelIndex]->ControllPoints.erase(crvs[DelIndex]->ControllPoints.begin() + ptInd);
    if(crvs[DelIndex]->getCurveType() == 0)crvs[DelIndex]->Bezier(3, crvPtNum);
    if(crvs[DelIndex]->getCurveType() == 1)crvs[DelIndex]->Bspline(3, crvPtNum);
    if(crvs[DelIndex]->ControllPoints.size() == 0){
        crvs.erase(crvs.begin() + DelIndex);
        refCrv.erase(refCrv.begin() + DelIndex);
    }
    crvs[DelIndex]->ControllPoints.shrink_to_fit();
    crvs.shrink_to_fit();
    refCrv.shrink_to_fit();
    if(outline->isClosed && crvs[DelIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){
        if(crvs[DelIndex]->getCurveType() == 0)crvs[DelIndex]->BezierRulings(outline, DivSize, crvPtNum);
        if(crvs[DelIndex]->getCurveType() == 1)crvs[DelIndex]->BsplineRulings(outline, DivSize, crvPtNum, curveDimention);
        addRulings();
    }
}

void Model::Check4Param(int curveDimention, std::vector<int>& deleteIndex){
    deleteIndex.clear();
    for(int i = (int)crvs.size() - 1; i >= 0; i--){
        if((crvs[i]->getCurveType() == 0 && (int)crvs[i]->ControllPoints.size() < curveDimention) ||
                (crvs[i]->getCurveType() == 1 && (int)crvs[i]->ControllPoints.size() <= curveDimention) ||
                (crvs[i]->getCurveType() == 2 && (int)crvs[i]->ControllPoints.size() < 2)) {
            deleteIndex.push_back(i);
            crvs.erase(crvs.begin() + i);
            refCrv.erase(refCrv.begin() + i);
        }

    }
    crvs.shrink_to_fit();
    refCrv.shrink_to_fit();
    linearrange[0] = LinearRange();
    linearrange[1] = LinearRange();
}

int Model::AddNewCurve(int curveType, int DivSize){
    CRV *crv = new CRV(crvPtNum, DivSize);
    crv->setCurveType(curveType);
    crvs.insert(crvs.begin(), crv);
    refCrv.insert(refCrv.begin(), 0);
    std::fill(refCrv.begin(), refCrv.end(), 0);
    refCrv[0] = 1;
    return 0;
}

void Model::AddControlPoint(QPointF pt, int curveDimention, int DivSize){

    int AddPtIndex = IsSelectedCurve();
    if(AddPtIndex == -1)return;
    if(crvs[AddPtIndex]->getCurveType() == 1){
        double dist = 5.0;
        int n = crvs[AddPtIndex]->movePtIndex(pt, dist);
        if(n == -1){
            crvs[AddPtIndex]->ControllPoints.push_back(glm::f64vec3{pt.x(), pt.y(), 0});
            crvs[AddPtIndex]->Bspline(curveDimention, crvPtNum);
        }
    }
    else if(crvs[AddPtIndex]->getCurveType() == 2){
        
        if(crvs[AddPtIndex]->ControllPoints.size() < 2)crvs[AddPtIndex]->ControllPoints.push_back(glm::f64vec3{pt.x(), pt.y(), 0});
        if(crvs[AddPtIndex]->ControllPoints.size() == 2)crvs[AddPtIndex]->Line();

    }
    if(outline->isClosed && crvs[AddPtIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){
        //if(crv->getCurveType() == 0)crvs[AddPtIndex]->BezierRulings(outline->vertices, DivSize, crvPtNum);
        if(crvs[AddPtIndex]->getCurveType() == 1)crvs[AddPtIndex]->BsplineRulings(outline, DivSize, crvPtNum, curveDimention);
        if(crvs[AddPtIndex]->getCurveType() == 2)crvs[AddPtIndex]->LineRulings(outline, DivSize);
        addRulings();
    }
}

void Model::MoveCurvePoint(QPointF pt, int MoveIndex, int ptInd, int curveDimention, int DivSize){
    if(MoveIndex == -1 || ptInd == -1)return;

    crvs[MoveIndex]->ControllPoints[ptInd] = glm::f64vec3{pt.x(), pt.y(),0};
    if(crvs[MoveIndex]->getCurveType() == 1)crvs[MoveIndex]->Bspline(curveDimention, crvPtNum);
    else if(crvs[MoveIndex]->getCurveType() == 2)crvs[MoveIndex]->Line();

    if(outline->isClosed && crvs[MoveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){
        //if(crvs[MoveIndex]->getCurveType() == 0)crvs[MoveIndex]->BezierRulings(outline->vertices, DivSize, crvPtNum);
        //if(crvs[MoveIndex]->getCurveType() == 1)crvs[MoveIndex]->BsplineRulings(outline->vertices, DivSize, crvPtNum, curveDimention);
        if(crvs[MoveIndex]->getCurveType() == 2)crvs[MoveIndex]->LineRulings(outline, DivSize);
        //addRulings();
    }
}

bool Model::CrossDection4AllCurve(){
    if(crvs.size() == 0)return false;
    for(int i = 0; i < (int)crvs.size(); i++){
        CRV *c1 = crvs[i];
        if(c1->isempty)continue;
        for(int j = i + 1; j < (int)crvs.size(); j++){
            CRV *c2 = crvs[j];
            if(c2->isempty)continue;
            for(auto& r1: c1->Rulings){
                for(auto& r2: c2->Rulings){
                    if(IsIntersect(std::get<0>(r1->r)->p, std::get<1>(r1->r)->p, std::get<0>(r2->r)->p, std::get<1>(r2->r)->p)){
                        r2->IsCrossed = 1;
                    }
                }
            }
        }
    }
    return true;
}

void Model::makePath(Face *start, Face *end, std::list<LinearRange>& path){
    path.clear();

    LinearRange l = LinearRange();
    l.face = start;
    path.push_back(l);

    HalfEdge *he; HalfEdge *NextFace;
    glm::f64vec3 EndPos = glm::f64vec3{0,0,0};
    he = end->halfedge;
    int cnt = 0;
    do{
        EndPos += he->vertex->p;
        cnt++; he = he->next;
    }while(end->halfedge != he);
    EndPos /= cnt;


    double dist;
    Face *face = start;
    while(face != end){
        he = face->halfedge;
        NextFace = nullptr;
        dist = 1000;
        do{
            if(he->pair != nullptr){
                //if(he->pair->face == end) break;
                glm::f64vec3 v = (he->vertex->p + he->next->vertex->p); v /= 2;
                if(glm::distance(v, EndPos) < dist){
                    NextFace = he;
                    dist = glm::distance(v, EndPos);
                }
            }
            he = he->next;
        }while(he != face->halfedge);
        if(NextFace != nullptr){
           face = NextFace->pair->face;
           l = LinearRange();
            l.face = NextFace->face;
            path.push_back(l);
        }else break;
    }

    return;
}
