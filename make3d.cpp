//#define PI 3.14159265359
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
    refFL.clear();
    Fgrad.clear();
    befFaceNum = 0;
}

void Model::clear(){
    vertices.clear();
    Edges.clear();
    Faces.clear();
    ol_vertices.clear();
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

glm::f64vec3 Model::SetOnGrid(QPointF& cursol, double gridsize){
    int x = (int)cursol.x() % (int)gridsize, y = (int)cursol.y() % (int)gridsize;
    x = (cursol.x() - x + gridsize/2);
    y = (cursol.y() - y + gridsize/2);
    return glm::f64vec3{x,y,0};
}

void Model::InsertVertex(Vertex *v){
    for(auto& f: Faces){
        HalfEdge *he = f->halfedge;
        do{
            if ((glm::dot(glm::normalize(he->vertex->p - v->p), glm::normalize(he->next->vertex->p - v->p))) <= -1 + FLT_EPSILON){
                HalfEdge *NewEdge = new HalfEdge(v, EdgeType::r);
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
        HalfEdge *he1 = new HalfEdge(v1, EdgeType::r);
        he1->r = r;
        HalfEdge *he2 = new HalfEdge(v2, EdgeType::r);
        he2->r = r;
        r->he[0] = he1; r->he[1] = he2;
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
        R = glm::rotate(he->r->Gradation * M_PI/255.0, axis);
        do{
            T = glm::translate(he->vertex->p - v1);          
            he->vertex->p3 = A * R * T * glm::f64vec4{0,0,0, 1};
            mesh.push_back(he->vertex->p3);
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

void Model::deform(){
    if(Faces.empty())return;
    glm::f64mat4x4 T, R, A;

    std::queue<glm::f64mat4x4> MatQue;
    std::queue<HalfEdge*> FacesQue;
    HalfEdge*Pos;
    glm::f64vec3 v1, v2;

    HalfEdge *he;
    for(auto& h: Edges)h->vertex->deformed = false;
    for(auto&f: Faces)f->bend = false;
    he = (*Faces.begin())->halfedge;
    do{
        T = glm::translate(he->vertex->p);
        he->vertex->p3 = T * glm::f64vec4{0,0,0, 1};
        if(he->pair != nullptr){
            FacesQue.push(he->pair);
            MatQue.push(T);
        }
        he = he->next;
        he->vertex->deformed = true;
    }while(he != (*Faces.begin())->halfedge);
    he->face->bend = true;

    while(!FacesQue.empty()){
        Pos = FacesQue.front(); A = MatQue.front();
        FacesQue.pop(); MatQue.pop();
        if(Pos->face->bend)continue;
        he = Pos;
        v1 = Pos->pair->vertex->p;
        v2 = Pos->vertex->p;
        glm::f64vec3 axis = glm::normalize(v2 - v1);
        double phi = he->r->Gradation * M_PI/255.0;
        R = glm::rotate(phi, axis);
        do{
            T = glm::translate(he->vertex->p - v1);
            he->vertex->p3 = A * R * T * glm::f64vec4{0,0,0, 1};
            he->vertex->deformed = true;
            if(he->pair != nullptr && !he->pair->face->bend){
                FacesQue.push(he->pair);
                MatQue.push(A * R * T);
            }
            he = he->next;
        }while(he != Pos);
        Pos->face->bend = true;
    }

}

void Model::applyFL(){
}

void modify2Druling(){

}

void Model::setOutline(){
    std::vector<Vertex*> _outline = outline->getVertices();
    int n = _outline.size();
    if(n < 3)return;
    clear();
    HalfEdge*he;
    for(auto& p: _outline){
        vertices.push_back(p);
        he = new HalfEdge(p, EdgeType::ol);
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

void Model:: addConstraint(QPointF& cursol, int type, int gridsize, glm::f64vec3 (&axis)[2]){
    if(outline->IsClosed()){
        std::cout<<"constraint can be applied only not closed outline"<<std::endl;
        return;
    }
    glm::f64vec3 p = SetOnGrid(cursol, gridsize);
    if(axis[0] == glm::f64vec3{-1,-1,0}){axis[0] = p; return;}
    else if(axis[1] == glm::f64vec3{-1,-1,0} && axis[0] != p)axis[1] = p;
    glm::f64vec3 V = glm::normalize(axis[1] - axis[0]);
    glm::f64vec3 N = glm::f64vec3{-V.y, V.x, 0};
    std::vector<Vertex*> SymPts;
    std::vector<Vertex*> Vertices = outline->getVertices();
    if(type == 0){
        for(auto&v: Vertices){
            double t = glm::length(glm::cross((v->p - axis[0]), (axis[0] - axis[1])))/glm::length(axis[0] - axis[1]);
            if(glm::dot(axis[0] - v->p, N) < 0) N *= -1;
            SymPts.push_back(new Vertex(v->p + 2 * t * N));
        }
    }
    //鏡映反転したことで作成した曲線の始点、終点が元の曲線の始点、終点のいずれかと十分近い場合接続する。
    //基準の距離 5px

    int n = Vertices.size();
    int IsConnected = 0;
    double dist = gridsize/2;
    if(glm::distance(SymPts[0]->p, Vertices[0]->p) < dist && glm::distance(SymPts[n-1]->p, Vertices[n-1]->p) < dist){
        for(int i = 1; i < (int)SymPts.size(); i++) outline->addVertex(SymPts[i], 0);
        IsConnected = 2;
    }
    else if(glm::distance(SymPts[0]->p, Vertices[0]->p) < dist && glm::distance(SymPts[n-1]->p, Vertices[n-1]->p) > dist){
        for(int i = 1; i <(int)SymPts.size(); i++) outline->addVertex(SymPts[i], 0);
        IsConnected = 1;
    }
    else if(glm::distance(SymPts[0]->p, Vertices[0]->p) > dist && glm::distance(SymPts[n-1]->p, Vertices[n-1]->p) < dist){
        for(int i = (int)SymPts.size() - 2; i >= 0; i--) outline->addVertex(SymPts[i], outline->getVertices().size());
        IsConnected = 1;
    }
    if(IsConnected == 0){
        ol_vertices.push_back(SymPts);
    }
    else if(IsConnected == 2){
        outline->ConnectEdges();
    }
}

void Model::drawOutline(QPointF& cursol, int drawtype, double gridsize, bool IsClicked){
    glm::f64vec3 p = SetOnGrid(cursol, gridsize);
    if(drawtype == 0 || drawtype == 1){
        if(outline->hasPtNum != 2)outline->addVertex(p);
    }else if(drawtype == 2) outline->drawPolygon(p, IsClicked);
    if(outline->IsClosed())setOutline();

}

void Model::editOutlineVertex(QPointF& cursol, double gridsize, int event){
    glm::f64vec3 p = SetOnGrid(cursol, gridsize);
    static int grabedOutlineVertex = -1;
    if(event == 0){
        float dist = 5;
        grabedOutlineVertex = -1;//0~: vertex
        std::vector<Vertex*> _vertices = outline->getVertices();
        for(int i = 0; i < (int)_vertices.size(); i++){
            if(glm::distance(_vertices[i]->p, p) < dist){
                grabedOutlineVertex = i; dist = glm::distance(_vertices[i]->p, p);
            }
        }
    }else if(event == 1){
        outline->MoveVertex(p, grabedOutlineVertex);
        if(outline->IsClosed())deform();
    }else if(event == 2)grabedOutlineVertex = -1;
}

void Model::ConnectOutline(QPointF& cursol, double gridsize){
    glm::f64vec3 p = SetOnGrid(cursol, gridsize);
    std::vector<Vertex*> V = outline->getVertices();
    for(auto&v: V){
        if(p == v->p){
            if(Connect2Vertices[0] == nullptr)Connect2Vertices[0] = v;
            else if(Connect2Vertices[1] == nullptr && Connect2Vertices[0] != v)Connect2Vertices[1] = v;
        }
    }
}

void Model::LinearInterPolation(std::vector<HalfEdge*>& path){
    if(path.size() < 2)return;

    glm::f64vec3 befcenter = glm::f64vec3{-1,-1,-1}, center;
    std::vector<glm::f64vec3> vertices;
    double len = 0.0;
    for(auto& l: path){
        if(befcenter == glm::f64vec3{-1,-1,-1}){         
            befcenter = (l->vertex->p + l->next->vertex->p)/2.0;
            continue;
        }     
        center = (l->vertex->p + l->next->vertex->p)/2.0;
        len += glm::distance(center, befcenter);
        befcenter = center;
    }
    double r = (GradationPoints[1]->r->Gradation - GradationPoints[0]->r->Gradation)/len;
    befcenter = glm::f64vec3{-1,-1,-1};
    HalfEdge *bef = nullptr;
    for(auto&l : path){
        if(bef == nullptr){
            bef = l;
            befcenter = (l->vertex->p + l->next->vertex->p)/2.0;
            continue;
        }
        center = (l->vertex->p + l->next->vertex->p)/2.0;
        l->r->Gradation = r * glm::distance(center, befcenter) + bef->r->Gradation;
        bef = l;
        befcenter = center;
    }
}

void Model::SplineInterPolation(std::vector<HalfEdge*>& path, std::vector<glm::f64vec2>& CurvePath){//凹型に関してはとりあえず虫
    if((int)GradationPoints.size() < 2)return;
    CurvePath.clear();
    int N =(int)GradationPoints.size() - 1;
    auto getCenter = [](HalfEdge *h) { return glm::f64vec3(h->vertex->p + h->next->vertex->p)/2.0;};
    auto getItr = [](std::vector<HalfEdge*>& _path, HalfEdge *h){ return (std::find(_path.begin(), _path.end(), h) != _path.end()) ? std::find(_path.begin(), _path.end(), h):  std::find(_path.begin(), _path.end(), h->pair);};
    Eigen::VectorXd v(N - 1);
    std::vector<double>h(N);
    for(int i = 1; i < N + 1; i++)h[i - 1] = glm::distance(getCenter(GradationPoints[i]),getCenter(GradationPoints[i - 1]));
    for(int i = 1; i < N; i++){
        double a = (h[i] != 0) ? (GradationPoints[i + 1]->r->Gradation - GradationPoints[i]->r->Gradation)/h[i] : 0, b = (h[i - 1] != 0) ? (GradationPoints[i]->r->Gradation - GradationPoints[i - 1]->r->Gradation)/h[i - 1]: 0;
        v(i - 1) = 6 * (a - b);
    }
    Eigen::VectorXd u;
    if(N == 1){
        u = Eigen::VectorXd::Zero(2);
    }
    if(N == 2){
        u = Eigen::VectorXd::Zero(3);
        double dx1, dx2, dx3;
        dx1 = glm::distance(getCenter(GradationPoints[2]), getCenter(GradationPoints[0]));
        dx2 = glm::distance(getCenter(GradationPoints[2]), getCenter(GradationPoints[1]));
        dx3 = glm::distance(getCenter(GradationPoints[1]), getCenter(GradationPoints[0]));
        if(abs(dx1) < FLT_EPSILON) u(1) = 0;
        else{
            double a = (dx2 == 0) ? 0: (GradationPoints[2]->r->Gradation - GradationPoints[1]->r->Gradation)/dx2;
            double b = (dx3 == 0) ? 0: (GradationPoints[1]->r->Gradation - GradationPoints[0]->r->Gradation)/dx3;
            u(1) = 3 * (a - b)/dx1;
        }
    }
    if(N > 2){
        Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N-1,N-1);
        for(int i = 0; i < N-1;i++){
            A(i,i) = 2 * (h[i] + h[i + 1]);
            if(i != 0){
                A(i-1,i) = A(i,i-1) = h[i];
            }

        }
        Eigen::VectorXd t = A.colPivHouseholderQr().solve(v);
        u = Eigen::VectorXd::Zero(N+1);
        for(int i = 0; i < (int)t.size(); i++)u(i+1) = t(i);
    }

    double x, a, b, c, d, y;

    int cnt = 0;
    for(int i = 0; i < N; i++){
        auto itr_cur = getItr(path, GradationPoints[i]);
        int cur = (itr_cur != path.end()) ? std::distance(path.begin(), itr_cur): -1;
        auto itr_next = getItr(path, GradationPoints[i + 1]);
        int next = (itr_next != path.end()) ? std::distance(path.begin(), itr_next): -1;
        if(cur == -1){
        std::cout<< i << "  cur"<<std::endl;
            return;
        }
        if(next == -1 ){
            std::cout<< i << "  next"<<std::endl;
            return;
        }
        glm::f64vec3 befcenter = getCenter(GradationPoints[i]);
        glm::f64vec3 center = getCenter(GradationPoints[i + 1]);
        double den = glm::distance(befcenter, center);
        a = (den != 0)? (u(i + 1) - u(i))/(6 * den) : 0;
        b = u(i)/2;
        c = - den * (2 * u(i) + u(i + 1))/6.0;
        c += (den != 0)? (GradationPoints[i + 1]->r->Gradation - GradationPoints[i]->r->Gradation)/den : 0;
        d = GradationPoints[i]->r->Gradation;
        for(int k = cur + 1; k < next; k++){
            x =  glm::distance(getCenter(path[k]), befcenter);
            y = a * std::pow(x, 3) + b * std::pow(x, 2) + c * x + d;
            path[k]->r->Gradation = y;
            CurvePath.push_back(glm::f64vec2{cnt++,y});
        }
    }
}

void Model::setGradationValue(int val, HalfEdge *refHE, int InterpolationType, std::vector<glm::f64vec2>& CurvePath){
    if(Faces.size() == 0 || refHE == nullptr || std::find(Edges.begin(), Edges.end(), refHE) == Edges.end()){std::cout<<"no selected" << std::endl; return;}
    refHE->r->Gradation += val;
    refHE->r->Gradation = (refHE->r->Gradation < -255)? -255 : (255 < refHE->r->Gradation)? 255 : refHE->r->Gradation;

    if(InterpolationType == 0){
        if(GradationPoints.size() == 0)GradationPoints.push_back(refHE);
        if(std::find(GradationPoints.begin(), GradationPoints.end(), refHE) == GradationPoints.end() && GradationPoints.size() < 2)GradationPoints.push_back(refHE);
        std::vector<HalfEdge*> path = makePath();
        LinearInterPolation(path);
    }
    if(InterpolationType == 1){
        if(GradationPoints.size() == 0)GradationPoints.push_back(refHE);
        if(std::find(GradationPoints.begin(), GradationPoints.end(), refHE) == GradationPoints.end())GradationPoints.push_back(refHE);

        std::vector<HalfEdge*> path = makePath();
        SplineInterPolation(path, CurvePath);
    }
    return;
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
                //f->rulings.push_back(r);
                //f->Gradation = r->Gradation = r->pt->color;
            }
        }
    }
    //ruling *tmp;
    //if(f->rulings.size() == 0){//safty的な感じ。面の中心とrulingの中点の距離が最小となるrulingを選択
        //std::cout<<"safty fail in LinkRulingAndGradationArea"<<std::endl;
        /*
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
        */
    //}
}

void Model::deleteHE(HalfEdge *he){
    HalfEdge *prev = he->prev, *next = he->next;
    prev->next = next;
    next->prev = prev;
    Vertex *v = he->vertex;
    std::vector<HalfEdge*>::iterator itr_v = std::find(v->halfedge.begin(), v->halfedge.end(), he);
    if(itr_v != v->halfedge.end())v->halfedge.erase(v->halfedge.begin() + std::distance(v->halfedge.begin(), itr_v));
    std::vector<HalfEdge*>::iterator itr = std::find(Edges.begin(), Edges.end(), he);
    delete he;
    if(itr != Edges.end())Edges.erase(itr);
}

void Model::ConnectEdge(HalfEdge *he){
    HalfEdge *prev = he->prev, *next = he->next;
    prev->next = next;
    next->prev = prev;
    std::vector<Vertex*>::iterator itr_v = std::find(vertices.begin(), vertices.end(), he->vertex);
    delete he->vertex;
    if(itr_v != vertices.end())vertices.erase(itr_v);
    delete he;
    std::vector<HalfEdge*>::iterator itr = std::find(Edges.begin(), Edges.end(), he);
    if(itr != Edges.end())Edges.erase(itr);

}

void Model::addRulings(){
    if(!outline->IsClosed())return;
    setOutline();
    for(auto&c: crvs)CrossDetection(outline->getFace(), c);
    CrossDection4AllCurve();

    for(auto& c: crvs){
        if(c->isempty)continue;
        for(auto rl: c->Rulings){
            if(rl->IsCrossed != -1)continue;
            InsertVertex(std::get<0>(rl->r));
            InsertVertex(std::get<1>(rl->r));
            HalfEdge *he1 = new HalfEdge(std::get<0>(rl->r), EdgeType::r);
            HalfEdge *he2 = new HalfEdge(std::get<1>(rl->r), EdgeType::r);
            he1->r = rl; he2->r = rl;
            rl->he[0] = he1; rl->he[1] = he2;
            devide(he1, he2, Faces);
            Edges.push_back(he1); Edges.push_back(he2);
            vertices.push_back(std::get<0>(rl->r));
            vertices.push_back(std::get<1>(rl->r));
            /*
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
            */
        }
    }
    /*
    for(auto&f : Faces)LinkRulingAndGradationArea(f);
    for(auto it = Faces.begin(); it != Faces.end();){
        if((*it)->rulings.empty()){
            ConnectFaces(*it, it);
            continue;
        }
        it++;
    }*/
    //Fgrad.clear();
    //for(auto&f: Faces){
        //Fgrad.push_back(FaceGradation(f->halfedge, &f->Gradation));
    //}

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
    if(outline->IsClosed() && crvs[DelIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){
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
                (crvs[i]->getCurveType() == 2 && (int)crvs[i]->ControllPoints.size() < 2) ||
                (crvs[i]->getCurveType() == 3 && (int)crvs[i]->ControllPoints.size() < 3)) {
            deleteIndex.push_back(i);
            crvs.erase(crvs.begin() + i);
            refCrv.erase(refCrv.begin() + i);
        }

    }
    crvs.shrink_to_fit();
    refCrv.shrink_to_fit();
    GradationPoints.clear();
    Axis4Const[0] = glm::f64vec3{-1,-1,0};
    Axis4Const[1] = glm::f64vec3{-1,-1,0};
    Connect2Vertices[0] = nullptr;
    Connect2Vertices[1] = nullptr;
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

void Model::AddControlPoint(glm::f64vec3& p, int curveDimention, int DivSize){
    int AddPtIndex = IsSelectedCurve();
    if(AddPtIndex == -1)return;
    if(crvs[AddPtIndex]->getCurveType() == 1){
        double dist = 5.0;
        int n = crvs[AddPtIndex]->movePtIndex(p, dist);
        if(n == -1){
            crvs[AddPtIndex]->ControllPoints.push_back(p);
            crvs[AddPtIndex]->Bspline(curveDimention, crvPtNum);
        }
    }
    else if(crvs[AddPtIndex]->getCurveType() == 2){        
        if(crvs[AddPtIndex]->ControllPoints.size() < 2)crvs[AddPtIndex]->ControllPoints.push_back(p);
        if(crvs[AddPtIndex]->ControllPoints.size() == 2)crvs[AddPtIndex]->Line();
    }
    else if(crvs[AddPtIndex]->getCurveType() == 3){
        if(crvs[AddPtIndex]->ControllPoints.size() == 0){
            if(outline->IsClosed()){
                Face *f = outline->getFace();
                bool PointOnLines = false;
                HalfEdge *h = f->halfedge;
                do{
                    if(glm::dot(glm::normalize(h->vertex->p - p),glm::normalize(h->next->vertex->p - p)) >= 1 - FLT_EPSILON)PointOnLines = true;
                    h = h->next;
                }while(h != f->halfedge);
                if(!cn(f, p) || PointOnLines)crvs[AddPtIndex]->ControllPoints.push_back(p);
            }
            else crvs[AddPtIndex]->ControllPoints.push_back(p);
        }
        else if(crvs[AddPtIndex]->ControllPoints.size() == 1)crvs[AddPtIndex]->ControllPoints.push_back(p);
        else if(crvs[AddPtIndex]->ControllPoints.size() == 2){
            double l = glm::distance(crvs[AddPtIndex]->ControllPoints[0], crvs[AddPtIndex]->ControllPoints[1]);
            double l2 = glm::distance(crvs[AddPtIndex]->ControllPoints[0], p);
            glm::f64vec3 point = l/l2 * (p - crvs[AddPtIndex]->ControllPoints[0]) + crvs[AddPtIndex]->ControllPoints[0];
            crvs[AddPtIndex]->ControllPoints.push_back(point);
        }
        if(crvs[AddPtIndex]->ControllPoints.size() == 3)crvs[AddPtIndex]->Arc(crvPtNum);

    }
    if(outline->IsClosed() && crvs[AddPtIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){
        //if(crv->getCurveType() == 0)crvs[AddPtIndex]->BezierRulings(outline->vertices, DivSize, crvPtNum);
        if(crvs[AddPtIndex]->getCurveType() == 1)crvs[AddPtIndex]->BsplineRulings(outline, DivSize, crvPtNum, curveDimention);
        if(crvs[AddPtIndex]->getCurveType() == 2)crvs[AddPtIndex]->LineRulings(outline, DivSize);
        if(crvs[AddPtIndex]->getCurveType() == 3)crvs[AddPtIndex]->ArcRulings(outline, DivSize);
        addRulings();
    }
}

void Model::MoveCurvePoint(glm::f64vec3& p, int MoveIndex, int ptInd, int curveDimention, int DivSize){
    if(MoveIndex == -1 || ptInd == -1)return;

    crvs[MoveIndex]->ControllPoints[ptInd] = p;
    if(crvs[MoveIndex]->getCurveType() == 1)crvs[MoveIndex]->Bspline(curveDimention, crvPtNum);
    else if(crvs[MoveIndex]->getCurveType() == 2)crvs[MoveIndex]->Line();
    else if(crvs[MoveIndex]->getCurveType() == 3)crvs[MoveIndex]->Arc(crvPtNum);

    if(outline->IsClosed() && crvs[MoveIndex]->CurvePoints[0].pt != glm::f64vec3{-1,-1,-1}){
        //if(crvs[MoveIndex]->getCurveType() == 0)crvs[MoveIndex]->BezierRulings(outline->vertices, DivSize, crvPtNum);
        //if(crvs[MoveIndex]->getCurveType() == 1)crvs[MoveIndex]->BsplineRulings(outline->vertices, DivSize, crvPtNum, curveDimention);
        if(crvs[MoveIndex]->getCurveType() == 2)crvs[MoveIndex]->LineRulings(outline, DivSize);
        if(crvs[MoveIndex]->getCurveType() == 3)crvs[MoveIndex]->ArcRulings(outline, DivSize);
        //addRulings();
    }
}

bool Model::AddControlPoint_FL(glm::f64vec3& p, int event, int curveDimention){
    bool res = false;
    if(event == 0){
        res = FL[0]->addCtrlPt(p, curveDimention, outline, Faces, Edges, vertices);
    }else if(event == 1){
        res = FL[0]->delCtrlPt(p, curveDimention, outline);
    }
    return res;
}

bool Model::CrossDection4AllCurve(){
    if(crvs.empty() || crvs[0]->isempty)return false;
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

std::vector<HalfEdge*> Model::makePath(){
    std::vector<HalfEdge*> path;
    if(GradationPoints.empty())return path;
    path.push_back(GradationPoints[0]);
    for(int i = 0; i < (int)GradationPoints.size() - 1; i++){
        Face *f = GradationPoints[i]->face;
        HalfEdge *he, *nextHE;
        glm::f64vec3 nextPos = (std::get<0>(GradationPoints[i+1]->r->r)->p + std::get<1>(GradationPoints[i+1]->r->r)->p)/2.0;
        do{
            he = f->halfedge;
            double dist = 1000;
            nextHE = nullptr;
            do{
                if(he->pair != nullptr){
                    double d = glm::length(glm::cross((nextPos - he->vertex->p), (he->vertex->p - he->pair->vertex->p)))/glm::length(he->vertex->p - he->pair->vertex->p);
                    if(d < dist){
                        dist = d;
                        nextHE = he;
                    }
                }
                he = he->next;
            }while(he != f->halfedge);
            if(nextHE != nullptr && std::find(Edges.begin(), Edges.end(), nextHE) != Edges.end()){
                bool haspath = false;
                for(auto& he: path){
                    if(he->r == nextHE->r){
                        haspath = true;
                        break;
                    }
                }
                if(!haspath)path.push_back(nextHE);
            }
            f = nextHE->pair->face;
        }while(nextHE->r != GradationPoints[i+1]->r);
    }
    return path;
}
