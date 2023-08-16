//#define PI 3.14159265359
#include "make3d.h"
using namespace MathTool;

Model::Model(){
    clear();
    outline = std::make_shared<OUTLINE>();
}

Model::Model(int _crvPtNum){
    crvPtNum = _crvPtNum;
    clear();
    outline = std::make_shared<OUTLINE>();
    refCrv.clear();
    refFL.clear();
    FoldCurveIndex = befFaceNum = 0;
    NTree_fl = NTree<std::shared_ptr<FoldLine>>();
}

void Model::clear(){
    //vertices.clear();
    Rulings.clear();
    //Faces.clear();
    ol_vertices.clear();
    ColorPt = ColorPoint(200, std::numbers::pi/2.0);
}

void Model::Initialize(){
    clear();
    outline = std::make_shared<OUTLINE>();
    refCrv.clear();
}

void Model::SetMaxFold(double val){
    ColorPt.angle = val * std::numbers::pi/180.0;
    if(!outline->IsClosed())return;
    for(auto& c: crvs){
        if(c->isempty){
            std::cout << "crvs has empty"<<std::endl;
            return;
        }
    }
    deform();
}


glm::f64vec3 Model::SetOnGrid(QPointF& cursol, double gridsize){
    int x = (int)cursol.x() % (int)gridsize, y = (int)cursol.y() % (int)gridsize;
    x = (cursol.x() - x + gridsize/2);
    y = (cursol.y() - y + gridsize/2);
    return glm::f64vec3{x,y,0};
}

void Model::deform(){
    auto Color2Angle = [](double a, ColorPoint CP){
        if(-CP.color <= a && a <= CP.color) return -a * CP.angle/CP.color;
        else if(CP.color < a) return -((std::numbers::pi - CP.angle)/(255.0 - CP.color)* (abs(a) - CP.color) + CP.angle);
        return ((std::numbers::pi - CP.angle)/(255.0 - CP.color)* (abs(a) - CP.color) + CP.angle);
    };

    if(Rulings.empty() || !outline->IsClosed())return;
    glm::f64mat4x4 T, R, M;
    std::vector<glm::f64mat4x4> TMs;
    std::vector<std::vector<std::shared_ptr<Vertex>>> Polygons;
    std::vector<std::shared_ptr<Vertex>> polygon;

    for(auto& l: outline->Lines) polygon.push_back(l->v);
    Polygons.push_back(polygon);
    for(auto itr_r = Rulings.begin(); itr_r != Rulings.end(); itr_r++){
        for(auto&P: Polygons){
            int vind = -1, oind = -1;
            for(int i = 0; i < (int)P.size(); i++){
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
    T = glm::translate(Rulings[0]->o->p);
    R = M = glm::f64mat4(1.0); TMs.push_back(M * R * T);
    for(int i = 1; i < (int)Rulings.size(); i++){
        glm::f64vec3 axis = glm::normalize(Rulings[i-1]->v->p - Rulings[i-1]->o->p);
        M = TMs.back();
        T = glm::translate(Rulings[i]->o->p - Rulings[i-1]->o->p);
        R = glm::rotate(Color2Angle(Rulings[i-1]->color, ColorPt), axis);
        Rulings[i]->o->p3_ori = Rulings[i]->o->p3 = M * R * T * glm::f64vec4{0,0,0, 1};
        TMs.push_back(M * R * T);
        T = glm::translate(Rulings[i]->v->p - Rulings[i-1]->o->p);
        Rulings[i]->v->p3_ori = Rulings[i]->v->p3 = M * R * T * glm::f64vec4{0,0,0, 1};
    }

    for(auto&l: outline->Lines){
        for(auto& P: Polygons){
            if(std::find(P.begin(), P.end(), l->v) == P.end())continue;
            for(int i = 0; i < (int)Rulings.size(); i++){
                if(std::find(P.begin(), P.end(), Rulings[i]->v) != P.end()){
                    if(i == 0){
                    }else if(i == Rulings.size() - 1){
                        glm::f64vec3 axis = glm::normalize(Rulings.back()->v->p - Rulings.back()->o->p);
                        M = TMs.back();
                        T = glm::translate(l->v->p - Rulings.back()->o->p);
                        R = glm::rotate(Color2Angle(Rulings.back()->color, ColorPt), axis);
                        l->v->p3_ori = l->v->p3 = M * R * T * glm::f64vec4{0,0,0, 1};
                        break;
                    }else{
                        T = glm::translate(l->v->p - Rulings[i]->o->p);
                        l->v->p3_ori = l->v->p3 = TMs[i] * T * glm::f64vec4{0,0,0, 1};
                        break;
                    }
                }
            }
        }
    }
}

void Model::ChangeFoldLineState(){
    if((int)FL.size() == 1)return;
    FoldCurveIndex++;
    for(auto itr = FL.begin(); itr != FL.end();){
        if((*itr)->FoldingCurve.empty() && std::distance(FL.begin(), itr) < FoldCurveIndex){
            itr = FL.erase(itr);
        }
        else itr++;
    }
    //vertices.shrink_to_fit();
    FoldCurveIndex = FL.size() - 1;
}

void Model::UpdateFLOrder(int dim){
    //下から上へとn分木での実装が必要かも
    std::vector<std::shared_ptr<FoldLine>> hasFoldingCurve;
    for(auto&fl: FL){
        if((int)fl->CtrlPts.size() > dim)hasFoldingCurve.push_back(fl);
    }

    std::shared_ptr<Line> btm = outline->Lines.front();//一番下の辺を探索
    for(auto&l: outline->Lines){
        if(((l->v->p + l->o->p)/2.0).y > ((btm->v->p + btm->o->p)/2.0).y)btm = l;
    }
    int btm_i = std::distance(outline->Lines.begin(), std::find(outline->Lines.begin(), outline->Lines.end(),btm));
    int i = btm_i;
    std::list<std::shared_ptr<FoldLine>> List_FL;//frontかbackのいずれかが探索するline上にある時、先頭と違っていれば先頭の要素を親に新しいものを子としてNTree_flに挿入、同じであればqueueから吐き出す
    //折曲線が交差しないという前提ならこの実装方法で正しいはず
    class LineOnFL{
    public:
        std::shared_ptr<FoldLine> FL;
        double t;
        LineOnFL(std::shared_ptr<FoldLine> _FL, double _t): FL(_FL), t(_t){}
    };
    NTree_fl.clear();
    do{
       std::vector<LineOnFL> LoF;
        for(auto&fl: hasFoldingCurve){
            if(outline->Lines[i]->is_on_line(fl->FoldingCurve.front().first->p))
                LoF.push_back(LineOnFL(fl, glm::length(fl->FoldingCurve.front().first->p - outline->Lines[i]->o->p)/glm::length(outline->Lines[i]->v->p - outline->Lines[i]->o->p)));
            if(outline->Lines[i]->is_on_line(fl->FoldingCurve.back().first->p))
                LoF.push_back(LineOnFL(fl, glm::length(fl->FoldingCurve.back().first->p - outline->Lines[i]->o->p)/glm::length(outline->Lines[i]->v->p - outline->Lines[i]->o->p)));
        }
        if(!LoF.empty()){
            std::sort(LoF.begin(), LoF.end(), [](const LineOnFL& a, const LineOnFL& b){return a.t < b.t;});
            for(auto&x: LoF){
                if(List_FL.empty()){
                    List_FL.push_back(x.FL);
                    if(NTree_fl.empty())NTree_fl = NTree(x.FL);
                }else{
                    auto res = NTree_fl.find(x.FL);
                    if(res == nullptr){
                        NTree_fl.insert(List_FL.front(), x.FL);
                        List_FL.push_front(x.FL);
                    }else{
                        auto it = std::find(List_FL.begin(), List_FL.end(), x.FL);
                        if(it != List_FL.end()) List_FL.erase(it);
                    }
                }
            }

        }
        i = (i + 1) % (int)outline->Lines.size();
    }while(i != btm_i);

    NTree_fl.print();
}

bool Model::BendingModel(double wb, double wp, int dim, bool ConstFunc){
    UpdateFLOrder(dim);
    std::shared_ptr<NTreeNode<std::shared_ptr<FoldLine>>> root = NTree_fl.GetRoot();
    if(root == nullptr)return false;
    std::queue<std::shared_ptr<NTreeNode<std::shared_ptr<FoldLine>>>> q;
    std::vector<std::shared_ptr<Vertex>> Poly_V = outline->getVertices();
    //root->data->Optimization_FlapAngle(Poly_V, wb, wp, ConstFunc);
    q.push(root);
    while(!q.empty()){
        auto cur = q.front(); q.pop();
        cur->data->Optimization_FlapAngle(Poly_V, wb, wp, ConstFunc);
        for (const auto& child : cur->children){
            if(child != nullptr){
                child->data->reassinruling(cur->data);
                q.push(child);
            }
        }
    }
    return true;
}

void Model::applyAAAMethod(double a){
    auto Poly_V = outline->getVertices();
    FL[FoldCurveIndex]->applyAAAMethod(Poly_V, a);
}

bool Model::RevisionCrosPtsPosition(){return FL[FoldCurveIndex]->RevisionCrosPtsPosition();}

bool Model::SplitRulings(int dim){
    using namespace MathTool;
    glm::f64vec3 UpVec{0,-1,0};
    auto getCrossPoint = [](std::vector<glm::f64vec3>& CtrlPts,  const std::shared_ptr<Vertex>& v, const std::shared_ptr<Vertex>& o, int dim){
        std::vector<double>arcT = BezierClipping(CtrlPts, v, o, dim);
        for(auto&t: arcT){
            if(t < 0 || 1 < t){std::cout<<"t is not correct value " << t << std::endl; continue;}
            glm::f64vec3 v2{0,0,0};
            for (int i = 0; i < int(CtrlPts.size()); i++) v2 += MathTool::BernsteinBasisFunc(dim, i, t) * CtrlPts[i];
            if(!MathTool::is_point_on_line(v2, v->p, o->p))continue;
            double sa = glm::distance(v2, o->p), sc = glm::distance(o->p, v->p);
            glm::f64vec3 v3 = sa/sc * (v->p3 - o->p3) + o->p3;
            std::shared_ptr<CrvPt_FL> P = std::make_shared<CrvPt_FL>(CrvPt_FL(v2, v3, t));
            P->set(v2, o, v);
            return P;
        }
        return std::shared_ptr<CrvPt_FL>(nullptr);
    };
    if(FL.empty() || FL[FoldCurveIndex]->CtrlPts.size() <= dim)return false;
    /*
    if(FL.size() == 1){
        for(auto& r: Rulings){
            CrvPt_FL *P = getCrossPoint(FL[0]->CtrlPts, r->v, r->o, dim);
            if(P!= nullptr){
                if(glm::dot(UpVec, glm::normalize(r->v->p - r->o->p)) > 0)FL[0]->FoldingCurve.push_back(Vertex4d(P, r->v, r->o));
                else FL[0]->FoldingCurve.push_back(Vertex4d(P, r->o, r->v));
            }
        }
        for(auto& l: outline->Lines){
            CrvPt_FL *P = getCrossPoint(FL[0]->CtrlPts, l->v, l->o, dim);
            if(P!= nullptr){
                if(glm::dot(UpVec, glm::normalize(l->v->p - l->o->p)) > 0)FL[0]->FoldingCurve.push_back(Vertex4d(P, l->v, l->o));
                else FL[0]->FoldingCurve.push_back(Vertex4d(P, l->o, l->v));
            }
        }
        FL[0]->SortCurve();
    }else{
        for(auto&fl: FL){

            if(fl == FL[FoldCurveIndex])continue;
            for(auto& r: fl->FoldingCurve){
                CrvPt_FL *P = getCrossPoint(FL[FoldCurveIndex]->CtrlPts, r.second, r.first, dim);
                if(P!= nullptr){
                    FL[FoldCurveIndex]->FoldingCurve.push_back(Vertex4d(P, r.second, r.first)); r.second = P;
                    continue;
                }
                P = getCrossPoint(FL[FoldCurveIndex]->CtrlPts, r.third, r.first, dim);
                if(P!= nullptr){
                    FL[FoldCurveIndex]->FoldingCurve.push_back(Vertex4d(P, r.third, r.first)); r.third = P;
                    continue;
                }
            }
            FL[FoldCurveIndex]->SortCurve();
        }
    }*/

    for(auto&fl: FL){
        fl->FoldingCurve.clear();
        for(auto& l: outline->Lines){
            std::shared_ptr<CrvPt_FL> P = getCrossPoint(fl->CtrlPts, l->v, l->o, dim);
            if(P!= nullptr){
                if(glm::dot(UpVec, glm::normalize(l->v->p - l->o->p)) > 0)fl->FoldingCurve.push_back(Vertex4d(P, l->v, l->o));
                else fl->FoldingCurve.push_back(Vertex4d(P, l->o, l->v));
            }
        }
        fl->SortCurve();
    }

    UpdateFLOrder(dim);
    auto root = NTree_fl.GetRoot();
    if(root == nullptr)return false;
    for(auto& r: Rulings){
        std::shared_ptr<CrvPt_FL> P = getCrossPoint(root->data->CtrlPts, r->v, r->o, dim);
        if(P!= nullptr){
            if(glm::dot(UpVec, glm::normalize(r->v->p - r->o->p)) > 0)root->data->FoldingCurve.push_back(Vertex4d(P, r->v, r->o));
            else root->data->FoldingCurve.push_back(Vertex4d(P, r->o, r->v));
        }
    }
    root->data->SortCurve();
    /*
    std::queue<std::shared_ptr<NTreeNode<std::shared_ptr<FoldLine>>>>q;
    q.push(root);
    while (!q.empty()) {
        auto cur = q.front(); q.pop();
        for (const auto& child : cur->children){
            if(child != nullptr){
                child->data->reassinruling(cur->data);
                q.push(child);
            }

        }
    }*/
    return true;
}

void Model::modifyFoldingCurvePositionOn3d(){
    for(auto&FC: FL){
        for(auto&fldCrv: FC->FoldingCurve){
            for(auto&r: Rulings){
                if(!MathTool::is_point_on_line(fldCrv.first->p, r->o->p, r->v->p))continue;
                double t = glm::length(fldCrv.first->p - r->o->p)/glm::length(r->o->p - r->v->p);
                fldCrv.first->p3_ori = fldCrv.first->p3 = t * (r->v->p3 - r->o->p3) + r->o->p3;
                break;
            }
            for(auto&l: outline->Lines){
                if(!MathTool::is_point_on_line(fldCrv.first->p, l->o->p, l->v->p))continue;
                double t = glm::length(fldCrv.first->p - l->o->p)/glm::length(l->o->p - l->v->p);
                fldCrv.first->p3_ori = fldCrv.first->p3 = t * (l->v->p3 - l->o->p3) + l->o->p3;
                break;
            }
        }
    }
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
    std::vector<std::shared_ptr<Vertex>> SymPts;
    std::vector<std::shared_ptr<Vertex>> Vertices = outline->getVertices();
    if(type == 0){
        for(auto&v: Vertices){
            double t = glm::length(glm::cross((v->p - axis[0]), (axis[0] - axis[1])))/glm::length(axis[0] - axis[1]);
            if(glm::dot(axis[0] - v->p, N) < 0) N *= -1;
            SymPts.push_back(std::make_shared<Vertex>(v->p + 2 * t * N));
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
    else if(IsConnected == 2){//元々あったedgeを削除してつなぎなおす
        std::cout<<"can't use now"<<std::endl;
    }
}

void Model::drawOutline(QPointF& cursol, int drawtype, double gridsize, bool IsClicked){
    glm::f64vec3 p = SetOnGrid(cursol, gridsize);
    if(drawtype == 0 || drawtype == 1){
        if(outline->hasPtNum != 2)outline->addVertex(p);
    }else if(drawtype == 2) outline->drawPolygon(p, IsClicked);

}

void Model::editOutlineVertex(QPointF& cursol, double gridsize, int event){
    glm::f64vec3 p = SetOnGrid(cursol, gridsize);
    static int grabedOutlineVertex = -1;
    if(event == 0){
        float dist = 5;
        grabedOutlineVertex = -1;//0~: vertex
        std::vector<std::shared_ptr<Vertex>> _vertices = outline->getVertices();
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
    std::vector<std::shared_ptr<Vertex>> V = outline->getVertices();
    for(auto&v: V){
        if(p == v->p){
            if(Connect2Vertices[0] == nullptr)Connect2Vertices[0] = v;
            else if(Connect2Vertices[1] == nullptr && Connect2Vertices[0] != v)Connect2Vertices[1] = v;
        }
    }
}

void Model::LinearInterPolation(const std::vector<std::shared_ptr<Line>>& path){
    if(path.size() < 2)return;

    glm::f64vec3 befcenter = glm::f64vec3{-1,-1,-1}, center;
    double len = 0.0;
    for(auto& l: path){
        if(befcenter == glm::f64vec3{-1,-1,-1}){         
            befcenter = (l->v->p + l->o->p)/2.0;
            continue;
        }     
        center = (l->v->p + l->o->p)/2.0;
        len += glm::distance(center, befcenter);
        befcenter = center;
    }
    double r = (GradationPoints[1]->color - GradationPoints[0]->color)/len;
    befcenter = glm::f64vec3{-1,-1,-1};
    std::shared_ptr<Line> bef = std::shared_ptr<Line>(nullptr);
    for(auto&l : path){
        if(bef == nullptr){
            bef = l;
            befcenter = (l->v->p + l->o->p)/2.0;
            continue;
        }
        center = (l->v->p + l->o->p)/2.0;
        l->color = r * glm::distance(center, befcenter) + bef->color;
        bef = l;
        befcenter = center;
    }
}

void Model::SplineInterPolation(const std::vector<std::shared_ptr<Line>>& path, std::vector<glm::f64vec2>& CurvePath){//凹型に関してはとりあえず虫
    if((int)GradationPoints.size() < 2)return;
    CurvePath.clear();
    int N =(int)GradationPoints.size() - 1;
    auto getCenter = [](const std::shared_ptr<Line>& L) { return glm::f64vec3(L->o->p + L->v->p)/2.0;};
    Eigen::VectorXd v(N - 1);
    std::vector<double>h(N);
    for(int i = 1; i < N + 1; i++)h[i - 1] = glm::distance(getCenter(GradationPoints[i]),getCenter(GradationPoints[i - 1]));
    for(int i = 1; i < N; i++){
        double a = (h[i] != 0) ? (GradationPoints[i + 1]->color - GradationPoints[i]->color)/h[i] : 0, b = (h[i - 1] != 0) ? (GradationPoints[i]->color - GradationPoints[i - 1]->color)/h[i - 1]: 0;
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
        if(abs(dx1) < DBL_EPSILON) u(1) = 0;
        else{
            double a = (dx2 == 0) ? 0: (GradationPoints[2]->color - GradationPoints[1]->color)/dx2;
            double b = (dx3 == 0) ? 0: (GradationPoints[1]->color - GradationPoints[0]->color)/dx3;
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
        auto itr_cur = std::find(path.begin(), path.end(), GradationPoints[i]);
        int cur = (itr_cur != path.end()) ? std::distance(path.begin(), itr_cur): -1;
        auto itr_next = std::find(path.begin(), path.end(), GradationPoints[i+1]);
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
        c += (den != 0)? (GradationPoints[i + 1]->color - GradationPoints[i]->color)/den : 0;
        d = GradationPoints[i]->color;
        for(int k = cur + 1; k < next; k++){
            x =  glm::distance(getCenter(path[k]), befcenter);
            y = a * std::pow(x, 3) + b * std::pow(x, 2) + c * x + d;
            path[k]->color = y;
            CurvePath.push_back(glm::f64vec2{cnt++,y});
        }
    }
}

void Model::setGradationValue(int val, const std::shared_ptr<Line>& refL, int InterpolationType, std::vector<glm::f64vec2>& CurvePath){
    //if(Faces.size() == 0 || std::find(Edges.begin(), Edges.end(), refHE) == Edges.end()){std::cout<<"no selected" << std::endl; return;
    if(Rulings.empty() || std::find(Rulings.begin(), Rulings.end(), refL) == Rulings.end()){std::cout<<"no selected" << std::endl; return;}
    if(refL->et != EdgeType::r)return;
    refL->color += val;
    refL->color = (refL->color < -255.0)? -255.0 : (255.0 < refL->color)? 255.0 : refL->color;

    if(InterpolationType == 0){
        if(GradationPoints.size() == 0)GradationPoints.push_back(refL);
        if(std::find(GradationPoints.begin(), GradationPoints.end(), refL) == GradationPoints.end() && GradationPoints.size() < 2)GradationPoints.push_back(refL);
        if(GradationPoints.empty())return;
        std::vector<std::shared_ptr<Line>> path;
        if(GradationPoints.size() == 1) path = {refL};
        else{
            int i0 = std::distance(Rulings.begin(), std::find(Rulings.begin(), Rulings.end(), GradationPoints[0]));
            int i1 = std::distance(Rulings.begin(), std::find(Rulings.begin(), Rulings.end(), GradationPoints[1]));
            path = {Rulings.begin() + std::min(i0, i1), Rulings.begin() + std::max(i0, i1)+1};
        }
        LinearInterPolation(path);//単一のruling制御曲線からのみ操作するものとする
    }
    if(InterpolationType == 1){
        //if(GradationPoints.size() == 0)GradationPoints.push_back(refL);
        //if(std::find(GradationPoints.begin(), GradationPoints.end(), refL) == GradationPoints.end())GradationPoints.push_back(refHE);
        //std::vector<Line*> path = makePath();
        SplineInterPolation(Rulings, CurvePath);
    }
    return;
}


void Model::addRulings(){
    if(!outline->IsClosed())return;
    Rulings.clear();
    for(auto&c: crvs)CrossDetection(outline, c);
    CrossDection4AllCurve();
    for(auto&c: crvs){
        for(auto&r: c->Rulings){
            if(r->IsCrossed == -1){
                this->Rulings.push_back(r);
            }
        }
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
                glm::f64vec2 cp = crvs[j]->CurvePoints[i];
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
    bool res = false;
    crvs[DelIndex]->ControllPoints.erase(crvs[DelIndex]->ControllPoints.begin() + ptInd);
    //if(crvs[DelIndex]->getCurveType() == CurveType::bezier3)crvs[DelIndex]->drawBezier(3, crvPtNum);
    if(crvs[DelIndex]->getCurveType() == CurveType::bsp3)res = crvs[DelIndex]->drawBspline(3, crvPtNum);
    if(crvs[DelIndex]->ControllPoints.size() == 0){
        crvs.erase(crvs.begin() + DelIndex);
        refCrv.erase(refCrv.begin() + DelIndex);
        res = false;
    }
    crvs[DelIndex]->ControllPoints.shrink_to_fit();
    crvs.shrink_to_fit();
    refCrv.shrink_to_fit();
    if(outline->IsClosed() && res){
        if(crvs[DelIndex]->getCurveType() == CurveType::bezier3)crvs[DelIndex]->BezierRulings(outline, DivSize, crvPtNum);
        if(crvs[DelIndex]->getCurveType() == CurveType::bsp3)crvs[DelIndex]->BsplineRulings(outline, DivSize, crvPtNum, curveDimention);
        addRulings();
    }
}

void Model::Check4Param(int curveDimention, std::vector<int>& deleteIndex){
    deleteIndex.clear();
    for(int i = (int)crvs.size() - 1; i >= 0; i--){
        if((crvs[i]->getCurveType() == CurveType::bezier3 && (int)crvs[i]->ControllPoints.size() < curveDimention) ||
                (crvs[i]->getCurveType() == CurveType::bsp3 && (int)crvs[i]->ControllPoints.size() <= curveDimention) ||
                (crvs[i]->getCurveType() == CurveType::line && (int)crvs[i]->ControllPoints.size() < 2) ||
                (crvs[i]->getCurveType() == CurveType::arc && (int)crvs[i]->ControllPoints.size() < 3)) {
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
    Connect2Vertices[0] = std::shared_ptr<Vertex>(nullptr);
    Connect2Vertices[1] = std::shared_ptr<Vertex>(nullptr);
}

int Model::AddNewCurve(CurveType curveType, int DivSize){
    std::shared_ptr<CRV> crv = std::make_shared<CRV>(crvPtNum, DivSize);
    crv->setCurveType(curveType);
    crvs.insert(crvs.begin(), crv);
    refCrv.insert(refCrv.begin(), 0);
    std::fill(refCrv.begin(), refCrv.end(), 0);
    refCrv[0] = 1;
    return 0;
}

bool Model::AddControlPoint(glm::f64vec3& p, int curveDimention, int DivSize){
    int AddPtIndex = IsSelectedCurve();
    if(AddPtIndex == -1)return false;
    bool res = false;
    if(crvs[AddPtIndex]->getCurveType() == CurveType::bsp3){
        double dist = 5.0;
        int n = crvs[AddPtIndex]->movePtIndex(p, dist);
        if(n == -1){
            crvs[AddPtIndex]->ControllPoints.push_back(p);
            res = crvs[AddPtIndex]->drawBspline(curveDimention, crvPtNum);
        }
    }
    else if(crvs[AddPtIndex]->getCurveType() == CurveType::line){
        if(crvs[AddPtIndex]->ControllPoints.size() < 2)crvs[AddPtIndex]->ControllPoints.push_back(p);
        if(crvs[AddPtIndex]->ControllPoints.size() == 2)res = crvs[AddPtIndex]->drawLine();
    }
    else if(crvs[AddPtIndex]->getCurveType() == CurveType::arc){
        if(crvs[AddPtIndex]->ControllPoints.size() == 0){
            if(outline->IsClosed()){
                bool PointInFace = outline->IsPointInFace(p);
                bool PointOnLines = false;
                for(auto&l : outline->Lines){if(l->is_on_line(p))PointOnLines = true;}
                if(!PointInFace || PointOnLines)crvs[AddPtIndex]->ControllPoints.push_back(p);
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
        if(crvs[AddPtIndex]->ControllPoints.size() == 3)res = crvs[AddPtIndex]->drawArc(crvPtNum);

    }
    if(outline->IsClosed() && res){
        //if(crv->getCurveType() == 0)crvs[AddPtIndex]->BezierRulings(outline->vertices, DivSize, crvPtNum);
        if(crvs[AddPtIndex]->getCurveType() == CurveType::bsp3)crvs[AddPtIndex]->BsplineRulings(outline, DivSize, crvPtNum, curveDimention);
        if(crvs[AddPtIndex]->getCurveType() == CurveType::line)crvs[AddPtIndex]->LineRulings(outline, DivSize);
        if(crvs[AddPtIndex]->getCurveType() == CurveType::arc)crvs[AddPtIndex]->ArcRulings(outline, DivSize);
        //addRulings();
        return true;
    }
    return false;
}

void Model::MoveCurvePoint(glm::f64vec3& p, int MoveIndex, int ptInd, int curveDimention, int DivSize){
    if(MoveIndex == -1 || ptInd == -1)return;
    bool res = false;
    if(crvs[MoveIndex]->getCurveType() == CurveType::arc && ptInd == 0){
        bool PointOnLines = false;
        bool PointInFace = outline->IsPointInFace(p);
        for(auto&l : outline->Lines){if(is_point_on_line(p, l->o->p, l->v->p))PointOnLines = true;}
        if(!PointInFace || PointOnLines) crvs[MoveIndex]->ControllPoints[ptInd] = p;
    }
    else crvs[MoveIndex]->ControllPoints[ptInd] = p;
    if(crvs[MoveIndex]->getCurveType() == CurveType::bsp3)res = crvs[MoveIndex]->drawBspline(curveDimention, crvPtNum);
    else if(crvs[MoveIndex]->getCurveType() == CurveType::line)res = crvs[MoveIndex]->drawLine();
    else if(crvs[MoveIndex]->getCurveType() == CurveType::arc)res = crvs[MoveIndex]->drawArc(crvPtNum);

    if(outline->IsClosed() && res){
        //if(crvs[MoveIndex]->getCurveType() == 0)crvs[MoveIndex]->BezierRulings(outline->vertices, DivSize, crvPtNum);
        //if(crvs[MoveIndex]->getCurveType() == 1)crvs[MoveIndex]->BsplineRulings(outline->vertices, DivSize, crvPtNum, curveDimention);
        if(crvs[MoveIndex]->getCurveType() == CurveType::line)crvs[MoveIndex]->LineRulings(outline, DivSize);
        if(crvs[MoveIndex]->getCurveType() == CurveType::arc)crvs[MoveIndex]->ArcRulings(outline, DivSize);
    }
}

bool Model::AddControlPoint_FL(glm::f64vec3& p, int event, int curveDimention){
    bool res = false;
    if(event == 0){
        res = FL[FoldCurveIndex]->addCtrlPt(p, curveDimention);
    }else if(event == 1){
        res = FL[FoldCurveIndex]->delCtrlPt(p, curveDimention, outline);
    }
    return res;
}

bool Model::CrossDection4AllCurve(){
    if(crvs.empty() || crvs[0]->isempty)return false;
    for(int i = 0; i < (int)crvs.size(); i++){
        std::shared_ptr<CRV> c1 = crvs[i];
        if(c1->isempty)continue;
        for(int j = i + 1; j < (int)crvs.size(); j++){
            std::shared_ptr<CRV> c2 = crvs[j];
            if(c2->isempty)continue;
            for(auto& r1: c1->Rulings){
                for(auto& r2: c2->Rulings){
                    if(IsIntersect(r1->v->p, r1->o->p, r2->v->p, r2->o->p)){
                        r2->IsCrossed = 1;
                    }
                }
            }
        }
    }
    return true;
}

/*
std::vector<Line*> Model::makePath(){
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
               if(nextHE != nullptr && std::find(Rulings.begin(), Edges.end(), nextHE) != Edges.end()){
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
}*/
