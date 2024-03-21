#include "make3d.h"

Model::Model(){
    clear();
    outline = std::make_shared<OUTLINE>();
}

Model::Model(int _crvPtNum){
    crvPtNum = _crvPtNum;
    clear();
    outline = std::make_shared<OUTLINE>();
    refCrv.clear();
    refCreases = nullptr;
    movingVertexIndex = -1;
    NTree_Creases = std::shared_ptr<FoldLine>();
}

std::shared_ptr<Model> Model::deepCopy(){
    std::shared_ptr<Model> m = std::make_shared<Model>();
    m->crvPtNum = crvPtNum;
    m->refCrv = refCrv;
    m->ColorPt = ColorPt;
    m->Rulings.resize(Rulings.size());
    for(int i = 0; i < (int)Rulings.size(); i++)m->Rulings[i] = Rulings[i]->deepCopy();
    m->outline = outline->deepCopy();
    m->RulingCurve.resize(RulingCurve.size());
    for(int i = 0; i < (int)RulingCurve.size(); i++)m->RulingCurve[i] = RulingCurve[i]->deepCopy();
    return m;
}

void Model::clear(){
    Rulings.clear();
    ColorPt = ColorPoint(200, std::numbers::pi/2.0);
}

void Model::Initialize(){
    clear();
    outline = std::make_shared<OUTLINE>();
    refCrv.clear();
}

//クリックで選択した曲線の決定
void Model::detectClickedObj(const QPointF& curPos){
    refCreases = nullptr;
    const double step = 0.002;
    double mindist = 7.0;
    Eigen::Vector3d pos{curPos.x(), curPos.y(), 0};
    std::queue<std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>>> q;
    for(auto&btmcrease: NTree_Creases.GetRoot()->children)q.push(btmcrease);
    while (!q.empty()) {
        std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>> cur = q.front(); q.pop();
        for (auto child : cur->children){if (child != nullptr)q.push(child);}
        if(cur->data->CtrlPts.size() <= 3)continue;
        for(double t = 0.0; t <= 1.0; t += step){
            Eigen::Vector3d p = MathTool::bezier(cur->data->CtrlPts, t, 3);//三次ベジエ曲線に限定しているため
            double d = (p - pos).norm();
            if(d < mindist){
                mindist = d; refCreases = cur->data;
            }
        }
    }
    qDebug()<<"clicked object  " << refCreases.get();
}

//アフィン変換の実行(中身はtransform.cppに記載)
void Model::AffinTranse_Crease(int type, const QPointF& befPos, const QPointF& curPos, const QPointF& basePos){
    if(refCreases == nullptr)return;
    switch (type) {
    case 0:
        AffinTrans(refCreases, befPos, curPos);
        break;
    case 1:
        AffinScale(refCreases, basePos, befPos, curPos);
    default:
    case 2:
        AffinRotate(refCreases, basePos, befPos, curPos);
        break;
    }
}

//色の割り当て時に色の0~255を折りの0~piに線形性のある対応にしていない
//ある色の値を閾値として、それよりも大きい値の色を指定すると曲げ角度が急激に変化する(deform関数のColor2Angleを参照)
//閾値の編集
void Model::SetMaxFold(double val){
    ColorPt.angle = val * std::numbers::pi/180.0;
    if(!outline->IsClosed())return;
    for(auto& c: RulingCurve){
        if(c->isempty){
            qDebug() << "crvs has empty";
            return;
        }
    }
    deform();
}

//rulingに割り当てられた色から離散可展面を作成(折り目を含めない)
//色 -> 曲げ角度では閾値をColorPoint変数で定義
//閾値を境界としてrulingに割り当てられた色と曲げ角度は折れ線のような形状となる
void Model::deform(){
    auto Color2Angle = [](double a, ColorPoint CP){
        if(-CP.color <= a && a <= CP.color) return -a * CP.angle/CP.color;
        else if(CP.color < a) return -((std::numbers::pi - CP.angle)/(255.0 - CP.color)* (abs(a) - CP.color) + CP.angle);
        return ((std::numbers::pi - CP.angle)/(255.0 - CP.color)* (abs(a) - CP.color) + CP.angle);
    };
    if(Rulings.empty() || !outline->IsClosed())return;
    std::vector<Eigen::Transform<double, 3, Eigen::Affine>> TMs;
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
   Eigen::Transform<double, 3, Eigen::Affine> M = Eigen::Transform<double, 3, Eigen::Affine>::Identity(), R, T;
   R = Eigen::Matrix3d::Identity();
   T = Eigen::Translation3d(Rulings[0]->o->p);
   TMs.push_back(M * R * T);
   for(int i = 1; i < (int)Rulings.size(); i++){
       Eigen::Vector3d axis = (Rulings[i-1]->v->p - Rulings[i-1]->o->p).normalized();
       M = TMs.back();
       T = Eigen::Translation3d(Rulings[i]->o->p - Rulings[i-1]->o->p);
       R = Eigen::AngleAxisd(Color2Angle(Rulings[i-1]->color, ColorPt), axis);
       Rulings[i]->o->p3_ori = Rulings[i]->o->p3 = M * R * T * Eigen::Vector3d(0,0,0);
       TMs.push_back(M * R * T);
       T = Eigen::Translation3d(Rulings[i]->v->p - Rulings[i-1]->o->p);
       Rulings[i]->v->p3_ori = Rulings[i]->v->p3 = M * R * T * Eigen::Vector3d(0,0,0);
   }


    for(auto&l: outline->Lines){
        for(auto& P: Polygons){
            if(std::find(P.begin(), P.end(), l->v) == P.end())continue;
            for(int i = 0; i < (int)Rulings.size(); i++){
                if(std::find(P.begin(), P.end(), Rulings[i]->v) != P.end()){
                    if(i == (int)Rulings.size() - 1){
                        Eigen::Vector3d axis = (Rulings.back()->v->p - Rulings.back()->o->p).normalized();
                        M = TMs.back();
                        T = Eigen::Translation3d(l->v->p - Rulings.back()->o->p);
                        R = Eigen::AngleAxisd(Color2Angle(Rulings.back()->color, ColorPt), axis);
                        l->v->p3_ori = l->v->p3 = (M * R * T) * Eigen::Vector3d(0,0,0);
                        break;
                    }else{
                        T = Eigen::Translation3d(l->v->p - Rulings[i]->o->p);
                        l->v->p3_ori = l->v->p3 = TMs[i] * T * Eigen::Vector3d(0,0,0);
                        break;
                    }
                }
            }
        }
    }
}

//折り目の追加と不要な折り目(制御点の数が足りずベジエ曲線を作成できないもの)の削除
void Model::RemoveUnable2GenCurve(){
    std::queue<std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>>> q;
    for(auto&btmcrease: NTree_Creases.GetRoot()->children)q.push(btmcrease);
    while (!q.empty()) {
        std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>> cur = q.front(); q.pop();
        for (auto child : cur->children){if (child != nullptr)q.push(child);}
        if(cur->data->CtrlPts.size() <= 3 && cur->data != refCreases)NTree_Creases.Erase(cur->data);
    }
}

//線分とベジエ曲線で表される折り目の交点を計算
std::vector<std::shared_ptr<CrvPt_FL>> getCrossPoint(std::vector<Eigen::Vector3d>& CtrlPts,  const std::shared_ptr<Vertex>& v, const std::shared_ptr<Vertex>& o, int dim){
    std::vector<double>arcT = BezierClipping(CtrlPts, v, o, dim);
    std::vector<std::shared_ptr<CrvPt_FL>> P;
    for(auto&t: arcT){
        if(t < 0 || 1 < t){qDebug()<<"t is not correct value " << t; continue;}
        Eigen::Vector3d v2(0,0,0);
        for (int i = 0; i < int(CtrlPts.size()); i++) v2 += MathTool::BernsteinBasisFunc(dim, i, t) * CtrlPts[i];
        if(!MathTool::is_point_on_line(v2, v->p, o->p))continue;
        double sa = (v2 - o->p).norm(), sc = (v->p - o->p).norm();
        Eigen::Vector3d v3 = sa/sc * (v->p3 - o->p3) + o->p3;
        P.push_back(std::make_shared<CrvPt_FL>(CrvPt_FL(v2, v3, t)));
    }
    return P;
}

void Model::initializeSurfaceVertices(){ for(auto&v: outline->vertices)v->p3 = v->p3_ori;}


//曲線同士の交差を想定しない
void Model::MakeTree(std::shared_ptr<FoldLine>& NewCrease){

    //C：木構造内の折り目, Cnew:追加する折り目
    //-2: 同じ折り目　-1:Cnewが親 0：兄弟, 1:Cnewが子
    auto parent_child_relation = []( const std::shared_ptr<FoldLine>& C, const std::shared_ptr<FoldLine>& Cnew,
                                    const std::vector<std::shared_ptr<Line>>& contour, int btm)->int{
        if(C == Cnew)return -2;
        int i = btm;
        int dim = 3;
        double sumlen = 0.0;
        double dist_l, dist_r, dist_cl, dist_cr;//曲面上の基準とした頂点からの全ての辺を一直線に伸ばしたときの距離
        std::vector<double> C_dist, Cnew_dist;
        std::vector<std::shared_ptr<CrvPt_FL>> C_endpoint, Cnew_endpoint;
        for(const auto&e: contour){
            std::vector<std::shared_ptr<CrvPt_FL>> P = getCrossPoint(C->CtrlPts, e->o, e->v, dim);
            for(const auto&p: P)C_endpoint.push_back(p);
            P = getCrossPoint(Cnew->CtrlPts, e->o, e->v, dim);
            for(const auto&p: P)Cnew_endpoint.push_back(p);
        }
        if(C_endpoint[0]->s > C_endpoint[1]->s)std::swap(C_endpoint[0], C_endpoint[1]);
        if(Cnew_endpoint[0]->s > Cnew_endpoint[1]->s)std::swap(Cnew_endpoint[0], Cnew_endpoint[1]);
        do{
            Eigen::Vector3d x = contour[i]->o->p, x2 = contour[(i + 1) % (int)contour.size()]->o->p;
            if(MathTool::is_point_on_line(C_endpoint[0]->p, x, x2))dist_l = sumlen + (C_endpoint[0]->p - x).norm();
            if(MathTool::is_point_on_line(C_endpoint[1]->p, x, x2))dist_r = sumlen + (C_endpoint[1]->p - x).norm();
            if(MathTool::is_point_on_line(Cnew_endpoint[0]->p, x, x2))dist_cl = sumlen + (Cnew_endpoint[0]->p - x).norm();
            if(MathTool::is_point_on_line(Cnew_endpoint[1]->p, x, x2))dist_cr = sumlen + (Cnew_endpoint[1]->p - x).norm();
            sumlen += (x - x2).norm();
            i = (i + 1) % (int)contour.size();
        }while(i != btm);
        if(dist_l < dist_cl && dist_cl < dist_cr && dist_cr < dist_r)return 1;
        if(dist_cl < dist_l && dist_l < dist_r && dist_r < dist_cr)return -1;
        return 0;
    };

    if((int)NewCrease->CtrlPts.size() <= 3)return;
    NewCrease->FoldingCurve.clear();

    std::queue<std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>>> q;
    q.push(NTree_Creases.GetRoot());
    while (!q.empty()) {
        std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>> cur = q.front(); q.pop();
        for (const auto& child : cur->children){if (child != nullptr)q.push(child);}
        if(cur == NTree_Creases.GetRoot())continue;
        cur->data->FoldingCurve.clear();
    }
    Eigen::Vector3d UpVec(0,-1,0);
    std::vector<std::shared_ptr<FoldLine>> hasFoldingCurve;
    std::shared_ptr<Line> btm = outline->Lines.front();//一番下の辺を探索

    for(auto&l: outline->Lines){
        if(((l->v->p + l->o->p)/2.0).y() > ((btm->v->p + btm->o->p)/2.0).y())btm = l;
    }
    int btm_i = std::distance(outline->Lines.begin(), std::find(outline->Lines.begin(), outline->Lines.end(),btm));
    std::list<std::shared_ptr<FoldLine>> List_Creases;

    class LineOnFL{
    public:
        std::shared_ptr<FoldLine> FL;
        double t;
        std::vector<int> indices;
        LineOnFL(std::shared_ptr<FoldLine> _FL, double _t): FL(_FL), t(_t){}
    };

    //各折曲線上と曲面の交点を導出
    //曲面の各辺と交点があるか探索する
    //曲線が交差した辺のインデックスを持っておく
    RemoveUnable2GenCurve();

    if(NTree_Creases.GetRoot()->children.empty()){
        NTree_Creases.Insert(NTree_Creases.GetRoot(), NewCrease);
    }else{
        auto node = NTree_Creases.GetRoot();
        std::queue<std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>>> q;
        q.push(node);
        while (!q.empty()) {
            std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>> cur = q.front(); q.pop();
            for (auto child : cur->children){if (child != nullptr)q.push(child);}
            if(cur == NTree_Creases.GetRoot() || (int)cur->data->CtrlPts.size() <= 3)continue;
            int res = parent_child_relation(cur->data, NewCrease, outline->Lines, btm_i);
            if(res == 1){
                if(cur->children.empty()){
                    NTree_Creases.Insert(cur->data, NewCrease);
                    return;
                }

            }else{
                auto parent = NTree_Creases.GetParent(cur->data);
                NTree_Creases.Insert(parent->data, NewCrease);
                std::vector<int> res_child;
                for(auto& child: parent->children)  res_child.push_back(parent_child_relation(child->data, NewCrease, outline->Lines, btm_i));
                for(int j = 0; j < (int)res_child.size(); j++){
                    if(res_child[j] == -1)NTree_Creases.Planting(NewCrease, parent->children[j]->data);
                }
            }
        }
    }
    NTree_Creases.print();
    qDebug() <<"NewCrease = " << NewCrease.get() << "   ,  refCrease = " << refCreases.get();
    return;
}

template <typename T>
bool IsOverlapped(std::shared_ptr<T>&p, std::shared_ptr<T>&q){return ((p->p - q->p).norm() < 1e-8)? true: false;};

//端点の空間座標の計算
void Model::SetEndPoint(std::shared_ptr<Vertex4d>&v4d, const std::vector<std::shared_ptr<Line>>& Surface, const std::vector<std::shared_ptr<Line>>& Rulings, bool IsupdateEndPt){
    int s = Surface.size();
    std::shared_ptr<Vertex> prev = nullptr, next = nullptr;
    std::shared_ptr<Vertex> downcast = std::dynamic_pointer_cast<Vertex>(v4d->first);
    for(int i = 0; i < s; i++){
        if(IsOverlapped(downcast, Surface[i]->o)){v4d->first->p3 = Surface[i]->o->p3; return;}
        if(MathTool::is_point_on_line(v4d->first->p, Surface[i]->o->p, Surface[i]->v->p)){prev = Surface[i]->o; next = Surface[i]->v;}
    }

    std::queue<std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>>> q;
    for(auto&btmcrease: NTree_Creases.GetRoot()->children)q.push(btmcrease);
    while (!q.empty()) {
        std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>> cur = q.front(); q.pop();
        for (auto child : cur->children){if (child != nullptr)q.push(child);}
        if(cur->data->FoldingCurve.empty())continue;
        for(auto it = cur->data->FoldingCurve.begin() + 1; it != cur->data->FoldingCurve.end()-1;it++){
            if(!(*it)->IsCalc)continue;
            if(IsOverlapped(downcast, (*it)->second)){v4d->first->p3 = (*it)->second->p3; return;}
            if(IsOverlapped(downcast, (*it)->third)){v4d->first->p3 = (*it)->third->p3; return;}
            if(MathTool::is_point_on_line(v4d->first->p, (*it)->second->p, next->p) && MathTool::is_point_on_line((*it)->second->p, prev->p, next->p)){prev = (*it)->second;continue;}
            if(MathTool::is_point_on_line(v4d->first->p, (*it)->second->p, prev->p) && MathTool::is_point_on_line((*it)->second->p, prev->p, next->p)){next = (*it)->second;continue;}
            if(MathTool::is_point_on_line(v4d->first->p, (*it)->third->p, next->p) && MathTool::is_point_on_line((*it)->third->p, prev->p, next->p)){prev = (*it)->third;continue;}
            if(MathTool::is_point_on_line(v4d->first->p, (*it)->third->p, prev->p) && MathTool::is_point_on_line((*it)->third->p, prev->p, next->p)){next = (*it)->third;continue;}

            std::shared_ptr<Vertex> downcast_f = std::dynamic_pointer_cast<Vertex>((*it)->first);
            //if((*it) == v4d->first)continue;
            if(IsOverlapped(downcast, downcast_f)){v4d->first->p3 = downcast_f->p3; return;}
            if(MathTool::is_point_on_line(v4d->first->p, (*it)->first->p, next->p) && MathTool::is_point_on_line((*it)->first->p, prev->p, next->p)){prev = (*it)->first;continue;}
            if(MathTool::is_point_on_line(v4d->first->p, (*it)->first->p, prev->p) && MathTool::is_point_on_line((*it)->first->p, prev->p, next->p)){next = (*it)->first;continue;}
        }
    }

    for(auto&r: Rulings){
        if(r->hasCrossPoint)continue;
        if(IsOverlapped(downcast, r->v)){v4d->first->p3 = r->v->p3; return;}
        if(IsOverlapped(downcast, r->o)){v4d->first->p3 = r->o->p3; return;}
        if(MathTool::is_point_on_line(v4d->first->p, r->v->p, prev->p) && MathTool::is_point_on_line(r->v->p, prev->p, next->p))next = r->v;
        if(MathTool::is_point_on_line(v4d->first->p, r->v->p, next->p) && MathTool::is_point_on_line(r->v->p, prev->p, next->p))prev = r->v;
        if(MathTool::is_point_on_line(v4d->first->p, r->o->p, prev->p) && MathTool::is_point_on_line(r->o->p, prev->p, next->p))next = r->o;
        if(MathTool::is_point_on_line(v4d->first->p, r->o->p, next->p) && MathTool::is_point_on_line(r->o->p, prev->p, next->p))prev = r->o;
    }

    auto EdgeSplit = [](std::vector<std::shared_ptr<Line>>& Lines, std::shared_ptr<Vertex>& v){
        for(auto&l: Lines){
            if(MathTool::is_point_on_line(v->p, l->o->p, l->v->p) && (l->o->p - v->p).norm() > 1e-5 && (l->v->p - v->p).norm() > 1e-5){
                std::shared_ptr<Line> l_split = std::make_shared<Line>(l->o, v, l->et);
                l->o = v;
                Lines.push_back(l_split);
                return;
            }
        }
    };
    std::vector<std::shared_ptr<Line>> polyline_surface;
    for(auto&l: Surface) polyline_surface.push_back(std::make_shared<Line>(l->o, l->v, l->et));
    for(auto&r: Rulings){
        if(r->hasCrossPoint)continue;
        EdgeSplit(polyline_surface, r->v);
        EdgeSplit(polyline_surface, r->o);
    }

    std::queue<std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>>> empty;
    std::swap(q, empty);
    for(auto&btmcrease: NTree_Creases.GetRoot()->children)q.push(btmcrease);
    while (!q.empty()) {
        std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>> cur = q.front(); q.pop();
        for (auto child : cur->children){if (child != nullptr)q.push(child);}
        if(cur->data->FoldingCurve.empty())continue;
        for(auto it = cur->data->FoldingCurve.begin(); it != cur->data->FoldingCurve.end();it++){
            if(!(*it)->IsCalc || (*it) == v4d)continue;
            std::shared_ptr<Vertex> dowscastV = std::dynamic_pointer_cast<Vertex>((*it)->first);
            EdgeSplit(polyline_surface, dowscastV);
            if(it == cur->data->FoldingCurve.begin() || it == cur->data->FoldingCurve.end()-1)continue;
            EdgeSplit(polyline_surface, (*it)->second);
            EdgeSplit(polyline_surface, (*it)->third);
        }
    }

    for(auto&l: polyline_surface){
        if(!IsOverlapped(l->o, downcast) && !IsOverlapped(l->v, downcast) && MathTool::is_point_on_line(v4d->first->p, l->v->p, l->o->p)){
            if(IsupdateEndPt){
                v4d->first->p3_ori = v4d->first->p3 = (l->v->p - v4d->first->p).norm()*(l->o->p3 - l->v->p3).normalized() + l->v->p3;
            }
            v4d->second = l->o; v4d->third = l->v;

            std::swap(q, empty);
            for(auto&btmcrease: NTree_Creases.GetRoot()->children)q.push(btmcrease);
            while (!q.empty()) {
                std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>> cur = q.front(); q.pop();
                for (auto child : cur->children){if (child != nullptr)q.push(child);}
                if(cur->data->FoldingCurve.empty())continue;
                cur->data->AlignmentVertex4dDirection();
            }
            return;
        }
    }
}

//曲面の輪郭がもつ三次元上の頂点を正しい位置に配置する
//端点の修正は曲面の頂点修正後に行う
void Model::SetOnVertices_outline(){
    class vertexinfo{
    public:
        std::shared_ptr<Vertex> v;
        int vtype;//0:頂点,1:first,2:second, 3:third, 4:ruling
        std::shared_ptr<Vertex> basePt;
        vertexinfo(const std::shared_ptr<Vertex>&_v, int t, const std::shared_ptr<Vertex>&b){v = _v; vtype=t; basePt = b;}
    };
    std::vector<vertexinfo> VerticesInfo;
    for(auto&v: outline->vertices)VerticesInfo.push_back(vertexinfo(v, 0, nullptr));
    initializeSurfaceVertices();
    auto EdgeSplit = [&VerticesInfo](const std::shared_ptr<Vertex>& v, int vtype, const std::shared_ptr<Vertex>&basePt){
        int s = VerticesInfo.size();
        for(int i = 0; i < s; i++){
            //if((v->p - VerticesInfo[i].v->p).norm() < 1e-10)return;
            if(MathTool::is_point_on_line(v->p, VerticesInfo[i].v->p, VerticesInfo[(i+1) % s].v->p)){
                VerticesInfo.insert(VerticesInfo.begin() + (i+1) % s, vertexinfo(v, vtype, basePt));
                return;
            }
        }
    };

    std::queue<std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>>> q;
    for(auto&btmcrease: NTree_Creases.GetRoot()->children)q.push(btmcrease);
    while (!q.empty()) {
        std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>> cur = q.front(); q.pop();
        for (auto child : cur->children){if (child != nullptr)q.push(child);}
        if((int)cur->data->FoldingCurve.size() <= 2)continue;
        std::vector<std::shared_ptr<Vertex4d>> ValidFC;
        for(auto&fc: cur->data->FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
        EdgeSplit(ValidFC.front()->first, 1, ValidFC[1]->first);
        EdgeSplit(ValidFC.back()->first, 1, ValidFC.end()[-2]->first);
        for(auto itr = ValidFC.begin()+1; itr != ValidFC.end()-1; itr++){
            if(!(*itr)->IsCalc)continue;
            EdgeSplit((*itr)->second, 2, (*itr)->first);
            EdgeSplit((*itr)->third, 3, (*itr)->first);
        }
    }

    for(auto&r: Rulings){
        if(r->hasCrossPoint)continue;
        EdgeSplit(r->v, 4, r->o);
        EdgeSplit(r->o, 4, r->v);
    }
    int visize = VerticesInfo.size();
    for(auto&p: outline->vertices){
        auto itr = std::find_if(VerticesInfo.begin(), VerticesInfo.end(), [&p](const vertexinfo& v){return v.v == p;});
        int i = std::distance(VerticesInfo.begin(), itr);
        int i_prev = (i - 1 + visize) % visize, i_next = (i + 1) % visize;
        std::shared_ptr<Vertex> prev, next, basePt;
        //thirdの方を基準として計算する
        if(VerticesInfo[i_prev].vtype == 2 && VerticesInfo[i_next].vtype == 3){
            next =  VerticesInfo[i_next].v; basePt =  VerticesInfo[i_next].basePt;
        }
        else if(VerticesInfo[i_prev].vtype == 3 && VerticesInfo[i_next].vtype == 2){

        }else{
            auto itr_overlapped = std::find_if(VerticesInfo.begin(), VerticesInfo.end(), [&](const vertexinfo& v){return (p != v.v && (p->p - v.v->p).norm() < 1e-9);});
            if(itr_overlapped != VerticesInfo.end()){
                p->p3 = (*itr_overlapped).v->p3;
                continue;
            }
            prev = VerticesInfo[i_prev].v; next = VerticesInfo[i_next].v; basePt = (VerticesInfo[i_prev].basePt != nullptr)? VerticesInfo[i_prev].basePt: VerticesInfo[i_next].basePt;
            if(basePt == nullptr)continue;
            Eigen::Vector3d po = p->p - basePt->p, a = prev->p - basePt->p, b = next->p - basePt->p;
            double s,t;
            if(a.norm() < 1e-15 || b.norm() < 1e-15){continue; }
            if(abs(a.y()) <  1e-12){
                s = po.y()/b.y();
                t = (po.x() - s*b.x())/a.x();
            }else if(abs(b.y()) < 1e-12){
                t = po.y()/a.y();
                s = (po.x() - t*a.x())/b.x();
            }else{
                s = (a.y() * po.x() - a.x()*po.y())/(a.y()*b.x() - a.x()*b.y());
                t = (po.y() - s*b.y())/a.y();
            }
            p->p3 = t * (prev->p3 - basePt->p3) + s * (next->p3 - basePt->p3) + basePt->p3;
        }
    }

    std::queue<std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>>> empty;
    std::swap(q, empty);
    for(auto&btmcrease: NTree_Creases.GetRoot()->children)q.push(btmcrease);
    while (!q.empty()) {
        std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>> cur = q.front(); q.pop();
        for (auto child : cur->children){if (child != nullptr)q.push(child);}
        if(cur->data->FoldingCurve.empty())continue;
        cur->data->SortCurve();
    }
    return;

    int s = outline->vertices.size();
    initializeSurfaceVertices();
    auto Creases = NTree_Creases.NTree2Array();
    for(auto&fl:Creases)fl->data->SortCurve();
    for(int i = 0; i < s; i++){
        std::shared_ptr<Vertex> p = outline->vertices[i];
        bool Isvertexlapped = false;
        std::shared_ptr<Vertex> prev = outline->vertices[(i - 1 + s)  % s], next = outline->vertices[(i + 1) % s];//前後の頂点で作られるlineのうち最も近いものを入れる
        std::shared_ptr<Vertex> basePt = nullptr;
        for(auto&fl: Creases){
            if((int)fl->data->FoldingCurve.size() <= 2)continue;
            for(auto&v4d: fl->data->FoldingCurve){
                if(!v4d->IsCalc)continue;
                std::shared_ptr<Vertex> downcast = std::dynamic_pointer_cast<Vertex>(v4d->first);
                if(IsOverlapped(downcast, outline->vertices[i])){
                    Isvertexlapped = true; basePt = v4d->first; break;
                }
                if(IsOverlapped(v4d->second, outline->vertices[i])){
                    if(fl->data->FoldingCurve.front() != v4d && fl->data->FoldingCurve.back() != v4d){
                        Isvertexlapped = true; basePt = v4d->second; break;
                    }

                }
                if(IsOverlapped(v4d->third, outline->vertices[i])){
                    if(fl->data->FoldingCurve.front() != v4d && fl->data->FoldingCurve.back() != v4d){
                        Isvertexlapped = true; basePt = v4d->third; break;
                    }
                }
                if(MathTool::is_point_on_line(v4d->first->p, prev->p, outline->vertices[i]->p)&& !IsOverlapped(downcast, prev) && !IsOverlapped(downcast, outline->vertices[i])){
                    prev = v4d->first;
                }
                if(MathTool::is_point_on_line(v4d->first->p, next->p, outline->vertices[i]->p)&& !IsOverlapped(downcast, next) && !IsOverlapped(downcast, outline->vertices[i])){
                    next = v4d->first;
                }
                if(MathTool::is_point_on_line(v4d->second->p, prev->p, outline->vertices[i]->p) && !IsOverlapped(v4d->second, prev) && !IsOverlapped(v4d->second, outline->vertices[i])){
                    prev = v4d->second;
                    if(!IsOverlapped(downcast, prev) && !IsOverlapped(downcast, next))basePt = v4d->first;
                }
                if(MathTool::is_point_on_line(v4d->second->p, next->p, outline->vertices[i]->p) && !IsOverlapped(v4d->second, next) && !IsOverlapped(v4d->second, outline->vertices[i])){
                    next = v4d->second;
                    if(!IsOverlapped(downcast, prev) && !IsOverlapped(downcast, next))basePt = v4d->first;
                }
                if(MathTool::is_point_on_line(v4d->third->p, prev->p, outline->vertices[i]->p) && !IsOverlapped(v4d->third, prev) && !IsOverlapped(v4d->third, outline->vertices[i])){
                    prev = v4d->third;
                    if(!IsOverlapped(downcast, prev) && !IsOverlapped(downcast, next))basePt = v4d->first;
                }
                if(MathTool::is_point_on_line(v4d->third->p, next->p, outline->vertices[i]->p) && !IsOverlapped(v4d->third, next) && !IsOverlapped(v4d->third, outline->vertices[i])){
                    next = v4d->third;
                    if(!IsOverlapped(downcast, prev) && !IsOverlapped(downcast, next))basePt = v4d->first;
                }
            }
        }
        for(auto&r: Rulings){
            if(r->hasCrossPoint)continue;
            if(Isvertexlapped)break;
            if(IsOverlapped(r->v, outline->vertices[i])){Isvertexlapped = true; basePt = r->v; break;}
            if(IsOverlapped(r->o, outline->vertices[i])){Isvertexlapped = true; basePt = r->o; break;}
            if(MathTool::is_point_on_line(r->v->p, prev->p, outline->vertices[i]->p) && MathTool::is_point_on_line(r->v->p, outline->vertices[i]->p, prev->p)){
                prev = r->v;
                basePt = r->o;
            }
            if(MathTool::is_point_on_line(r->o->p, prev->p, outline->vertices[i]->p) && MathTool::is_point_on_line(r->o->p, outline->vertices[i]->p, prev->p)){
                prev = r->o;
                basePt = r->v;
            }
            if(MathTool::is_point_on_line(r->v->p, next->p, outline->vertices[i]->p) && MathTool::is_point_on_line(r->v->p, next->p, outline->vertices[i]->p)){
                next = r->v;
                basePt = r->o;
            }
            if(MathTool::is_point_on_line(r->o->p, next->p, outline->vertices[i]->p) && MathTool::is_point_on_line(r->o->p, next->p, outline->vertices[i]->p)){
                next = r->o;
                basePt = r->v;
            }
        }
        if(Isvertexlapped){
            outline->vertices[i]->p3 = basePt->p3;
        }else{
            auto itr = std::find_if(Rulings.begin(), Rulings.end(), [&](const std::shared_ptr<Line>& r){return (prev->p == r->o->p && next->p == r->v->p) || (prev->p == r->v->p && next->p == r->o->p);});
            if(itr != Rulings.end() || basePt == nullptr){}
            else{
                Eigen::Vector3d po = outline->vertices[i]->p - basePt->p, a = prev->p - basePt->p, b = next->p - basePt->p;
                double s,t;
                if(a.norm() < 1e-15 || b.norm() < 1e-15) continue;

                if(abs(a.y()) <  1e-12){
                    s = po.y()/b.y();
                    t = (po.x() - s*b.x())/a.x();
                }else if(abs(b.y()) < 1e-12){
                    t = po.y()/a.y();
                    s = (po.x() - t*a.x())/b.x();
                }else{
                    s = (a.y() * po.x() - a.x()*po.y())/(a.y()*b.x() - a.x()*b.y());
                    t = (po.y() - s*b.y())/a.y();
                }
                outline->vertices[i]->p3 = t * (prev->p3 - basePt->p3) + s * (next->p3 - basePt->p3) + basePt->p3;
            }
        }
    }
    for(auto&fl: Creases){
        if(fl->data->FoldingCurve.empty())continue;
        fl->data->SortCurve();
    }
}

//曲線折りの実行
bool Model::BendingModel(double wp, double wsim, int dim, int alg){
    AddNewCrease(refCreases);
    if(NTree_Creases.GetRoot()->children.empty())return false;
    std::vector<std::shared_ptr<Vertex>> Poly_V = outline->getVertices();
    std::queue<std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>>> q; 
    for(auto&btmcrease: NTree_Creases.GetRoot()->children)q.push(btmcrease);
    int rank = 0;
    while(!q.empty()){
        auto cur = q.front(); q.pop();
        if(alg == -1){
        auto itr = std::find(NTree_Creases.GetRoot()->children.begin(), NTree_Creases.GetRoot()->children.end(), cur);
        if(itr != NTree_Creases.GetRoot()->children.end())cur->data->initialize_foldstate(Poly_V);
        }else{
            cur->data->ReassignColor();
            cur->data->Optimization_FlapAngle(Poly_V, wp, wsim, rank);
            cur->data->PropagateOptimization_Vertex(Poly_V, wp, wsim);
            cur->data->applyAAAMethod(Poly_V, cur->data->a_flap);
        }
        for (const auto& child : cur->children){
            if(child != nullptr){
                child->data->reassignruling(cur->data, outline->Lines, Rulings);
                rank++;
                q.push(child);
            }
        }     
    }
    SetOnVertices_outline();
    return true;
}

//曲面の輪郭上にのっている頂点の情報をかえす
std::vector<vertexinfo> Model::MappingVertex(bool IsRemovingOverlapping){
    std::vector<vertexinfo> VerticesInfo;
    for(auto&v: outline->vertices)VerticesInfo.push_back(vertexinfo(v, 0));
    auto EdgeSplit = [&](const std::shared_ptr<Vertex>& v, int vtype){
        int s = VerticesInfo.size();
        for(int i = 0; i < s; i++){
            if((v->p - VerticesInfo[i].v->p).norm() < 1e-10){//頂点が重なる場合
                if(IsRemovingOverlapping)return;
                VerticesInfo[i].v = v; VerticesInfo[i].vtype = vtype;
                return;
            }
            if(MathTool::is_point_on_line(v->p, VerticesInfo[i].v->p, VerticesInfo[(i+1) % s].v->p)){
                if((v->p - VerticesInfo[(i+1)%s].v->p).norm() < 1e-10 && IsRemovingOverlapping)return;
                VerticesInfo.insert(VerticesInfo.begin() + (i+1) % s, vertexinfo(v, vtype));
                return;
            }
        }
    };
    auto Creases = NTree_Creases.NTree2Array();
    for(auto&fl:Creases){
        if((int)fl->data->FoldingCurve.size() <= 2)continue;
        EdgeSplit(fl->data->FoldingCurve.front()->first, 1);EdgeSplit(fl->data->FoldingCurve.back()->first, 1);

        for(auto itr = fl->data->FoldingCurve.begin()+1; itr != fl->data->FoldingCurve.end()-1; itr++){
            EdgeSplit((*itr)->second, 2); EdgeSplit((*itr)->third, 3);
        }
    }
    for(auto&r: Rulings){
        if(r->hasCrossPoint)continue;
        EdgeSplit(r->v, 4); EdgeSplit(r->o, 4);
    }
    return VerticesInfo;
}

//折り目を追加したときの木構造の再構成
bool Model::AddNewCrease(std::shared_ptr<FoldLine>& NewCrease){
    if(NewCrease == nullptr)return false;
    MakeTree(NewCrease);
    int dim = 3;

    std::vector<std::shared_ptr<Vertex>> Poly_V = outline->vertices;
    std::queue<std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>>> q;
    for(auto&btmcrease: NTree_Creases.GetRoot()->children)q.push(btmcrease);
    while (!q.empty()) {
        std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>> cur = q.front(); q.pop();
        for (const auto& child : cur->children){if (child != nullptr)q.push(child);}
        auto parentNode = NTree_Creases.GetParent(cur);
        if(parentNode == NTree_Creases.GetRoot()) SplitRulings(dim);
        else cur->data->reassignruling(parentNode->data, outline->Lines, Rulings);
    }
    return true;
}

//rulingと折り目の交点を計算し三次元区間上の折り目の頂点座標の計算
bool Model::SplitRulings(int dim){
    Eigen::Vector3d UpVec(0,-1,0);
    if(NTree_Creases.GetRoot()->children.empty())return false;

    for(auto&v: outline->Lines){
        v->o->p3 = v->o->p3_ori;
        v->v->p3 = v->v->p3_ori;
    }
    for(const auto& crease: NTree_Creases.GetRoot()->children){
        if((int)crease->data->CtrlPts.size() <= dim)continue;
        crease->data->FoldingCurve.clear();
        for(auto& r: Rulings){
            std::vector<std::shared_ptr<CrvPt_FL>> P = getCrossPoint(crease->data->CtrlPts, r->v, r->o, dim);
            if(!P.empty()){
                for(auto&p: P){
                    r->hasCrossPoint = true;
                    std::shared_ptr<Vertex> sec, thi;
                    if(UpVec.dot((r->v->p - r->o->p).normalized()) > 0){ sec = r->v->deepCopy(); thi = r->o->deepCopy();}
                    else{sec = r->o->deepCopy(); thi = r->v->deepCopy(); }
                    std::shared_ptr<Vertex4d> v4d = std::make_shared<Vertex4d>(p, sec, thi); v4d->addline(thi, sec);
                    crease->data->FoldingCurve.push_back(v4d);
                }
            }
        }
        for(auto& line: outline->Lines){
            std::vector<std::shared_ptr<CrvPt_FL>> P = getCrossPoint(crease->data->CtrlPts, line->v, line->o, dim);
            if(!P.empty()){
                for(auto&p: P){
                    std::shared_ptr<Vertex> sec, thi;
                    if(UpVec.dot((line->v->p - line->o->p).normalized()) > 0){ sec = line->v->deepCopy(); thi = line->o->deepCopy();}
                    else{sec = line->o->deepCopy(); thi = line->v->deepCopy(); }
                    std::shared_ptr<Vertex4d> v4d = std::make_shared<Vertex4d>(p, sec, thi); v4d->addline(thi, sec);
                    crease->data->FoldingCurve.push_back(v4d);
                }
            }
        }
        crease->data->SortCurve();
        crease->data->AlignmentVertex4dDirection();
        crease->data->validsize = crease->data->FoldingCurve.size();
    }
    return true;
}

//折り目を最小二乗平面上に乗せる
void Model::flatten_lsp(std::shared_ptr<FoldLine>& FldLine){
    if(!NTree_Creases.find(FldLine))return;
    _flatten_lsp(FldLine);
}

//折り目とrulingの交点を間引く
void Model::SimplifyModel(int iselim){
    if(refCreases == nullptr || refCreases->FoldingCurve.empty())return;
    bool isroot = false;
    for(const auto& node: NTree_Creases.GetRoot()->children){
        if(node->data == refCreases)isroot = true;
    }

    int validsize = refCreases->validsize - iselim;
    refCreases->SimplifyModel(validsize, isroot);
}

//間引いた後にrulingの向きを修正して滑らかさの補間(可展性を考慮していない)
bool Model::Smoothing(){
    auto root = NTree_Creases.GetRoot();
    if(root == nullptr)return false;
    std::queue<std::shared_ptr<TreeNode<std::shared_ptr<FoldLine>>>>q;
    q.push(root);
    auto Poly_V = outline->getVertices();
    while (!q.empty()) {
        auto cur = q.front(); q.pop();
        cur->data->SimpleSmooothSrf(Poly_V);
        for (const auto& child : cur->children){
            if(child != nullptr){
                child->data->reassignruling(cur->data, outline->Lines, Rulings);
                q.push(child);
            }
        }
    }
    return true;
}

void Model::modifyFoldingCurvePositionOn3d(){
    auto Creases = NTree_Creases.NTree2Array();
    for(auto&crease: Creases){
        for(auto&fldCrv: crease->data->FoldingCurve){
            for(auto&r: Rulings){
                if(!MathTool::is_point_on_line(fldCrv->first->p, r->o->p, r->v->p))continue;
                double t = (fldCrv->first->p - r->o->p).norm()/(r->o->p - r->v->p).norm();
                fldCrv->first->p3_ori = fldCrv->first->p3 = t * (r->v->p3 - r->o->p3) + r->o->p3;
                break;
            }
            for(auto&l: outline->Lines){
                if(!MathTool::is_point_on_line(fldCrv->first->p, l->o->p, l->v->p))continue;
                double t = (fldCrv->first->p - l->o->p).norm()/(l->o->p - l->v->p).norm();
                fldCrv->first->p3_ori = fldCrv->first->p3 = t * (l->v->p3 - l->o->p3) + l->o->p3;
                break;
            }
        }
    }
}

//曲面の輪郭が凹型でないと仮定
void Model::drawOutline(QPointF& cursol, int drawtype, double gridsize, const QSize& S, bool IsClicked){
    QPointF _p = SetOnGrid(cursol, gridsize, S);
    Eigen::Vector3d p(_p.x(), _p.y(), 0);
    if(drawtype == 0 || drawtype == 1){
        if(outline->hasPtNum != 2)outline->addVertex(p);
    }else if(drawtype == 2) outline->drawPolygon(p, IsClicked);

    //頂点の順番を時計回りに修正する
    if(outline->IsClosed()){
        Eigen::Vector3d center(0,0,0);
        for(auto&p: outline->vertices)center += p->p;
        center /= static_cast<double>(outline->vertices.size());
        Eigen::Vector3d N = (outline->vertices[0]->p - center).cross(outline->vertices[1]->p - center);
        if(N.dot(Eigen::Vector3d(0,0,1)) < 0){
            std::reverse(outline->vertices.begin(), outline->vertices.end());
            for(auto&l: outline->Lines) std::swap(l->o, l->v);
            std::reverse(outline->Lines.begin(), outline->Lines.end());
        }
    }
}

void Model::editOutlineVertex(QPointF& cursol, double gridsize, const QSize& S, int event){
    QPointF _p = SetOnGrid(cursol, gridsize, S);
    Eigen::Vector3d p(_p.x(), _p.y(), 0);
    if(event == 0){
        double dist = 5;
        movingVertexIndex = -1;//0~: vertex
        for(int i = 0; i < (int)outline->vertices.size(); i++){
            if((outline->vertices[i]->p - p).norm() < dist){
                movingVertexIndex = i; dist = (outline->vertices[i]->p - p).norm();
            }
        }
    }else if(event == 1){
        outline->MoveVertex(p, movingVertexIndex);
        if(outline->IsClosed())deform();
    }else if(event == 2)movingVertexIndex = -1;
}

void Model::LinearInterPolation(const std::vector<std::shared_ptr<Line>>& path){
    if(path.size() < 2)return;

    Eigen::Vector3d befcenter = Eigen::Vector3d(-1,-1,-1), center;
    double len = 0.0;
    for(auto& l: path){
        if(befcenter == Eigen::Vector3d(-1,-1,-1)){
            befcenter = (l->v->p + l->o->p)/2.0;
            continue;
        }     
        center = (l->v->p + l->o->p)/2.0;
        len += (center - befcenter).norm();
        befcenter = center;
    }
    double r = (GradationPoints[1]->color - GradationPoints[0]->color)/len;
    befcenter = Eigen::Vector3d(-1,-1,-1);
    std::shared_ptr<Line> bef = std::shared_ptr<Line>(nullptr);
    for(auto&l : path){
        if(bef == nullptr){
            bef = l;
            befcenter = (l->v->p + l->o->p)/2.0;
            continue;
        }
        center = (l->v->p + l->o->p)/2.0;
        l->color = r * (center -  befcenter).norm() + bef->color;
        bef = l;
        befcenter = center;
    }
}

void Model::SetGradationValue(int val, const std::shared_ptr<Line>& refL){
    if(Rulings.empty() || std::find(Rulings.begin(), Rulings.end(), refL) == Rulings.end())return;
    if(refL->et != EdgeType::r)return;
    refL->color += val;
    refL->color = (refL->color < -255.0)? -255.0 : (255.0 < refL->color)? 255.0 : refL->color;

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
    LinearInterPolation(path);
    return;
}


void Model::addRulings(){
    if(!outline->IsClosed())return;
    Rulings.clear();
    for(auto&c: RulingCurve)CrossDetection(outline, c);
    CrossDection4AllCurve();
    for(auto&c: RulingCurve){
        for(auto&r: c->Rulings){
            if(r->IsCrossed == -1 || r->IsCrossed == 0) Rulings.push_back(r);
        }
    }
}

void Model::SelectCurve(QPointF pt){
    int ptInd;
    auto crvInd = searchPointIndex(pt, ptInd, 1);
    std::fill(refCrv.begin(), refCrv.end(), 0);
    if(crvInd[0] != -1)refCrv[crvInd[0]] = 1;
}

int Model::IsSelectedCurve(){
    auto itr = std::find(refCrv.begin(), refCrv.end(), 1);
    if(itr == refCrv.end())return -1;
    return std::distance(refCrv.begin(), itr);
}

std::array<int, 2> Model::getSelectedCurveIndex(QPointF pt){
    int ptInd;
    return searchPointIndex(pt, ptInd, 1);
}

//type = 0 -> Control Point, 1: Curve Point
std::array<int, 2> Model::searchPointIndex(QPointF _p, int& ptInd, int type){
    Eigen::Vector3d p(_p.x(), _p.y(), 0);
    double dist = 5.0;
    ptInd = -1;
    std::array<int, 2> CrvInds{-1,-1};
    if(type == 0){
        for(int j = 0; j < (int)RulingCurve.size(); j++){
            for(int i = 0; i < (int)RulingCurve[j]->ControllPoints.size(); i++){
                Eigen::Vector3d cp = RulingCurve[j]->ControllPoints[i];
                if(dist  > (cp - p).norm()){
                    dist = (cp - p).norm();
                    ptInd = i;
                    CrvInds[0] = j;
                }
            }
        }

    }else{
        for(int j = 0; j < (int)RulingCurve.size(); j++){
            for(int i = 0; i < (int)RulingCurve[j]->CurvePoints.size(); i++){
                Eigen::Vector3d cp = RulingCurve[j]->CurvePoints[i];
                if(dist > (cp - p).norm()){
                    dist = (cp - p).norm();
                    ptInd = i;
                    CrvInds[0] = j;
                }
            }
        }
    }

    for(const auto&crease: NTree_Creases.NTree2Array()){
        if(crease->data == nullptr)continue;
        for(int j = 0; j < (int)crease->data->CtrlPts.size(); j++){
            double d = (crease->data->CtrlPts[j] - p).norm();
            if(dist > d){
                dist = (crease->data->CtrlPts[j] - p).norm();
                ptInd = j; refCreases = crease->data;
            }
        }
    }

    return CrvInds;
}

int Model::DeleteCurve(){
    int DelIndex = IsSelectedCurve();
    if(DelIndex == -1)return -1;
    RulingCurve.erase(RulingCurve.begin() + DelIndex);
    refCrv.erase(refCrv.begin() + DelIndex);
    RulingCurve.shrink_to_fit();
    refCrv.shrink_to_fit();
    return DelIndex;
}

void Model::Check4Param(int curveDimention, std::vector<int>& deleteIndex){
    deleteIndex.clear();
    for(int i = (int)RulingCurve.size() - 1; i >= 0; i--){
        if((RulingCurve[i]->getCurveType() == CurveType::bezier3 && (int)RulingCurve[i]->ControllPoints.size() < curveDimention) ||
                (RulingCurve[i]->getCurveType() == CurveType::bsp3 && (int)RulingCurve[i]->ControllPoints.size() <= curveDimention) ||
                (RulingCurve[i]->getCurveType() == CurveType::line && (int)RulingCurve[i]->ControllPoints.size() < 2) ||
                (RulingCurve[i]->getCurveType() == CurveType::arc && (int)RulingCurve[i]->ControllPoints.size() < 3)) {
            deleteIndex.push_back(i);
            RulingCurve.erase(RulingCurve.begin() + i);
            refCrv.erase(refCrv.begin() + i);
        }

    }
    RulingCurve.shrink_to_fit();
    refCrv.shrink_to_fit();
    GradationPoints.clear();
}

int Model::AddNewCurve(CurveType curveType, int DivSize){
    std::shared_ptr<CRV> NewCurve = std::make_shared<CRV>(crvPtNum, DivSize);
    NewCurve->setCurveType(curveType);
    RulingCurve.insert(RulingCurve.begin(), NewCurve);
    refCrv.insert(refCrv.begin(), 0);
    std::fill(refCrv.begin(), refCrv.end(), 0);
    refCrv[0] = 1;
    return 0;
}

bool Model::ControlPoint(PaintTool dtype, int event, QPointF p, int curveDimention, int DivSize, std::shared_ptr<FoldLine>& newCrease){
    Eigen::Vector3d _p(p.x(), p.y(), 0);
    if(dtype == PaintTool::Crease){
        if(refCreases == nullptr)return false;
        if(event == 0){
            refCreases->addCtrlPt(_p, curveDimention);
            if((int)refCreases->CtrlPts.size() <= 3)return false;
            else return true;

        }else{
            return refCreases->delCtrlPt(_p, curveDimention, outline);
        }
    }else{
        if(event == 0){
            int AddPtIndex = IsSelectedCurve();
            if(AddPtIndex == -1)return false;
            bool res = false;
            if(RulingCurve[AddPtIndex]->getCurveType() == CurveType::bsp3){
                double dist = 5.0;
                int n = RulingCurve[AddPtIndex]->movePtIndex(_p, dist);
                if(n == -1){
                    RulingCurve[AddPtIndex]->ControllPoints.push_back(_p);
                    res = RulingCurve[AddPtIndex]->drawBspline(curveDimention, crvPtNum);
                }
            }
            else if(RulingCurve[AddPtIndex]->getCurveType() == CurveType::line){
                if(RulingCurve[AddPtIndex]->ControllPoints.size() < 2)RulingCurve[AddPtIndex]->ControllPoints.push_back(_p);
                if(RulingCurve[AddPtIndex]->ControllPoints.size() == 2)res = RulingCurve[AddPtIndex]->drawLine();
            }
            else if(RulingCurve[AddPtIndex]->getCurveType() == CurveType::arc){
                if(RulingCurve[AddPtIndex]->ControllPoints.size() == 0){
                    if(outline->IsClosed()){
                        bool PointInFace = outline->IsPointInFace(_p);
                        bool PointOnLines = false;
                        for(auto&l : outline->Lines){
                            if(l->is_on_line(_p))PointOnLines = true;
                        }
                        if(!PointInFace || PointOnLines)RulingCurve[AddPtIndex]->ControllPoints.push_back(_p);
                    }
                    else RulingCurve[AddPtIndex]->ControllPoints.push_back(_p);
                }
                else if(RulingCurve[AddPtIndex]->ControllPoints.size() == 1)RulingCurve[AddPtIndex]->ControllPoints.push_back(_p);
                else if(RulingCurve[AddPtIndex]->ControllPoints.size() == 2){
                    double l = (RulingCurve[AddPtIndex]->ControllPoints[0] - RulingCurve[AddPtIndex]->ControllPoints[1]).norm();
                    double l2 = (RulingCurve[AddPtIndex]->ControllPoints[0] - _p).norm();
                    Eigen::Vector3d point = l/l2 * (_p - RulingCurve[AddPtIndex]->ControllPoints[0]) + RulingCurve[AddPtIndex]->ControllPoints[0];
                    RulingCurve[AddPtIndex]->ControllPoints.push_back(point);
                }
                if(RulingCurve[AddPtIndex]->ControllPoints.size() == 3)res = RulingCurve[AddPtIndex]->drawArc(crvPtNum);

            }
            if(outline->IsClosed() && res){
                if(RulingCurve[AddPtIndex]->getCurveType() == CurveType::bsp3)RulingCurve[AddPtIndex]->BsplineRulings(outline, DivSize, crvPtNum, curveDimention);
                if(RulingCurve[AddPtIndex]->getCurveType() == CurveType::line)RulingCurve[AddPtIndex]->LineRulings(outline, DivSize);
                if(RulingCurve[AddPtIndex]->getCurveType() == CurveType::arc)RulingCurve[AddPtIndex]->ArcRulings(outline, DivSize);
                return true;
            }
            return false;
        }else if(event == 1){
            int ptInd;
            auto DelIndex = searchPointIndex(p, ptInd, 0);
            if(DelIndex[0] == -1 && refCreases == nullptr)return false;
            bool res = false;
            if(DelIndex[0] != -1){
                RulingCurve[DelIndex[0]]->ControllPoints.erase(RulingCurve[DelIndex[0]]->ControllPoints.begin() + ptInd);
                if(RulingCurve[DelIndex[0]]->getCurveType() == CurveType::bsp3)res = RulingCurve[DelIndex[0]]->drawBspline(3, crvPtNum);
                if(RulingCurve[DelIndex[0]]->ControllPoints.size() == 0){
                    RulingCurve.erase(RulingCurve.begin() + DelIndex[0]);
                    refCrv.erase(refCrv.begin() + DelIndex[0]);
                    res = false;
                }
                RulingCurve[DelIndex[0]]->ControllPoints.shrink_to_fit();
                RulingCurve.shrink_to_fit();
                refCrv.shrink_to_fit();
                if(outline->IsClosed() && res){
                    if(RulingCurve[DelIndex[0]]->getCurveType() == CurveType::bsp3)RulingCurve[DelIndex[0]]->BsplineRulings(outline, DivSize, crvPtNum, curveDimention);
                    addRulings();
                }
            }else{
                if(refCreases != nullptr){
                    refCreases->CtrlPts.erase(refCreases->CtrlPts.begin() + ptInd);
                    if((int)refCreases->CtrlPts.size() <= curveDimention){
                        NTree_Creases.Erase(refCreases);
                        refCreases = nullptr;
                    }
                }
            }
        }
    }
    return false;
}

bool Model::MoveCurvePoint(QPointF _p, int MoveIndex, int ptInd, int curveDimention, int DivSize){

    if((MoveIndex == -1 && refCreases == nullptr) || ptInd == -1)return false;
    Eigen::Vector3d p(_p.x(), _p.y(), 0);
    bool res = false;
    if(MoveIndex != -1){
        if(RulingCurve[MoveIndex]->getCurveType() == CurveType::arc && ptInd == 0){
            bool PointOnLines = false;
            bool PointInFace = outline->IsPointInFace(p);
            for(auto&l : outline->Lines){if(MathTool::is_point_on_line(p, l->o->p, l->v->p))PointOnLines = true;}
            if(!PointInFace || PointOnLines) RulingCurve[MoveIndex]->ControllPoints[ptInd] = p;
        }
        else RulingCurve[MoveIndex]->ControllPoints[ptInd] = p;
        if(RulingCurve[MoveIndex]->getCurveType() == CurveType::bsp3)res = RulingCurve[MoveIndex]->drawBspline(curveDimention, crvPtNum);
        else if(RulingCurve[MoveIndex]->getCurveType() == CurveType::line)res = RulingCurve[MoveIndex]->drawLine();
        else if(RulingCurve[MoveIndex]->getCurveType() == CurveType::arc)res = RulingCurve[MoveIndex]->drawArc(crvPtNum);

        if(outline->IsClosed() && res){
            if(RulingCurve[MoveIndex]->getCurveType() == CurveType::line)RulingCurve[MoveIndex]->LineRulings(outline, DivSize);
            if(RulingCurve[MoveIndex]->getCurveType() == CurveType::arc)RulingCurve[MoveIndex]->ArcRulings(outline, DivSize);
        }
    }else{
        return refCreases->moveCtrlPt(p, ptInd, curveDimention);
    }
    return true;
}


bool Model::CrossDection4AllCurve(){
    if(RulingCurve.empty() || RulingCurve[0]->isempty)return false;
    for(int i = 0; i < (int)RulingCurve.size(); i++){
        std::shared_ptr<CRV> c1 = RulingCurve[i];
        if(c1->isempty)continue;
        for(int j = i + 1; j < (int)RulingCurve.size(); j++){
            std::shared_ptr<CRV> c2 = RulingCurve[j];
            if(c2->isempty)continue;
            Eigen::Vector3d q;
            for(auto& r1: c1->Rulings){
                for(auto& r2: c2->Rulings){
                    if(MathTool::IsIntersect(r1->v->p, r1->o->p, r2->v->p, r2->o->p, q)){
                        r2->IsCrossed = 1;
                    }
                }
            }
        }
    }
    return true;
}
