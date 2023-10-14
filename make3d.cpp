#include "make3d.h"
#include <QDebug>
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
    Rulings.clear();
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
            qDebug() << "crvs has empty";
            return;
        }
    }
    deform();
}


Eigen::Vector3d Model::SetOnGrid(QPointF& cursol, double gridsize, const QSize& S){
    int x = (int)cursol.x() % (int)gridsize, y = (int)cursol.y() % (int)gridsize;
    x = (cursol.x() - x + gridsize/2);
    y = (cursol.y() - y + gridsize/2);
    return Eigen::Vector3d(x,y,0);
}

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

void Model::ChangeFoldLineState(){
    if((int)FL.size() == 1){FoldCurveIndex = 0; return;}
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

void Model::UpdateFLOrder(int dim){
    Eigen::Vector3d UpVec(0,-1,0);
    std::vector<std::shared_ptr<FoldLine>> hasFoldingCurve;
    std::shared_ptr<Line> btm = outline->Lines.front();//一番下の辺を探索

    initializeSurfaceVertices();

    for(auto&l: outline->Lines){
        if(((l->v->p + l->o->p)/2.0).y() > ((btm->v->p + btm->o->p)/2.0).y())btm = l;
    }
    int btm_i = std::distance(outline->Lines.begin(), std::find(outline->Lines.begin(), outline->Lines.end(),btm));
    int i = btm_i;
    std::list<std::shared_ptr<FoldLine>> List_FL;

    class LineOnFL{
    public:
        std::shared_ptr<FoldLine> FL;
        double t;
        LineOnFL(std::shared_ptr<FoldLine> _FL, double _t): FL(_FL), t(_t){}
    };

    for(auto&fl: FL){
        fl->FoldingCurve.clear();
        if((int)fl->CtrlPts.size() > dim){
            hasFoldingCurve.push_back(fl);
            for(auto& l: outline->Lines){
                std::vector<std::shared_ptr<CrvPt_FL>> P = getCrossPoint(fl->CtrlPts, l->v, l->o, dim);
                if(!P.empty()){
                    for(auto&p: P){
                        std::shared_ptr<Vertex> sec, thi;
                        if(UpVec.dot((l->v->p - l->o->p).normalized()) > 0){sec = l->v->deepCopy(); thi = l->o->deepCopy();}
                        else{sec = l->o->deepCopy(); thi = l->v->deepCopy();}
                        std::shared_ptr<Vertex4d> v4d = std::make_shared<Vertex4d>(p, sec, thi); v4d->addline(l);
                        fl->FoldingCurve.push_back(v4d);
                    }

                }
            }
            fl->SortCurve();
        }
    }

    NTree_fl.clear();
    do{
       std::vector<LineOnFL> LoF;
        for(auto&fl: hasFoldingCurve){
            if(outline->Lines[i]->is_on_line(fl->FoldingCurve.front()->first->p))
                LoF.push_back(LineOnFL(fl, (fl->FoldingCurve.front()->first->p - outline->Lines[i]->o->p).norm()/(outline->Lines[i]->v->p - outline->Lines[i]->o->p).norm()));
            if(outline->Lines[i]->is_on_line(fl->FoldingCurve.back()->first->p))
                LoF.push_back(LineOnFL(fl, (fl->FoldingCurve.back()->first->p - outline->Lines[i]->o->p).norm()/(outline->Lines[i]->v->p - outline->Lines[i]->o->p).norm()));
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
    if(DebugMode::Singleton::getInstance().isdebug())qDebug() <<"order of FoldLines";
    //NTree_fl.print();
    for(auto&r: Rulings)r->hasCrossPoint = false;
    return;
}

int Model::getLayerNum(){return NTree_fl.getLayerNum();}

//曲面の輪郭がもつ三次元上の頂点を正しい位置に配置する
//端点の修正は曲面の頂点修正後に行う
void Model::SetOnVertices_outline(bool IsupdateEndPt){
    auto IsOverlapped = [](std::shared_ptr<Vertex>&p, std::shared_ptr<Vertex>&q)->bool{return ((p->p - q->p).norm() < 1e-8)? true: false;};
    auto SetEndPoint = [&](std::shared_ptr<Vertex4d>&v4d, const std::vector<std::shared_ptr<Line>>& Surface,
                           const std::vector<std::shared_ptr<Line>>& Rulings, int ind, bool IsupdateEndPt){
        int s = Surface.size();
        std::shared_ptr<Vertex> prev = nullptr, next = nullptr;
        std::shared_ptr<Vertex> downcast = std::dynamic_pointer_cast<Vertex>(v4d->first);
        for(int i = 0; i < s; i++){
            if(IsOverlapped(downcast, Surface[i]->o)){v4d->first->p3 = Surface[i]->o->p3; return;}
            if(MathTool::is_point_on_line(v4d->first->p, Surface[i]->o->p, Surface[i]->v->p)){prev = Surface[i]->o; next = Surface[i]->v;}
        }

        for(int i = 0; i < (int)FL.size(); i++){
            if(FL[i]->FoldingCurve.empty())continue;
            for(auto it = FL[i]->FoldingCurve.begin() + 1; it != FL[i]->FoldingCurve.end()-1;it++){
                if(IsOverlapped(downcast, (*it)->second)){v4d->first->p3 = (*it)->second->p3; return;}
                if(IsOverlapped(downcast, (*it)->third)){v4d->first->p3 = (*it)->third->p3; return;}
                if(MathTool::is_point_on_line(v4d->first->p, (*it)->second->p, next->p) && MathTool::is_point_on_line((*it)->second->p, prev->p, next->p)){prev = (*it)->second;continue;}
                if(MathTool::is_point_on_line(v4d->first->p, (*it)->second->p, prev->p) && MathTool::is_point_on_line((*it)->second->p, prev->p, next->p)){next = (*it)->second;continue;}
                if(MathTool::is_point_on_line(v4d->first->p, (*it)->third->p, next->p) && MathTool::is_point_on_line((*it)->third->p, prev->p, next->p)){prev = (*it)->third;continue;}
                if(MathTool::is_point_on_line(v4d->first->p, (*it)->third->p, prev->p) && MathTool::is_point_on_line((*it)->third->p, prev->p, next->p)){next = (*it)->third;continue;}
                if(i == ind)continue;
                std::shared_ptr<Vertex> downcast_f = std::dynamic_pointer_cast<Vertex>((*it)->first);
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
        if(prev == nullptr || next == nullptr)return;
        if(IsupdateEndPt)v4d->first->p3 = (prev->p - v4d->first->p).norm()*(next->p3 - prev->p3).normalized() + prev->p3;
        v4d->second = next; v4d->third = prev;
    };
    int s = outline->vertices.size();
    initializeSurfaceVertices();
    for(auto&fl:FL)fl->SortCurve();
    for(int i = 0; i < s; i++){
        std::shared_ptr<Vertex> p = outline->vertices[i];
        bool Isvertexlapped = false;
        std::shared_ptr<Vertex> prev = outline->vertices[(i - 1 + s)  % s], next = outline->vertices[(i + 1) % s];//前後の頂点で作られるlineのうち最も近いものを入れる
        std::shared_ptr<Vertex> basePt = nullptr;
        for(auto&fl: FL){
            for(auto&v4d: fl->FoldingCurve){
                if(!v4d->IsCalc)continue;
                std::shared_ptr<Vertex> downcast = std::dynamic_pointer_cast<Vertex>(v4d->first);
                if(IsOverlapped(downcast, outline->vertices[i])){Isvertexlapped = true; basePt = v4d->first; break;}
                if(IsOverlapped(v4d->second, outline->vertices[i])){Isvertexlapped = true; basePt = v4d->second; break;}
                if(IsOverlapped(v4d->third, outline->vertices[i])){Isvertexlapped = true; basePt = v4d->third; break;}
                if(MathTool::is_point_on_line(v4d->first->p, prev->p, outline->vertices[i]->p) && MathTool::is_point_on_line(v4d->first->p, outline->vertices[i]->p, prev->p)){ prev = v4d->first;}
                if(MathTool::is_point_on_line(v4d->first->p, next->p, outline->vertices[i]->p) && MathTool::is_point_on_line(v4d->first->p, next->p, outline->vertices[i]->p)){ next = v4d->first;}
                if(MathTool::is_point_on_line(v4d->second->p, prev->p, outline->vertices[i]->p) && MathTool::is_point_on_line(v4d->second->p, outline->vertices[i]->p, prev->p)){
                    prev = v4d->second;
                    if(!IsOverlapped(downcast, prev) && !IsOverlapped(downcast, next))basePt = v4d->first;
                }
                if(MathTool::is_point_on_line(v4d->second->p, next->p, outline->vertices[i]->p) && MathTool::is_point_on_line(v4d->second->p, next->p, outline->vertices[i]->p)){
                    next = v4d->second;
                    if(!IsOverlapped(downcast, prev) && !IsOverlapped(downcast, next))basePt = v4d->first;
                }
                if(MathTool::is_point_on_line(v4d->third->p, prev->p, outline->vertices[i]->p) && MathTool::is_point_on_line(v4d->third->p, outline->vertices[i]->p, prev->p)){
                    prev = v4d->third;
                    if(!IsOverlapped(downcast, prev) && !IsOverlapped(downcast, next))basePt = v4d->first;
                }
                if(MathTool::is_point_on_line(v4d->third->p, next->p, outline->vertices[i]->p) && MathTool::is_point_on_line(v4d->third->p, next->p, outline->vertices[i]->p)){
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

                auto p = outline->vertices[i]->p;
                Eigen::Vector3d po = outline->vertices[i]->p - basePt->p, a = prev->p - basePt->p, b = next->p - basePt->p;
                double s,t;
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

    for(int i = 0; i < (int)FL.size(); i++){
        if(FL[i]->FoldingCurve.empty())continue;
        SetEndPoint(FL[i]->FoldingCurve.front(), outline->Lines, Rulings, i, IsupdateEndPt);
        SetEndPoint(FL[i]->FoldingCurve.back(), outline->Lines, Rulings, i, IsupdateEndPt);
        FL[i]->SortCurve();
    }

}

bool Model::BendingModel(double wb, double wp, int dim, double tol, int bendrank, int alg, bool ConstFunc){
    UpdateFLOrder(dim);
    SplitRulings(dim);
    std::shared_ptr<NTreeNode<std::shared_ptr<FoldLine>>> root = NTree_fl.GetRoot();
    if(root == nullptr)return false;
    if(bendrank == 0){
        AssignRuling(dim, tol, false);
        return true;
    }
    std::queue<std::shared_ptr<NTreeNode<std::shared_ptr<FoldLine>>>> q;
    std::vector<std::shared_ptr<Vertex>> Poly_V = outline->getVertices();
    q.push(root);
    int rank = 0;
    int bendnum = 0;
    while(!q.empty()){
        if(bendnum == bendrank){
            AssignRuling(dim, tol, false);
            return true;
        }
        auto cur = q.front(); q.pop();
        cur->data->ReassignColor();
        //cur->data->RevisionCrosPtsPosition();//端点の修正
        if(ConstFunc){
            qDebug() << "trim each iteration";
            bool res = cur->data->Optimization_FlapAngle(Poly_V, wb, wp, rank, 1, ConstFunc);//正しい第5引数はalgだけど検証用に1
            bool res_true = res;
            res = true;
            while(!res){
                bool isroot = (cur == root)? true: false;
                int validsize = cur->data->validsize - 1;
                cur->data->SimplifyModel(validsize, isroot);
                //cur->data->RevisionCrosPtsPosition();//端点の修正
                res = cur->data->Optimization_FlapAngle(Poly_V, wb, wp, rank, alg, ConstFunc);
                //if(DebugMode::Singleton::getInstance().isdebug())
                qDebug() << "optimization result " << res << "  ,  tol = " << tol << ", ruling num = " << cur->data->validsize;
            }
            //if(!res_true)cur->data->Optimization_FlapAngle(Poly_V, wb, wp, rank, 1, ConstFunc);//一応消しておく(flap angleの最適化に制約関数を入れてないため)
            qDebug() <<"bending result : tol = " << cur->data->tol << " valid ruling num = " << cur->data->validsize  << " , a_flap = " << cur->data->a_flap ;
            bool isroot = (root == cur)? true: false;
            cur->data->applyAAAMethod(Poly_V, false, cur->data->a_flap, cur->data->tol, isroot);
            //cur->data->SimpleSmooothSrf(Poly_V);
        }else{
            //qDebug() << "remove cross ruling and modify";
            //bool res = cur->data->Optimization_FlapAngle(Poly_V, wb, wp, rank, alg, ConstFunc);
            //cur->data->revisecrossedruling(Poly_V);
        }
        if(alg == 2) cur->data->Optimization_Vertex(Poly_V);

        SetOnVertices_outline(false);
        for (const auto& child : cur->children){
            if(child != nullptr){
                child->data->reassignruling(cur->data, outline->Lines, Rulings);
                rank++;
                q.push(child);
            }
        }     
        bendnum++;
    }
    SetOnVertices_outline(false);
    return true;
}

void Model::applyAAAMethod(double a, double tol, bool begincenter){
    auto Poly_V = outline->getVertices();
    auto root = NTree_fl.GetRoot();
    bool isroot = (root->data == FL[FoldCurveIndex])? true: false;
    FL[FoldCurveIndex]->applyAAAMethod(Poly_V, begincenter, a, tol, isroot);
}

bool Model::RevisionCrosPtsPosition(){return FL[FoldCurveIndex]->RevisionCrosPtsPosition();}

bool Model::AssignRuling(int dim, double tol, bool begincenter){
    UpdateFLOrder(dim);
    SplitRulings(dim);
    auto Poly_V = outline->getVertices();
    auto root = NTree_fl.GetRoot();
    if(root == nullptr)return false;
    std::queue<std::shared_ptr<NTreeNode<std::shared_ptr<FoldLine>>>>q;
    q.push(root);
    while (!q.empty()) {
        auto cur = q.front(); q.pop();
        if(cur->data->isbend()){
            bool isroot = (root == cur)? true: false;
            //cur->data->RevisionCrosPtsPosition();//端点の修正
            cur->data->applyAAAMethod(Poly_V, begincenter, cur->data->a_flap, cur->data->tol, isroot);
            cur->data->SimpleSmooothSrf(Poly_V);
        }
        for (const auto& child : cur->children){
            if(child != nullptr){
                child->data->reassignruling(cur->data, outline->Lines, Rulings);
                SetOnVertices_outline(true);
                q.push(child);
            }
        }
    }
    return true;
}

bool Model::SplitRulings(int dim){
    Eigen::Vector3d UpVec(0,-1,0);
    if(FL.empty() || (int)FL[FoldCurveIndex]->CtrlPts.size() <= dim)return false;
    auto root = NTree_fl.GetRoot();
    if(root == nullptr)return false;
    for(auto& r: Rulings){
        std::vector<std::shared_ptr<CrvPt_FL>> P = getCrossPoint(root->data->CtrlPts, r->v, r->o, dim);
        if(!P.empty()){
            for(auto&p: P){
                r->hasCrossPoint = true;
                std::shared_ptr<Vertex> sec, thi;
                if(UpVec.dot((r->v->p - r->o->p).normalized()) > 0){ sec = r->v->deepCopy(); thi = r->o->deepCopy();}
                else{sec = r->o->deepCopy(); thi = r->v->deepCopy(); }
                std::shared_ptr<Vertex4d> v4d = std::make_shared<Vertex4d>(p, sec, thi); v4d->addline(r);
                root->data->FoldingCurve.push_back(v4d);
            }
        }
    }
    root->data->SortCurve();
    SetOnVertices_outline(true);
    root->data->validsize = root->data->FoldingCurve.size();
    return true;
}

void Model::SimplifyModel(int iselim){
    if(!(0 <= FoldCurveIndex && FoldCurveIndex < (int)FL.size()) || FL[FoldCurveIndex]->FoldingCurve.empty())return;
    auto root = NTree_fl.GetRoot();
    bool isroot = (root->data == FL[FoldCurveIndex])? true: false;
    int validsize = FL[FoldCurveIndex]->validsize - iselim;
    FL[FoldCurveIndex]->SimplifyModel(validsize, isroot);
}

void Model::SimplifyModel(double tol){
    if(!(0 <= FoldCurveIndex && FoldCurveIndex < (int)FL.size()) || FL[FoldCurveIndex]->FoldingCurve.empty())return;
    auto root = NTree_fl.GetRoot();
    bool isroot = (root->data == FL[FoldCurveIndex])? true: false;
    FL[FoldCurveIndex]->SimplifyModel(tol, isroot);
}

bool Model::Smoothing(){
    auto root = NTree_fl.GetRoot();
    if(root == nullptr)return false;
    std::queue<std::shared_ptr<NTreeNode<std::shared_ptr<FoldLine>>>>q;
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
    for(auto&FC: FL){
        for(auto&fldCrv: FC->FoldingCurve){
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

void Model:: addConstraint(QPointF& cursol, int type, int gridsize, Eigen::Vector3d (&axis)[2], const QSize& S){
    if(outline->IsClosed()){
        qDebug()<<"constraint can be applied only not closed outline";
        return;
    }
    Eigen::Vector3d p = SetOnGrid(cursol, gridsize, S);
    if(axis[0] == Eigen::Vector3d(-1,-1,0)){axis[0] = p; return;}
    else if(axis[1] == Eigen::Vector3d(-1,-1,0) && axis[0] != p)axis[1] = p;
    Eigen::Vector3d V = (axis[1] - axis[0]).normalized();
    Eigen::Vector3d N(-V.y(), V.x(), 0);
    std::vector<std::shared_ptr<Vertex>> SymPts;
    std::vector<std::shared_ptr<Vertex>> Vertices = outline->getVertices();
    if(type == 0){
        for(auto&v: Vertices){
            double t = ((v->p - axis[0]).cross(axis[0] - axis[1])).norm()/(axis[0] - axis[1]).norm();
            if(N.dot(axis[0] - v->p) < 0) N *= -1;
            SymPts.push_back(std::make_shared<Vertex>(v->p + 2 * t * N));
        }
    }
    //鏡映反転したことで作成した曲線の始点、終点が元の曲線の始点、終点のいずれかと十分近い場合接続する。
    //基準の距離 5px

    int n = Vertices.size();
    int IsConnected = 0;
    double dist = gridsize/2;
    if((SymPts[0]->p - Vertices[0]->p).norm() < dist && (SymPts[n-1]->p - Vertices[n-1]->p).norm() < dist){
        for(int i = 1; i < (int)SymPts.size(); i++) outline->addVertex(SymPts[i], 0);
        IsConnected = 2;
    }
    else if((SymPts[0]->p - Vertices[0]->p).norm() < dist && (SymPts[n-1]->p - Vertices[n-1]->p).norm() > dist){
        for(int i = 1; i <(int)SymPts.size(); i++) outline->addVertex(SymPts[i], 0);
        IsConnected = 1;
    }
    else if((SymPts[0]->p - Vertices[0]->p).norm() > dist && (SymPts[n-1]->p - Vertices[n-1]->p).norm() < dist){
        for(int i = (int)SymPts.size() - 2; i >= 0; i--) outline->addVertex(SymPts[i], outline->getVertices().size());
        IsConnected = 1;
    }
    if(IsConnected == 0){
        ol_vertices.push_back(SymPts);
    }
    else if(IsConnected == 2){//元々あったedgeを削除してつなぎなおす
        qDebug()<<"can't use now";
    }
}

void Model::drawOutline(QPointF& cursol, int drawtype, double gridsize, const QSize& S, bool IsClicked){
    Eigen::Vector3d p = SetOnGrid(cursol, gridsize, S);
    if(drawtype == 0 || drawtype == 1){
        if(outline->hasPtNum != 2)outline->addVertex(p);
    }else if(drawtype == 2) outline->drawPolygon(p, IsClicked);

}

void Model::editOutlineVertex(QPointF& cursol, double gridsize, const QSize& S, int event){
    Eigen::Vector3d p = SetOnGrid(cursol, gridsize, S);
    static int grabedOutlineVertex = -1;
    if(event == 0){
        float dist = 5;
        grabedOutlineVertex = -1;//0~: vertex
        std::vector<std::shared_ptr<Vertex>> _vertices = outline->getVertices();
        for(int i = 0; i < (int)_vertices.size(); i++){
            if((_vertices[i]->p - p).norm() < dist){
                grabedOutlineVertex = i; dist = (_vertices[i]->p - p).norm();
            }
        }
    }else if(event == 1){
        outline->MoveVertex(p, grabedOutlineVertex);
        if(outline->IsClosed())deform();
    }else if(event == 2)grabedOutlineVertex = -1;
}

void Model::ConnectOutline(QPointF& cursol, double gridsize, const QSize& S){
    Eigen::Vector3d p = SetOnGrid(cursol, gridsize, S);
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

void Model::SplineInterPolation(const std::vector<std::shared_ptr<Line>>& path, std::vector<Eigen::Vector2d>& CurvePath){//凹型に関してはとりあえず虫
    if((int)GradationPoints.size() < 2)return;
    CurvePath.clear();
    int N =(int)GradationPoints.size() - 1;
    auto getCenter = [](const std::shared_ptr<Line>& L) { return Eigen::Vector3d(L->o->p + L->v->p)/2.0;};
    Eigen::VectorXd v(N - 1);
    std::vector<double>h(N);
    for(int i = 1; i < N + 1; i++)h[i - 1] = (getCenter(GradationPoints[i]) - getCenter(GradationPoints[i - 1])).norm();
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
        dx1 = (getCenter(GradationPoints[2]) - getCenter(GradationPoints[0])).norm();
        dx2 = (getCenter(GradationPoints[2]) -  getCenter(GradationPoints[1])).norm();
        dx3 = (getCenter(GradationPoints[1]) - getCenter(GradationPoints[0])).norm();
        if(abs(dx1) < 1e-9) u(1) = 0;
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
        qDebug()<< i << "  cur";
            return;
        }
        if(next == -1 ){
            qDebug()<< i << "  next";
            return;
        }
        Eigen::Vector3d befcenter = getCenter(GradationPoints[i]);
        Eigen::Vector3d center = getCenter(GradationPoints[i + 1]);
        double den = (befcenter - center).norm();
        a = (den != 0)? (u(i + 1) - u(i))/(6 * den) : 0;
        b = u(i)/2;
        c = - den * (2 * u(i) + u(i + 1))/6.0;
        c += (den != 0)? (GradationPoints[i + 1]->color - GradationPoints[i]->color)/den : 0;
        d = GradationPoints[i]->color;
        for(int k = cur + 1; k < next; k++){
            x =  (getCenter(path[k]) - befcenter).norm();
            y = a * std::pow(x, 3) + b * std::pow(x, 2) + c * x + d;
            path[k]->color = y;
            CurvePath.push_back(Eigen::Vector2d(cnt++,y));
        }
    }
}

void Model::setGradationValue(int val, const std::shared_ptr<Line>& refL, int InterpolationType, std::vector<Eigen::Vector2d>& CurvePath){
    if(Rulings.empty() || std::find(Rulings.begin(), Rulings.end(), refL) == Rulings.end()){qDebug()<<"no selected"; return;}
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
std::array<int, 2> Model::searchPointIndex(QPointF pt, int& ptInd, int type){
    double dist = 5.0;
    ptInd = -1;
    std::array<int, 2> CrvInds{-1,-1};
    Eigen::Vector3d p(pt.x(), pt.y(), 0);
    if(type == 0){
        for(int j = 0; j < (int)crvs.size(); j++){
            for(int i = 0; i < (int)crvs[j]->ControllPoints.size(); i++){
                Eigen::Vector3d cp = crvs[j]->ControllPoints[i];
                if(dist  > (cp - p).norm()){
                    dist = (cp - p).norm();
                    ptInd = i;
                    CrvInds[0] = j;
                }
            }
        }

    }else{
        for(int j = 0; j < (int)crvs.size(); j++){
            for(int i = 0; i < (int)crvs[j]->CurvePoints.size(); i++){
                Eigen::Vector3d cp = crvs[j]->CurvePoints[i];
                if(dist > (cp - p).norm()){
                    dist = (cp - p).norm();
                    ptInd = i;
                    CrvInds[0] = j;
                }
            }
        }
    }
    for(int i = 0; i < (int)FL.size(); i++){
        for(int j = 0; j < (int)FL[i]->CtrlPts.size(); j++){
            if(dist > (FL[i]->CtrlPts[j] - p).norm()){
                dist = (FL[i]->CtrlPts[j] - p).norm();
                ptInd = j; FoldCurveIndex = i;
            }
        }
    }
    return CrvInds;
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
    auto DelIndex = searchPointIndex(pt, ptInd, 0);
    if(DelIndex[0] == -1 && FoldCurveIndex == -1)return;
    bool res = false;
    if(DelIndex[0] != -1){
        crvs[DelIndex[0]]->ControllPoints.erase(crvs[DelIndex[0]]->ControllPoints.begin() + ptInd);
        //if(crvs[DelIndex]->getCurveType() == CurveType::bezier3)crvs[DelIndex]->drawBezier(3, crvPtNum);
        if(crvs[DelIndex[0]]->getCurveType() == CurveType::bsp3)res = crvs[DelIndex[0]]->drawBspline(3, crvPtNum);
        if(crvs[DelIndex[0]]->ControllPoints.size() == 0){
            crvs.erase(crvs.begin() + DelIndex[0]);
            refCrv.erase(refCrv.begin() + DelIndex[0]);
            res = false;
        }
        crvs[DelIndex[0]]->ControllPoints.shrink_to_fit();
        crvs.shrink_to_fit();
        refCrv.shrink_to_fit();
        if(outline->IsClosed() && res){
            if(crvs[DelIndex[0]]->getCurveType() == CurveType::bezier3)crvs[DelIndex[0]]->BezierRulings(outline, DivSize, crvPtNum);
            if(crvs[DelIndex[0]]->getCurveType() == CurveType::bsp3)crvs[DelIndex[0]]->BsplineRulings(outline, DivSize, crvPtNum, curveDimention);
            addRulings();
        }
    }else{
        if(FoldCurveIndex != -1){
            FL[FoldCurveIndex]->CtrlPts.erase(FL[FoldCurveIndex]->CtrlPts.begin() + ptInd);
            if((int)FL[FoldCurveIndex]->CtrlPts.size() <= curveDimention){
                FL.erase(FL.begin() + FoldCurveIndex);
                FoldCurveIndex = FL.size() - 1;
            }
        }
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
    Axis4Const[0] = Eigen::Vector3d(-1,-1,0);
    Axis4Const[1] = Eigen::Vector3d(-1,-1,0);
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

bool Model::AddControlPoint(Eigen::Vector3d& p, int curveDimention, int DivSize){
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
            double l = (crvs[AddPtIndex]->ControllPoints[0] - crvs[AddPtIndex]->ControllPoints[1]).norm();
            double l2 = (crvs[AddPtIndex]->ControllPoints[0] - p).norm();
            Eigen::Vector3d point = l/l2 * (p - crvs[AddPtIndex]->ControllPoints[0]) + crvs[AddPtIndex]->ControllPoints[0];
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

bool Model::MoveCurvePoint(Eigen::Vector3d& p, int MoveIndex, int ptInd, int curveDimention, int DivSize){
    if((MoveIndex == -1 && FoldCurveIndex == -1) || ptInd == -1)return false;
    bool res = false;
    if(MoveIndex != -1){
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
    }else{
        return FL[FoldCurveIndex]->moveCtrlPt(p, ptInd, curveDimention);
    }
    return true;
}

bool Model::AddControlPoint_FL(Eigen::Vector3d& p, int event, int curveDimention){
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
           Eigen::Vector3d nextPos = (std::get<0>(GradationPoints[i+1]->r->r)->p + std::get<1>(GradationPoints[i+1]->r->r)->p)/2.0;
           do{
               he = f->halfedge;
               double dist = 1000;
               nextHE = nullptr;
               do{
                   if(he->pair != nullptr){
                       double d = ((nextPos - he->vertex->p).cross((he->vertex->p - he->pair->vertex->p))).norm()/(he->vertex->p - he->pair->vertex->p).norm();
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
