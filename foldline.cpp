
#include "foldline.h"

const double eps = 1e-5;
using namespace MathTool;

void _FoldingAAAMethod(std::vector<std::pair<CrvPt_FL*, Vertex*>>& FoldingCurve, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Poly_V,
                                 double a, double phi02, double phim1, double& ff, double& feq, double& fcon, int type);
Vertex* findAxisVertex(std::pair<CrvPt_FL*, Vertex*> R, Vertex *e1, Vertex *e2);

bool IsRulingCrossed(glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint,  std::vector<Vertex*>& Poly_V);

FoldLine::FoldLine(int crvNum, int rsize, PaintTool _type)
{
    CtrlPts.clear();
    CurvePts.clear();
    Rulings_2dL.clear();
    Rulings_2dR.clear();
    Rulings_3dL.clear();
    Rulings_3dR.clear();

    CurvePts.resize(crvNum);
    maxRsize = rsize;
    color = 0;
    type = _type;
    he = he2 = nullptr;
    CtrlPts_res.clear();
    Curve_res.clear();
}

std::vector<glm::f64vec3> FoldLine::getCtrlPt(){return CtrlPts;}
std::vector<glm::f64vec3> FoldLine::getCtrlPt2d(){return CtrlPts_res2d;}
double FoldLine::getColor(){return color;}

bool FoldLine::delCtrlPt(glm::f64vec3& p, int dim, OUTLINE *outline){
    if(!outline->IsClosed()) return false;
    if(!CtrlPts.empty())CtrlPts.pop_back();
    double dist = 10;
    int ind = -1;
    for(int i = 0; i < (int)CtrlPts.size(); i++){
        if(glm::distance(p, CtrlPts[i]) < dist){
            dist = glm::distance(p, CtrlPts[i]);
            ind = i;
        }
    }
    if(ind != -1){
        CtrlPts.erase(CtrlPts.begin() + ind);
    }
    bool hasCrv = setCurve(dim);
    if(outline->IsClosed() && hasCrv){
        //bool res = applyCurvedFolding(dim, outline);
        //return res;
    }
}

bool FoldLine::moveCtrlPt(glm::f64vec3& p, int movePtIndex){
    if((int)CtrlPts.size() < movePtIndex || movePtIndex < 0)return false;
    CtrlPts[movePtIndex] = p;
    if(FoldingCurve.empty())return false;
    return true;

}

bool FoldLine::addCtrlPt(glm::f64vec3& p, int dim){
    //if((int)CtrlPts.size() <= dim)
    CtrlPts.push_back(p);
    bool hasCrv = setCurve(dim);

    return hasCrv;
}


//type == 0 line
//type == 1 b-spline
bool FoldLine::setCurve(int dim){
    using namespace MathTool;
    int crvPtNum = CurvePts.size();
    if(type == PaintTool::FoldLine_line){
        if(CtrlPts.size() < 2)return false;
        while(CtrlPts.size() != 2){
            CtrlPts.erase(CtrlPts.end() - 1);
            CtrlPts.shrink_to_fit();
        }
        glm::f64vec3 V = CtrlPts[1] - CtrlPts[0];
        for(int i = 0; i < crvPtNum; i++){
            CurvePts[i] = (double)i/(double)(crvPtNum - 1) * V + CtrlPts[0];
        }
    }
    else if(type == PaintTool::FoldLine_arc){
        if((int)CtrlPts.size() < 3) return false;

        double l = glm::distance(CtrlPts[0], CtrlPts[1]);
        double l2 = glm::distance(CtrlPts[0], CtrlPts[2]);
        CtrlPts[2] = l/l2 * (CtrlPts[2] - CtrlPts[0]) + CtrlPts[0];
        double phi = acos(glm::dot(glm::normalize(CtrlPts[1] - CtrlPts[0]), glm::normalize(CtrlPts[2] - CtrlPts[0])));
        glm::f64vec3 axis = glm::cross(glm::normalize(CtrlPts[1] - CtrlPts[0]), glm::normalize(CtrlPts[2] - CtrlPts[0]));
        glm::f64mat4x4  T = glm::translate(CtrlPts[0]);
        glm::f64mat4x4 invT = glm::translate(-CtrlPts[0]);
        glm::f64mat4x4  R;
        for(int i = 0; i < crvPtNum; i++){
            R = glm::rotate(phi * (double)i/(double)crvPtNum, axis);
            CurvePts[i] = T * R * invT * glm::f64vec4{CtrlPts[1],1};
        }

    }
    else if(type == PaintTool::FoldLine_bezier){
        if((int)CtrlPts.size() <= dim)return false;
        if((int)CtrlPts.size() > dim + 1)CtrlPts.pop_back();
        double t = 0.0;
        glm::f64vec3 v;
        for(int n = 0; n < crvPtNum; n++){
            v = {0.f,0.f, 0.f};
            for (int i = 0; i < int(CtrlPts.size()); i++) {
                v += MathTool::BernsteinBasisFunc(dim, i, t) * CtrlPts[i];
            }
            CurvePts[n] = v;
            t += 1/(double)crvPtNum;
        }
    }

    return true;
}

bool IsRulingCrossed (glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint,  std::vector<Vertex*>& Poly_V){
    double l = 1000, minDist = 1000;
    bool IsIntersected = false;
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    N = glm::normalize(N);
    int n = Poly_V.size();
    for(int i = 0; i < n; i++){
        glm::f64vec3 v = Poly_V[(i+1) % n]->p - Poly_V[i]->p;
        A(0,0) = v.x; A(0,1) = -l*N.x; A(1,0) = v.y; A(1,1) = -l * N.y;
        b(0) = cp.x - Poly_V[i]->p.x; b(1) = cp.y - Poly_V[i]->p.y;
        if(abs(glm::dot(glm::normalize(v), N))>= 1-DBL_EPSILON)continue;
        Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
        if (0 <= x(1) && x(1) <= 1) {
            glm::f64vec3 p = x(1)*l*N + cp;
            if(glm::distance(p, cp) < minDist){
                crossPoint = p;
                minDist = glm::distance(p, cp);
            }
            IsIntersected = true;
        }
    }
    return IsIntersected;
}


double Fbend(std::vector<Vertex*>& Poly_V,std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices){
    auto Cot_jk = [](glm::f64vec3 ej, glm::f64vec3 ek){
        double a = glm::angle(ej, ek);
        return 1.0/std::tan(a);
    };

    double f = 0.0;
    //三角形分割
    //すべてのエッジに対して曲げエネルギーの計算
    //相和を返す
    std::vector<HalfEdge*> _edges = EdgeCopy(Edges, Vertices);
    std::vector<Face*> _faces; EdgeRecconection(Poly_V,_faces, _edges);
    /*
    for(int i = 0; i < (int)_faces.size(); i++)_faces[i]->TrianglationSplit(_edges,_faces);
    for(int i = 0; i < (int)_edges.size(); i++){
        HalfEdge *h = _edges[i];
        if(h->edgetype == EdgeType::ol || h->pair == nullptr || (h->face->edgeNum(false) != 3 && h->pair->face->edgeNum(false) != 3))continue;
        glm::f64vec3 _e0 = (h->next->vertex->p3 - h->vertex->p3);
        glm::f64vec3 _e1 = (h->vertex->p3 - h->prev->vertex->p3);
        glm::f64vec3 _e3 = (h->next->next->vertex->p3 - h->next->vertex->p3);
        glm::f64vec3 _e2 = (h->pair->next->next->vertex->p3 - h->vertex->p3);
        glm::f64vec3 _e4 = (h->pair->vertex->p3 - h->pair->prev->vertex->p3);

        double A0 = glm::length(glm::cross(_e0, -_e1))/2.0, A1 = glm::length(glm::cross(_e0, _e2))/2.0;
        _e0 = glm::normalize(_e0); _e1 = glm::normalize(_e1); _e2 = glm::normalize(_e2); _e3 = glm::normalize(_e3); _e4 = glm::normalize(_e4);
        glm::f64vec4 K{Cot_jk(-_e0,_e3)+ Cot_jk(-_e0,-_e4), Cot_jk(_e0,-_e1)+ Cot_jk(_e0,_e2), -Cot_jk(_e0,-_e1)- Cot_jk(-_e0,_e3), -Cot_jk(_e0,_e2)- Cot_jk(-_e0,-_e4)};

        f += 3.0/(2.0*(A0 + A1))*glm::dot(K,K)*(glm::dot(h->vertex->p3,h->vertex->p3) + glm::dot(h->next->vertex->p3,h->next->vertex->p3) +
                                                glm::dot(h->prev->vertex->p3,h->prev->vertex->p3) +glm::dot(h->pair->prev->vertex->p3,h->pair->prev->vertex->p3));
    }*/

    for(const auto&e: _edges){
        if(e->pair == nullptr)continue;
        glm::f64vec3 N = e->face->getNormalVec();
        glm::f64vec3 Np = e->pair->face->getNormalVec();
        double phi = (abs(glm::dot(N,Np)) > 1)? std::numbers::pi: std::numbers::pi - std::acos(glm::dot(N,Np));
        f += 1.0/(phi*phi);
    }

    return f;
}

double FoldLine::Fruling(){
    double f = 0.0;
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    //int PSize = Poly_V.size();
    for(int i = 2; i < (int)FoldingCurve.size() - 3; i++){
        A(0,0) = FoldingCurve[i].second->p.x - FoldingCurve[i].first->p.x; A(0,1) = -(FoldingCurve[i-1].second->p.x - FoldingCurve[i-1].first->p.x);
        A(1,0) = FoldingCurve[i].second->p.y - FoldingCurve[i].first->p.y; A(1,1) = -(FoldingCurve[i-1].second->p.y - FoldingCurve[i-1].first->p.y);
        b(0) = FoldingCurve[i-1].first->p.x - FoldingCurve[i].first->p.x;
        b(1) = FoldingCurve[i-1].first->p.y - FoldingCurve[i].first->p.y;
        Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
        if(0 < ts(0) && ts(0) < 1) f += 1.0/ts(0);

        A(0,0) = FoldingCurve[i].second->p.x - FoldingCurve[i].first->p.x; A(0,1) = -(FoldingCurve[i+1].second->p.x - FoldingCurve[i+1].first->p.x);
        A(1,0) = FoldingCurve[i].second->p.y - FoldingCurve[i].first->p.y; A(1,1) = -(FoldingCurve[i+1].second->p.y - FoldingCurve[i+1].first->p.y);
        b(0) = FoldingCurve[i+1].first->p.x - FoldingCurve[i].first->p.x;
        b(1) = FoldingCurve[i+1].first->p.y - FoldingCurve[i].first->p.y;
        ts = A.colPivHouseholderQr().solve(b);
        if(0 < ts(0) && ts(0) < 1) f += 1.0/ts(0);
    }
    return f;
}

//https://programming.pc-note.net/c/struct.html
typedef struct{
    std::vector<std::pair<CrvPt_FL*, Vertex*>> FoldingCurve;
    std::vector<HalfEdge*> Edges;
    std::vector<Vertex*> Vertices, Poly_V;
    double phi02, phim1;
    int type;
}ObjData;//初期化の変数代入の順番が大事

double ObjFunc(const std::vector<double> &a, std::vector<double> &grad, void* f_data){
    double h = 1e-7;
    double ff, feq, fcon;
    ObjData *od = (ObjData *)f_data;
    double phi02 = od->phi02, phim1 = od->phim1;
    int type = od->type;
    std::vector<std::pair<CrvPt_FL*, Vertex*>> FoldingCurve = od->FoldingCurve;
    std::vector<HalfEdge*> Edges = od->Edges;
    std::vector<Vertex*> Vertices = od->Vertices, Poly_V = od->Poly_V;
    if(!grad.empty()){
        _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a[0] + h, phi02, phim1, ff, feq, fcon, type);
       double fp = Fbend(Poly_V, Edges, Vertices);
        _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a[0] - h, phi02, phim1, ff, feq, fcon, type);
       double fm = Fbend(Poly_V, Edges, Vertices);
       grad[0] = (fp - fm)/(2.0 * h);
    }

    return Fbend(Poly_V, Edges, Vertices);
}

bool FoldLine::Optimization(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, int type){
    if(FoldingCurve.empty())return false;

    struct CrossRegion{
        double first, second;
        bool IsContinuation;
        CrossRegion(double a, bool _IsContinuation): first(a), second(-1), IsContinuation(_IsContinuation){}
    };

    double phi02, phim1;
    phi02 = glm::angle(glm::normalize(FoldingCurve[0].second->p - FoldingCurve[0].first->p),glm::normalize(FoldingCurve[1].first->p - FoldingCurve[0].first->p));
    phim1 = glm::angle(glm::normalize(FoldingCurve.end()[-2].first->p - FoldingCurve.back().first->p),glm::normalize(FoldingCurve.back().second->p - FoldingCurve.back().first->p));
    double a = 0, a2 = 0.0;
    std::ofstream ofs2;
    std::filesystem::create_directory("./Optimization");

    std::string AngleFile = "./Optimization/ChangeAngle.csv";
    ofs2.open(AngleFile, std::ios::out);
    ofs2 << "a, Fruling, f_ff, f_eq, fcon, f_con'"<<std::endl;

    int maxitr = 100;
    double wb = 1e-3, wr = 1, wa = 1e-2;
    double h = 1e-6, th = 1e-6;
    double ff, feq, fcon, fcon_bef;

    std::vector<CrossRegion> NoCrossRegion;
    bool IsCrossedBefore = false;
    while(a2 <= 2.0 * std::numbers::pi){
        a = a2;
        double fcon = 0, ff = 0, feq = 0;
        _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a, phi02, phim1, ff, feq, fcon, type);
        double fr = Fruling();
        if(fr <= DBL_EPSILON){
            if(IsCrossedBefore) NoCrossRegion.push_back(CrossRegion{a2, false});
            else NoCrossRegion.back().IsContinuation = (fcon < 1e-2)? true: NoCrossRegion.back().IsContinuation;
            IsCrossedBefore = false;
        }else{
            if(!IsCrossedBefore && !NoCrossRegion.empty()) NoCrossRegion.back().second = a2;
            IsCrossedBefore = true;
        }

        double diff_fcon = (a2 == 0.0)? 0.0: (fcon - fcon_bef)/(std::numbers::pi/(double)maxitr);
        fcon_bef = fcon;
        ofs2 << glm::degrees(a2) << ", "  << fr << ", " << ff << ", " << feq << ", " << fcon << ", "  << diff_fcon << std::endl;
        a2 += 1.0/(double)maxitr;
    }
    ofs2.close();
    for(auto&L: NoCrossRegion){
        if(L.IsContinuation){
            std::cout << "Continuation : " <<  glm::degrees(L.first) << " <=  a <= " << glm::degrees(L.second) << std::endl;
        }else{
            std::cout << "No Continuation : " <<  glm::degrees(L.first) << " <=  a <= " << glm::degrees(L.second) << std::endl;
            a = (L.first + L.second)/2.0;
            {

                std::ofstream ofs;
                std::string OptFile =  "./Optimization/Optimization.csv";
                ofs.open(OptFile, std::ios::out);
                ofs << "a, Fbend, Fruling, f" << std::endl;
                ObjData od = {FoldingCurve, Edges, Vertices, Poly_V, phi02, phim1, type};
                nlopt::opt opt(nlopt::LD_MMA, 1);
                opt.set_lower_bounds(L.first);
                opt.set_upper_bounds(L.second);
                opt.set_min_objective(ObjFunc, &od);
                opt.set_xtol_rel(1e-9);
                double minf;
                  try {
                    std::vector<double> _a{a};
                      nlopt::result result = opt.optimize(_a, minf);
                      std::cout << "found minimum at f(" << glm::degrees(_a[0]) << ") = "
                          << std::setprecision(10) << minf << std::endl;
                  }
                  catch (std::exception& e) {
                      std::cout << "nlopt failed: " << e.what() << std::endl;
                  }
                 /*
                auto hasBoundary = [](double a, double a_min, double a_max){
                    double h = 1e-7;
                    double b = (a - a_min > h)? std::log(a - a_min): std::log(h), b2 = (a_max - a  > h)? std::log(a_max - a): std::log(1e-7);
                    return -(b + b2);
                };
                double mu = 1;
                a = (L.first + L.second)/2.0;
                for(int i = 0; i < maxitr; i++){
                    mu = 1.0/(double)(i+ 1);
                    for(int j = 0; j < maxitr; j++){

                         _FoldingAAAMethod(Edges, Vertices, Poly_V, a + h, phi02, phim1, ff, feq, fcon, type);
                        double fp = wb * Fbend(Poly_V, Edges, Vertices) +  mu * hasBoundary(a + h, L.first, L.second);
                         _FoldingAAAMethod(Edges, Vertices, Poly_V, a - h, phi02, phim1, ff, feq, fcon, type);
                        double fm = wb * Fbend(Poly_V, Edges, Vertices) +  mu * hasBoundary(a - h, L.first, L.second);
                        double da = (fp - fm)/(2.0 * h);
                        if(abs(da) < th)break;
                        a2 = a - wa * da;
                        a2 = (a2 < 0) ? 2.0*std::numbers::pi - 2.0 * std::numbers::pi: (a2 > 2.0 * std::numbers::pi)? a2 - (double)(floor(a2 / (2.0 * std::numbers::pi))) * 2.0 * std::numbers::pi: a2;
                        a = a2;
                    }
                     _FoldingAAAMethod(Edges, Vertices, Poly_V, a, phi02, phim1, ff, feq, fcon, type);
                    double f = wb * Fbend(Poly_V, Edges, Vertices) +  mu * hasBoundary(a, L.first, L.second);
                    a = a2;
                    if(f  < th)break;
                    ofs << glm::degrees(a2) << ", " << wb * Fbend(Poly_V, Edges, Vertices) <<", " <<  wr * Fruling() << ", " << f << std::endl;
                }*/
                ofs.close();
            }


            _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a2, phi02, phim1, ff, feq, fcon, type);

        }
    }


    std::cout << "finish"<<std::endl;
    return true;
}


bool FoldLine::modify2DRulings(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_v, int dim, int t_type){
    using namespace MathTool;
    auto FindAnotherPoint = [](const std::vector<Vertex*>& Poly_v, std::vector<HalfEdge*>& Edges, Vertex *v, glm::f64vec3 n){
        int N = Poly_v.size();
        for(int i = 1; i <= N; i++){
            if(MathTool::is_point_on_line(v->p, Poly_v[i % N]->p, Poly_v[(i - 1) % N]->p)){
                for(const auto& e: Edges){
                    if(e->vertex->p == v->p){
                        if(glm::dot(glm::normalize(Poly_v[(i - 1) % N]->p - v->p), glm::normalize(n)) > glm::dot(glm::normalize(Poly_v[i % N]->p - v->p), glm::normalize(n)))return Poly_v[(i - 1) % N];
                        return Poly_v[i % N];
                    }
                }
            }
        }
        return static_cast<Vertex*>(nullptr);
    };
    std::vector<CrvPt_FL*> T_crs;
    std::vector<glm::f64vec3>edges_ol;
    std::vector<HalfEdge*> edges_he;
    double t_max = -1, t_min = 1;
    for(auto&e: Poly_v){
        glm::f64vec3 p = e->p;
        edges_he.push_back(new HalfEdge(new Vertex(p), EdgeType::ol));
        edges_ol.push_back(p);
    }
    for(int i = 0; i < (int)edges_he.size(); i++){
        edges_he[i]->prev = edges_he[(i + 1) % edges_he.size()];
        edges_he[i]->next = edges_he[(i - 1) % edges_he.size()];
    }
    Face *face_ol = new Face(edges_he[0]);

    std::vector<HalfEdge*> SearchedEdge;  
    if(type == PaintTool::FoldLine_line){
        int divSize = 1;

    }else if(type == PaintTool::FoldLine_arc){
        for(auto& e: Edges){
            if((e->edgetype == EdgeType::r && (std::find(SearchedEdge.begin(), SearchedEdge.end(), e) != SearchedEdge.end() || std::find(SearchedEdge.begin(), SearchedEdge.end(), e->pair) != SearchedEdge.end()))
                    || e->edgetype != EdgeType::r)continue;
            SearchedEdge.push_back(e);
            glm::f64vec3 v{0,0,0};
            double x1 = e->vertex->p.x, y1 = e->vertex->p.y, x2 = e->next->vertex->p.x, y2 = e->next->vertex->p.y;
            double xd = x2 - x1, yd = y2 - y1, X = x1 - CtrlPts[0].x, Y = y1 - CtrlPts[0].y;
            double a = xd * xd + yd * yd, b = xd * X + yd * Y, c = X * X + Y * Y - pow(glm::distance(CtrlPts[1], CtrlPts[0]),2);
            double D = b * b - a*c;
            if(D < 0) continue;
            double s1 = (-b + sqrt(D)) / a;
            double s2 = (-b - sqrt(D)) / a;
            if(0 <= s1 && s1 <= 1 && face_ol->IsPointInFace(glm::f64vec3{x1 + xd*s1, y1 + yd*s1, 0}))v = glm::f64vec3{x1 + xd*s1, y1 + yd*s1, 0};
            if(0 <= s2 && s2 <= 1 && face_ol->IsPointInFace(glm::f64vec3{x1 + xd*s1, y1 + yd*s1, 0}))v = glm::f64vec3{x1 + xd*s2, y1 + yd*s2, 0};

            if(!is_point_on_line(v, e->vertex->p, e->next->vertex->p)) continue;
            double sa = glm::distance(v, e->vertex->p), sc = glm::distance(e->vertex->p, e->next->vertex->p);
            glm::f64vec3 v3 = sa/sc * (e->next->vertex->p3 - e->vertex->p3) + e->vertex->p3;
            CrvPt_FL *P = new CrvPt_FL{v, v3, 0};
            T_crs.push_back(P);

        }
    }else if(type == PaintTool::FoldLine_bezier){
        for(auto& e: Edges){
            if((e->edgetype == EdgeType::r && (std::find(SearchedEdge.begin(), SearchedEdge.end(), e) != SearchedEdge.end() || std::find(SearchedEdge.begin(), SearchedEdge.end(), e->pair) != SearchedEdge.end()))
                    )continue;
            SearchedEdge.push_back(e);
            std::vector<double>arcT = BezierClipping(CtrlPts, e, dim);

            for(auto&t: arcT){
                if(t < 0 || 1 < t){std::cout<<"t is not correct value " << t << std::endl; continue;}
                glm::f64vec3 v2{0,0,0};
                for (int i = 0; i < int(CtrlPts.size()); i++) v2 += MathTool::BernsteinBasisFunc(dim, i, t) * CtrlPts[i];
                if(!is_point_on_line(v2, e->vertex->p, e->next->vertex->p)) continue;
                double sa = glm::distance(v2, e->vertex->p), sc = glm::distance(e->vertex->p, e->next->vertex->p);
                glm::f64vec3 v3 = sa/sc * (e->next->vertex->p3 - e->vertex->p3) + e->vertex->p3;
                //VectorDigAlign(v2); VectorDigAlign(v3);
                t_max = std::max(t_max, t); t_min = std::min(t_min, t);
                CrvPt_FL *P = new CrvPt_FL(v2, v3, t);
                T_crs.push_back(P);
            }
        }
    }else if(type == PaintTool::FoldLine_test){
        int faceVerticesSize = edges_ol.size();
        for(int i = 0; i < faceVerticesSize; i++){
            bool hasPoint = MathTool::IsIntersect(edges_ol[i], edges_ol[(i+1)%faceVerticesSize], CtrlPts[0], CtrlPts[1], true);
            if(hasPoint){
                glm::f64vec3 NewPoint = MathTool::getIntersectionPoint(edges_ol[i], edges_ol[(i+1)%faceVerticesSize], CtrlPts[0], CtrlPts[1]);
                CrvPt_FL *P = new CrvPt_FL(NewPoint, i);
                T_crs.push_back(P);
            }
        }
    }
    std::sort(T_crs.begin(), T_crs.end(), [](CrvPt_FL* T, CrvPt_FL* T2){
        return T->s > T2->s;
    });
    double tol = 1.0*std::numbers::pi/ 180.0;
    TrimPoints(Edges, Faces, Vertices, T_crs, tol);

    for(int i = 0; i < Faces.size(); i++){
        std::cout << i << " , " << Faces[i]->edgeNum() << std::endl;
    }
    std::cout <<"-----------------------"<<std::endl;

    for(auto&t: T_crs){
        Vertices.push_back(t);
    }
    for(auto&t: T_crs){
        std::vector<HalfEdge*> H_new;
        for(auto&he: Edges){
            H_new = he->Split(t, Edges);
            if(H_new.size() > 0)break;          
        }
    }

    std::vector<Face*> Faces_new;
    FoldingCurve.clear();
    double tbef = -1;
    int faceNum = Faces.size();
    std::pair<CrvPt_FL*, Vertex*> R;
    HalfEdge_FL *h1, *h2;
    for(int i = 0; i < faceNum; i++){
        Face *f = Faces[i];
        HalfEdge *h = f->halfedge;
        tbef = -1;
        std::vector<HalfEdge*> InsertEdges;
        std::vector<CrvPt_FL*> InsertPoints;
        do{
            for(auto&t: T_crs){//indexが小さいほうがtの値が小さい
                if(t->p == h->vertex->p){
                if(tbef == -1){
                    InsertEdges.push_back(h);
                    InsertPoints.push_back(t);
                }
                else if(tbef != -1 && tbef < t->s){
                    InsertEdges.push_back(h);
                    InsertPoints.push_back(t);
                }
                else if(tbef != -1 && tbef > t->s){
                    InsertEdges.insert(InsertEdges.begin() + 0, h);
                    InsertPoints.insert(InsertPoints.begin()+0, t);
                }
                tbef = t->s;
                break;
                }
            }
            h = h->next;
        }while(h != f->halfedge);

        if(InsertEdges.size() == 2){

            h1 = new HalfEdge_FL(InsertPoints[0], EdgeType::fl);
            h2 = new HalfEdge_FL(InsertPoints[1], EdgeType::fl);
            h1->pair = h2; h2->pair = h1;
            h1->prev = InsertEdges[0]->prev; h2->prev = InsertEdges[1]->prev;
            h1->next = InsertEdges[1]; h2->next = InsertEdges[0];
            InsertEdges[1]->prev->next = h2; InsertEdges[0]->prev->next = h1;
            InsertEdges[0]->prev = h2; InsertEdges[1]->prev = h1;

            Edges.push_back(h1); Edges.push_back(h2);
            Face *fn = new Face(h1);
            f->ReConnect(h2);
            fn->ReConnect(h1);
            Faces_new.push_back(fn);
            if(h1->v->s > h2->v->s){ R = {InsertPoints[1], he->prev->vertex};FoldingCurve.push_back(R);}
            else {R = {InsertPoints[0], h1->prev->vertex};FoldingCurve.push_back(R);}
        }
    }
    glm::f64vec3 befN = glm::normalize(FoldingCurve.front().second->p - FoldingCurve.front().first->p);


    for(int i = 0; i < (int)FoldingCurve.size(); i++){
        if(FoldingCurve[i].first == T_crs.front() || FoldingCurve[i].first == T_crs.back() )FoldingCurve.erase(FoldingCurve.begin() + i);
    }

    std::sort(FoldingCurve.begin(), FoldingCurve.end(), [](std::pair<CrvPt_FL*, Vertex*>& F, std::pair<CrvPt_FL*, Vertex*>& F2){
        return F.first->s > F2.first->s;
    });
    /*
    for(int i = 0; i < (int)FoldingCurve.size(); i++){
        for(int j = i + 1; j < (int)FoldingCurve.size(); j++){
            if(FoldingCurve[i].first->s < FoldingCurve[j].first->s){
               auto tmp = FoldingCurve[i];
                FoldingCurve[i] = FoldingCurve[j];
               FoldingCurve[j] = tmp;
                i = 0;
           }
        }
    }*/

    R = {T_crs.front(), FindAnotherPoint(Poly_v, Edges, T_crs.front(), befN)}; FoldingCurve.insert(FoldingCurve.begin(), R);
    R = {T_crs.back(), FindAnotherPoint(Poly_v, Edges, T_crs.back(), befN)}; FoldingCurve.push_back(R);
    //FoldingCurve.erase(std::unique(FoldingCurve.begin(), FoldingCurve.end()), FoldingCurve.end());
    Faces.insert(Faces.end(), Faces_new.begin(), Faces_new.end()); // 連結
    for(int i = 0; i < Faces.size(); i++){
        std::cout << i << " , " << Faces[i]->edgeNum() << std::endl;
    }
    return true;
}

double ObjFunc2(const std::vector<double> &Phi, std::vector<double> &grad, void* f_data){
    double h = 1e-7;
    double ff, feq, fcon;
    ObjData *od = (ObjData *)f_data;
    double phi02 = od->phi02, phim1 = od->phim1;
    int type = od->type;
    std::vector<std::pair<CrvPt_FL*, Vertex*>> FoldingCurve = od->FoldingCurve;
    std::vector<HalfEdge*> Edges = od->Edges;
    std::vector<Vertex*> Vertices = od->Vertices, Poly_V = od->Poly_V;
    if(!grad.empty()){

    }

    return Fbend(Poly_V, Edges, Vertices);
}

bool FoldLine::SetNewPointsOnCurve(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_v, int dim, int divSize){
    using namespace MathTool;

    if(divSize < 1 || Poly_v.empty() || (int)CtrlPts.size() <= dim)return false;

    double t_max = -1, t_min = 2;
    std::vector<HalfEdge*> SearchedEdge;
    for(auto& e: Edges){
        if(std::find(SearchedEdge.begin(), SearchedEdge.end(), e) != SearchedEdge.end() || (e->pair != nullptr && std::find(SearchedEdge.begin(), SearchedEdge.end(), e->pair) != SearchedEdge.end()))continue;
        SearchedEdge.push_back(e);
        std::vector<double>arcT = BezierClipping(CtrlPts, e, dim);
        for(auto&t: arcT){
            if(t < 0 || 1 < t){std::cout<<"t is not correct value " << t << std::endl; continue;}
            t_max = std::max(t_max, t); t_min = std::min(t_min, t);
        }
    }

    Points_On_Curve.clear();
    for(int i = 0; i < divSize; i++){
        double t = (double)(i + 1) * (t_max - t_min)/(double)(divSize + 1) + t_min;
        glm::f64vec3 v{0,0,0};
        for (int i = 0; i < int(CtrlPts.size()); i++) v += MathTool::BernsteinBasisFunc(dim, i, t) * CtrlPts[i];
        CrvPt_FL *vertex = new CrvPt_FL(v,v,t); Vertices.push_back(vertex);
        Points_On_Curve.push_back(vertex);
    }
    std::sort(Points_On_Curve.begin(), Points_On_Curve.end(), [](CrvPt_FL* F, CrvPt_FL* F2){
        return F->s > F2->s;
    });

    glm::f64vec3 v{0,0,0};
    for (int i = 0; i < int(CtrlPts.size()); i++) v += MathTool::BernsteinBasisFunc(dim, i, t_max) * CtrlPts[i];
    CrvPt_FL *v_max = new CrvPt_FL(v, v, t_max); Vertices.push_back(v_max);

    v = glm::f64vec3{0,0,0};
    for (int i = 0; i < int(CtrlPts.size()); i++) v += MathTool::BernsteinBasisFunc(dim, i, t_min) * CtrlPts[i];
    CrvPt_FL *v_min = new CrvPt_FL(v, v, t_min); Vertices.push_back(v_min);

    return true;
}

void FoldLine::Optimization2(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_v, double a){
    if(Points_On_Curve.empty() || a < DBL_EPSILON || abs(a - std::numbers::pi) < DBL_EPSILON || abs(a - 2.0*std::numbers::pi) < DBL_EPSILON)return;
    std::vector<double> Phi4;
    for(int i = 1; i < (int)Points_On_Curve.size() - 1; i++){
        glm::f64vec3 x = Points_On_Curve[i]->p;
        glm::f64vec3 e = glm::normalize(Points_On_Curve[i-1]->p - x);
        glm::f64vec3 e2 = glm::normalize(Points_On_Curve[i+1]->p - x);
        double phi = (glm::dot(glm::cross(e, e2), glm::f64vec3{0,0,1}) < 0)? std::acos(glm::dot(e, e2)): 2.0*std::numbers::pi - std::acos(glm::dot(e, e2));
        Phi4.push_back(phi/2.0);
    }
    nlopt::opt opt(nlopt::LD_MMA, (int)Phi4.size());

    opt.set_min_objective(ObjFunc, NULL);
    opt.set_xtol_rel(1e-9);
    double minf;
      try {
          nlopt::result result = opt.optimize(Phi4, minf);

      }
      catch (std::exception& e) {
          std::cout << "nlopt failed: " << e.what() << std::endl;
      }

}

void FoldLine::TrimPoints(std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces, std::vector<Vertex*>& Vertices, std::vector<CrvPt_FL*>& T_crs, double tol){
    //すべての交点に対してkを導出
    int Ind;
    double mink;
    do{
        Ind = -1;
        mink = tol;
        for(int i = 1; i < (int)T_crs.size() - 1; i++){
            double k = std::acos(glm::dot(glm::normalize(T_crs[i-1]->p - T_crs[i]->p), glm::normalize(T_crs[i+1]->p - T_crs[i]->p)));
            double k_diff = (std::numbers::pi - k);
            if(k_diff == std::min(mink, k_diff)){//最もpiに近くてtolよりも差が小さいkをもつ交点を選択
                Ind = i;
                mink = k_diff;
            }
        }
        if(Ind != -1){
            //交点とruling、Edgeを削除
            for(auto&h: Edges){
                if(MathTool::is_point_on_line(T_crs[Ind]->p, h->vertex->p, h->next->vertex->p)){
                    std::vector<HalfEdge*>::iterator itr_hpn = std::find(Edges.begin(), Edges.end(), h->pair->next);
                    std::vector<HalfEdge*>::iterator itr_hn = std::find(Edges.begin(), Edges.end(), h->next);
                    std::vector<HalfEdge*>::iterator itr_hp = std::find(Edges.begin(), Edges.end(), h->pair);
                    h->prev->next = h->pair->next->next; h->pair->prev->next = h->next->next;
                    h->pair->next->next->prev = h->prev; h->next->next->prev = h->pair->prev;
                    std::vector<Face*>::iterator itr_f = std::find(Faces.begin(), Faces.end(), h->pair->face);
                    h->face->ReConnect(h->prev);
                    //h->pair->vertex->halfedge.erase(std::find(h->pair->vertex->halfedge.begin(), h->pair->vertex->halfedge.end(), h->pair));
                    //h->vertex->halfedge.erase(std::find(h->vertex->halfedge.begin(), h->vertex->halfedge.end(), h));

                    Vertices.erase(std::find(Vertices.begin(), Vertices.end(), h->pair->vertex));
                    Vertices.erase(std::find(Vertices.begin(), Vertices.end(), h->vertex));

                    Edges.erase(itr_hn);
                    Edges.erase(itr_hpn);
                    Edges.erase(itr_hp);
                    Edges.erase(std::find(Edges.begin(), Edges.end(), h));
                    Faces.erase(itr_f);
                    T_crs.erase(T_crs.begin() + Ind);
                    break;
                }
            }
        }
    }while(Ind != -1);
}

Vertex* findAxisVertex(std::pair<CrvPt_FL*, Vertex*> R, Vertex *e1, Vertex *e2){
    for(auto& h: R.first->halfedge){
        if(h->next->vertex->p != R.second->p && h->next->vertex->p != e1->p && h->next->vertex->p != e2->p)return h->next->vertex;
    }
    return nullptr;
}

void FoldLine::drawRulingInAllAngles(std::vector<std::array<glm::f64vec3, 2>>& _Rulings){

    double a2 = 0;
    glm::f64vec3 e, e2, x;
    glm::f64vec3 N, r, befN;
    double veclen = 40;
    std::array<glm::f64vec3, 2> newVec;
    _Rulings.clear();
    double angle = 0.0;
    int flsize = FoldingCurve.size();
    while(angle <= 2*std::numbers::pi){
        double a = angle;

        double tau, dir;
        for(int i = 1; i < flsize - 1; i++){
            x = FoldingCurve[i].first->p3;
            e = (FoldingCurve[i-1].first->p3 - x)/glm::length((FoldingCurve[i-1].first->p3 - x));
            e2 = (FoldingCurve[i+1].first->p3 - x)/glm::length((FoldingCurve[i+1].first->p3 - x));

            double beta = std::acos(glm::dot(e,e2));
            if(i != 1){
                glm::f64vec3 srfN = MathTool::ProjectionVector(glm::cross(e, e2), -e, true);
                dir = glm::dot(-e, glm::cross(befN, srfN));
                tau = std::acos(glm::dot(srfN, befN));
                if(dir <= 0)tau = 2.0*std::numbers::pi - tau;
                a = (a2 - tau < 0) ? a2 - tau + 2.0 * std::numbers::pi: a2 - tau;
            }
            Vertex *AxisV = findAxisVertex(FoldingCurve[i], FoldingCurve[i-1].first, FoldingCurve[i+1].first);
            if(AxisV == nullptr)continue;
            glm::f64vec3 Axis = (AxisV->p - FoldingCurve[i].first->p)/glm::length((AxisV->p - FoldingCurve[i].first->p));
            double phi3 = std::acos(glm::dot((FoldingCurve[i+1].first->p - FoldingCurve[i].first->p)/glm::length((FoldingCurve[i+1].first->p - FoldingCurve[i].first->p)),
                    Axis));
            double phi4 = std::acos(glm::dot((FoldingCurve[i-1].first->p - FoldingCurve[i].first->p)/glm::length((FoldingCurve[i-1].first->p - FoldingCurve[i].first->p)),
                                    Axis));
            double k = 2.0 * std::numbers::pi - phi3 - phi4;
            double _phi1 = atan((std::cos(k) - std::cos(beta))/(std::sin(beta)*std::cos(a)- std::sin(k)));
            double phi1 = (_phi1 < 0)? _phi1+std::numbers::pi: _phi1;
            double phi2 = k - phi1;
            N = glm::normalize(glm::rotate (a, -e)  * glm::f64vec4{glm::cross(e2, e),1.0});
            r = glm::rotate(phi1, -N) * glm::f64vec4{e, 1.0};//新しいruling方向
            newVec = std::array<glm::f64vec3, 2>{veclen * r + x, x};_Rulings.push_back(newVec);

            double sin_a = sin(phi1)*sin(a)/sin(phi2);
            if(sin_a > 1)sin_a = 1;
            else if(sin_a < -1)sin_a = -1;
            double cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
            if(cos_a > 1)cos_a = 1;
            else if(cos_a < -1)cos_a = -1;
            a2 = (sin_a >= 0)? acos(cos_a): (sin_a < 0 && cos_a < 0)? 2.0*std::numbers::pi - acos(cos_a): 2.0*std::numbers::pi + asin(sin_a);
            if(i == (int)FoldingCurve.size() -2)break;
            befN = MathTool::ProjectionVector(glm::cross(e, e2), e2, true);
        }
        angle += 1e-3;
    }

}

void FoldLine::applyAAAMethod(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, double a, int type){
    if(FoldingCurve.empty())return;
    double phi02 = glm::angle(glm::normalize(FoldingCurve[0].second->p - FoldingCurve[0].first->p),glm::normalize(FoldingCurve[1].first->p - FoldingCurve[0].first->p));
    double phim1 = glm::angle(glm::normalize(FoldingCurve.end()[-2].first->p - FoldingCurve.back().first->p),glm::normalize(FoldingCurve.back().second->p - FoldingCurve.back().first->p));
    double ff, feq, fcon;
     _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a, phi02, phim1, ff, feq, fcon, type);
    //std::cout << a << " : "<<  Fbend(Poly_v, edges, Vertices) << " , " << Fruling(Poly_v) <<  std::endl;
}

void FoldLine::devide(Vertex *v1, Vertex *v2, std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, EdgeType _type){
    for(auto&f: Faces){
        HalfEdge *h = f->halfedge, *h1 = nullptr, *h2 = nullptr;
        do{
            if(h->vertex == v1)h1 = h;
            if(h->vertex == v2)h2 = h;
            h = h->next;
        }while(h != f->halfedge);
        if(h1 != nullptr && h2 != nullptr){
            HalfEdge *h1_new = new HalfEdge(h1->vertex, _type), *h2_new = new HalfEdge(h2->vertex, _type);
            h2_new->pair = h1_new; h1_new->pair = h2_new;

            h1_new->prev = h1->prev; h1_new->next = h2;
            h2_new->prev = h2->prev; h2_new->next = h1;
            h1->prev->next = h1_new; h2->prev->next = h2_new;
            h1->prev = h2_new; h2->prev = h1_new;

            HalfEdge *h = h1;
            f->halfedge = h1;
            do{
                h->face = f;
                h = h->next;
            }while(h != h1);
            h = h2;
            Face* face2 = new Face(h2);
            Faces.push_back(face2);
            face2->halfedge = h2;
            do{
                h2->face = face2;
                h = h->next;
            }while(h != h2);
            Edges.push_back(h1_new); Edges.push_back(h2_new);
            return;
        }
    }
}

void _FoldingAAAMethod(std::vector<std::pair<CrvPt_FL*, Vertex*>>& FoldingCurve, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Poly_V,
                                 double a, double phi02, double phim1, double& ff, double& feq, double& fcon, int type){

    auto calcTargetDistanceOnPlane = [](Vertex *p, Vertex *o, Vertex *v1, Vertex *v2){
        Eigen::Matrix2d A; Eigen::Vector2d b;
        b(0) = p->p.x - o->p.x; b(1) = p->p.y - o->p.y;
        A(0,0) = v1->p.x - o->p.x; A(0,1) = v2->p.x - o->p.x;
        A(1,0) = v1->p.y - o->p.y; A(1,1) = v2->p.y - o->p.y;
        Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
        return  x(0) * (v1->p3 - o->p3) + x(1) * (v2->p3 - o->p3) + o->p3;
    };

    auto getClosestVertex = [](Vertex *v, Vertex* o,  std::vector<HalfEdge*> Edges){
       Vertex *V_max = nullptr;
       double t_max = -1;
       for(auto&e: Edges){
           if((glm::length(e->vertex->p - v->p) < 1e-7|| glm::length(e->vertex->p - o->p) < 1e-7) || e->edgetype != EdgeType::ol)continue;
           if(MathTool::is_point_on_line(e->vertex->p, v->p, o->p)){
               double t = glm::distance(e->vertex->p, o->p)/glm::distance(v->p, o->p);
               if(t > t_max){
                   t_max = t; V_max = e->vertex;
               }
           }
       }
       return V_max;
   };

    ff = feq = fcon = 0.0;
    double a2 = 0, l;
    glm::f64vec3 e, e2, x, e_bef;
    glm::f64vec3 N, r, n, befN;
    //double veclen = 40;
    //std::array<glm::f64vec3, 2> newVec;
    //SingleRuling.clear();NewRuling2d.clear();
    int flsize = FoldingCurve.size();
    double tau, dir, _tau;
    for(int i = 1; i < flsize - 1; i++){
        x = FoldingCurve[i].first->p3;
        e = (FoldingCurve[i-1].first->p3 - x)/glm::length((FoldingCurve[i-1].first->p3 - x));
        e2 = (FoldingCurve[i+1].first->p3 - x)/glm::length((FoldingCurve[i+1].first->p3 - x));

        double beta = std::acos(glm::dot(e,e2));
        if(i != 1){
            glm::f64vec3 e_next = MathTool::ProjectionVector(e2, -e, true);
            glm::f64vec3 srfN = MathTool::ProjectionVector(glm::cross(e, e2), -e, true);
            //tau = abs(std::numbers::pi - std::acos(glm::dot(glm::normalize(glm::cross(e, e2)), befN)));
            //if(glm::dot(glm::cross(befN, glm::normalize(glm::cross(e, e2) ) ), e ) > 0)tau = 2.0 * std::numbers::pi - tau;
            dir = glm::dot(-e, glm::cross(befN, srfN));
            if(type == 0)tau = std::acos(glm::dot(e_bef, e_next));
            if(type == 1)tau = std::acos(glm::dot(srfN, befN));
            _tau = std::acos(glm::dot(srfN, befN));
            if(dir <= 0)_tau = 2.0*std::numbers::pi - _tau;
            //if(_tau < 0) _tau += 2.0 * std::numbers::pi;
            //if(tau > 2.0 * std::numbers::pi)tau -= 2.0*std::numbers::pi;
            tau = _tau;
            a = (a2 - tau < 0) ? a2 - tau + 2.0 * std::numbers::pi: a2 - tau;
        }
        Vertex *AxisV = findAxisVertex(FoldingCurve[i], FoldingCurve[i-1].first, FoldingCurve[i+1].first);
        if(AxisV == nullptr)continue;
        glm::f64vec3 Axis = (AxisV->p - FoldingCurve[i].first->p)/glm::length((AxisV->p - FoldingCurve[i].first->p));
        double phi3 = std::acos(glm::dot((FoldingCurve[i+1].first->p - FoldingCurve[i].first->p)/glm::length((FoldingCurve[i+1].first->p - FoldingCurve[i].first->p)),
                Axis));
        double tmp = glm::dot(glm::cross(glm::normalize(FoldingCurve[i+1].first->p - FoldingCurve[i].first->p), glm::normalize(FoldingCurve[i-1].first->p - FoldingCurve[i].first->p)), glm::f64vec3{0,0,1});
        double phi4 = std::acos(glm::dot((FoldingCurve[i-1].first->p - FoldingCurve[i].first->p)/glm::length((FoldingCurve[i-1].first->p - FoldingCurve[i].first->p)),
                                Axis));
        double k = 2.0 * std::numbers::pi - phi3 - phi4;
        double _phi1 = atan((std::cos(k) - std::cos(beta))/(std::sin(beta)*std::cos(a)- std::sin(k)));
        double phi1 = (_phi1 < 0)? _phi1+std::numbers::pi: _phi1;
        double phi2 = k - phi1;

        ff += abs(std::numbers::pi - phi1 - phi3);
        feq += abs(k/2.0 - phi1);
        fcon += abs(phi1 + phi4 - std::numbers::pi);
        // std::cout <<"angle "<< i << " :  " << "phi1 = " <<  std::setprecision(10) <<  glm::degrees(phi1) << ", phi2 = " << glm::degrees(phi2) << ", phi3 = " <<
         //        glm::degrees(phi3) << ", phi4 = " << glm::degrees(phi4)<< " : k = " <<  glm::degrees(k) << " , beta = " << glm::degrees(beta) << " , sum = " << phi1+phi2+phi3+phi4 - 2.0*std::numbers::pi<<  std::endl;
        N = glm::normalize(glm::rotate (a, -e)  * glm::f64vec4{glm::cross(e2, e),1.0});
        r = glm::rotate(phi1, -N) * glm::f64vec4{e, 1.0};r = glm::normalize(r);//新しいruling方向
        n = glm::rotate(phi1, glm::f64vec3{0,0,-1.0})* glm::f64vec4{(glm::normalize(FoldingCurve[i-1].first->p- FoldingCurve[i].first->p)), 1.0};n = glm::normalize(n);//展開図のruling方向
        //newVec = std::array<glm::f64vec3, 2>{veclen * e + x, x};SingleRuling.push_back(newVec);
        //newVec = std::array<glm::f64vec3, 2>{veclen * N + x, x};SingleRuling.push_back(newVec);
        //newVec = std::array<glm::f64vec3, 2>{veclen * r + x, x};SingleRuling.push_back(newVec);
        //newVec = std::array<glm::f64vec3, 2>{veclen * (e + e2) + x, x};SingleRuling.push_back(newVec);

        glm::f64vec3 crossPoint;
        bool hasPointOnEdge = IsRulingCrossed(n, FoldingCurve[i].first->p, crossPoint, Poly_V);
        if(hasPointOnEdge){
            l = glm::distance(crossPoint, FoldingCurve[i].first->p);
            FoldingCurve[i].second->p = crossPoint;
            FoldingCurve[i].second->p3 = l * r + FoldingCurve[i].first->p3;
        }else std::cout << "cross point could not be founded in " << i << std::endl;
        double sin_a = sin(phi1)*sin(a)/sin(phi2);
        if(sin_a > 1)sin_a = 1;
        else if(sin_a < -1)sin_a = -1;
        double cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
        if(cos_a > 1)cos_a = 1;
        else if(cos_a < -1)cos_a = -1;
        a2 = (sin_a >= 0)? acos(cos_a): (sin_a < 0 && cos_a < 0)? 2.0*std::numbers::pi - acos(cos_a): 2.0*std::numbers::pi + asin(sin_a);
        if(i == (int)FoldingCurve.size() -2)break;
        e_bef = MathTool::ProjectionVector(e, e2, true);
        befN = MathTool::ProjectionVector(glm::cross(e, e2), e2, true);
    }

    {//i = 0(端の面)

        Vertex *v_clst = getClosestVertex(FoldingCurve[0].second , FoldingCurve[0].first, Edges);
        if(v_clst == nullptr){
            x = FoldingCurve[0].first->p3;
            e = glm::normalize(FoldingCurve[1].first->p3 - x);
            glm::f64vec3 ee = glm::normalize(FoldingCurve[1].second->p3- FoldingCurve[1].first->p3);
            N = glm::normalize(glm::cross(ee,-e));
            r = glm::rotate(phi02, N) * glm::f64vec4{e,1};r = glm::normalize(r);//新しいruling方向
            l = glm::length(FoldingCurve[0].second->p - FoldingCurve[0].first->p);
            FoldingCurve[0].second->p3 = l * r + FoldingCurve[0].first->p3;
        }else{
            int j;
            for(j = 0; j < (int)FoldingCurve.size(); j++){
                if(FoldingCurve[j].second == v_clst)break;
            }
            if(j == (int)FoldingCurve.size())std::cout <<"can't find correct edge" << std::endl;
            else{
                auto p = calcTargetDistanceOnPlane(FoldingCurve[0].second, v_clst,  FoldingCurve[j].first, FoldingCurve[j+1].first);
                FoldingCurve[0].second->p3 = p;
            }
            //HalfEdge *h_clst = nullptr;
            //for(auto&h: v_clst->halfedge){if(h->edgetype == EdgeType::ol)h_clst = h;}
            //if(h_clst != nullptr){
                //auto p = calcTargetDistanceOnPlane(FoldingCurve[0].second, v_clst,  h_clst->prev->vertex, h_clst->prev->prev->vertex);
                //FoldingCurve[0].second->p3 = p;
            //}else {
                //std::cout <<"can't find correct edge" << std::endl;
            //}
        }
    }

    {
        Vertex *v_clst = getClosestVertex(FoldingCurve.back().second , FoldingCurve.back().first, Edges);
        if(v_clst == nullptr){
            x = FoldingCurve.back().first->p3;
            e = glm::normalize(FoldingCurve.end()[-2].first->p3 - x);
            glm::f64vec3 ee = glm::normalize(FoldingCurve.end()[-2].second->p3- FoldingCurve.end()[-2].first->p3);
            N = glm::normalize(glm::cross(ee,e));
            r = (glm::rotate(-phim1, N) * glm::f64vec4{e,1}); r = glm::normalize(r);//新しいruling方向
            l = glm::length(FoldingCurve.back().second->p - FoldingCurve.back().first->p);
            FoldingCurve.back().second->p3 = l * r + x;

            //newVec = std::array<glm::f64vec3, 2>{veclen * e + x, x};SingleRuling.push_back(newVec);
            //newVec = std::array<glm::f64vec3, 2>{veclen * N + x, x};SingleRuling.push_back(newVec);
            //newVec = std::array<glm::f64vec3, 2>{veclen * r + x, x};SingleRuling.push_back(newVec);
           // newVec = std::array<glm::f64vec3, 2>{veclen * e + x, x};SingleRuling.push_back(newVec);
        }else{
            HalfEdge *h_clst = nullptr;
            for(auto&h: v_clst->halfedge){if(h->edgetype != EdgeType::ol)h_clst = h;}
            if(h_clst != nullptr){
                auto p = calcTargetDistanceOnPlane(FoldingCurve.back().second, v_clst,  h_clst->next->vertex, h_clst->next->next->vertex);
                FoldingCurve.back().second->p3 = p;
            }else {
                std::cout <<"can't find correct edge" << std::endl;
            }
        }
    }

    //std::cout <<"||||||||||||||||||||||||||||||" << std::endl;
}

