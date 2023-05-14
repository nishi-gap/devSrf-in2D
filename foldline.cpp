
#include "foldline.h"

const double eps = 1e-5;
using namespace MathTool;

Vertex4d::Vertex4d(CrvPt_FL *v, Vertex *v2){
    first = v; second = v2; IsCalc = true;
}
Vertex4d::Vertex4d(const Vertex4d& V4d){
    first = V4d.first; second = V4d.second; IsCalc = V4d.IsCalc;
}

void _FoldingAAAMethod(std::vector<Vertex4d>& FoldingCurve, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Poly_V,
                                 double a, double phi02, double phim1, int type);
Vertex* findAxisVertex(Vertex4d& R, Vertex *e1, Vertex *e2);
std::vector<Vertex4d> TrimPoints2(std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces, std::vector<Vertex*>& Vertices, std::vector<Vertex4d>& FoldingCurve, double tol);
bool IsRulingCrossed(glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint,  std::vector<Vertex*>& Poly_V);
std::vector<CrvPt_FL*> SetPointsOnCurve(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_v, std::vector<glm::f64vec3>& CtrlPts, int dim, int divSize);
inline glm::f64vec3 calcCrossPoint_2Vector(Vertex* p1, Vertex* q1, Vertex* p2, Vertex* q2);
inline glm::f64vec3 calcTargetDistanceOnPlane(Vertex *p, Vertex *o, Vertex *v1, Vertex *v2);
void Douglas_Peucker_algorithm(std::vector<Vertex4d>& FoldingCurve, std::vector<Vertex4d>& res, double tol = std::numbers::pi/9.0);

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

//https://programming.pc-note.net/c/struct.html
typedef struct{
    std::vector<Vertex4d> FoldingCurve;
    std::vector<HalfEdge*> Edges;
    std::vector<Vertex*> Vertices, Poly_V;
    double phi02, phim1;
    int type;
}ObjData;//初期化の変数代入の順番が大事

double Fruling(const std::vector<double> &a, std::vector<double> &grad, void* f_data)
{
    auto f_r =[](std::vector<Vertex4d>& FoldingCurve)->double{
        double f = 0.0;
        Eigen::Matrix2d A;
        Eigen::Vector2d b;
        /*
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
        }*/
        for(int i = 1; i < (int)FoldingCurve.size() - 2; i++){

            A(0,0) = FoldingCurve[i].second->p.x - FoldingCurve[i].first->p.x; A(0,1) = -(FoldingCurve[i+1].second->p.x - FoldingCurve[i+1].first->p.x);
            A(1,0) = FoldingCurve[i].second->p.y - FoldingCurve[i].first->p.y; A(1,1) = -(FoldingCurve[i+1].second->p.y - FoldingCurve[i+1].first->p.y);
            b(0) = FoldingCurve[i+1].first->p.x - FoldingCurve[i].first->p.x;
            b(1) = FoldingCurve[i+1].first->p.y - FoldingCurve[i].first->p.y;
            Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
            if(0 < ts(0) && ts(0) < 1) f += 1.0/ts(0);
        }
        return f;
    };
    ObjData *od = (ObjData *)f_data;
    double h = 1e-7;
    double phi02 = od->phi02, phim1 = od->phim1;
    int type = od->type;
    std::vector<Vertex4d> FoldingCurve = od->FoldingCurve;
    double f = f_r(FoldingCurve);
    if(!grad.empty()){
        std::vector<HalfEdge*> Edges = od->Edges;
        std::vector<Vertex*> Poly_V = od->Poly_V;
        std::vector<Vertex> tmp;
        for(auto v: FoldingCurve)tmp.push_back(v.second);
        _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a[0] + h, phi02, phim1, type);
       double fp = f_r(FoldingCurve);
        _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a[0] - h, phi02, phim1, type);
       double fm = f_r(FoldingCurve);
       grad[0] = (fp - fm)/(2.0 * h);
       for(int j = 0; j < (int)FoldingCurve.size(); j++)FoldingCurve[j].second->p3 = tmp[j].p3;
    }
    std::cout <<"Fruling " << glm::degrees(a[0]) << ", " <<  f << std::endl;
    return f;
 }

double ObjFunc(const std::vector<double> &a, std::vector<double> &grad, void* f_data){
    double h = 1e-7;
    ObjData *od = (ObjData *)f_data;
    double phi02 = od->phi02, phim1 = od->phim1;
    int type = od->type;
    std::vector<Vertex4d> FoldingCurve = od->FoldingCurve;
    std::vector<HalfEdge*> Edges = od->Edges;
    std::vector<Vertex*> Vertices = od->Vertices, Poly_V = od->Poly_V;
    if(!grad.empty()){
        _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a[0] + h, phi02, phim1, type);
       double fp = Fbend(Poly_V, Edges, Vertices);
        _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a[0] - h, phi02, phim1, type);
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
        double fcon;
        CrossRegion(double a, bool _IsContinuation): first(a), second(-1), IsContinuation(_IsContinuation){}
    };

    double phi02, phim1;
    phi02 = glm::angle(glm::normalize(FoldingCurve[0].second->p - FoldingCurve[0].first->p),glm::normalize(FoldingCurve[1].first->p - FoldingCurve[0].first->p));
    phim1 = glm::angle(glm::normalize(FoldingCurve.end()[-2].first->p - FoldingCurve.back().first->p),glm::normalize(FoldingCurve.back().second->p - FoldingCurve.back().first->p));
    double a = 0, a2 = 0.0;
    /*
    std::ofstream ofs2;
    std::filesystem::create_directory("./Optimization");

    std::string AngleFile = "./Optimization/ChangeAngle.csv";
    ofs2.open(AngleFile, std::ios::out);
    ofs2 << "a, Fruling'"<<std::endl;
    */

    Vertex *AxisV = findAxisVertex(FoldingCurve[1], FoldingCurve[0].first, FoldingCurve[2].first);
    glm::f64vec3 Axis = (AxisV->p3 - FoldingCurve[1].first->p3)/glm::length((AxisV->p3 - FoldingCurve[1].first->p3));
    double phi3 = std::acos(glm::dot(glm::normalize(FoldingCurve[2].first->p3 - FoldingCurve[1].first->p3),Axis));
    double phi4 = std::acos(glm::dot(glm::normalize(FoldingCurve[0].first->p3 - FoldingCurve[1].first->p3),Axis));
    double k = 2.0 * std::numbers::pi - phi3 - phi4;
    double a_min, a_max;
    bool IsMount;
    for(const auto&e: Edges){
        if((e->vertex == FoldingCurve[1].first && e->next->vertex == FoldingCurve[1].second)){
            if(e->r->Gradation < 0)IsMount = false;
            else IsMount = true;
            break;
        }
    }
    double beta = std::acos(glm::dot(glm::normalize(FoldingCurve[2].first->p3 - FoldingCurve[1].first->p3), glm::normalize(FoldingCurve[0].first->p3 - FoldingCurve[1].first->p3)));
    double a_con = std::acos((sin(k) + (cos(beta) - cos(k))/tan(phi4))/sin(beta));
    std::cout << "continuatin " << a_con << " , degree = " << glm::degrees(a_con) << std::endl;
    if(k < std::numbers::pi && IsMount){a_min = a_con - std::numbers::pi;a_max = a_con - eps;}
    if(k >= std::numbers::pi && IsMount){a_min = a_con + eps; a_max = std::numbers::pi + a_con;}
    if(k < std::numbers::pi && !IsMount){a_min = a_con - std::numbers::pi;a_max = a_con - eps;}
    if(k >= std::numbers::pi && !IsMount){a_min = a_con + eps; a_max = std::numbers::pi + a_con;}

    ObjData od = {FoldingCurve, Edges, Vertices, Poly_V, phi02, phim1, type};
    nlopt::opt opt(nlopt::LD_MMA, 1);
    opt.set_lower_bounds(a_min);
    opt.set_upper_bounds(a_max);
    opt.set_min_objective(ObjFunc, &od);
    opt.add_inequality_constraint(Fruling, &od);
    opt.set_param("inner_maxeval", 300);
    opt.set_xtol_rel(1e-9);
    a = (a_min + a_max)/2.0;
    std::cout << "area " << glm::degrees(a_min) << " < " << glm::degrees(a) << " < " << glm::degrees(a_max) << std::endl;
    double minf;
      try {
        std::vector<double> _a{a};
          nlopt::result result = opt.optimize(_a, minf);
          std::cout <<"result :  " <<result << std::endl;
          std::cout << "found minimum at f(" << glm::degrees(_a[0]) << ") = "
              << std::setprecision(10) << minf << std::endl;
          a2 = _a[0];
      }
      catch (std::exception& e) {
          std::cout << "nlopt failed: " << e.what() << std::endl;
      }
    std::cout << "finish"<<std::endl;
    return true;
}

namespace RevisionVertices{
    struct EndEdge{
        glm::f64vec3 r;
        glm::f64vec3 o;
        bool empty;
        EndEdge(){empty = true;}
    };

    typedef struct{
        //std::vector<glm::f64vec3> R2d;//ruling
        //std::vector<glm::f64vec3> O2d;//始点
        std::vector<glm::f64vec3> Vertices;
        EndEdge Right;
        EndEdge Left;
    }data;
    inline double getK(const glm::f64vec3& o, const glm::f64vec3& x, const glm::f64vec3& x2){
        double k = std::acos(glm::dot(glm::normalize(x - o), glm::normalize(x2 - o)));
        if(glm::dot(glm::f64vec3{0,0,1}, glm::cross(glm::normalize(x - o), glm::normalize(x2 - o))) > 0)k = 2.0*std::numbers::pi - k;
        return k;
    }
    double ObjFunc(const std::vector<glm::f64vec3>& V){
        double f = 0.0;
        std::vector<double> K((int)V.size() - 2), K_avg((int)V.size() - 2);
        for(int i = 1; i < (int)V.size() - 1; i++)K[i - 1] = RevisionVertices::getK(V[i], V[i-1], V[i+1]);
        //移動平均(区間2)
        for(int i = 1; i < (int)K.size(); i++){
           K_avg[i] = (K[i] + K[i-1])/2.0;
        }
        if((int)K_avg.size() > 1)K_avg[0] = (K_avg[1] - K_avg[2]) + K_avg[1];
        for(int i = 0; i < (int)K.size(); i++)f += abs(K_avg[i] - K[i]);
        return f;
    }

    double Minimize(const std::vector<double> &T, std::vector<double> &grad, void* f_data){
        auto Fright = [](glm::f64vec3 v0, glm::f64vec3&v, glm::f64vec3 v2, glm::f64vec3& v3, glm::f64vec3& v4, double t, glm::f64vec3& R, glm::f64vec3& O){
          std::vector<double> F(2);
          double k0 = RevisionVertices::getK(v, v0, v2);
          double k1 = RevisionVertices::getK(v2, v, v3);
          double k2 = RevisionVertices::getK(v3, v2, v4);
          F[0] = abs(k0 - (2. * k1 - k2));
          glm::f64vec3 vp = (t + eps) * R + O;
          glm::f64vec3 vm = (t - eps) * R + O;
          double f0p = abs(RevisionVertices::getK(v, vp, v2) - (2.*k1 - k2)), f0m = abs(RevisionVertices::getK(v, vm, v2) - (2.*k1 - k2));
          F[1] = (f0p - f0m)/(2.0*eps);
          return F;
        };
        auto Fleft = [](glm::f64vec3 v0, glm::f64vec3&v, glm::f64vec3 v2, glm::f64vec3& v3, glm::f64vec3& v4, double t, glm::f64vec3& R, glm::f64vec3& O){
          std::vector<double> F(2);
          double k0 = RevisionVertices::getK(v, v2, v0);
          double k1 = RevisionVertices::getK(v2, v3, v);
          double k2 = RevisionVertices::getK(v3, v4, v2);
          F[0] = abs(k0 - (2. * k1 - k2));
          glm::f64vec3 vp = (t + eps) * R + O;
          glm::f64vec3 vm = (t - eps) * R + O;
          double f0p = abs(RevisionVertices::getK(v, v2, vp) - (2.*k1 - k2)), f0m = abs(RevisionVertices::getK(v, v2, vm) - (2.*k1 - k2));
          F[1] = (f0p - f0m)/(2.0*eps);
          return F;
        };
        data *od =  reinterpret_cast<data*>(f_data);
        //std::vector<glm::f64vec3> R2d = od->R2d;
        //std::vector<glm::f64vec3> O2d = od->O2d;
        std::vector<glm::f64vec3> Vertices = od->Vertices;
        EndEdge Left = od->Left, Right = od->Right;
        double f = 0.0;
        std::vector<double> f_end;
        if((int)grad.size() == 2){
            f_end = Fright(T[0] * Right.r + Right.o, Vertices[0], Vertices[1], Vertices[2], Vertices[3], T[0], Right.r, Right.o);
            f += f_end[0];  grad[0] = f_end[1];
            f_end = Fleft(T[1] * Left.r + Left.o, Vertices.back(), Vertices.end()[-2], Vertices.end()[-3], Vertices.end()[-4], T[1], Left.r, Left.o);
            f = f_end[0]; grad[0] = f_end[1];
        }else if((int)grad.size() == 1){
            if(!Left.empty)f_end = Fleft(T[0] * Left.r + Left.o, Vertices.back(), Vertices.end()[-2], Vertices.end()[-3], Vertices.end()[-4], T[0], Left.r, Left.o);
            if(!Right.empty)f_end = Fright(T[0] * Right.r + Right.o, Vertices[0], Vertices[1], Vertices[2], Vertices[3], T[0], Right.r, Right.o);
            f = f_end[0]; grad[0] = f_end[1];
        }
        std::cout <<"f = " << f << std::endl;
        return f;
    }
}

inline glm::f64vec3 calcCrossPoint_2Vector(Vertex* p1, Vertex* q1, Vertex* p2, Vertex* q2){
    Eigen::Matrix2d A; Eigen::Vector2d b;
    glm::f64vec3 v1 = p1->p - q1->p, v2 = p2->p - q2->p;
    b(0) = p2->p.x - p1->p.x; b(1) = p2->p.y - p1->p.y;
    A(0,0) = v1.x; A(0,1) = -v2.x;
    A(1,0) = v1.y; A(1,1) = -v2.y;
    Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
    return  x(0) * v1 + p1->p;
}

inline glm::f64vec3 calcTargetDistanceOnPlane(Vertex *p, Vertex *o, Vertex *v1, Vertex *v2){
    Eigen::Matrix2d A; Eigen::Vector2d b;
    b(0) = p->p.x - o->p.x; b(1) = p->p.y - o->p.y;
    A(0,0) = v1->p.x - o->p.x; A(0,1) = v2->p.x - o->p.x;
    A(1,0) = v1->p.y - o->p.y; A(1,1) = v2->p.y - o->p.y;
    Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
    return  x(0) * (v1->p3 - o->p3) + x(1) * (v2->p3 - o->p3) + o->p3;
}

void FoldLine::ApproximatePolyLine(){
    int ind_bef;
    double ay, by, cy;
    auto AnsQuadEq =[](double a, double b, double c) {
        std::vector<double> x;
       if(abs(a) > DBL_EPSILON){
           double _x1 = (-b - std::sqrt(b*b - 4.0*a*c))/(2.0*a), _x2 = (-b + std::sqrt(b*b - 4.0*a*c))/(2.0*a);
           x.push_back(_x1);
           if((_x2 - _x1) < 1e-3)x.insert(x.begin(), _x2);
           else if(_x1 - _x2 < 1e-3)x.push_back(_x2);
       }else x.push_back(-c/b);
       return x;
    };
    if(CtrlPts.size() != 4){std::cout <<"input curve is only cubic bezier curve " << std::endl; return;}
    ay = (CtrlPts[0].y - 3.*CtrlPts[1].y + 3.*CtrlPts[2].y - CtrlPts[3].y);
    by = (-2.*CtrlPts[0].y + 4.*CtrlPts[1].y - 2.*CtrlPts[2].y);
    cy = (CtrlPts[0].y - CtrlPts[1].y);
    std::vector<double>Ty = AnsQuadEq(ay, by, cy);
    if((int)Ty.size() == 1 || abs(Ty[0] - Ty[1]) < 1e-3)return;

    for(int i = 0; i < (int)Ty.size(); i++){
        double t_min = 10000;
        int ind = -1;
        for(auto&_v: FoldingCurve){
            if(abs(_v.first->s - Ty[i]) < t_min){
                ind = std::distance(FoldingCurve.begin(), std::find(FoldingCurve.begin(), FoldingCurve.end(), _v));
                t_min = abs(_v.first->s - Ty[i]);
            }
        }

        if(ind == -1)continue;
        if(i == 0){
            CrvPt_FL *p = FoldingCurve.back().first;
            for(int j = ind + 1; j < (int)FoldingCurve.size() - 1; j++){

                for(auto&h: FoldingCurve[j].second->halfedge){
                    if(h->edgetype == EdgeType::r){
                        FoldingCurve[j].first->p = calcCrossPoint_2Vector(std::get<0>(h->r->r), std::get<1>(h->r->r), p, FoldingCurve[ind].first);
                    }
                }
                double s = glm::length(FoldingCurve[j].first->p - p->p)/glm::length(p->p - FoldingCurve[ind].first->p);
                FoldingCurve[j].first->p3 = s * (FoldingCurve[ind].first->p3 - p->p3) + p->p3;
                for(auto&h: FoldingCurve[j].first->halfedge){
                    if(h->edgetype == EdgeType::fl)continue;
                    h->next->vertex->p3 = calcTargetDistanceOnPlane(h->next->vertex, FoldingCurve[j].first, FoldingCurve.back().second, FoldingCurve[ind].first);
                }
            }
            ind_bef = ind;
        }else if(i == (int)Ty.size() - 1){
            CrvPt_FL *p = FoldingCurve.front().first;
            for(int j = 1; j < ind; j++){
                for(auto&h: FoldingCurve[j].second->halfedge){
                    if(h->edgetype == EdgeType::r)FoldingCurve[j].first->p = calcCrossPoint_2Vector(std::get<0>(h->r->r), std::get<1>(h->r->r), p, FoldingCurve[ind].first);
                }
                double s = glm::length(FoldingCurve[j].first->p - p->p)/glm::length(p->p - FoldingCurve[ind].first->p);
                FoldingCurve[j].first->p3 = s * (FoldingCurve[ind].first->p3 - p->p3) + p->p3;
                for(auto&h: FoldingCurve[j].first->halfedge){
                    if(h->edgetype == EdgeType::fl)continue;
                    h->next->vertex->p3 = calcTargetDistanceOnPlane(h->next->vertex, FoldingCurve[j].first, FoldingCurve.front().second, FoldingCurve[ind].first);
                }
            }
            if(ind > ind_bef)std::swap(ind, ind_bef);
            for(int j = ind + 1; j < ind_bef; j++){
                for(auto&h: FoldingCurve[j].second->halfedge){
                    if(h->edgetype == EdgeType::r)FoldingCurve[j].first->p = calcCrossPoint_2Vector(std::get<0>(h->r->r), std::get<1>(h->r->r), FoldingCurve[ind_bef].first, FoldingCurve[ind].first);
                }
                double s = glm::length(FoldingCurve[j].first->p - FoldingCurve[ind_bef].first->p)/glm::length(FoldingCurve[ind_bef].first->p - FoldingCurve[ind].first->p);
                FoldingCurve[j].first->p3 = s * (FoldingCurve[ind].first->p3 - FoldingCurve[ind_bef].first->p3) + FoldingCurve[ind_bef].first->p3;
                for(auto&h: FoldingCurve[j].first->halfedge){
                    if(h->edgetype == EdgeType::fl)continue;
                    h->next->vertex->p3 = calcTargetDistanceOnPlane(h->next->vertex, FoldingCurve[j].first, FoldingCurve[ind_bef].second, FoldingCurve[ind].second);
                }
            }
        }else{

        }
    }
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
    if(type == PaintTool::FoldLine_arc){
        /*
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

        }*/
    }else if(type == PaintTool::FoldLine_bezier){
        for(auto& e: Edges){
            if((e->edgetype == EdgeType::r && (std::find(SearchedEdge.begin(), SearchedEdge.end(), e) != SearchedEdge.end() ||
                                               std::find(SearchedEdge.begin(), SearchedEdge.end(), e->pair) != SearchedEdge.end())))continue;
            SearchedEdge.push_back(e);
            std::vector<double>arcT = BezierClipping(CtrlPts, e, dim);

            for(auto&t: arcT){
                if(t < 0 || 1 < t){std::cout<<"t is not correct value " << t << std::endl; continue;}
                glm::f64vec3 v2{0,0,0};
                for (int i = 0; i < int(CtrlPts.size()); i++) v2 += MathTool::BernsteinBasisFunc(dim, i, t) * CtrlPts[i];
                if(!is_point_on_line(v2, e->vertex->p, e->next->vertex->p)) continue;
                double sa = glm::distance(v2, e->vertex->p), sc = glm::distance(e->vertex->p, e->next->vertex->p);
                glm::f64vec3 v3 = sa/sc * (e->next->vertex->p3 - e->vertex->p3) + e->vertex->p3;
                t_max = std::max(t_max, t); t_min = std::min(t_min, t);
                CrvPt_FL *P = new CrvPt_FL(v2, v3, t);
                P->set(v2, e->vertex, e->next->vertex);
                T_crs.push_back(P);
            }
        }
    }
    std::sort(T_crs.begin(), T_crs.end(), [](CrvPt_FL* T, CrvPt_FL* T2){
        return T->s > T2->s;
    });
    double tol = 1.0*std::numbers::pi/ 180.0;
    //std::cout << "before k" << T_crs.size() << std::endl;
    for(int i = 1; i < T_crs.size()-1; i++){
        double k = std::acos(glm::dot(glm::normalize(T_crs[i-1]->p - T_crs[i]->p), glm::normalize(T_crs[i+1]->p - T_crs[i]->p)));
        if(glm::dot(glm::f64vec3{0,0,1}, glm::cross(glm::normalize(T_crs[i-1]->p - T_crs[i]->p), glm::normalize(T_crs[i+1]->p - T_crs[i]->p))) > 0)k = 2.0*std::numbers::pi - k;
        //std::cout << k << std::endl;
    }

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
            if(h1->v->s > h2->v->s){FoldingCurve.push_back(Vertex4d(InsertPoints[1], he->prev->vertex));}
            else {FoldingCurve.push_back(Vertex4d(InsertPoints[0], h1->prev->vertex));}
        }
    }
    glm::f64vec3 befN = glm::normalize(FoldingCurve.front().second->p - FoldingCurve.front().first->p);


    for(int i = 0; i < (int)FoldingCurve.size(); i++){
        if(FoldingCurve[i].first == T_crs.front() || FoldingCurve[i].first == T_crs.back() )FoldingCurve.erase(FoldingCurve.begin() + i);
    }

    std::sort(FoldingCurve.begin(), FoldingCurve.end(), [](Vertex4d& F, Vertex4d& F2){
        return F.first->s > F2.first->s;
    });

    FoldingCurve.insert(FoldingCurve.begin(), Vertex4d(T_crs.front(), FindAnotherPoint(Poly_v, Edges, T_crs.front(), befN)));
    FoldingCurve.push_back(Vertex4d(T_crs.back(), FindAnotherPoint(Poly_v, Edges, T_crs.back(), befN)));
    Faces.insert(Faces.end(), Faces_new.begin(), Faces_new.end()); // 連結
    return true;
}

bool FoldLine::RevisionCrosPtsPosition(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, int type, bool TrimMode){
    if(FoldingCurve.empty())return false;

    //Y字のとき折り線は山であると仮定(谷の場合は色を反転)
    typedef struct {
        ruling *line;
        int type_mvk;//0:mountain && k < pi, 1: mountain && k >= pi, 2: vally && k < pi, 3: vally && k >= pi, -1: gradation = 0
        //0の数 + 3の数 == rulingの数 || 1の数 + 2の数 == rulingの数 -> 山谷とkの割り当てが正しい
        //それ以外
         //0の数 > 2の数 -> 　2の色を変換 (逆もまた然り)(1と3でも同様に)
    }MVK;

    //ApproximatePolyLine();
    //先にpiに近い値を削除してもう一度挑戦
    double tol = 0 * std::numbers::pi/180.0;
    auto tmp = TrimPoints2(Edges, Faces, Vertices, FoldingCurve, tol);
    for(auto&V4d: FoldingCurve){
        if(std::find(tmp.begin(), tmp.end(), V4d) != tmp.end()){         
        }else V4d.IsCalc = false;
    }

    std::vector<MVK> InitState;
    std::vector<int>MV(5,0);
    for(int i = 1; i < (int)tmp.size()-1; i++){
        double k = std::acos(glm::dot(glm::normalize(tmp[i-1].first->p - tmp[i].first->p), glm::normalize(tmp[i+1].first->p - tmp[i].first->p)));
        if(glm::dot(glm::f64vec3{0,0,1}, glm::cross(glm::normalize(tmp[i-1].first->p - tmp[i].first->p), glm::normalize(tmp[i+1].first->p - tmp[i].first->p))) > 0)k = 2.0*std::numbers::pi - k;
        for(auto&e: Edges){
            if(e->edgetype == EdgeType::r && MathTool::is_point_on_line(tmp[i].first->p, e->next->vertex->p, e->vertex->p)){
                int mv = (e->r->Gradation == 0)? -1: (e->r->Gradation > 0)? 0: 1;
                int type_mvk = (mv == 0 && k < std::numbers::pi)? 0: (mv == 0 && k >= std::numbers::pi)? 1: (mv == 1 && k < std::numbers::pi)? 2: (mv == 1 && k >= std::numbers::pi)? 3: -1;
                InitState.push_back(MVK(e->r, type_mvk));
                MV[type_mvk + 1] += 1;
                break;
            }
        }
    }
    if(MV[0] != 0){std::cout << "no folding pattern was found"<< std::endl; return false;}
    while(!(MV[1] + MV[4] == (int)tmp.size() - 2 || MV[2] + MV[3] == (int)tmp.size() - 2)){
        if(MV[1] == MV[2] && MV[3] == MV[4] && MV[1] == MV[4]){
        }else{
            for(auto& IS: InitState){
                if((MV[1] > MV[3] || (MV[1] == MV[3] && MV[2] < MV[4]))&& IS.type_mvk == 2){IS.line->Gradation *= -1; IS.type_mvk = 0; MV[1]++; MV[3]--; }
                else if((MV[1] < MV[3] || (MV[1] == MV[3] && MV[2] > MV[4])) && IS.type_mvk == 0){IS.line->Gradation *= -1; IS.type_mvk = 2; MV[1]--; MV[3]++;}
                else if((MV[2] > MV[4] || (MV[2] == MV[4] && MV[1] < MV[3])) && IS.type_mvk == 3){IS.line->Gradation *= -1; IS.type_mvk = 1; MV[2]++; MV[4]--;}
                else if((MV[2] < MV[4] || (MV[2] == MV[4] && MV[1] > MV[3])) && IS.type_mvk == 1){IS.line->Gradation *= -1; IS.type_mvk = 3; MV[2]--; MV[4]++;}
            }
        }
    }
    std::cout << "before revision" << std::endl;
    for(int i = 1; i < (int)tmp.size()-1; i++){
        double k = std::acos(glm::dot(glm::normalize(tmp[i-1].first->p - tmp[i].first->p), glm::normalize(tmp[i+1].first->p - tmp[i].first->p)));
        if(glm::dot(glm::f64vec3{0,0,1}, glm::cross(glm::normalize(tmp[i-1].first->p - tmp[i].first->p), glm::normalize(tmp[i+1].first->p - tmp[i].first->p))) > 0)k = 2.0*std::numbers::pi - k;
        std::cout << k << std::endl;
    }
    RevisionVertices::EndEdge Edge_fr, Edge_bc;
    int paramSize = 0;
    std::vector<double>lb, ub, T_rlen;
    double k0 = RevisionVertices::getK(tmp[1].first->p, tmp[0].first->p, tmp[2].first->p);
    double k1 = RevisionVertices::getK(tmp[2].first->p, tmp[1].first->p, tmp[3].first->p);
    double k2 = RevisionVertices::getK(tmp[3].first->p, tmp[2].first->p, tmp[4].first->p);
     std::cout <<"front " << abs(k0 - k1) << ", " << abs(k1 - k2) << std::endl;
    if(abs(k0 - k1) > abs(k1 - k2)){
        Edge_fr.r = FoldingCurve.front().first->ve->p - FoldingCurve.front().first->vo->p;
        Edge_fr.o = FoldingCurve.front().first->vo->p; Edge_fr.empty = false;
        T_rlen.push_back(FoldingCurve.front().first->rt);
        paramSize++; lb.push_back(0); ub.push_back(1);
    }
    k0 = RevisionVertices::getK(tmp.end()[-2].first->p, tmp.end()[-3].first->p, tmp.back().first->p);
    k1 = RevisionVertices::getK(tmp.end()[-3].first->p, tmp.end()[-4].first->p, tmp.end()[-2].first->p);
    k2 = RevisionVertices::getK(tmp.end()[-4].first->p, tmp.end()[-5].first->p, tmp.end()[-3].first->p);
    std::cout <<"back " << abs(k0 - k1) << ", " << abs(k1 - k2) << ", " << 2.*k1 - k2 <<  std::endl;
    if(abs(k0 - k1) > abs(k1 - k2)){
        Edge_bc.r = FoldingCurve.back().first->ve->p - FoldingCurve.back().first->vo->p;
        Edge_bc.o = FoldingCurve.back().first->vo->p; Edge_bc.empty = false;
        T_rlen.push_back(FoldingCurve.back().first->rt);
        paramSize++; lb.push_back(0); ub.push_back(1);
    }

    //std::vector<glm::f64vec3> T_R2d = {FoldingCurve.front().first->ve->p - FoldingCurve.front().first->vo->p, FoldingCurve.back().first->ve->p - FoldingCurve.back().first->vo->p};
    //std::vector<glm::f64vec3> T_O2d = {FoldingCurve.front().first->vo->p, FoldingCurve.back().first->vo->p};
    std::vector<glm::f64vec3> Vertices_CrosPt;
    for(int i = 1; i < (int)tmp.size() - 1; i++){
        Vertices_CrosPt.push_back(tmp[i].first->p);
    }

    RevisionVertices::data od = {Vertices_CrosPt, Edge_fr, Edge_bc};
    nlopt::opt opt(nlopt::LD_LBFGS, paramSize);

    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);
    opt.set_param("inner_maxeval", 100);
    opt.set_min_objective(RevisionVertices::Minimize, &od);
    opt.set_xtol_rel(1e-7);
    double minf;
      try {
          nlopt::result result = opt.optimize(T_rlen, minf);
          if(result == nlopt::SUCCESS){
              std::cout <<"revision finished  "<<std::endl;
              for(auto&t: T_rlen)std::cout <<t <<"  ";
              std::cout <<std::endl;
          }
      }
      catch (std::exception& e) {
          std::cout << "nlopt failed: " << e.what() << std::endl;
      }
    if(!Edge_fr.empty && !Edge_bc.empty){
        FoldingCurve.front().first->p = T_rlen[0] * Edge_fr.r + Edge_fr.o;
        FoldingCurve.front().first->p3 = T_rlen[0] * (FoldingCurve.front().first->ve->p_test - FoldingCurve.front().first->vo->p_test) + FoldingCurve.front().first->vo->p_test;
        FoldingCurve.back().first->p = T_rlen[1] * Edge_bc.r + Edge_bc.o;
        FoldingCurve.back().first->p3 = T_rlen[1] * (FoldingCurve.back().first->ve->p_test - FoldingCurve.back().first->vo->p_test) + FoldingCurve.back().first->vo->p_test;
    }else{
        if(!Edge_fr.empty){
            FoldingCurve.front().first->p = T_rlen[0] * Edge_fr.r + Edge_fr.o;
            FoldingCurve.front().first->p3 = T_rlen[0] * (FoldingCurve.front().first->ve->p_test - FoldingCurve.front().first->vo->p_test) + FoldingCurve.front().first->vo->p_test;
        }
        if(!Edge_bc.empty){
            FoldingCurve.back().first->p = T_rlen[0] * Edge_bc.r + Edge_bc.o;
            FoldingCurve.back().first->p3 = T_rlen[0] * (FoldingCurve.back().first->ve->p_test - FoldingCurve.back().first->vo->p_test) + FoldingCurve.back().first->vo->p_test;
        }
    }

    std::cout << "after revision" << std::endl;
    for(int i = 1; i < (int)tmp.size()-1; i++){
        double k = std::acos(glm::dot(glm::normalize(tmp[i-1].first->p - tmp[i].first->p), glm::normalize(tmp[i+1].first->p - tmp[i].first->p)));
        if(glm::dot(glm::f64vec3{0,0,1}, glm::cross(glm::normalize(tmp[i-1].first->p - tmp[i].first->p), glm::normalize(tmp[i+1].first->p - tmp[i].first->p))) > 0)k = 2.0*std::numbers::pi - k;
        std::cout << k << std::endl;
    }
    return Optimization(Edges, Vertices, Poly_V, type);
    return true;
}

std::vector<CrvPt_FL*> SetPointsOnCurve(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_v, std::vector<glm::f64vec3>& CtrlPts, int dim, int divSize){
    using namespace MathTool;

    std::vector<CrvPt_FL*>Points_On_Curve;
    if(divSize < 1 || Poly_v.empty() || (int)CtrlPts.size() <= dim)return Points_On_Curve;

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

    for(int i = 0; i < divSize; i++){
        double t = (double)(i + 1) * (t_max - t_min)/(double)(divSize + 1) + t_min;
        glm::f64vec3 v{0,0,0};
        for (int i = 0; i < int(CtrlPts.size()); i++) v += MathTool::BernsteinBasisFunc(dim, i, t) * CtrlPts[i];
        CrvPt_FL *vertex = new CrvPt_FL(v,v,t); //Vertices.push_back(vertex);
        Points_On_Curve.push_back(vertex);
    }
    std::sort(Points_On_Curve.begin(), Points_On_Curve.end(), [](CrvPt_FL* F, CrvPt_FL* F2){
        return F->s > F2->s;
    });

    glm::f64vec3 v{0,0,0};
    for (int i = 0; i < int(CtrlPts.size()); i++) v += MathTool::BernsteinBasisFunc(dim, i, t_max) * CtrlPts[i];
    CrvPt_FL *v_max = new CrvPt_FL(v, v, t_max);
    v = glm::f64vec3{0,0,0};
    for (int i = 0; i < int(CtrlPts.size()); i++) v += MathTool::BernsteinBasisFunc(dim, i, t_min) * CtrlPts[i];
    CrvPt_FL *v_min = new CrvPt_FL(v, v, t_min);
    Points_On_Curve.insert(Points_On_Curve.begin(), v_max);
    Points_On_Curve.push_back(v_min);
    return Points_On_Curve;
}

namespace OptimRulingDir{

    typedef struct{
        std::vector<HalfEdge*> Edges;
        std::vector<Vertex*> Vertices;
        std::vector<Vertex*> Poly_v;
        int divSize;
        int BezierDim;
        int BezierNum;

    }data;
    std::vector<glm::f64vec3> Projection2d(std::vector<glm::f64vec3>& Vec3d);

    void GenerateRulings(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_v, std::vector<glm::f64vec3>& CtrlPts,
                         std::vector<double>& Phi4, double a,int Bezierdim, int divSize){
        auto Points_On_Curve = SetPointsOnCurve(Edges, Vertices, Poly_v, CtrlPts, Bezierdim,  divSize);
        glm::f64vec3 befN;
        double tau, dir, a2;
        for(int i = 1; i < (int)Points_On_Curve.size() - 1; i++){

            glm::f64vec3 x = Points_On_Curve[i]->p3;
            glm::f64vec3 e = glm::normalize(Points_On_Curve[i-1]->p3 - x), e2 = glm::normalize(Points_On_Curve[i+1]->p3 - x);
            glm::f64vec3 e_2d = glm::normalize(Points_On_Curve[i-1]->p - Points_On_Curve[i]->p), e2_2d = glm::normalize(Points_On_Curve[i+1]->p - Points_On_Curve[i]->p);
            double beta = std::acos(glm::dot(e,e2));
            if(i != 1){
                glm::f64vec3 srfN = MathTool::ProjectionVector(glm::cross(e, e2), -e, true);
                dir = glm::dot(-e, glm::cross(befN, srfN));
                tau = std::acos(glm::dot(srfN, befN));
                if(dir <= 0)tau = 2.0*std::numbers::pi - tau;
                a = (a2 - tau <= 0) ? a2 - tau + 2.0 * std::numbers::pi: a2 - tau;
            }
            double k, phi4;
            glm::f64vec3 Nl, Nr;
            if(glm::dot(glm::cross(e_2d, e2_2d),glm::f64vec3{0,0,1}) < 0){
                k = 2.*std::numbers::pi - std::acos(glm::dot(e_2d, e2_2d));
                phi4 = 2.*std::numbers::pi - std::acos(glm::dot(e_2d, e2_2d));
                Nl = glm::normalize(glm::rotate (a, -e)  * glm::f64vec4{glm::cross(e2, e),1.0});
                Nr = glm::normalize(glm::rotate (-a, -e)  * glm::f64vec4{glm::cross(e2, e),1.0});
            }else{
                phi4 = std::acos(glm::dot(e_2d, e2_2d));
                k = std::acos(glm::dot(e_2d, e2_2d));
                Nl = glm::normalize(glm::rotate (a, -e)  * glm::f64vec4{glm::cross(e, e2),1.0});
                Nr = glm::normalize(glm::rotate (-a, -e)  * glm::f64vec4{glm::cross(e, e2),1.0});
            }
            phi4 /= 2;
            double phi1;
            if(abs(std::sin(beta)*std::cos(a)- std::sin(k)) < 1e-7) phi1 = std::numbers::pi/2.0;
            else{
                phi1 = atan2((std::cos(k) - std::cos(beta)),(std::sin(beta)*std::cos(a)- std::sin(k)));
                if(phi1 < 0)phi1 += std::numbers::pi;
            }
            double phi2 = k - phi1;

            glm::f64vec3 ruling_l2d = glm::normalize(glm::rotate(phi1, glm::f64vec3{0,0,-1.0})* glm::f64vec4{e_2d, 1.0});//展開図のruling方向
            glm::f64vec3 ruling_l3d = glm::normalize(glm::rotate(phi1, -Nl) * glm::f64vec4{e, 1.0});//新しいruling方向
            glm::f64vec3 ruling_r2d = glm::normalize(glm::rotate(phi4, glm::f64vec3{0,0,1}) * glm::f64vec4{e_2d, 1.0});
            glm::f64vec3 ruling_r3d = glm::normalize(glm::rotate(phi4, -Nr) * glm::f64vec4{e, 1.0});

            double sin_a = sin(phi1)*sin(a)/sin(phi2);
            if(sin_a > 1)sin_a = 1;
            else if(sin_a < -1)sin_a = -1;
            double cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
            if(cos_a > 1)cos_a = 1;
            else if(cos_a < -1)cos_a = -1;
            a2 = (sin_a >= 0)? acos(cos_a): (sin_a < 0 && cos_a < 0)? 2.0*std::numbers::pi - acos(cos_a): 2.0*std::numbers::pi + asin(sin_a);
            if(i == (int)Points_On_Curve.size() -2)break;
            befN = MathTool::ProjectionVector(glm::cross(e, e2), e2, true);
        }

    }

    double F_prj(std::vector<glm::f64vec3>& CtrlPts3d, std::vector<glm::f64vec3>& CtrlPts2d){
        double f = 0.0;
        if(CtrlPts2d.size() != CtrlPts3d.size()){
            std::cout <<"size is different"<<std::endl; return -1;
        }
        for(int k = 0; k < (int)CtrlPts2d.size() - 1; k++){
            glm::f64vec3 _CtrlPts3d = glm::f64vec3{CtrlPts3d[k].x, CtrlPts3d[k].y, 0};
            double delta = 1.0/3.0, sigma = 1e-2, d = glm::length(CtrlPts2d[k+1] - CtrlPts2d[k]);
            double w = std::exp(-d*d/(2.0*delta*delta)) + sigma;
            f += glm::length(_CtrlPts3d[k] - CtrlPts2d[k]) + w * (glm::length(glm::f64vec3{CtrlPts3d[k+1].x, CtrlPts3d[k+1].y, 0} - _CtrlPts3d) - glm::length(CtrlPts2d[k+1] - CtrlPts2d[k]));
        }
        return f;
    }
    double F_vari(){
        double f = 0.0;
        return f;
    }

    std::vector<glm::f64vec3> Projection2d(std::vector<glm::f64vec3>& Vec3d){
        std::vector<glm::f64vec3> Vec2d;
        for(auto&v: Vec3d)Vec2d.push_back(glm::f64vec3{v.x, v.y, 0});
        return Vec2d;
    }

    double Minimize(const std::vector<double> &X, std::vector<double> &grad, void* f_data){
        data *od =  reinterpret_cast<data*>(f_data);
        double f;
        std::vector<HalfEdge*> Edges = od->Edges;
        std::vector<Vertex*> Vertices = od->Vertices, Poly_v = od->Poly_v;
        int divSize = od->divSize, BezierDim = od->BezierDim, BezierNum = od->BezierNum;

        std::vector<double> Phi4(divSize);
        std::vector<glm::f64vec3> Bezier_ControlPts(BezierNum);

        double a = X.front();
        for(int i = 0; i < divSize; i++)Phi4[i] = X[i+1];
        for(int i = 0; i < BezierNum; i++)Bezier_ControlPts[i] = glm::f64vec3{X[i + divSize + 1], X[i + divSize + 2], X[i + divSize + 3]};

        GenerateRulings(Edges, Vertices, Poly_v, Bezier_ControlPts, Phi4, a, BezierDim, divSize);
        if(!grad.empty()){

        }
        return f;
    }
}

void FoldLine::Optimization2(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_v, double a){
    //if(a < DBL_EPSILON || abs(a - std::numbers::pi) < DBL_EPSILON || abs(a - 2.0*std::numbers::pi) < DBL_EPSILON)return;
    Points_On_Curve = SetPointsOnCurve(Edges, Vertices, Poly_v, CtrlPts, 3, 10);
    std::vector<double> X = {a};//最適化を行うパラメータ, a：折角度(1つ)、phi4(分割数)、B：ベジエ曲線の3次元空間上の制御点(4つ←今後4つ以上でもできるように)
    for(int i = 1; i < (int)Points_On_Curve.size() - 1; i++){
        glm::f64vec3 x = Points_On_Curve[i]->p;
        glm::f64vec3 e = glm::normalize(Points_On_Curve[i-1]->p - x);
        glm::f64vec3 e2 = glm::normalize(Points_On_Curve[i+1]->p - x);
        double phi = (glm::dot(glm::cross(e, e2), glm::f64vec3{0,0,1}) < 0)? std::acos(glm::dot(e, e2)): 2.0*std::numbers::pi - std::acos(glm::dot(e, e2));
        X.push_back(phi/2.0);
    }
    for(auto&P: CtrlPts){
        X.push_back(P.x);
        X.push_back(P.y);
        X.push_back(P.z);
    }
    int BezierDim = 3, BezierNum = 1;
    OptimRulingDir::data mydata = {Edges, Vertices, Poly_v, (int)Points_On_Curve.size() - 2, BezierDim, BezierNum};
    nlopt::opt opt(nlopt::LD_MMA, (int)X.size());
    opt.set_min_objective(OptimRulingDir::Minimize, &mydata);
    opt.set_xtol_rel(1e-9);

    double veclen = 60;
    glm::f64vec3 befN;
    NewRuling2d.clear(); AllRulings.clear();
    double tau, dir, a2;
    for(int i = 1; i < (int)Points_On_Curve.size() - 1; i++){
        glm::f64vec3 x = Points_On_Curve[i]->p3;
        glm::f64vec3 e = glm::normalize(Points_On_Curve[i-1]->p3 - x), e2 = glm::normalize(Points_On_Curve[i+1]->p3 - x);
        glm::f64vec3 e_2d = glm::normalize(Points_On_Curve[i-1]->p - Points_On_Curve[i]->p), e2_2d = glm::normalize(Points_On_Curve[i+1]->p - Points_On_Curve[i]->p);
        double beta = std::acos(glm::dot(e,e2));
        if(i != 1){
            glm::f64vec3 srfN = MathTool::ProjectionVector(glm::cross(e, e2), -e, true);
            dir = glm::dot(-e, glm::cross(befN, srfN));
            tau = std::acos(glm::dot(srfN, befN));
            if(dir <= 0)tau = 2.0*std::numbers::pi - tau;
            a = (a2 - tau <= 0) ? a2 - tau + 2.0 * std::numbers::pi: a2 - tau;
        }
        double k, phi4;
        glm::f64vec3 Nl, Nr;
        if(glm::dot(glm::cross(e_2d, e2_2d),glm::f64vec3{0,0,1}) < 0){
            k = 2.*std::numbers::pi - std::acos(glm::dot(e_2d, e2_2d));
            phi4 = 2.*std::numbers::pi - std::acos(glm::dot(e_2d, e2_2d));
            Nl = glm::normalize(glm::rotate (a, -e)  * glm::f64vec4{glm::cross(e2, e),1.0});
            Nr = glm::normalize(glm::rotate (-a, -e)  * glm::f64vec4{glm::cross(e2, e),1.0});
        }else{
            phi4 = std::acos(glm::dot(e_2d, e2_2d));
            k = std::acos(glm::dot(e_2d, e2_2d));
            Nl = glm::normalize(glm::rotate (a, -e)  * glm::f64vec4{glm::cross(e, e2),1.0});
            Nr = glm::normalize(glm::rotate (-a, -e)  * glm::f64vec4{glm::cross(e, e2),1.0});
        }
        phi4 /= 2;
        double phi1;
        if(abs(std::sin(beta)*std::cos(a)- std::sin(k)) < 1e-7) phi1 = std::numbers::pi/2.0;
        else{
            phi1 = atan2((std::cos(k) - std::cos(beta)),(std::sin(beta)*std::cos(a)- std::sin(k)));
            if(phi1 < 0)phi1 += std::numbers::pi;
        }
        double phi2 = k - phi1;

        glm::f64vec3 ruling_l2d = glm::normalize(glm::rotate(phi1, glm::f64vec3{0,0,-1.0})* glm::f64vec4{e_2d, 1.0});//展開図のruling方向
        glm::f64vec3 ruling_l3d = glm::normalize(glm::rotate(phi1, Nl) * glm::f64vec4{e, 1.0});//新しいruling方向
        glm::f64vec3 ruling_r2d = glm::normalize(glm::rotate(phi4, glm::f64vec3{0,0,1}) * glm::f64vec4{e_2d, 1.0});
        glm::f64vec3 ruling_r3d = glm::normalize(glm::rotate(phi4, Nr) * glm::f64vec4{e, 1.0});
        NewRuling2d.push_back({Points_On_Curve[i]->p, veclen * ruling_l2d + Points_On_Curve[i]->p});
        NewRuling2d.push_back({Points_On_Curve[i]->p, veclen * ruling_r2d + Points_On_Curve[i]->p});
        NewRuling2d.push_back({Points_On_Curve[i]->p, 20.0 * e_2d + Points_On_Curve[i]->p});
        std::cout << "i " << i << " , a = " << a << " , k = " << k  << " , phi1 = " << phi1 << ", phi4 = " << phi4 << std::endl;
        AllRulings.push_back({x, veclen * ruling_l3d + x});
        AllRulings.push_back({x, veclen * ruling_r3d + x});
        double sin_a = sin(phi1)*sin(a)/sin(phi2);
        if(sin_a > 1)sin_a = 1;
        else if(sin_a < -1)sin_a = -1;
        double cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
        if(cos_a > 1)cos_a = 1;
        else if(cos_a < -1)cos_a = -1;
        a2 = (sin_a >= 0)? acos(cos_a): (sin_a < 0 && cos_a < 0)? 2.0*std::numbers::pi - acos(cos_a): 2.0*std::numbers::pi + asin(sin_a);
        if(i == (int)Points_On_Curve.size() -2)break;
        befN = MathTool::ProjectionVector(glm::cross(e, e2), e2, true);
    }//std::cout <<"--------------"<<std::endl;

    double minf;
      try {
          //nlopt::result result = opt.optimize(X, minf);

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
            for(auto &h: Edges){
                if(MathTool::is_point_on_line(T_crs[Ind]->p, h->vertex->p, h->next->vertex->p)){

                    HalfEdge *hpn = h->pair->next, *hn = h->next, *hp = h->pair, *_h = h;
                    h->prev->next = h->pair->next->next; h->pair->prev->next = h->next->next;
                    h->pair->next->next->prev = h->prev; h->next->next->prev = h->pair->prev;
                    std::vector<Face*>::iterator itr_f = std::find(Faces.begin(), Faces.end(), h->pair->face);
                    h->face->ReConnect(h->prev);

                    //h->pair->vertex->halfedge.erase(std::find(h->pair->vertex->halfedge.begin(), h->pair->vertex->halfedge.end(), h->pair));
                    //h->vertex->halfedge.erase(std::find(h->vertex->halfedge.begin(), h->vertex->halfedge.end(), h));

                    Vertices.erase(std::find(Vertices.begin(), Vertices.end(), h->pair->vertex));
                    Vertices.erase(std::find(Vertices.begin(), Vertices.end(), h->vertex));

                    Edges.erase(std::find(Edges.begin(), Edges.end(), hn));
                    Edges.erase(std::find(Edges.begin(), Edges.end(), hpn));
                    Edges.erase(std::find(Edges.begin(), Edges.end(), hp));
                    Edges.erase(std::find(Edges.begin(), Edges.end(), _h));
                    Faces.erase(itr_f);
                    T_crs.erase(T_crs.begin() + Ind);
                    break;
                }
            }
        }
    }while(Ind != -1);
}

Vertex* findAxisVertex(Vertex4d& R, Vertex *e1, Vertex *e2){
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
            glm::f64vec3 Axis = (AxisV->p3 - FoldingCurve[i].first->p3)/glm::length((AxisV->p3 - FoldingCurve[i].first->p3));
            double phi3 = std::acos(glm::dot((FoldingCurve[i+1].first->p3 - FoldingCurve[i].first->p3)/glm::length((FoldingCurve[i+1].first->p3 - FoldingCurve[i].first->p3)),
                    Axis));
            double phi4 = std::acos(glm::dot((FoldingCurve[i-1].first->p3 - FoldingCurve[i].first->p3)/glm::length((FoldingCurve[i-1].first->p3 - FoldingCurve[i].first->p3)),
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
     _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a, phi02, phim1, type);
    //std::cout << a << " : "<<  Fbend(Poly_v, edges, Vertices) << " , " << Fruling(Poly_v) <<  std::endl;
}

void _FoldingAAAMethod(std::vector<Vertex4d>& FoldingCurve, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Poly_V,
                                 double a, double phi02, double phim1, int type){

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

    a = a - 2.0 * std::numbers::pi * (int)(a / (2.0 * std::numbers::pi));
    if(a < 0)a += 2.0 * std::numbers::pi;

    double a2 = 0, l;
    glm::f64vec3 e, e2, x, e_bef;
    glm::f64vec3 N, r, n, befN;
    //double veclen = 40;
    //std::array<glm::f64vec3, 2> newVec;
    //SingleRuling.clear();NewRuling2d.clear();
    std::vector<int> Vertices_Ind;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){
        if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);
    }

    double tau, dir, _tau;
    for(int ind = 1; ind < (int)Vertices_Ind.size() - 1; ind++){
        Vertex4d fc = FoldingCurve[Vertices_Ind[ind]];
        Vertex4d fc_bef = FoldingCurve[Vertices_Ind[ind - 1]];
        Vertex4d fc_next = FoldingCurve[Vertices_Ind[ind + 1]];
        x = fc.first->p3;
        e = (fc_bef.first->p3 - x)/glm::length((fc_bef.first->p3 - x));
        e2 = (fc_next.first->p3 - x)/glm::length((fc_next.first->p3 - x));

        double beta = std::acos(glm::dot(e,e2));
        if(ind != 1){
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
            a = (a2 - tau <= 0) ? a2 - tau + 2.0 * std::numbers::pi: a2 - tau;
        }
        Vertex *AxisV = findAxisVertex(fc, fc_bef.first, fc_next.first);
        if(AxisV == nullptr)continue;
        glm::f64vec3 Axis = (AxisV->p3 - fc.first->p3)/glm::length((AxisV->p3 - fc.first->p3));
        double phi3 = std::acos(glm::dot((fc_next.first->p3 - fc.first->p3)/glm::length((fc_next.first->p3 - fc.first->p3)),Axis));
        double phi4 = std::acos(glm::dot((fc_bef.first->p3 - fc.first->p3)/glm::length((fc_bef.first->p3 - fc.first->p3)), Axis));
        double k = 2.0 * std::numbers::pi - phi3 - phi4;
        double _phi1 = atan((std::cos(k) - std::cos(beta))/(std::sin(beta)*std::cos(a)- std::sin(k)));
        double phi1 = (_phi1 < 0)? _phi1+std::numbers::pi: _phi1;
        double phi2 = k - phi1;

        N = glm::normalize(glm::rotate (a, -e)  * glm::f64vec4{glm::cross(e2, e),1.0});
        r = glm::rotate(phi1, -N) * glm::f64vec4{e, 1.0};r = glm::normalize(r);//新しいruling方向
        n = glm::rotate(phi1, glm::f64vec3{0,0,-1.0})* glm::f64vec4{(glm::normalize(fc_bef.first->p- fc.first->p)), 1.0};n = glm::normalize(n);//展開図のruling方向

        glm::f64vec3 crossPoint;
        bool hasPointOnEdge = IsRulingCrossed(n, fc.first->p, crossPoint, Poly_V);
        if(hasPointOnEdge){
            l = glm::distance(crossPoint, fc.first->p);
            fc.second->p = crossPoint;
            fc.second->p3 = l * r + fc.first->p3;
        }else std::cout << "cross point could not be founded in " << Vertices_Ind[ind] << std::endl;
        double sin_a = sin(phi1)*sin(a)/sin(phi2);
        if(sin_a > 1)sin_a = 1;
        else if(sin_a < -1)sin_a = -1;
        double cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
        if(cos_a > 1)cos_a = 1;
        else if(cos_a < -1)cos_a = -1;
        a2 = (sin_a >= 0)? acos(cos_a): (sin_a < 0 && cos_a < 0)? 2.0*std::numbers::pi - acos(cos_a): 2.0*std::numbers::pi + asin(sin_a);
        if(ind == (int)FoldingCurve.size() -2)break;
        e_bef = MathTool::ProjectionVector(e, e2, true);
        befN = MathTool::ProjectionVector(glm::cross(e, e2), e2, true);
    }

    {//i = 0(端の面)

        Vertex *v_clst = getClosestVertex(FoldingCurve[Vertices_Ind[0]].second , FoldingCurve[Vertices_Ind[0]].first, Edges);
        if(v_clst == nullptr){
            x = FoldingCurve[Vertices_Ind[0]].first->p3;
            e = glm::normalize(FoldingCurve[Vertices_Ind[1]].first->p3 - x);
            glm::f64vec3 ee = glm::normalize(FoldingCurve[Vertices_Ind[1]].second->p3- FoldingCurve[Vertices_Ind[1]].first->p3);
            N = glm::normalize(glm::cross(ee,-e));
            r = glm::rotate(phi02, N) * glm::f64vec4{e,1};r = glm::normalize(r);//新しいruling方向
            l = glm::length(FoldingCurve[Vertices_Ind[0]].second->p - FoldingCurve[Vertices_Ind[0]].first->p);
            FoldingCurve[Vertices_Ind[0]].second->p3 = l * r + FoldingCurve[Vertices_Ind[0]].first->p3;
        }else{
            int j;
            for(j = 0; j < (int)Vertices_Ind.size(); j++){
                if(FoldingCurve[Vertices_Ind[j]].second == v_clst)break;
            }
            if(j == (int)Vertices_Ind.size())std::cout <<"can't find correct edge" << std::endl;
            else{
                auto p = calcTargetDistanceOnPlane(FoldingCurve[Vertices_Ind[0]].second, v_clst,  FoldingCurve[Vertices_Ind[j]].first, FoldingCurve[Vertices_Ind[j+1]].first);
                FoldingCurve[Vertices_Ind[0]].second->p3 = p;
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
        Vertex *v_clst = getClosestVertex(FoldingCurve[Vertices_Ind.back()].second , FoldingCurve[Vertices_Ind.back()].first, Edges);
        if(v_clst == nullptr){
            x = FoldingCurve[Vertices_Ind.back()].first->p3;
            e = glm::normalize(FoldingCurve[Vertices_Ind.end()[-2]].first->p3 - x);
            glm::f64vec3 ee = glm::normalize(FoldingCurve[Vertices_Ind.end()[-2]].second->p3 - FoldingCurve[Vertices_Ind.end()[-2]].first->p3);
            N = glm::normalize(glm::cross(ee,e));
            r = (glm::rotate(-phim1, N) * glm::f64vec4{e,1}); r = glm::normalize(r);//新しいruling方向
            l = glm::length(FoldingCurve[Vertices_Ind.back()].second->p - FoldingCurve[Vertices_Ind.back()].first->p);
            FoldingCurve[Vertices_Ind.back()].second->p3 = l * r + x;
        }else{
            HalfEdge *h_clst = nullptr;
            for(auto&h: v_clst->halfedge){if(h->edgetype != EdgeType::ol)h_clst = h;}
            if(h_clst != nullptr){
                auto p = calcTargetDistanceOnPlane(FoldingCurve[Vertices_Ind.back()].second, v_clst,  h_clst->next->vertex, h_clst->next->next->vertex);
                FoldingCurve[Vertices_Ind.back()].second->p3 = p;
            }else {
                std::cout <<"can't find correct edge" << std::endl;
            }
        }
    }

    for(int j = 0; j < (int)Vertices_Ind.size() - 1; j++){
        glm::f64vec3 Q = calcCrossPoint_2Vector(FoldingCurve[Vertices_Ind[j]].first, FoldingCurve[Vertices_Ind[j]].second, FoldingCurve[Vertices_Ind[j+1]].first, FoldingCurve[Vertices_Ind[j+1]].second);
        Vertex *Qv = new Vertex(Q);
        for(int i = Vertices_Ind[j] + 1; i < Vertices_Ind[j + 1]; i++){
            auto axis = findAxisVertex(FoldingCurve[i], FoldingCurve[i-1].first, FoldingCurve[i+1].first);
            if(glm::dot(glm::normalize(axis->p - FoldingCurve[i].first->p), glm::normalize(Q - FoldingCurve[i].first->p)) > 0)
                Qv->p = 1000. * (FoldingCurve[i].first->p - Q) + FoldingCurve[i].first->p;

            for(int k = 0; k < (int)Poly_V.size(); k++){
                glm::f64vec3 q = calcCrossPoint_2Vector(FoldingCurve[i].first, Qv, Poly_V[k], Poly_V[(k + 1) % (int)Poly_V.size()]);
                if(MathTool::is_point_on_line(q, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) &&
                        MathTool::is_point_on_line(q, FoldingCurve[i].first->p, Qv->p) ){

                    FoldingCurve[i].second->p = q;
                    double s = glm::length(FoldingCurve[i].first->p - FoldingCurve[Vertices_Ind[j]].first->p)/glm::length(FoldingCurve[Vertices_Ind[j]].first->p - FoldingCurve[Vertices_Ind[j+1]].first->p);
                    FoldingCurve[i].first->p3 = s * (FoldingCurve[Vertices_Ind[j+1]].first->p3 - FoldingCurve[Vertices_Ind[j]].first->p3) + FoldingCurve[Vertices_Ind[j]].first->p3;
                    for(auto&h: FoldingCurve[i].first->halfedge){
                        if(h->edgetype == EdgeType::fl)continue;
                        h->next->vertex->p3 = calcTargetDistanceOnPlane(h->next->vertex, FoldingCurve[i].first, FoldingCurve[Vertices_Ind[j]].first, FoldingCurve[Vertices_Ind[j+1]].first);
                    }
                }
            }
        }
    }

}

std::vector<Vertex4d> TrimPoints2(std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces, std::vector<Vertex*>& Vertices, std::vector<Vertex4d>& FoldingCurve, double tol){
    int Ind;
    double mink;
    std::vector<HalfEdge*>::iterator itr, itr_h;
    std::vector<Vertex*>::iterator itr_v;
    std::vector<Face*>::iterator itr_f;
    std::vector<Vertex4d> res = FoldingCurve;
    do{
       Ind = -1;
       mink = tol;
       for(int i = 1; i < (int)res.size() - 1; i++){
           double k = std::acos(glm::dot(glm::normalize(res[i-1].first->p - res[i].first->p), glm::normalize(res[i+1].first->p - res[i].first->p)));
           double k_diff = (std::numbers::pi - k);
           if(k_diff == std::min(mink, k_diff)){//最もpiに近くてtolよりも差が小さいkをもつ交点を選択
           Ind = i;
           mink = k_diff;
            }
        }
       if(Ind != -1){
           res.erase(res.begin() + Ind);
       }
       /*
        if(Ind != -1){
            for(auto&e: Edges){
                if(e->vertex != FoldingCurve[Ind].first) continue;
                std::vector<HalfEdge*> H = {e->vertex->halfedge.front(), nullptr, nullptr, nullptr};
                for(auto&h: e->vertex->halfedge){
                    if(h == H[0])continue;
                    if(h == H[0]->pair->next)H[1] = h;
                    else if(h == H[0]->prev->pair)H[3] = h;
                    else H[2] = h;
                }
                H[0]->prev->pair = H[2]->prev; H[2]->prev->pair = H[0]->prev;
                H[0]->prev->next = H[1]->next; H[1]->next->prev = H[0]->prev;
                H[0]->next->next->prev = H[1]->next->next; H[1]->next->next->next = H[0]->next->next;
                H[2]->prev->next = H[3]->next; H[3]->next->prev = H[2]->prev;
                H[3]->next->next->next = H[2]->next->next; H[2]->next->next->prev = H[3]->next->next;


                H[0]->prev->face->ReConnect(H[0]->prev);
                H[2]->prev->face->ReConnect(H[2]->prev);             
                itr = H.begin();
                Edges.erase(std::find(Edges.begin(), Edges.end(), H[0]->pair));
                Edges.erase(std::find(Edges.begin(), Edges.end(), H[2]->pair));
                Edges.erase(std::find(Edges.begin(), Edges.end(), H[0]->next));
                Edges.erase(std::find(Edges.begin(), Edges.end(), H[2]->next));
                while(itr != H.end()){
                    Edges.erase(std::find(Edges.begin(), Edges.end(), *itr));
                    itr = H.erase(H.begin());
                }

                itr_v = std::find(Vertices.begin(), Vertices.end(), FoldingCurve[Ind].first); if(itr_v != Vertices.end())Vertices.erase(itr_v);else std::cout <<"not found vertex"<<std::endl;
                FoldingCurve.erase(FoldingCurve.begin() + Ind);
                break;

            }
        }*/
    }while(Ind != -1);
    return res;
}


//https://issekinichou.wordpress.com/2022/12/28/douglas-peucker-algorithm/
void Douglas_Peucker_algorithm(std::vector<Vertex4d>& FoldingCurve, std::vector<Vertex4d>& res, double tol){
    auto perpendicularDistance = [](glm::f64vec3& p, glm::f64vec3& l_start, glm::f64vec3& l_end)->double{
        glm::f64vec3 line = (l_start - l_end);
        glm::f64vec3 OP = p - l_start;
        double d = glm::dot(OP, line)/glm::length(line);
        glm::f64vec3 H = l_start + d * line;
        return glm::distance(H, p);
    };
    if((int)FoldingCurve.size() < 2)return;
    double dmax = 0.0;
    size_t index = 0;
    size_t end = FoldingCurve.size()-1;
    for(size_t i = 1; i < end; i++){
        double d = perpendicularDistance(FoldingCurve[i].first->p, FoldingCurve[0].first->p, FoldingCurve.back().first->p);
        if (d > dmax){
            index = i; dmax = d;
        }
    }

    // If max distance is greater than epsilon, recursively simplify
    if(dmax > tol){
        // Recursive call
        std::vector<Vertex4d> res_first, res_last;
        std::vector<Vertex4d> firstLine{FoldingCurve.begin(), FoldingCurve.begin()+index+1};
        std::vector<Vertex4d> lastLine{FoldingCurve.begin() + index, FoldingCurve.end()};
        Douglas_Peucker_algorithm(firstLine,res_first, tol);
        Douglas_Peucker_algorithm(lastLine,res_last, tol);

        // Build the result list
        res.assign(res_first.begin(), res_first.end()-1);
        res.insert(res.end(), res_last.begin(), res_last.end());
        if(res.size()<2)return;
    } else {
        //Just return start and end points
        res.clear();
        res.push_back(FoldingCurve[0]);
        res.push_back(FoldingCurve.back());
    }
}


