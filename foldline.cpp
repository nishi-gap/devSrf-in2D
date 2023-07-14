
#include "foldline.h"

const double eps = 1e-7;
using namespace MathTool;

std::string File_Ebend = "./Optimization/Ebend.csv";
std::string File_Eruling = "./Optimization/Eruling.csv";
std::ofstream ofs_Ebend, ofs_Eruling;

inline glm::f64vec3 _calcruling3d(const double& a, glm::f64vec3 e, glm::f64vec3 e2, glm::f64vec3 axis, double& beta, std::vector<double>& Phi);
void CalcRuling(double a, Vertex4d& xbef, Vertex4d& x, Vertex4d& xnext, std::vector<Vertex*>& Poly_V, double& a2, glm::f64vec3& SrfN);
inline double update_flapangle(double a, const glm::f64vec3& befN, const glm::f64vec3& SrfN, const glm::f64vec3& e);
void _FoldingAAAMethod(std::vector<Vertex4d>& FoldingCurve, std::vector<Vertex*>& Poly_V,double a);
std::vector<Vertex4d> TrimPoints2(std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces, std::vector<Vertex*>& Vertices, std::vector<Vertex4d>& FoldingCurve, double tol);
bool IsRulingCrossed(glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint,  std::vector<Vertex*>& Poly_V);
std::vector<CrvPt_FL*> SetPointsOnCurve(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_v, std::vector<glm::f64vec3>& CtrlPts, int dim, int divSize);
void Douglas_Peucker_algorithm(std::vector<Vertex4d>& FoldingCurve, std::vector<Vertex4d>& res, double tol = std::numbers::pi/9.0);
void Douglas_Peucker_algorithm2(std::vector<Vertex4d>& FoldingCurve, std::vector<Vertex4d>& res, int size);
Vertex* getClosestVertex(Vertex *v, Vertex* o,  const std::vector<Vertex4d>& FoldingCurve, bool SkipTrimedPoint = true);

inline glm::f64vec3 calcCrossPoint_2Vector(Vertex* p1, Vertex* q1, Vertex* p2, Vertex* q2){
    Eigen::Matrix2d A; Eigen::Vector2d b;
    glm::f64vec3 v1 = q1->p - p1->p, v2 = q2->p - p2->p;
    b(0) = p2->p.x - p1->p.x; b(1) = p2->p.y - p1->p.y;
    A(0,0) = v1.x; A(0,1) = -v2.x;
    A(1,0) = v1.y; A(1,1) = -v2.y;
    double t = ((p2->p.x - p1->p.x)*(p2->p.y - q2->p.y) - (p2->p.x - q2->p.x)*(p2->p.y - p1->p.y))/((q1->p.x - p1->p.x) * (p2->p.y - q2->p.y) - (p2->p.x - q2->p.x)*(q1->p.y - p1->p.y));
    return glm::f64vec3{t * (q1->p.x - p1->p.x) + p1->p.x, t * (q1->p.y - p1->p.y) + p1->p.y, 0};
    Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
    return  x(0) * v1 + p1->p;
}
inline glm::f64vec3 calcCrossPoint_2Vector(glm::f64vec3 p1, glm::f64vec3 q1, glm::f64vec3 p2, glm::f64vec3 q2){

    double t = ((p2.x - p1.x)*(p2.y - q2.y) - (p2.x - q2.x)*(p2.y - p1.y))/((q1.x - p1.x) * (p2.y - q2.y) - (p2.x - q2.x)*(q1.y - p1.y));
    return glm::f64vec3{t * (q1.x - p1.x) + p1.x, t * (q1.y - p1.y) + p1.y, 0};

    Eigen::Matrix2d A; Eigen::Vector2d b;
    b(0) = p2.x - p1.x; b(1) = p2.y - p1.y;
    A(0,0) = (q1 - p1).x; A(0,1) = -(q2 - p2).x;
    A(1,0) = (q1 - p1).y; A(1,1) = -(q2 - p2).y;
    Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
    return  x(0) * (q1 - p1) + p1;
}
inline bool IsParallel(Vertex* p1, Vertex* q1, Vertex* p2, Vertex* q2){
    glm::f64vec3 v1 = glm::normalize(q1->p - p1->p), v2 = glm::normalize(q2->p - p2->p);
    if(abs(glm::dot(v1, v2)) >= 1.0 - 1e-5)return true;
    return false;
}

inline glm::f64vec3 calcTargetDistanceOnPlane(glm::f64vec3 p, Vertex *o, Vertex *v1, Vertex *v2){
    Eigen::Matrix2d A; Eigen::Vector2d b;
    b(0) = p.x - o->p.x; b(1) = p.y - o->p.y;
    A(0,0) = v1->p.x - o->p.x; A(0,1) = v2->p.x - o->p.x;
    A(1,0) = v1->p.y - o->p.y; A(1,1) = v2->p.y - o->p.y;
    Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
    return  x(0) * (v1->p3 - o->p3) + x(1) * (v2->p3 - o->p3) + o->p3;
}

namespace RevisionVertices{
    using FoldLine3d = std::vector<Vertex4d>;
    struct OptimizeParam{
        FoldLine3d FC;
        std::vector<Vertex*> Vertices, Poly_V;
        std::vector<double*> res_Fbend, res_Fruling, res_a;
        double wb, wp;
        void AddWeight(double _wb, double _wp){ wb = _wb; wp = _wp;}

        OptimizeParam(FoldLine3d& _FC, std::vector<Vertex*>& _Vertices, std::vector<Vertex*>& _Poly_V):FC{_FC}, Vertices{_Vertices}, Poly_V{_Poly_V}{}
        ~OptimizeParam(){}
    private:

    };
    struct OptimizeParam_v: public OptimizeParam{
        double a;
        OptimizeParam_v(double _a, FoldLine3d& _FC,  std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V):
            a{_a}, OptimizeParam::OptimizeParam( _FC, Vertices,  Poly_V){}
        ~OptimizeParam_v(){}
    };
    struct SmoothingArea{
        Vertex4d stP, lastP;
        std::vector<Vertex*> OriginalVertices;
        int st_ind, last_ind;
        Vertex *qt, *qb;
        SmoothingArea(Vertex4d& a, Vertex4d& _end, int si, int li, std::vector<Vertex*>& OV): stP{a}, lastP{_end}, st_ind{si}, last_ind{li}, OriginalVertices{OV}{
            auto getV = [](Vertex *o, Vertex *x, Vertex *o2, Vertex *x2)->Vertex*{
                if(IsParallel(o, x, o2, x2))return nullptr;
                glm::f64vec3 p2d = calcCrossPoint_2Vector(o, x, o2, x2);
                glm::f64vec3 p3d = calcTargetDistanceOnPlane(p2d, o,  x, x2);
                return  new Vertex(p2d, p3d);
            };
            qt = getV(stP.first, stP.second, lastP.first, lastP.second);
            qb = getV(stP.first, stP.third, lastP.first, lastP.third);
        }
        double Edev(std::vector<glm::f64vec3>& P, bool IsConnectEndPoint, const double th = 1e-5){
            auto Angle = [](glm::f64vec3& e, glm::f64vec3& e2)->double{return (glm::dot(e, e2) >= 1)? 0: (glm::dot(e, e2) <= -1)? std::numbers::pi: std::acos(glm::dot(e, e2));  };
            double f = 0.0;
            glm::f64vec3 el, er, et, eb;
            for(int i = 0; i < (int)P.size(); i++){
                if(i == 0)er = glm::normalize(stP.first->p3 - P[0]); else er = glm::normalize(P[i-1] - P[i]);
                if(i == (int)P.size() - 1 || IsConnectEndPoint) el = glm::normalize(lastP.first->p3 - P[i]);
                 else  el = glm::normalize(P[i+1] - P[i]);
                if(qt != nullptr)et = glm::normalize(qt->p3 - P[i]);
                else {
                    et = (st_ind != 0)? glm::normalize(stP.second->p3 - stP.first->p3): glm::normalize(lastP.second->p3 - lastP.first->p3);
                }
                if(qb != nullptr)eb = glm::normalize(qb->p3 - P[i]);
                else {
                    eb = (st_ind != 0)? glm::normalize(stP.third->p3 - stP.first->p3): glm::normalize(lastP.third->p3 - lastP.first->p3);
                }
                f += (std::abs(2.0 * std::numbers::pi - (Angle(el, et) + Angle(et, er) + Angle(er, eb) + Angle(eb, el))) - th);
            }
            return f;
        }

        double E_density(const std::vector<glm::f64vec3>& Qt, const std::vector<glm::f64vec3>& Qb){
            double f = .0;
            return f;
        }

        ~SmoothingArea(){}
    };
    struct SmoothSurface{
        std::vector<SmoothingArea> SA;
        bool IsConnectEndPoint;
    };

    using ObjData = OptimizeParam;
    using ObjData_v = OptimizeParam_v;
    using ObjData_smooth = SmoothSurface;

    inline double getK(const glm::f64vec3 o, const glm::f64vec3 x, const glm::f64vec3 x2);
    inline glm::f64vec3 getVertex(const glm::f64vec3& o, const glm::f64vec3& p, const double t);
    double F_vpos(FoldLine3d &FC, const std::vector<double>& T, std::vector<double>& grad);
    double F_k(FoldLine3d &FC, const std::vector<double>& T, std::vector<double>& grad);
    double F_ruling(const std::vector<double> &T, std::vector<double> &grad, void* f_data);

    double Const_Edev(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    double E_fair(std::vector<std::vector<glm::f64vec3>>& P, std::vector<double>& grad);
    double E_sim(std::vector<std::vector<glm::f64vec3>>& P, std::vector<std::vector<glm::f64vec3>>& Pori, std::vector<double>& grad);
    double E_iso(std::vector<std::vector<glm::f64vec3>>& P, std::vector<std::vector<glm::f64vec3>>& Pori, std::vector<double>& grad);
    double E_planerity(const std::vector<Vertex4d>& FC, bool IsLeft4pt, bool IsRight4pt);


    double Minimize_SmoothSrf(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    double Minimize_Vpos(const std::vector<double> &T, std::vector<double> &grad, void* f_data);
}

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
    //he = he2 = nullptr;
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

double Fparallel(std::vector<Vertex4d>& FoldingCurve){
    std::vector<int>Vertices_Ind;
    double f = 0.0;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}
    for(int i = 1; i < (int)Vertices_Ind.size() - 2; i++){
        f += glm::dot(glm::normalize(FoldingCurve[i].second->p - FoldingCurve[i].first->p), glm::normalize(FoldingCurve[i+1].second->p - FoldingCurve[i+1].first->p));
    }
    return f;
}

double Fbend2(std::vector<Vertex4d>& FoldingCurve){
    std::vector<int> Vertices_Ind;
    double f = 0.0;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}
    glm::f64vec3 Nt = glm::normalize(glm::cross(FoldingCurve[0].first->p3 - FoldingCurve[Vertices_Ind[1]].first->p3, FoldingCurve[Vertices_Ind[1]].second->p3 - FoldingCurve[Vertices_Ind[1]].first->p3));
    glm::f64vec3 Nb = glm::normalize(glm::cross(FoldingCurve[Vertices_Ind[1]].third->p3 - FoldingCurve[Vertices_Ind[1]].first->p3, FoldingCurve[0].first->p3 - FoldingCurve[Vertices_Ind[1]].first->p3));
    double phi = ((glm::dot(Nt,Nb)) > 1)? std::numbers::pi: ((glm::dot(Nt,Nb)) < -1)? 0:  std::numbers::pi - abs(std::acos(glm::dot(Nt,Nb)));
    f += 1.0/(phi*phi);
    for(int i = 1; i < (int)Vertices_Ind.size() - 1; i++){
        glm::f64vec3 Ntp = glm::normalize(glm::cross(FoldingCurve[Vertices_Ind[i]].first->p3 - FoldingCurve[Vertices_Ind[i+1]].first->p3,
                FoldingCurve[Vertices_Ind[i+1]].second->p3 - FoldingCurve[Vertices_Ind[i+1]].first->p3));
        glm::f64vec3 Nbp = glm::normalize(glm::cross(FoldingCurve[Vertices_Ind[i+1]].third->p3 - FoldingCurve[Vertices_Ind[i+1]].first->p3,
                FoldingCurve[Vertices_Ind[i]].first->p3 - FoldingCurve[Vertices_Ind[i+1]].first->p3));

        phi = ((glm::dot(Ntp,Nbp)) > 1)? std::numbers::pi: ((glm::dot(Ntp,Nbp)) < -1)? 0:  std::numbers::pi - abs(std::acos(glm::dot(Ntp,Nbp)));
        f += 1.0/(phi*phi);
        phi = ((glm::dot(Ntp,Nt)) > 1)? std::numbers::pi: ((glm::dot(Ntp,Nt)) < -1)? 0:  std::numbers::pi - std::acos(glm::dot(Ntp,Nt));
        f += 1.0/(phi*phi);
        Nt = Ntp; Nb = Nbp;
    }
    return f;
}

double RulingsCrossed(std::vector<Vertex4d>& FoldingCurve){
    double f = 0.0;
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    std::vector<int> Vertices_Ind;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}
    for(int i = 1; i < (int)Vertices_Ind.size() -1; i++){
        for(int j = -1; j <= 1; j +=2){
            int i2 = i + j;
            if(i2 <= 0 || i2 >= (int)Vertices_Ind.size() - 1)continue;
            glm::f64vec3 p1 = FoldingCurve[Vertices_Ind[i]].first->p, q1 = FoldingCurve[Vertices_Ind[i]].second->p, p2 = FoldingCurve[Vertices_Ind[i2]].first->p, q2 = FoldingCurve[Vertices_Ind[i2]].second->p;
            double t = ((p2.x - p1.x)*(p2.y - q2.y) - (p2.x - q2.x)*(p2.y - p1.y))/((q1.x - p1.x) * (p2.y - q2.y) - (p2.x - q2.x)*(q1.y - p1.y));
            if(0 < t && t < 1) f += 1.0/t;
            continue;
            A(0,0) = FoldingCurve[Vertices_Ind[i]].first->p.x - FoldingCurve[Vertices_Ind[i]].second->p.x;
            A(0,1) = -(FoldingCurve[Vertices_Ind[i2]].first->p.x - FoldingCurve[Vertices_Ind[i2]].second->p.x);
            A(1,0) = FoldingCurve[Vertices_Ind[i]].first->p.y - FoldingCurve[Vertices_Ind[i]].second->p.y;
            A(1,1) = -(FoldingCurve[Vertices_Ind[i2]].first->p.y - FoldingCurve[Vertices_Ind[i2]].second->p.y);
            b(0) = FoldingCurve[Vertices_Ind[i]].first->p.x - FoldingCurve[Vertices_Ind[i2]].first->p.x;
            b(1) = FoldingCurve[Vertices_Ind[i]].first->p.y - FoldingCurve[Vertices_Ind[i2]].first->p.y;
            Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
            if(0 < ts(0) && ts(0) < 1) f += 1.0/ts(0);
        }
    }
    return f;
}

double Fruling(const std::vector<double> &a, std::vector<double> &grad, void* f_data)
{
    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    RevisionVertices::FoldLine3d FoldingCurve = od->FC;
    double f = RulingsCrossed(FoldingCurve);
    if(!grad.empty()){
        std::vector<Vertex*> Poly_V = od->Poly_V;
        std::vector<Vertex> tmp;
        for(auto v: FoldingCurve)tmp.push_back(v.second);
        _FoldingAAAMethod(FoldingCurve, Poly_V, a[0] + eps);
       double fp = RulingsCrossed(FoldingCurve);
        _FoldingAAAMethod(FoldingCurve, Poly_V, a[0] - eps);
       double fm = RulingsCrossed(FoldingCurve);
       grad[0] = (fp - fm)/(2.0 * eps);
       //for(int j = 0; j < (int)FoldingCurve.size(); j++)FoldingCurve[j].second->p3 = tmp[j].p3;
    }
    if(ofs_Eruling.is_open()) ofs_Eruling << a[0] <<", " << glm::degrees(a[0]) << ", " <<  f <<  " ,  " << grad[0] << std::endl;
    return f;
 }

double ObjFunc_FlapAngle(const std::vector<double> &a, std::vector<double> &grad, void* f_data){
    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    RevisionVertices::FoldLine3d FoldingCurve = od->FC;
    std::vector<Vertex*>Poly_V = od->Poly_V;
    _FoldingAAAMethod(FoldingCurve, Poly_V, a[0]);
    double fb = od->wb * Fbend2(FoldingCurve), fp = od->wp * Fparallel(FoldingCurve);
    if(!grad.empty()){
        _FoldingAAAMethod(FoldingCurve, Poly_V, a[0] + eps);
       double fp =  od->wb * Fbend2(FoldingCurve) + od->wp * Fparallel(FoldingCurve);
        _FoldingAAAMethod(FoldingCurve, Poly_V, a[0] - eps);
       double fm =  od->wb * Fbend2(FoldingCurve) + od->wp * Fparallel(FoldingCurve);
       grad[0] = (fp - fm)/(2.0 * eps);
    }
    //std::cout <<"Fbend(" << glm::degrees(a[0]) << ")  =  " <<  fb <<  " ,  " << grad[0] << "  ,  Fparalell = " << fp <<  std::endl;
    if(ofs_Ebend.is_open())ofs_Ebend << a[0] <<", " << glm::degrees(a[0]) << ", " <<  fb <<  " ,  " << grad[0] << std::endl;
    return fb + fp;
}

double RevisionVertices::Const_Edev(const std::vector<double>& X, std::vector<double> &grad, void *f_data){
    SmoothSurface *od = (SmoothSurface*)f_data;
    std::vector<SmoothingArea> SA = od->SA;
    bool ICE = od->IsConnectEndPoint;
    double th = 1e-6;
    int ind = 0;
    std::vector<std::vector<glm::f64vec3>> P;
    for(auto&sa: SA){
        std::vector<glm::f64vec3> _P;
        for(int n = 0; n < (int)sa.OriginalVertices.size(); n++){
            _P.push_back(glm::f64vec3{X[3 * ind], X[3 * ind + 1], X[3 * ind + 2]});
            ind++;
        }
        P.push_back(_P);
    }
    double f = 0.0;
    ind = 0;
    if(!grad.empty()){
        for(int i = 0; i < (int)SA.size(); i++){
            f += SA[i].Edev(P[i], ICE, th);
            double fp, fm;
            for(int n = 0; n < (int)P[i].size(); n++){
                P[i][n].x += eps; fp = SA[i].Edev(P[i], ICE, th); P[i][n].x -= 2.0*eps; fm = SA[i].Edev(P[i], ICE, th);
                grad[ind] = (fp - fm)/(2.0*eps); P[i][n].x = X[ind++];
                P[i][n].y += eps; fp = SA[i].Edev(P[i], ICE, th); P[i][n].y -= 2.0*eps; fm = SA[i].Edev(P[i], ICE, th);
                grad[ind] = (fp - fm)/(2.0*eps); P[i][n].y = X[ind++];
                P[i][n].z += eps; fp = SA[i].Edev(P[i], ICE, th); P[i][n].z -= 2.0*eps; fm = SA[i].Edev(P[i], ICE, th);
                grad[ind] = (fp - fm)/(2.0*eps); P[i][n].z = X[ind++];
            }
        }
    }
    std::cout<<"developability  = " << f << std::endl;
    return f;
}

double RevisionVertices::E_planerity(const std::vector<Vertex4d>& FC, bool IsLeft4pt, bool IsRight4pt){
    double f = .0;
    auto CalcPlanarity = [](glm::f64vec3& p, glm::f64vec3& q, glm::f64vec3& p2, glm::f64vec3& q2){
        double l_avg = (glm::distance(p, q2) + glm::distance(q, p2))/2.0;
        glm::f64vec3 u1 = glm::normalize(p - q2), u2 = glm::normalize(q -  p2);
        if(glm::length(glm::cross(u1, u2)) < DBL_EPSILON)return glm::distance(p2 + glm::dot(q - p2, u2)*u2, q)/l_avg;
        else return glm::length(glm::dot(glm::cross(u1,u2),  q2 - p2))/glm::length(glm::cross(u1, u2))/l_avg;
    };

    for(int i = 1; i < (int)FC.size() - 2; i++)f += CalcPlanarity(FC[i].first->p3, FC[i].second->p3, FC[i+1].first->p3, FC[i+1].second->p3);

    if(IsLeft4pt)f += CalcPlanarity(FC[0].first->p3, FC[0].second->p3, FC[1].first->p3, FC[1].second->p3);
    if(IsRight4pt)f += CalcPlanarity(FC.back().first->p3, FC.back().second->p3, FC.end()[-2].first->p3, FC.end()[-2].second->p3);

    return f;
}

double RevisionVertices::E_iso(std::vector<std::vector<glm::f64vec3>>& P, std::vector<std::vector<glm::f64vec3>>& Pori, std::vector<double>& grad){
    double e = 0.0;
    auto f = [](std::vector<glm::f64vec3>& P, std::vector<glm::f64vec3>& Pori)->double{
        double val;
         for(int i = 1; i < (int)P.size(); i++)val += std::abs(glm::length(P[i] - P[i-1]) - glm::length(Pori[i] - Pori[i-1]));
        return val;
    };
    for(int i = 0; i < (int)P.size(); i++)e += f(P[i], Pori[i]);
    int i = 0;
    double fp, fm;
    for(int n = 0; n < (int)P.size(); n++){
        for(int j = 1; j < (int)P[n].size() - 1; j++){
            P[n][j].x += eps; fp = f(P[n], Pori[n]); P[n][j].x -= 2.0 * eps; fm = f(P[n], Pori[n]);  P[n][j].x += eps; grad[i++] = (fp - fm)/(2.0*eps);
            P[n][j].y += eps; fp = f(P[n], Pori[n]); P[n][j].y -= 2.0 * eps; fm = f(P[n], Pori[n]);  P[n][j].y += eps; grad[i++] = (fp - fm)/(2.0*eps);
            P[n][j].z += eps; fp = f(P[n], Pori[n]); P[n][j].z -= 2.0 * eps; fm = f(P[n], Pori[n]);  P[n][j].z += eps; grad[i++] = (fp - fm)/(2.0*eps);
        }
    }
    return e;
}

double RevisionVertices::E_fair(std::vector<std::vector<glm::f64vec3>>& P, std::vector<double>& grad){

    auto f = [](std::vector<glm::f64vec3>& P)->double{
        double val;
        for(int i = 1; i < (int)P.size() - 1; i++)val += std::pow(glm::length(P[i-1] - 2.0 * P[i] + P[i+1]), 2);
        return val;
    };
    double val = 0.0; for(auto& p: P)val += f(p);
    int i = 0;
    double fp, fm;
    for(auto&p : P){
        for(int j = 1; j < (int)p.size() - 1; j++){
            p[j].x += eps; fp = f(p); p[j].x -= 2.0 * eps; fm = f(p);  p[j].x += eps; grad[i++] = (fp - fm)/(2.0*eps);
            p[j].y += eps; fp = f(p); p[j].y -= 2.0 * eps; fm = f(p);  p[j].y += eps; grad[i++] = (fp - fm)/(2.0*eps);
            p[j].z += eps; fp = f(p); p[j].z -= 2.0 * eps; fm = f(p);  p[j].z += eps; grad[i++] = (fp - fm)/(2.0*eps);
        }
    }
    return val;
}

double RevisionVertices::E_sim(std::vector<std::vector<glm::f64vec3>>& P, std::vector<std::vector<glm::f64vec3>>& Pori, std::vector<double>& grad){
    double e = 0.0;
    auto f = [](std::vector<glm::f64vec3>& P, std::vector<glm::f64vec3>& Pori)->double{
        double val;
         for(int i = 0; i < (int)P.size(); i++)val += glm::length(P[i] - Pori[i]);
        return val;
    };
    for(int i = 0; i < (int)P.size(); i++)e += f(P[i], Pori[i]);
    int i = 0;
    double fp, fm;
    for(int n = 0; n < (int)P.size(); n++){
        for(int j = 1; j < (int)P[n].size() - 1; j++){
            P[n][j].x += eps; fp = f(P[n], Pori[n]); P[n][j].x -= 2.0 * eps; fm = f(P[n], Pori[n]);  P[n][j].x += eps; grad[i++] = (fp - fm)/(2.0*eps);
            P[n][j].y += eps; fp = f(P[n], Pori[n]); P[n][j].y -= 2.0 * eps; fm = f(P[n], Pori[n]);  P[n][j].y += eps; grad[i++] = (fp - fm)/(2.0*eps);
            P[n][j].z += eps; fp = f(P[n], Pori[n]); P[n][j].z -= 2.0 * eps; fm = f(P[n], Pori[n]);  P[n][j].z += eps; grad[i++] = (fp - fm)/(2.0*eps);
        }
    }
    return e;
}

double RevisionVertices::Minimize_SmoothSrf(const std::vector<double>& X, std::vector<double> &grad, void *f_data){
    SmoothSurface *od = (SmoothSurface*)f_data;
    std::vector<SmoothingArea> SA = od->SA;
    std::vector<std::vector<glm::f64vec3>> P,Pori;
    int i = 0;
    for(auto&sa: SA){
        std::vector<glm::f64vec3> tmp = {sa.stP.first->p3}, tmp_ori = {sa.stP.first->p3_ori};

        for(auto&p: sa.OriginalVertices){
            tmp_ori.push_back(p->p3_ori);
            tmp.push_back(glm::f64vec3{X[i], X[i+1], X[i+2]}); i += 3;
        }tmp.push_back(sa.lastP.first->p3); tmp_ori.push_back(sa.lastP.first->p3_ori);
        P.push_back(tmp); Pori.push_back(tmp_ori);
    }
    double f_sim, f_fair, f_iso;
    double w_sim = 100, w_fair = 100, w_iso = 0.1;
    if(!grad.empty()){
        std::vector<double> dE_sim(grad.size()), dE_fair(grad.size()), dE_iso(grad.size());
        f_sim = RevisionVertices::E_sim(P, Pori, dE_sim); f_fair = RevisionVertices::E_fair(P, dE_fair), f_iso = RevisionVertices::E_iso(P, Pori, dE_iso);
        for(int i = 0; i < (int)grad.size(); i++)grad[i] = w_sim * dE_sim[i] + w_fair * dE_fair[i] + w_iso * dE_iso[i];

    }
    std::cout<< "similarilty = " << f_sim << "  ,  fairness = " << f_fair << ", isometric = " << f_iso << std::endl;
    return w_sim * f_sim + w_fair * f_fair + w_iso * f_iso;

}

bool FoldLine::SimpleSmooothSrf(const std::vector<Vertex*>& Poly_v){

    auto getV = [](Vertex *o, Vertex *x, Vertex *o2, Vertex *x2)->Vertex*{
        if(IsParallel(o, x, o2, x2))return nullptr;
        glm::f64vec3 p2d = calcCrossPoint_2Vector(o, x, o2, x2);
        glm::f64vec3 p3d = calcTargetDistanceOnPlane(p2d, o,  x, x2);
        return  new Vertex(p2d, p3d);
    };

    glm::f64vec3 r3d, r2d, _r3d, _r2d;
    std::vector<int>Vertices_Ind;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}

    Vertex *qt;
    for(int j = 0; j < (int)Vertices_Ind.size() - 1; j++){
        qt = nullptr;
        if(j != 0 && j != (int)Vertices_Ind.size() - 2){
            qt = getV(FoldingCurve[Vertices_Ind[j]].first, FoldingCurve[Vertices_Ind[j]].second, FoldingCurve[Vertices_Ind[j+1]].first, FoldingCurve[Vertices_Ind[j+1]].second);
             _r3d = glm::normalize(FoldingCurve[Vertices_Ind[j]].second->p3 - FoldingCurve[Vertices_Ind[j]].first->p3), _r2d = glm::normalize(FoldingCurve[Vertices_Ind[j]].second->p - FoldingCurve[Vertices_Ind[j]].first->p);
        }
        for(int i = Vertices_Ind[j] + 1; i < Vertices_Ind[j+1]; i++){
            FoldingCurve[i].first->p3 = FoldingCurve[i].first->p3_ori; FoldingCurve[i].first->p = FoldingCurve[i].first->p2_ori;
            if(j == 0){
                r2d = glm::normalize(FoldingCurve[Vertices_Ind[1]].second->p - FoldingCurve[Vertices_Ind[1]].first->p);
                r3d = glm::normalize(FoldingCurve[Vertices_Ind[1]].second->p3 - FoldingCurve[Vertices_Ind[1]].first->p3);
            }else if(j == (int)Vertices_Ind.size() - 2){
                r2d = glm::normalize(FoldingCurve[Vertices_Ind[j]].second->p - FoldingCurve[Vertices_Ind[j]].first->p);
                r3d = glm::normalize(FoldingCurve[Vertices_Ind[j]].second->p3 - FoldingCurve[Vertices_Ind[j]].first->p3);
            }else{
                if(qt != nullptr){
                    r2d = glm::normalize(qt->p - FoldingCurve[i].first->p); r3d = glm::normalize(qt->p3 - FoldingCurve[i].first->p3);
                    double k = glm::dot(glm::normalize(FoldingCurve[i-1].first->p - FoldingCurve[i].first->p), glm::normalize(FoldingCurve[i+1].first->p - FoldingCurve[i].first->p));
                    k = (k < -1)? std::numbers::pi: (k > 1)? 0: std::acos(k);
                    glm::f64vec3 Nori = glm::cross(FoldingCurve[i-1].first->p - FoldingCurve[i].first->p, FoldingCurve[i+1].first->p - FoldingCurve[i].first->p);
                    if(Nori.z > 0)k = 2.0*std::numbers::pi - k;
                    double phi = glm::dot(glm::normalize(FoldingCurve[i-1].first->p - FoldingCurve[i].first->p), r2d);
                    phi = (phi < -1)? std::numbers::pi: (phi > 1)? 0: std::acos(phi);
                    glm::f64vec3 N = glm::cross(FoldingCurve[i-1].first->p - FoldingCurve[i].first->p, r2d);
                    if(N.z > 0)phi = 2.0*std::numbers::pi - phi;
                    if(phi > k)r2d *= -1;
                    if(glm::dot(r3d, _r3d) < 0)r3d *= -1;
                }else{
                    r2d = glm::normalize(FoldingCurve[Vertices_Ind[j]].second->p - FoldingCurve[Vertices_Ind[j]].first->p);
                    r3d = glm::normalize(FoldingCurve[Vertices_Ind[j]].second->p3 - FoldingCurve[Vertices_Ind[j]].first->p3);
                }
            }

            glm::f64vec3 p;
            for(int k = 0; k < (int)Poly_v.size(); k++){
                p = calcCrossPoint_2Vector(FoldingCurve[i].first->p, 1000.0 * r2d + FoldingCurve[i].first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
                if(MathTool::is_point_on_line(p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p) && glm::dot(p - FoldingCurve[i].first->p, r2d) > 0)break;
            }
            FoldingCurve[i].second->p = p;
            FoldingCurve[i].second->p3 = glm::length(FoldingCurve[i].second->p - FoldingCurve[i].first->p) * r3d + FoldingCurve[i].first->p3;
            FoldingCurve[i].third->p = FoldingCurve[i].third->p2_ori;
            FoldingCurve[i].third->p3 = FoldingCurve[i].third->p3_ori;
        }
    }
    Vertex* v_clst = getClosestVertex(FoldingCurve[0].second, FoldingCurve[0].first, FoldingCurve, false);
    if(v_clst != nullptr){
        int j;
        for(j = 0; j < (int)FoldingCurve.size(); j++){if(v_clst == FoldingCurve[j].second)break;}
        FoldingCurve[0].second->p3 = calcTargetDistanceOnPlane(FoldingCurve[0].second->p, FoldingCurve[j].first, FoldingCurve[j].second, FoldingCurve[j+1].second);
    }
    v_clst = getClosestVertex(FoldingCurve.back().second, FoldingCurve.back().first, FoldingCurve, false);
        if(v_clst != nullptr){
            int j;
            for(j = 0; j < (int)FoldingCurve.size(); j++){if(v_clst == FoldingCurve[j].second)break;}
            FoldingCurve.back().second->p3 = calcTargetDistanceOnPlane(FoldingCurve.back().second->p, FoldingCurve[j].first, FoldingCurve[j].second, FoldingCurve[j-1].second);
        }
    delete qt;
    return true;
}

std::vector<std::vector<glm::f64vec3>> FoldLine::Optimization_SmooothSrf(const std::vector<Vertex*>& Poly_v, bool IsConnectEndPoint){
    std::vector<double> X;
    std::vector<std::vector<glm::f64vec3>> res_qt;
    Vertex4d bef = FoldingCurve.front();
    std::vector<RevisionVertices::SmoothingArea> SA;
    std::vector<Vertex*> OriginalVertices;
    int bef_ind;
    for(auto&FC : FoldingCurve){
        if(!FC.IsCalc){
            X.push_back(FC.first->p3_ori.x); X.push_back(FC.first->p3_ori.y); X.push_back(FC.first->p3_ori.z);
            OriginalVertices.push_back(FC.first);
        }
        else{
            int ind = std::distance(FoldingCurve.begin(), std::find(FoldingCurve.begin(), FoldingCurve.end(), FC));
            if(FC == bef){bef_ind = ind; continue;}
            RevisionVertices::SmoothingArea sm = {bef, FC, bef_ind, ind, OriginalVertices};
            if(bef_ind == 0 || ind == (int)FoldingCurve.size() - 1){sm.qt = nullptr; sm.qb = nullptr;}
            SA.push_back(sm);
            bef_ind = ind;
            bef = FC;
            OriginalVertices.clear();
        }
    }
    RevisionVertices::ObjData_smooth od = {SA, IsConnectEndPoint};
    if(X.empty())return res_qt;
    nlopt::opt opt;
    opt = nlopt::opt(nlopt::LD_MMA, X.size());
    opt.set_min_objective(RevisionVertices::Minimize_SmoothSrf, &od);
    opt.add_inequality_constraint(RevisionVertices::Const_Edev, &od);
    opt.set_xtol_rel(1e-13);

    double minf;
    try {
        nlopt::result result = opt.optimize(X, minf);
        std::cout <<"result :  smoothing  " <<result << std::endl;
        std::cout << "found minimum at f = "  << std::setprecision(10) << minf << std::endl;
        int i = 0;
        for(auto&FC: FoldingCurve){
            if(!FC.IsCalc){
                FC.first->p3 = glm::f64vec3{X[i], X[i+1], X[i+2]}; i += 3;
            }
        }
        glm::f64vec3 p, r3d, r2d;
        for(auto sm: od.SA){
            glm::f64vec3 vt = glm::normalize(sm.stP.second->p3 - sm.stP.first->p3), vt2d = glm::normalize(sm.stP.second->p - sm.stP.first->p);
            glm::f64vec3 SpinAxis = (vt2d.y < 0)? glm::f64vec3{0,0,-1}: glm::f64vec3{0,0,1};
            for(int j = sm.st_ind + 1; j < sm.last_ind; j++){
                glm::f64vec3 e = FoldingCurve[j].first->p3 - FoldingCurve[j-1].first->p3;
                double l = glm::length(e), phi = std::acos(glm::dot(vt, glm::normalize(e)));
                FoldingCurve[j].first->p = l * glm::f64vec3{glm::rotate(phi, SpinAxis)* glm::f64vec4{ vt2d, 1.0}} + FoldingCurve[j-1].first->p;

                if(sm.qt == nullptr){
                    if(sm.st_ind != 0){r2d = glm::normalize(sm.stP.second->p - sm.stP.first->p); r3d = glm::normalize(sm.stP.second->p3 - sm.stP.first->p3);}
                    else {r2d = glm::normalize(sm.lastP.second->p - sm.lastP.first->p); r3d = glm::normalize(sm.lastP.second->p3 - sm.lastP.first->p3);}
                }else{
                    res_qt.push_back({sm.stP.first->p3, sm.qt->p3, sm.lastP.first->p3});

                    r2d = glm::normalize(sm.qt->p - FoldingCurve[j].first->p);
                    if(glm::dot(r2d, FoldingCurve[sm.st_ind].second->p - FoldingCurve[sm.st_ind].first->p) < 0)r2d *= -1;
                    FoldingCurve[j].second->p = calcCrossPoint_2Vector(FoldingCurve[j].first->p, 1000.0 * r2d + FoldingCurve[j].first->p, sm.stP.second->p, sm.lastP.second->p);
                    r3d = glm::normalize(sm.qt->p3 - FoldingCurve[j].first->p3);
                    if(glm::dot(r3d, glm::normalize(sm.stP.second->p3 - sm.stP.first->p3)) < 0)r3d *= -1;
                    vt =glm::normalize(sm.qt->p3 - FoldingCurve[j].first->p3); vt2d =glm::normalize(sm.qt->p - FoldingCurve[j].first->p);
                    SpinAxis = (vt2d.y < 0)? glm::f64vec3{0,0,-1}: glm::f64vec3{0,0,1};
                }
                for(int k = 0; k < (int)Poly_v.size(); k++){
                    p = calcCrossPoint_2Vector(FoldingCurve[j].first->p, 1000.0 * r2d + FoldingCurve[j].first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
                    if(MathTool::is_point_on_line(p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p) && glm::dot(r2d, p - FoldingCurve[j].first->p) > 0)break;
                }
                FoldingCurve[j].second->p = p;
                FoldingCurve[j].second->p3 = glm::length(FoldingCurve[j].second->p - FoldingCurve[j].first->p) * r3d + FoldingCurve[j].first->p3;

                if(sm.qb == nullptr){
                    if(sm.st_ind != 0){r2d = glm::normalize(sm.stP.third->p - sm.stP.first->p); r3d = glm::normalize(sm.stP.third->p3 - sm.stP.first->p3);}
                    else {r2d = glm::normalize(sm.lastP.third->p - sm.lastP.first->p); r3d = glm::normalize(sm.lastP.third->p3 - sm.lastP.first->p3);}
                }else{
                    r2d = glm::normalize(sm.qb->p - FoldingCurve[j].first->p);r3d = glm::normalize(sm.qb->p3 - FoldingCurve[j].first->p3);
                    if(glm::dot(r2d, FoldingCurve[sm.st_ind].third->p - FoldingCurve[sm.st_ind].first->p) < 0)r2d *= -1;
                }
                for(int k = 0; k < (int)Poly_v.size(); k++){
                    p = calcCrossPoint_2Vector(FoldingCurve[j].first->p, 1000.0 * r2d + FoldingCurve[j].first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
                    if(MathTool::is_point_on_line(p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p) && glm::dot(r2d, p - FoldingCurve[j].first->p) > 0)break;
                }
                FoldingCurve[j].third->p = p;
                FoldingCurve[j].third->p3 = glm::length(p - FoldingCurve[j].first->p) * r3d + FoldingCurve[j].first->p3;
            }
        }
        Vertex* v_clst = getClosestVertex(FoldingCurve[0].second, FoldingCurve[0].first, FoldingCurve, false);
        if(v_clst != nullptr){
            int j; for(j = 0; j < (int)FoldingCurve.size(); j++){if(v_clst == FoldingCurve[j].second)break;}
            FoldingCurve[0].second->p3 = calcTargetDistanceOnPlane(FoldingCurve[0].second->p, FoldingCurve[j].first, FoldingCurve[j].second, FoldingCurve[j+1].second);
        }
        v_clst = getClosestVertex(FoldingCurve.back().second, FoldingCurve.back().first, FoldingCurve, false);
        if(v_clst != nullptr){
            int j; for(j = 0; j < (int)FoldingCurve.size(); j++){if(v_clst == FoldingCurve[j].second)break;}
            FoldingCurve.back().second->p3 = calcTargetDistanceOnPlane(FoldingCurve.back().second->p, FoldingCurve[j].first, FoldingCurve[j].second, FoldingCurve[j-1].second);
        }
        std::cout<<"smoothing finish"<<std::endl;
    }catch (std::exception& e) {std::cout << "nlopt failed: " << e.what() << std::endl; }

    for(auto&sm: od.SA){
        delete sm.qb; delete sm.qt;
    }

    return res_qt;
}

bool FoldLine::Optimization_FlapAngle(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, double wb, double wp){
    if(FoldingCurve.empty())return false;

    std::vector<int> Vertices_Ind;for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}

    double a = 0, a2 = 0.0;
    glm::f64vec3 e = glm::normalize(FoldingCurve[Vertices_Ind[0]].first->p3 - FoldingCurve[Vertices_Ind[1]].first->p3),
            e2 = glm::normalize(FoldingCurve[Vertices_Ind[2]].first->p3 - FoldingCurve[Vertices_Ind[1]].first->p3);
    glm::f64vec3 Nb = -glm::normalize(glm::cross(e, e2)), N4;
    N4 = -glm::normalize(glm::cross(FoldingCurve[Vertices_Ind[1]].third->p3 - FoldingCurve[Vertices_Ind[1]].first->p3, FoldingCurve[Vertices_Ind[0]].first->p3 - FoldingCurve[Vertices_Ind[1]].first->p3));
    double a_ll = std::atan2(glm::dot(glm::cross(Nb,N4), e),glm::dot(-Nb, N4));
    if(a_ll < 0)a_ll += 2.0*std::numbers::pi;
    double a_con = a_ll + std::numbers::pi;
    if(a_con > 2.0*std::numbers::pi)a_con -= 2.0*std::numbers::pi;
    if(a_con < 0)a_con +=2.0*std::numbers::pi;

    double _a_con = a_con + std::numbers::pi;
    if(_a_con < 0)_a_con += 2.0*std::numbers::pi;else if(_a_con > 2.0* std::numbers::pi)_a_con -= 2.0*std::numbers::pi;

    std::cout <<glm::degrees(a_ll) << " , " << glm::degrees(a_con)<<std::endl;
    std::ofstream ofs2;
    std::filesystem::create_directory("./Optimization");

    std::string AngleFile = "./Optimization/ChangeAngle.csv" ;
    ofs2.open(AngleFile, std::ios::out); ofs2 << "a(radian),a(degree),Eruling, Ebend"<<std::endl;
    while(a <= 2.0*std::numbers::pi){
         _FoldingAAAMethod(FoldingCurve, Poly_V, a);
         double f = RulingsCrossed(FoldingCurve);
         double fb = Fbend2(FoldingCurve);
        double val = (f == 0)? fb: 1;
         ofs2 << a << "," << a * 180.0/std::numbers::pi << ", " << f << ", " << fb << ", " << val <<  std::endl;
        a += 1e-3;
    }ofs2.close();

    glm::f64vec3 Axis = (FoldingCurve[Vertices_Ind[1]].third->p3 - FoldingCurve[Vertices_Ind[1]].first->p3)/glm::length((FoldingCurve[Vertices_Ind[1]].third->p3 - FoldingCurve[Vertices_Ind[1]].first->p3));
    double phi3 = std::acos(glm::dot(glm::normalize(FoldingCurve[Vertices_Ind[2]].first->p3 - FoldingCurve[Vertices_Ind[1]].first->p3),Axis));
    double phi4 = std::acos(glm::dot(glm::normalize(FoldingCurve[Vertices_Ind[0]].first->p3 - FoldingCurve[Vertices_Ind[1]].first->p3),Axis));
    double k = 2.0 * std::numbers::pi - phi3 - phi4;
    double a_min, a_max;
    bool IsMount;
    for(const auto&e: Edges){
        if((e->vertex == FoldingCurve[Vertices_Ind[1]].first && e->next->vertex == FoldingCurve[Vertices_Ind[1]].second)){
            if(e->r->Gradation < 0)IsMount = false;
            else IsMount = true;
            break;
        }
    }

    if(k < std::numbers::pi && IsMount){a_min = a_con + std::numbers::pi; a_max = 2.0 * std::numbers::pi;}
    if(k >= std::numbers::pi && IsMount){a_min = std::numbers::pi; a_max = std::numbers::pi + a_con;}
    if(k < std::numbers::pi && !IsMount){a_min = 0.0;a_max = a_con - std::numbers::pi;}
    if(k >= std::numbers::pi && !IsMount){a_min = a_con - std::numbers::pi; a_max = std::numbers::pi;}
    a = (a_min + a_max)/2.0;
    //another version
    //if(k < std::numbers::pi && IsMount){a_min = a_con - std::numbers::pi; a_max = a_con; a = a_min + 0.1;}
    //if(k >= std::numbers::pi && IsMount){a_min = a_con - 2.0*std::numbers::pi; a_max = a_con - std::numbers::pi; a = a_max - 0.1;}
    //if(k < std::numbers::pi && !IsMount){a_min = a_con;a_max = a_con + std::numbers::pi; a = a_max - 0.1;}
    //if(k >= std::numbers::pi && !IsMount){a_min = a_con - std::numbers::pi; a_max = a_con; a = a_min + 0.1;}
    //a_min = std::min(a_con, _a_con); a_max = std::max(a_con, _a_con);
    RevisionVertices::ObjData od = {FoldingCurve, Vertices, Poly_V}; od.AddWeight(wb, wp);
    nlopt::opt opt;
    opt = nlopt::opt(nlopt::LD_MMA, 1);
    opt.set_min_objective(ObjFunc_FlapAngle, &od);
    opt.set_lower_bounds(a_min);
    opt.set_upper_bounds(a_max);
    opt.add_inequality_constraint(Fruling, &od);
    //opt.set_param("inner_maxeval", 100);
    opt.set_maxtime(2.0);//stop over this time
    //opt.set_ftol_rel(1e-10);
    opt.set_xtol_rel(1e-13);
    //
    //a = a_min;
    std::cout << "area " << glm::degrees(a_min) << " < " << glm::degrees(a) << " < " << glm::degrees(a_max) << std::endl;

    AngleFile = "./Optimization/OptimArea.csv";
    ofs2.open(AngleFile, std::ios::out);
    ofs2 << "a(radian), a(degree) , Eruling, Ebend, Ebend"<<std::endl;
    for(double _a = a_min; _a <= a_max; _a+= 1e-3){
         _FoldingAAAMethod(FoldingCurve, Poly_V, _a);
         double f = RulingsCrossed(FoldingCurve);
         double fb = Fbend2(FoldingCurve);
         double val = (f == 0)? fb: 1;
        ofs2 << _a << ", " << glm::degrees(_a) << " , " << f << ", " << fb << ", "<< val << std::endl;
    }ofs2.close();

    double minf_amin, minf_amax, f_ruling_amin, f_ruling_amax;
    double res_amin, res_amax;
    ofs_Ebend.open(File_Ebend, std::ios::out); ofs_Eruling.open(File_Eruling, std::ios::out);
      try {
        std::vector<double> _a{a_min + 0.5};
          nlopt::result result = opt.optimize(_a, minf_amin);

          std::cout <<"result :  lower bound" <<result << std::endl;
          f_ruling_amin = RulingsCrossed(FoldingCurve);
          res_amin = _a[0];
          std::cout << "found minimum at f(" << glm::degrees(_a[0]) << ") = " << std::setprecision(10) << minf_amin << "  ,  " << f_ruling_amin <<  std::endl;          
      }
      catch (std::exception& e) {
          std::cout << "nlopt failed: " << e.what() << std::endl;
      }ofs_Ebend.close(); ofs_Eruling.close();

    try {
      std::vector<double> _a{a_max - 0.5};
        nlopt::result result = opt.optimize(_a, minf_amax);
        f_ruling_amax = RulingsCrossed(FoldingCurve);
        std::cout <<"result :  upper bound" <<result << std::endl;
        std::cout << "found minimum at f(" << glm::degrees(_a[0]) << ") = " << std::setprecision(10) << minf_amax << "  ,  "  << f_ruling_amax << std::endl;
        res_amax = _a[0];

    }
    catch (std::exception& e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
    double th_ruling = 1e-9;
    if(f_ruling_amin < th_ruling && f_ruling_amax < th_ruling){
        if(minf_amax > minf_amin)a2 = res_amin;
        else a2 = res_amax;
    }else if(f_ruling_amin < th_ruling && f_ruling_amax >= th_ruling){a2 = res_amin; }
    else if(f_ruling_amin >= th_ruling && f_ruling_amax < th_ruling){ a2 = res_amax; }else{
        std::cout<<"could not avoid ruling cross"<<std::endl;
        return false;
    }
   _FoldingAAAMethod(FoldingCurve, Poly_V, a2);
   std::cout << "result : smaller = " << glm::degrees(a2) << ",  f_bend = " << wb * Fbend2(FoldingCurve) << "  ,  f_paralell = " << wp * Fparallel(FoldingCurve) << std::endl;
    std::cout << "finish"<<std::endl;
    return true;
}

void FoldLine::Optimization_Vertices(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V){
    if(FoldingCurve.empty())return;
    glm::f64vec3 Nb = glm::normalize(glm::cross(FoldingCurve[2].first->p3 - FoldingCurve[1].first->p3, FoldingCurve[0].first->p3 - FoldingCurve[1].first->p3));
    glm::f64vec3 Nf = glm::normalize(glm::cross(FoldingCurve[1].first->p3 - FoldingCurve[0].first->p3, FoldingCurve[1].second->p3 - FoldingCurve[1].first->p3));
    glm::f64vec3 SpinAxis = glm::normalize(FoldingCurve[0].first->p3 - FoldingCurve[1].first->p3);
    double a = (glm::dot(Nb, Nf) < -1)? std::numbers::pi: (glm::dot(Nb, Nf) > 1)? 0: std::acos(glm::dot(Nb, Nf));
    if(glm::dot(glm::cross(Nb, Nf), SpinAxis) > 0)a = 2.0 * std::numbers::pi - a;

    glm::f64vec3 N4;
    for(auto&he: Edges){
        if(he->vertex->p3 == FoldingCurve[0].first->p3 && he->next->vertex->p3 == FoldingCurve[1].first->p3){N4 = he->face->getNormalVec();break;}
    }

    double a_ll = std::atan2(glm::dot(glm::cross(Nb,N4), glm::normalize(FoldingCurve[0].first->p3 - FoldingCurve[1].first->p3)),glm::dot(-Nb, N4));
    if(a_ll < 0)a_ll += 2.0*std::numbers::pi;
    double a_con = a_ll + std::numbers::pi;
    if(a_con > 2.0*std::numbers::pi)a_con -= 2.0*std::numbers::pi;
    if(a_con < 0)a_con +=2.0*std::numbers::pi;

    glm::f64vec3 Axis = (FoldingCurve[1].third->p3 - FoldingCurve[1].first->p3)/glm::length((FoldingCurve[1].third->p3 - FoldingCurve[1].first->p3));
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

    if(k < std::numbers::pi && IsMount){a_min = a_con + std::numbers::pi; a_max = 2.0 * std::numbers::pi;}
    if(k >= std::numbers::pi && IsMount){a_min = std::numbers::pi; a_max = std::numbers::pi + a_con;}
    if(k < std::numbers::pi && !IsMount){a_min = 0.0;a_max = a_con - std::numbers::pi;}
    if(k >= std::numbers::pi && !IsMount){a_min = a_con - std::numbers::pi; a_max = std::numbers::pi;}

    RevisionVertices::ObjData_v od_v = RevisionVertices::ObjData_v{a, FoldingCurve, Vertices, Poly_V};
    std::vector<double> X(FoldingCurve.size() + 1, 1.0);
    std::vector<double>UpperBnds(FoldingCurve.size() + 1, 0), LowerBnds(FoldingCurve.size() + 1,0);
    for(int i = 0; i < (int)FoldingCurve.size(); i++)
        UpperBnds[i] = glm::distance(FoldingCurve[i].second->p2_ori, FoldingCurve[i].third->p2_ori)/
                glm::distance(FoldingCurve[i].first->p2_ori, FoldingCurve[i].third->p2_ori);
    UpperBnds.back() = a_max; LowerBnds.back() = a_min; X.back() = a;
    nlopt::opt opt_v(nlopt::LD_MMA, X.size());

    opt_v.set_lower_bounds(LowerBnds);
    opt_v.set_upper_bounds(UpperBnds);

    opt_v.set_min_objective(RevisionVertices::Minimize_Vpos, &od_v);
    opt_v.add_inequality_constraint(RevisionVertices::F_ruling, &od_v);
    opt_v.set_param("inner_maxeval", 300);
    opt_v.set_xtol_rel(1e-8);
    double minf;
    try {
          nlopt::result result = opt_v.optimize(X, minf);
          std::cout <<"result :  " <<result << std::endl;
          std::cout << "found minimum at f(" << glm::degrees(X.back()) << ") = " << std::setprecision(10) << minf << std::endl;
          for(int i = 0; i < (int)FoldingCurve.size(); i++){
              FoldingCurve[i].first->p = X[i] * (FoldingCurve[i].first->p2_ori - FoldingCurve[i].third->p2_ori) + FoldingCurve[i].third->p2_ori;
              FoldingCurve[i].first->p3 = X[i] * (FoldingCurve[i].first->p3_ori - FoldingCurve[i].third->p3_ori) + FoldingCurve[i].third->p3_ori;
          }
         _FoldingAAAMethod(FoldingCurve, Poly_V, X.back());
    }catch (std::exception& e) {std::cout << "nlopt failed: " << e.what() << std::endl; }
}

inline double RevisionVertices::getK(const glm::f64vec3 o, const glm::f64vec3 x, const glm::f64vec3 x2){
    double k = std::acos(glm::dot(glm::normalize(x - o), glm::normalize(x2 - o)));
    if(glm::dot(glm::f64vec3{0,0,1}, glm::cross(glm::normalize(x - o), glm::normalize(x2 - o))) > 0)k = 2.0*std::numbers::pi - k;
    return k;
}
inline glm::f64vec3 RevisionVertices::getVertex(const glm::f64vec3& o, const glm::f64vec3& p, const double t){return t * (p - o) + o;}

double RevisionVertices::F_vpos(FoldLine3d &FC, const std::vector<double>& T, std::vector<double>& grad){
    double f = 0.0;
    for(int i = 0; i < (int)FC.size(); i++){
        f += glm::distance(FC[i].first->p2_ori, getVertex(FC[i].third->p2_ori, FC[i].first->p2_ori, T[i]));;
        double fp = glm::distance(FC[i].first->p2_ori, (T[i] + eps) * (FC[i].first->p2_ori - FC[i].third->p2_ori) + FC[i].third->p2_ori);
        double fm = glm::distance(FC[i].first->p2_ori, (T[i] - eps) * (FC[i].first->p2_ori - FC[i].third->p2_ori) + FC[i].third->p2_ori);
        grad[i] = (fp - fm)/(2.0*eps);
    }
    if(FC.size() != grad.size())grad.back() = 0;
    return f;
}

double RevisionVertices::F_k(FoldLine3d &FC, const std::vector<double>& T, std::vector<double>& grad){
    double f = 0.0, f_ori, fp, fm;

    if(!grad.empty()){
        glm::f64vec3 v, vbef, vnext;
        vbef = getVertex(FC[0].third->p2_ori, FC[0].first->p2_ori, T[0]); v = getVertex(FC[1].third->p2_ori, FC[1].first->p2_ori, T[1]);
        vnext = getVertex(FC[2].third->p2_ori, FC[2].first->p2_ori, T[2]);
        f_ori = getK(FC[1].first->p2_ori, FC[0].first->p2_ori, FC[2].first->p2_ori);
        fp = abs(f_ori - getK(v, getVertex(FC[0].third->p2_ori, FC[0].first->p2_ori, T[0] + eps), vnext));
        fm = abs(f_ori - getK(v, getVertex(FC[0].third->p2_ori, FC[0].first->p2_ori, T[0] - eps), vnext));
        grad[0] = (fp - fm)/(2.0 * eps);
        f += abs(f_ori - getK(v, vbef, vnext));

        for(int i = 1; i < (int)FC.size() - 1; i++){
            v = getVertex(FC[i].third->p2_ori, FC[i].first->p2_ori, T[i]); vbef = getVertex(FC[i-1].third->p2_ori, FC[i-1].first->p2_ori, T[i-1]);
            vnext = getVertex(FC[i+1].third->p2_ori, FC[i+1].first->p2_ori, T[i+1]);
            f_ori = getK(FC[i].first->p2_ori, FC[i-1].first->p2_ori, FC[i+1].first->p2_ori);
            fp = abs(f_ori - getK(getVertex(FC[i].third->p2_ori, FC[i].first->p2_ori, T[i] + eps), vbef, vnext));
            fm = abs(f_ori - getK(getVertex(FC[i].third->p2_ori, FC[i].first->p2_ori, T[i] - eps), vbef, vnext));
            f += abs(getK(v, vbef, vnext) - f_ori);
            grad[i] = (fp - fm)/(2.0*eps);
        }
        int lastInd = (grad.size() != FC.size())? -1: 0;
        v = getVertex(FC.end()[-2].third->p2_ori, FC.end()[-2].first->p2_ori, T.end()[-2 +lastInd]); vbef = getVertex(FC.end()[-3].third->p2_ori, FC.end()[-3].first->p2_ori, T.end()[-3+lastInd]);
        vnext = getVertex(FC.back().third->p2_ori, FC.back().first->p2_ori, T.end()[-1 + lastInd]);
        f_ori = getK(FC.end()[-2].first->p2_ori, FC.end()[-3].first->p2_ori, FC.back().first->p2_ori);
        fp = abs(f_ori - getK(v, vbef, getVertex(FC.back().third->p2_ori, FC.back().first->p2_ori, T.end()[-1 + lastInd] + eps)));
        fm = abs(f_ori - getK(v, vbef, getVertex(FC.back().third->p2_ori, FC.back().first->p2_ori, T.end()[-1 * lastInd] - eps)));
        grad.end()[-1 + lastInd] = (fp - fm)/(2.0 * eps);
        f += abs(f_ori - getK(v, vbef, vnext));
        if(lastInd == -1)grad.back() = 0;
    }
    return  f;
}

double RevisionVertices::F_ruling(const std::vector<double> &X, std::vector<double> &grad, void* f_data){
    auto fr =[](FoldLine3d& FoldingCurve)->double{
        double f = 0.0;
        Eigen::Matrix2d A;
        Eigen::Vector2d b;
        for(int i = 1; i < (int)FoldingCurve.size() - 1; i++){
            int cnt = 0;
            double f2 = 0;
            for(int i2 = 1; i2 < FoldingCurve.size() - 1; i2 += 1){
                if(i == i2)continue;
                glm::f64vec3 p1 = FoldingCurve[i].first->p, q1 = FoldingCurve[i].second->p, p2 = FoldingCurve[i2].first->p, q2 = FoldingCurve[i2].second->p;
                double t = ((p2.x - p1.x)*(p2.y - q2.y) - (p2.x - q2.x)*(p2.y - p1.y))/((q1.x - p1.x) * (p2.y - q2.y) - (p2.x - q2.x)*(q1.y - p1.y));
                if(0 < t && t < 1) f += 1.0/t;
                continue;
                A(0,0) = FoldingCurve[i].second->p.x - FoldingCurve[i].first->p.x; A(0,1) = -(FoldingCurve[i2].second->p.x - FoldingCurve[i2].first->p.x);
                A(1,0) = FoldingCurve[i].second->p.y - FoldingCurve[i].first->p.y; A(1,1) = -(FoldingCurve[i2].second->p.y - FoldingCurve[i2].first->p.y);
                b(0) = FoldingCurve[i2].first->p.x - FoldingCurve[i].first->p.x;
                b(1) = FoldingCurve[i2].first->p.y - FoldingCurve[i].first->p.y;
                Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
                if(0 < ts(0) && ts(0) < 1){ f2 += 1.0/ts(0); cnt++; }
            }
            f += cnt*f2;
        }
        return f;
    };

    double f = 0.0;
    RevisionVertices::ObjData_v *data = static_cast<RevisionVertices::ObjData_v*>(f_data);
    FoldLine3d FC = data->FC;
    double a = X.back();

    if(!grad.empty()){
        std::vector<std::array<glm::f64vec3, 2>> tmp;
        for(int i = 0; i < (int)FC.size(); i++ ){
            tmp.push_back({FC[i].first->p, FC[i].first->p3});
            FC[i].first->p = getVertex(FC[i].third->p2_ori, FC[i].first->p2_ori, X[i]);
            FC[i].first->p3 = getVertex(FC[i].third->p3_ori, FC[i].first->p3_ori, X[i]);
        } _FoldingAAAMethod(FC, data->Poly_V, a);

        f = fr(FC);
        for(int i = 0; i < (int)FC.size(); i++){
            FC[i].first->p = getVertex(FC[i].third->p2_ori, FC[i].first->p2_ori, X[i] + eps); FC[i].first->p3 = getVertex(FC[i].third->p3_ori, FC[i].first->p3_ori, X[i] + eps);
            _FoldingAAAMethod(FC, data->Poly_V, a);
            double fp = fr(FC);
            FC[i].first->p = getVertex(FC[i].third->p2_ori, FC[i].first->p2_ori, X[i] - eps); FC[i].first->p3 = getVertex(FC[i].third->p3_ori, FC[i].first->p3_ori, X[i] - eps);
            _FoldingAAAMethod(FC, data->Poly_V, a);
            double fm = fr(FC);
            grad[i] = (fp - fm)/(2.0*eps);
            FC[i].first->p = tmp[i][0]; FC[i].first->p3 = tmp[i][1];
        }
        if(grad.size() != FC.size()){
            _FoldingAAAMethod(FC, data->Poly_V, a + eps);
            double fp = fr(FC);
            _FoldingAAAMethod(FC, data->Poly_V, a - eps);
            double fm = fr(FC);
            grad.back() = (fp - fm)/(2.0*eps);
        }
    }
    std::cout <<"F_ruling(vertices) = " << f << std::endl;
    return f;
}

double RevisionVertices::Minimize_Vpos(const std::vector<double> &X, std::vector<double> &grad, void* f_data){
    ObjData_v *od = (ObjData_v *)f_data;
    FoldLine3d FC = od->FC;
    double fk, fvpos, fa;

    if(!grad.empty()){
        double a = (grad.size() != FC.size())? X.back(): od->a;
        std::vector<double> grad_vpos(grad.size()), grad_k(grad.size()), grad_a(grad.size(), 0);
        std::vector<Vertex*> Poly_V = od->Poly_V;
        for(int i = 0; i < (int)FC.size(); i++ ){
            FC[i].first->p = getVertex(FC[i].third->p2_ori, FC[i].first->p2_ori, X[i]);
            FC[i].first->p3 = getVertex(FC[i].third->p3_ori, FC[i].first->p3_ori, X[i]);
        }  _FoldingAAAMethod(FC, Poly_V, a);
        fk = F_k(FC, X, grad_k); fvpos = F_vpos(FC, X, grad_vpos);
        fa = abs(a - od->a); grad.back() = (abs(a + eps - od->a) - abs(a - eps - od->a)/(2.0*eps));
        for(int i = 0; i < (int)X.size(); i++){
            grad[i] = grad_k[i] + grad_vpos[i] + grad_a[i];
        }
    }
    std::cout <<"fk = " << fk <<", fvpos = " << fvpos <<", fa = " << fa <<  "  fk + fvpos = " << fk + fvpos + fa << std::endl;
    return fk + fvpos;
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
            //if(std::find_if(T_crs.begin(), T_crs.end(), [&e](CrvPt_FL* T){return T->p == e->vertex->p || T->p == e->next->vertex->p;}) != T_crs.end())continue;
            std::vector<double>arcT = BezierClipping(CtrlPts, e, dim);

            for(auto&t: arcT){
                if(t < 0 || 1 < t){std::cout<<"t is not correct value " << t << std::endl; continue;}
                if(!T_crs.empty() && std::find_if(T_crs.begin(), T_crs.end(), [&t](CrvPt_FL* T){return abs(T->s - t) < 1e-9;}) != T_crs.end()){std::cout<<"has" <<std::endl; continue;}
                glm::f64vec3 v2{0,0,0};
                for (int i = 0; i < int(CtrlPts.size()); i++) v2 += MathTool::BernsteinBasisFunc(dim, i, t) * CtrlPts[i];
                //if(!is_point_on_line(v2, e->vertex->p, e->next->vertex->p)){std::cout<<"not on line  "<< t << std::endl; continue;}
                double sa = glm::distance(v2, e->vertex->p), sc = glm::distance(e->vertex->p, e->next->vertex->p);
                glm::f64vec3 v3 = sa/sc * (e->next->vertex->p3 - e->vertex->p3) + e->vertex->p3;
                t_max = std::max(t_max, t); t_min = std::min(t_min, t);
                CrvPt_FL *P = new CrvPt_FL(v2, v3, t);
                P->set(v2, e->vertex, e->next->vertex);
                T_crs.push_back(P);
            }
        }
    }
    FoldingCurve.clear();
    std::sort(T_crs.begin(), T_crs.end(), [](CrvPt_FL* T, CrvPt_FL* T2){return T->s > T2->s;});
    for(auto&t: T_crs){
        Vertices.push_back(t);
    }
    glm::f64vec3 UpVec{0,-1,0};
    for(auto&t: T_crs){
        if(t == T_crs.front() || t == T_crs.back()){
            for(int i = 0; i < (int)Poly_v.size(); i++){
                if(!MathTool::is_point_on_line(t->p, Poly_v[i]->p, Poly_v[(i + 1) % (int)Poly_v.size()]->p))continue;
                Vertex *top, *bottom;
                if(glm::dot(UpVec, glm::normalize(Poly_v[i]->p - Poly_v[(i + 1) % (int)Poly_v.size()]->p)) > 0){
                    for(auto&he: Edges){
                        if(he->vertex->p == Poly_v[i]->p)top = he->vertex;
                        if(he->vertex->p == Poly_v[(i + 1) % (int)Poly_v.size()]->p)bottom = he->vertex;
                    }
                }
                else {
                    for(auto&he: Edges){
                        if(he->vertex->p == Poly_v[i]->p)bottom = he->vertex;
                        if(he->vertex->p == Poly_v[(i + 1) % (int)Poly_v.size()]->p)top = he->vertex;
                    }
                }
                FoldingCurve.push_back(Vertex4d(t, top, bottom));
                break;
            }
        }else{
            for(auto&he: Edges){
                if(!MathTool::is_point_on_line(t->p, he->vertex->p, he->next->vertex->p))continue;
                if(glm::dot(UpVec, glm::normalize(he->next->vertex->p - he->vertex->p)) > 0)FoldingCurve.push_back(Vertex4d(t, he->next->vertex, he->vertex));
                else FoldingCurve.push_back(Vertex4d(t, he->vertex, he->next->vertex));
                break;
            }
        }

    }
    return true;


    for(auto&t: T_crs){
        std::vector<HalfEdge*> H_new;
        for(auto&he: Edges){
            H_new = he->Split(t, Edges);
            if(H_new.size() > 0)break;
        }
    }

    std::vector<Face*> Faces_new;

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
            if(h1->v->s > h2->v->s){FoldingCurve.push_back(Vertex4d(InsertPoints[1], h2->prev->vertex, h1->next->next->vertex));}
            else {FoldingCurve.push_back(Vertex4d(InsertPoints[0], h1->prev->vertex, h2->next->next->vertex));}
        }
    }
    glm::f64vec3 befN = glm::normalize(FoldingCurve.front().second->p - FoldingCurve.front().first->p);


    for(int i = 0; i < (int)FoldingCurve.size(); i++){
        if(FoldingCurve[i].first == T_crs.front() || FoldingCurve[i].first == T_crs.back() )FoldingCurve.erase(FoldingCurve.begin() + i);
    }

    std::sort(FoldingCurve.begin(), FoldingCurve.end(), [](Vertex4d& F, Vertex4d& F2){return F.first->s > F2.first->s;});

    FoldingCurve.insert(FoldingCurve.begin(), Vertex4d(T_crs.front(), FindAnotherPoint(Poly_v, Edges, T_crs.front(), befN), FindAnotherPoint(Poly_v, Edges, T_crs.front(), -befN)));
    FoldingCurve.push_back(Vertex4d(T_crs.back(), FindAnotherPoint(Poly_v, Edges, T_crs.back(), befN), FindAnotherPoint(Poly_v, Edges, T_crs.back(), -befN)));
    Faces.insert(Faces.end(), Faces_new.begin(), Faces_new.end()); // 連結

    return true;
}

void FoldLine::SimplifyModel(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, double tol){
    for(auto&fc: FoldingCurve){
        fc.IsCalc = true;
        fc.first->p = fc.first->p2_ori; fc.first->p3 = fc.first->p3_ori;
        fc.second->p = fc.second->p2_ori; fc.second->p3 = fc.second->p3_ori;
        fc.third->p = fc.third->p2_ori; fc.third->p3 = fc.third->p3_ori;
    }
    auto tmp = TrimPoints2(Edges, Faces, Vertices, FoldingCurve, tol);
    for(auto&V4d: FoldingCurve){
        if(std::find(tmp.begin(), tmp.end(), V4d) != tmp.end()){
        }else V4d.IsCalc = false;
    }
    std::vector<int>Vertices_Ind;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}
    for(int i = 0; i < (int)Vertices_Ind.size() - 1; i++){
    int ind_branch = Vertices_Ind[i + 1], ind_root = Vertices_Ind[i];
    if(i == 0){
        if(getClosestVertex(FoldingCurve[Vertices_Ind[0]].second , FoldingCurve[Vertices_Ind[0]].first, FoldingCurve) != nullptr){
            ind_root = Vertices_Ind[1]; ind_branch = Vertices_Ind[2];
        }
    }
    if(i == (int)Vertices_Ind.size() -2 ){
        if(getClosestVertex(FoldingCurve[Vertices_Ind.back()].second , FoldingCurve[Vertices_Ind.back()].first, FoldingCurve) != nullptr)
        {
            ind_root = Vertices_Ind[i]; ind_branch = Vertices_Ind[i-1];
        }

    }
        for(int j = Vertices_Ind[i] + 1; j < Vertices_Ind[i+1]; j++){
            FoldingCurve[j].first->p =  calcCrossPoint_2Vector(FoldingCurve[j].second, FoldingCurve[j].third, FoldingCurve[Vertices_Ind[i]].first, FoldingCurve[Vertices_Ind[i+1]].first);
            FoldingCurve[j].first->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].first->p, FoldingCurve[ind_root].first, FoldingCurve[ind_branch].third, FoldingCurve[ind_branch].second);
            FoldingCurve[j].third->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].third->p, FoldingCurve[ind_root].first, FoldingCurve[ind_branch].first, FoldingCurve[ind_branch].second);
            FoldingCurve[j].second->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].second->p, FoldingCurve[ind_root].first, FoldingCurve[ind_branch].first, FoldingCurve[ind_branch].second);

        }
    }
}

bool FoldLine::RevisionCrosPtsPosition(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, int type, bool TrimMode){
    if(FoldingCurve.empty())return false;

    //ApproximatePolyLine();
    //先にpiに近い値を削除してもう一度挑戦
    if(FoldingCurve.size() > 5){
        std::cout << "before revision" << std::endl;
        for(int i = 1; i < (int)FoldingCurve.size()-1; i++){
            double k = std::acos(glm::dot(glm::normalize(FoldingCurve[i-1].first->p - FoldingCurve[i].first->p), glm::normalize(FoldingCurve[i+1].first->p - FoldingCurve[i].first->p)));
            if(glm::dot(glm::f64vec3{0,0,1}, glm::cross(glm::normalize(FoldingCurve[i-1].first->p - FoldingCurve[i].first->p), glm::normalize(FoldingCurve[i+1].first->p - FoldingCurve[i].first->p))) > 0)k = 2.0*std::numbers::pi - k;
            //std::cout << k << std::endl;
        }
        auto ReviseEndPosition = [](glm::f64vec3& o, glm::f64vec3 e, glm::f64vec3 &x_ori,  double k)->glm::f64vec3{
            return glm::normalize(glm::rotate (k, glm::f64vec3{0,0,-1})  * glm::f64vec4{e,1.0});
        };
        auto getCrossPosition = [](Vertex* v, Vertex* o, Vertex *p, Vertex *q){
            Eigen::Matrix2d A;
            Eigen::Vector2d b;
            glm::f64vec3 v2 = q->p - p->p;
            A(0,0) = v->p.x; A(0,1) = -v2.x; A(1,0) = v->p.y; A(1,1) = -v2.y;
            b(0) = p->p.x - o->p.x; b(1) = p->p.y - o->p.y;
            Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
            return  o->p + ts(0) * v->p;
        };

        double k0 = RevisionVertices::getK(FoldingCurve[1].first->p, FoldingCurve[0].first->p, FoldingCurve[2].first->p);
        double k1 = RevisionVertices::getK(FoldingCurve[2].first->p, FoldingCurve[1].first->p, FoldingCurve[3].first->p);
        double k2 = RevisionVertices::getK(FoldingCurve[3].first->p, FoldingCurve[2].first->p, FoldingCurve[4].first->p);
         std::cout <<"front " << abs(k0 - k1) << ", " << abs(k1 - k2) << " , " << 2*k1 - k2 << std::endl;
        //if(abs(k0 - k1) > abs(k1 - k2)){
            std::cout <<"front"<<std::endl;
            FoldingCurve[0].first->p = ReviseEndPosition(FoldingCurve[1].first->p, FoldingCurve[2].first->p - FoldingCurve[1].first->p, FoldingCurve[0].first->p, -(2.*k1 - k2));
            glm::f64vec3 newPos = getCrossPosition(FoldingCurve[0].first, FoldingCurve[1].first, FoldingCurve[0].first->vo, FoldingCurve[0].first->ve);
            FoldingCurve[0].first->p = newPos;
            FoldingCurve[0].first->p3 = FoldingCurve[0].first->vo->p3 + glm::length(newPos - FoldingCurve[0].first->vo->p)/glm::length(FoldingCurve[0].first->vo->p - FoldingCurve[0].first->ve->p)
                    * (FoldingCurve[0].first->ve->p3 - FoldingCurve[0].first->vo->p3);
       // }
        k0 = RevisionVertices::getK(FoldingCurve.end()[-2].first->p, FoldingCurve.end()[-3].first->p, FoldingCurve.back().first->p);
        k1 = RevisionVertices::getK(FoldingCurve.end()[-3].first->p, FoldingCurve.end()[-4].first->p, FoldingCurve.end()[-2].first->p);
        k2 = RevisionVertices::getK(FoldingCurve.end()[-4].first->p, FoldingCurve.end()[-5].first->p, FoldingCurve.end()[-3].first->p);
        std::cout <<"back " << abs(k0 - k1) << ", " << abs(k1 - k2) << ", " << 2.*k1 - k2 <<  std::endl;
        //if(abs(k0 - k1) > abs(k1 - k2)){
            std::cout << "back"<<std::endl;
            FoldingCurve.back().first->p = ReviseEndPosition(FoldingCurve.end()[-2].first->p, FoldingCurve.end()[-3].first->p - FoldingCurve.end()[-2].first->p, FoldingCurve.back().first->p, 2.*k1 - k2);
            glm::f64vec3 newPos2 = getCrossPosition(FoldingCurve.back().first, FoldingCurve.end()[-2].first, FoldingCurve.back().first->vo, FoldingCurve.back().first->ve);
            FoldingCurve.back().first->p = newPos2;
            FoldingCurve.back().first->p3 = FoldingCurve.back().first->vo->p3 + glm::length(newPos2 - FoldingCurve.back().first->vo->p)/
                    glm::length(FoldingCurve.back().first->vo->p - FoldingCurve.back().first->ve->p) * (FoldingCurve.back().first->ve->p3 - FoldingCurve.back().first->vo->p3);
        //}

    }

    double tol = 0 * std::numbers::pi/180.0;
    auto tmp = TrimPoints2(Edges, Faces, Vertices, FoldingCurve, tol);

    std::vector<int>Vertices_Ind;
    for(auto&V4d: FoldingCurve){
        if(std::find(tmp.begin(), tmp.end(), V4d) == tmp.end())V4d.IsCalc = false;
    }
    for(int i = 0; i < FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}
    for(int i = 0; i < (int)Vertices_Ind.size() - 1; i++){
        for(int j = Vertices_Ind[i] + 1; j < Vertices_Ind[i+1]; j++){
            FoldingCurve[j].first->p =  calcCrossPoint_2Vector(FoldingCurve[j].first, FoldingCurve[j].third, FoldingCurve[Vertices_Ind[i]].first, FoldingCurve[Vertices_Ind[i+1]].first);
            FoldingCurve[j].first->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].first->p, FoldingCurve[Vertices_Ind[i]].first, FoldingCurve[Vertices_Ind[i]].third, FoldingCurve[Vertices_Ind[i+1]].third);
            FoldingCurve[j].third->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].third->p, FoldingCurve[Vertices_Ind[i]].first, FoldingCurve[Vertices_Ind[i]].third, FoldingCurve[Vertices_Ind[i+1]].third);
        }
    }
    return true;

}

void FoldLine::ReassignColor(std::vector<HalfEdge*>& Edges, ColorPoint& CP){
    //Y字のとき折り線は山であると仮定(谷の場合は色を反転)
    typedef struct {
        ruling *line;
        int type_mvk;//0:mountain && k < pi, 1: mountain && k >= pi, 2: vally && k < pi, 3: vally && k >= pi, -1: gradation = 0
        //0の数 + 3の数 == rulingの数 || 1の数 + 2の数 == rulingの数 -> 山谷とkの割り当てが正しい
        //それ以外
         //0の数 > 2の数 -> 　2の色を変換 (逆もまた然り)(1と3でも同様に)
    }MVK;
    std::vector<MVK> InitState;
    std::vector<int> Vertices_Ind;for(int i = 0; i < (int)FoldingCurve.size();i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}
    std::vector<int>MV(5,0);
    for(int i = 1; i < (int)Vertices_Ind.size()-1; i++){
        double k = std::acos(glm::dot(glm::normalize(FoldingCurve[Vertices_Ind[i-1]].first->p - FoldingCurve[Vertices_Ind[i]].first->p),
                glm::normalize(FoldingCurve[Vertices_Ind[i+1]].first->p - FoldingCurve[Vertices_Ind[i]].first->p)));
        if(glm::dot(glm::f64vec3{0,0,1}, glm::cross(glm::normalize(FoldingCurve[Vertices_Ind[i-1]].first->p - FoldingCurve[Vertices_Ind[i]].first->p),
                                                    glm::normalize(FoldingCurve[Vertices_Ind[i+1]].first->p - FoldingCurve[Vertices_Ind[i]].first->p))) > 0)k = 2.0*std::numbers::pi - k;
       glm::f64vec3 N = glm::normalize(glm::cross(FoldingCurve[Vertices_Ind[i]].third->p3 - FoldingCurve[Vertices_Ind[i]].first->p3,
               FoldingCurve[Vertices_Ind[i-1]].first->p3 - FoldingCurve[Vertices_Ind[i]].first->p3));
       glm::f64vec3 Np = glm::normalize(glm::cross(FoldingCurve[Vertices_Ind[i+1]].first->p3 - FoldingCurve[Vertices_Ind[i]].first->p3,
               FoldingCurve[Vertices_Ind[i]].third->p3 - FoldingCurve[Vertices_Ind[i]].first->p3));
       double phi = std::acos(glm::dot(N, Np));

        double color;
        if(phi < CP.angle)color = CP.color/CP.angle *phi;
        else color = ((255.0 - CP.color)*(phi - CP.angle)/(std::numbers::pi - CP.angle) + CP.color);
        for(auto&e: Edges){
            if(e->edgetype == EdgeType::r && MathTool::is_point_on_line(FoldingCurve[Vertices_Ind[i]].first->p, e->next->vertex->p, e->vertex->p)){
                int mv = (e->r->Gradation == 0)? -1: (e->r->Gradation > 0)? 0: 1;
                if(mv == 0)e->r->Gradation = color; else e->r->Gradation = -color;
                int type_mvk = (mv == 0 && k < std::numbers::pi)? 0: (mv == 0 && k >= std::numbers::pi)? 1: (mv == 1 && k < std::numbers::pi)? 2: (mv == 1 && k >= std::numbers::pi)? 3: -1;
                InitState.push_back(MVK(e->r, type_mvk));
                MV[type_mvk + 1] += 1;
                break;
            }
        }
    }

    if(MV[0] != 0){std::cout << "no folding pattern was found"<< std::endl; return;}
    while(!(MV[1] + MV[4] == (int)Vertices_Ind.size() - 2 || MV[2] + MV[3] == (int)Vertices_Ind.size() - 2)){
        if(MV[1] == MV[2] && MV[3] == MV[4] && MV[1] == MV[4]){
        }else{
            for(auto& IS: InitState){
                if(MV[1] + MV[4] > MV[2] + MV[3] ){
                    if(IS.type_mvk == 1){IS.line->Gradation *= -1; IS.type_mvk = 3; MV[2]--; MV[4]++;}
                    if(IS.type_mvk == 2){IS.line->Gradation *= -1; IS.type_mvk = 0; MV[3]--; MV[1]++;}
                }else{
                    if(IS.type_mvk == 0){IS.line->Gradation *= -1; IS.type_mvk = 2; MV[1]--; MV[3]++;}
                    if(IS.type_mvk == 3){IS.line->Gradation *= -1; IS.type_mvk = 1; MV[4]--; MV[2]++;}
                }
            }
        }
    }
    std::cout <<"color changed"<<std::endl;
}

void FoldLine::modifyFoldingCurvePositionOn3d(const std::vector<HalfEdge*>& Edges){
    for(auto&fc: FoldingCurve){
        for(auto&h: Edges){
            if(!MathTool::is_point_on_line(fc.first->p, h->vertex->p, h->next->vertex->p))continue;
            double t = glm::length(fc.first->p - h->vertex->p)/glm::length(h->vertex->p - h->next->vertex->p);
            fc.first->p3_ori = fc.first->p3 = t * (h->next->vertex->p3 - h->vertex->p3) + h->vertex->p3;
            break;
        }
    }
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
            glm::f64vec3 Axis = (FoldingCurve[i].third->p3 - FoldingCurve[i].first->p3)/glm::length((FoldingCurve[i].third->p3 - FoldingCurve[i].first->p3));
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

void FoldLine::applyAAAMethod(std::vector<Vertex*>& Poly_V, double a){
    if(FoldingCurve.empty())return;
     _FoldingAAAMethod(FoldingCurve, Poly_V, a);
}

inline double update_flapangle(double a, const glm::f64vec3& befN, const glm::f64vec3& SrfN, const glm::f64vec3& e){
    double dir = glm::dot(-e, glm::cross(befN, SrfN));
    double tau = std::acos(glm::dot(SrfN, befN));
    if(dir <= 0)tau = 2.0*std::numbers::pi - tau;
    if(a - tau <= 0)return a - tau + 2.0 * std::numbers::pi;
    return a - tau;
}
inline glm::f64vec3 _calcruling3d(const double& a, glm::f64vec3 e, glm::f64vec3 e2, glm::f64vec3 axis, double& beta, std::vector<double>& Phi){
    e = glm::normalize(e); e2 = glm::normalize(e2); axis = glm::normalize(axis);
    beta = std::acos(glm::dot(e, e2));
    Phi.resize(4);
    Phi[2] = std::acos(glm::dot(e2,axis));
    Phi[3] = std::acos(glm::dot(e, axis));

    double k = 2.0 * std::numbers::pi - Phi[2] - Phi[3];
    double _phi1 = std::atan((std::cos(k) - std::cos(beta))/(std::sin(beta)*std::cos(a)- std::sin(k)));
    Phi[0] = (_phi1 < 0)? _phi1+std::numbers::pi: _phi1;
    Phi[1] = k - Phi[0];
    glm::f64vec3 N = glm::normalize(glm::rotate (a, -e)  * glm::f64vec4{glm::cross(e2, e),1.0});
    glm::f64vec3 r3d = glm::rotate(Phi[0], -N) * glm::f64vec4{e, 1.0};
    return glm::normalize(r3d);//新しいruling方向
}

void CalcRuling(double a, Vertex4d& xbef, Vertex4d& x, Vertex4d& xnext, std::vector<Vertex*>& Poly_V,  double& a2, glm::f64vec3& SrfN){

    glm::f64vec3 e = glm::normalize(xbef.first->p3 - x.first->p3), e2 = glm::normalize(xnext.first->p3 - x.first->p3);
    double beta;
    std::vector<double> Phi;
    glm::f64vec3 Axis = glm::normalize(x.third->p3 - x.first->p3);
    glm::f64vec3 r3d = _calcruling3d(a, e, e2, Axis, beta, Phi);
    glm::f64vec3 r2d = glm::rotate(Phi[0], glm::f64vec3{0,0,-1.0})* glm::f64vec4{(glm::normalize(xbef.first->p- x.first->p)), 1.0};r2d = glm::normalize(r2d);//展開図のruling方向

    glm::f64vec3 crossPoint;
    bool hasPointOnEdge = IsRulingCrossed(r2d, x.first->p, crossPoint, Poly_V);
    if(hasPointOnEdge){
        x.second->p = crossPoint;
        x.second->p3 = glm::distance(crossPoint, x.first->p) * r3d + x.first->p3;
    }
    double sin_a = sin(Phi[0])*sin(a)/sin(Phi[1]);
    if(sin_a > 1)sin_a = 1;
    else if(sin_a < -1)sin_a = -1;
    double cos_a = (cos(Phi[0]) - cos(Phi[1])*cos(beta))/(sin(Phi[1])*sin(beta));
    if(cos_a > 1)cos_a = 1;
    else if(cos_a < -1)cos_a = -1;
    a2 = std::atan2(sin_a, cos_a);
    SrfN = MathTool::ProjectionVector(glm::cross(e, e2), e2, true);
}

Vertex* getClosestVertex(Vertex *v, Vertex* o,  const std::vector<Vertex4d>& FoldingCurve, bool SkipTrimedPoint){
   Vertex *V_max = nullptr;
   double t_max = -1;
   for(auto&fc: FoldingCurve){
       if(SkipTrimedPoint && !fc.IsCalc)continue;
       if((glm::length(fc.second->p - v->p) < 1e-7|| glm::length(fc.second->p - o->p) < 1e-7))continue;
       if(MathTool::is_point_on_line(fc.second->p, v->p, o->p)){
           double t = glm::distance(fc.second->p, o->p)/glm::distance(v->p, o->p);
           if(t > t_max){
               t_max = t; V_max = fc.second;
           }
       }
   }
   return V_max;
}

void _FoldingAAAMethod(std::vector<Vertex4d>& FoldingCurve, std::vector<Vertex*>& Poly_V, double a){

    auto SetOnPlane = [&](Vertex4d& V, std::vector<Vertex*>& Poly_V, Vertex4d& p, Vertex4d& p2, Vertex *q, Vertex *q2, int IsEnd){
        //p, q: 左、p2, q2：右
        glm::f64vec3 vec;
        if(IsEnd == -1){vec = glm::normalize(q2->p - p2.first->p);}//左が端
        else if(IsEnd == 1){vec = glm::normalize(q->p - p.first->p);}//右が端
        else{
            if(!IsParallel(p.first, q, p2.first, q2))vec = glm::normalize(calcCrossPoint_2Vector(p.first, q, p2.first, q2) - V.first->p);
            else vec = glm::normalize(q->p - p.first->p);

            double k = glm::dot(glm::normalize(p2.first->p - V.first->p), glm::normalize(p.first->p - V.first->p)); k = (k < -1)? std::numbers::pi: (k > 1)? 0: std::acos(k);
            glm::f64vec3 N = glm::cross(p2.first->p - V.first->p, p.first->p - V.first->p);
            if(N.z > 0)k = 2.0*std::numbers::pi - k;
            double phi = glm::dot(glm::normalize(p2.first->p - V.first->p), vec); phi = (phi < -1)? std::numbers::pi: (phi > 1)? 0: std::acos(phi);
            N = glm::cross(p2.first->p - V.first->p, vec);if(N.z > 0)phi = 2.0*std::numbers::pi - phi;
            if(phi > k)vec *= -1;

        }
        vec = 1000. * vec + V.first->p;
        Vertex *Qv = new Vertex(vec);

        for(int k = 0; k < (int)Poly_V.size(); k++){
            glm::f64vec3 q;
            if(IsParallel(V.first, Qv, Poly_V[k], Poly_V[(k + 1) % (int)Poly_V.size()]))continue;
            q = calcCrossPoint_2Vector(V.first, Qv, Poly_V[k], Poly_V[(k + 1) % (int)Poly_V.size()]);
            if(MathTool::is_point_on_line(q, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) && MathTool::is_point_on_line(q, V.first->p, Qv->p) ){
                V.second->p = q;
                V.second->p3 = calcTargetDistanceOnPlane(V.second->p, q2, p.first, p2.first);
            }
        }
    };

    double a2 = 0, l;
    glm::f64vec3 e, e2, x;
    glm::f64vec3 N, r, befN;
    std::vector<int> Vertices_Ind;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}

    double phi02 = glm::angle(glm::normalize(FoldingCurve[Vertices_Ind[0]].second->p - FoldingCurve[Vertices_Ind[0]].first->p),
            glm::normalize(FoldingCurve[Vertices_Ind[1]].first->p - FoldingCurve[Vertices_Ind[0]].first->p));
    double phim1 = glm::angle(glm::normalize(FoldingCurve[Vertices_Ind.end()[-2]].first->p - FoldingCurve.back().first->p),
            glm::normalize(FoldingCurve.back().second->p - FoldingCurve.back().first->p));

    a = (a < 0.0)? a + 2.0 * std::numbers::pi: (a > 2.0*std::numbers::pi)? a - 2.0*std::numbers::pi: a;
    for(int ind = 1; ind < (int)Vertices_Ind.size() - 1; ind++){

        Vertex4d fc = FoldingCurve[Vertices_Ind[ind]];
        Vertex4d fc_bef = FoldingCurve[Vertices_Ind[ind - 1]];
        Vertex4d fc_next = FoldingCurve[Vertices_Ind[ind + 1]];

        x = fc.first->p3;
        e = (fc_bef.first->p3 - x)/glm::length((fc_bef.first->p3 - x));
        e2 = (fc_next.first->p3 - x)/glm::length((fc_next.first->p3 - x));
        if(ind != 1){glm::f64vec3 SrfN = MathTool::ProjectionVector(glm::cross(e, e2), -e, true);a = update_flapangle(a2, befN, SrfN, e);}
        CalcRuling(a, fc_bef, fc, fc_next, Poly_V, a2, befN);
        if(ind == (int)FoldingCurve.size() -2)break;
        befN = MathTool::ProjectionVector(glm::cross(e, e2), e2, true);
        if(ind != 1){
            for(int i = Vertices_Ind[ind-1] + 1; i < Vertices_Ind[ind]; i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind[ind]],FoldingCurve[Vertices_Ind[ind-1]], FoldingCurve[Vertices_Ind[ind]].second,  FoldingCurve[Vertices_Ind[ind-1]].second, 0);
        }
    }
    {//i = 0(端の面)
        Vertex *v_clst = getClosestVertex(FoldingCurve[Vertices_Ind[0]].second , FoldingCurve[Vertices_Ind[0]].first, FoldingCurve);
        if(v_clst == nullptr){
            x = FoldingCurve[Vertices_Ind[0]].first->p3;
            e = glm::normalize(FoldingCurve[Vertices_Ind[1]].first->p3 - x);
            glm::f64vec3 ee = glm::normalize(FoldingCurve[Vertices_Ind[1]].second->p3- FoldingCurve[Vertices_Ind[1]].first->p3);
            N = glm::normalize(glm::cross(ee,-e));
            r = glm::rotate(phi02, N) * glm::f64vec4{e,1};r = glm::normalize(r);//新しいruling方向
            l = glm::length(FoldingCurve[Vertices_Ind[0]].second->p - FoldingCurve[Vertices_Ind[0]].first->p);
            FoldingCurve[Vertices_Ind[0]].second->p3 = l * r + FoldingCurve[Vertices_Ind[0]].first->p3;
            for(int i = Vertices_Ind[0] + 1; i < Vertices_Ind[1]; i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind[1]], FoldingCurve[Vertices_Ind[0]], FoldingCurve[Vertices_Ind[1]].second, FoldingCurve[Vertices_Ind[0]].second, 1);
        }else{
            int j;
            for(j = 0; j < (int)Vertices_Ind.size(); j++){
                if(FoldingCurve[Vertices_Ind[j]].second == v_clst)break;
            }
            if(j == (int)Vertices_Ind.size())std::cout <<"can't find correct edge right" << std::endl;
            else{
                auto p = calcTargetDistanceOnPlane(FoldingCurve[Vertices_Ind[0]].second->p, v_clst,  FoldingCurve[Vertices_Ind[j]].first, FoldingCurve[Vertices_Ind[j+1]].first);
                FoldingCurve[Vertices_Ind[0]].second->p3 = p;
            }
            for(int i = Vertices_Ind[0] + 1; i < Vertices_Ind[1]; i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind[1]], FoldingCurve[Vertices_Ind[0]], FoldingCurve[Vertices_Ind[1]].second, v_clst, 1);
        }

    }

    {
        Vertex *v_clst = getClosestVertex(FoldingCurve[Vertices_Ind.back()].second , FoldingCurve[Vertices_Ind.back()].first, FoldingCurve);
        if(v_clst == nullptr){
            x = FoldingCurve[Vertices_Ind.back()].first->p3;
            e = glm::normalize(FoldingCurve[Vertices_Ind.end()[-2]].first->p3 - x);
            glm::f64vec3 ee = glm::normalize(FoldingCurve[Vertices_Ind.end()[-2]].second->p3 - FoldingCurve[Vertices_Ind.end()[-2]].first->p3);
            N = glm::normalize(glm::cross(ee,e));
            r = (glm::rotate(-phim1, N) * glm::f64vec4{e,1}); r = glm::normalize(r);//新しいruling方向
            l = glm::length(FoldingCurve[Vertices_Ind.back()].second->p - FoldingCurve[Vertices_Ind.back()].first->p);
            FoldingCurve[Vertices_Ind.back()].second->p3 = l * r + x;
            for(int i = Vertices_Ind.end()[-2] + 1; i < Vertices_Ind.back(); i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind.back()], FoldingCurve[Vertices_Ind.end()[-2]], FoldingCurve[Vertices_Ind.back()].second, FoldingCurve[Vertices_Ind.end()[-2]].second, -1);
        }else{
            int j;
            for(j = 0; j < (int)Vertices_Ind.size(); j++){
                if(FoldingCurve[Vertices_Ind[j]].second == v_clst)break;
            }
            if(j == (int)Vertices_Ind.size())std::cout <<"can't find correct edge right" << std::endl;
            else{
                auto p = calcTargetDistanceOnPlane(FoldingCurve[Vertices_Ind.back()].second->p, v_clst,  FoldingCurve[Vertices_Ind[j]].first, FoldingCurve[Vertices_Ind[j-1]].first);
                FoldingCurve[Vertices_Ind.back()].second->p3 = p;
            }
            for(int i = Vertices_Ind.end()[-2] + 1; i < Vertices_Ind.back(); i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind.back()], FoldingCurve[Vertices_Ind.end()[-2]], v_clst, FoldingCurve[Vertices_Ind.end()[-2]].second, -1);
        }

    }
}

//kokoga genninn
std::vector<Vertex4d> TrimPoints2(std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces, std::vector<Vertex*>& Vertices, std::vector<Vertex4d>& FoldingCurve, double tol){
    int Ind;
    double mink;
    std::vector<HalfEdge*>::iterator itr;
    std::vector<Vertex*>::iterator itr_v;
    std::vector<Vertex4d> res = FoldingCurve;
    res.clear();
    size_t st = FoldingCurve.size();

    while(1){
        Douglas_Peucker_algorithm(FoldingCurve,res, tol);
        if(res.size() == st)break;
        st = res.size();
    }

    return res;

    for(auto fc = FoldingCurve.begin(); fc != FoldingCurve.end();){
        if(std::find(res.begin(), res.end(), *fc) != res.end()){fc++; continue;}
            for(auto&e: Edges){
                if(e->vertex != (*fc).first) continue;
                std::vector<HalfEdge*> H = {nullptr, nullptr, nullptr, nullptr};
                std::vector<int>H_PolygonSize(4,0);
                for(auto&h: e->vertex->halfedge){
                    if(h->prev->vertex == (*(fc + 1)).first)H[0] = h;
                    if(h->next->vertex == (*(fc - 1)).first)H[1] = h;
                    else if(h->next->vertex == (*(fc + 1)).first)H[3] = h;
                    else if(h->prev->vertex == (*(fc - 1)).first) H[2] = h;
                    //tmp1 = (*(fc - 1)).first; tmp2 =(*(fc + 1)).first;
                }
                for(int i = 0; i < 4; i++){
                    HalfEdge *h = H[i];
                    do{
                        h = h->next; H_PolygonSize[i]++;
                    }while(h != H[i]);
                }

                H[0]->prev->pair = H[2]->prev; H[2]->prev->pair = H[0]->prev;
                H[2]->prev->next = H[3]->next; H[3]->next->prev = H[2]->prev;
                H[0]->prev->next = H[1]->next; H[1]->next->prev = H[0]->prev;
                if(H_PolygonSize[0] == 3){//hidariue ga 3
                   H[1]->next->next->next = H[0]->prev; H[0]->prev->prev = H[1]->next->next;
                   H[3]->prev->prev->next = H[2]->next->next; H[2]->next->next->prev = H[3]->prev->prev;
                }else if(H_PolygonSize[1] == 3){//migiue ga 3
                     H[1]->next->next = H[0]->prev->prev; H[0]->prev->prev->prev = H[1]->next;
                     H[2]->next->next->prev = H[3]->next->next; H[3]->next->next->next = H[2]->next->next;
                }else if(H_PolygonSize[2] == 3){//migisita ga 3
                    H[0]->next->next->prev = H[1]->next->next; H[1]->next->next->next = H[0]->next->next;
                    H[3]->next->next->next = H[2]->next->next; H[2]->next->next->prev = H[3]->next->next;
                }else if(H_PolygonSize[3] == 3){//hidarisita ga 3
                    H[1]->next->next->next = H[0]->next->next; H[0]->next->next->prev = H[1]->next->next;
                    H[3]->next->next = H[2]->prev->prev; H[2]->prev->prev->next = H[3]->next;
                }else{
                    H[0]->next->next->prev = H[1]->next->next; H[1]->next->next->next = H[0]->next->next;
                    H[3]->next->next->next = H[2]->next->next; H[2]->next->next->prev = H[3]->next->next;
                }


                H[0]->prev->face->ReConnect(H[0]->prev);
                H[2]->prev->face->ReConnect(H[2]->prev);

                auto itr_f = std::find(Faces.begin(), Faces.end(), H[1]->face);
                if(itr_f != Faces.end()){delete *itr_f; Faces.erase(itr_f);}else std::cout<<"can't find face"<<std::endl;
                itr_f = std::find(Faces.begin(), Faces.end(), H[3]->face);
                if(itr_f != Faces.end()){delete *itr_f; Faces.erase(itr_f);}else std::cout<<"can't find face"<<std::endl;

                itr_v = std::find(Vertices.begin(), Vertices.end(), (*fc).second);if(itr_v != Vertices.end()){delete (*fc).second; Vertices.erase(itr_v);}else std::cout <<"not found vertex"<<std::endl;
                itr_v = std::find(Vertices.begin(), Vertices.end(), (*fc).third);if(itr_v != Vertices.end()){delete (*fc).third; Vertices.erase(itr_v);}else std::cout <<"not found vertex"<<std::endl;
                itr_v = std::find(Vertices.begin(), Vertices.end(), (*fc).first);if(itr_v != Vertices.end()){delete (*fc).first; Vertices.erase(itr_v);}else std::cout <<"not found vertex"<<std::endl;

                itr = H.begin();
                auto itr_h0p = std::find(Edges.begin(), Edges.end(), H[1]->prev); if(itr_h0p != Edges.end()){delete *itr_h0p; Edges.erase(itr_h0p); }else std::cout<<"not found itr_h0p"<<std::endl;
                auto itr_h2p = std::find(Edges.begin(), Edges.end(), H[3]->prev); if(itr_h2p != Edges.end()){delete *itr_h2p; Edges.erase(itr_h2p); }else std::cout<<"not found itr_h2p"<<std::endl;
                auto itr_h0n = std::find(Edges.begin(), Edges.end(), H[0]->next); if(itr_h0n != Edges.end()){delete *itr_h0n; Edges.erase(itr_h0n); }else std::cout<<"not found itr_h0n"<<std::endl;
                auto itr_h2n = std::find(Edges.begin(), Edges.end(), H[2]->next); if(itr_h2n != Edges.end()){delete *itr_h2n; Edges.erase(itr_h2n); }else std::cout<<"not found itr_h2n"<<std::endl;
                for(auto&h: H){
                    itr = std::find(Edges.begin(), Edges.end(), h);
                    delete *itr; Edges.erase(itr); //
                }
                fc = FoldingCurve.erase(fc);
                break;
            }
    }

    return res;

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
    }while(Ind != -1);
    return res;
}

//間引いた後の数を指定
//候補をVlistにいれる、確定をresにいれる
void Douglas_Peucker_algorithm2(std::vector<Vertex4d>& FoldingCurve, std::vector<Vertex4d>& res, int size){
    struct Point{
        Vertex4d start;
        Vertex4d p;
        Vertex4d end;
        int index;
        double dist;
        Point(){}
        Point(Vertex4d s, Vertex4d _p, Vertex4d e, size_t i, double d): start(s), p(_p), end(e), index(i), dist(d) {}
        //bool operator =(const Point& P){start = P.start; p = P.p; end = P.end; index = P.index; dist = P.dist;}
    };

    auto perpendicularDistance = [](glm::f64vec3& p, glm::f64vec3& l_start, glm::f64vec3& l_end)->double{
        glm::f64vec3 line = glm::normalize(l_end - l_start);
        glm::f64vec3 OP = glm::normalize(p - l_start);
        glm::f64vec3 H = l_start + glm::dot(line, OP) * (p - l_start);
        return glm::distance(H, p);
    };
    auto DPA = [&](std::vector<Vertex4d>& FoldingCurve, const std::vector<Vertex4d>& res){
        double dmax = 0.0;
        size_t index = -1, end = FoldingCurve.size()-1;
        if((int)FoldingCurve.size() < 2)return Point(FoldingCurve.front(), FoldingCurve.front(), FoldingCurve.back(), index, dmax);
        for(size_t i = 1; i < end; i++){
            if(std::find(res.begin(), res.end(), FoldingCurve[i]) != res.end())continue;
            double d = perpendicularDistance(FoldingCurve[i].first->p, FoldingCurve[0].first->p, FoldingCurve.back().first->p);
            if (d > dmax){
                index = i; dmax = d;
            }
        }
        if(index != -1)return Point(FoldingCurve.front(), FoldingCurve[index], FoldingCurve.back(), index, dmax);
        else return Point(FoldingCurve.front(), FoldingCurve.front(), FoldingCurve.back(), index, dmax);
    };

    int index = 0;
    res.push_back(FoldingCurve.front());res.push_back(FoldingCurve.back());
    do{
        std::vector<Vertex4d> firstLine{FoldingCurve.begin(), FoldingCurve.begin()+index+1};
        std::vector<Vertex4d> lastLine{FoldingCurve.begin() + index, FoldingCurve.end()};
        Point p_first = DPA(firstLine, res), p_second = DPA(lastLine,res);
        Point p;//最も遠い点
        if(p_first.index == -1 && p_second.index == -1)return;
        if(p_first.index == -1 && p_second.index != -1)p = p_second;
        else if(p_first.index != -1 && p_second.index == -1)p = p_first;
        else {
            if(p_first.dist > p_second.dist)p = p_first; else p = p_second;
        }
        for(int i = 1; i < res.size(); i++){
            if(res[i].first->s < p.p.first->s && p.p.first->s < res[i-1].first->s){res.insert(res.begin() + i, p.p);break;}
        }
        index = p.index;
    }while(res.size() < size);
}


//https://issekinichou.wordpress.com/2022/12/28/douglas-peucker-algorithm/
void Douglas_Peucker_algorithm(std::vector<Vertex4d>& FoldingCurve, std::vector<Vertex4d>& res, double tol){
    auto perpendicularDistance = [](glm::f64vec3& p, glm::f64vec3& l_start, glm::f64vec3& l_end)->double{
        glm::f64vec3 line = glm::normalize(l_end - l_start);
        glm::f64vec3 OP = glm::normalize(p - l_start);
        glm::f64vec3 H = l_start + glm::dot(line, OP) * (p - l_start);
        return glm::distance(H, p);
    };
    if((int)FoldingCurve.size() < 2)return;
    double dmax = 0.0;
    size_t index = 0;
    size_t end = FoldingCurve.size()-1;
    for(size_t i = 1; i < end; i++){
        double d = perpendicularDistance(FoldingCurve[i].first->p2_ori, FoldingCurve[0].first->p2_ori, FoldingCurve.back().first->p2_ori);
        if (d > dmax){index = i; dmax = d;}
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
