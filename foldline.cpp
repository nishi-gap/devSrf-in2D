
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
std::vector<Vertex4d> TrimPoints2(std::vector<Vertex4d>& FoldingCurve, double tol);
bool IsRulingCrossed(glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint,  std::vector<Vertex*>& Poly_V);
void Douglas_Peucker_algorithm(std::vector<Vertex4d>& FoldingCurve, std::vector<Vertex4d>& res, double tol = std::numbers::pi/9.0);
void Douglas_Peucker_algorithm2(std::vector<Vertex4d>& FoldingCurve, std::vector<Vertex4d>& res, int size);
Vertex* getClosestVertex(Vertex *v, Vertex* o,  const std::vector<Vertex4d>& FoldingCurve, bool SkipTrimedPoint = true);

inline glm::f64vec3 calcCrossPoint_2Vertex(Vertex* p1, Vertex* q1, Vertex* p2, Vertex* q2){
    Eigen::Matrix2d A; Eigen::Vector2d b;
    glm::f64vec3 v1 = q1->p - p1->p, v2 = q2->p - p2->p;
    b(0) = p2->p.x - p1->p.x; b(1) = p2->p.y - p1->p.y;
    A(0,0) = v1.x; A(0,1) = -v2.x;
    A(1,0) = v1.y; A(1,1) = -v2.y;
    return MathTool::calcCrossPoint_2Vector(p1->p, q1->p, p2->p, q2->p);

    double t = ((p2->p.x - p1->p.x)*(p2->p.y - q2->p.y) - (p2->p.x - q2->p.x)*(p2->p.y - p1->p.y))/((q1->p.x - p1->p.x) * (p2->p.y - q2->p.y) - (p2->p.x - q2->p.x)*(q1->p.y - p1->p.y));
    return glm::f64vec3{t * (q1->p.x - p1->p.x) + p1->p.x, t * (q1->p.y - p1->p.y) + p1->p.y, 0};
    Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
    return  x(0) * v1 + p1->p;
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
        double wb, wp ,wr;
        void AddWeight(double _wb, double _wp, double _wr){ wb = _wb; wp = _wp; wr = _wr;}

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
    public:
        Vertex4d stP, lastP;
        std::vector<Vertex*> OriginalVertices;
        int st_ind, last_ind;
        Vertex *qt, *qb;

        SmoothingArea(Vertex4d& a, Vertex4d& _end, int si, int li, std::vector<Vertex*>& OV): stP{a}, lastP{_end}, st_ind{si}, last_ind{li}, OriginalVertices{OV}{
            qt = getV(stP.first, stP.second, lastP.first, lastP.second);
            qb = getV(stP.first, stP.third, lastP.first, lastP.third);
        }

        SmoothingArea(Vertex4d& a, Vertex4d& _end, int si, int li): stP{a}, lastP{_end}, st_ind{si}, last_ind{li}{
            qt = getV(stP.first, stP.second, lastP.first, lastP.second);
            qb = getV(stP.first, stP.third, lastP.first, lastP.third);
        }

        ~SmoothingArea(){}
    private:
        Vertex *getV(Vertex *o, Vertex *x, Vertex *o2, Vertex *x2){
            if(IsParallel(o, x, o2, x2))return nullptr;
            glm::f64vec3 p2d = calcCrossPoint_2Vertex(o, x, o2, x2);
            glm::f64vec3 p3d = calcTargetDistanceOnPlane(p2d, o,  x, x2);
            return  new Vertex(p2d, p3d);
        }


    };
    struct SmoothSurface{
        std::vector<SmoothingArea> SA;
        bool IsConnectEndPoint;

        double Edev(std::vector<glm::f64vec3>& P, bool IsConnectEndPoint, const double th = 1e-5){
            auto Angle = [](glm::f64vec3& e, glm::f64vec3& e2)->double{return (glm::dot(e, e2) >= 1)? 0: (glm::dot(e, e2) <= -1)? std::numbers::pi: std::acos(glm::dot(e, e2));  };
            double f = 0.0;
            int j = 0;
            std::vector<std::array<glm::f64vec3, 3>> X;
            glm::f64vec3 rt, rb;
            for(const auto&sa: SA){
                rt = glm::normalize(sa.stP.second->p3 - sa.stP.first->p3);
                rb = glm::normalize(sa.stP.third->p3 - sa.stP.first->p3);
                X.push_back(std::array{sa.stP.first->p3, rt, rb});
                for(int i = 0; i < (int)sa.OriginalVertices.size(); i++){
                    rt = (sa.qt != nullptr)? glm::normalize(sa.qt->p3 - P[j]): glm::normalize(sa.stP.second->p3 - sa.stP.first->p3);
                    rb = (sa.qb != nullptr)? glm::normalize(sa.qb->p3 - P[j]): glm::normalize(sa.stP.third->p3 - sa.stP.first->p3);
                    X.push_back(std::array{P[j], rt, rb});
                    j++;
                }
            }
            rt = glm::normalize(SA.back().lastP.second->p3 - SA.back().lastP.first->p3);
            rb = glm::normalize(SA.back().lastP.third->p3 - SA.back().lastP.first->p3);
            X.push_back(std::array{SA.back().lastP.first->p3, rt, rb});

            glm::f64vec3 el, er;
            for(int i = 1; i < X.size() - 1; i++){
                el = glm::normalize(X[i+1][0] - X[i][0]); er = glm::normalize(X[i-1][0] - X[i][0]);
                f += std::abs((2.0 * std::numbers::pi - (Angle(el, X[i][1]) + Angle(X[i][1], er) + Angle(er, X[i][2]) + Angle(X[i][2], el))) - th);
            }

            return f;
        }

        double Econv(std::vector<glm::f64vec3>& P, std::vector<glm::f64vec3>& Pori){
            auto SgdArea = [](glm::f64vec3& A, glm::f64vec3& B, glm::f64vec3& C)->double{
              return ((A.x*B.y - B.x*A.y) + (B.x*C.y - C.x*B.y) + (C.x*A.y - A.x*C.y))/2.0;
            };
            double f = 0.0;
            double sum = 0.0;
            if(P.empty())return 0.0;
            int j = 0;
            std::vector<glm::f64vec3> X, Xori;
            for(const auto&sa: SA){
                X.push_back(sa.stP.first->p); Xori.push_back(sa.stP.first->p);
                for(int i = 0; i < (int)sa.OriginalVertices.size(); i++){
                    X.push_back(P[j]); Xori.push_back(Pori[j]);
                    j++;
                }
            }
            X.push_back(SA.back().lastP.first->p); Xori.push_back(SA.back().lastP.first->p);
            for(int i = 1; i < (int)X.size() - 1; i++){
                f =  -SgdArea(X[i-1],X[i],X[i+1]) * SgdArea(Xori[i-1],Xori[i],Xori[i+1]);
                sum += (f > 0)? f: 0;
            }
            return sum;
        }
    };

    using ObjData = OptimizeParam;
    using ObjData_v = OptimizeParam_v;
    using ObjData_smooth = SmoothSurface;
    using FoldCrv = std::vector<Vertex4d>;

    inline double getK(const glm::f64vec3 o, const glm::f64vec3 x, const glm::f64vec3 x2);

    double Const_Edev(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    double Const_Econv(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    double Const_Eplanarity(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    double E_fair(std::vector<std::vector<glm::f64vec3>>& P, std::vector<double>& grad);
    double E_sim(std::vector<std::vector<glm::f64vec3>>& P, std::vector<std::vector<glm::f64vec3>>& Pori, std::vector<double>& grad);
    double E_iso(std::vector<std::vector<glm::f64vec3>>& P, std::vector<std::vector<glm::f64vec3>>& Pori, std::vector<double>& grad);
    double E_normal(std::vector<glm::f64vec3>& R, const std::vector<Vertex4d>& FC);
    double Minimize_PlanaritySrf(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    double Minimize_SmoothSrf(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    glm::f64vec3 decideRulingDirectionOn3d(glm::f64vec3 e, glm::f64vec3 N, double a, double phi);

}

FoldLine::FoldLine(PaintTool _type)
{
    CtrlPts.clear();
    Rulings_2dL.clear();
    Rulings_2dR.clear();
    Rulings_3dL.clear();
    Rulings_3dR.clear();
    int crvNum = 300;
    CurvePts.resize(crvNum);
    color = 0;
    type = _type;
}

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
        glm::f64vec3 v = glm::normalize(FoldingCurve[i].second->p - FoldingCurve[i].first->p), v2 = glm::normalize(FoldingCurve[i+1].second->p - FoldingCurve[i+1].first->p);
        double _f = glm::dot(v, v2);
        f += 1.0 - _f;
    }
    return f;
}

double Fbend2(std::vector<Vertex4d>& FoldingCurve){
    std::vector<int> Vertices_Ind;
    double f = 0.0;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}
    glm::f64vec3 Nt = glm::normalize(glm::cross(FoldingCurve[0].first->p3 - FoldingCurve[Vertices_Ind[1]].first->p3, FoldingCurve[Vertices_Ind[1]].second->p3 - FoldingCurve[Vertices_Ind[1]].first->p3));
    for(int i = 1; i < (int)Vertices_Ind.size() - 1; i++){
        glm::f64vec3 Ntp = glm::normalize(glm::cross(FoldingCurve[Vertices_Ind[i]].first->p3 - FoldingCurve[Vertices_Ind[i+1]].first->p3,
                FoldingCurve[Vertices_Ind[i+1]].second->p3 - FoldingCurve[Vertices_Ind[i+1]].first->p3));
        double phi = ((glm::dot(Ntp,Nt)) > 1)? std::numbers::pi: ((glm::dot(Ntp,Nt)) < -1)? 0:  std::numbers::pi - std::acos(glm::dot(Ntp,Nt));
        f += 1.0/(phi*phi);
        Nt = Ntp;
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

double RulingsCrossed2(std::vector<Vertex4d>& FoldingCurve){
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
            f += 1.0/ts(0);
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
    double fb =  Fbend2(FoldingCurve), fp = Fparallel(FoldingCurve);
    double fr = (od->wr != -1)?RulingsCrossed2(FoldingCurve): 0;
    if(!grad.empty()){
        _FoldingAAAMethod(FoldingCurve, Poly_V, a[0] + eps);
       double fp = Fbend2(FoldingCurve);
       double fp2 = Fparallel(FoldingCurve);
       double fp3 = (od->wr != -1)?RulingsCrossed2(FoldingCurve): 0;
        _FoldingAAAMethod(FoldingCurve, Poly_V, a[0] - eps);
       double fm = Fbend2(FoldingCurve);
       double fm2 = Fparallel(FoldingCurve);
       double fm3 = (od->wr != -1)?RulingsCrossed2(FoldingCurve): 0;
       grad[0] = od->wb * (fp - fm)/(2.0 * eps) + od->wp * (fp2 - fm2)/(2.0 * eps) + 100.0*(fp3 - fm3)/(2.0*eps);
    }
    //std::cout <<"Fbend(" << glm::degrees(a[0]) << ")  =  " <<  fb <<  " ,  " << grad[0] << "  ,  Fparalell = " << fp << ", Fruling2 = "<< fr <<  std::endl;
    if(ofs_Ebend.is_open())ofs_Ebend << a[0] <<", " << glm::degrees(a[0]) << ", " <<  fb << ", " << fp <<   ",  " << grad[0] << std::endl;
    return od->wb * fb + od->wp * fp + 100.0 * fr;
}

glm::f64vec3 RevisionVertices::decideRulingDirectionOn3d(glm::f64vec3 e, glm::f64vec3 N, double a, double phi){
    e = glm::normalize(e);
    N = glm::normalize(glm::rotate(a, -e) * glm::f64vec4{glm::normalize(N), 1.0});
    glm::f64vec3 r = (glm::rotate(phi, N) * glm::f64vec4{e, 1.0}); r /= glm::length(r);
    return r;
}

double RevisionVertices::Const_Edev(const std::vector<double>& X, std::vector<double> &grad, void *f_data){
    SmoothSurface *od = (SmoothSurface*)f_data;
    std::vector<SmoothingArea> SA = od->SA;
    bool ICE = od->IsConnectEndPoint;
    double th = 1e-6;
    int ind = 0;
    std::vector<glm::f64vec3> P;
    for(auto&sa: SA){
        for(int n = 0; n < (int)sa.OriginalVertices.size(); n++){
            P.push_back(glm::f64vec3{X[3 * ind], X[3 * ind + 1], X[3 * ind + 2]});
            ind++;
        }
    }

    double f = 0.0;
    ind = 0;
    for(int i = 0; i < (int)SA.size(); i++)f += od->Edev(P, ICE, th);
    if(!grad.empty()){
        for(auto& p: P){
            p.x += eps; double fp = od->Edev(P, ICE, th); p.x -= 2.0*eps; double fm = od->Edev(P, ICE, th);
            grad[ind] = (fp - fm)/(2.0*eps); p.x = X[ind++];
            p.y += eps; fp = od->Edev(P, ICE, th); p.y -= 2.0*eps; fm = od->Edev(P, ICE, th);
            grad[ind] = (fp - fm)/(2.0*eps); p.y = X[ind++];
            p.z += eps; fp = od->Edev(P, ICE, th); p.z -= 2.0*eps; fm = od->Edev(P, ICE, th);
            grad[ind] = (fp - fm)/(2.0*eps); p.z = X[ind++];
        }
    }
    std::cout<<"developability  = " << f << std::endl;
    return f;
}

double RevisionVertices::Const_Econv(const std::vector<double>& X, std::vector<double> &grad, void *f_data){
    SmoothSurface *od = (SmoothSurface*)f_data;
    std::vector<SmoothingArea> SA = od->SA;
    int ind = 0;
    std::vector<glm::f64vec3> P, Pori;

    for(auto sm: SA){
        glm::f64vec3 vt = glm::normalize(sm.stP.second->p3 - sm.stP.first->p3), vt2d = glm::normalize(sm.stP.second->p - sm.stP.first->p);
        glm::f64vec3 SpinAxis = (vt2d.y < 0)? glm::f64vec3{0,0,-1}: glm::f64vec3{0,0,1};
        glm::f64vec3 befP3 = sm.stP.first->p3, befP2 = sm.stP.first->p;
        for(int n = 0; n < (int)sm.OriginalVertices.size(); n++){
            glm::f64vec3 p = glm::f64vec3{X[3 * ind], X[3 * ind + 1], X[3 * ind + 2]};
            glm::f64vec3 e = p - befP3;
            double l = glm::length(e), phi = std::acos(glm::dot(vt, glm::normalize(e)));
            glm::f64vec3 p2d = l * glm::f64vec3{glm::rotate(phi, SpinAxis)* glm::f64vec4{ vt2d, 1.0}} + befP2;
            P.push_back(p2d); Pori.push_back(sm.OriginalVertices[n]->p);
            ind++;
            befP2 = p2d; befP3 = p;
            if(sm.qt != nullptr){
                vt =glm::normalize(sm.qt->p3 - p); vt2d =glm::normalize(sm.qt->p - p2d);
                SpinAxis = (vt2d.y < 0)? glm::f64vec3{0,0,-1}: glm::f64vec3{0,0,1};
            }
        }
    }
    double f = 0.0;
    f = od->Econv(P, Pori);
    if(!grad.empty()){
        for(int i = 0; i < (int)P.size(); i++){
            double fp, fm;
            P[i].x += eps; fp = od->Econv(P, Pori); P[i].x -= 2.0*eps; fm = od->Econv(P, Pori);
            grad[i] = (fp - fm)/(2.0*eps); P[i].x = X[i];
            P[i].y += eps; fp = od->Econv(P, Pori); P[i].y -= 2.0*eps; fm = od->Econv(P, Pori);
            grad[i] = (fp - fm)/(2.0*eps); P[i].y = X[i];
            P[i].z += eps; fp = od->Econv(P, Pori); P[i].z -= 2.0*eps; fm = od->Econv(P, Pori);
            grad[i] = (fp - fm)/(2.0*eps); P[i].z = X[i];
        }
    }
    std::cout<<"convexity  = " << f << std::endl;
    return f;
}

double RevisionVertices::Const_Eplanarity(const std::vector<double>& X, std::vector<double> &grad, void *f_data){
    double f = .0, fp, fm;
    auto Normal_Sim = [](std::vector<Vertex4d>& FC, std::vector<glm::f64vec3>& R){
        double f = 0.0;
        int ind = 0;
        for(auto itr = FC.begin(); itr != FC.end(); itr++){
            if(!itr->IsCalc){
                auto _next = itr + 1;
                glm::f64vec3 e = _next->first->p3 - itr->first->p3; e /= glm::length(e);
                glm::f64vec3 N = glm::cross(R[ind], e); N /= glm::length(N);
                glm::f64vec3 N2 = glm::cross(-e, R[ind+1]); N2 /= glm::length(N2);
                f += 1.0 - glm::dot(N, N2);
                ind++;
            }
            if(ind == (int)R.size() - 1)break;
        }
        return f;
    };
    std::vector<Vertex4d> *FC = (std::vector<Vertex4d>*)f_data;
    std::vector<glm::f64vec3> Rulings;

    int x_ind = 0;
    for(int j = 1; j < (int)FC[0].size(); j++){
        if(!FC[0][j].IsCalc){
            glm::f64vec3 e = (FC[0][j].first->p3 - FC[0][j-1].first->p3);
            glm::f64vec3 N = (glm::cross(e, FC[0][j+1].first->p3 - FC[0][j-1].first->p3));
            glm::f64vec3 r = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2*x_ind], X[2*x_ind + 1]);
            Rulings.push_back(r);
            x_ind++;
        }
    }
    f = Normal_Sim(FC[0], Rulings);
    if(!grad.empty()){
        x_ind = 0;
        for(int j = 1; j < (int)FC[0].size(); j++){
            if(!FC[0][j].IsCalc){
                glm::f64vec3 e = (FC[0][j].first->p3 - FC[0][j-1].first->p3);
                glm::f64vec3 N = (glm::cross(e, FC[0][j+1].first->p3 - FC[0][j-1].first->p3));

                Rulings[x_ind] = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2*x_ind] + eps, X[2 * x_ind + 1]);
                fp = Normal_Sim(FC[0], Rulings);
                Rulings[x_ind] = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2*x_ind] - eps, X[2 * x_ind + 1]);
                fm = Normal_Sim(FC[0], Rulings);
                grad[2 * x_ind] = (fp - fm)/(2.0*eps);

                Rulings[x_ind] = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2 * x_ind], X[2 * x_ind + 1] + eps);
                fp = Normal_Sim(FC[0], Rulings);
                Rulings[x_ind] = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2 * x_ind], X[2 * x_ind + 1] - eps);
                fm = Normal_Sim(FC[0], Rulings);
                grad[2 * x_ind + 1] = (fp - fm)/(2.0*eps);
                Rulings[x_ind] = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2 * x_ind], X[2 * x_ind + 1]);
                x_ind++;
            }
        }
    }
    //std::cout <<"E_planarity " << f << std::endl;
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

double RevisionVertices::E_normal(std::vector<glm::f64vec3>& R, const std::vector<Vertex4d>& FC){
    std::vector<int> Vertices_Indicator;
    for(int i = 0; i < (int)FC.size(); i++){if(FC[i].IsCalc)Vertices_Indicator.push_back(i);}
    double f = .0;
    int j = 0;
    for(int k = 0; k < (int)Vertices_Indicator.size() - 1; k++){
        for(int i = Vertices_Indicator[k] + 1; i < Vertices_Indicator[k+1]; i++){
            glm::f64vec3 N = glm::normalize(glm::cross(glm::normalize(FC[i-1].first->p3 - FC[i].first->p3), R[j]));
            glm::f64vec3 N_ori = glm::normalize(glm::cross(glm::normalize(FC[i-1].first->p3 - FC[i].first->p3), glm::normalize(FC[i].second->p3 - FC[i].first->p3)));
            f += 1.0 - glm::dot(N, N_ori);
            N = glm::normalize(glm::cross(R[j], glm::normalize(FC[i+1].first->p3 - FC[i].first->p3)));
            N_ori = glm::normalize(glm::cross(glm::normalize(FC[i].second->p3 - FC[i].first->p3), glm::normalize(FC[i+1].first->p3 - FC[i].first->p3)));
            f += 1.0 - glm::dot(N, N_ori);
            j++;
        }
    }

    return f;
}

double RevisionVertices::Minimize_SmoothSrf(const std::vector<double>& X, std::vector<double> &grad, void *f_data){
    SmoothSurface *od = (SmoothSurface*)f_data;
    std::vector<SmoothingArea> SA = od->SA;
    std::vector<std::vector<glm::f64vec3>> P,Pori, P2d;
    int i = 0;
    for(auto&sa: SA){
        std::vector<glm::f64vec3> tmp = {sa.stP.first->p3}, tmp_ori = {sa.stP.first->p3_ori};
        glm::f64vec3 vt = glm::normalize(sa.stP.second->p3 - sa.stP.first->p3), vt2d = glm::normalize(sa.stP.second->p - sa.stP.first->p);
        glm::f64vec3 SpinAxis = (vt2d.y < 0)? glm::f64vec3{0,0,-1}: glm::f64vec3{0,0,1};
        std::vector<glm::f64vec3> _P2d;
        glm::f64vec3 befP3 = sa.stP.first->p3, befP2 = sa.stP.first->p;

        for(auto&p: sa.OriginalVertices){
            glm::f64vec3 p3d = glm::f64vec3{X[i], X[i+1], X[i+2]};
            tmp_ori.push_back(p->p3_ori);
            tmp.push_back(p3d);
            i += 3;
            glm::f64vec3 e = p3d - befP3;
            double l = glm::length(e), phi = std::acos(glm::dot(vt, glm::normalize(e)));
            glm::f64vec3 p2d = l * glm::f64vec3{glm::rotate(phi, SpinAxis)* glm::f64vec4{ vt2d, 1.0}} + befP2;
            _P2d.push_back(p2d);
            befP2 = p2d; befP3 = p3d;
            if(sa.qt != nullptr){
                vt =glm::normalize(sa.qt->p3 - p3d); vt2d =glm::normalize(sa.qt->p - p2d);
                SpinAxis = (vt2d.y < 0)? glm::f64vec3{0,0,-1}: glm::f64vec3{0,0,1};
            }
        }tmp.push_back(sa.lastP.first->p3); tmp_ori.push_back(sa.lastP.first->p3_ori);
        P.push_back(tmp); Pori.push_back(tmp_ori); P2d.push_back(_P2d);
    }

    double f_sim, f_fair, f_iso;
    double w_sim = 10, w_fair = 0.1, w_iso = 0.;
    if(!grad.empty()){
        std::vector<double> dE_sim(grad.size()), dE_fair(grad.size()), dE_iso(grad.size());
        f_sim = RevisionVertices::E_sim(P, Pori, dE_sim); f_fair = RevisionVertices::E_fair(P2d, dE_fair), f_iso = RevisionVertices::E_iso(P, Pori, dE_iso);
        for(int i = 0; i < (int)grad.size(); i++)grad[i] = w_sim * dE_sim[i] + w_fair * dE_fair[i] + w_iso * dE_iso[i];

    }
    std::cout<< "similarilty = " << f_sim << "  ,  fairness = " << f_fair << ", isometric = " << f_iso << std::endl;
    return w_sim * f_sim + w_fair * f_fair + w_iso * f_iso;
}

double RevisionVertices::Minimize_PlanaritySrf(const std::vector<double>& X, std::vector<double> &grad, void *f_data){

    std::vector<Vertex4d> *FC = (std::vector<Vertex4d>*)f_data;
    std::vector<glm::f64vec3> Rulings;

    int x_ind = 0;
    for(int j = 1; j < (int)FC[0].size(); j++){

        if(!FC[0][j].IsCalc){
            glm::f64vec3 e = (FC[0][j].first->p3 - FC[0][j-1].first->p3);
            glm::f64vec3 N = (glm::cross(e, FC[0][j+1].first->p3 - FC[0][j-1].first->p3));
            glm::f64vec3 r = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2*x_ind], X[2*x_ind + 1]);
            Rulings.push_back(r);
            x_ind++;
        }
    }

    double f = 0.0, fp, fm;
    f = RevisionVertices::E_normal(Rulings, FC[0]);
    if(!grad.empty()){
        x_ind = 0;
        for(int j = 1; j < (int)FC[0].size(); j++){
            if(!FC[0][j].IsCalc){
                glm::f64vec3 e = (FC[0][j].first->p3 - FC[0][j-1].first->p3);
                glm::f64vec3 N = (glm::cross(e, FC[0][j+1].first->p3 - FC[0][j-1].first->p3));

                Rulings[x_ind] = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2*x_ind] + eps, X[2 * x_ind + 1]);
                fp = RevisionVertices::E_normal(Rulings, FC[0]);
                Rulings[x_ind] = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2*x_ind] - eps, X[2 * x_ind + 1]);
                fm = RevisionVertices::E_normal(Rulings, FC[0]);
                grad[2 * x_ind] = (fp - fm)/(2.0*eps);

                Rulings[x_ind] = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2 * x_ind], X[2 * x_ind + 1] + eps);
                fp = RevisionVertices::E_normal(Rulings, FC[0]);
                Rulings[x_ind] = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2 * x_ind], X[2 * x_ind + 1] - eps);
                fm = RevisionVertices::E_normal(Rulings, FC[0]);
                grad[2 * x_ind + 1] = (fp - fm)/(2.0*eps);
                Rulings[x_ind] = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2 * x_ind], X[2 * x_ind + 1]);
                x_ind++;
            }
        }
    }
    //std::cout <<"E_normal " << f << std::endl;
    return f;
}

bool FoldLine::SimpleSmooothSrf(const std::vector<Vertex*>& Poly_v){

    auto getV = [](Vertex *o, Vertex *x, Vertex *o2, Vertex *x2)->Vertex*{
        if(IsParallel(o, x, o2, x2))return nullptr;
        glm::f64vec3 p2d = calcCrossPoint_2Vertex(o, x, o2, x2);
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
                p = MathTool::calcCrossPoint_2Vector(FoldingCurve[i].first->p, 1000.0 * r2d + FoldingCurve[i].first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
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
    }else FoldingCurve[0].second->p3 = calcTargetDistanceOnPlane(FoldingCurve[0].second->p, FoldingCurve[1].first, FoldingCurve[0].first, FoldingCurve[1].second);

    v_clst = getClosestVertex(FoldingCurve.back().second, FoldingCurve.back().first, FoldingCurve, false);
    if(v_clst != nullptr){
        int j;
        for(j = 0; j < (int)FoldingCurve.size(); j++){if(v_clst == FoldingCurve[j].second)break;}
        FoldingCurve.back().second->p3 = calcTargetDistanceOnPlane(FoldingCurve.back().second->p, FoldingCurve[j].first, FoldingCurve[j].second, FoldingCurve[j-1].second);
    }else FoldingCurve.back().second->p3 = calcTargetDistanceOnPlane(FoldingCurve.back().second->p, FoldingCurve.end()[-2].first, FoldingCurve.back().first, FoldingCurve.end()[-2].second);

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
    opt.add_inequality_constraint(RevisionVertices::Const_Econv, &od);
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
                    FoldingCurve[j].second->p = MathTool::calcCrossPoint_2Vector(FoldingCurve[j].first->p, 1000.0 * r2d + FoldingCurve[j].first->p, sm.stP.second->p, sm.lastP.second->p);
                    r3d = glm::normalize(sm.qt->p3 - FoldingCurve[j].first->p3);
                    if(glm::dot(r3d, glm::normalize(sm.stP.second->p3 - sm.stP.first->p3)) < 0)r3d *= -1;
                    vt =glm::normalize(sm.qt->p3 - FoldingCurve[j].first->p3); vt2d =glm::normalize(sm.qt->p - FoldingCurve[j].first->p);
                    SpinAxis = (vt2d.y < 0)? glm::f64vec3{0,0,-1}: glm::f64vec3{0,0,1};
                }
                for(int k = 0; k < (int)Poly_v.size(); k++){
                    p = MathTool::calcCrossPoint_2Vector(FoldingCurve[j].first->p, 1000.0 * r2d + FoldingCurve[j].first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
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
                    p = MathTool::calcCrossPoint_2Vector(FoldingCurve[j].first->p, 1000.0 * r2d + FoldingCurve[j].first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
                    if(MathTool::is_point_on_line(p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p) && glm::dot(r2d, p - FoldingCurve[j].first->p) > 0)break;
                }
                FoldingCurve[j].third->p = p;
                FoldingCurve[j].third->p3 = glm::length(p - FoldingCurve[j].first->p) * r3d + FoldingCurve[j].first->p3;

                auto v4d = FoldingCurve[j];
                glm::f64vec3 et = v4d.second->p3 -v4d.first->p3, er = FoldingCurve[j-1].first->p3 -v4d.first->p3, eb = v4d.third->p3 -v4d.first->p3, el = FoldingCurve[j+1].first->p3 -v4d.first->p3;
                et /= glm::length(et); er /= glm::length(er); eb /= glm::length(eb); el /= glm::length(el);
                double phi1 = std::acos(glm::dot(et, er)), phi2 = std::acos(glm::dot(et, el)), phi3 = std::acos(glm::dot(eb, el)), phi4 = std::acos(glm::dot(eb, er));
                std::cout << j << " : " << abs(2.0*std::numbers::pi - phi1 - phi2 - phi3 - phi4)  << " , " << glm::degrees(phi1) << " , " << glm::degrees(phi2) << " , " << glm::degrees(phi3) << ", " << glm::degrees(phi4) << std::endl;
            }
        }
        Vertex* v_clst = getClosestVertex(FoldingCurve[0].second, FoldingCurve[0].first, FoldingCurve, false);
        if(v_clst != nullptr){
            int j; for(j = 0; j < (int)FoldingCurve.size(); j++){if(v_clst == FoldingCurve[j].second)break;}
            FoldingCurve[0].second->p3 = calcTargetDistanceOnPlane(FoldingCurve[0].second->p, FoldingCurve[j].first, FoldingCurve[j].second, FoldingCurve[j+1].second);
        } else FoldingCurve[0].second->p3 = calcTargetDistanceOnPlane(FoldingCurve[0].second->p, FoldingCurve[1].second, FoldingCurve[0].first, FoldingCurve[1].first);
        v_clst = getClosestVertex(FoldingCurve.back().second, FoldingCurve.back().first, FoldingCurve, false);
        if(v_clst != nullptr){
            int j; for(j = 0; j < (int)FoldingCurve.size(); j++){if(v_clst == FoldingCurve[j].second)break;}
            FoldingCurve.back().second->p3 = calcTargetDistanceOnPlane(FoldingCurve.back().second->p, FoldingCurve[j].first, FoldingCurve[j].second, FoldingCurve[j-1].second);
        }else FoldingCurve.back().second->p3 = calcTargetDistanceOnPlane(FoldingCurve.back().second->p, FoldingCurve.end()[-2].second, FoldingCurve.end()[-2].first, FoldingCurve.back().first);
        std::cout<<"smoothing finish"<<std::endl;
    }catch (std::exception& e) {std::cout << "nlopt failed: " << e.what() << std::endl; }

    for(auto&sm: od.SA){
        delete sm.qb; delete sm.qt;
    }

    return res_qt;
}

/*
std::vector<std::vector<glm::f64vec3>> FoldLine::Optimization_SmooothSrf2(const std::vector<Vertex*>& Poly_v, bool IsConnectEndPoint){

    std::vector<std::vector<glm::f64vec3>> res_qt;
    Vertex4d bef = FoldingCurve.front();
    std::vector<RevisionVertices::SmoothingArea> SA;
    std::vector<Vertex*> OriginalVertices;
    int bef_ind;
    for(auto itr = FoldingCurve.begin()+1; itr != FoldingCurve.end() - 1; itr++){
        if(itr->IsCalc)continue;
        std::vector<double> X;
        X.push_back(itr->first->p3_ori.x); X.push_back(itr->first->p3_ori.y); X.push_back(itr->first->p3_ori.z);
        auto itr2 = itr + 1;
        nlopt::opt opt;
        while(!itr2->IsCalc)itr2++;
        opt = nlopt::opt(nlopt::LD_MMA, X.size());
        opt.set_min_objective(RevisionVertices::Minimize_SmoothSrf, &od);
        opt.add_inequality_constraint(RevisionVertices::Const_Edev, &od);
        opt.add_inequality_constraint(RevisionVertices::Const_Econv, &od);
        opt.set_xtol_rel(1e-13);
    }
    RevisionVertices::ObjData_smooth od = {SA, IsConnectEndPoint};
    if(X.empty())return res_qt;



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
                    FoldingCurve[j].second->p = MathTool::calcCrossPoint_2Vector(FoldingCurve[j].first->p, 1000.0 * r2d + FoldingCurve[j].first->p, sm.stP.second->p, sm.lastP.second->p);
                    r3d = glm::normalize(sm.qt->p3 - FoldingCurve[j].first->p3);
                    if(glm::dot(r3d, glm::normalize(sm.stP.second->p3 - sm.stP.first->p3)) < 0)r3d *= -1;
                    vt =glm::normalize(sm.qt->p3 - FoldingCurve[j].first->p3); vt2d =glm::normalize(sm.qt->p - FoldingCurve[j].first->p);
                    SpinAxis = (vt2d.y < 0)? glm::f64vec3{0,0,-1}: glm::f64vec3{0,0,1};
                }
                for(int k = 0; k < (int)Poly_v.size(); k++){
                    p = MathTool::calcCrossPoint_2Vector(FoldingCurve[j].first->p, 1000.0 * r2d + FoldingCurve[j].first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
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
                    p = MathTool::calcCrossPoint_2Vector(FoldingCurve[j].first->p, 1000.0 * r2d + FoldingCurve[j].first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
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
        } else FoldingCurve[0].second->p3 = calcTargetDistanceOnPlane(FoldingCurve[0].second->p, FoldingCurve[1].second, FoldingCurve[0].first, FoldingCurve[1].first);
        v_clst = getClosestVertex(FoldingCurve.back().second, FoldingCurve.back().first, FoldingCurve, false);
        if(v_clst != nullptr){
            int j; for(j = 0; j < (int)FoldingCurve.size(); j++){if(v_clst == FoldingCurve[j].second)break;}
            FoldingCurve.back().second->p3 = calcTargetDistanceOnPlane(FoldingCurve.back().second->p, FoldingCurve[j].first, FoldingCurve[j].second, FoldingCurve[j-1].second);
        }else FoldingCurve.back().second->p3 = calcTargetDistanceOnPlane(FoldingCurve.back().second->p, FoldingCurve.end()[-2].second, FoldingCurve.end()[-2].first, FoldingCurve.back().first);
        std::cout<<"smoothing finish"<<std::endl;
    }catch (std::exception& e) {std::cout << "nlopt failed: " << e.what() << std::endl; }

    for(auto&sm: od.SA){
        delete sm.qb; delete sm.qt;
    }

    return res_qt;
}
*/
std::vector<std::vector<glm::f64vec3>> FoldLine::Optimization_PlanaritySrf(const std::vector<Vertex*>& Poly_v){
    std::vector<double> X;
    std::vector<std::vector<glm::f64vec3>> res_qt;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){
        if(!FoldingCurve[i].IsCalc){
           glm::f64vec3 e = glm::normalize(FoldingCurve[i-1].first->p3 - FoldingCurve[i].first->p3);
           glm::f64vec3 N = glm::normalize(glm::cross(e, FoldingCurve[i+1].first->p3 - FoldingCurve[i].first->p3));
           glm::f64vec3 r = glm::normalize(FoldingCurve[i].second->p3 - FoldingCurve[i].first->p3);
           double phi = (glm::dot(e,r) <= -1)? std::numbers::pi: (glm::dot(e,r) >= 1)? 0: std::acos(glm::dot(e, r));
           r = MathTool::ProjectionVector(r, e, true); N = MathTool::ProjectionVector(N, e, true);
           double a = (glm::dot(N,r) <= -1)? std::numbers::pi: (glm::dot(N,r) >= 1)? 0: std::acos(glm::dot(N, r));
           if(glm::dot(glm::cross(r, N), e) < 0){
               a = std::numbers::pi/2.0 + a;
           }
           else{
               a = (std::numbers::pi/2.0 - a > 0)? 2.0*std::numbers::pi - a: std::numbers::pi/2.0 - a;
           }
           //a = 2.0 * std::numbers::pi - a;
           std::cout << glm::degrees(a)<<std::endl;
           X.push_back(a); X.push_back(phi);
        }
    }

    if(X.empty())return res_qt;
    nlopt::opt opt;
    opt = nlopt::opt(nlopt::LD_MMA, X.size());
    opt.set_min_objective(RevisionVertices::Minimize_PlanaritySrf, &(FoldingCurve));
    opt.add_inequality_constraint(RevisionVertices::Const_Eplanarity, &FoldingCurve);
    opt.set_xtol_rel(1e-13);

    double minf;
    try {
        nlopt::result result = opt.optimize(X, minf);
        std::cout <<"result :  planarity  " <<result << std::endl;
        std::cout << "found minimum at f = "  << std::setprecision(10) << minf << std::endl;
        int x_ind = 0;
        glm::f64vec3 p;
        for(auto itr = FoldingCurve.begin(); itr != FoldingCurve.end(); itr++){
            if(!itr->IsCalc){
                auto _prev = itr-1, _next = itr + 1;
                std::cout << std::distance(FoldingCurve.begin(), itr) << ": "  << glm::degrees(X[2 * x_ind]) << " , " << glm::degrees(X[2 * x_ind + 1]) << std::endl;
                glm::f64vec3 e = _prev->first->p3 - itr->first->p3;
                glm::f64vec3 N = glm::cross(e, _next->first->p3 - itr->first->p3);
                glm::f64vec3 r3d = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2 * x_ind], X[2 * x_ind + 1]);
                glm::f64vec3 r2d = glm::normalize(glm::rotate(X[2 * x_ind + 1], glm::f64vec3{0,0,-1}) * glm::f64vec4{_prev->first->p - itr->first->p, 1});
                for(int k = 0; k < (int)Poly_v.size(); k++){
                    p = MathTool::calcCrossPoint_2Vector(itr->first->p, 1000.0 * r2d + itr->first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
                    if(MathTool::is_point_on_line(p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p) && glm::dot(r2d, p - itr->first->p) > 0)break;
                }
                itr->second->p = p;
                std::cout << glm::length(itr->second->p3 - glm::length(p - itr->first->p) * r3d + itr->first->p3) << std::endl;
                itr->second->p3 = glm::length(p - itr->first->p) * r3d + itr->first->p3;
                x_ind++;
            }
        }


    }catch (std::exception& e) {std::cout << "nlopt failed: " << e.what() << std::endl; }
    std::cout <<"planarity optimization finish"<<std::endl;
    return res_qt;
}

bool FoldLine::Optimization_FlapAngle(std::vector<Line*>& Rulings, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, double wb, double wp, bool ConstFunc){
    if(FoldingCurve.empty())return false;

    std::vector<int> Vertices_Ind;for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}

    double a = 0, a2 = 0.0;
    glm::f64vec3 e = glm::normalize(FoldingCurve[Vertices_Ind[0]].first->p3 - FoldingCurve[Vertices_Ind[1]].first->p3),
            e2 = glm::normalize(FoldingCurve[Vertices_Ind[2]].first->p3 - FoldingCurve[Vertices_Ind[1]].first->p3);
    glm::f64vec3 SpinAxis = glm::normalize(FoldingCurve[Vertices_Ind[1]].third->p3 - FoldingCurve[Vertices_Ind[1]].first->p3);
    glm::f64vec3 Nb = -glm::normalize(glm::cross(e, e2)), N4 = glm::normalize(glm::cross(e, SpinAxis));

    double a_ll = std::atan2(glm::dot(glm::cross(Nb,N4), e),glm::dot(-Nb, N4));
    if(a_ll < 0)a_ll += 2.0*std::numbers::pi;
    double a_con = a_ll + std::numbers::pi;
    if(a_con > 2.0*std::numbers::pi)a_con -= 2.0*std::numbers::pi;
    if(a_con < 0)a_con +=2.0*std::numbers::pi;

    double phi3 = std::acos(glm::dot(e2,SpinAxis)), phi4 = std::acos(glm::dot(e,SpinAxis));
    double k = 2.0 * std::numbers::pi - phi3 - phi4;
    double a_min, a_max;

    glm::f64vec3 FaceNp = glm::normalize(glm::cross(SpinAxis, e2)), CrossV = glm::normalize(glm::cross(N4, FaceNp));
    double val = glm::dot(SpinAxis, CrossV);
    bool IsMount = (glm::dot(SpinAxis, CrossV) < 0)? true: false;

    if(k < std::numbers::pi && IsMount){a_min = a_con + std::numbers::pi; a_max = 2.0 * std::numbers::pi;}
    if(k >= std::numbers::pi && IsMount){a_min = std::numbers::pi; a_max = std::numbers::pi + a_con;}
    if(k < std::numbers::pi && !IsMount){a_min = 0.0;a_max = a_con - std::numbers::pi;}
    if(k >= std::numbers::pi && !IsMount){a_min = a_con - std::numbers::pi; a_max = std::numbers::pi;}
    a = (a_min + a_max)/2.0;

    std::cout <<glm::degrees(a_ll) << " , " << glm::degrees(a_con)<<std::endl;
    std::ofstream ofs2;
    std::filesystem::create_directory("./Optimization");

    std::string AngleFile = "./Optimization/ChangeAngle.csv" ;
    ofs2.open(AngleFile, std::ios::out); ofs2 << "a(radian),a(degree),Eruling, Ebend, Eparalell, Eruling2"<<std::endl;
    while(a <= 2.0*std::numbers::pi){
         _FoldingAAAMethod(FoldingCurve, Poly_V, a);
         double f = RulingsCrossed(FoldingCurve);
         double fb = Fbend2(FoldingCurve);
         double fp =  Fparallel(FoldingCurve);
         double f2 = RulingsCrossed2(FoldingCurve);
         ofs2 << a << "," << a * 180.0/std::numbers::pi << ", " << f << ", " << fb << ", " << fp <<", " << f2 <<  std::endl;
        a += 1e-3;
    }ofs2.close();

    RevisionVertices::ObjData od = {FoldingCurve, Vertices, Poly_V};
    if(ConstFunc)od.AddWeight(wb, wp, -1); else od.AddWeight(wb, wp, 100.0);
    nlopt::opt opt;
    opt = nlopt::opt(nlopt::LD_MMA, 1);
    opt.set_min_objective(ObjFunc_FlapAngle, &od);
    opt.set_lower_bounds(a_min);
    opt.set_upper_bounds(a_max);
    if(ConstFunc)opt.add_inequality_constraint(Fruling, &od);
    //opt.set_param("inner_maxeval", 100);
    opt.set_maxtime(2.0);//stop over this time
    //opt.set_ftol_rel(1e-10);
    opt.set_xtol_rel(1e-13);
    //
    //a = a_min;
    std::cout << "area " << glm::degrees(a_min) << " < " << glm::degrees(a) << " < " << glm::degrees(a_max) << std::endl;

    AngleFile = "./Optimization/OptimArea.csv";
    ofs2.open(AngleFile, std::ios::out);
    ofs2 << "a(radian), a(degree) , Eruling, Ebend, Eparalell, Ebend, Eruling2"<<std::endl;
    for(double _a = a_min; _a <= a_max; _a+= 1e-3){
         _FoldingAAAMethod(FoldingCurve, Poly_V, _a);
         double f = RulingsCrossed(FoldingCurve);
         double fb = Fbend2(FoldingCurve);
          double fp =  Fparallel(FoldingCurve);
         double fr = RulingsCrossed2(FoldingCurve);
        ofs2 << _a << ", " << glm::degrees(_a) << " , " << f << ", " << fb << ", " <<fp <<", "<< fr << std::endl;
    }ofs2.close();

    double minf_amin, minf_amax, f_ruling_amin, f_ruling_amax;
    double res_amin, res_amax;
    ofs_Ebend.open(File_Ebend, std::ios::out); ofs_Eruling.open(File_Eruling, std::ios::out);
      try {
        std::vector<double> _a{a_min + 0.5};
          nlopt::result result = opt.optimize(_a, minf_amin);

          std::cout <<"result :  lower bound" <<result << std::endl;
          f_ruling_amin = (ConstFunc)? RulingsCrossed(FoldingCurve): RulingsCrossed2(FoldingCurve);
          res_amin = _a[0];
          std::cout << "found minimum at f(" << glm::degrees(_a[0]) << ") = " << std::setprecision(10) << minf_amin << "  ,  " << f_ruling_amin <<  std::endl;          
      }
      catch (std::exception& e) {
          std::cout << "nlopt failed: " << e.what() << std::endl;
      }ofs_Ebend.close(); ofs_Eruling.close();

    try {
      std::vector<double> _a{a_max - 0.5};
        nlopt::result result = opt.optimize(_a, minf_amax);
        f_ruling_amax = (ConstFunc)? RulingsCrossed(FoldingCurve): RulingsCrossed2(FoldingCurve);
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
   std::cout << "result : smaller = " << glm::degrees(a2) << ",  f_bend = " << Fbend2(FoldingCurve) << "  ,  f_paralell = " << Fparallel(FoldingCurve) << std::endl;
    std::cout << "finish"<<std::endl;
    return true;
}

inline double RevisionVertices::getK(const glm::f64vec3 o, const glm::f64vec3 x, const glm::f64vec3 x2){
    double k = std::acos(glm::dot(glm::normalize(x - o), glm::normalize(x2 - o)));
    if(glm::dot(glm::f64vec3{0,0,1}, glm::cross(glm::normalize(x - o), glm::normalize(x2 - o))) > 0)k = 2.0*std::numbers::pi - k;
    return k;
}

void FoldLine::SimplifyModel(double tol){
    for(auto&fc: FoldingCurve){
        fc.IsCalc = true;
        fc.first->p = fc.first->p2_ori; fc.first->p3 = fc.first->p3_ori;
        fc.second->p = fc.second->p2_ori; fc.second->p3 = fc.second->p3_ori;
        fc.third->p = fc.third->p2_ori; fc.third->p3 = fc.third->p3_ori;
    }
    auto tmp = TrimPoints2(FoldingCurve, tol);
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
            FoldingCurve[j].first->p = calcCrossPoint_2Vertex(FoldingCurve[j].second, FoldingCurve[j].third, FoldingCurve[Vertices_Ind[i]].first, FoldingCurve[Vertices_Ind[i+1]].first);
            FoldingCurve[j].first->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].first->p, FoldingCurve[ind_root].first, FoldingCurve[ind_branch].third, FoldingCurve[ind_branch].second);
            FoldingCurve[j].third->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].third->p, FoldingCurve[ind_root].first, FoldingCurve[ind_branch].first, FoldingCurve[ind_branch].second);
            FoldingCurve[j].second->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].second->p, FoldingCurve[ind_root].first, FoldingCurve[ind_branch].first, FoldingCurve[ind_branch].second);

        }
    }
}

bool FoldLine::RevisionCrosPtsPosition(){
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
        auto ReviseEndPosition = [](glm::f64vec3 e,  double k)->glm::f64vec3{
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
        if(abs(k0 - k1) > abs(k1 - k2)){
            std::cout <<"front"<<std::endl;
            FoldingCurve[0].first->p = ReviseEndPosition(FoldingCurve[2].first->p - FoldingCurve[1].first->p,  -(2.*k1 - k2));
            glm::f64vec3 newPos = getCrossPosition(FoldingCurve[0].first, FoldingCurve[1].first, FoldingCurve[0].first->vo, FoldingCurve[0].first->ve);
            FoldingCurve[0].first->p = newPos;
            FoldingCurve[0].first->p3 = FoldingCurve[0].first->vo->p3 + glm::length(newPos - FoldingCurve[0].first->vo->p)/glm::length(FoldingCurve[0].first->vo->p - FoldingCurve[0].first->ve->p)
                    * (FoldingCurve[0].first->ve->p3 - FoldingCurve[0].first->vo->p3);
        }
        k0 = RevisionVertices::getK(FoldingCurve.end()[-2].first->p, FoldingCurve.end()[-3].first->p, FoldingCurve.back().first->p);
        k1 = RevisionVertices::getK(FoldingCurve.end()[-3].first->p, FoldingCurve.end()[-4].first->p, FoldingCurve.end()[-2].first->p);
        k2 = RevisionVertices::getK(FoldingCurve.end()[-4].first->p, FoldingCurve.end()[-5].first->p, FoldingCurve.end()[-3].first->p);
        std::cout <<"back " << abs(k0 - k1) << ", " << abs(k1 - k2) << ", " << 2.*k1 - k2 <<  std::endl;
        if(abs(k0 - k1) > abs(k1 - k2)){
            std::cout << "back"<<std::endl;
            FoldingCurve.back().first->p = ReviseEndPosition(FoldingCurve.end()[-3].first->p - FoldingCurve.end()[-2].first->p, 2.*k1 - k2);
            glm::f64vec3 newPos2 = getCrossPosition(FoldingCurve.back().first, FoldingCurve.end()[-2].first, FoldingCurve.back().first->vo, FoldingCurve.back().first->ve);
            FoldingCurve.back().first->p = newPos2;
            FoldingCurve.back().first->p3 = FoldingCurve.back().first->vo->p3 + glm::length(newPos2 - FoldingCurve.back().first->vo->p)/
                    glm::length(FoldingCurve.back().first->vo->p - FoldingCurve.back().first->ve->p) * (FoldingCurve.back().first->ve->p3 - FoldingCurve.back().first->vo->p3);
        }

    }
    return true;

}

void FoldLine::ReassignColor(std::vector<Line*>& Rulings, ColorPoint& CP){
    //Y字のとき折り線は山であると仮定(谷の場合は色を反転)
    typedef struct {
        Line *r;
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
        for(auto&r: Rulings){
            if(r->et == EdgeType::r && MathTool::is_point_on_line(FoldingCurve[Vertices_Ind[i]].first->p, r->v->p, r->o->p)){
                int mv = (r->color == 0)? -1: (r->color > 0)? 0: 1;
                int type_mvk = (mv == 0 && k < std::numbers::pi)? 0: (mv == 0 && k >= std::numbers::pi)? 1: (mv == 1 && k < std::numbers::pi)? 2: (mv == 1 && k >= std::numbers::pi)? 3: -1;
                InitState.push_back(MVK(r, type_mvk));
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
                    if(IS.type_mvk == 1){IS.r->color *= -1; IS.type_mvk = 3; MV[2]--; MV[4]++;}
                    if(IS.type_mvk == 2){IS.r->color *= -1; IS.type_mvk = 0; MV[3]--; MV[1]++;}
                }else{
                    if(IS.type_mvk == 0){IS.r->color *= -1; IS.type_mvk = 2; MV[1]--; MV[3]++;}
                    if(IS.type_mvk == 3){IS.r->color *= -1; IS.type_mvk = 1; MV[4]--; MV[2]++;}
                }
            }
        }
    }
    std::cout <<"color changed"<<std::endl;
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
     return;
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
    return RevisionVertices::decideRulingDirectionOn3d(e, glm::cross(e,e2), a, Phi[0]);//新しいruling方向
}

void CalcRuling(double a, Vertex4d& xbef, Vertex4d& x, Vertex4d& xnext, std::vector<Vertex*>& Poly_V,  double& a2, glm::f64vec3& SrfN){

    glm::f64vec3 e = glm::normalize(xbef.first->p3 - x.first->p3), e2 = glm::normalize(xnext.first->p3 - x.first->p3);
    double beta;
    std::vector<double> Phi;
    glm::f64vec3 Axis = glm::normalize(x.third->p3 - x.first->p3);
    glm::f64vec3 r3d = _calcruling3d(a, e, e2, Axis, beta, Phi);
    glm::f64vec3 r2d = glm::rotate(Phi[0], glm::f64vec3{0,0,-1.0})* glm::f64vec4{(glm::normalize(xbef.first->p- x.first->p)), 1.0};r2d = glm::normalize(r2d);//展開図のruling方向

    glm::f64vec3 crossPoint;
    for(int k = 0; k < (int)Poly_V.size(); k++){
        crossPoint = MathTool::calcCrossPoint_2Vector(x.first->p, 1000.0 * r2d + x.first->p, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p);
        if(MathTool::is_point_on_line(crossPoint, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) && glm::dot(crossPoint - x.first->p, r2d) > 0)break;
    }
    x.second->p = crossPoint;
    x.second->p3 = glm::length(crossPoint - x.first->p) * r3d + x.first->p3;

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
            if(!IsParallel(p.first, q, p2.first, q2))vec = glm::normalize(calcCrossPoint_2Vertex(p.first, q, p2.first, q2) - V.first->p);
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
            q = calcCrossPoint_2Vertex(V.first, Qv, Poly_V[k], Poly_V[(k + 1) % (int)Poly_V.size()]);
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
                auto p = calcTargetDistanceOnPlane(FoldingCurve[Vertices_Ind[0]].second->p, v_clst, FoldingCurve[Vertices_Ind[j]].first, FoldingCurve[Vertices_Ind[j+1]].first);
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
std::vector<Vertex4d> TrimPoints2(std::vector<Vertex4d>& FoldingCurve, double tol){
    std::vector<Vertex4d> res = FoldingCurve;
    res.clear();
    size_t st = FoldingCurve.size();

    while(1){
        Douglas_Peucker_algorithm(FoldingCurve,res, tol);
        if(res.size() == st)break;
        st = res.size();
    }
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
