#include "foldline.h"

const double eps = 1e-7;
std::string File_Ebend = "./Optimization/Ebend.csv";
std::string File_Eruling = "./Optimization/Eruling.csv";
std::ofstream ofs_Ebend, ofs_Eruling;

inline Eigen::Vector3d _calcruling3d(const double& a, Eigen::Vector3d e, Eigen::Vector3d e2, Eigen::Vector3d axis, double& beta, std::vector<double>& Phi);
void CalcRuling(double a, std::shared_ptr<Vertex4d>& xbef, std::shared_ptr<Vertex4d>& x, std::shared_ptr<Vertex4d>& xnext, const std::vector<std::shared_ptr<Vertex>>& Poly_V, double& a2, Eigen::Vector3d& SrfN, const Eigen::Vector3d& SpinAixs);
inline double update_flapangle(double a, const Eigen::Vector3d& befN, const Eigen::Vector3d& SrfN, const Eigen::Vector3d& e);
void _FoldingAAAMethod(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle);
void _FoldingAAAMethod_center(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle);
std::vector<std::shared_ptr<Vertex4d>> TrimPoints(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, double tol);
bool IsRulingCrossed(Eigen::Vector3d N, Eigen::Vector3d& cp, Eigen::Vector3d& crossPoint,  const std::vector<std::shared_ptr<Vertex>>& Poly_V);
void Douglas_Peucker_algorithm(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, std::vector<std::shared_ptr<Vertex4d>>& res, double tol = std::numbers::pi/9.0);
std::shared_ptr<Vertex> getClosestVertex(const std::shared_ptr<Vertex>& v, const std::shared_ptr<Vertex>& o,  const std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, bool SkipTrimedPoint = true);

inline Eigen::Vector3d calcCrossPoint_2Vertex(const std::shared_ptr<Vertex>& p1, const std::shared_ptr<Vertex>& q1, const std::shared_ptr<Vertex>& p2, const std::shared_ptr<Vertex>& q2){
    Eigen::Matrix2d A; Eigen::Vector2d b;
    Eigen::Vector3d v1 = q1->p - p1->p, v2 = q2->p - p2->p;
    b(0) = p2->p.x() - p1->p.x(); b(1) = p2->p.y() - p1->p.y();
    A(0,0) = v1.x(); A(0,1) = -v2.x();
    A(1,0) = v1.y(); A(1,1) = -v2.y();
    return MathTool::calcCrossPoint_2Vector(p1->p, q1->p, p2->p, q2->p);
}

inline bool IsParallel(const std::shared_ptr<Vertex>& p1, const std::shared_ptr<Vertex>& q1, const std::shared_ptr<Vertex>& p2, const std::shared_ptr<Vertex>& q2){
    auto v1 = (q1->p - p1->p).normalized(), v2 = (q2->p - p2->p).normalized();
    return (abs(v1.dot(v2)) >= 1.0 - 1e-9)? true: false;
}

inline Eigen::Vector3d calcTargetDistanceOnPlane(Eigen::Vector3d p, const std::shared_ptr<Vertex>& o, const std::shared_ptr<Vertex>& v1, const std::shared_ptr<Vertex>& v2){
    Eigen::Matrix2d A; Eigen::Vector2d b;
    b(0) = p.x() - o->p.x(); b(1) = p.y() - o->p.y();
    A(0,0) = v1->p.x() - o->p.x(); A(0,1) = v2->p.x() - o->p.x();
    A(1,0) = v1->p.y() - o->p.y(); A(1,1) = v2->p.y() - o->p.y();
    Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
    return  x(0) * (v1->p3 - o->p3) + x(1) * (v2->p3 - o->p3) + o->p3;
}

namespace RevisionVertices{
    using FoldLine3d = std::vector<std::shared_ptr<Vertex4d>>;
    struct OptimizeParam{
        FoldLine3d FC;
        std::vector<std::shared_ptr<Vertex>> Poly_V;
        std::vector<double*> res_Fbend, res_Fruling, res_a;
        double wb, wp ,wr;
        void AddWeight(double _wb, double _wp, double _wr){ wb = _wb; wp = _wp; wr = _wr;}
        OptimizeParam(FoldLine3d& _FC,  const std::vector<std::shared_ptr<Vertex>>& _Poly_V):FC{_FC}, Poly_V{_Poly_V}{}
        ~OptimizeParam(){}
    private:
    };

    class OptimizeParam_v: public OptimizeParam{
    public:
        OptimizeParam_v() = delete;
        OptimizeParam_v(double _a, FoldLine3d& _FC, const std::vector<std::shared_ptr<Vertex>>& Poly_V): a(_a), OptimizeParam::OptimizeParam( _FC,  Poly_V){}
        ~OptimizeParam_v(){}
        double a;
    private:
    };

    class SmoothingArea{
    public:
        SmoothingArea() = delete;
        SmoothingArea(std::shared_ptr<Vertex4d>& a, std::shared_ptr<Vertex4d>& _end, int si, int li, std::vector<std::shared_ptr<Vertex>>& OV): stP{a}, lastP{_end}, st_ind{si}, last_ind{li}, OriginalVertices{OV}{
            qt = getV(stP->first, stP->second, lastP->first, lastP->second);
            qb = getV(stP->first, stP->third, lastP->first, lastP->third);
        }
        SmoothingArea(std::shared_ptr<Vertex4d>& a, std::shared_ptr<Vertex4d>& _end, int si, int li): stP{a}, lastP{_end}, st_ind{si}, last_ind{li}{
            qt = getV(stP->first, stP->second, lastP->first, lastP->second);
            qb = getV(stP->first, stP->third, lastP->first, lastP->third);
        }
        ~SmoothingArea(){}
        std::shared_ptr<Vertex4d> stP, lastP;
        std::vector<std::shared_ptr<Vertex>> OriginalVertices;
        int st_ind, last_ind;
        std::shared_ptr<Vertex> qt,qb;
    private:
        std::shared_ptr<Vertex> getV(const std::shared_ptr<Vertex>& o, const std::shared_ptr<Vertex>& x, const std::shared_ptr<Vertex>& o2, const std::shared_ptr<Vertex>& x2){
            if(IsParallel(o, x, o2, x2))return nullptr;
            Eigen::Vector3d p2d = calcCrossPoint_2Vertex(o, x, o2, x2);
            Eigen::Vector3d p3d = calcTargetDistanceOnPlane(p2d, o,  x, x2);
            return  std::make_shared<Vertex>(Vertex(p2d, p3d));
        }
    };
    struct SmoothSurface{
        std::vector<SmoothingArea> SA;
        bool IsConnectEndPoint;

        double Edev(std::vector<Eigen::Vector3d>& P, bool IsConnectEndPoint, const double th = 1e-5){
            auto Angle = [](Eigen::Vector3d& e, Eigen::Vector3d& e2)->double{return (e.dot(e2) >= 1)? 0: (e.dot(e2) <= -1)? std::numbers::pi: std::acos(e.dot(e2));  };
            double f = 0.0;
            int j = 0;
            std::vector<std::array<Eigen::Vector3d, 3>> X;
            Eigen::Vector3d rt, rb;
            for(const auto&sa: SA){
                rt = (sa.stP->second->p3 - sa.stP->first->p3).normalized();
                rb = (sa.stP->third->p3 - sa.stP->first->p3).normalized();
                X.push_back(std::array{sa.stP->first->p3, rt, rb});
                for(int i = 0; i < (int)sa.OriginalVertices.size(); i++){
                    rt = (sa.qt != nullptr)? (sa.qt->p3 - P[j]).normalized(): (sa.stP->second->p3 - sa.stP->first->p3).normalized();
                    rb = (sa.qb != nullptr)? (sa.qb->p3 - P[j]).normalized(): (sa.stP->third->p3 - sa.stP->first->p3).normalized();
                    X.push_back(std::array{P[j], rt, rb});
                    j++;
                }
            }
            rt = (SA.back().lastP->second->p3 - SA.back().lastP->first->p3).normalized();
            rb = (SA.back().lastP->third->p3 - SA.back().lastP->first->p3).normalized();
            X.push_back(std::array{SA.back().lastP->first->p3, rt, rb});

            Eigen::Vector3d el, er;
            for(int i = 1; i < (int)X.size() - 1; i++){
                el = (X[i+1][0] - X[i][0]).normalized(); er = (X[i-1][0] - X[i][0]).normalized();
                f += std::abs((2.0 * std::numbers::pi - (Angle(el, X[i][1]) + Angle(X[i][1], er) + Angle(er, X[i][2]) + Angle(X[i][2], el))) - th);
            }
            return f;
        }

        double Econv(std::vector<Eigen::Vector3d>& P, std::vector<Eigen::Vector3d>& Pori){
            auto SgdArea = [](Eigen::Vector3d& A, Eigen::Vector3d& B, Eigen::Vector3d& C)->double{
              return ((A.x() * B.y() - B.x() * A.y()) + (B.x() * C.y() - C.x() * B.y()) + (C.x() * A.y() - A.x() * C.y()))/2.0;
            };
            double f = 0.0;
            double sum = 0.0;
            if(P.empty())return 0.0;
            int j = 0;
            std::vector<Eigen::Vector3d> X, Xori;
            for(const auto&sa: SA){
                X.push_back(sa.stP->first->p); Xori.push_back(sa.stP->first->p);
                for(int i = 0; i < (int)sa.OriginalVertices.size(); i++){
                    X.push_back(P[j]); Xori.push_back(Pori[j]);
                    j++;
                }
            }
            X.push_back(SA.back().lastP->first->p); Xori.push_back(SA.back().lastP->first->p);
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
    using FoldCrv = std::vector<std::shared_ptr<Vertex4d>>;

    inline double getK(const Eigen::Vector3d o, const Eigen::Vector3d x, const Eigen::Vector3d x2);

    double Const_Edev(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    double Const_Econv(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    double Const_Eplanarity(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    double E_fair(std::vector<std::vector<Eigen::Vector3d>>& P, std::vector<double>& grad);
    double E_sim(std::vector<std::vector<Eigen::Vector3d>>& P, std::vector<std::vector<Eigen::Vector3d>>& Pori, std::vector<double>& grad);
    double E_iso(std::vector<std::vector<Eigen::Vector3d>>& P, std::vector<std::vector<Eigen::Vector3d>>& Pori, std::vector<double>& grad);
    double E_normal(std::vector<Eigen::Vector3d>& R, const std::vector<std::shared_ptr<Vertex4d>>& FC);
    double Minimize_PlanaritySrf(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    double Minimize_SmoothSrf(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    Eigen::Vector3d decideRulingDirectionOn3d(Eigen::Vector3d e, Eigen::Vector3d N, double a, double phi);

}

FoldLine::FoldLine(PaintTool _type)
{
    CtrlPts.clear();
    curveNum = 300;
    color = 0;
    type = _type;
    CurvePts.clear();
    a_flap = -1;
    tol = 0;
    validsize = 0;
}

void FoldLine::SortCurve(bool ascending){
    if(FoldingCurve.empty())return;
    if(!ascending)std::sort(FoldingCurve.begin(), FoldingCurve.end(), [](const std::shared_ptr<Vertex4d>& V1, const std::shared_ptr<Vertex4d>& V2){return V1->first->s > V2->first->s;});//左から右への曲線の流れにしたい
    else std::sort(FoldingCurve.begin(), FoldingCurve.end(), [](const std::shared_ptr<Vertex4d>& V1, const std::shared_ptr<Vertex4d>& V2){return V1->first->s < V2->first->s;});//左から右への曲線の流れにしたい
}

double FoldLine::getColor(){return color;}

bool FoldLine::delCtrlPt(Eigen::Vector3d& p, int dim, std::shared_ptr<OUTLINE>& outline){
    if(!outline->IsClosed()) return false;
    if(!CtrlPts.empty())CtrlPts.pop_back();
    double dist = 10;
    int ind = -1;
    for(int i = 0; i < (int)CtrlPts.size(); i++){
        if((p - CtrlPts[i]).norm() < dist){
            dist = (p - CtrlPts[i]).norm();
            ind = i;
        }
    }
    if(ind != -1) CtrlPts.erase(CtrlPts.begin() + ind);
    bool hasCrv = setCurve(dim);
    if(outline->IsClosed() && hasCrv){
        //bool res = applyCurvedFolding(dim, outline);
        //return res;
    }
}

bool FoldLine::moveCtrlPt(Eigen::Vector3d& p, int movePtIndex, int dim){
    if((int)CtrlPts.size() < movePtIndex || movePtIndex < 0)return false;
    CtrlPts[movePtIndex] = p;
    setCurve(dim);
    if(FoldingCurve.empty())return false;
    return true;

}

bool FoldLine::addCtrlPt(Eigen::Vector3d& p, int dim){
    //if((int)CtrlPts.size() <= dim)
    CtrlPts.push_back(p);
    bool hasCrv = setCurve(dim);
    return hasCrv;
}

//type == 0 line
//type == 1 b-spline
bool FoldLine::setCurve(int dim){
    using namespace MathTool;
    CurvePts.clear();
    if(type == PaintTool::FoldLine_line){
        if((int)CtrlPts.size() < 2)return false;
        while((int)CtrlPts.size() != 2){
            CtrlPts.erase(CtrlPts.end() - 1);
            CtrlPts.shrink_to_fit();
        }
        Eigen::Vector3d V = CtrlPts[1] - CtrlPts[0];
        for(int i = 0; i < curveNum; i++){
            CurvePts.push_back((double)i/(double)(curveNum - 1) * V + CtrlPts[0]);
        }
    }
    else if(type == PaintTool::FoldLine_arc){
        if((int)CtrlPts.size() < 3) return false;

        double l = (CtrlPts[0] - CtrlPts[1]).norm(), l2 = (CtrlPts[0], CtrlPts[2]).norm();
        CtrlPts[2] = l/l2 * (CtrlPts[2] - CtrlPts[0]) + CtrlPts[0];
        Eigen::Vector3d v = (CtrlPts[1] - CtrlPts[0]).normalized(), v2 = (CtrlPts[2] - CtrlPts[0]).normalized();
        double phi = acos(v.dot(v2));
        Eigen::Vector3d axis = v.cross(v2);
        Eigen::Translation3d T(CtrlPts[0]), invT(-CtrlPts[0]);// 平行移動ベクトルを作成
        for(int i = 0; i < curveNum; i++){
            Eigen::AngleAxisd R = Eigen::AngleAxisd(phi * (double)i/(double)curveNum, axis); // 回転行列を作成
            Eigen::Transform<double, 3, Eigen::Affine> transform = T * R * invT;// Transform行列を作成
            CurvePts.push_back(transform * CtrlPts[1]);
        }

    }
    else if(type == PaintTool::FoldLine_bezier){
        if((int)CtrlPts.size() <= dim)return false;
        if((int)CtrlPts.size() > dim + 1)CtrlPts.pop_back();
        double t = 0.0;
        Eigen::Vector3d v;
        for(int n = 0; n < curveNum; n++){
            v = Eigen::Vector3d(0,0,0);
            for (int i = 0; i < int(CtrlPts.size()); i++)  v += MathTool::BernsteinBasisFunc(dim, i, t) * CtrlPts[i];  
            CurvePts.push_back(v);
            t += 1/(double)curveNum;
        }
    }

    return true;
}

bool IsRulingCrossed (Eigen::Vector3d N, Eigen::Vector3d& cp, Eigen::Vector3d& crossPoint,  const std::vector<std::shared_ptr<Vertex>>& Poly_V){
    double l = 1000, minDist = 1000;
    bool IsIntersected = false;
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    N = N.normalized();
    int n = Poly_V.size();
    for(int i = 0; i < n; i++){
        Eigen::Vector3d v = Poly_V[(i+1) % n]->p - Poly_V[i]->p;
        A(0,0) = v.x(); A(0,1) = -l*N.x(); A(1,0) = v.y(); A(1,1) = -l * N.y();
        b(0) = cp.x() - Poly_V[i]->p.x(); b(1) = cp.y() - Poly_V[i]->p.y();
        if(abs(N.dot(v.normalized()) >= 1 - 1e-9))continue;
        Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
        if (0 <= x(1) && x(1) <= 1) {
            Eigen::Vector3d p = x(1) * l * N + cp;
            if((p - cp).norm() < minDist){
                crossPoint = p;
                minDist = (p - cp).norm();
            }
            IsIntersected = true;
        }
    }
    return IsIntersected;
}

double Fparallel(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve){
    std::vector<std::shared_ptr<Vertex4d>>ValidFC;
    double f = 0.0;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    for(int i = 1; i < (int)ValidFC.size() - 2; i++){
        Eigen::Vector3d v = (ValidFC[i]->second->p - ValidFC[i]->first->p).normalized(), v2 = (ValidFC[i+1]->second->p - ValidFC[i+1]->first->p).normalized();
        f += 1.0 - v.dot(v2);
    }
    return f;
}

double Fbend2(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve){
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    double f = 0.0;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    auto v1 = (ValidFC[0]->first->p3 - ValidFC[1]->first->p3).normalized(), v2 = (ValidFC[1]->second->p3 - ValidFC[1]->first->p3).normalized();
    Eigen::Vector3d Nt = (v1.cross(v2)).normalized();
    for(int i = 1; i < (int)ValidFC.size() - 1; i++){
        auto n = (ValidFC[i]->first->p3 - ValidFC[i+1]->first->p3);
        Eigen::Vector3d Ntp = (n.cross(ValidFC[i+1]->second->p3 - ValidFC[i+1]->first->p3)).normalized();
        double phi = ((Ntp.dot(Nt)) > 1)? std::numbers::pi: ((Ntp.dot(Nt)) < -1)? 0:  std::numbers::pi - std::acos(Ntp.dot(Nt));
        f += 1.0/(phi*phi);
        Nt = Ntp;
    }
    return f;
}

double RulingsCrossed(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve){
    double f = 0.0;
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    for(int i = 1; i < (int)ValidFC.size() -1; i++){
        for(int j = i; j < (int)ValidFC.size() -1; j ++){
            int i2 = j;
            if(i == j)continue;
            Eigen::Vector3d p1 = ValidFC[i]->first->p, q1 = ValidFC[i]->second->p, p2 = ValidFC[i2]->first->p, q2 = ValidFC[i2]->second->p;
            double t = ((p2.x() - p1.x())*(p2.y() - q2.y()) - (p2.x() - q2.x())*(p2.y() - p1.y()))/((q1.x() - p1.x()) * (p2.y() - q2.y()) - (p2.x() - q2.x())*(q1.y() - p1.y()));
            if(0 < t && t < 1) f += abs(1.0/t);
            continue;
            A(0,0) = ValidFC[i]->first->p.x() - ValidFC[i]->second->p.x();
            A(0,1) = -(ValidFC[i2]->first->p.x() - ValidFC[i2]->second->p.x());
            A(1,0) = ValidFC[i]->first->p.y() - ValidFC[i]->second->p.y();
            A(1,1) = -(ValidFC[i2]->first->p.y() - ValidFC[i2]->second->p.y());
            b(0) = ValidFC[i]->first->p.x() - ValidFC[i2]->first->p.x();
            b(1) = ValidFC[i]->first->p.y() - ValidFC[i2]->first->p.y();
            Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
            if(0 < ts(0) && ts(0) < 1)f += 1.0/ts(0);
        }
    }
    return f;
}

double RulingsCrossed_NoOutline(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve){
    double f = 0.0;
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    for(int i = 1; i < (int)ValidFC.size() -1; i++){
        for(int j = -1; j <= 1; j +=2){
            int i2 = i + j;
            if(i2 <= 0 || i2 >= (int)ValidFC.size() - 1)continue;
            Eigen::Vector3d p1 = ValidFC[i]->first->p, q1 = ValidFC[i]->second->p, p2 = ValidFC[i2]->first->p, q2 = ValidFC[i2]->second->p;
            double t = ((p2.x() - p1.x())*(p2.y() - q2.y()) - (p2.x() - q2.x())*(p2.y() - p1.y()))/((q1.x() - p1.x()) * (p2.y() - q2.y()) - (p2.x() - q2.x())*(q1.y() - p1.y()));
            //f += abs(1.0/t);
            //continue;
            A(0,0) = ValidFC[i]->first->p.x() - ValidFC[i]->second->p.x();
            A(0,1) = -(ValidFC[i2]->first->p.x() - ValidFC[i2]->second->p.x());
            A(1,0) = ValidFC[i]->first->p.y() - ValidFC[i]->second->p.y();
            A(1,1) = -(ValidFC[i2]->first->p.y() - ValidFC[i2]->second->p.y());
            b(0) = ValidFC[i]->first->p.x() - ValidFC[i2]->first->p.x();
            b(1) = ValidFC[i]->first->p.y() - ValidFC[i2]->first->p.y();
            Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
            f += abs(1.0/ts(0));
        }
    }
    return f;
}

bool FoldLine::isbend(){
    double fr = RulingsCrossed(FoldingCurve);
    return (fr < 1e-9 && a_flap != -1)? true: false;
}

double Fruling(const std::vector<double> &a, std::vector<double> &grad, void* f_data)
{
    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    RevisionVertices::FoldLine3d FoldingCurve = od->FC;
    double f = RulingsCrossed(FoldingCurve);
    if(!grad.empty()){
        std::vector<std::shared_ptr<Vertex>> Poly_V = od->Poly_V;
        _FoldingAAAMethod(FoldingCurve, Poly_V, a[0] + eps);
       double fp = RulingsCrossed(FoldingCurve);
        _FoldingAAAMethod(FoldingCurve, Poly_V, a[0] - eps);
       double fm = RulingsCrossed(FoldingCurve);
       grad[0] = (fp - fm)/(2.0 * eps);
       //for(int j = 0; j < (int)FoldingCurve.size(); j++)FoldingCurve[j]->second->p3 = tmp[j].p3;
    }
    if(ofs_Eruling.is_open()) ofs_Eruling << a[0] <<", " << MathTool::rad2deg(a[0]) << ", " <<  f <<  " ,  " << grad[0] << std::endl;
    if(DebugMode::Singleton::getInstance().isdebug())std::cout <<"constraint function = " << MathTool::rad2deg(a[0]) << "(" << a[0] << ")  , " << f << std::endl;
    return f;
 }

double ObjFunc_FlapAngle(const std::vector<double> &a, std::vector<double> &grad, void* f_data){
    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    RevisionVertices::FoldLine3d FoldingCurve = od->FC;
    std::vector<std::shared_ptr<Vertex>>Poly_V = od->Poly_V;
    _FoldingAAAMethod(FoldingCurve, Poly_V, a[0]);
    double fb =  Fbend2(FoldingCurve), fp = Fparallel(FoldingCurve);
    double fr = (od->wr != -1)?RulingsCrossed_NoOutline(FoldingCurve): 0;
    if(!grad.empty()){
        _FoldingAAAMethod(FoldingCurve, Poly_V, a[0] + eps);
       double fp = Fbend2(FoldingCurve);
       double fp2 = Fparallel(FoldingCurve);
       double fp3 = (od->wr != -1)?RulingsCrossed_NoOutline(FoldingCurve): 0;
        _FoldingAAAMethod(FoldingCurve, Poly_V, a[0] - eps);
       double fm = Fbend2(FoldingCurve);
       double fm2 = Fparallel(FoldingCurve);
       double fm3 = (od->wr != -1)?RulingsCrossed_NoOutline(FoldingCurve): 0;
       grad[0] = od->wb * (fp - fm)/(2.0 * eps) + od->wp * (fp2 - fm2)/(2.0 * eps) + 100.0*(fp3 - fm3)/(2.0*eps);
    }
    //std::cout <<"Fbend(" << glm::degrees(a[0]) << ")  =  " <<  fb <<  " ,  " << grad[0] << "  ,  Fparalell = " << fp << ", Fruling2 = "<< fr <<  std::endl;
    if(ofs_Ebend.is_open())ofs_Ebend << a[0] <<", " << MathTool::rad2deg(a[0]) << ", " <<  fb << ", " << fp <<   ",  " << grad[0] << std::endl;
    return od->wb * fb + od->wp * fp + 100.0 * fr;
}

Eigen::Vector3d RevisionVertices::decideRulingDirectionOn3d(Eigen::Vector3d e, Eigen::Vector3d N, double a, double phi){
    e = e.normalized();
    N = (Eigen::AngleAxisd(a, -e) * N).normalized(); // 回転行列を作成
    return (Eigen::AngleAxisd(phi, N) * e).normalized(); 
}

double RevisionVertices::Const_Edev(const std::vector<double>& X, std::vector<double> &grad, void *f_data){
    SmoothSurface *od = (SmoothSurface*)f_data;
    std::vector<SmoothingArea> SA = od->SA;
    bool ICE = od->IsConnectEndPoint;
    double th = 1e-6;
    int ind = 0;
    std::vector<Eigen::Vector3d> P;
    for(auto&sa: SA){
        for(int n = 0; n < (int)sa.OriginalVertices.size(); n++){
            P.push_back(Eigen::Vector3d(X[3 * ind], X[3 * ind + 1], X[3 * ind + 2]));
            ind++;
        }
    }

    double f = 0.0;
    ind = 0;
    for(int i = 0; i < (int)SA.size(); i++)f += od->Edev(P, ICE, th);
    if(!grad.empty()){
        for(auto& p: P){
            p.x() += eps; double fp = od->Edev(P, ICE, th); p.x() -= 2.0*eps; double fm = od->Edev(P, ICE, th);
            grad[ind] = (fp - fm)/(2.0*eps); p(0) = X[ind++];
            p.y() += eps; fp = od->Edev(P, ICE, th); p.y() -= 2.0*eps; fm = od->Edev(P, ICE, th);
            grad[ind] = (fp - fm)/(2.0*eps); p.y() = X[ind++];
            p.z() += eps; fp = od->Edev(P, ICE, th); p.z() -= 2.0*eps; fm = od->Edev(P, ICE, th);
            grad[ind] = (fp - fm)/(2.0*eps); p.z() = X[ind++];
        }
    }
    if(DebugMode::Singleton::getInstance().isdebug())std::cout<<"developability  = " << f << std::endl;
    return f;
}

double RevisionVertices::Const_Econv(const std::vector<double>& X, std::vector<double> &grad, void *f_data){
    SmoothSurface *od = (SmoothSurface*)f_data;
    std::vector<SmoothingArea> SA = od->SA;
    int ind = 0;
    std::vector<Eigen::Vector3d> P, Pori;

    for(auto sm: SA){
        Eigen::Vector3d vt = (sm.stP->second->p3 - sm.stP->first->p3).normalized(), vt2d = (sm.stP->second->p - sm.stP->first->p).normalized();
        Eigen::Vector3d SpinAxis = (vt2d.y() < 0) ? Eigen::Vector3d(0,0,-1): Eigen::Vector3d(0,0,1);
        Eigen::Vector3d befP3 = sm.stP->first->p3, befP2 = sm.stP->first->p;
        for(int n = 0; n < (int)sm.OriginalVertices.size(); n++){
            Eigen::Vector3d p = Eigen::Vector3d(X[3 * ind], X[3 * ind + 1], X[3 * ind + 2]);
            Eigen::Vector3d e = p - befP3;
            double l = e.norm(), phi = std::acos(vt.dot(e.normalized()));
            Eigen::Vector3d p2d = l * (Eigen::AngleAxisd(phi, SpinAxis)* vt2d) + befP2;
            P.push_back(p2d); Pori.push_back(sm.OriginalVertices[n]->p);
            ind++;
            befP2 = p2d; befP3 = p;
            if(sm.qt != nullptr){
                vt = (sm.qt->p3 - p).normalized(); vt2d = (sm.qt->p - p2d).normalized();
                SpinAxis = (vt2d.y() < 0)? Eigen::Vector3d(0,0,-1): Eigen::Vector3d(0,0,1);
            }
        }
    }
    double f = 0.0;
    f = od->Econv(P, Pori);
    if(!grad.empty()){
        for(int i = 0; i < (int)P.size(); i++){
            double fp, fm;
            P[i].x() += eps; fp = od->Econv(P, Pori); P[i].x() -= 2.0*eps; fm = od->Econv(P, Pori);
            grad[i] = (fp - fm)/(2.0*eps); P[i].x() = X[i];
            P[i].y() += eps; fp = od->Econv(P, Pori); P[i].y() -= 2.0*eps; fm = od->Econv(P, Pori);
            grad[i] = (fp - fm)/(2.0*eps); P[i].y() = X[i];
            P[i].z() += eps; fp = od->Econv(P, Pori); P[i].z() -= 2.0*eps; fm = od->Econv(P, Pori);
            grad[i] = (fp - fm)/(2.0*eps); P[i].z() = X[i];
        }
    }
    if(DebugMode::Singleton::getInstance().isdebug())std::cout<<"convexity  = " << f << std::endl;
    return f;
}

double RevisionVertices::Const_Eplanarity(const std::vector<double>& X, std::vector<double> &grad, void *f_data){
    double f = .0, fp, fm;
    auto Normal_Sim = [](std::vector<std::shared_ptr<Vertex4d>>& FC, std::vector<Eigen::Vector3d>& R){
        double f = 0.0;
        int ind = 0;
        for(auto itr = FC.begin(); itr != FC.end(); itr++){
            if(!(*itr)->IsCalc){
                auto _next = itr + 1;
                Eigen::Vector3d e = ((*_next)->first->p3 - (*itr)->first->p3).normalized();
                Eigen::Vector3d N = (R[ind].cross(e)).normalized();
                Eigen::Vector3d N2 = ((-e).cross(R[ind+1])).normalized(); 
                f += 1.0 - N.dot(N2);
                ind++;
            }
            if(ind == (int)R.size() - 1)break;
        }
        return f;
    };
    std::vector<std::shared_ptr<Vertex4d>> *FC = (std::vector<std::shared_ptr<Vertex4d>>*)f_data;
    std::vector<Eigen::Vector3d> Rulings;
    int x_ind = 0;
    for(int j = 1; j < (int)FC[0].size(); j++){
        if(!FC[0][j]->IsCalc){
            Eigen::Vector3d e = (FC[0][j]->first->p3 - FC[0][j-1]->first->p3);
            Eigen::Vector3d N = (e.cross(FC[0][j+1]->first->p3 - FC[0][j-1]->first->p3));
            Rulings.push_back(RevisionVertices::decideRulingDirectionOn3d(e, N, X[2*x_ind], X[2*x_ind + 1]));
            x_ind++;
        }
    }
    f = Normal_Sim(FC[0], Rulings);
    if(!grad.empty()){
        x_ind = 0;
        for(int j = 1; j < (int)FC[0].size(); j++){
            if(!FC[0][j]->IsCalc){
                Eigen::Vector3d e = FC[0][j]->first->p3 - FC[0][j-1]->first->p3;
                Eigen::Vector3d N = e.cross(FC[0][j+1]->first->p3 - FC[0][j-1]->first->p3);

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

double RevisionVertices::E_iso(std::vector<std::vector<Eigen::Vector3d>>& P, std::vector<std::vector<Eigen::Vector3d>>& Pori, std::vector<double>& grad){
    double e = 0.0;
    auto f = [](const std::vector<Eigen::Vector3d>& P, const std::vector<Eigen::Vector3d>& Pori)->double{
        double val;
         for(int i = 1; i < (int)P.size(); i++)val += std::abs((P[i] - P[i-1]).norm() - (Pori[i] - Pori[i-1]).norm());
        return val;
    };
    for(int i = 0; i < (int)P.size(); i++)e += f(P[i], Pori[i]);
    int i = 0;
    double fp, fm;
    for(int n = 0; n < (int)P.size(); n++){
        for(int j = 1; j < (int)P[n].size() - 1; j++){
            P[n][j].x() += eps; fp = f(P[n], Pori[n]); P[n][j].x() -= 2.0 * eps; fm = f(P[n], Pori[n]);  P[n][j].x() += eps; grad[i++] = (fp - fm)/(2.0*eps);
            P[n][j].y() += eps; fp = f(P[n], Pori[n]); P[n][j].y() -= 2.0 * eps; fm = f(P[n], Pori[n]);  P[n][j].y() += eps; grad[i++] = (fp - fm)/(2.0*eps);
            P[n][j].z() += eps; fp = f(P[n], Pori[n]); P[n][j].z() -= 2.0 * eps; fm = f(P[n], Pori[n]);  P[n][j].z() += eps; grad[i++] = (fp - fm)/(2.0*eps);
        }
    }
    return e;
}

double RevisionVertices::E_fair(std::vector<std::vector<Eigen::Vector3d>>& P, std::vector<double>& grad){

    auto f = [](std::vector<Eigen::Vector3d>& P)->double{
        double val;
        for(int i = 1; i < (int)P.size() - 1; i++)val += std::pow((P[i-1] - 2.0 * P[i] + P[i+1]).norm(), 2);
        return val;
    };
    double val = 0.0; for(auto& p: P)val += f(p);
    int i = 0;
    double fp, fm;
    for(auto&p : P){
        for(int j = 1; j < (int)p.size() - 1; j++){
            p[j].x() += eps; fp = f(p); p[j].x() -= 2.0 * eps; fm = f(p);  p[j].x() += eps; grad[i++] = (fp - fm)/(2.0*eps);
            p[j].y() += eps; fp = f(p); p[j].y() -= 2.0 * eps; fm = f(p);  p[j].y() += eps; grad[i++] = (fp - fm)/(2.0*eps);
            p[j].z() += eps; fp = f(p); p[j].z() -= 2.0 * eps; fm = f(p);  p[j].z() += eps; grad[i++] = (fp - fm)/(2.0*eps);
        }
    }
    return val;
}

double RevisionVertices::E_sim(std::vector<std::vector<Eigen::Vector3d>>& P, std::vector<std::vector<Eigen::Vector3d>>& Pori, std::vector<double>& grad){
    double e = 0.0;
    auto f = [](std::vector<Eigen::Vector3d>& P, std::vector<Eigen::Vector3d>& Pori)->double{
        double val;
         for(int i = 0; i < (int)P.size(); i++)val += (P[i] - Pori[i]).norm();
        return val;
    };
    for(int i = 0; i < (int)P.size(); i++)e += f(P[i], Pori[i]);
    int i = 0;
    double fp, fm;
    for(int n = 0; n < (int)P.size(); n++){
        for(int j = 1; j < (int)P[n].size() - 1; j++){
            P[n][j].x() += eps; fp = f(P[n], Pori[n]); P[n][j].x() -= 2.0 * eps; fm = f(P[n], Pori[n]);  P[n][j].x() += eps; grad[i++] = (fp - fm)/(2.0*eps);
            P[n][j].y() += eps; fp = f(P[n], Pori[n]); P[n][j].y() -= 2.0 * eps; fm = f(P[n], Pori[n]);  P[n][j].y() += eps; grad[i++] = (fp - fm)/(2.0*eps);
            P[n][j].z() += eps; fp = f(P[n], Pori[n]); P[n][j].z() -= 2.0 * eps; fm = f(P[n], Pori[n]);  P[n][j].z() += eps; grad[i++] = (fp - fm)/(2.0*eps);
        }
    }
    return e;
}

double RevisionVertices::E_normal(std::vector<Eigen::Vector3d>& R, const std::vector<std::shared_ptr<Vertex4d>>& FC){
    std::vector<int> Vertices_Indicator;
    for(int i = 0; i < (int)FC.size(); i++){if(FC[i]->IsCalc)Vertices_Indicator.push_back(i);}
    double f = .0;
    int j = 0;
    for(int k = 0; k < (int)Vertices_Indicator.size() - 1; k++){
        for(int i = Vertices_Indicator[k] + 1; i < Vertices_Indicator[k+1]; i++){
            Eigen::Vector3d N = ((FC[i-1]->first->p3 - FC[i]->first->p3).cross(R[j])).normalized();
            Eigen::Vector3d N_ori = ((FC[i-1]->first->p3 - FC[i]->first->p3).cross(FC[i]->second->p3 - FC[i]->first->p3)).normalized();
            f += 1.0 - N.dot(N_ori);
            N = (R[j].cross(FC[i+1]->first->p3 - FC[i]->first->p3)).normalized();
            N_ori = ((FC[i]->second->p3 - FC[i]->first->p3).cross(FC[i+1]->first->p3 - FC[i]->first->p3)).normalized();
            f += 1.0 - N.dot(N_ori);
            j++;
        }
    }

    return f;
}

double RevisionVertices::Minimize_SmoothSrf(const std::vector<double>& X, std::vector<double> &grad, void *f_data){
    SmoothSurface *od = (SmoothSurface*)f_data;
    std::vector<SmoothingArea> SA = od->SA;
    std::vector<std::vector<Eigen::Vector3d>> P,Pori, P2d;
    int i = 0;
    for(auto&sa: SA){
        std::vector<Eigen::Vector3d> tmp = {sa.stP->first->p3}, tmp_ori = {sa.stP->first->p3_ori};
        Eigen::Vector3d vt = (sa.stP->second->p3 - sa.stP->first->p3).normalized(), vt2d = (sa.stP->second->p - sa.stP->first->p).normalized();
        Eigen::Vector3d SpinAxis = (vt2d.y() < 0)? -Eigen::Vector3d(0,0,-1): Eigen::Vector3d(0,0,1);
        std::vector<Eigen::Vector3d> _P2d;
        Eigen::Vector3d befP3 = sa.stP->first->p3, befP2 = sa.stP->first->p;

        for(auto&p: sa.OriginalVertices){
            Eigen::Vector3d p3d(X[i], X[i+1], X[i+2]);
            tmp_ori.push_back(p->p3_ori);
            tmp.push_back(p3d);
            i += 3;
            Eigen::Vector3d e = p3d - befP3;
            double l = e.norm(), phi = std::acos(vt.dot(e.normalized()));
            Eigen::Vector3d p2d = l * (Eigen::AngleAxisd(phi, SpinAxis) * vt2d) + befP2;
            _P2d.push_back(p2d);
            befP2 = p2d; befP3 = p3d;
            if(sa.qt != nullptr){
                vt = (sa.qt->p3 - p3d).normalized(); vt2d = (sa.qt->p - p2d).normalized();
                SpinAxis = (vt2d.y() < 0)? Eigen::Vector3d(0,0,-1): Eigen::Vector3d(0,0,1);
            }
        }tmp.push_back(sa.lastP->first->p3); tmp_ori.push_back(sa.lastP->first->p3_ori);
        P.push_back(tmp); Pori.push_back(tmp_ori); P2d.push_back(_P2d);
    }

    double f_sim, f_fair, f_iso;
    double w_sim = 10, w_fair = 0.1, w_iso = 0.;
    if(!grad.empty()){
        std::vector<double> dE_sim(grad.size()), dE_fair(grad.size()), dE_iso(grad.size());
        f_sim = RevisionVertices::E_sim(P, Pori, dE_sim); f_fair = RevisionVertices::E_fair(P2d, dE_fair), f_iso = RevisionVertices::E_iso(P, Pori, dE_iso);
        for(int i = 0; i < (int)grad.size(); i++)grad[i] = w_sim * dE_sim[i] + w_fair * dE_fair[i] + w_iso * dE_iso[i];

    }
    if(DebugMode::Singleton::getInstance().isdebug())std::cout<< "similarilty = " << f_sim << "  ,  fairness = " << f_fair << ", isometric = " << f_iso << std::endl;
    return w_sim * f_sim + w_fair * f_fair + w_iso * f_iso;
}

double RevisionVertices::Minimize_PlanaritySrf(const std::vector<double>& X, std::vector<double> &grad, void *f_data){

    std::vector<std::shared_ptr<Vertex4d>> *FC = (std::vector<std::shared_ptr<Vertex4d>>*)f_data;
    std::vector<Eigen::Vector3d> Rulings;

    int x_ind = 0;
    for(int j = 1; j < (int)FC[0].size(); j++){

        if(!FC[0][j]->IsCalc){
            Eigen::Vector3d e = (FC[0][j]->first->p3 - FC[0][j-1]->first->p3);
            Eigen::Vector3d N = e.cross(FC[0][j+1]->first->p3 - FC[0][j-1]->first->p3);
            Rulings.push_back(RevisionVertices::decideRulingDirectionOn3d(e, N, X[2*x_ind], X[2*x_ind + 1]));
            x_ind++;
        }
    }

    double f = 0.0, fp, fm;
    f = RevisionVertices::E_normal(Rulings, FC[0]);
    if(!grad.empty()){
        x_ind = 0;
        for(int j = 1; j < (int)FC[0].size(); j++){
            if(!FC[0][j]->IsCalc){
                Eigen::Vector3d e = FC[0][j]->first->p3 - FC[0][j-1]->first->p3;
                Eigen::Vector3d N = e.cross(FC[0][j+1]->first->p3 - FC[0][j-1]->first->p3);

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

bool FoldLine::SimpleSmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v){

    auto getV = [](const std::shared_ptr<Vertex>& o, const std::shared_ptr<Vertex>& x, const std::shared_ptr<Vertex>& o2, const std::shared_ptr<Vertex>& x2)->std::shared_ptr<Vertex>{
        if(IsParallel(o, x, o2, x2))return std::shared_ptr<Vertex>(nullptr);
        Eigen::Vector3d p2d = calcCrossPoint_2Vertex(o, x, o2, x2);
        Eigen::Vector3d p3d = calcTargetDistanceOnPlane(p2d, o,  x, x2);
        return std::make_shared<Vertex>(p2d, p3d);
    };

    Eigen::Vector3d r3d, r2d, _r3d, _r2d;
    std::vector<int>Vertices_Ind;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i]->IsCalc)Vertices_Ind.push_back(i);}

    std::shared_ptr<Vertex> qt;
    for(int j = 0; j < (int)Vertices_Ind.size() - 1; j++){
        qt = nullptr;
        if(j != 0 && j != (int)Vertices_Ind.size() - 2){
            qt = getV(FoldingCurve[Vertices_Ind[j]]->first, FoldingCurve[Vertices_Ind[j]]->second, FoldingCurve[Vertices_Ind[j+1]]->first, FoldingCurve[Vertices_Ind[j+1]]->second);
             _r3d = (FoldingCurve[Vertices_Ind[j]]->second->p3 - FoldingCurve[Vertices_Ind[j]]->first->p3).normalized(), _r2d = (FoldingCurve[Vertices_Ind[j]]->second->p - FoldingCurve[Vertices_Ind[j]]->first->p).normalized();
        }
        for(int i = Vertices_Ind[j] + 1; i < Vertices_Ind[j+1]; i++){
            FoldingCurve[i]->first->p3 = FoldingCurve[i]->first->p3_ori; FoldingCurve[i]->first->p = FoldingCurve[i]->first->p2_ori;
            if(j == 0){
                r2d = (FoldingCurve[Vertices_Ind[1]]->second->p - FoldingCurve[Vertices_Ind[1]]->first->p).normalized();
                r3d = (FoldingCurve[Vertices_Ind[1]]->second->p3 - FoldingCurve[Vertices_Ind[1]]->first->p3).normalized();
            }else if(j == (int)Vertices_Ind.size() - 2){
                r2d = (FoldingCurve[Vertices_Ind[j]]->second->p - FoldingCurve[Vertices_Ind[j]]->first->p).normalized();
                r3d = (FoldingCurve[Vertices_Ind[j]]->second->p3 - FoldingCurve[Vertices_Ind[j]]->first->p3).normalized();
            }else{
                if(qt != nullptr){
                    auto e = (FoldingCurve[i-1]->first->p - FoldingCurve[i]->first->p).normalized(), e2 = (FoldingCurve[i+1]->first->p - FoldingCurve[i]->first->p).normalized();
                    r2d = (qt->p - FoldingCurve[i]->first->p).normalized(); r3d = (qt->p3 - FoldingCurve[i]->first->p3).normalized();
                    double k = e.dot(e2);
                    k = (k < -1)? std::numbers::pi: (k > 1)? 0: std::acos(k);
                    if((e.cross(e2)).z() > 0)k = 2.0*std::numbers::pi - k;
                    double phi = e.dot(r2d); phi = (phi < -1)? std::numbers::pi: (phi > 1)? 0: std::acos(phi);
                    if((e.cross(r2d)).z() > 0)phi = 2.0*std::numbers::pi - phi;
                    if(phi > k)r2d *= -1;
                    if(r3d.dot(_r3d) < 0)r3d *= -1;
                }else{
                    r2d = (FoldingCurve[Vertices_Ind[j]]->second->p - FoldingCurve[Vertices_Ind[j]]->first->p).normalized();
                    r3d = (FoldingCurve[Vertices_Ind[j]]->second->p3 - FoldingCurve[Vertices_Ind[j]]->first->p3).normalized();
                }
            }
            
            Eigen::Vector3d v2, p_clst;
            for(int k = 0; k < (int)Poly_v.size(); k++){
                v2 = MathTool::calcCrossPoint_2Vector(FoldingCurve[i]->first->p, 1000.0 * r2d + FoldingCurve[i]->first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
                if(MathTool::is_point_on_line(v2, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p) && r2d.dot(v2 - FoldingCurve[i]->first->p) > 0) p_clst = v2;
            }
            FoldingCurve[i]->IsCalc = true;
            FoldingCurve[i]->second->p = p_clst;
            FoldingCurve[i]->second->p3 = (p_clst - FoldingCurve[i]->first->p).norm() * r3d + FoldingCurve[i]->first->p3;
            FoldingCurve[i]->third->p = FoldingCurve[i]->third->p2_ori;
            FoldingCurve[i]->third->p3 = FoldingCurve[i]->third->p3_ori;
        }
    }
    std::shared_ptr<Vertex> v_clst = getClosestVertex(FoldingCurve[0]->second, FoldingCurve[0]->first, FoldingCurve, false);
    if(v_clst != nullptr){
        int j;
        for(j = 0; j < (int)FoldingCurve.size(); j++){if(v_clst == FoldingCurve[j]->second)break;}
        FoldingCurve[0]->second->p3 = calcTargetDistanceOnPlane(FoldingCurve[0]->second->p, FoldingCurve[j]->first, FoldingCurve[j]->second, FoldingCurve[j+1]->second);
    }else FoldingCurve[0]->second->p3 = calcTargetDistanceOnPlane(FoldingCurve[0]->second->p, FoldingCurve[1]->first, FoldingCurve[0]->first, FoldingCurve[1]->second);

    v_clst = getClosestVertex(FoldingCurve.back()->second, FoldingCurve.back()->first, FoldingCurve, false);
    if(v_clst != nullptr){
        int j;
        for(j = 0; j < (int)FoldingCurve.size(); j++){if(v_clst == FoldingCurve[j]->second)break;}
        FoldingCurve.back()->second->p3 = calcTargetDistanceOnPlane(FoldingCurve.back()->second->p, FoldingCurve[j]->first, FoldingCurve[j]->second, FoldingCurve[j-1]->second);
    }else FoldingCurve.back()->second->p3 = calcTargetDistanceOnPlane(FoldingCurve.back()->second->p, FoldingCurve.end()[-2]->first, FoldingCurve.back()->first, FoldingCurve.end()[-2]->second);
    return true;
}

std::vector<std::vector<Eigen::Vector3d>> FoldLine::Optimization_SmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v, bool IsConnectEndPoint){
    std::vector<double> X;
    std::vector<std::vector<Eigen::Vector3d>> res_qt;
    std::shared_ptr<Vertex4d> bef = FoldingCurve.front();
    std::vector<RevisionVertices::SmoothingArea> SA;
    std::vector<std::shared_ptr<Vertex>> OriginalVertices;
    int bef_ind;
    for(auto&FC : FoldingCurve){
        if(!FC->IsCalc){
            X.push_back(FC->first->p3_ori.x()); X.push_back(FC->first->p3_ori.y()); X.push_back(FC->first->p3_ori.z());
            OriginalVertices.push_back(FC->first);
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
        if(DebugMode::Singleton::getInstance().isdebug()){
            std::cout <<"result :  smoothing  " <<result << std::endl;
            std::cout << "found minimum at f = "  << std::setprecision(10) << minf << std::endl;
        }
        int i = 0;
        for(auto&FC: FoldingCurve){
            if(!FC->IsCalc){
                FC->first->p3 = Eigen::Vector3d(X[i], X[i+1], X[i+2]); i += 3;
            }
        }
        Eigen::Vector3d p, r3d, r2d;
        for(auto sm: od.SA){
            Eigen::Vector3d vt = (sm.stP->second->p3 - sm.stP->first->p3).normalized(), vt2d = (sm.stP->second->p - sm.stP->first->p).normalized();
            Eigen::Vector3d SpinAxis = (vt2d.y() < 0)? Eigen::Vector3d(0,0,-1): Eigen::Vector3d(0,0,1);
            for(int j = sm.st_ind + 1; j < sm.last_ind; j++){
                Eigen::Vector3d e = FoldingCurve[j]->first->p3 - FoldingCurve[j-1]->first->p3;
                double l = e.norm(), phi = std::acos(vt.dot(e.normalized()));
                FoldingCurve[j]->first->p = l * (Eigen::AngleAxisd(phi, SpinAxis) * vt2d) + FoldingCurve[j-1]->first->p;

                if(sm.qt == nullptr){
                    if(sm.st_ind != 0){r2d = (sm.stP->second->p - sm.stP->first->p).normalized(); r3d = (sm.stP->second->p3 - sm.stP->first->p3).normalized();}
                    else {r2d = (sm.lastP->second->p - sm.lastP->first->p).normalized(); r3d = (sm.lastP->second->p3 - sm.lastP->first->p3).normalized();}
                }else{
                    res_qt.push_back({sm.stP->first->p3, sm.qt->p3, sm.lastP->first->p3});

                    r2d = (sm.qt->p - FoldingCurve[j]->first->p).normalized();
                    if(r2d.dot(FoldingCurve[sm.st_ind]->second->p - FoldingCurve[sm.st_ind]->first->p) < 0)r2d *= -1;
                    FoldingCurve[j]->second->p = MathTool::calcCrossPoint_2Vector(FoldingCurve[j]->first->p, 1000.0 * r2d + FoldingCurve[j]->first->p, sm.stP->second->p, sm.lastP->second->p);
                    r3d = (sm.qt->p3 - FoldingCurve[j]->first->p3).normalized();
                    if(r3d.dot((sm.stP->second->p3 - sm.stP->first->p3).normalized()) < 0)r3d *= -1;
                    vt = (sm.qt->p3 - FoldingCurve[j]->first->p3).normalized(); vt2d = (sm.qt->p - FoldingCurve[j]->first->p).normalized();
                    SpinAxis = (vt2d.y() < 0)? Eigen::Vector3d(0,0,-1): Eigen::Vector3d(0,0,1);
                }
                for(int k = 0; k < (int)Poly_v.size(); k++){
                    p = MathTool::calcCrossPoint_2Vector(FoldingCurve[j]->first->p, 1000.0 * r2d + FoldingCurve[j]->first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
                    if(MathTool::is_point_on_line(p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p) && r2d.dot(p - FoldingCurve[j]->first->p) > 0)break;
                }
                FoldingCurve[j]->second->p = p;
                FoldingCurve[j]->second->p3 = (FoldingCurve[j]->second->p - FoldingCurve[j]->first->p).norm() * r3d + FoldingCurve[j]->first->p3;

                if(sm.qb == nullptr){
                    if(sm.st_ind != 0){r2d = (sm.stP->third->p - sm.stP->first->p).normalized(); r3d = (sm.stP->third->p3 - sm.stP->first->p3).normalized();}
                    else {r2d = (sm.lastP->third->p - sm.lastP->first->p).normalized(); r3d = (sm.lastP->third->p3 - sm.lastP->first->p3).normalized();}
                }else{
                    r2d = (sm.qb->p - FoldingCurve[j]->first->p).normalized();r3d = (sm.qb->p3 - FoldingCurve[j]->first->p3).normalized();
                    if(r2d.dot(FoldingCurve[sm.st_ind]->third->p - FoldingCurve[sm.st_ind]->first->p) < 0)r2d *= -1;
                }
                for(int k = 0; k < (int)Poly_v.size(); k++){
                    p = MathTool::calcCrossPoint_2Vector(FoldingCurve[j]->first->p, 1000.0 * r2d + FoldingCurve[j]->first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
                    if(MathTool::is_point_on_line(p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p) && r2d.dot(p - FoldingCurve[j]->first->p) > 0)break;
                }
                FoldingCurve[j]->third->p = p;
                FoldingCurve[j]->third->p3 = (p - FoldingCurve[j]->first->p).norm() * r3d + FoldingCurve[j]->first->p3;

                auto v4d = FoldingCurve[j];
                Eigen::Vector3d et = (v4d->second->p3 -v4d->first->p3).normalized(), er = (FoldingCurve[j-1]->first->p3 -v4d->first->p3).normalized(), eb = (v4d->third->p3 -v4d->first->p3).normalized(), el = (FoldingCurve[j+1]->first->p3 -v4d->first->p3).normalized();
                double phi1 = std::acos(et.dot(er)), phi2 = std::acos(et.dot(el)), phi3 = std::acos(eb.dot(el)), phi4 = std::acos(eb.dot(er));
                if(DebugMode::Singleton::getInstance().isdebug())
                    std::cout << j << " : " << abs(2.0*std::numbers::pi - phi1 - phi2 - phi3 - phi4)  << " , " << MathTool::rad2deg(phi1) << " , " << MathTool::rad2deg(phi2) << " , " << MathTool::rad2deg(phi3) << ", " << MathTool::rad2deg(phi4) << std::endl;
            }
        }
        std::shared_ptr<Vertex> v_clst = getClosestVertex(FoldingCurve[0]->second, FoldingCurve[0]->first, FoldingCurve, false);
        if(v_clst != nullptr){
            int j; for(j = 0; j < (int)FoldingCurve.size(); j++){if(v_clst == FoldingCurve[j]->second)break;}
            FoldingCurve[0]->second->p3 = calcTargetDistanceOnPlane(FoldingCurve[0]->second->p, FoldingCurve[j]->first, FoldingCurve[j]->second, FoldingCurve[j+1]->second);
        } else FoldingCurve[0]->second->p3 = calcTargetDistanceOnPlane(FoldingCurve[0]->second->p, FoldingCurve[1]->second, FoldingCurve[0]->first, FoldingCurve[1]->first);
        v_clst = getClosestVertex(FoldingCurve.back()->second, FoldingCurve.back()->first, FoldingCurve, false);
        if(v_clst != nullptr){
            int j; for(j = 0; j < (int)FoldingCurve.size(); j++){if(v_clst == FoldingCurve[j]->second)break;}
            FoldingCurve.back()->second->p3 = calcTargetDistanceOnPlane(FoldingCurve.back()->second->p, FoldingCurve[j]->first, FoldingCurve[j]->second, FoldingCurve[j-1]->second);
        }else FoldingCurve.back()->second->p3 = calcTargetDistanceOnPlane(FoldingCurve.back()->second->p, FoldingCurve.end()[-2]->second, FoldingCurve.end()[-2]->first, FoldingCurve.back()->first);
        std::cout<<"smoothing finish"<<std::endl;
    }catch (std::exception& e) {std::cout << "nlopt failed: " << e.what() << std::endl; }

    return res_qt;
}

std::vector<std::vector<Eigen::Vector3d>> FoldLine::Optimization_PlanaritySrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v){
    std::vector<double> X;
    std::vector<std::vector<Eigen::Vector3d>> res_qt;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){
        if(!FoldingCurve[i]->IsCalc){
           Eigen::Vector3d e = (FoldingCurve[i-1]->first->p3 - FoldingCurve[i]->first->p3).normalized();
           Eigen::Vector3d N = (e.cross(FoldingCurve[i+1]->first->p3 - FoldingCurve[i]->first->p3)).normalized();
           Eigen::Vector3d r = (FoldingCurve[i]->second->p3 - FoldingCurve[i]->first->p3).normalized();
           double phi = (e.dot(r) <= -1)? std::numbers::pi: (e.dot(r) >= 1)? 0: std::acos(e.dot(r));
           r = MathTool::ProjectionVector(r, e, true); N = MathTool::ProjectionVector(N, e, true);
           double a = (N.dot(r) <= -1)? std::numbers::pi: (N.dot(r) >= 1)? 0: std::acos(N.dot(r));
           if(e.dot(r.cross(N)) < 0) a = std::numbers::pi/2.0 + a;  
           else a = (std::numbers::pi/2.0 - a > 0)? 2.0*std::numbers::pi - a: std::numbers::pi/2.0 - a;
           
           //a = 2.0 * std::numbers::pi - a;
           if(DebugMode::Singleton::getInstance().isdebug()) std::cout << MathTool::rad2deg(a)<<std::endl;
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
        if(DebugMode::Singleton::getInstance().isdebug()){
            std::cout <<"result :  planarity  " <<result << std::endl;
            std::cout << "found minimum at f = "  << std::setprecision(10) << minf << std::endl;
        }
        int x_ind = 0;
        Eigen::Vector3d p;
        for(auto itr = FoldingCurve.begin(); itr != FoldingCurve.end(); itr++){
            if(!(*itr)->IsCalc){
                auto _prev = itr-1, _next = itr + 1;
                if(DebugMode::Singleton::getInstance().isdebug())
                    std::cout << std::distance(FoldingCurve.begin(), itr) << ": "  << MathTool::rad2deg(X[2 * x_ind]) << " , " << MathTool::rad2deg(X[2 * x_ind + 1]) << std::endl;
                Eigen::Vector3d e = (*_prev)->first->p3 - (*itr)->first->p3;
                Eigen::Vector3d N = e.cross((*_next)->first->p3 - (*itr)->first->p3);
                Eigen::Vector3d r3d = RevisionVertices::decideRulingDirectionOn3d(e, N, X[2 * x_ind], X[2 * x_ind + 1]);
                Eigen::Vector3d r2d = (Eigen::AngleAxisd(X[2 * x_ind + 1], Eigen::Vector3d(0,0,-1)) * ((*_prev)->first->p - (*itr)->first->p)).normalized();
                for(int k = 0; k < (int)Poly_v.size(); k++){
                    p = MathTool::calcCrossPoint_2Vector((*itr)->first->p, 1000.0 * r2d + (*itr)->first->p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p);
                    if(MathTool::is_point_on_line(p, Poly_v[k]->p, Poly_v[(k + 1) % (int)Poly_v.size()]->p) && r2d.dot(p - (*itr)->first->p) > 0)break;
                }
                (*itr)->second->p = p;
                (*itr)->second->p3 = (p - (*itr)->first->p).norm() * r3d + (*itr)->first->p3;
                x_ind++;
            }
        }
    }catch (std::exception& e) {std::cout << "nlopt failed: " << e.what() << std::endl; }
    std::cout <<"planarity optimization finish"<<std::endl;
    return res_qt;
}

bool FoldLine::Optimization_FlapAngle(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double wb, double wp, bool ConstFunc){
    if(FoldingCurve.empty())return false;

    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    double a = 0, a2 = 0.0;
    Eigen::Vector3d e = (ValidFC[0]->first->p3 - ValidFC[1]->first->p3).normalized(),e2 = (ValidFC[2]->first->p3 - ValidFC[1]->first->p3).normalized();
    Eigen::Vector3d SpinAxis = (ValidFC[1]->third->p3 - ValidFC[1]->first->p3).normalized();
    Eigen::Vector3d Nb = -(e.cross(e2)).normalized(), N4 = (e.cross(SpinAxis)).normalized();

    double a_ll = std::atan2(e.dot(Nb.cross(N4)), N4.dot(-Nb));
    if(a_ll < 0)a_ll += 2.0*std::numbers::pi;
    double a_con = a_ll + std::numbers::pi;
    if(a_con > 2.0*std::numbers::pi)a_con -= 2.0*std::numbers::pi;
    if(a_con < 0)a_con +=2.0*std::numbers::pi;

    double phi3 = std::acos(e2.dot(SpinAxis)), phi4 = std::acos(e.dot(SpinAxis));
    double k = 2.0 * std::numbers::pi - phi3 - phi4;
    double a_min, a_max;

    Eigen::Vector3d FaceNp = (SpinAxis.cross(e2)).normalized(), CrossV = (N4.cross(FaceNp)).normalized();
    bool IsMount = (SpinAxis.dot(CrossV) < 0)? true: false;

    if(k < std::numbers::pi && IsMount){a_min = a_con + std::numbers::pi; a_max = 2.0 * std::numbers::pi;}
    if(k >= std::numbers::pi && IsMount){a_min = std::numbers::pi; a_max = std::numbers::pi + a_con;}
    if(k < std::numbers::pi && !IsMount){a_min = 0.0;a_max = a_con - std::numbers::pi;}
    if(k >= std::numbers::pi && !IsMount){a_min = a_con - std::numbers::pi; a_max = std::numbers::pi;}
    a = (a_min + a_max)/2.0;

    if(DebugMode::Singleton::getInstance().isdebug()) std::cout <<MathTool::rad2deg(a_ll) << " , " << MathTool::rad2deg(a_con)<<std::endl;
    if(DebugMode::Singleton::getInstance().isdebug()){
        std::ofstream ofs2;
        std::filesystem::create_directory("./Optimization");

        std::string AngleFile = "./Optimization/ChangeAngle.csv" ;
        ofs2.open(AngleFile, std::ios::out); ofs2 << "a(radian),a(degree),Eruling, Ebend, Eparalell, Eruling2"<<std::endl;
        while(a <= 2.0*std::numbers::pi){
             _FoldingAAAMethod(FoldingCurve, Poly_V, a);
             double f = RulingsCrossed(FoldingCurve);
             double fb = Fbend2(FoldingCurve);
             double fp =  Fparallel(FoldingCurve);
             double f2 =RulingsCrossed_NoOutline(FoldingCurve);
             ofs2 << a << "," << MathTool::rad2deg(a) << ", " << f << ", " << fb << ", " << fp <<", " << f2 <<  std::endl;
            a += 1e-3;
        }ofs2.close();

        AngleFile = "./Optimization/OptimArea.csv";
        ofs2.open(AngleFile, std::ios::out);
        ofs2 << "a(radian), a(degree) , Eruling, Ebend, Eparalell, Ebend, Eruling2"<<std::endl;
        for(double _a = a_min; _a <= a_max; _a+= 1e-3){
             _FoldingAAAMethod(FoldingCurve, Poly_V, _a);
             double f = RulingsCrossed(FoldingCurve);
             double fb = Fbend2(FoldingCurve);
              double fp =  Fparallel(FoldingCurve);
             double fr =RulingsCrossed_NoOutline(FoldingCurve);
            ofs2 << _a << ", " << MathTool::rad2deg(_a) << " , " << f << ", " << fb << ", " <<fp <<", "<< fr << std::endl;
        }ofs2.close();
    }

    RevisionVertices::ObjData od = {FoldingCurve, Poly_V};
    ConstFunc = false;
    if(!ConstFunc)od.AddWeight(wb, wp, -1); else od.AddWeight(wb, wp, 100.0);
    //od.AddWeight(wb, wp, 100.0);
    nlopt::opt opt;
    opt = nlopt::opt(nlopt::LD_MMA, 1);
    opt.set_min_objective(ObjFunc_FlapAngle, &od);
    opt.set_lower_bounds(a_min);
    opt.set_upper_bounds(a_max);
    opt.add_inequality_constraint(Fruling, &od);
    //opt.set_param("inner_maxeval", 100);
    opt.set_maxtime(2.0);//stop over this time
    opt.set_xtol_rel(1e-13);
    if(DebugMode::Singleton::getInstance().isdebug())
        std::cout << "area " << MathTool::rad2deg(a_min) << " < " << MathTool::rad2deg(a) << " < " << MathTool::rad2deg(a_max) << std::endl;

    double minf_amin, minf_amax, f_ruling_amin, f_ruling_amax;
    double res_amin, res_amax;
    ofs_Ebend.open(File_Ebend, std::ios::out); ofs_Eruling.open(File_Eruling, std::ios::out);
      try {
        std::vector<double> _a{a_min + 0.5};
          nlopt::result result = opt.optimize(_a, minf_amin);
          _FoldingAAAMethod(FoldingCurve, Poly_V, _a[0]);
          std::cout <<"result :  lower bound" <<result << std::endl;
          //f_ruling_amin = (ConstFunc)? RulingsCrossed(FoldingCurve):RulingsCrossed_NoOutline(FoldingCurve);

          f_ruling_amin = RulingsCrossed(FoldingCurve);
          res_amin = _a[0];
          std::cout << "found minimum at f(" << MathTool::rad2deg(_a[0]) << ") = " << std::setprecision(10) << minf_amin << "  ,  " << f_ruling_amin <<  std::endl;          
      }
      catch (std::exception& e) {
          std::cout << "nlopt failed: " << e.what() << std::endl;
      }ofs_Ebend.close(); ofs_Eruling.close();

    try {
      std::vector<double> _a{a_max - 0.5};
        nlopt::result result = opt.optimize(_a, minf_amax);
        //f_ruling_amax = (ConstFunc)? RulingsCrossed(FoldingCurve):RulingsCrossed_NoOutline(FoldingCurve);
         _FoldingAAAMethod(FoldingCurve, Poly_V, _a[0]);
        f_ruling_amax = RulingsCrossed(FoldingCurve);
        std::cout <<"result :  upper bound" <<result << std::endl;
        std::cout << "found minimum at f(" << MathTool::rad2deg(_a[0]) << ") = " << std::setprecision(10) << minf_amax << "  ,  "  << f_ruling_amax << std::endl;
        res_amax = _a[0];

    }
    catch (std::exception& e){std::cout << "nlopt failed: " << e.what() << std::endl;}
    
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
   a_flap = a2;
    std::cout << "result : smaller = " << MathTool::rad2deg(a2)  << "(" << a2 << "), tol = " << tol <<  "  f_bend = " << Fbend2(FoldingCurve) << "  ,  f_ruling = " <<RulingsCrossed(FoldingCurve) << std::endl;
    std::cout << "finish"<<std::endl;
    return true;
}

inline double RevisionVertices::getK(const Eigen::Vector3d o, const Eigen::Vector3d x, const Eigen::Vector3d x2){
    auto v = (x - o).normalized(), v2 = (x2 - o).normalized();
    double k = std::acos(v.dot(v2));
    if(Eigen::Vector3d::UnitZ().dot(v.cross(v2)) > 0)k = 2.0*std::numbers::pi - k;
    return k;
}

void FoldLine::SimplifyModel(int iselim, bool isroot){
    auto elim_rulings = [&]()->std::vector<std::shared_ptr<Vertex4d>>{
        struct Point{
            std::shared_ptr<Vertex4d> p;
            int index;
            double dist;
            Point(){}
            Point(std::shared_ptr<Vertex4d> _p, int i, double d): p(_p), index(i), dist(d) {}
        };

        auto perpendicularDistance = [](Eigen::Vector3d& p, Eigen::Vector3d& l_start, Eigen::Vector3d& l_end)->double{
            Eigen::Vector3d line = (l_end - l_start).normalized(), OP = (p - l_start).normalized();
            return (l_start + line.dot(OP) * (p - l_start) - p).norm();
        };
        auto DPA = [&](std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex4d>>& res){
            double dmax = 0.0;
            int index = -1, end = (int)FoldingCurve.size()-1;
            if((int)FoldingCurve.size() < 2)return Point(FoldingCurve.front(), index, dmax);
            for(int i = 1; i < end; i++){
                if(std::find(res.begin(), res.end(), FoldingCurve[i]) != res.end())continue;
                double d = perpendicularDistance(FoldingCurve[i]->first->p, FoldingCurve[0]->first->p, FoldingCurve.back()->first->p);
                if (d > dmax){ index = i; dmax = d;}
            }
            return (index != -1) ? Point(FoldingCurve[index], index, dmax): Point(FoldingCurve.front(), index, dmax);
        };
        int index = 0;
        std::vector<std::shared_ptr<Vertex4d>> res{FoldingCurve.front(), FoldingCurve.back()};
        do{
            std::vector<std::shared_ptr<Vertex4d>> firstLine{FoldingCurve.begin(), FoldingCurve.begin()+index+1};
            std::vector<std::shared_ptr<Vertex4d>> lastLine{FoldingCurve.begin() + index, FoldingCurve.end()};
            Point p_first = DPA(firstLine, res), p_second = DPA(lastLine,res);
            Point p;//最も遠い点
            if(p_first.index == -1 && p_second.index == -1)return res;
            if(p_first.index == -1 && p_second.index != -1)p = p_second;
            else if(p_first.index != -1 && p_second.index == -1)p = p_first;
            else p = (p_first.dist > p_second.dist) ? p_first : p_second;

            for(int i = 1; i < (int)res.size(); i++){
                if(res[i]->first->s < p.p->first->s && p.p->first->s < res[i-1]->first->s){res.insert(res.begin() + i, p.p);break;}
            }
            index = p.index;
        }while((int)res.size() < validsize);
        return res;
    };

    int n = (isroot)? 0: 1;
    validsize -= iselim;
    for(auto it = FoldingCurve.begin() + n; it != FoldingCurve.end() - n; it++){
        (*it)->IsCalc = true;
        (*it)->first->p = (*it)->first->p2_ori; (*it)->first->p3 = (*it)->first->p3_ori;
        (*it)->second->p = (*it)->second->p2_ori; (*it)->second->p3 = (*it)->second->p3_ori;
        (*it)->third->p = (*it)->third->p2_ori; (*it)->third->p3 = (*it)->third->p3_ori;
    }
    auto res = elim_rulings();
    for(auto&V4d: FoldingCurve){
        if(std::find_if(res.begin(), res.end(), [&V4d](const std::shared_ptr<Vertex4d>& V){return V4d->first == V->first;}) != res.end()){
        }else V4d->IsCalc = false;
    }
    return;
}

void FoldLine::SimplifyModel(double tol, bool isroot){
    int n = (isroot)? 0: 1;
    for(auto it = FoldingCurve.begin() + n; it != FoldingCurve.end() - n; it++){
        (*it)->IsCalc = true;
        (*it)->first->p = (*it)->first->p2_ori; (*it)->first->p3 = (*it)->first->p3_ori;
        (*it)->second->p = (*it)->second->p2_ori; (*it)->second->p3 = (*it)->second->p3_ori;
        (*it)->third->p = (*it)->third->p2_ori; (*it)->third->p3 = (*it)->third->p3_ori;
    }
    auto tmp = TrimPoints(FoldingCurve, tol);
    for(auto&V4d: FoldingCurve){
        if(std::find_if(tmp.begin(), tmp.end(), [&V4d](const std::shared_ptr<Vertex4d>& V){return V4d->first == V->first;}) != tmp.end()){
        }else V4d->IsCalc = false;
    }
    this->tol = tol;
    return;
    std::vector<int>Vertices_Ind;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i]->IsCalc)Vertices_Ind.push_back(i);}
    for(int i = 0; i < (int)Vertices_Ind.size() - 1; i++){
    int ind_branch = Vertices_Ind[i + 1], ind_root = Vertices_Ind[i];
    if(i == 0){
        if(getClosestVertex(FoldingCurve[Vertices_Ind[0]]->second , FoldingCurve[Vertices_Ind[0]]->first, FoldingCurve) != nullptr){
            ind_root = Vertices_Ind[1]; ind_branch = Vertices_Ind[2];
        }
    }
    if(i == (int)Vertices_Ind.size() -2 ){
        if(getClosestVertex(FoldingCurve[Vertices_Ind.back()]->second , FoldingCurve[Vertices_Ind.back()]->first, FoldingCurve) != nullptr)
        {
            ind_root = Vertices_Ind[i]; ind_branch = Vertices_Ind[i-1];
        }

    }
        for(int j = Vertices_Ind[i] + 1; j < Vertices_Ind[i+1]; j++){
            FoldingCurve[j]->first->p = calcCrossPoint_2Vertex(FoldingCurve[j]->second, FoldingCurve[j]->third, FoldingCurve[Vertices_Ind[i]]->first, FoldingCurve[Vertices_Ind[i+1]]->first);
            //FoldingCurve[j]->first->p3 = calcTargetDistanceOnPlane(FoldingCurve[j]->first->p, FoldingCurve[ind_root]->first, FoldingCurve[ind_branch]->third, FoldingCurve[ind_branch]->second);
            double s = (FoldingCurve[j]->first->p - FoldingCurve[Vertices_Ind[i]]->first->p).norm()/(FoldingCurve[Vertices_Ind[i+1]]->first->p - FoldingCurve[Vertices_Ind[i]]->first->p).norm();
            FoldingCurve[j]->first->p3 = s * (FoldingCurve[Vertices_Ind[i+1]]->first->p3 - FoldingCurve[Vertices_Ind[i]]->first->p3)  + FoldingCurve[Vertices_Ind[i]]->first->p3;
            FoldingCurve[j]->third->p3 = calcTargetDistanceOnPlane(FoldingCurve[j]->third->p, FoldingCurve[ind_root]->first, FoldingCurve[ind_branch]->first, FoldingCurve[ind_branch]->second);
            FoldingCurve[j]->second->p3 = calcTargetDistanceOnPlane(FoldingCurve[j]->second->p, FoldingCurve[ind_root]->first, FoldingCurve[ind_branch]->first, FoldingCurve[ind_branch]->second);

        }
    }
}

bool FoldLine::RevisionCrosPtsPosition(){
    if(FoldingCurve.empty())return false;
    FoldingCurve.front()->first->p3 = FoldingCurve.front()->first->p3_ori;
    FoldingCurve.back()->first->p3 = FoldingCurve.back()->first->p3_ori;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve)if(fc->IsCalc)ValidFC.push_back(fc);
    if(FoldingCurve.size() > 5){
        std::cout << "before revision" << std::endl;
        for(int i = 1; i < (int)FoldingCurve.size()-1; i++){
            auto v = (FoldingCurve[i-1]->first->p - FoldingCurve[i]->first->p).normalized(), v2 = (FoldingCurve[i+1]->first->p - FoldingCurve[i]->first->p).normalized();
            double k = std::acos(v.dot(v2));
            if(Eigen::Vector3d::UnitZ().dot(v.cross(v2)) > 0)k = 2.0*std::numbers::pi - k;
        }
        auto ReviseEndPosition = [](Eigen::Vector3d e,  double k)->Eigen::Vector3d{
             return (Eigen::AngleAxisd(k, -Eigen::Vector3d::UnitZ()) * e).normalized();
        };
        auto getCrossPosition = [](const std::shared_ptr<Vertex>& v, const std::shared_ptr<Vertex>& o, const std::shared_ptr<Vertex>& p, const std::shared_ptr<Vertex>& q){
            Eigen::Matrix2d A;
            Eigen::Vector2d b;
            Eigen::Vector3d v2 = q->p - p->p;
            A(0,0) = v->p.x(); A(0,1) = -v2.x(); A(1,0) = v->p.y(); A(1,1) = -v2.y();
            b(0) = p->p.x() - o->p.x(); b(1) = p->p.y() - o->p.y();
            Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
            return  o->p + ts(0) * v->p;
        };

        double k0 = RevisionVertices::getK(FoldingCurve[1]->first->p, FoldingCurve[0]->first->p, FoldingCurve[2]->first->p);
        double k1 = RevisionVertices::getK(FoldingCurve[2]->first->p, FoldingCurve[1]->first->p, FoldingCurve[3]->first->p);
        double k2 = RevisionVertices::getK(FoldingCurve[3]->first->p, FoldingCurve[2]->first->p, FoldingCurve[4]->first->p);
        //std::cout <<"front " << abs(k0 - k1) << ", " << abs(k1 - k2) << " , " << 2*k1 - k2 << std::endl;
        if(abs(k0 - k1) > abs(k1 - k2)){
            std::cout <<"front"<<std::endl;
            FoldingCurve[0]->first->p = ReviseEndPosition(FoldingCurve[2]->first->p - FoldingCurve[1]->first->p,  -(2.*k1 - k2));
            Eigen::Vector3d newPos = getCrossPosition(FoldingCurve[0]->first, FoldingCurve[1]->first, FoldingCurve[0]->third, FoldingCurve[0]->second);
            FoldingCurve[0]->first->p = newPos;
            FoldingCurve[0]->first->p3 = (newPos - FoldingCurve[0]->third->p).norm()/(FoldingCurve[0]->third->p - FoldingCurve[0]->second->p).norm()
                    * (FoldingCurve[0]->second->p3 - FoldingCurve[0]->third->p3) + FoldingCurve[0]->third->p3;
        }
        k0 = RevisionVertices::getK(FoldingCurve.end()[-2]->first->p, FoldingCurve.end()[-3]->first->p, FoldingCurve.back()->first->p);
        k1 = RevisionVertices::getK(FoldingCurve.end()[-3]->first->p, FoldingCurve.end()[-4]->first->p, FoldingCurve.end()[-2]->first->p);
        k2 = RevisionVertices::getK(FoldingCurve.end()[-4]->first->p, FoldingCurve.end()[-5]->first->p, FoldingCurve.end()[-3]->first->p);
        //std::cout <<"back " << abs(k0 - k1) << ", " << abs(k1 - k2) << ", " << 2.*k1 - k2 <<  std::endl;
        if(abs(k0 - k1) > abs(k1 - k2)){
            std::cout << "back"<<std::endl;
            FoldingCurve.back()->first->p = ReviseEndPosition(FoldingCurve.end()[-3]->first->p - FoldingCurve.end()[-2]->first->p, 2.*k1 - k2);
            Eigen::Vector3d newPos2 = getCrossPosition(FoldingCurve.back()->first, FoldingCurve.end()[-2]->first, FoldingCurve.back()->third, FoldingCurve.back()->second);
            FoldingCurve.back()->first->p = newPos2;
            FoldingCurve.back()->first->p3 = (newPos2 - FoldingCurve.back()->third->p).norm()/
                    (FoldingCurve.back()->third->p - FoldingCurve.back()->second->p).norm() * (FoldingCurve.back()->second->p3 - FoldingCurve.back()->third->p3) + FoldingCurve.back()->third->p3;
        }

    }
    return true;

}

void FoldLine::ReassignColor(std::vector<std::shared_ptr<Line>>& Rulings, ColorPoint& CP){
    //Y字のとき折り線は山であると仮定(谷の場合は色を反転)
    typedef struct {
        std::shared_ptr<Line> r;
        int type_mvk;//0:mountain && k < pi, 1: mountain && k >= pi, 2: vally && k < pi, 3: vally && k >= pi, -1: gradation = 0
        //0の数 + 3の数 == rulingの数 || 1の数 + 2の数 == rulingの数 -> 山谷とkの割り当てが正しい
        //それ以外
         //0の数 > 2の数 -> 　2の色を変換 (逆もまた然り)(1と3でも同様に)
    }MVK;
    std::vector<MVK> InitState;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}

    std::vector<int>MV(5,0);
    for(int i = 1; i < (int)ValidFC.size()-1; i++){
        auto e = (ValidFC[i-1]->first->p - ValidFC[i]->first->p).normalized(),
        e2 = (ValidFC[i+1]->first->p - ValidFC[i]->first->p).normalized();
        double k = std::acos(e.dot(e2));
        if(Eigen::Vector3d::UnitZ().dot(e.cross(e2)) > 0)k = 2.0*std::numbers::pi - k;
       Eigen::Vector3d N = ((ValidFC[i]->third->p3 - ValidFC[i]->first->p3).cross(e)).normalized();
       Eigen::Vector3d Np = (e2.cross(ValidFC[i]->third->p3 - ValidFC[i]->first->p3)).normalized();
       double phi = std::acos(N.dot(Np));

        double color;
        if(phi < CP.angle)color = CP.color/CP.angle *phi;
        else color = ((255.0 - CP.color)*(phi - CP.angle)/(std::numbers::pi - CP.angle) + CP.color);
        for(auto&r: Rulings){
            if(r->et == EdgeType::r && MathTool::is_point_on_line(ValidFC[i]->first->p, r->v->p, r->o->p)){
                int mv = (r->color == 0)? -1: (r->color > 0)? 0: 1;
                int type_mvk = (mv == 0 && k < std::numbers::pi)? 0: (mv == 0 && k >= std::numbers::pi)? 1: (mv == 1 && k < std::numbers::pi)? 2: (mv == 1 && k >= std::numbers::pi)? 3: -1;
                InitState.push_back(MVK(r, type_mvk));
                MV[type_mvk + 1] += 1;
                break;
            }
        }
    }

    if(MV[0] != 0){std::cout << "no folding pattern was found"<< std::endl; return;}
    while(!(MV[1] + MV[4] == (int)ValidFC.size() - 2 || MV[2] + MV[3] == (int)ValidFC.size() - 2)){
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

template <typename T>
class PointOnEndEdge{
    public:
    std::shared_ptr<T> v;
    double t;
    PointOnEndEdge(std::shared_ptr<T>&_v, double _t): v(_v), t(_t){}
};

void FoldLine::reassignruling(std::shared_ptr<FoldLine>& parent){
    auto getCrossPoint = [](std::vector<Eigen::Vector3d>& CtrlPts, std::shared_ptr<Vertex> v, std::shared_ptr<Vertex> o, int dim){
        std::vector<double>arcT = BezierClipping(CtrlPts, v, o, dim);
        for(auto&t: arcT){
            if(t < 0 || 1 < t){std::cout<<"t is not correct value " << t << std::endl; continue;}
            Eigen::Vector3d v2(0,0,0);
            for (int i = 0; i < int(CtrlPts.size()); i++) v2 += MathTool::BernsteinBasisFunc(dim, i, t) * CtrlPts[i];
            if(!MathTool::is_point_on_line(v2, v->p, o->p))continue;
            double s = (v2 - o->p).norm()/(o->p - v->p).norm();
            Eigen::Vector3d v3 = s * (v->p3 - o->p3) + o->p3;
            std::shared_ptr<CrvPt_FL> P = std::make_shared<CrvPt_FL>(CrvPt_FL(v2, v3, t));
            return P;
        }
        return std::shared_ptr<CrvPt_FL>(nullptr);
    };
    if(parent->FoldingCurve.empty() || FoldingCurve.empty())return;
    int dim = CtrlPts.size() - 1;
    auto SplitOnEndLine = [&](std::shared_ptr<Vertex4d>& V, std::shared_ptr<Vertex4d>& line_v, std::vector<std::shared_ptr<Vertex4d>>& FC, int dim){
        std::shared_ptr<CrvPt_FL> p = getCrossPoint(CtrlPts, line_v->first, line_v->second, dim);
        if(p == nullptr)return;
        std::vector<PointOnEndEdge<Vertex>> PointOnEndEdges;
        std::shared_ptr<Vertex> dowscastV = std::dynamic_pointer_cast<Vertex>(line_v->first);
        PointOnEndEdges.push_back(PointOnEndEdge(dowscastV, 0.0));
        PointOnEndEdges.push_back(PointOnEndEdge(line_v->second, (line_v->second->p - line_v->first->p).norm()));
        for(auto it = FC.begin() + 1; it != FC.end() - 1; it++){
            if(MathTool::is_point_on_line((*it)->second->p, line_v->first->p, line_v->second->p)){
                PointOnEndEdges.push_back(PointOnEndEdge((*it)->second, ((*it)->second->p - line_v->first->p).norm()));
            }
        }
        std::sort(PointOnEndEdges.begin(), PointOnEndEdges.end(), [](const PointOnEndEdge<Vertex>& V1, const PointOnEndEdge<Vertex>& V2){return V1.t < V2.t;});
        for(auto it = PointOnEndEdges.begin() + 1; it != PointOnEndEdges.end(); it++){
            if(MathTool::is_point_on_line(p->p, (it - 1)->v->p, it->v->p)){
                V->first = p; V->third = (it - 1)->v;
                V->second = line_v->second;//V->second = it->v;
                line_v->second = p;
                return;
            }
        }
    };

    SplitOnEndLine(FoldingCurve.front(), parent->FoldingCurve.front(), parent->FoldingCurve, dim);
    SplitOnEndLine(FoldingCurve.back(),  parent->FoldingCurve.back(), parent->FoldingCurve, dim);
    for(auto it = parent->FoldingCurve.begin() + 1; it != parent->FoldingCurve.end() - 1; it++){
        if(!(*it)->IsCalc)continue;
        std::shared_ptr<CrvPt_FL> p = getCrossPoint(CtrlPts, (*it)->first, (*it)->second, dim);
        if(p != nullptr){
            FoldingCurve.push_back(std::make_shared<Vertex4d>(p, (*it)->second, (*it)->first));
            (*it)->second = p;
        }
    }
    validsize = FoldingCurve.size();
    SortCurve();
    parent->SortCurve();
}

void FoldLine::drawRulingInAllAngles(std::vector<std::array<Eigen::Vector3d, 2>>& _Rulings){

    double a2 = 0;
    Eigen::Vector3d e, e2, x;
    Eigen::Vector3d N, r, befN;
    double veclen = 40;
    std::array<Eigen::Vector3d, 2> newVec;
    _Rulings.clear();
    double angle = 0.0;
    int flsize = FoldingCurve.size();
    while(angle <= 2*std::numbers::pi){
        double a = angle;
        double tau, dir;
        for(int i = 1; i < flsize - 1; i++){
            x = FoldingCurve[i]->first->p3;
            e = (FoldingCurve[i-1]->first->p3 - x).normalized();
            e2 = (FoldingCurve[i+1]->first->p3 - x).normalized();

            double beta = std::acos(e.dot(e2));
            if(i != 1){
                Eigen::Vector3d srfN = MathTool::ProjectionVector(e.cross(e2), -e, true);
                dir = (-e).dot(befN.cross(srfN));
                tau = std::acos(srfN.dot(befN));
                if(dir <= 0)tau = 2.0*std::numbers::pi - tau;
                a = (a2 - tau < 0) ? a2 - tau + 2.0 * std::numbers::pi: a2 - tau;
            }
            Eigen::Vector3d Axis = (FoldingCurve[i]->third->p3 - FoldingCurve[i]->first->p3).normalized();
            double phi3 = std::acos(Axis.dot((FoldingCurve[i+1]->first->p3 - FoldingCurve[i]->first->p3).normalized()));
            double phi4 = std::acos(Axis.dot((FoldingCurve[i-1]->first->p3 - FoldingCurve[i]->first->p3).normalized()));
            double k = 2.0 * std::numbers::pi - phi3 - phi4;
            double _phi1 = atan((std::cos(k) - std::cos(beta))/(std::sin(beta)*std::cos(a)- std::sin(k)));
            double phi1 = (_phi1 < 0)? _phi1+std::numbers::pi: _phi1;
            double phi2 = k - phi1;
            N = (Eigen::AngleAxisd(a, -e)  * e2.cross(e)).normalized();
            r = Eigen::AngleAxisd(phi1, -N) * e;//新しいruling方向
            newVec = std::array<Eigen::Vector3d, 2>{veclen * r + x, x};_Rulings.push_back(newVec);

            double sin_a = sin(phi1)*sin(a)/sin(phi2);
            if(sin_a > 1)sin_a = 1;
            else if(sin_a < -1)sin_a = -1;
            double cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
            if(cos_a > 1)cos_a = 1;
            else if(cos_a < -1)cos_a = -1;
            a2 = (sin_a >= 0)? acos(cos_a): (sin_a < 0 && cos_a < 0)? 2.0*std::numbers::pi - acos(cos_a): 2.0*std::numbers::pi + asin(sin_a);
            if(i == (int)FoldingCurve.size() -2)break;
            befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
        }
        angle += 1e-3;
    }

}

void FoldLine::applyAAAMethod(const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool begincener, double a, double _tol, bool isroot){
    if(FoldingCurve.empty() && a_flap == -1)return;
    SimplifyModel(validsize, isroot);
    //if(!begincener)_FoldingAAAMethod(FoldingCurve, Poly_V, a);
    //else _FoldingAAAMethod_center(FoldingCurve, Poly_V, a);
    _FoldingAAAMethod(FoldingCurve, Poly_V, a);
    if(DebugMode::Singleton::getInstance().isdebug())std::cout << a << "  ,,  " << RulingsCrossed(FoldingCurve) << " , " <<RulingsCrossed_NoOutline(FoldingCurve) << std::endl;
    if(RulingsCrossed(FoldingCurve) < 1e-9){
        a_flap = a; tol = _tol;
    }
    return;
}

inline double update_flapangle(double a, const Eigen::Vector3d& befN, const Eigen::Vector3d& SrfN, const Eigen::Vector3d& e){
    double dir = (-e).dot(befN.cross(SrfN));
    double tau = std::acos(SrfN.dot(befN));
    if(dir <= 0)tau = 2.0*std::numbers::pi - tau;
    if(a - tau <= 0)return a - tau + 2.0 * std::numbers::pi;
    return a - tau;
}

inline Eigen::Vector3d _calcruling3d(const double& a, Eigen::Vector3d e, Eigen::Vector3d e2, Eigen::Vector3d axis, double& beta, std::vector<double>& Phi){
    e = e.normalized(); e2 = e2.normalized(); axis = axis.normalized();
    beta = std::acos(e.dot(e2));
    Phi.resize(4);
    Phi[2] = std::acos(e2.dot(axis)); Phi[3] = std::acos(e.dot(axis));

    double k = 2.0 * std::numbers::pi - Phi[2] - Phi[3];
    double _phi1 = std::atan((std::cos(k) - std::cos(beta))/(std::sin(beta)*std::cos(a)- std::sin(k)));
    Phi[0] = (_phi1 < 0)? _phi1+std::numbers::pi: _phi1;
    Phi[1] = k - Phi[0];
    return RevisionVertices::decideRulingDirectionOn3d(e, e.cross(e2), a, Phi[0]);//新しいruling方向
}

void CalcRuling(double a, std::shared_ptr<Vertex4d>& xbef, std::shared_ptr<Vertex4d>& x, std::shared_ptr<Vertex4d>& xnext, const std::vector<std::shared_ptr<Vertex>>& Poly_V,  double& a2, Eigen::Vector3d& SrfN, const Eigen::Vector3d &SpinAxis ){

    Eigen::Vector3d e = (xbef->first->p3 - x->first->p3).normalized(), e2 = (xnext->first->p3 - x->first->p3).normalized();
    double beta;
    std::vector<double> Phi;
    Eigen::Vector3d Axis = (x->third->p3 - x->first->p3).normalized();
    Eigen::Vector3d r3d = _calcruling3d(a, e, e2, Axis, beta, Phi);
    Eigen::Vector3d r2d = (Eigen::AngleAxisd(Phi[0], SpinAxis)* (xbef->first->p- x->first->p)).normalized();//展開図のruling方向

    Eigen::Vector3d crossPoint;
    for(int k = 0; k < (int)Poly_V.size(); k++){
        crossPoint = MathTool::calcCrossPoint_2Vector(x->first->p, 1000.0 * r2d + x->first->p, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p);
        if(MathTool::is_point_on_line(crossPoint, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) && r2d.dot(crossPoint - x->first->p) > 0)break;
    }
    x->second->p = crossPoint;
    x->second->p3 = (crossPoint - x->first->p).norm() * r3d + x->first->p3;

    double sin_a = sin(Phi[0])*sin(a)/sin(Phi[1]);
    if(sin_a > 1)sin_a = 1;
    else if(sin_a < -1)sin_a = -1;
    double cos_a = (cos(Phi[0]) - cos(Phi[1])*cos(beta))/(sin(Phi[1])*sin(beta));
    if(cos_a > 1)cos_a = 1;
    else if(cos_a < -1)cos_a = -1;
    a2 = std::atan2(sin_a, cos_a);
    SrfN = MathTool::ProjectionVector(e.cross(e2), e2, true);
}

std::shared_ptr<Vertex> getClosestVertex(const std::shared_ptr<Vertex>& v, const std::shared_ptr<Vertex>& o,  const std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, bool SkipTrimedPoint){
   std::shared_ptr<Vertex> V_max = nullptr;
   double t_max = -1;
   for(auto&fc: FoldingCurve){
       if(SkipTrimedPoint && !fc->IsCalc)continue;
       if(((fc->second->p - v->p).norm() < 1e-7|| (fc->second->p - o->p).norm() < 1e-7))continue;
       if(MathTool::is_point_on_line(fc->second->p, v->p, o->p)){
           double t = (fc->second->p - o->p).norm()/(v->p - o->p).norm();
           if(t > t_max){
               t_max = t; V_max = fc->second;
           }
       }
   }
   return V_max;
}

void _FoldingAAAMethod_center(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle){

    auto SetOnPlane = [&](std::shared_ptr<Vertex4d>& V, const std::vector<std::shared_ptr<Vertex>>& Poly_V, std::shared_ptr<Vertex4d>& p, std::shared_ptr<Vertex4d>& p2, const std::shared_ptr<Vertex>& q, const std::shared_ptr<Vertex>& q2, int IsEnd){
        //p, q: 左、p2, q2：右
        Eigen::Vector3d vec;
        if(IsEnd == -1){vec = (q2->p - p2->first->p).normalized();}//左が端
        else if(IsEnd == 1){vec = (q->p - p->first->p).normalized();}//右が端
        else{
            if(!IsParallel(p->first, q, p2->first, q2))vec = (calcCrossPoint_2Vertex(p->first, q, p2->first, q2) - V->first->p).normalized();
            else vec = (q->p - p->first->p).normalized();
            auto e = (p2->first->p - V->first->p).normalized(), e2 = (p->first->p - V->first->p).normalized();
            double k = e.dot(e2); k = (k < -1)? std::numbers::pi: (k > 1)? 0: std::acos(k);
            Eigen::Vector3d N = e.cross(e2);
            if(N.z() > 0)k = 2.0*std::numbers::pi - k;
            double phi = e.dot(vec); phi = (phi < -1)? std::numbers::pi: (phi > 1)? 0: std::acos(phi);
            N = e.cross(vec);if(N.z() > 0)phi = 2.0*std::numbers::pi - phi;
            if(phi > k)vec *= -1;
        }
        vec = 1000. * vec + V->first->p;
        std::shared_ptr<Vertex> Qv = std::make_shared<Vertex>(Vertex(vec));

        for(int k = 0; k < (int)Poly_V.size(); k++){
            Eigen::Vector3d q;
            if(IsParallel(V->first, Qv, Poly_V[k], Poly_V[(k + 1) % (int)Poly_V.size()]))continue;
            q = calcCrossPoint_2Vertex(V->first, Qv, Poly_V[k], Poly_V[(k + 1) % (int)Poly_V.size()]);
            if(MathTool::is_point_on_line(q, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) && MathTool::is_point_on_line(q, V->first->p, Qv->p) ){
                V->second->p = q;
                V->second->p3 = calcTargetDistanceOnPlane(V->second->p, q2, p->first, p2->first);
            }
        }
    };

    double a2 = 0, l;
    Eigen::Vector3d e, e2, x;
    Eigen::Vector3d N, r, befN;
    std::vector<int> Vertices_Ind;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i]->IsCalc)Vertices_Ind.push_back(i);}

    double phi02 = std::acos(((FoldingCurve[Vertices_Ind[0]]->second->p - FoldingCurve[Vertices_Ind[0]]->first->p).normalized()).dot(
        (FoldingCurve[Vertices_Ind[1]]->first->p - FoldingCurve[Vertices_Ind[0]]->first->p).normalized()));
    double phim1 = std::acos(((FoldingCurve[Vertices_Ind.end()[-2]]->first->p - FoldingCurve.back()->first->p).normalized()).dot(
        (FoldingCurve.back()->second->p - FoldingCurve.back()->first->p).normalized()));

    double a = (angle < 0.0)? angle + 2.0 * std::numbers::pi: (angle > 2.0*std::numbers::pi)? angle - 2.0*std::numbers::pi: angle;
    std::shared_ptr<Vertex4d> fc, fc_bef, fc_next;
    for(int ind = Vertices_Ind.size()/2; ind < (int)Vertices_Ind.size() - 1; ind++){
        fc = FoldingCurve[Vertices_Ind[ind]]; fc_bef = FoldingCurve[Vertices_Ind[ind - 1]]; fc_next = FoldingCurve[Vertices_Ind[ind + 1]];
        x = fc->first->p3;
        e = (fc_bef->first->p3 - x).normalized();
        e2 = (fc_next->first->p3 - x).normalized();
        if(ind != (int)Vertices_Ind.size()/2){
            Eigen::Vector3d SrfN = MathTool::ProjectionVector(e.cross(e2), -e, true);
            a = update_flapangle(a2, befN, SrfN, e);
        }
        CalcRuling(a, fc_bef, fc, fc_next, Poly_V, a2, befN, -Eigen::Vector3d::UnitZ());
        if(ind == (int)FoldingCurve.size() -2)break;
        befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
        if(ind != 1){
            for(int i = Vertices_Ind[ind-1] + 1; i < Vertices_Ind[ind]; i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind[ind]],FoldingCurve[Vertices_Ind[ind-1]], FoldingCurve[Vertices_Ind[ind]]->second,  FoldingCurve[Vertices_Ind[ind-1]]->second, 0);
        }
    }
    a = (angle < 0.0)? angle + 2.0 * std::numbers::pi: (angle > 2.0*std::numbers::pi)? angle - 2.0*std::numbers::pi: angle;
    for(int ind = Vertices_Ind.size()/2 -1; ind > 0; ind--){
        fc = FoldingCurve[Vertices_Ind[ind]]; fc_bef = FoldingCurve[Vertices_Ind[ind + 1]]; fc_next = FoldingCurve[Vertices_Ind[ind - 1]];
        x = fc->first->p3;
        e = (fc_bef->first->p3 - x).normalized();
        e2 = (fc_next->first->p3 - x).normalized();
        if(ind != (int)Vertices_Ind.size()/2 -1){Eigen::Vector3d SrfN = MathTool::ProjectionVector(e.cross(e2), -e, true);a = update_flapangle(a2, befN, SrfN, e);}
        CalcRuling(a, fc_bef, fc, fc_next, Poly_V, a2, befN, Eigen::Vector3d::UnitZ());
        if(ind == (int)FoldingCurve.size() -2)break;
        befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
        if(ind != 1){
            for(int i = Vertices_Ind[ind-1] + 1; i < Vertices_Ind[ind]; i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind[ind]],FoldingCurve[Vertices_Ind[ind-1]], FoldingCurve[Vertices_Ind[ind]]->second,  FoldingCurve[Vertices_Ind[ind-1]]->second, 0);
        }
    }
}


void _FoldingAAAMethod(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle){

    auto SetOnPlane = [&](std::shared_ptr<Vertex4d>& V, const std::vector<std::shared_ptr<Vertex>>& Poly_V, std::shared_ptr<Vertex4d>& p, std::shared_ptr<Vertex4d>& p2, const std::shared_ptr<Vertex>& q, const std::shared_ptr<Vertex>& q2, int IsEnd){
        //p, q: 左、p2, q2：右
        Eigen::Vector3d vec;
        if(IsEnd == -1){vec = (q2->p - p2->first->p).normalized();}//左が端
        else if(IsEnd == 1){vec = (q->p - p->first->p).normalized();}//右が端
        else{
            if(!IsParallel(p->first, q, p2->first, q2))vec = (calcCrossPoint_2Vertex(p->first, q, p2->first, q2) - V->first->p).normalized();
            else vec = (q->p - p->first->p).normalized();
            auto e = (p2->first->p - V->first->p).normalized(), e2 = (p->first->p - V->first->p).normalized();
            double k = e.dot(e2); k = (k < -1)? std::numbers::pi: (k > 1)? 0: std::acos(k);
            auto N = e.cross(e2);
            if(N.z() > 0)k = 2.0*std::numbers::pi - k;
            double phi = e.dot(vec); phi = (phi < -1)? std::numbers::pi: (phi > 1)? 0: std::acos(phi);
            N = e.cross(vec);if(N.z() > 0)phi = 2.0*std::numbers::pi - phi;
            if(phi > k)vec *= -1;
        }
        vec = 1000. * vec + V->first->p;
        std::shared_ptr<Vertex> Qv = std::make_shared<Vertex>(Vertex(vec));

        for(int k = 0; k < (int)Poly_V.size(); k++){
            Eigen::Vector3d q;
            if(IsParallel(V->first, Qv, Poly_V[k], Poly_V[(k + 1) % (int)Poly_V.size()]))continue;
            q = calcCrossPoint_2Vertex(V->first, Qv, Poly_V[k], Poly_V[(k + 1) % (int)Poly_V.size()]);
            if(MathTool::is_point_on_line(q, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) && MathTool::is_point_on_line(q, V->first->p, Qv->p) ){
                V->second->p = q;
                V->second->p3 = calcTargetDistanceOnPlane(V->second->p, q2, p->first, p2->first);
            }
        }
    };

    double a2 = 0, l;
    Eigen::Vector3d e, e2, x;
    Eigen::Vector3d N, r, befN;
    std::vector<int> Vertices_Ind;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i]->IsCalc)Vertices_Ind.push_back(i);}

    double phi02 = std::acos(((FoldingCurve[Vertices_Ind[0]]->second->p - FoldingCurve[Vertices_Ind[0]]->first->p).normalized()).dot(
        (FoldingCurve[Vertices_Ind[1]]->first->p - FoldingCurve[Vertices_Ind[0]]->first->p).normalized()));
    double phim1 = std::acos(((FoldingCurve[Vertices_Ind.end()[-2]]->first->p - FoldingCurve.back()->first->p).normalized()).dot(
        (FoldingCurve.back()->second->p - FoldingCurve.back()->first->p).normalized()));

    double a = (angle < 0.0)? angle + 2.0 * std::numbers::pi: (angle > 2.0*std::numbers::pi)? angle - 2.0*std::numbers::pi: angle;
    for(int ind = 1; ind < (int)Vertices_Ind.size() - 1; ind++){
        std::shared_ptr<Vertex4d> fc = FoldingCurve[Vertices_Ind[ind]];
        std::shared_ptr<Vertex4d> fc_bef = FoldingCurve[Vertices_Ind[ind - 1]];
        std::shared_ptr<Vertex4d> fc_next = FoldingCurve[Vertices_Ind[ind + 1]];

        x = fc->first->p3;
        e = (fc_bef->first->p3 - x).normalized();
        e2 = (fc_next->first->p3 - x).normalized();
        if(ind != 1){
            Eigen::Vector3d SrfN = MathTool::ProjectionVector(e.cross(e2), -e, true);
            a = update_flapangle(a2, befN, SrfN, e);
        }
        CalcRuling(a, fc_bef, fc, fc_next, Poly_V, a2, befN, -Eigen::Vector3d::UnitZ());
        if(ind == (int)FoldingCurve.size() -2)break;
        befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
        if(ind != 1){
            for(int i = Vertices_Ind[ind-1] + 1; i < Vertices_Ind[ind]; i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind[ind]],FoldingCurve[Vertices_Ind[ind-1]], FoldingCurve[Vertices_Ind[ind]]->second,  FoldingCurve[Vertices_Ind[ind-1]]->second, 0);
        }
    }
    {//i = 0(端の面)
        std::shared_ptr<Vertex> v_clst = getClosestVertex(FoldingCurve[Vertices_Ind[0]]->second , FoldingCurve[Vertices_Ind[0]]->first, FoldingCurve);
        if(v_clst == nullptr){
            FoldingCurve[Vertices_Ind[0]]->second->p3 = calcTargetDistanceOnPlane(FoldingCurve[Vertices_Ind[0]]->second->p, FoldingCurve[Vertices_Ind[1]]->first, FoldingCurve[Vertices_Ind[0]]->first, FoldingCurve[Vertices_Ind[1]]->second);
            for(int i = Vertices_Ind[0] + 1; i < Vertices_Ind[1]; i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind[1]], FoldingCurve[Vertices_Ind[0]], FoldingCurve[Vertices_Ind[1]]->second, FoldingCurve[Vertices_Ind[0]]->second, 1);
        }else{
            int j;
            for(j = 0; j < (int)Vertices_Ind.size(); j++){
                if(FoldingCurve[Vertices_Ind[j]]->second->p == v_clst->p)break;
            }
            if(j == (int)Vertices_Ind.size())std::cout <<"can't find correct edge right" << std::endl;
            else{
                auto p = calcTargetDistanceOnPlane(FoldingCurve[Vertices_Ind[0]]->second->p, v_clst, FoldingCurve[Vertices_Ind[j]]->first, FoldingCurve[Vertices_Ind[j+1]]->first);
                FoldingCurve[Vertices_Ind[0]]->second->p3 = p;
            }
            for(int i = Vertices_Ind[0] + 1; i < Vertices_Ind[1]; i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind[1]], FoldingCurve[Vertices_Ind[0]], FoldingCurve[Vertices_Ind[1]]->second, v_clst, 1);
        }

    }

    {
        std::shared_ptr<Vertex> v_clst = getClosestVertex(FoldingCurve[Vertices_Ind.back()]->second , FoldingCurve[Vertices_Ind.back()]->first, FoldingCurve);
        if(v_clst == nullptr){
            FoldingCurve[Vertices_Ind.back()]->second->p3 = calcTargetDistanceOnPlane(FoldingCurve[Vertices_Ind.back()]->second->p, FoldingCurve[Vertices_Ind.end()[-2]]->first, FoldingCurve[Vertices_Ind.back()]->first, FoldingCurve[Vertices_Ind.end()[-2]]->second);
            for(int i = Vertices_Ind.end()[-2] + 1; i < Vertices_Ind.back(); i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind.back()], FoldingCurve[Vertices_Ind.end()[-2]], FoldingCurve[Vertices_Ind.back()]->second, FoldingCurve[Vertices_Ind.end()[-2]]->second, -1);
        }else{
            int j;
            for(j = 0; j < (int)Vertices_Ind.size(); j++){
                if(FoldingCurve[Vertices_Ind[j]]->second->p == v_clst->p)break;
            }
            if(j == (int)Vertices_Ind.size())std::cout <<"can't find correct edge right" << std::endl;
            else{
                auto p = calcTargetDistanceOnPlane(FoldingCurve[Vertices_Ind.back()]->second->p, v_clst,  FoldingCurve[Vertices_Ind[j]]->first, FoldingCurve[Vertices_Ind[j-1]]->first);
                FoldingCurve[Vertices_Ind.back()]->second->p3 = p;
            }
            for(int i = Vertices_Ind.end()[-2] + 1; i < Vertices_Ind.back(); i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind.back()], FoldingCurve[Vertices_Ind.end()[-2]], v_clst, FoldingCurve[Vertices_Ind.end()[-2]]->second, -1);
        }

    }
}

//kokoga genninn
std::vector<std::shared_ptr<Vertex4d>> TrimPoints(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, double tol){
    std::vector<std::shared_ptr<Vertex4d>> res;
    size_t st = FoldingCurve.size();

    while(1){
        Douglas_Peucker_algorithm(FoldingCurve,res, tol);
        if(res.size() == st)break;
        st = res.size();
    }
    return res;

}

//https://issekinichou.wordpress.com/2022/12/28/douglas-peucker-algorithm/
void Douglas_Peucker_algorithm(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, std::vector<std::shared_ptr<Vertex4d>>& res, double tol){
    auto perpendicularDistance = [](Eigen::Vector3d& p, Eigen::Vector3d& l_start, Eigen::Vector3d& l_end)->double{
        Eigen::Vector3d line = (l_end - l_start).normalized(), OP = (p - l_start).normalized();
        Eigen::Vector3d H = l_start + line.dot(OP) * (p - l_start);
        return (H - p).norm();
    };
    if((int)FoldingCurve.size() < 2)return;
    double dmax = 0.0;
    size_t index = 0;
    size_t end = FoldingCurve.size()-1;
    for(size_t i = 1; i < end; i++){
        double d = perpendicularDistance(FoldingCurve[i]->first->p2_ori, FoldingCurve[0]->first->p2_ori, FoldingCurve.back()->first->p2_ori);
        if (d > dmax){index = i; dmax = d;}
    }

    // If max distance is greater than epsilon, recursively simplify
    if(dmax > tol){
        // Recursive call
        std::vector<std::shared_ptr<Vertex4d>> res_first, res_last;
        std::vector<std::shared_ptr<Vertex4d>> firstLine{FoldingCurve.begin(), FoldingCurve.begin()+index+1};
        std::vector<std::shared_ptr<Vertex4d>> lastLine{FoldingCurve.begin() + index, FoldingCurve.end()};
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
