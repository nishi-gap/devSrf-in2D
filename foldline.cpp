#include "foldline.h"

const double w_at = 1e-4;
const double eps = 1e-7;
inline Eigen::Vector3d _calcruling3d(const double& a, Eigen::Vector3d e, Eigen::Vector3d e2, Eigen::Vector3d axis, double& beta, std::vector<double>& Phi);
void CalcRuling(double a, std::shared_ptr<Vertex4d>& xbef, std::shared_ptr<Vertex4d>& x, std::shared_ptr<Vertex4d>& xnext, const std::vector<std::shared_ptr<Vertex>>& Poly_V, double& a2, Eigen::Vector3d& SrfN, const Eigen::Vector3d& SpinAixs);
void CalcRuling_back(double a, std::shared_ptr<Vertex4d>& xbef, std::shared_ptr<Vertex4d>& x, std::shared_ptr<Vertex4d>& xnext, const std::vector<std::shared_ptr<Vertex>>& Poly_V, double& a2, Eigen::Vector3d& SrfN, const Eigen::Vector3d& SpinAixs);
inline double update_flapangle(double a, const Eigen::Vector3d& befN, const Eigen::Vector3d& SrfN, const Eigen::Vector3d& e);
inline double update_flapangle2(double a, const Eigen::Vector3d& befN, const Eigen::Vector3d& SrfN, const Eigen::Vector3d& e);
void _FoldingAAAMethod(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle);
void _FoldingAAAMethod_center(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle);
std::vector<std::shared_ptr<Vertex4d>> TrimPoints(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, double tol);
bool IsRulingCrossed(Eigen::Vector3d N, Eigen::Vector3d& cp, Eigen::Vector3d& crossPoint,  const std::vector<std::shared_ptr<Vertex>>& Poly_V);
void Douglas_Peucker_algorithm(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, std::vector<std::shared_ptr<Vertex4d>>& res, double tol = std::numbers::pi/9.0);
std::shared_ptr<Vertex> getClosestVertex(const std::shared_ptr<Vertex>& v, const std::shared_ptr<Vertex>& o,  const std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, bool SkipTrimedPoint = true);

inline double logistic(double x, double a = 1){return 1.0/(1+exp(-a*x));}

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

    class FrenetFrame{
        public:
        std::vector<Eigen::Vector3d> T, N, B, R;
        std::vector<double> Beta, Tau, Kappa;
        FrenetFrame(std::vector<Eigen::Vector3d>& _T, std::vector<Eigen::Vector3d>& _N, std::vector<Eigen::Vector3d>& _B, std::vector<Eigen::Vector3d>& _R,
                    std::vector<double>& _beta, std::vector<double>_tau, std::vector<double>& _kappa): T(_T), N(_N), B(_B), Beta(_beta), Tau(_tau), Kappa(_kappa), R(_R) {}
    };

    struct OptimizeParam{
        FoldLine3d FC;
        std::vector<std::shared_ptr<Vertex>> Poly_V;
        std::vector<double*> res_Fbend, res_Fruling, res_a;
        FoldLine3d BasePt;
        double wb, wp ,wr;
        void AddWeight(double _wb, double _wp, double _wr){ wb = _wb; wp = _wp; wr = _wr;}
        void AddBasePt(std::vector<std::shared_ptr<Vertex4d>>& FC){BasePt = FC;}
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
    Eigen::Vector3d decideRulingDirectionOn3d(Eigen::Vector3d e, Eigen::Vector3d N, double a, double phi);

    Eigen::Vector3d getIntersection(const std::shared_ptr<Vertex4d>& v, const std::shared_ptr<Vertex4d>& v2){//v2-vのつながりとして与えないとうまくいかないかも
        Eigen::Vector3d e = (v->first->p - v2->first->p).normalized(), r = (v->second->p - v->first->p).normalized(), r2 = (v2->second->p - v2->first->p).normalized();
        double phi = std::acos(r.dot(-e)), phi2 = std::acos(r2.dot(e));
        double l = sin(phi2)/sin(phi+phi2)*(v->first->p - v2->first->p).norm();
        return l*(v->second->p - v->first->p).normalized() + v->first->p;
    }
    std::vector<Eigen::Vector3d> getRegCrvPt(const std::vector<std::shared_ptr<Vertex4d>>&FC){
        std::vector<Eigen::Vector3d> Pts;
        std::vector<std::shared_ptr<Vertex4d>> ValidFC;
        for(auto&fc: FC){if(fc->IsCalc)ValidFC.push_back(fc);}
        for(int i = 0; i < (int)ValidFC.size() - 1; i++){
            Eigen::Vector3d p = getIntersection(ValidFC[i], ValidFC[i+1]);
            Pts.push_back(p);
        }
        return Pts;
    }
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

    //erase closed or overlapped vertex
    FoldingCurve.erase(std::unique(FoldingCurve.begin(), FoldingCurve.end(), [&](std::shared_ptr<Vertex4d>&a, std::shared_ptr<Vertex4d>& b){return (a->first->p - b->first->p).norm() < 1e-9;}), FoldingCurve.end());

    if(FoldingCurve.size() >= 3){
        if((FoldingCurve[0]->second->p - FoldingCurve[0]->first->p).normalized().dot((FoldingCurve[1]->second->p - FoldingCurve[1]->first->p).normalized()) < 0)
            std::swap(FoldingCurve[0]->second, FoldingCurve[0]->third);
        if((FoldingCurve.back()->second->p - FoldingCurve.back()->first->p).normalized().dot((FoldingCurve.end()[-2]->second->p - FoldingCurve.end()[-2]->first->p).normalized()) < 0)
            std::swap(FoldingCurve.back()->second, FoldingCurve.back()->third);

    }
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
        double w = 0.0;
        double _f = 0.0;
        for(int j = -1; j <= 1; j +=2){
            int i2 = i + j;
            if(i2 <= 0 || i2 >= (int)ValidFC.size() - 1)continue;
            //Eigen::Vector3d p1 = ValidFC[i]->first->p, q1 = ValidFC[i]->second->p, p2 = ValidFC[i2]->first->p, q2 = ValidFC[i2]->second->p;
            //double t = ((p2.x() - p1.x())*(p2.y() - q2.y()) - (p2.x() - q2.x())*(p2.y() - p1.y()))/((q1.x() - p1.x()) * (p2.y() - q2.y()) - (p2.x() - q2.x())*(q1.y() - p1.y()));
            //if(0 < t && t < 1) f += abs(1.0/t);
            //continue;
            A(0,0) = ValidFC[i]->first->p.x() - ValidFC[i]->second->p.x();
            A(0,1) = -(ValidFC[i2]->first->p.x() - ValidFC[i2]->second->p.x());
            A(1,0) = ValidFC[i]->first->p.y() - ValidFC[i]->second->p.y();
            A(1,1) = -(ValidFC[i2]->first->p.y() - ValidFC[i2]->second->p.y());
            b(0) = ValidFC[i]->first->p.x() - ValidFC[i2]->first->p.x();
            b(1) = ValidFC[i]->first->p.y() - ValidFC[i2]->first->p.y();
            Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
            if(0 < ts(0) && ts(0) < 1){_f += 1.0/ts(0); w++;}
        }
        f += w * _f;
    }
    return f;
}

double RulingCrossed_all(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve){
    double f = 0.0;
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    for(int i = 1; i < (int)ValidFC.size() -1; i++){
        double w = 0.0;
        double _f = 0.0;
        for(int j = 1; j < (int)ValidFC.size() -1; j++){
            if(i == j)continue;
            A(0,0) = ValidFC[i]->first->p.x() - ValidFC[i]->second->p.x();
            A(0,1) = -(ValidFC[j]->first->p.x() - ValidFC[j]->second->p.x());
            A(1,0) = ValidFC[i]->first->p.y() - ValidFC[i]->second->p.y();
            A(1,1) = -(ValidFC[j]->first->p.y() - ValidFC[j]->second->p.y());
            b(0) = ValidFC[i]->first->p.x() - ValidFC[j]->first->p.x();
            b(1) = ValidFC[i]->first->p.y() - ValidFC[j]->first->p.y();
            Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
            if(0 < ts(0) && ts(0) < 1){_f += 1.0/ts(0); w++;}
        }
        f += _f;
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

Eigen::Vector3d RevisionVertices::decideRulingDirectionOn3d(Eigen::Vector3d e, Eigen::Vector3d N, double a, double phi){
    e = e.normalized();
    N = (Eigen::AngleAxisd(a, -e) * N).normalized(); // 回転行列を作成
    return (Eigen::AngleAxisd(phi, N) * e).normalized();
}

double Fmin_TriangleArea(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, bool IsSameWeight){
    double f = 0.0;
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    //triangle area
    for(int i = 1; i < (int)ValidFC.size()-2; i++)
    if(((ValidFC[i]->second->p - ValidFC[i]->first->p).normalized()).dot((ValidFC[i+1]->second->p - ValidFC[i+1]->first->p).normalized()) > 1.0 - 1e-13){
        f += 0.0;
    }else{
        Eigen::Vector3d e = (ValidFC[i+1]->first->p - ValidFC[i]->first->p), r = (ValidFC[i]->second->p - ValidFC[i]->first->p).normalized(), r2 = (ValidFC[i+1]->second->p - ValidFC[i+1]->first->p).normalized();
        double phi1 = std::acos(r.dot(e.normalized())), phi2 = std::acos(r2.dot(-e.normalized()));
        double l = sin(phi2)/sin(phi1+phi2)*(ValidFC[i+1]->first->p - ValidFC[i]->first->p).norm();
        double s = 2.0/(e.cross(l*r)).norm();
        if(IsSameWeight)f += s;
        else f += (std::sin(std::numbers::pi - (ValidFC.size()/2  - i)*std::numbers::pi) + 0.1)*s;
    }
    if(!IsSameWeight){
        Eigen::Vector3d e = (ValidFC[2]->first->p - ValidFC[1]->first->p), r = (ValidFC[1]->second->p - ValidFC[1]->first->p).normalized(), r2 = (ValidFC[2]->second->p - ValidFC[2]->first->p).normalized();
        double phi1 = std::acos(r.dot(e.normalized())), phi2 = std::acos(r2.dot(-e.normalized()));
        double l = sin(phi2)/sin(phi1+phi2)*(ValidFC[2]->first->p - ValidFC[1]->first->p).norm();
        double s = 2.0/(e.cross(l*r)).norm();
        f += (std::sin(std::numbers::pi - ValidFC.size()/2*std::numbers::pi) + 0.1)*s;

        e = (ValidFC.end()[-2]->first->p - ValidFC.end()[-3]->first->p), r = (ValidFC.end()[-3]->second->p - ValidFC[1]->first->p).normalized(), r2 = (ValidFC.end()[-2]->second->p - ValidFC.end()[-2]->first->p).normalized();
        phi1 = std::acos(r.dot(e.normalized())), phi2 = std::acos(r2.dot(-e.normalized()));
        l = sin(phi2)/sin(phi1+phi2)*(ValidFC.end()[-2]->first->p - ValidFC.end()[-3]->first->p).norm();
        s = 2.0/(e.cross(l*r)).norm();
        f += (std::sin(std::numbers::pi - ValidFC.size()/2*std::numbers::pi) + 0.1)*s;
    }
    return f;
}

//Xの先頭はflap angle
double Fmin_simularity(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, std::vector<std::shared_ptr<Vertex4d>>& BasePts){
    double f = 0.0;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){
        f += (FoldingCurve[i]->first->p - BasePts[i]->first->p).norm();
    }
    return f;
}

double Fconst_conv(const std::vector<double> &X, std::vector<double> &grad, void* f_data){
    auto k_conv = [](std::shared_ptr<CrvPt_FL>& v, std::shared_ptr<Vertex4d>& v2, std::shared_ptr<CrvPt_FL>& v3)->double{
        Eigen::Vector3d e = (v->p - v2->first->p).normalized(), e2 = (v3->p - v2->first->p).normalized(), axis = (v2->third->p - v2->first->p).normalized();
        return 2*std::numbers::pi - std::acos(e.dot(axis)) - std::acos(e2.dot(axis));
    };

    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    for(int i = 1; i < static_cast<int>(X.size()); i++){
        od->FC[i-1]->first->p =  X[i]*(od->BasePt[i-1]->first->p - od->BasePt[i-1]->third->p) + od->FC[i-1]->third->p;
    }
    double f = 0.0;
    std::vector<std::shared_ptr<Vertex>> Poly_V = od->Poly_V;
     _FoldingAAAMethod(od->FC, Poly_V, X[0]);
    std::vector<double> K_base;
    for(int i = 1; i < (int)od->FC.size() -1;i++){
        double k = k_conv(od->BasePt[i-1]->first, od->BasePt[i], od->BasePt[i+1]->first), k2 = k_conv(od->FC[i-1]->first, od->FC[i], od->FC[i+1]->first);
        f += -(std::numbers::pi - k)*(std::numbers::pi - k2);
        K_base.push_back(k);
    }

    std::vector<double> a = X;
    if(!grad.empty()){      
        for(int i = 0; i < (int)X.size(); i++){
            double fp = 0.0, fm  = 0.0;
            a[i] = X[i] + eps;
            if(i != 0){
                od->FC[i-1]->first->p =  a[i]*(od->BasePt[i-1]->first->p - od->BasePt[i-1]->third->p) + od->FC[i-1]->third->p;
            }
            for(int i = 1; i < (int)od->FC.size() -1;i++)
                fp += -(std::numbers::pi - K_base[i])*(std::numbers::pi - k_conv(od->FC[i-1]->first, od->FC[i], od->FC[i+1]->first));

            a[i] = X[i] - eps;
            if(i != 0){
                od->FC[i-1]->first->p =  a[i]*(od->BasePt[i-1]->first->p - od->BasePt[i-1]->third->p) + od->FC[i-1]->third->p;
            }
            for(int i = 1; i < (int)od->FC.size() -1;i++)
                fm += -(std::numbers::pi - K_base[i])*(std::numbers::pi - k_conv(od->FC[i-1]->first, od->FC[i], od->FC[i+1]->first));

            if(i != 0)od->FC[i-1]->first->p =  X[i]*(od->BasePt[i-1]->first->p - od->BasePt[i-1]->third->p) + od->FC[i-1]->third->p;
            grad[i] = (fp - fm)/(2.0 * eps);
            if(i != 0)grad[i] *= w_at;
            a[i] = X[i];
        }
    }
    if(DebugMode::Singleton::getInstance().isdebug())qDebug() <<"constraint function = " << MathTool::rad2deg(X[0]) << "(" << X[0] << ")  , " << f ;
    return f;
}

bool FoldLine::isbend(){
    double fr = RulingsCrossed(FoldingCurve);
    return (fr < 1e-9 && a_flap != -1)? true: false;
}

double Fruling(const std::vector<double> &X, std::vector<double> &grad, void* f_data)
{
    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    std::vector<std::shared_ptr<Vertex>> Poly_V = od->Poly_V;
    _FoldingAAAMethod(od->FC, Poly_V, X[0]);
    double f = RulingsCrossed(od->FC);
    std::vector<double> a = X;
    if(!grad.empty()){
        double fp, fm;
        for(int i = 0; i < (int)X.size(); i++){
            a[i] = X[i] + eps;
            _FoldingAAAMethod(od->FC, Poly_V, a[0]); fp = RulingsCrossed(od->FC);
            a[i] = X[i] - eps;
            _FoldingAAAMethod(od->FC, Poly_V, a[0]);  fm = RulingsCrossed(od->FC);

            grad[i] = (fp - fm)/(2.0 * eps);
        }
    }

    if(DebugMode::Singleton::getInstance().isdebug())qDebug() <<"constraint function = " << MathTool::rad2deg(a[0]) << "(" << a[0] << ")  , " << f ;
    return f;
 }


double Fruling4Vertex(const std::vector<double> &X, std::vector<double> &grad, void* f_data)
{
    auto f_intersection = [&X](const std::vector<std::shared_ptr<Vertex4d>>& _FC){
        auto calcCrossPt = [](const std::shared_ptr<Vertex4d>& x, const std::shared_ptr<Vertex4d>& x2)->double{
            Eigen::Vector3d e = (x2->first->p - x->first->p), r = (x->second->p - x->first->p).normalized(), r2 = (x2->second->p - x2->first->p).normalized();
            double phi1 = std::acos(r.dot(e.normalized())), phi2 = std::acos(r2.dot(-e.normalized()));
            return  sin(phi2)/sin(phi1+phi2)*(x2->first->p - x->first->p).norm();
        };

        if((int)_FC.size() <= 3)return 0.0;
        std::vector<std::shared_ptr<Vertex4d>> FC;
        double f = 0.0;
        for(auto&fc: _FC){if(fc->IsCalc)FC.push_back(fc);}
        std::vector<Eigen::Vector3d> RegCrv = RevisionVertices::getRegCrvPt(FC);
        double sum = 0.0;
        Eigen::Vector3d center = getCenter(RegCrv, sum);

        for(int i = 0; i< (int)FC.size(); i++){
            if(i == 0){
                double lmax = (FC[1]->second->p - FC[1]->first->p).norm(), l = calcCrossPt(FC[1],FC[2]);
                if(0 < l && l < lmax){
                    double w = std::sin(std::numbers::pi - (FC.size()/2  - i)*std::numbers::pi) + 0.1;
                    f += (RegCrv[i]-center).norm()/sum * l/lmax;
                }
            }else if(i == (int)FC.size() - 1){
                double lmax = (FC[i-3]->second->p - FC[i-3]->first->p).norm(), l = calcCrossPt(FC[i-3],FC[i-2]);
                if(0 < l && l < lmax){
                    double w = std::sin(std::numbers::pi - (FC.size()/2  - i)*std::numbers::pi) + 0.1;
                    f += (RegCrv[i]-center).norm()/sum * l/lmax;
                }
            }else{
                double lmax = (FC[i]->second->p - FC[i]->first->p).norm(), l = calcCrossPt(FC[i],FC[i+1]);
                if(0 < l && l < lmax){
                    double w = std::sin(std::numbers::pi - (FC.size()/2  - i)*std::numbers::pi) + 0.1;
                    f += (RegCrv[i]-center).norm()/sum * l/lmax;
                }
            }
        }
        return f;
    };
    RevisionVertices::ObjData_v *od = (RevisionVertices::ObjData_v *)f_data;
    std::vector<std::shared_ptr<Vertex>> Poly_V = od->Poly_V;
    for(int i = 0; i < static_cast<int>(X.size()); i++){
        od->FC[i]->first->p = X[i] * (od->FC[i]->line_parent->v->p - od->FC[i]->first->p).normalized() + od->FC[i]->first->p;
        od->FC[i]->first->p3 = X[i] * (od->FC[i]->line_parent->v->p3 - od->FC[i]->first->p3).normalized() + od->FC[i]->first->p3;
    }
    _FoldingAAAMethod(od->FC, Poly_V, od->a);
    double f = RulingsCrossed(od->FC);
    for(int i = 0; i < (int)od->FC.size(); i++){
        if(i == 0)f += f_intersection(od->FC);
    }
    std::vector<double> a = X;
    if(!grad.empty()){
        double fp, fm;
        for(int i = 0; i < (int)X.size(); i++){
            a[i] = X[i] + eps;
            if(i != 0){
                od->FC[i]->first->p = a[i] * (od->FC[i]->line_parent->v->p - od->FC[i]->first->p).normalized() + od->FC[i]->first->p;
                od->FC[i]->first->p3 = a[i] * (od->FC[i]->line_parent->v->p3 - od->FC[i]->first->p3).normalized() + od->FC[i]->first->p3;
            }
             _FoldingAAAMethod(od->FC, Poly_V, od->a);fp = f_intersection(od->FC);
            a[i] = X[i] - eps;
            if(i != 0){
                od->FC[i]->first->p = a[i] * (od->FC[i]->line_parent->v->p - od->FC[i]->first->p).normalized() + od->FC[i]->first->p;
                od->FC[i]->first->p3 = a[i] * (od->FC[i]->line_parent->v->p3 - od->FC[i]->first->p3).normalized() + od->FC[i]->first->p3;
            }
            _FoldingAAAMethod(od->FC, Poly_V, od->a);  fm = f_intersection(od->FC);
            if(i != 0){
                od->FC[i]->first->p = X[i] * (od->FC[i]->line_parent->v->p - od->FC[i]->first->p).normalized() + od->FC[i]->first->p;
                od->FC[i]->first->p3 = X[i] * (od->FC[i]->line_parent->v->p3 - od->FC[i]->first->p3).normalized() + od->FC[i]->first->p3;
            }
            grad[i] = (fp - fm)/(2.0 * eps);
            a[i] = X[i];
        }
    }

    qDebug() <<"constraint function = "  << f ;
    return f;
}

double ObjFunc_RulingIntersection(const std::vector<double> &X, std::vector<double> &grad, void* f_data){
    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    std::vector<std::shared_ptr<Vertex>> Poly_V = od->Poly_V;
    for(int i = 1; i < static_cast<int>(X.size()); i++){
        od->FC[i-1]->first->p =  X[i]*(od->BasePt[i-1]->first->p - od->BasePt[i-1]->third->p) + od->FC[i-1]->third->p;
        od->FC[i-1]->first->p3 =  X[i]*(od->BasePt[i-1]->first->p3 - od->BasePt[i-1]->third->p3) + od->FC[i-1]->third->p3;
    }
    _FoldingAAAMethod(od->FC, Poly_V, X[0]);
    double fb =  Fbend2(od->FC), fp = Fparallel(od->FC),fr = (od->wr != -1)?RulingsCrossed_NoOutline(od->FC): 0;

    std::vector<double> a = X;
    if(!grad.empty()){
       for(int i = 0; i < (int)X.size(); i++){
            a[i] = X[i] + eps;
            if(i != 0){
                od->FC[i-1]->first->p =  a[i]*(od->BasePt[i-1]->first->p - od->BasePt[i-1]->third->p) + od->FC[i-1]->third->p;
                od->FC[i-1]->first->p3 =  a[i]*(od->BasePt[i-1]->first->p3 - od->BasePt[i-1]->third->p3) + od->FC[i-1]->third->p3;
            }
            _FoldingAAAMethod(od->FC, Poly_V, a[0]);
            double fp = Fbend2(od->FC), fp2 = Fparallel(od->FC), fp3 = (od->wr != -1)?RulingsCrossed_NoOutline(od->FC): 0;

            a[i] = X[i] - eps;
            if(i != 0){
                od->FC[i-1]->first->p =  a[i]*(od->BasePt[i-1]->first->p - od->BasePt[i-1]->third->p) + od->FC[i-1]->third->p;
                od->FC[i-1]->first->p3 =  a[i]*(od->BasePt[i-1]->first->p3 - od->BasePt[i-1]->third->p3) + od->FC[i-1]->third->p3;
            }
            _FoldingAAAMethod(od->FC, Poly_V, a[0]);
            double fm = Fbend2(od->FC), fm2 = Fparallel(od->FC), fm3 = (od->wr != -1)?RulingsCrossed_NoOutline(od->FC): 0;
            grad[i] = od->wb * (fp - fm)/(2.0 * eps) + od->wp * (fp2 - fm2)/(2.0 * eps) + 100.0*(fp3 - fm3)/(2.0*eps);
            if(i != 0)grad[i] *= w_at;
            a[i] = X[i];
            if(i != 0){
                od->FC[i-1]->first->p =  X[i]*(od->BasePt[i-1]->first->p - od->BasePt[i-1]->third->p) + od->FC[i-1]->third->p;
                od->FC[i-1]->first->p3 =  X[i]*(od->BasePt[i-1]->first->p3 - od->BasePt[i-1]->third->p3) + od->FC[i-1]->third->p3;
            }
       }
    }
    return od->wb * fb + od->wp * fp + 100.0 * fr;
}

double ObjFunc_RegressionCurve(const std::vector<double> &a, std::vector<double> &grad, void* f_data){
    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    std::vector<std::shared_ptr<Vertex>>Poly_V = od->Poly_V;
    double fb = 0.0, fp = 0.0, fr = 0.0;
    _FoldingAAAMethod(od->FC, Poly_V, a[0]);
    fb =  Fbend2(od->FC); fp = Fparallel(od->FC); fr = Fmin_TriangleArea(od->FC, true);

    if(!grad.empty()){
       std::vector<double> X = a;
        for(int i = 0; i < (int)a.size(); i++){
        X[i] = a[i] + eps;

        _FoldingAAAMethod(od->FC, Poly_V, X[0]);
        double fp = 0*Fbend2(od->FC), fp2 = 0*Fparallel(od->FC), fp3 = Fmin_TriangleArea(od->FC, true);
        X[i] = a[i] - eps;
        _FoldingAAAMethod(od->FC, Poly_V, X[0]);
        double fm = 0*Fbend2(od->FC), fm2 = 0*Fparallel(od->FC), fm3 = Fmin_TriangleArea(od->FC, true);
        double dfb = od->wb * (fp - fm)/(2.0 * eps), dfp = 0 * (fp2 - fm2)/(2.0 * eps), dfa = 1.0*(fp3 - fm3)/(2.0*eps);
        grad[i] = dfb + dfp + dfa;
       }

    }
    if(DebugMode::Singleton::getInstance().isdebug())qDebug() <<"Triangle Area(" << MathTool::rad2deg(a[0]) << "(" << a[0] << ")  =  " <<  fb <<  " , grad = " << grad[0] << ", Area = "<< fr;
    return od->wb * fb + 0 * fp + 1.0 * fr;
}

double ObjFunc_Vertex(const std::vector<double> &X, std::vector<double> &grad, void* f_data){
    const double wr = 1e-3, wa = 1e-1;
    RevisionVertices::ObjData_v *od = (RevisionVertices::ObjData_v *)f_data;
    std::vector<std::shared_ptr<Vertex>> Poly_V = od->Poly_V;
    for(int i = 0; i < static_cast<int>(X.size()); i++){
       od->FC[i]->first->p = X[i] * (od->FC[i]->line_parent->v->p - od->FC[i]->first->p).normalized() + od->FC[i]->first->p;
       od->FC[i]->first->p3 = X[i] * (od->FC[i]->line_parent->v->p3 - od->FC[i]->first->p3).normalized() + od->FC[i]->first->p3;
    }
    _FoldingAAAMethod(od->FC, Poly_V, od->a);
    std::vector<Eigen::Vector3d> RegCrv = RevisionVertices::getRegCrvPt(od->FC);
    double sum = 0.0;
    Eigen::Vector3d center = getCenter(RegCrv, sum);
    double freg = 0.0, fsim = 0.0, ftri = 0.0;
    for(auto&x: RegCrv)freg += wr*(center - x).norm();
    for(int i = 0; i < (int)X.size(); i++){fsim += (RegCrv[i]-center).norm()/sum * X[i];}
    ftri = wa*Fmin_TriangleArea(od->FC, false);
    std::vector<double> a = X;
    if(!grad.empty()){
       for(int i = 0; i < (int)X.size(); i++){
        a[i] = X[i] + eps;
        if(i != 0){
            od->FC[i]->first->p = a[i] * (od->FC[i]->line_parent->v->p - od->FC[i]->first->p).normalized() + od->FC[i]->first->p;
            od->FC[i]->first->p3 = a[i] * (od->FC[i]->line_parent->v->p3 - od->FC[i]->first->p3).normalized() + od->FC[i]->first->p3;
        }
        _FoldingAAAMethod(od->FC, Poly_V, od->a);
        RegCrv = RevisionVertices::getRegCrvPt(od->FC);
        center = getCenter(RegCrv, sum);
        double frp = 0.0, fsimp = 0.0;
        //for(auto&x: RegCrv)frp += wr*(center - x).norm();
        for(int i = 0; i < (int)X.size(); i++){fsimp += (RegCrv[i]-center).norm()/sum * X[i];}
        double ftrip = wa*Fmin_TriangleArea(od->FC, false);
        a[i] = X[i] - eps;
        if(i != 0){
            od->FC[i]->first->p = a[i] * (od->FC[i]->line_parent->v->p - od->FC[i]->first->p).normalized() + od->FC[i]->first->p;
            od->FC[i]->first->p3 = a[i] * (od->FC[i]->line_parent->v->p3 - od->FC[i]->first->p3).normalized() + od->FC[i]->first->p3;
        }
        _FoldingAAAMethod(od->FC, Poly_V, a[0]);
        RegCrv = RevisionVertices::getRegCrvPt(od->FC);
        center = getCenter(RegCrv, sum);
        double frm = 0.0, fsimm = 0.0;
        //for(auto&x: RegCrv)frm += wr*(center - x).norm();
        for(int i = 0; i < (int)X.size(); i++){fsimm += (RegCrv[i]-center).norm()/sum * X[i];}
        double ftrim = wa*Fmin_TriangleArea(od->FC, false);
        grad[i] = ((frp + fsimp + ftrip) - (frm + fsimm + ftrim))/(2.0 * eps);
        a[i] = X[i];
        if(i != 0){
            od->FC[i]->first->p = X[i] * (od->FC[i]->line_parent->v->p - od->FC[i]->first->p).normalized() + od->FC[i]->first->p;
            od->FC[i]->first->p3 = X[i] * (od->FC[i]->line_parent->v->p3 - od->FC[i]->first->p3).normalized() + od->FC[i]->first->p3;
        }
       }
    }
    qDebug() <<"simlarity = " <<fsim << " ,  density of intersection = " << freg << " , sum of triangle area = " << ftri;
    return freg + fsim;
}

bool FoldLine::Optimization_Vertex(const std::vector<std::shared_ptr<Vertex>>& Poly_V){
    nlopt::opt opt;
    std::vector<double> X, bnd_lower, bnd_upper;
    std::vector<std::shared_ptr<Vertex4d>> FC;
    for(auto&fc: FoldingCurve){
        std::shared_ptr<Vertex4d> v4d = fc->deepCopy(); v4d->addline(fc->line_parent);
        FC.push_back(v4d);
        X.push_back(0);
        bnd_lower.push_back(-(fc->first->p - fc->line_parent->o->p).norm());
        bnd_upper.push_back((fc->first->p - fc->line_parent->v->p).norm());
    }
    RevisionVertices::ObjData_v od = {a_flap, FC, Poly_V};
    opt = nlopt::opt(nlopt::LD_MMA, X.size());
    opt.set_min_objective(ObjFunc_Vertex, &od);
    opt.add_inequality_constraint(Fruling4Vertex, &od);
    opt.set_lower_bounds(bnd_lower);
    opt.set_upper_bounds(bnd_upper);
    //opt.set_maxtime(10.0);//stop over this time
    //opt.set_xtol_rel(1e-9);
    double minf;
    try {
        nlopt::result result = opt.optimize(X, minf);
        for(int i = 0; i < (int)X.size(); i++){
            FC[i]->first->p = X[i] * (FoldingCurve[i]->line_parent->v->p - FoldingCurve[i]->first->p).normalized() + FoldingCurve[i]->first->p;
            FC[i]->first->p3 = X[i] * (FoldingCurve[i]->line_parent->v->p3 - FoldingCurve[i]->first->p3).normalized() + FoldingCurve[i]->first->p3;
        }
        _FoldingAAAMethod(FC, Poly_V, a_flap);
        double fc = RulingsCrossed(FC);
        qDebug() <<"result = " << fc<< " , error " << minf <<result;
        for(int i = 0; i< (int)FC.size(); i++){FoldingCurve[i]->first->p = FC[i]->first->p; FoldingCurve[i]->first->p3 = FC[i]->first->p3;}

    }catch (std::exception& e){qDebug() << "nlopt failed: " << e.what();}
    return false;
}

bool FoldLine::Optimization_FlapAngle(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double wb, double wp, int rank, int alg, bool ConstFunc){
    if(FoldingCurve.empty())return false;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    double a = 0;
    Eigen::Vector3d e = (ValidFC[0]->first->p3 - ValidFC[1]->first->p3).normalized(),e2 = (ValidFC[2]->first->p3 - ValidFC[1]->first->p3).normalized();
    Eigen::Vector3d SpinAxis = (ValidFC[1]->third->p3 - ValidFC[1]->first->p3).normalized();
    Eigen::Vector3d Nb = -(e.cross(e2)).normalized(), N4 = (e.cross(SpinAxis)).normalized();

    double a_ll = std::atan2(e.dot(Nb.cross(N4)), N4.dot(-Nb));
    if(a_ll < 0)a_ll += 2.0*std::numbers::pi;
    double a_con = a_ll + std::numbers::pi;
    if(a_con > 2.0*std::numbers::pi)a_con -= 2.0*std::numbers::pi;
    if(a_con < 0)a_con +=2.0*std::numbers::pi;
    qDebug() << "a_con = " << MathTool::rad2deg(a_con);
    double phi3 = std::acos(e2.dot(SpinAxis)), phi4 = std::acos(e.dot(SpinAxis));
    double k = 2.0 * std::numbers::pi - phi3 - phi4;
    double a_min, a_max;
    std::vector<double>bnd_lower,bnd_upper;
    Eigen::Vector3d FaceNp = (SpinAxis.cross(e2)).normalized(), CrossV = (N4.cross(FaceNp)).normalized();
    bool IsMount = (SpinAxis.dot(CrossV) < 0)? true: false;

    if(k < std::numbers::pi && IsMount){a_min = a_con + std::numbers::pi; a_max = 2.0 * std::numbers::pi;}
    if(k >= std::numbers::pi && IsMount){a_min = std::numbers::pi; a_max = std::numbers::pi + a_con;}
    if(k < std::numbers::pi && !IsMount){a_min = 0.0;a_max = a_con - std::numbers::pi;}
    if(k >= std::numbers::pi && !IsMount){a_min = a_con - std::numbers::pi; a_max = std::numbers::pi;}

    if(DebugMode::Singleton::getInstance().isdebug()) qDebug() <<MathTool::rad2deg(a_ll) << " , " << MathTool::rad2deg(a_con);
    if(DebugMode::Singleton::getInstance().isdebug()){
        std::ofstream ofs2;
        std::filesystem::create_directory("./Optimization");
        std::string AngleFile;
        AngleFile = "./Optimization/OptimArea_" + std::to_string(rank) + "_" + std::to_string(validsize) + ".csv";
        ofs2.open(AngleFile, std::ios::out);
        ofs2 << "a(radian), a(degree) , Eruling, Ebend, Eparalell, Eruling(N ooutline), Eruling(all) , Earea\n";
        for(double _a = a_min; _a <= a_max; _a+= 1e-3){
            _FoldingAAAMethod(FoldingCurve, Poly_V, _a);
            double f = RulingsCrossed(FoldingCurve);
            double fb = Fbend2(FoldingCurve);
            double fp =  Fparallel(FoldingCurve);
            double fr = RulingsCrossed_NoOutline(FoldingCurve);
            double fr2 = RulingCrossed_all(FoldingCurve);
            double fa = Fmin_TriangleArea(FoldingCurve, false);
            ofs2 << _a << ", " << MathTool::rad2deg(_a) << " , " << f << ", " << fb << ", " <<fp <<", "<< fr << ", " << fr2 << ", " << fa << std::endl ;
        }ofs2.close();
    }
    //std::vector<std::shared_ptr<Vertex4d>> FC_amin, FC_amax;
    //for(auto&fc: FoldingCurve){FC_amin.push_back(fc->deepCopy()); FC_amax.push_back(fc->deepCopy());}
    RevisionVertices::ObjData od_amin = {FoldingCurve, Poly_V}, od_amax = {FoldingCurve, Poly_V};
    od_amin.AddWeight(wb, 0.1*wp, -1); od_amax.AddWeight(wb, 0.1*wp, -1);
    //if(!ConstFunc)od.AddWeight(wb, wp, -1); else od.AddWeight(wb, wp, 100.0);
    //od.AddWeight(wb, wp, 100.0);
    nlopt::opt opt_amin, opt_amax;
    std::vector<double> X_amin{a_min + 0.5}, X_amax{a_max - 0.5};
    bnd_lower.push_back(a_min); bnd_upper.push_back(a_max);
    if(alg == 0){//ruling intersection as a constraint
        opt_amin = nlopt::opt(nlopt::LD_MMA, X_amin.size()); opt_amax = nlopt::opt(nlopt::LD_MMA, X_amax.size());
        opt_amin.set_min_objective(ObjFunc_RulingIntersection, &od_amin);
        opt_amax.set_min_objective(ObjFunc_RulingIntersection, &od_amax);
    }else if(alg == 1){//maximize triangle area
        od_amin.AddBasePt(FoldingCurve);
        opt_amin = nlopt::opt(nlopt::LD_MMA, 1); opt_amax = nlopt::opt(nlopt::LD_MMA, 1);
        opt_amin.set_min_objective(ObjFunc_RegressionCurve, &od_amin);
        opt_amin.add_inequality_constraint(Fruling, &od_amin);

        od_amax.AddBasePt(FoldingCurve);
        opt_amax.set_min_objective(ObjFunc_RegressionCurve, &od_amax);
        opt_amax.add_inequality_constraint(Fruling, &od_amax);
    }

    //qDebug()<<"now optimization of flap angle is not applied constraint function of ruling intersection";
    //if(ConstFunc){opt_amin.add_inequality_constraint(Fruling, &od_amin); opt_amax.add_inequality_constraint(Fruling, &od_amax);}
    //opt_amin.add_inequality_constraint(Fruling, &od_amin); opt_amax.add_inequality_constraint(Fruling, &od_amax);
    opt_amin.set_lower_bounds(bnd_lower);   opt_amax.set_lower_bounds(bnd_lower);
    opt_amin.set_upper_bounds(bnd_upper);   opt_amax.set_upper_bounds(bnd_upper);

    //opt.set_param("inner_maxeval", 100);
    opt_amin.set_maxtime(5.0);//stop over this time
    opt_amin.set_xtol_rel(1e-7);

    opt_amax.set_maxtime(5.0);//stop over this time
    opt_amax.set_xtol_rel(1e-7);

    qDebug() << "area " << MathTool::rad2deg(a_min) << " < " << MathTool::rad2deg(a) << " < " << MathTool::rad2deg(a_max) ;

    double minf_amin, minf_amax, f_ruling_amin, f_ruling_amax;
    try {
        qDebug() <<"small a" ;
        nlopt::result result = opt_amin.optimize(X_amin, minf_amin);
        _FoldingAAAMethod(FoldingCurve, Poly_V, X_amin[0]);
        qDebug() <<"result :  lower bound" <<result ;
        f_ruling_amin = RulingsCrossed(FoldingCurve);
        qDebug() << "found minimum at f(" << MathTool::rad2deg(X_amin[0]) << ") = " << minf_amin << "  , fruling = " << f_ruling_amin << ", ruling num = " << validsize << "\n";
    }catch (std::exception& e) { qDebug() << "nlopt failed: " << e.what() ; }

    try {
        qDebug() << "large a" ;
        nlopt::result result = opt_amax.optimize(X_amax, minf_amax);
        _FoldingAAAMethod(FoldingCurve, Poly_V, X_amax[0]);
        f_ruling_amax = RulingsCrossed(FoldingCurve);
        qDebug() <<"result :  upper bound" <<result;
        qDebug() << "found minimum at f(" << MathTool::rad2deg(X_amax[0]) << ") = " << minf_amax << "  ,  fruling = "  << f_ruling_amax << " , ruling num = " << validsize << "\n";
    }catch (std::exception& e){qDebug() << "nlopt failed: " << e.what();}

    auto movePt = [](std::vector<std::shared_ptr<Vertex4d>>&dest, std::vector<std::shared_ptr<Vertex4d>>& src, bool IsPrintError){
        double er = 0.0;
        for(int i = 0; i < (int)src.size(); i++){
            er += (dest[i]->first->p - src[i]->first->p).norm();
            if(IsPrintError)
                qDebug() << i <<  "  :  distance of point 2d = " << dest[i]->first->p.x() << " , " << dest[i]->first->p.y() << ", " << dest[i]->first->p.z() << " , new point = " <<
                            src[i]->first->p.x() << " , " << src[i]->first->p.y() << ", " << src[i]->first->p.z();
            dest[i]->first->p2_ori = dest[i]->first->p = src[i]->first->p; dest[i]->first->p3_ori = dest[i]->first->p3 = src[i]->first->p3;
        }
        if(IsPrintError)qDebug()<<"error of distance = " << er;
    };
    bool IsPrintError = false;

    double th_ruling = 1e-9;
    if(f_ruling_amin >= th_ruling && f_ruling_amax > th_ruling){

        if(minf_amin < minf_amax){
            //movePt(FoldingCurve, FC_amin, IsPrintError);
            a_flap = X_amin[0];
        }else{
            //movePt(FoldingCurve, FC_amax, IsPrintError);
            a_flap = X_amax[0];
        }
        _FoldingAAAMethod(FoldingCurve, Poly_V, a_flap);
        qDebug()<<"could not avoid ruling cross";
        qDebug() << "result : smaller = " << MathTool::rad2deg(a_flap)  << "(" << a_flap << ") " <<  "  f_bend = " << Fbend2(FoldingCurve) << "  ,  f_ruling = " <<RulingsCrossed(FoldingCurve) ;
        return false;
    }

    if(f_ruling_amin < th_ruling && f_ruling_amax < th_ruling){
        if(minf_amin < minf_amax){
            //movePt(FoldingCurve, FC_amin, IsPrintError);
            a_flap = X_amin[0];
        }else{
            //movePt(FoldingCurve, FC_amax, IsPrintError);
            a_flap = X_amax[0];
        }
    }else if(f_ruling_amin < th_ruling && f_ruling_amax >= th_ruling){
       // movePt(FoldingCurve, FC_amin, IsPrintError);
        a_flap = X_amin[0];
    }else if(f_ruling_amin >= th_ruling && f_ruling_amax < th_ruling){
        //movePt(FoldingCurve, FC_amax, IsPrintError);
        a_flap = X_amax[0];
    }
    _FoldingAAAMethod(FoldingCurve, Poly_V, a_flap);
    qDebug() << "result : smaller = " << MathTool::rad2deg(a_flap)  << "(" << a_flap << "), tol = " << tol <<  "  f_bend = " << Fbend2(FoldingCurve) << "  ,  f_ruling = " <<RulingsCrossed(FoldingCurve) ;
    qDebug() << "finish";
    return true;
}

void FoldLine::revisecrossedruling(const std::vector<std::shared_ptr<Vertex>>& Poly_v){
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    Eigen::VectorXd E_crossdst;
    while(1){
        E_crossdst = Eigen::VectorXd(FoldingCurve.size() - 2);
        for(int i = 1; i < (int)ValidFC.size() -1; i++){
            for(int j = -1; j <= 1; j +=2){
                int i2 = i + j;
                if(i2 <= 0 || i2 >= (int)ValidFC.size() - 1)continue;
                A(0,0) = ValidFC[i]->first->p.x() - ValidFC[i]->second->p.x();
                A(0,1) = -(ValidFC[i2]->first->p.x() - ValidFC[i2]->second->p.x());
                A(1,0) = ValidFC[i]->first->p.y() - ValidFC[i]->second->p.y();
                A(1,1) = -(ValidFC[i2]->first->p.y() - ValidFC[i2]->second->p.y());
                b(0) = ValidFC[i]->first->p.x() - ValidFC[i2]->first->p.x();
                b(1) = ValidFC[i]->first->p.y() - ValidFC[i2]->first->p.y();
                Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
                if(0 < ts(0) && ts(0) < 1){
                    E_crossdst(i-1) += 1.0/ts(0);
                    ValidFC[i]->IsCalc = false;
                }
            }
        }
        int cnt = 0;
        for(auto&fc: FoldingCurve)if(fc->IsCalc)cnt++;
        qDebug() << "valid size = " << cnt;
        SimpleSmooothSrf(Poly_v);
        return;
    }
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
    /*
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
    //validsize = FoldingCurve.size();
    */
    return true;
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
    validsize = iselim;
    validsize = (validsize < 0)? 0: (validsize > (int)FoldingCurve.size())? (int)FoldingCurve.size(): validsize;
    for(auto it = FoldingCurve.begin() + n; it != FoldingCurve.end() - n; it++){
        (*it)->IsCalc = true;
        (*it)->first->p = (*it)->first->p2_ori; (*it)->first->p3 = (*it)->first->p3_ori;
        (*it)->second->p = (*it)->second->p2_ori; (*it)->second->p3 = (*it)->second->p3_ori;
        (*it)->third->p = (*it)->third->p2_ori; (*it)->third->p3 = (*it)->third->p3_ori;
    }
    auto res = elim_rulings();
    for(auto&V4d: FoldingCurve){
        if(std::find_if(res.begin(), res.end(), [&V4d](const std::shared_ptr<Vertex4d>& V){return V4d->first->p == V->first->p;}) != res.end()){
        }else V4d->IsCalc = false;
    }
    return;
}

/*
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
*/
bool FoldLine::RevisionCrosPtsPosition(){
    if(FoldingCurve.empty() || (int)FoldingCurve.size() <= 4)return false;
    FoldingCurve[1]->IsCalc = false;
    FoldingCurve.end()[-2]->IsCalc = false;
    validsize -= 2;
    return true;
    FoldingCurve.front()->first->p3 = FoldingCurve.front()->first->p3_ori;
    FoldingCurve.back()->first->p3 = FoldingCurve.back()->first->p3_ori;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve)if(fc->IsCalc)ValidFC.push_back(fc);
    if(FoldingCurve.size() > 5){
        qDebug() << "before revision";
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
        double k0,k1,k2;
        //double k0 = RevisionVertices::getK(FoldingCurve[1]->first->p, FoldingCurve[0]->first->p, FoldingCurve[2]->first->p);
        //double k1 = RevisionVertices::getK(FoldingCurve[2]->first->p, FoldingCurve[1]->first->p, FoldingCurve[3]->first->p);
        //double k2 = RevisionVertices::getK(FoldingCurve[3]->first->p, FoldingCurve[2]->first->p, FoldingCurve[4]->first->p);
        //qDebug() <<"front " << abs(k0 - k1) << ", " << abs(k1 - k2) << " , " << 2*k1 - k2 ;
        if(abs(k0 - k1) > abs(k1 - k2)){
            qDebug() <<"front";
            FoldingCurve[0]->first->p = ReviseEndPosition(FoldingCurve[2]->first->p - FoldingCurve[1]->first->p,  -(2.*k1 - k2));
            Eigen::Vector3d newPos = getCrossPosition(FoldingCurve[0]->first, FoldingCurve[1]->first, FoldingCurve[0]->third, FoldingCurve[0]->second);
            FoldingCurve[0]->first->p = newPos;
            FoldingCurve[0]->first->p3 = (newPos - FoldingCurve[0]->third->p).norm()/(FoldingCurve[0]->third->p - FoldingCurve[0]->second->p).norm()
                    * (FoldingCurve[0]->second->p3 - FoldingCurve[0]->third->p3) + FoldingCurve[0]->third->p3;
        }
        //k0 = RevisionVertices::getK(FoldingCurve.end()[-2]->first->p, FoldingCurve.end()[-3]->first->p, FoldingCurve.back()->first->p);
        //k1 = RevisionVertices::getK(FoldingCurve.end()[-3]->first->p, FoldingCurve.end()[-4]->first->p, FoldingCurve.end()[-2]->first->p);
        //k2 = RevisionVertices::getK(FoldingCurve.end()[-4]->first->p, FoldingCurve.end()[-5]->first->p, FoldingCurve.end()[-3]->first->p);
        //qDebug() <<"back " << abs(k0 - k1) << ", " << abs(k1 - k2) << ", " << 2.*k1 - k2 <<  std::endl;
        if(abs(k0 - k1) > abs(k1 - k2)){
            qDebug() << "back";
            FoldingCurve.back()->first->p = ReviseEndPosition(FoldingCurve.end()[-3]->first->p - FoldingCurve.end()[-2]->first->p, 2.*k1 - k2);
            Eigen::Vector3d newPos2 = getCrossPosition(FoldingCurve.back()->first, FoldingCurve.end()[-2]->first, FoldingCurve.back()->third, FoldingCurve.back()->second);
            FoldingCurve.back()->first->p = newPos2;
            FoldingCurve.back()->first->p3 = (newPos2 - FoldingCurve.back()->third->p).norm()/
                    (FoldingCurve.back()->third->p - FoldingCurve.back()->second->p).norm() * (FoldingCurve.back()->second->p3 - FoldingCurve.back()->third->p3) + FoldingCurve.back()->third->p3;
        }

    }
    return true;

}

void FoldLine::ReassignColor(){
    //Y字のとき折り線は山であると仮定(谷の場合は色を反転)
    typedef struct {
        std::shared_ptr<Vertex4d> v;
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
        auto e = (ValidFC[i-1]->first->p - ValidFC[i]->first->p).normalized(), e2 = (ValidFC[i+1]->first->p - ValidFC[i]->first->p).normalized(), axis = (ValidFC[i]->third->p - ValidFC[i]->first->p).normalized();
        double k = 2.0*std::numbers::pi - std::acos(e.dot(axis)) - std::acos(e2.dot(axis));
        Eigen::Vector3d f_nv = ((ValidFC[i]->third->p3 - ValidFC[i]->first->p3).cross(ValidFC[i-1]->first->p3 - ValidFC[i]->first->p3)).normalized();
        Eigen::Vector3d fp_nv = ((ValidFC[i+1]->first->p3 - ValidFC[i]->first->p3).cross(ValidFC[i]->third->p3 - ValidFC[i]->first->p3)).normalized();
        Eigen::Vector3d SpinAxis = (ValidFC[i]->third->p3 - ValidFC[i]->first->p3).normalized();

        //std::string s = (SpinAxis.dot(f_nv.cross(fp_nv)) <-1e-5)? "red": "blue";
        //qDebug() << "color = " << s ;
        int mv = (abs(SpinAxis.dot(f_nv.cross(fp_nv))) < 1e-9 )? -1: (SpinAxis.dot(f_nv.cross(fp_nv)) < 1e-9)? 0: 1;
        int type_mvk = (mv == 0 && k < std::numbers::pi)? 0: (mv == 0 && k >= std::numbers::pi)? 1: (mv == 1 && k < std::numbers::pi)? 2: (mv == 1 && k >= std::numbers::pi)? 3: -1;
        InitState.push_back(MVK(ValidFC[i], type_mvk));
        MV[type_mvk + 1] += 1;
    }

    if(MV[0] != 0){qDebug() << "no folding pattern was found"; return;}
    if(MV[1] == MV[2] && MV[3] == MV[4] && MV[1] == MV[4]){
    }else{
        for(auto& IS: InitState){
            if(MV[1] + MV[4] > MV[2] + MV[3] ){
                if(IS.type_mvk == 1 || IS.type_mvk == 2){IS.v->IsCalc = false;}
            }else{
                if(IS.type_mvk == 0 || IS.type_mvk == 3){IS.v->IsCalc = false;}
            }
        }
    }
    qDebug() <<"color changed";
}

template class PointOnEndEdge<std::shared_ptr<Vertex>>;

void FoldLine::reassignruling(std::shared_ptr<FoldLine>& parent, const std::vector<std::shared_ptr<Line>>& Surface, const std::vector<std::shared_ptr<Line>>& Rulings){
    auto getCrossPoint = [](std::vector<Eigen::Vector3d>& CtrlPts, std::shared_ptr<Vertex> v, std::shared_ptr<Vertex> o, int dim){
        std::vector<double>arcT = BezierClipping(CtrlPts, v, o, dim);
        for(auto&t: arcT){
            if(t < 0 || 1 < t){qDebug()<<"t is not correct value " << t ; continue;}
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

    if(parent->FoldingCurve.empty() || FoldingCurve.empty())return;
    int dim = CtrlPts.size() - 1;


    std::vector<std::shared_ptr<Line>> polyline_surface;
    for(auto&l: Surface) polyline_surface.push_back(std::make_shared<Line>(l->o, l->v, l->et));

    for(auto it = parent->FoldingCurve.begin(); it != parent->FoldingCurve.end(); it++){
        if(!(*it)->IsCalc)continue;
        std::shared_ptr<Vertex> dowscastV = std::dynamic_pointer_cast<Vertex>((*it)->first);
        EdgeSplit(polyline_surface, dowscastV);
        EdgeSplit(polyline_surface, (*it)->second);
        EdgeSplit(polyline_surface, (*it)->third);
    }
    for(auto&r: Rulings){
        if(r->hasCrossPoint)continue;
        EdgeSplit(polyline_surface, r->v);
        EdgeSplit(polyline_surface, r->o);
    }


    //Surfaceの各Lineの要素に影響がないか調べること
    for(auto&l: polyline_surface){
        std::shared_ptr<CrvPt_FL> p = getCrossPoint(CtrlPts, l->o, l->v, dim);
        if(p == nullptr)continue;
        if(std::abs(p->s - FoldingCurve.front()->first->s) < std::abs(p->s - FoldingCurve.back()->first->s)){
            FoldingCurve.front()->third = l->o; FoldingCurve.front()->second = l->v; FoldingCurve.front()->first = p;
        }else{
            FoldingCurve.back()->third = l->o; FoldingCurve.back()->second = l->v; FoldingCurve.back()->first = p;
        }
    }

    //SplitOnEndLine(FoldingCurve.front(), parent->FoldingCurve.front(), parent->FoldingCurve, dim);
    //SplitOnEndLine(FoldingCurve.back(),  parent->FoldingCurve.back(), parent->FoldingCurve, dim);
    for(auto it = parent->FoldingCurve.begin() + 1; it != parent->FoldingCurve.end() - 1; it++){
        if(!(*it)->IsCalc)continue;
        std::shared_ptr<CrvPt_FL> p = getCrossPoint(CtrlPts, (*it)->first, (*it)->second, dim);
        if(p != nullptr){
            std::shared_ptr<Vertex4d> v4d = std::make_shared<Vertex4d>(p, (*it)->second, (*it)->first);
            v4d->addline((*it));
            FoldingCurve.push_back(v4d);
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

std::vector<std::vector<std::shared_ptr<Vertex>>> FoldLine::CalclateRegressionCurve(double a, const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool IsWriteCSV, std::vector<std::vector<std::shared_ptr<Vertex>>>& Tri_fixside){
    auto getCrossPoint = [](double l, std::shared_ptr<Vertex4d>& x, bool IsFixed){
        Eigen::Vector3d r2d = (!IsFixed)?(x->second->p - x->first->p).normalized(): (x->third->p - x->first->p).normalized();
        Eigen::Vector3d r3d = (!IsFixed)?(x->second->p3 - x->first->p3).normalized(): (x->third->p3 - x->first->p3).normalized();
        return std::make_shared<Vertex>(Vertex(l*r2d + x->first->p, l*r3d + x->first->p3));
    };

    if(FoldingCurve.size() < 4)return{};
    Tri_fixside.clear();
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(const auto&FC: FoldingCurve){if(FC->IsCalc)ValidFC.push_back(FC);}
    a = MathTool::deg2rad(a);
    std::vector<std::array<double,3>> ABKs;
    double a2;
    Eigen::Vector3d befedge, befN;
    for(int ind = 1; ind < (int)ValidFC.size() - 1; ind++){
        Eigen::Vector3d e = (ValidFC[ind-1]->first->p3 - ValidFC[ind]->first->p3).normalized(), e2 = (ValidFC[ind+1]->first->p3 - ValidFC[ind]->first->p3).normalized(), axis = (ValidFC[ind]->third->p3 - ValidFC[ind]->first->p3).normalized();
        if(ind != 1){
            Eigen::Vector3d edge = MathTool::ProjectionVector(e2, -e, true);
            Eigen::Vector3d SrfN = MathTool::ProjectionVector(e.cross(e2), -e, true);
            a = update_flapangle(a2, befN, SrfN, e);
        }
        CalcRuling(a, ValidFC[ind-1], ValidFC[ind], ValidFC[ind+1], Poly_V, a2, befN, -Eigen::Vector3d::UnitZ());
        double b = std::acos(e.dot(e2)), k = 2.0*std::numbers::pi - (std::acos(axis.dot(e)) + std::acos(axis.dot(e2)));
        ABKs.push_back({a,b,k});
        if(ind == (int)FoldingCurve.size() -2)break;
        befN = MathTool::ProjectionVector(e.cross(e2), e2, true);befedge = MathTool::ProjectionVector(e, e2, true);
    }

    std::vector<std::vector<std::shared_ptr<Vertex>>> Triangles;
    if(IsWriteCSV){

        Eigen::Vector3d e = (ValidFC[0]->first->p3 - ValidFC[1]->first->p3).normalized(),e2 = (ValidFC[2]->first->p3 - ValidFC[1]->first->p3).normalized();
        Eigen::Vector3d SpinAxis = (ValidFC[1]->third->p3 - ValidFC[1]->first->p3).normalized();
        Eigen::Vector3d Nb = -(e.cross(e2)).normalized(), N4 = (e.cross(SpinAxis)).normalized();

        double a_ll = std::atan2(e.dot(Nb.cross(N4)), N4.dot(-Nb));
        if(a_ll < 0)a_ll += 2.0*std::numbers::pi;
        double acon = a_ll + std::numbers::pi;
        if(acon > 2.0*std::numbers::pi)acon -= 2.0*std::numbers::pi;
        if(acon < 0)acon +=2.0*std::numbers::pi;

        std::ofstream ofs;
        std::filesystem::create_directory("./Optimization");
        std::string File = "./Optimization/RegressionCurve_";
        File += std::to_string(ValidFC.size() - 3) + ".csv" ;
        ofs.open(File, std::ios::out);
        ofs << "a(radian), a(degree), area, Ebend, Eruling, |a - a_con|"<<std::endl;
        for(double _a = 0.0; _a <= 2.0*std::numbers::pi; _a+= 1e-3){
            a = _a;
            for(int ind = 1; ind < (int)ValidFC.size() - 1; ind++){
                Eigen::Vector3d e = (ValidFC[ind-1]->first->p3 - ValidFC[ind]->first->p3).normalized(), e2 = (ValidFC[ind+1]->first->p3 - ValidFC[ind]->first->p3).normalized(), axis = (ValidFC[ind]->third->p3 - ValidFC[ind]->first->p3).normalized();
                if(ind != 1){
                    Eigen::Vector3d edge = MathTool::ProjectionVector(e2, -e, true);
                    a = update_flapangle2(a2, befedge, edge, e);
                }
                CalcRuling(a, ValidFC[ind-1], ValidFC[ind], ValidFC[ind+1], Poly_V, a2, befN, -Eigen::Vector3d::UnitZ());
                if(ind == (int)FoldingCurve.size() -2)break;
                befN = MathTool::ProjectionVector(e.cross(e2), e2, true);befedge = MathTool::ProjectionVector(e, e2, true);
            }
            double area = 0.0, Ebend = 0.0, Eruling = 0.0;
            for(int i = 1; i < (int)ValidFC.size()-2; i++){
                //triangle area
                if(((ValidFC[i]->second->p - ValidFC[i]->first->p).normalized()).dot((ValidFC[i+1]->second->p - ValidFC[i+1]->first->p).normalized()) > 1.0 - 1e-7){
                    area = 0.0; break;
                }else{
                    auto crosspoint = calcCrossPoint_2Vertex(ValidFC[i]->first, ValidFC[i]->second, ValidFC[i+1]->first, ValidFC[i+1]->second);   
                    Eigen::Vector3d e = ValidFC[i+1]->first->p - ValidFC[i]->first->p, e2 = crosspoint - ValidFC[i]->first->p;
                    double s = 2.0/(e.cross(e2)).norm();
                    area += s;
                }

                //bending energy
                Eigen::Vector3d N = (ValidFC[i-1]->first->p3 - ValidFC[i]->first->p3).cross(ValidFC[i]->second->p3 - ValidFC[i]->first->p3).normalized();
                Eigen::Vector3d N2 = (ValidFC[i]->second->p3 - ValidFC[i]->first->p3).cross(ValidFC[i+1]->first->p3 - ValidFC[i]->first->p3).normalized();
                Ebend += 1.0/std::pow(std::acos(N.dot(N2)),2);

                //Is ruling crossed
                Eigen::Matrix2d A; Eigen::Vector2d b(ValidFC[i]->first->p.x() - ValidFC[i+1]->first->p.x(), ValidFC[i]->first->p.y() - ValidFC[i+1]->first->p.y());
                A(0,0) = ValidFC[i]->first->p.x() - ValidFC[i]->second->p.x();
                A(0,1) = -(ValidFC[i+1]->first->p.x() - ValidFC[i+1]->second->p.x());
                A(1,0) = ValidFC[i]->first->p.y() - ValidFC[i]->second->p.y();
                A(1,1) = -(ValidFC[i+1]->first->p.y() - ValidFC[i+1]->second->p.y());
                Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
                Eruling +=  1.0/ts(0);
            }
            ofs << _a << " , " << MathTool::rad2deg(_a) << " , " << area << " , " << Ebend << " , " << Eruling << " , " << std::abs(acon - _a) <<  std::endl;
        }
        ofs.close();
    }

    for(int i = 1; i < (int)ValidFC.size()-2; i++){
        Eigen::Vector3d e = (ValidFC[i+1]->first->p - ValidFC[i]->first->p).normalized(), r = (ValidFC[i]->second->p - ValidFC[i]->first->p).normalized(), r2 = (ValidFC[i+1]->second->p - ValidFC[i+1]->first->p).normalized();   
        double phi1 = std::acos(e.dot(r)), phi2 = std::acos(r2.dot(-e));
        std::shared_ptr<Vertex> p = getCrossPoint(sin(phi2)/sin(phi1+phi2)*(ValidFC[i+1]->first->p - ValidFC[i]->first->p).norm(), ValidFC[i], false);
        if(p != nullptr)Triangles.push_back({ValidFC[i]->first, ValidFC[i+1]->first, p});

        r = (ValidFC[i]->third->p - ValidFC[i]->first->p).normalized(), r2 = (ValidFC[i+1]->third->p - ValidFC[i+1]->first->p).normalized();
        phi1 = std::acos(e.dot(r)); phi2 = std::acos(r2.dot(-e));
        p = getCrossPoint(sin(phi2)/sin(phi1+phi2)*(ValidFC[i+1]->first->p - ValidFC[i]->first->p).norm(), ValidFC[i], true);
        if(p != nullptr)Tri_fixside.push_back({ValidFC[i]->first, ValidFC[i+1]->first, p});
    }
    return Triangles;
}

void FoldLine::applyAAAMethod(const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool begincener, double a, double _tol, bool isroot){

    if(FoldingCurve.empty() && a_flap == -1)return;
    SimplifyModel(validsize, isroot);
    //if(!begincener)_FoldingAAAMethod(FoldingCurve, Poly_V, a);
    //else _FoldingAAAMethod_center(FoldingCurve, Poly_V, a);
    _FoldingAAAMethod(FoldingCurve, Poly_V, a);

    if(RulingsCrossed(FoldingCurve) < 1e-9){
        a_flap = a; //tol = _tol;
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

inline double update_flapangle2(double a, const Eigen::Vector3d& befN, const Eigen::Vector3d& SrfN, const Eigen::Vector3d& e){
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
    auto N = (Eigen::AngleAxisd(a, -e) * e.cross(e2)).normalized();
    //x->second->p3 = 20*N + x->first->p3;
    //qDebug()<<"front " << N.x() << " , " << N.y() <<" , " << N.z();
    double sin_a = sin(Phi[0])*sin(a)/sin(Phi[1]);
    if(sin_a > 1)sin_a = 1;
    else if(sin_a < -1)sin_a = -1;
    double cos_a = (cos(Phi[0]) - cos(Phi[1])*cos(beta))/(sin(Phi[1])*sin(beta));
    if(cos_a > 1)cos_a = 1;
    else if(cos_a < -1)cos_a = -1;
    a2 = std::atan2(sin_a, cos_a);
    if(a2 < 0)a2 += 2*std::numbers::pi;//koko chuui
    SrfN = MathTool::ProjectionVector(e.cross(e2), e2, true);
}

void CalcRuling_back(double a, std::shared_ptr<Vertex4d>& xbef, std::shared_ptr<Vertex4d>& x, std::shared_ptr<Vertex4d>& xnext, const std::vector<std::shared_ptr<Vertex>>& Poly_V,  double& a2, Eigen::Vector3d& SrfN, const Eigen::Vector3d &SpinAxis ){

    Eigen::Vector3d e = (xbef->first->p3 - x->first->p3).normalized(), e2 = (xnext->first->p3 - x->first->p3).normalized();
    Eigen::Vector3d axis = (x->third->p3 - x->first->p3).normalized();

    double beta = std::acos(e.dot(e2));
    double k = 2.0 * std::numbers::pi - std::acos(e2.dot(axis)) - std::acos(e.dot(axis));
    double _phi1 = std::atan((std::cos(k) - std::cos(beta))/(std::sin(beta)*std::cos(a)- std::sin(k)));
    double phi1 = (_phi1 < 0)? _phi1+std::numbers::pi: _phi1;
    double phi2 = k - phi1;
    Eigen::Vector3d N = (Eigen::AngleAxisd(a, e) * e.cross(e2)).normalized(); // 回転行列を作成
    Eigen::Vector3d r3d = (Eigen::AngleAxisd(phi1, -N) * e).normalized();
    Eigen::Vector3d r2d = (Eigen::AngleAxisd(phi1, SpinAxis)* (xbef->first->p- x->first->p)).normalized();//展開図のruling方向
    qDebug() <<"phi1 " << MathTool::rad2deg(phi1) << " , a = " << MathTool::rad2deg(a);
    Eigen::Vector3d crossPoint;
    for(int k = 0; k < (int)Poly_V.size(); k++){
        crossPoint = MathTool::calcCrossPoint_2Vector(x->first->p, 1000.0 * r2d + x->first->p, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p);
        if(MathTool::is_point_on_line(crossPoint, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) && r2d.dot(crossPoint - x->first->p) > 0)break;
    }
    x->second->p = crossPoint;
    x->second->p3 = (crossPoint - x->first->p).norm() * r3d + x->first->p3;
    //x->second->p3 = 20*N + x->first->p3;
    double sin_a = sin(phi1)*sin(a)/sin(phi2);
    if(sin_a > 1)sin_a = 1; else if(sin_a < -1)sin_a = -1;
    double cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
    if(cos_a > 1)cos_a = 1; else if(cos_a < -1)cos_a = -1;
    a2 = std::atan2(sin_a, cos_a);
    if(a2 < 0)a2 += 2*std::numbers::pi;//koko chuui
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

   double a2 = 0;
    Eigen::Vector3d e, e2;
    Eigen::Vector3d N, r, befN, befedge;
    std::vector<std::shared_ptr<Vertex4d>> FC;
    for(auto&fc:FoldingCurve){if(fc->IsCalc)FC.push_back(fc);}
    int mid = FC.size()/2;
    double a = (angle < 0.0)? angle + 2.0 * std::numbers::pi: (angle > 2.0*std::numbers::pi)? angle - 2.0*std::numbers::pi: angle;
    std::shared_ptr<Vertex4d> fc, fc_bef, fc_next;

    for(int ind = mid; ind < (int)FC.size() - 1; ind++){
        e = (FC[ind-1]->first->p3 - FC[ind]->first->p3).normalized();
        e2 = (FC[ind+1]->first->p3 - FC[ind]->first->p3).normalized();
        if(ind != mid){
            Eigen::Vector3d SrfN = MathTool::ProjectionVector(e.cross(e2), -e, true);
            a = update_flapangle(a2, befN, SrfN, e);
        }
        CalcRuling(a, FC[ind-1], FC[ind], FC[ind+1], Poly_V, a2, befN, -Eigen::Vector3d::UnitZ());
        if(ind == (int)FoldingCurve.size() -2)break;
        befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
    }

    a = (angle < 0.0)? angle + 2.0 * std::numbers::pi: (angle > 2.0*std::numbers::pi)? angle - 2.0*std::numbers::pi: angle;
    Eigen::Vector3d SrfN = MathTool::ProjectionVector(e.cross(e2), -e, true);
    befN = MathTool::ProjectionVector((FC[mid - 2]->first->p3 - FC[mid-1]->first->p3).normalized().cross(-e), -e, true);
    double dir = e.dot(befN.cross(SrfN));
    double tau = std::acos(SrfN.dot(befN));
    if(dir > 0)tau = std::numbers::pi - tau;
    if(a + tau > 2*std::numbers::pi)a= a + tau - 2.0 * std::numbers::pi;  else a = a + tau;
    befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
    for(int ind = mid -1; ind > 0; ind--){
        e = (FC[ind+1]->first->p3 - FC[ind]->first->p3).normalized();
        e2 = (FC[ind-1]->first->p3 - FC[ind]->first->p3).normalized();
        SrfN = MathTool::ProjectionVector(e.cross(e2), -e, true);
        if(ind != mid-1){
            dir = e.dot(befN.cross(SrfN));
            tau = std::acos(SrfN.dot(befN));
            if(dir <= 0)tau = 2.0*std::numbers::pi - tau;
            if(a2 - tau <= 0)a= a2 - tau + 2.0 * std::numbers::pi;  else a = a2 - tau;
        }
        CalcRuling_back(a, FC[ind+1], FC[ind], FC[ind-1], Poly_V, a2, befN, Eigen::Vector3d::UnitZ());
        Eigen::Vector3d r = (FC[ind]->second->p3 - FC[ind]->first->p3).normalized();
        Eigen::Vector3d axis = (FC[ind]->third->p3 - FC[ind]->first->p3).normalized();
        double phi1 = std::acos(r.dot(e)), phi2 = std::acos(r.dot(e2)), phi3 = std::acos(axis.dot(e)), phi4 = std::acos(axis.dot(e2));
        qDebug() <<ind << "  tau = " << MathTool::rad2deg(tau) << "  , a = " << MathTool::rad2deg(a) << " , phi1 = " << MathTool::rad2deg(phi1) <<
            ", phi2 = " << MathTool::rad2deg(phi2) << " ,  " << 2*std::numbers::pi - (phi1 + phi2 + phi3 + phi4);

        if(ind == 1)break;
        befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
    }
    qDebug() <<"====================";
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
    Eigen::Vector3d N, r, befN, befedge;
    std::vector<int> Vertices_Ind;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i]->IsCalc)Vertices_Ind.push_back(i);}

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
            Eigen::Vector3d edge = MathTool::ProjectionVector(e2, -e, true);
            a = update_flapangle(a2, befN, SrfN, e);
        }
        CalcRuling(a, fc_bef, fc, fc_next, Poly_V, a2, befN, -Eigen::Vector3d::UnitZ());
        if(ind == (int)FoldingCurve.size() -2)break;
        befN = MathTool::ProjectionVector(e.cross(e2), e2, true);befedge = MathTool::ProjectionVector(e, e2, true);
        if(ind != 1){
            for(int i = Vertices_Ind[ind-1] + 1; i < Vertices_Ind[ind]; i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind[ind]],FoldingCurve[Vertices_Ind[ind-1]], FoldingCurve[Vertices_Ind[ind]]->second,  FoldingCurve[Vertices_Ind[ind-1]]->second, 0);
        }
    }
    return;
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
            if(j == (int)Vertices_Ind.size())qDebug() <<"can't find correct edge right" ;
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
            FoldingCurve[Vertices_Ind.back()]->second->p3 = calcTargetDistanceOnPlane(FoldingCurve[Vertices_Ind.back()]->second->p, FoldingCurve[Vertices_Ind.end()[-2]]->second, FoldingCurve[Vertices_Ind.back()]->first, FoldingCurve[Vertices_Ind.end()[-2]]->first);
            for(int i = Vertices_Ind.end()[-2] + 1; i < Vertices_Ind.back(); i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind.back()], FoldingCurve[Vertices_Ind.end()[-2]], FoldingCurve[Vertices_Ind.back()]->second, FoldingCurve[Vertices_Ind.end()[-2]]->second, -1);
        }else{
            int j;
            for(j = 0; j < (int)Vertices_Ind.size(); j++){
                if(FoldingCurve[Vertices_Ind[j]]->second->p == v_clst->p)break;
            }
            if(j == (int)Vertices_Ind.size())qDebug() <<"can't find correct edge right" ;
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

