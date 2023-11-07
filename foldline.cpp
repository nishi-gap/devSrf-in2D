#include "foldline.h"

const double w_at = 1e-4;
const double eps = 1e-7;
double GetFlapAngle(int ind, const std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, double a_init, int StartingIndex);
inline Eigen::Vector3d _calcruling3d(const double& a, Eigen::Vector3d e, Eigen::Vector3d e2, Eigen::Vector3d axis, double& beta, std::vector<double>& Phi);
void CalcRuling(double a, std::shared_ptr<Vertex4d>& xbef, std::shared_ptr<Vertex4d>& x, std::shared_ptr<Vertex4d>& xnext, const std::vector<std::shared_ptr<Vertex>>& Poly_V, double& a2, Eigen::Vector3d& SrfN, const Eigen::Vector3d& SpinAixs);
void CalcRuling_back(double a, std::shared_ptr<Vertex4d>& xbef, std::shared_ptr<Vertex4d>& x, std::shared_ptr<Vertex4d>& xnext, const std::vector<std::shared_ptr<Vertex>>& Poly_V, double& a2, Eigen::Vector3d& SrfN, const Eigen::Vector3d& SpinAixs);
inline double update_flapangle(double a, const Eigen::Vector3d& befN, const Eigen::Vector3d& SrfN, const Eigen::Vector3d& e);
inline double update_flapangle2(double a, const Eigen::Vector3d& befN, const Eigen::Vector3d& SrfN, const Eigen::Vector3d& e);
void _FoldingAAAMethod_left(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle);
void _FoldingAAAMethod_center(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle, int StartingIndex);
void _FoldingAAAMethod_right(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle, int StartingIndex);
std::vector<std::shared_ptr<Vertex4d>> TrimPoints(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, double tol);
bool IsRulingCrossed(Eigen::Vector3d N, Eigen::Vector3d& cp, Eigen::Vector3d& crossPoint,  const std::vector<std::shared_ptr<Vertex>>& Poly_V);
void Douglas_Peucker_algorithm(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, std::vector<std::shared_ptr<Vertex4d>>& res, double tol = std::numbers::pi/9.0);

inline double logistic(double x, double a = 1){return 1.0/(1+exp(-a*x));}

void cb_Folding(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle, int StartingPoint){
    if(StartingPoint == (int)FoldingCurve.size()){_FoldingAAAMethod_right(FoldingCurve, Poly_V, angle, StartingPoint - 1);}//右手系のみを使った回転
    else _FoldingAAAMethod_center(FoldingCurve, Poly_V, angle, StartingPoint);
}

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

inline void update_Vertex(double t, std::shared_ptr<Vertex4d>&tar, const std::shared_ptr<Vertex4d>& src){
    tar->first->p = t * (src->first->p - src->third->p).normalized() + src->third->p;
    tar->first->p3 = t * (src->first->p3 - src->third->p3).normalized() + src->third->p3;
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
        FoldLine3d BasePt;
        int StartingIndex;
        double wb, wp ,wr;
        double mu;
        void AddWeight(double _wb, double _wp, double _wr, double _mu){ wb = _wb; wp = _wp; wr = _wr; mu = _mu;}
        void AddBasePt(std::vector<std::shared_ptr<Vertex4d>>& FC){BasePt = FC;}
        OptimizeParam(FoldLine3d& _FC,  const std::vector<std::shared_ptr<Vertex>>& _Poly_V, int _StartingIndex):FC{_FC}, Poly_V{_Poly_V}, StartingIndex(_StartingIndex){}
        ~OptimizeParam(){}
    private:
    };

    class OptimizeParam_v: public OptimizeParam{
    public:
        OptimizeParam_v() = delete;
        OptimizeParam_v(double _a, FoldLine3d& _FC, const std::vector<std::shared_ptr<Vertex>>& Poly_V, int _StartingIndex): a(_a), OptimizeParam::OptimizeParam( _FC,  Poly_V, _StartingIndex){}
        ~OptimizeParam_v(){}
        double a;
        int ind_l, ind_r;
        void SetCrossedIndex(int l, int r){ind_l = l; ind_r = r;}
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
    double getK(const std::shared_ptr<Vertex>&xbef, const std::shared_ptr<Vertex4d>&x ,const std::shared_ptr<Vertex>&xnext){
        Eigen::Vector3d e = (xbef->p - x->first->p).normalized(), e2 = (xnext->p - x->first->p).normalized(), axis = (x->third->p - x->first->p).normalized();
        return 2*std::numbers::pi - std::acos(e.dot(axis)) - std::acos(e2.dot(axis));
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
}

void FoldLine::AlignmentVertex4dDirection(){
    if(FoldingCurve.size() >= 3){
            if((FoldingCurve[0]->second->p - FoldingCurve[0]->third->p).normalized().dot((FoldingCurve[1]->second->p - FoldingCurve[1]->third->p).normalized()) < 0)
                std::swap(FoldingCurve[0]->second, FoldingCurve[0]->third);
            if((FoldingCurve.back()->second->p - FoldingCurve.back()->third->p).normalized().dot((FoldingCurve.end()[-2]->second->p - FoldingCurve.end()[-2]->third->p).normalized()) < 0)
                std::swap(FoldingCurve.back()->second, FoldingCurve.back()->third);
    }
}

void FoldLine::initialize_foldstate(bool IsStartEnd, const std::vector<std::shared_ptr<Vertex>>& Poly_V){
    if((int)FoldingCurve.size() < 3)return;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    int StartingPoint = (IsStartEnd)? 0: FoldingCurve.size()/2;
    Eigen::Vector3d e, e2, SpinAxis, Nb, N4;
    if(IsStartEnd){
            e = (ValidFC[0]->first->p3 - ValidFC[1]->first->p3).normalized(),e2 = (ValidFC[2]->first->p3 - ValidFC[1]->first->p3).normalized();
            SpinAxis = (ValidFC[1]->third->p3 - ValidFC[1]->first->p3).normalized();
            Nb = -(e.cross(e2)).normalized(), N4 = (e.cross(SpinAxis)).normalized();
    }else{
            int mid = ValidFC.size()/2;
            e = (ValidFC[mid-1]->first->p3 - ValidFC[mid]->first->p3).normalized(),e2 = (ValidFC[mid+1]->first->p3 - ValidFC[mid]->first->p3).normalized();
            SpinAxis = (ValidFC[mid]->third->p3 - ValidFC[mid]->first->p3).normalized();
            Nb = -(e.cross(e2)).normalized(), N4 = (e.cross(SpinAxis)).normalized();
    }

    double a_ll = std::atan2(e.dot(Nb.cross(N4)), N4.dot(-Nb));
    if(a_ll < 0)a_ll += 2.0*std::numbers::pi;
    double a_con = a_ll + std::numbers::pi;
    if(a_con > 2.0*std::numbers::pi)a_con -= 2.0*std::numbers::pi;
    if(a_con < 0)a_con +=2.0*std::numbers::pi;
    qDebug() << "a_con = " << MathTool::rad2deg(a_con);
    cb_Folding(FoldingCurve, Poly_V, a_con, StartingPoint);
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
            Eigen::Vector3d e = (ValidFC[i2]->first->p - ValidFC[i]->first->p), r = (ValidFC[i]->second->p - ValidFC[i]->first->p).normalized(), r2 = (ValidFC[i2]->second->p - ValidFC[i2]->first->p).normalized();
            double phi1 = std::acos(r.dot(e.normalized())), phi2 = std::acos(r2.dot(-e.normalized()));
            double l = sin(phi2)/sin(phi1+phi2)*(ValidFC[i2]->first->p - ValidFC[i]->first->p).norm();
            double lmax = (ValidFC[i]->second->p - ValidFC[i]->first->p).norm();
            if(0 <= l && l < lmax){_f += lmax/l; w++;}
            continue;
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
        if(IsSameWeight)f += abs(s);
        else f += abs((std::sin(std::numbers::pi - (ValidFC.size()/2  - i)*std::numbers::pi) + 0.1)*s);
    }
    if(!IsSameWeight){
        Eigen::Vector3d e = (ValidFC[2]->first->p - ValidFC[1]->first->p), r = (ValidFC[1]->second->p - ValidFC[1]->first->p).normalized(), r2 = (ValidFC[2]->second->p - ValidFC[2]->first->p).normalized();
        double phi1 = std::acos(r.dot(e.normalized())), phi2 = std::acos(r2.dot(-e.normalized()));
        double l = sin(phi2)/sin(phi1+phi2)*(ValidFC[2]->first->p - ValidFC[1]->first->p).norm();
        double s = 2.0/(e.cross(l*r)).norm();
        f += abs((std::sin(std::numbers::pi - ValidFC.size()/2*std::numbers::pi) + 0.1)*s);

        e = (ValidFC.end()[-2]->first->p - ValidFC.end()[-3]->first->p), r = (ValidFC.end()[-3]->second->p - ValidFC[1]->first->p).normalized(), r2 = (ValidFC.end()[-2]->second->p - ValidFC.end()[-2]->first->p).normalized();
        phi1 = std::acos(r.dot(e.normalized())), phi2 = std::acos(r2.dot(-e.normalized()));
        l = sin(phi2)/sin(phi1+phi2)*(ValidFC.end()[-2]->first->p - ValidFC.end()[-3]->first->p).norm();
        s = 2.0/(e.cross(l*r)).norm();
        f += abs((std::sin(std::numbers::pi - ValidFC.size()/2*std::numbers::pi) + 0.1)*s);
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

bool FoldLine::isbend(){
    double fr = RulingsCrossed(FoldingCurve);
    return (fr < 1e-9 && a_flap != -1)? true: false;
}

double Fruling(const std::vector<double> &X, std::vector<double> &grad, void* f_data)
{
    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    std::vector<std::shared_ptr<Vertex>> Poly_V = od->Poly_V;
    cb_Folding(od->FC, Poly_V, X[0], od->StartingIndex);
    double f = RulingsCrossed(od->FC);
    std::vector<double> a = X;
    if(!grad.empty()){
        double fp, fm;
        for(int i = 0; i < (int)X.size(); i++){
            a[i] = X[i] + eps;
            cb_Folding(od->FC, Poly_V, a[0], od->StartingIndex); fp = RulingsCrossed(od->FC);
            a[i] = X[i] - eps;
            cb_Folding(od->FC, Poly_V, a[0], od->StartingIndex);  fm = RulingsCrossed(od->FC);

            grad[i] = (fp - fm)/(2.0 * eps);
        }
    }

    if(DebugMode::Singleton::getInstance().isdebug())qDebug() <<"constraint function = " << MathTool::rad2deg(a[0]) << "(" << a[0] << ")  , " << f ;
    return f;
 }

inline double calcCrossPt4Constraint(const std::shared_ptr<Vertex4d>& x, const std::shared_ptr<Vertex4d>& x2){
    Eigen::Matrix2d A; Eigen::Vector2d b;
    A(0,0) = x->first->p.x() - x->second->p.x(); A(0,1) = -(x2->first->p.x() - x2->second->p.x());
    A(1,0) = x->first->p.y() - x->second->p.y(); A(1,1) = -(x2->first->p.y() - x2->second->p.y());
    b(0) = x->first->p.x() - x2->first->p.x(); b(1) = x->first->p.y() - x2->first->p.y();
    Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
    return (0 < ts(0) && ts(0) < 1) ? 1.0/ts(0): 0;
}

std::vector<double> getDynamicWeight(const std::vector<std::shared_ptr<Vertex4d>>& FC ,int StartingIndex){
    //動的な重みづけ
    //xiが直接影響するrulingの交差がない→重み0, 交差あり→中心から離れるほど重みを小さくして与える
    double mid = ((int)FC.size() % 2 == 1) ? static_cast<double>(StartingIndex): static_cast<double>(StartingIndex) - 0.5;
    std::vector<double> W(FC.size(), 0);
    for(int i = StartingIndex + 1; i < (int)FC.size(); i++){
        if(i > (int)FC.size() - 3) W[i] += (calcCrossPt4Constraint(FC.end()[-2],FC.end()[-3]) == 0)? 0: abs(std::pow(mid - (double)i/(mid), 3));
        else W[i] += (calcCrossPt4Constraint(FC[i],FC[i-1]) == 0)? 0: abs(std::pow(mid - (double)i/(mid), 3));
    }

    for(int i = 0; i < StartingIndex; i++){
        if(i < 2)W[i] += (calcCrossPt4Constraint(FC[1],FC[2]) == 0) ? 0 : abs(std::pow(mid - (double)i/(mid), 3));
        else W[i] += (calcCrossPt4Constraint(FC[i],FC[i+1]) == 0)? 0: abs(std::pow(mid - (double)i/(mid), 3));
    }
    return W;
}

double f_intersections(const std::vector<std::shared_ptr<Vertex4d>>& FC, int StartingIndex){
    double f = 0.0;

    for(int i = StartingIndex + 1; i < (int)FC.size(); i++){
        if(i > (int)FC.size() - 3) f += abs(std::pow((double)StartingIndex/(double)(StartingIndex - i), 3)) * calcCrossPt4Constraint(FC.end()[-2],FC.end()[-3]);
        else f += abs(std::pow((double)StartingIndex/(double)(StartingIndex - i), 3))* calcCrossPt4Constraint(FC[i],FC[i-1]);
    }
    for(int i = 0; i < StartingIndex - 1; i++){
        if(i < 2)f += abs(std::pow((double)StartingIndex/(double)(StartingIndex - i), 3))* calcCrossPt4Constraint(FC[1],FC[2]);
        else f += abs(std::pow((double)(StartingIndex)/(double)(abs(StartingIndex - i) + 1), 3))* calcCrossPt4Constraint(FC[i],FC[i+1]);
    }
    return f;
}

double f_intersections(const std::vector<std::shared_ptr<Vertex4d>>& FC, const std::vector<double>& W, int StartingIndex){
    double f = 0.0;
    for(int i = StartingIndex + 1; i < (int)FC.size(); i++){
        if(i > (int)FC.size() - 3) f += W[i]* calcCrossPt4Constraint(FC.end()[-2],FC.end()[-3]);
        else f += W[i]* calcCrossPt4Constraint(FC[i],FC[i-1]);
    }
    for(int i = 0; i < StartingIndex - 1; i++){
        if(i < 2)f += W[i]* calcCrossPt4Constraint(FC[1],FC[2]);
        else f += W[i]* calcCrossPt4Constraint(FC[i],FC[i+1]);
    }
    return f;
}

double Fconst_mixed(const std::vector<double> &X, std::vector<double> &grad, void* f_data){
    RevisionVertices::ObjData_v *od = (RevisionVertices::ObjData_v *)f_data;

    auto Fconv = [&](const std::vector<std::shared_ptr<Vertex4d>>& FC, const std::vector<std::shared_ptr<Vertex4d>>& BasePt, int s, int e){
        double f = 0.0;
        for(int i = s; i <= e; i++){
            double kori = RevisionVertices::getK(BasePt[i-1]->first, BasePt[i], BasePt[i+1]->first);
            double k = RevisionVertices::getK(FC[i-1]->first, FC[i], FC[i+1]->first);
            double tmp = (k - std::numbers::pi)*(kori - std::numbers::pi);
            if(tmp < 0)f += abs(tmp);
        }
        return f;
    };

    auto Fruling = [&](const std::vector<std::shared_ptr<Vertex4d>>& FC, int s){
        double f = 0.0;
        int mid = FC.size()/2;
        if(mid < s){//左側
            for(int i = mid + 1; i <= s; i++)f += calcCrossPt4Constraint(FC[i],FC[i-1]);
        }else{//右側
            for(int i = mid - 1; i >= s; i--)f += calcCrossPt4Constraint(FC[i],FC[i+1]);
        }

        return f;
    };

    if((int)X.size() == 2){
        update_Vertex(X[0], od->FC[od->ind_r - 1], od->BasePt[od->ind_r - 1]);
        update_Vertex(X[1], od->FC[od->ind_l + 1], od->BasePt[od->ind_l + 1]);
        cb_Folding(od->FC, od->Poly_V, od->a, od->StartingIndex);
        double fl_r = Fruling(od->FC, od->ind_l);
        double fl_c = Fconv(od->FC, od->BasePt, od->StartingIndex, od->ind_l);
        double fr_r = Fruling(od->FC, od->ind_r);
        double fr_c = Fconv(od->FC, od->BasePt, od->ind_r, od->StartingIndex );
        return fl_r + fl_c + fr_r + fr_c;
    }
    if((int)X.size() == 1){
        double fr, fc;
        if(od->ind_l != -1){
            update_Vertex(X[0], od->FC[od->ind_l + 1], od->BasePt[od->ind_l + 1]);
            cb_Folding(od->FC, od->Poly_V, od->a, od->StartingIndex);
            fr = Fruling(od->FC, od->ind_l);
            fc = Fconv(od->FC, od->BasePt, od->StartingIndex, od->ind_l);
        }else{
            update_Vertex(X[0], od->FC[od->ind_r - 1], od->BasePt[od->ind_r - 1]);
            cb_Folding(od->FC, od->Poly_V, od->a, od->StartingIndex);
            fr = Fruling(od->FC, od->ind_r);
            fc = Fconv(od->FC, od->BasePt, od->ind_r, od->StartingIndex);
        }
        return fr + fc;
    }

    //交点の導出
    for(int i = 0; i < static_cast<int>(X.size()); i++){
        update_Vertex(X[i], od->FC[i], od->BasePt[i]);
    }
    cb_Folding(od->FC, od->Poly_V, od->a, od->StartingIndex);

    const double wruling = 1e-0, wconv = 1e+4;//rulingの交差の値が大きいため抑制する
    std::vector<double> W = getDynamicWeight(od->FC, od->StartingIndex);
    double fc = 0.0;
    for(int i = 1; i < (int)od->FC.size() - 1; i++){
        double kori = RevisionVertices::getK(od->BasePt[i-1]->first, od->BasePt[i], od->BasePt[i+1]->first);
        double k = RevisionVertices::getK(od->FC[i-1]->first, od->FC[i], od->FC[i+1]->first);
        double tmp = (k - std::numbers::pi)*(kori - std::numbers::pi);
        if(tmp < 0)fc += wconv*abs(tmp);
    }
    double fr = wruling * f_intersections(od->FC, od->StartingIndex);

    std::vector<double> a = X;
    if(!grad.empty()){
        for(int i = 0; i < (int)X.size(); i++){
            double frp = 0.0, frm  = 0.0, fcp = 0.0, fcm = 0.0;
            update_Vertex(X[i] + eps, od->FC[i], od->BasePt[i]);
            cb_Folding(od->FC, od->Poly_V, od->a, od->StartingIndex);
            for(int j = 1; j < (int)od->FC.size() -1;j++){
                double k = RevisionVertices::getK(od->FC[j-1]->first, od->FC[j], od->FC[j+1]->first);
                double _k = RevisionVertices::getK(od->BasePt[j-1]->first, od->BasePt[j], od->BasePt[j+1]->first);
                double tmp = (k - std::numbers::pi)*(_k - std::numbers::pi);
                if(tmp < 0)fcp += wconv*abs(tmp);
            }
            frp = wruling * f_intersections(od->FC, od->StartingIndex);

            update_Vertex(X[i] - eps, od->FC[i], od->BasePt[i]);
            cb_Folding(od->FC, od->Poly_V, od->a, od->StartingIndex);
            for(int j = 1; j < (int)od->FC.size() -1;j++){
                double k = RevisionVertices::getK(od->FC[j-1]->first, od->FC[j], od->FC[j+1]->first);
                double _k = RevisionVertices::getK(od->BasePt[j-1]->first, od->BasePt[j], od->BasePt[j+1]->first);
                double tmp = (k - std::numbers::pi)*(_k - std::numbers::pi);
                if(tmp < 0)fcm += wconv*abs(tmp);
            }
            frm = wruling * f_intersections(od->FC, od->StartingIndex);
            update_Vertex(X[i], od->FC[i], od->BasePt[i]);
            grad[i] = W[i]*((fcp + frp) - (fcm + frm))/(2.0 * eps);

            a[i] = X[i];
        }
    }

    qDebug() <<"convex function = "  << fc << " , ruling function = " << fr;
    return fc + fr;
}

//比較に使う2つの交点の凹凸性が一致していることが前提
//凹凸が反転(変曲点をもつ)場合はひとまず無視
double F_sgnArea(const std::vector<std::shared_ptr<Vertex4d>>& FC, int StartingIndex){
    double f = 0.0;

    auto getphi = [](const std::shared_ptr<Vertex4d>& x, const std::shared_ptr<Vertex>& x2){
        Eigen::Vector3d e = (x2->p - x->first->p).normalized(), r = (x->second->p - x->first->p).normalized();
        return std::acos(e.dot(r));
    };

    double kbef = RevisionVertices::getK(FC[StartingIndex-1]->first, FC[StartingIndex], FC[StartingIndex+1]->first);
    for(int i = StartingIndex + 1; i < (int)FC.size() - 1; i++){
        double k = RevisionVertices::getK(FC[i-1]->first, FC[i], FC[i+1]->first);
        double phi1 = getphi(FC[i-1], FC[i]->first), phi2 = getphi(FC[i], FC[i-1]->first);
        if(kbef < std::numbers::pi && k < std::numbers::pi){
            if(phi1 + phi2 > std::numbers::pi) f += abs((phi1 + phi2) - std::numbers::pi);
        }
        else if(kbef > std::numbers::pi && k > std::numbers::pi){
            if(phi1 + phi2 < std::numbers::pi) f += abs((phi1 + phi2) - std::numbers::pi);
        }else{
            qDebug()<<"convex and concave is mixed";
        }
        kbef = k;
    }
    kbef = RevisionVertices::getK(FC[StartingIndex-1]->first, FC[StartingIndex], FC[StartingIndex+1]->first);
    for(int i = StartingIndex - 1; i > 0; i--){
        double k = RevisionVertices::getK(FC[i-1]->first, FC[i], FC[i+1]->first);
        double phi1 = getphi(FC[i+1], FC[i]->first), phi2 = getphi(FC[i], FC[i+1]->first);
        if(kbef < std::numbers::pi && k < std::numbers::pi){
            if(phi1 + phi2 > std::numbers::pi) f += abs((phi1 + phi2) - std::numbers::pi);
        }
        else if(kbef > std::numbers::pi && k > std::numbers::pi){
            if(phi1 + phi2 < std::numbers::pi) f += abs((phi1 + phi2) - std::numbers::pi);
        }else{
            qDebug()<<"convex and concave is mixed";
        }
        kbef = k;
    }

    return f;
}

double Fconst_conv(const std::vector<double> &X, std::vector<double> &grad, void* f_data){

    RevisionVertices::ObjData_v *od = (RevisionVertices::ObjData_v *)f_data;

    //元の折曲線のkを導出
    for(int i = 0; i < static_cast<int>(X.size()); i++){
        update_Vertex(X[i], od->FC[i], od->BasePt[i]);
    }

    cb_Folding(od->FC, od->Poly_V, od->a, od->StartingIndex);
    double f = 0.0;
    for(int i = 1; i < (int)od->FC.size() - 1; i++){
        double kori = RevisionVertices::getK(od->BasePt[i-1]->first, od->BasePt[i], od->BasePt[i+1]->first);
        double k = RevisionVertices::getK(od->FC[i-1]->first, od->FC[i], od->FC[i+1]->first);
        double tmp = (k - std::numbers::pi)*(kori - std::numbers::pi);
        if(tmp < 0)f += abs(tmp);
    }

    std::vector<double> W = getDynamicWeight(od->FC, od->StartingIndex);//交点位置に対する重みづけ

    std::vector<double> a = X;
    if(!grad.empty()){
        for(int i = 0; i < (int)X.size(); i++){
            double fp = 0.0, fm  = 0.0;
            a[i] = X[i] + eps;
            od->FC[i]->first->p =  a[i]*(od->BasePt[i]->first->p - od->BasePt[i]->third->p).normalized() + od->BasePt[i]->third->p;
            for(int j = 1; j < (int)od->FC.size() -1;j++){
                double k = RevisionVertices::getK(od->FC[j-1]->first, od->FC[j], od->FC[j+1]->first);
                double _k = RevisionVertices::getK(od->BasePt[j-1]->first, od->BasePt[j], od->BasePt[j+1]->first);
                double tmp = (k - std::numbers::pi)*(_k - std::numbers::pi);
                if(tmp < 0)fp += abs(tmp);
            }

            a[i] = X[i] - eps;
            od->FC[i]->first->p =  a[i]*(od->BasePt[i]->first->p - od->BasePt[i]->third->p).normalized() + od->BasePt[i]->third->p;
            for(int j = 1; j < (int)od->FC.size() -1;j++){
                double k = RevisionVertices::getK(od->FC[j-1]->first, od->FC[j], od->FC[j+1]->first);
                double _k = RevisionVertices::getK(od->BasePt[j-1]->first, od->BasePt[j], od->BasePt[j+1]->first);
                double tmp = (k - std::numbers::pi)*(_k - std::numbers::pi);
                if(tmp < 0)fm += abs(tmp);
            }

            od->FC[i]->first->p =  X[i]*(od->BasePt[i]->first->p - od->BasePt[i]->third->p).normalized() + od->BasePt[i]->third->p;
            grad[i] = W[i]*(fp - fm)/(2.0 * eps);

            a[i] = X[i];
        }
    }
    qDebug() <<"convex function = "  << f ;
    return f;
}

double Fruling4Vertex(const std::vector<double> &X, std::vector<double> &grad, void* f_data)
{
    //rulingの交点を導出
    //偏微分による勾配ベクトルでは交点xiを変化させて直接影響のあるruling ri, ri-1, ri+1の交点結果についてのみ計算する
    //直接影響のある場所のみ修正する
    //ステップサイズは中心ほど大きく、端に近いほど小さくする
    //本来だとxi以降のrulingの交点にも影響が出るためそれ以降のruling計算についても非線形な重みをかけて値を足すべきだと思うけど実装を簡単にするため

    RevisionVertices::ObjData_v *od = (RevisionVertices::ObjData_v *)f_data;
    std::vector<std::shared_ptr<Vertex>> Poly_V = od->Poly_V;
    for(int i = 0; i < static_cast<int>(X.size()); i++){
        od->FC[i]->first->p = X[i] * (od->BasePt[i]->first->p - od->BasePt[i]->third->p).normalized() + od->BasePt[i]->third->p;
        od->FC[i]->first->p3 = X[i] * (od->BasePt[i]->first->p3 - od->BasePt[i]->third->p3).normalized() + od->BasePt[i]->third->p3;
    }
    cb_Folding(od->FC, Poly_V, od->a, od->StartingIndex);

    std::vector<double> W = getDynamicWeight(od->FC, od->StartingIndex);
    double f = 1e-2*f_intersections(od->FC, W, od->StartingIndex);

    std::vector<double> a = X;
    if(!grad.empty()){
        double fp, fm;
        for(int i = 0; i < (int)X.size(); i++){
            a[i] = X[i] + eps;
            od->FC[i]->first->p = a[i] * (od->BasePt[i]->first->p - od->BasePt[i]->third->p).normalized() + od->BasePt[i]->third->p;
            od->FC[i]->first->p3 = a[i] * (od->BasePt[i]->first->p3 - od->BasePt[i]->third->p3).normalized() + od->BasePt[i]->third->p3;
            cb_Folding(od->FC, Poly_V, od->a, od->StartingIndex);fp = f_intersections(od->FC, od->StartingIndex);

            a[i] = X[i] - eps;
            od->FC[i]->first->p = a[i] * (od->BasePt[i]->first->p - od->BasePt[i]->third->p).normalized() + od->BasePt[i]->third->p;
            od->FC[i]->first->p3 = a[i] * (od->BasePt[i]->first->p3 - od->BasePt[i]->third->p3).normalized() + od->BasePt[i]->third->p3;
            cb_Folding(od->FC, Poly_V, od->a, od->StartingIndex);  fm = f_intersections(od->FC, od->StartingIndex);

            od->FC[i]->first->p = X[i] * (od->BasePt[i]->first->p - od->BasePt[i]->third->p).normalized() + od->BasePt[i]->third->p;
            od->FC[i]->first->p3 = X[i] * (od->BasePt[i]->first->p3 - od->BasePt[i]->third->p3).normalized() + od->BasePt[i]->third->p3;

            grad[i] = W[i] * (fp - fm)/(2.0 * eps);
            a[i] = X[i];
        }
    }

    qDebug() <<"ruling function = "  << f ;
    return f;
}

double ObjFunc_RulingIntersection(const std::vector<double> &X, std::vector<double> &grad, void* f_data){
    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    std::vector<std::shared_ptr<Vertex>> Poly_V = od->Poly_V;
    for(int i = 1; i < static_cast<int>(X.size()); i++){
        od->FC[i-1]->first->p =  X[i]*(od->BasePt[i-1]->first->p - od->BasePt[i-1]->third->p) + od->FC[i-1]->third->p;
        od->FC[i-1]->first->p3 =  X[i]*(od->BasePt[i-1]->first->p3 - od->BasePt[i-1]->third->p3) + od->FC[i-1]->third->p3;
    }
    cb_Folding(od->FC, Poly_V, X[0], od->StartingIndex);
    double fb =  Fbend2(od->FC), fp = Fparallel(od->FC),fr = (od->wr != -1)?RulingsCrossed_NoOutline(od->FC): 0;

    std::vector<double> a = X;
    if(!grad.empty()){
       for(int i = 0; i < (int)X.size(); i++){
            a[i] = X[i] + eps;
            if(i != 0){
                od->FC[i-1]->first->p =  a[i]*(od->BasePt[i-1]->first->p - od->BasePt[i-1]->third->p) + od->FC[i-1]->third->p;
                od->FC[i-1]->first->p3 =  a[i]*(od->BasePt[i-1]->first->p3 - od->BasePt[i-1]->third->p3) + od->FC[i-1]->third->p3;
            }
            cb_Folding(od->FC, Poly_V, a[0], od->StartingIndex);
            double fp = Fbend2(od->FC), fp2 = Fparallel(od->FC), fp3 = (od->wr != -1)?RulingsCrossed_NoOutline(od->FC): 0;

            a[i] = X[i] - eps;
            if(i != 0){
                od->FC[i-1]->first->p =  a[i]*(od->BasePt[i-1]->first->p - od->BasePt[i-1]->third->p) + od->FC[i-1]->third->p;
                od->FC[i-1]->first->p3 =  a[i]*(od->BasePt[i-1]->first->p3 - od->BasePt[i-1]->third->p3) + od->FC[i-1]->third->p3;
            }
            cb_Folding(od->FC, Poly_V, a[0], od->StartingIndex);
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
    cb_Folding(od->FC, Poly_V, a[0], od->StartingIndex);
    fb =  Fbend2(od->FC); fp = Fparallel(od->FC); fr = Fmin_TriangleArea(od->FC, true);

    if(!grad.empty()){
       std::vector<double> X = a;
        for(int i = 0; i < (int)a.size(); i++){
        X[i] = a[i] + eps;

        cb_Folding(od->FC, Poly_V, X[0], od->StartingIndex);
        double fp3 = Fmin_TriangleArea(od->FC, true);
        X[i] = a[i] - eps;
        cb_Folding(od->FC, Poly_V, X[0], od->StartingIndex);
        double fm3 = Fmin_TriangleArea(od->FC, true);
        double dfa = 1.0*(fp3 - fm3)/(2.0*eps);
        grad[i] = dfa;
       }

    }
    if(DebugMode::Singleton::getInstance().isdebug())qDebug() <<"Triangle Area(" << MathTool::rad2deg(a[0]) << "(" << a[0] << ")  =  " <<  fb  << ", Area = "<< fr;
    return od->wb * fb + 0 * fp + 1.0 * fr;
}

double ObjFunc_Vertex(const std::vector<double> &X, std::vector<double> &grad, void* f_data){
    const double w = 1e-5, wsim = 10, warea = 1e-2;

    RevisionVertices::ObjData_v *od = (RevisionVertices::ObjData_v *)f_data;
    std::vector<std::shared_ptr<Vertex>> Poly_V = od->Poly_V;
    int ind_l, ind_r;
    if((int)X.size() == 2){
       update_Vertex(X[0], od->FC[od->ind_r - 1], od->BasePt[od->ind_r - 1]);
       update_Vertex(X[1], od->FC[od->ind_l + 1], od->BasePt[od->ind_l + 1]);
       ind_l = std::max((int)od->FC.size(), od->ind_l + 1);
       ind_r = std::min(0, od->ind_r - 1);
    }
    if((int)X.size() == 1){
       if(od->ind_l != -1){
        ind_l = std::max((int)od->FC.size(), od->ind_l + 1); ind_r = 0;//ind_rの値が正しいか要検証
        update_Vertex(X[0], od->FC[od->ind_l + 1], od->BasePt[od->ind_l + 1]);
       }
       else {
        update_Vertex(X[0], od->FC[od->ind_r - 1], od->BasePt[od->ind_r - 1]);
        ind_r = std::min(0, od->ind_r - 1); ind_l = od->FC.size();//ind_lの値が正しいのか要検証
       }
    }

    cb_Folding(od->FC, Poly_V, od->a, od->StartingIndex);
    double fiso = 0.0;
    if((int)X.size() == 2){
       fiso += (od->BasePt[od->ind_r - 1]->first->p - od->FC[od->ind_r - 1]->first->p).norm();
       fiso += (od->BasePt[od->ind_l + 1]->first->p - od->FC[od->ind_l + 1]->first->p).norm();
    }else{
       if(od->ind_l != -1) fiso += (od->BasePt[od->ind_l + 1]->first->p - od->FC[od->ind_l + 1]->first->p).norm();
       else  fiso += (od->BasePt[od->ind_r - 1]->first->p -  od->FC[od->ind_r - 1]->first->p).norm();
    }
    std::vector<std::shared_ptr<Vertex4d>> subFC;
    for(int i = ind_r; i < ind_l; i++)subFC.push_back(od->FC[i]);
    double farea = Fmin_TriangleArea(subFC, true);
    return wsim * fiso + warea * farea;

    for(int i = 0; i < static_cast<int>(X.size()); i++){
       update_Vertex(X[i], od->FC[i], od->BasePt[i]);
    }
    std::vector<double> W = getDynamicWeight(od->FC, od->StartingIndex);

    double fsim = 0.0;
    //double fconv = wconv * Func_softconst(od->FC), fruling = wruilng * f_intersections(od->FC);
    for(int i = 0; i < (int)X.size(); i++){fsim +=  wsim * abs(X[i] - (od->BasePt[i]->first->p - od->BasePt[i]->third->p).norm());}
    //for(int i = 0; i < (int)X.size(); i++){fsim +=  wsim * abs(std::pow((double)(mid - i)/(double)mid, 3))*(RegCrv[i]-center).norm()/sum * abs(X[i]);}

    std::vector<double> a = X;
    if(!grad.empty()){
       for(int i = 0; i < (int)X.size(); i++){
        update_Vertex(X[i] + eps, od->FC[i], od->BasePt[i]);
        double fsimp = 0.0;
        for(int j = 0; j < (int)X.size(); j++){fsimp +=  wsim * abs(X[j] - (od->BasePt[j]->first->p - od->BasePt[j]->third->p).norm());}

        update_Vertex(X[i] - eps, od->FC[i], od->BasePt[i]);
        double fsimm = 0.0;
        for(int j = 0; j < (int)X.size(); j++){fsimm +=  wsim *abs(X[j] - (od->BasePt[j]->first->p - od->BasePt[j]->third->p).norm());}

        grad[i] = W[i] * w*((fsimp ) - (fsimm ))/(2.0 * eps);
        update_Vertex(X[i], od->FC[i], od->BasePt[i]);
       }
    }
    qDebug() <<"simlarity = " <<fsim ;
    return  fsim ;
}

double ObjFunc_AngleK(const std::vector<double> &X, std::vector<double> &grad, void* f_data){
    auto similarity_k = [&](const std::vector<std::shared_ptr<Vertex4d>>& FC, const std::vector<std::shared_ptr<Vertex4d>>& BasePt, std::vector<double>& W)->double{
        double f = 0.0, k, k2;
        for(int i = 0; i < (int)FC.size(); i++){
            if(i == 0){
                k = RevisionVertices::getK(BasePt[0]->first, BasePt[1], BasePt[2]->first);
                k2 = RevisionVertices::getK(FC[0]->first, FC[1], FC[2]->first);
            }else if(i == (int)X.size() - 1){
                k = RevisionVertices::getK(BasePt.end()[-1]->first, BasePt.end()[-2], BasePt.end()[-3]->first);
                k2 = RevisionVertices::getK(FC.end()[-1]->first, FC.end()[-2], FC.end()[-3]->first);
            }else{
                k = RevisionVertices::getK(BasePt[i-1]->first, BasePt[i], BasePt[i+1]->first);
                k2 = RevisionVertices::getK(FC[i-1]->first, FC[i], FC[i+1]->first);
            }
            f += W[i] * abs(k - k2);
        }
        return f;
    };

    RevisionVertices::ObjData_v *od = (RevisionVertices::ObjData_v *)f_data;
    std::vector<std::shared_ptr<Vertex>> Poly_V = od->Poly_V;
    for(int i = 0; i < static_cast<int>(X.size()); i++){
       update_Vertex(X[i], od->FC[i], od->BasePt[i]);
    }
    cb_Folding(od->FC, Poly_V, od->a, od->StartingIndex);
    std::vector<double> W = getDynamicWeight(od->FC, od->StartingIndex);

    double fsim = similarity_k(od->FC, od->BasePt, W);

    std::vector<double> a = X;
    if(!grad.empty()){
       for(int i = 0; i < (int)X.size(); i++){
        update_Vertex(X[i] + eps, od->FC[i], od->BasePt[i]);
        double fsimp = similarity_k(od->FC, od->BasePt, W);

        update_Vertex(X[i] - eps, od->FC[i], od->BasePt[i]);
        double fsimm = similarity_k(od->FC, od->BasePt, W);

        grad[i] = W[i] * ((fsimp ) - (fsimm ))/(2.0 * eps);
        update_Vertex(X[i], od->FC[i], od->BasePt[i]);
       }
    }
    qDebug() <<"simlarity of angle k = " <<fsim ;
    qDebug()<<"Weight = " << W[0] << " , " << W[1] << " , " << W[2] << " , " << W[3];
    return  fsim ;
}

void SetBounds(const std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V,
               int i, double& bnd_upper, double&bnd_lower, double& x){
    Eigen::Vector3d point;
    double range = 2.0;
    for(int k = 0; k < (int)Poly_V.size(); k++){
        point = MathTool::calcCrossPoint_2Vector(FoldingCurve[i]->first->p, 1000.0 * (FoldingCurve[i]->first->p - FoldingCurve[i]->third->p).normalized() + FoldingCurve[i]->first->p, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p);
        if(MathTool::is_point_on_line(point, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) && (FoldingCurve[i]->first->p - FoldingCurve[i]->third->p).normalized().dot(point - FoldingCurve[i]->first->p) > 0)break;
    }
    x = (FoldingCurve[i]->first->p - FoldingCurve[i]->third->p).norm();
    bnd_upper = std::min((point - FoldingCurve[i]->third->p).norm(), x + range);
    bnd_lower = std::max(0.0, x - range);
};

void update_rightside(std::vector<std::shared_ptr<Vertex4d>>& FC, const std::vector<std::shared_ptr<Vertex>>& Poly_V,
                      int ind_r, double x){
    double _d = (x - (FC[ind_r]->first->p - FC[ind_r]->third->p).norm());//後ろの交点も変化させる値(動かした交点同士の平均としておく)
    update_Vertex(x, FC[ind_r], FC[ind_r]);
    double tmp, up, low;
    for(int k = ind_r - 1; k >= 0; k--){
       double d = _d + (FC[k]->first->p -  FC[k]->third->p).norm();
       SetBounds(FC, Poly_V, k, up, low, tmp);
       d = (d < low)? low: (up < d)? up: d;
       update_Vertex(d, FC[k], FC[k]);
    }
}

void update_leftside(std::vector<std::shared_ptr<Vertex4d>>& FC, const std::vector<std::shared_ptr<Vertex>>& Poly_V,
                     int ind_l, double x){
    double _d = (x - (FC[ind_l]->first->p - FC[ind_l]->third->p).norm());//後ろの交点も変化させる値(動かした交点同士の平均としておく)
    update_Vertex(x, FC[ind_l], FC[ind_l]);
    double tmp, up, low;
    for(int k = ind_l + 1; k < (int)FC.size(); k++){
       double d = _d + (FC[k]->first->p -  FC[k]->third->p).norm();
       SetBounds(FC, Poly_V, k, up, low, tmp);
       d = (d < low)? low: (up < d)? up: d;
       update_Vertex(d, FC[k], FC[k]);
    }
}

//全探索による頂点最適化
//外側の頂点について振動させながら凹凸性とrulingの交差がない位置があれば終了
//外側の頂点で探索範囲内で見つからなければ内側の頂点も少しづつ振動させて見つける
//片側についてのみ探索するものとする←実装を簡単にするため＋検証目的
void FindFirstLocation(std::vector<double> &X, RevisionVertices::ObjData_v& od,
                const std::vector<double>& bnd_lower, const std::vector<double>& bnd_upper, bool& IsCrossed_l, bool& IsCrossed_r, int ind_l, int ind_r){

    auto Fconv = [&](const std::vector<std::shared_ptr<Vertex4d>>& FC, const std::vector<std::shared_ptr<Vertex4d>>& BasePt, int s, int e){
        double f = 0.0;
        for(int i = s; i <= e; i++){
            double kori = RevisionVertices::getK(BasePt[i-1]->first, BasePt[i], BasePt[i+1]->first);
            double k = RevisionVertices::getK(FC[i-1]->first, FC[i], FC[i+1]->first);
            double tmp = (k - std::numbers::pi)*(kori - std::numbers::pi);
            if(tmp < 0)f += abs(tmp);
        }
        return f;
    };

    auto Fruling = [&](const std::vector<std::shared_ptr<Vertex4d>>& FC, int s){
        double f = 0.0;
        int mid = FC.size()/2;
        if(mid < s){//左側
            for(int i = mid + 1; i <= s; i++)f += calcCrossPt4Constraint(FC[i],FC[i-1]);
        }else{//右側
            for(int i = mid - 1; i >= s; i--)f += calcCrossPt4Constraint(FC[i],FC[i+1]);
        }

        return f;
    };

    for(int i = 0; i < static_cast<int>(X.size()); i++) update_Vertex(X[i], od.FC[i], od.BasePt[i]);

    const double step = 1e-4;

    std::vector<double> minf_l{1e+10, 0,0}, minf_r{1e+10,0,0};
    double fl_r, fl_c, fr_r, fr_c;

    auto initloopCounter = [](const std::vector<double>& x, const std::vector<double>& bnd_lower, const std::vector<double>& bnd_upper, int i, double step){
        std::vector<double> X{x[i]};
        double tl = x[i] - step, tu = x[i] + step;
        do{
            if(tu < bnd_upper[i]){X.push_back(tu); tu += step;}
            if(bnd_lower[i] < tl){ X.push_back(tl); tl -=step;}
        }while(bnd_lower[i] <= tl || tu <= bnd_upper[i]);
        return X;
    };

    if(IsCrossed_l){
       std::vector<double> Xli = initloopCounter(X, bnd_lower, bnd_upper, ind_l, 1e-1), Xlo = initloopCounter(X, bnd_lower, bnd_upper, ind_l+1, step);
       for(int i = 0; i < 1 && IsCrossed_l; i++){
        //update_Vertex(Xli[i], od.FC[ind_l], od.BasePt[ind_l]);
            for(int j = 0; j < (int)Xlo.size() && IsCrossed_l; j++){
                update_Vertex(Xlo[j], od.FC[ind_l + 1], od.BasePt[ind_l + 1]);
                cb_Folding(od.FC, od.Poly_V, od.a, od.StartingIndex);
                fl_r = Fruling(od.FC, ind_l); fl_c = Fconv(od.FC, od.BasePt, od.StartingIndex, ind_l);
                if(minf_l[0] > fl_r){minf_l[0] = fl_r; minf_l[1] = fl_c;}
                if(fl_r == 0.0 && fl_c == 0.0){
                    //update_Vertex(Xli[i], od.BasePt[ind_l], od.BasePt[ind_l]);
                    update_Vertex(Xlo[j], od.BasePt[ind_l + 1], od.BasePt[ind_l + 1]);
                    double _d = (Xlo[j] - X[ind_l + 1]);//後ろの交点も変化させる値(動かした交点同士の平均としておく)
                    qDebug()<<"left side :  inner = " << Xli[i] - X[ind_l] << " , outer = " << Xlo[j] - X[ind_l + 1] << " , d = " <<_d;
                    for(int k = ind_l + 2; k < (int)X.size(); k++){
                        double d = _d + X[k];
                        d = (d < bnd_lower[k])? bnd_lower[k]: (bnd_upper[k] < d)? bnd_upper[k]: d;
                        update_Vertex(d, od.BasePt[k], od.BasePt[k]);
                    }
                    IsCrossed_l = false;
                    break;
                }
            }
       }
    }

    if(IsCrossed_r){
       std::vector<double> Xri = initloopCounter(X, bnd_lower, bnd_upper, ind_r, 1e-1), Xro = initloopCounter(X, bnd_lower, bnd_upper, ind_r - 1, step);
       for(int i = 0; i < 1 && IsCrossed_r; i++){
            //update_Vertex(Xri[i], od.FC[ind_r], od.BasePt[ind_r]);
            for(int j = 0; j < (int)Xro.size() && IsCrossed_r; j++){
                update_Vertex(Xro[j], od.FC[ind_r - 1], od.BasePt[ind_r - 1]);
                cb_Folding(od.FC, od.Poly_V, od.a, od.StartingIndex);
                fr_r = Fruling(od.FC, ind_r); fr_c = Fconv(od.FC, od.BasePt, ind_r, od.StartingIndex );
                if(minf_r[0] > fr_r){minf_r[0] = fr_r; minf_r[1] = fr_c;}
                if(fr_r == 0.0 && fr_c == 0.0){
                    //update_Vertex(Xri[i], od.BasePt[ind_r], od.BasePt[ind_r]);
                    update_Vertex(Xro[j], od.BasePt[ind_r - 1], od.BasePt[ind_r - 1]);
                    double _d = (Xro[j] - X[ind_r - 1]);//後ろの交点も変化させる値(動かした交点同士の平均としておく)
                    qDebug()<<"right side :  inner = " << Xri[i] - X[ind_r] << " , outer = " << Xro[j] - X[ind_r-1]<< " , d = " <<_d;
                    for(int k = ind_r - 2; k >= 0; k--){
                        double d = _d + X[k];
                        d = (d < bnd_lower[k])? bnd_lower[k]: (bnd_upper[k] < d)? bnd_upper[k]: d;
                        update_Vertex(d, od.BasePt[k], od.BasePt[k]);
                    }

                    IsCrossed_r = false;
                    break;
                }
            }
       }
    }
    qDebug() <<"result : minf_l ruling = " << minf_l[0] << " , conv = " << minf_l[1] << ", area = " << minf_l[2] << "  ,  minf_r ruling =  " << minf_r[0] << ", conv = " << minf_r[1] <<", area = " << minf_r[2];
    return;
}

void FullSearch(std::vector<double> &X, RevisionVertices::ObjData_v& od,
                const std::vector<double>& bnd_lower, const std::vector<double>& bnd_upper, bool& IsCrossed_l, bool& IsCrossed_r, int ind_l, int ind_r){
    auto Fconv = [&](const std::vector<std::shared_ptr<Vertex4d>>& FC, const std::vector<std::shared_ptr<Vertex4d>>& BasePt, int s, int e){
        double f = 0.0;
        for(int i = s; i <= e; i++){
            double kori = RevisionVertices::getK(BasePt[i-1]->first, BasePt[i], BasePt[i+1]->first);
            double k = RevisionVertices::getK(FC[i-1]->first, FC[i], FC[i+1]->first);
            double tmp = (k - std::numbers::pi)*(kori - std::numbers::pi);
            if(tmp < 0)f += abs(tmp);
        }
        return f;
    };

    auto Fruling = [&](const std::vector<std::shared_ptr<Vertex4d>>& FC, int s){
        double f = 0.0;
        int mid = FC.size()/2;
        if(mid < s){//左側
            for(int i = mid + 1; i <= s; i++)f += calcCrossPt4Constraint(FC[i],FC[i-1]);
        }else{//右側
            for(int i = mid - 1; i >= s; i--)f += calcCrossPt4Constraint(FC[i],FC[i+1]);
        }

        return f;
    };

    for(int i = 0; i < static_cast<int>(X.size()); i++) update_Vertex(X[i], od.FC[i], od.BasePt[i]);

    const double step = 1e-4;
    std::vector<double> minf_l{0, 1e+10, 1e+10, 1e+10}, minf_r{0,1e+10, 1e+10,1e+10};//変化量, rulingの交差判定, 凹凸性, 三角形の面積
    double fr, fc, fa;
    if(IsCrossed_l){
       for(double t = bnd_lower[ind_l+1]; t <= bnd_upper[ind_l+1]; t += step){
            update_Vertex(t, od.FC[ind_l+1], od.BasePt[ind_l+1]);
            cb_Folding(od.FC, od.Poly_V, od.a, od.StartingIndex);
            fr = Fruling(od.FC, ind_l); fc = Fconv(od.FC, od.BasePt, od.StartingIndex, ind_l);
            fa = Fmin_TriangleArea(od.FC, true);
            if(fr == 0.0 && fc == 0.0){//0であれば交差していないor凹凸性が失われていない
                IsCrossed_l = false;
                if(fa < minf_l[3]){//三角形の面積が大きい方を保持する
                    minf_l[0] = t; minf_l[1] = fr; minf_l[2] = fc; minf_l[3] = fa;
                }
            }
        }

       qDebug()<<"left side t = " << minf_l[0] << " , ruling = " << minf_l[1] << " , convesity = " << minf_l[2] << " , Farea = " << minf_l[3];
        double _d = minf_l[0] - X[ind_l+1];
        update_Vertex(minf_l[0], od.BasePt[ind_l+1], od.BasePt[ind_l+1]);
       for(int k = ind_l + 2; k < (int)X.size(); k++){
            double d = _d + X[k];
            d = (d < bnd_lower[k])? bnd_lower[k]: (bnd_upper[k] < d)? bnd_upper[k]: d;
            update_Vertex(d, od.BasePt[k], od.BasePt[k]);
       }

    }

    if(IsCrossed_r){
       for(double t = bnd_lower[ind_r-1]; t <= bnd_upper[ind_r-1]; t += step){
            update_Vertex(t, od.FC[ind_r-1], od.BasePt[ind_r-1]);
            cb_Folding(od.FC, od.Poly_V, od.a, od.StartingIndex);
            fr = Fruling(od.FC, ind_r); fc = Fconv(od.FC, od.BasePt, ind_r, od.StartingIndex);
            fa = Fmin_TriangleArea(od.FC, true);
            if(fr == 0.0 && fc == 0.0){
                IsCrossed_r = false;
                if(fa < minf_r[3]){//三角形の面積が大きい方を保持する
                    minf_r[0] = t; minf_r[1] = fr; minf_r[2] = fc; minf_r[3] = fa;
                }
            }
       }

       qDebug()<<"right side t = " << minf_r[0] << " , ruling = " << minf_r[1] << " , convesity = " << minf_r[2] << " , Farea = " << minf_r[3];
       double _d = minf_r[0] - X[ind_r-1];
       update_Vertex(minf_r[0], od.BasePt[ind_r-1], od.BasePt[ind_r-1]);
       for(int k = ind_r - 2; k >= 0; k--){
            double d = _d + X[k];
            d = (d < bnd_lower[k])? bnd_lower[k]: (bnd_upper[k] < d)? bnd_upper[k]: d;
            update_Vertex(d, od.BasePt[k], od.BasePt[k]);
       }
    }
}

void SetOptimizationParameter(nlopt::opt& opt, const std::vector<double>&bnd_lower, const std::vector<double>&bnd_upper, double ftol, double xtol, double maxtime){
    opt.set_lower_bounds(bnd_lower);
    opt.set_upper_bounds(bnd_upper);
    opt.set_xtol_rel(xtol);
    opt.set_ftol_rel(ftol);
    opt.set_maxtime(maxtime);//stop over this time
}

//片側の全探索による交点位置最適化
//i番目のrulingが交差した場合、i,i+1番目の交点位置を修正する(すべての交点位置を少しづつ動かして検証するのが正しいのかもしれないが、動かすのはi,i+1番目だけでよいという仮定)
//全探索により見つからない場合は境界条件を少しずらして再挑戦させる
bool FoldLine::ReviseVertexPos(const std::vector<std::shared_ptr<Vertex>>& Poly_V, int EndIndex_left, int EndIndex_right, int AlgOptim){
    int mid = FoldingCurve.size()/2;

    std::vector<double> bnd_lower((int)FoldingCurve.size(), -1), bnd_upper((int)FoldingCurve.size(), 1.0);
    std::vector<double> X(FoldingCurve.size());
    std::vector<std::shared_ptr<Vertex4d>> FC;
    double init_range = 3.0;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){
       std::shared_ptr<Vertex4d> v4d = FoldingCurve[i]->deepCopy(); v4d->addline(FoldingCurve[i]->line_parent); FC.push_back(v4d);
       X[i] = (FoldingCurve[i]->first->p - FoldingCurve[i]->third->p).norm();
       bnd_lower[i] = std::max(X[i] - init_range, 0.0);
       Eigen::Vector3d point;
       for(int k = 0; k < (int)Poly_V.size(); k++){
            point = MathTool::calcCrossPoint_2Vector(v4d->first->p, 1000.0 * (v4d->first->p - v4d->third->p).normalized() + v4d->first->p, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p);
            if(MathTool::is_point_on_line(point, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) && (v4d->first->p - v4d->third->p).normalized().dot(point - v4d->first->p) > 0)break;
       }
       bnd_upper[i] = std::min(X[i] + init_range, (point - v4d->third->p).norm());
    }

    RevisionVertices::ObjData_v od = {a_flap, FC, Poly_V, (int)FoldingCurve.size()/2};
    od.AddBasePt(FoldingCurve);
    bool IsCrossed_l = (EndIndex_left == (int)FoldingCurve.size() - 1)? false: true, IsCrossed_r= (EndIndex_right == 0)? false: true;
    int maxitr = 2;//最大5回境界条件を変えて検証する
    for(int i = 0; i < maxitr; i++){
       if(AlgOptim == 0)FindFirstLocation(X, od, bnd_lower,bnd_upper, IsCrossed_l, IsCrossed_r, EndIndex_left, EndIndex_right);
       else if(AlgOptim == 1)FullSearch(X, od, bnd_lower,bnd_upper, IsCrossed_l, IsCrossed_r, EndIndex_left, EndIndex_right);
       double fc_bef = RulingsCrossed(FoldingCurve);
       if(!IsCrossed_l || !IsCrossed_r){//片方だけでも交差除去できていたらtrueを返す
        cb_Folding(FoldingCurve, Poly_V, a_flap, mid);
        double fc = RulingsCrossed(FoldingCurve);
        qDebug() << "before optimization Fruling = " << fc_bef << "  , after fullsearch =  " << fc;
        return true;
       }else qDebug() << "not found in boundary condition";

       //境界条件のupdate
       //std::copy(bnd_upper.begin(),bnd_upper.end(), bnd_lower.begin());
       for(auto&b: bnd_upper)b+= 1.0;
    }

    return false;
}

bool FoldLine::PropagateOptimization_Vertex(const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool IsStartEnd, int OrderConstFunc, int OptimizationAlgorithm){
    if(IsStartEnd){
       qDebug()<<"this method is applied only starting from center";
       return false;
    }
      
    int bef_l = -1, bef_r = -1;
    int mid = FoldingCurve.size()/2;
    int ind_l = FoldingCurve.size() - 1, ind_r = 1;
    
    //交点を動かしたときの制約関数の値の変化をcsvファイルへ書き出し
    {
       int StartingIndex = FoldingCurve.size()/2;
       for(int i = mid+1; i < (int)FoldingCurve.size() - 1 && ind_l == (int)FoldingCurve.size() - 1; i++) ind_l = (calcCrossPt4Constraint(FoldingCurve[i], FoldingCurve[i-1]) != 0)? i: ind_l;
       for(int i = mid - 1; i > 0 && ind_r == 1; i--)  ind_r = (calcCrossPt4Constraint(FoldingCurve[i], FoldingCurve[i+1]) != 0)? i: ind_r;
       auto WriteCSV = [&](){
           //交点を動かしたときrulingの交差判定がどう変わるのか検証
           std::vector<std::shared_ptr<Vertex4d>> writeFC(ind_l - ind_r + 2), FC4csv;
           std::copy(&FoldingCurve[ind_r - 1], &FoldingCurve[ind_l + 1], writeFC.begin());//最適化を行う範囲(実際に最適化で変化させるのは両端の２つの交点)
           std::vector<double>X;
           for(int i = 0; i < (int)writeFC.size(); i++){
               std::shared_ptr<Vertex4d> v4d = writeFC[i]->deepCopy(); v4d->addline(writeFC[i]->line_parent); FC4csv.push_back(v4d);
               X.push_back((writeFC[i]->first->p - writeFC[i]->third->p).norm());
           }
           
           if(OptimizationAlgorithm == 0 && ind_l - ind_r + 2 == 4){
               std::ofstream ofs;
               std::string file = "./Optimization/OptimizationVertex.csv";
               ofs.open(file, std::ios::out); ofs << "t , s ,f_intersections,fruling,fconv, \n";
               for(double t = X[0] - 1; t <= X[0] + 1; t += 1e-2){
                   for(double s = X[1] - 2; s <= X[1] + 2; s += 1e-2){
                       FC4csv[0]->first->p  = t * (FC4csv[0]->first->p - FC4csv[0]->third->p).normalized() + FC4csv[0]->third->p;
                       FC4csv[0]->first->p3 = t * (FC4csv[0]->first->p3 - FC4csv[0]->third->p3).normalized() + FC4csv[0]->third->p3;
                       FC4csv[1]->first->p  = s * (FC4csv[1]->first->p - FC4csv[1]->third->p).normalized() + FC4csv[1]->third->p;
                       FC4csv[1]->first->p3 = s * (FC4csv[1]->first->p3 - FC4csv[1]->third->p3).normalized() + FC4csv[1]->third->p3;
                       cb_Folding(FC4csv, Poly_V, a_flap, StartingIndex);
                       double fr = f_intersections(FC4csv, StartingIndex);
                       double fr2 = calcCrossPt4Constraint(FC4csv[1],FC4csv[2]);
                       double k = RevisionVertices::getK(FC4csv[0]->first, FC4csv[1], FC4csv[2]->first),_k = RevisionVertices::getK(writeFC[0]->first, writeFC[1], writeFC[2]->first);
                       double tmp = (k - std::numbers::pi)*(_k - std::numbers::pi);
                       double fconv = (tmp < 0)? abs(tmp): 0;
                       ofs << t << " , " << s << " , " << fr <<" , " << fr2 << " , " << fconv << std::endl;
                   }   
               }
               ofs.close(); qDebug()<<"export csv file";
           }
       };
       //マルチスレッドのためにvoid関数でまとめておく
       std::mutex mtx_; // 排他制御用ミューテックス
       std::thread th_csv(WriteCSV);
       th_csv.join();
    }
    
    
    std::vector<Eigen::Vector3d> initX(FoldingCurve.size());
    for(int i = 0; i < (int)initX.size(); i++)initX[i] = FoldingCurve[i]->first->p;
    double fr;
    std::ofstream ofs;
    
    if(OptimizationAlgorithm == 0){//最適化による交点修正
       std::vector<std::shared_ptr<Vertex4d>> FC(FoldingCurve.size());
       do{
            std::vector<double> X, bnd_lower, bnd_upper;
            nlopt::opt opt;
            ind_l = FoldingCurve.size() - 1; ind_r = 0;//ind_l, ind_rがこの値から変わらない→片側でrulingの交差が起きていない
            for(int i = mid+1; i < (int)FoldingCurve.size() - 1 && ind_l == (int)FoldingCurve.size() - 1; i++) ind_l = (calcCrossPt4Constraint(FoldingCurve[i], FoldingCurve[i-1]) != 0)? i: ind_l;
            for(int i = mid - 1; i > 0 && ind_r == 0; i--)  ind_r = (calcCrossPt4Constraint(FoldingCurve[i], FoldingCurve[i+1]) != 0)? i: ind_r;

            qDebug()<<"ind_l = " << ind_l << " , ind_r = " << ind_r;
            for(int i = 0; i < (int)FoldingCurve.size(); i++){
                std::shared_ptr<Vertex4d> v4d = FoldingCurve[i]->deepCopy(); v4d->addline(FoldingCurve[i]->line_parent); FC[i] = v4d;
            }

            //境界条件と初期値を与える
            double up, down, x;
            if(ind_r != 0){
                    SetBounds(FoldingCurve, Poly_V, ind_r - 1, up, down, x); bnd_upper.push_back(up); bnd_lower.push_back(down); X.push_back(x);
            }else ind_r = -1;
            if(ind_l != (int)FoldingCurve.size()-1){
                SetBounds(FoldingCurve, Poly_V, ind_l + 1, up, down, x);  bnd_upper.push_back(up); bnd_lower.push_back(down); X.push_back(x);
            }else ind_l = -1;
            if(X.empty())return true;
            for(auto&x: X)qDebug()<<"before  x = " << x;

            RevisionVertices::ObjData_v od = {a_flap, FC, Poly_V, mid};
            od.SetCrossedIndex(ind_l, ind_r);
            od.AddBasePt(FoldingCurve);
            opt = nlopt::opt(nlopt::GN_ISRES, X.size());
            opt.add_inequality_constraint(Fconst_mixed, &od);
            opt.set_min_objective(ObjFunc_Vertex, &od);
            double ftol = 1e-9, xtol = 1e-9, maxtime = 30;
            SetOptimizationParameter(opt, bnd_lower, bnd_upper, ftol, xtol, maxtime);
            std::string file = "./Optimization/iteration_optimizer.csv";
            ofs.open(file, std::ios::out);
            double minf;
            try {
                    nlopt::result result = opt.optimize(X, minf);
                    qDebug() <<"result :  lower bound" <<result ;
            }catch (std::exception& e) { qDebug() << "nlopt failed: " << e.what() ; }
            if((int)X.size() == 2){
                    update_rightside(FoldingCurve, Poly_V, ind_r - 1, X[0]);
                    update_leftside(FoldingCurve, Poly_V, ind_l + 1, X[1]);
            }
            if((int)X.size() == 1){
                if(ind_r != -1)update_rightside(FoldingCurve, Poly_V, ind_r - 1, X[0]);
                else update_leftside(FoldingCurve, Poly_V, ind_l + 1, X[0]);
            }

            for(auto&x: X)qDebug()<<"after  x = " << x;

            for(int i = 0; i < (int)initX.size(); i++){
                double sgn = (MathTool::is_point_on_line(initX[i], FoldingCurve[i]->first->p, FoldingCurve[i]->third->p))?1: -1;
                qDebug() << i << " : error  = " << sgn * (initX[i] - FoldingCurve[i]->first->p).norm();
                ofs << sgn * (initX[i] - FoldingCurve[i]->first->p).norm() << " ,";
            }

            cb_Folding(FoldingCurve, Poly_V, a_flap, mid);
            fr = RulingsCrossed(FoldingCurve);
            ofs << "\n";
            if((bef_l != -1 && bef_r != -1) && (bef_l == ind_l && bef_r == ind_r))break;//同じ場所をloopしていたら修正を強制終了
            bef_r = ind_r; bef_l = ind_l;

       }while(fr != 0.0);

    }else{//逐次的な探索による修正
       bool res;
       std::string file = "./Optimization/iteration_";
       file += (OrderConstFunc == 0)? "FindFirstPoint.csv": "Fullsearch.csv";
       ofs.open(file, std::ios::out);
       do{
        ind_l = FoldingCurve.size() - 1; ind_r = 0;//ind_l, ind_rがこの値から変わらない→片側でrulingの交差が起きていない
        for(int i = mid+1; i < (int)FoldingCurve.size() - 1 && ind_l == (int)FoldingCurve.size() - 1; i++) ind_l = (calcCrossPt4Constraint(FoldingCurve[i], FoldingCurve[i-1]) != 0)? i: ind_l;
        for(int i = mid - 1; i > 0 && ind_r == 0; i--)  ind_r = (calcCrossPt4Constraint(FoldingCurve[i], FoldingCurve[i+1]) != 0)? i: ind_r;
        if(OrderConstFunc == 0)res = ReviseVertexPos(Poly_V, ind_l, ind_r, 0);//全探索ではないほう(最初に条件を満たす場所を見つける)
        else if(OrderConstFunc == 1)res = ReviseVertexPos(Poly_V, ind_l, ind_r, 1);//全探索
        cb_Folding(FoldingCurve, Poly_V, a_flap, mid);
        qDebug()<<"movement vertex ";

        for(int i = 0; i < (int)initX.size(); i++){
                double sgn = (MathTool::is_point_on_line(initX[i], FoldingCurve[i]->first->p, FoldingCurve[i]->third->p))?1: -1;
                qDebug() << i << " : error  = " << sgn * (initX[i] - FoldingCurve[i]->first->p).norm();
                ofs << sgn * (initX[i] - FoldingCurve[i]->first->p).norm() << " ,";
        }
        ofs << "\n";
        qDebug()<<"ind_l = " << ind_l   << " , ind_r = " << ind_r << " , result = " << res;
        fr = RulingsCrossed(FoldingCurve);
        if(!res || ((bef_l != -1 && bef_r != -1) && (bef_l == ind_l && bef_r == ind_r)))break;//同じ場所をloopしていたら修正を強制終了
        bef_l = ind_l; bef_r = ind_r;
       }while(fr != 0.0 && res);
    }
    
    ofs.close();
    return (fr == 0.0) ? true: false;
}

bool FoldLine::Optimization_Vertex(const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool IsStartEnd, int OrderConstFunc, int OptimizationAlgorithm){
    qDebug() <<"this optimization has to be started center of folding curve. End Point will not be mistaken result.";
    qDebug()<<"before optimization ";
    double fc_bef = RulingsCrossed(FoldingCurve);
    qDebug()<<"ruling intersection  : " << fc_bef;
    nlopt::opt opt;
    std::vector<double> X, bnd_lower, bnd_upper;
    std::vector<std::shared_ptr<Vertex4d>> FC;
    for(auto&fc: FoldingCurve){
        std::shared_ptr<Vertex4d> v4d = fc->deepCopy(); v4d->addline(fc->line_parent);
        FC.push_back(v4d);
        double init_l;
        if(fc == FoldingCurve.front() || fc == FoldingCurve.back()){
            init_l = (fc->first->p - fc->third->p).norm();
        }else{
            init_l = (v4d->first->p - v4d->third->p).norm();
            Eigen::Vector3d crossPoint;
            for(int k = 0; k < (int)Poly_V.size(); k++){
                crossPoint = MathTool::calcCrossPoint_2Vector(v4d->first->p, 1000.0 * (v4d->first->p - v4d->third->p).normalized() + v4d->first->p, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p);
                if(MathTool::is_point_on_line(crossPoint, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) && (v4d->first->p - v4d->third->p).normalized().dot(crossPoint - v4d->first->p) > 0)break;
            }
        }

        X.push_back(init_l);
        bnd_lower.push_back(init_l - 20);
        bnd_upper.push_back(init_l + 20);
        //bnd_lower.push_back(0);
        //bnd_upper.push_back((crossPoint - v4d->third->p).norm());
    }
    int StartingIndex = FoldingCurve.size()/2;
    RevisionVertices::ObjData_v od = {a_flap, FC, Poly_V, StartingIndex};
    od.AddBasePt(FoldingCurve);
    /*
    if(OptimizationAlgorithm == 0){
        qDebug()<<"optimization algorithm is MMA";
        opt = nlopt::opt(nlopt::LD_MMA, X.size());
    }else if(OptimizationAlgorithm == 1){
        qDebug()<<"optimization algorithm is CCSAQ";
        opt = nlopt::opt(nlopt::LD_CCSAQ, X.size());
    }else if(OptimizationAlgorithm == 2){
        qDebug()<<"optimization algorithm is SLSQP";
        opt = nlopt::opt(nlopt::LD_SLSQP, X.size());
    }else{
        qDebug()<<"not specified optimization algorithm";
        return false;
    }*/

    opt = nlopt::opt(nlopt::LD_SLSQP, X.size());
    opt.set_min_objective(ObjFunc_Vertex, &od);
    if(OrderConstFunc == 0){
        qDebug() << "order of constraint : convex -> ruling";
        opt.add_inequality_constraint(Fconst_conv, &od);
        opt.add_inequality_constraint(Fruling4Vertex, &od);
    }
    else if(OrderConstFunc == 1){
        qDebug() << "order of constraint : ruling -> convex";
        opt.add_inequality_constraint(Fruling4Vertex, &od);
        opt.add_inequality_constraint(Fconst_conv, &od);
    }else if(OrderConstFunc == 2){
        qDebug() << "order of constraint : mixed";
        //opt.add_inequality_constraint(Fconst_mixed, &od);

    }

    double ftol = 1e-9, xtol = 1e-9, maxtime = 100;
    SetOptimizationParameter(opt, bnd_lower, bnd_upper, ftol, xtol, maxtime);
    //opt.set_param("inner_maxeval", 10);
    double minf;
    try {
        nlopt::result result = opt.optimize(X, minf);
        qDebug()<<"result " << result;
    }catch (std::exception& e){qDebug() << "nlopt failed: " << e.what();}

    for(int i = 0; i < (int)X.size(); i++) qDebug() <<"norm error  " << i << "  ,  " << (FoldingCurve[i]->first->p - FoldingCurve[i]->third->p).norm() - X[i];

    for(int i = 0; i < (int)X.size(); i++){
        FoldingCurve[i]->first->p = X[i] * (FoldingCurve[i]->first->p - FoldingCurve[i]->third->p).normalized() + FoldingCurve[i]->third->p;
        FoldingCurve[i]->first->p3 = X[i] * (FoldingCurve[i]->first->p3 - FoldingCurve[i]->third->p3).normalized() + FoldingCurve[i]->third->p3;
    }
    cb_Folding(FoldingCurve, Poly_V, a_flap, StartingIndex);
    double fc = RulingsCrossed(FoldingCurve);
    qDebug() <<"before ruling" << fc_bef <<"  ,  result = " << fc<< " , error " << minf;


    return false;
}

bool FoldLine::Optimization_FlapAngle(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double wb, double wp, int rank, int alg, bool IsStartEnd, int OptimizationAlgorithm){

    auto apply_optimization = [&](nlopt::opt& opt, std::vector<double>& X, std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V,
                            int StartingIndex,  double& minf, double& fruling){
        try {
            qDebug() <<"small a" ;
            nlopt::result result = opt.optimize(X, minf);

            cb_Folding(FoldingCurve, Poly_V, X[0], StartingIndex);
            qDebug() <<"result :  lower bound" <<result ;
            fruling = RulingsCrossed(FoldingCurve);
            qDebug() << "found minimum at f(" << MathTool::rad2deg(X[0]) << ") = " << minf << "  , fruling = " << fruling;
        }catch (std::exception& e) { qDebug() << "nlopt failed: " << e.what() ; }
    };

    if(FoldingCurve.empty())return false;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    double a = 0, a_min, a_max;
    double ftol = 1e-9, xtol = 1e-9, maxtime = 1;
    const double th_ruling = 1e-9;

    Eigen::Vector3d e, e2, SpinAxis, Nb, N4;
    if(IsStartEnd){
        e = (ValidFC[0]->first->p3 - ValidFC[1]->first->p3).normalized(),e2 = (ValidFC[2]->first->p3 - ValidFC[1]->first->p3).normalized();
        SpinAxis = (ValidFC[1]->third->p3 - ValidFC[1]->first->p3).normalized();
        Nb = -(e.cross(e2)).normalized(), N4 = (e.cross(SpinAxis)).normalized();
    }else{
        int mid = ValidFC.size()/2;
        e = (ValidFC[mid-1]->first->p3 - ValidFC[mid]->first->p3).normalized(),e2 = (ValidFC[mid+1]->first->p3 - ValidFC[mid]->first->p3).normalized();
        SpinAxis = (ValidFC[mid]->third->p3 - ValidFC[mid]->first->p3).normalized();
        Nb = -(e.cross(e2)).normalized(), N4 = (e.cross(SpinAxis)).normalized();
    }

    double a_ll = std::atan2(e.dot(Nb.cross(N4)), N4.dot(-Nb));
    if(a_ll < 0)a_ll += 2.0*std::numbers::pi;
    double a_con = a_ll + std::numbers::pi;
    if(a_con > 2.0*std::numbers::pi)a_con -= 2.0*std::numbers::pi;
    if(a_con < 0)a_con +=2.0*std::numbers::pi;
    qDebug() << "a_con = " << MathTool::rad2deg(a_con);
    double phi3 = std::acos(e2.dot(SpinAxis)), phi4 = std::acos(e.dot(SpinAxis));
    double k = 2.0 * std::numbers::pi - phi3 - phi4;

    std::vector<double>bnd_lower,bnd_upper;
    Eigen::Vector3d FaceNp = (SpinAxis.cross(e2)).normalized(), CrossV = (N4.cross(FaceNp)).normalized();
    bool IsMount = (SpinAxis.dot(CrossV) < 0)? true: false;

    if(k < std::numbers::pi && IsMount){a_min = a_con + std::numbers::pi; a_max = 2.0 * std::numbers::pi;}
    if(k >= std::numbers::pi && IsMount){a_min = std::numbers::pi; a_max = std::numbers::pi + a_con;}
    if(k < std::numbers::pi && !IsMount){a_min = 0.0;a_max = a_con - std::numbers::pi;}
    if(k >= std::numbers::pi && !IsMount){a_min = a_con - std::numbers::pi; a_max = std::numbers::pi;}

    if(DebugMode::Singleton::getInstance().isdebug()) qDebug() <<MathTool::rad2deg(a_ll) << " , " << MathTool::rad2deg(a_con);
    qDebug() << "area " << MathTool::rad2deg(a_min) << " < " << MathTool::rad2deg(a) << " < " << MathTool::rad2deg(a_max) ;
    int StartingPoint = (IsStartEnd)? 0: FoldingCurve.size()/2;

    //csvファイルへの書き出し
    if(DebugMode::Singleton::getInstance().isdebug()){
        std::ofstream ofs2;
        std::filesystem::create_directory("./Optimization");
        std::string AngleFile;
        AngleFile = "./Optimization/OptimArea_" + std::to_string(rank) + "_" + std::to_string(validsize) + ".csv";
        ofs2.open(AngleFile, std::ios::out);
        ofs2 << "a(radian),a(degree) ,Eruling ,Ebend ,Eruling(N ooutline) ,Earea\n";
        for(double _a = a_min; _a <= a_max; _a+= 1e-3){
            cb_Folding(FoldingCurve, Poly_V, _a, StartingPoint);
            double f = RulingsCrossed(FoldingCurve);
            double fb = Fbend2(FoldingCurve);
            double fr = RulingsCrossed_NoOutline(FoldingCurve);
            double fa = Fmin_TriangleArea(FoldingCurve, false);
            ofs2 << _a << ", " << MathTool::rad2deg(_a) << " , " << f << ", " << fb << ", "<< fr  << ", " << fa << std::endl ;
        }ofs2.close();
    }

    RevisionVertices::ObjData od = {FoldingCurve, Poly_V, StartingPoint};
    double mu = 1e-1;
    od.AddWeight(wb, 0.1*wp, -1, mu);
    std::vector<double> X, minf, fruling;
    bnd_lower.push_back(a_min); bnd_upper.push_back(a_max);
    if(alg == 0){//ruling intersection as a constraint  
        od.AddBasePt(FoldingCurve);

        X.push_back(a_min + 0.5); minf.resize(2); fruling.resize(2);
        nlopt::opt opt;
        opt = nlopt::opt(nlopt::LD_MMA, 1);
        opt.set_min_objective(ObjFunc_RulingIntersection, &od);
        opt.add_inequality_constraint(Fruling, &od);
        SetOptimizationParameter(opt, bnd_lower, bnd_upper, ftol, xtol, maxtime);
        apply_optimization(opt, X, FoldingCurve, Poly_V, IsStartEnd, minf[0], fruling[0]);

        X[0] = a_max - 0.5;
        opt = nlopt::opt(nlopt::LD_MMA, 1);
        opt.set_min_objective(ObjFunc_RulingIntersection, &od);
        opt.add_inequality_constraint(Fruling, &od);
        SetOptimizationParameter(opt, bnd_lower, bnd_upper, ftol, xtol, maxtime);
        apply_optimization(opt, X, FoldingCurve, Poly_V, StartingPoint, minf[1], fruling[1]);

        //曲げエネルギーを使った最適化
        if(fruling[0] >= th_ruling && fruling[1] > th_ruling){
            if(minf[0] < minf[1]) a_flap = X[0];
            else a_flap = X[1];
            cb_Folding(FoldingCurve, Poly_V, a_flap, StartingPoint);
            qDebug()<<"could not avoid ruling cross";
            qDebug() << "result : smaller = " << MathTool::rad2deg(a_flap)  << "(" << a_flap << ") " <<  "  f_bend = " << Fbend2(FoldingCurve) << "  ,  f_ruling = " <<RulingsCrossed(FoldingCurve) ;
            return false;
        }
        if(fruling[0] < th_ruling && fruling[1] < th_ruling){
            if(minf[0] < minf[1]) a_flap = X[0];
            else a_flap = X[1];
        }else if(fruling[0] < th_ruling && fruling[0] >= th_ruling) a_flap = X[0];
        else if(fruling[0] >= th_ruling && fruling[1] < th_ruling)  a_flap = X[1];
        cb_Folding(FoldingCurve, Poly_V, a_flap, StartingPoint);
        qDebug() << "result : smaller = " << MathTool::rad2deg(a_flap)  << "(" << a_flap << "), tol = " << tol <<  "  f_bend = " << Fbend2(FoldingCurve) << "  ,  f_ruling = " <<RulingsCrossed(FoldingCurve) ;
        qDebug() << "finish";
        return true;

    }else if(alg == 1){//maximize triangle area
        od.AddBasePt(FoldingCurve);
        nlopt::opt opt;
        X.push_back((a_min + a_max)/2.0); minf.push_back(0); fruling.push_back(0);

        opt = nlopt::opt(nlopt::GN_ISRES, 1);
        opt.set_min_objective(ObjFunc_RegressionCurve, &od);
        opt.add_inequality_constraint(Fruling, &od);
        SetOptimizationParameter(opt, bnd_lower, bnd_upper, ftol, xtol, maxtime);

        apply_optimization(opt, X, FoldingCurve, Poly_V, StartingPoint, minf[0], fruling[0]);
        a_flap = X[0];
        cb_Folding(FoldingCurve, Poly_V, a_flap, StartingPoint);
        qDebug() << "result : " << MathTool::rad2deg(a_flap)  << "(" << a_flap << ") " <<  "  f_bend = " << Fbend2(FoldingCurve) << "  ,  f_ruling = " <<fruling[0] ;
        if(fruling[0] >= th_ruling){
            qDebug()<<"could not avoid ruling cross";
            return false;
        }
        qDebug() << "flap angle optimization is sucucessful";
        return true;
    }

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

void FoldLine::CheckIsCrossedRulings(){
    if((int)FoldingCurve.size() < 3)return;
    int mid = FoldingCurve.size()/2;
    for(int i = mid + 1; i < (int)FoldingCurve.size()-1; i++){
        if(calcCrossPt4Constraint(FoldingCurve[i], FoldingCurve[i-1]) != 0.0){
            for(int j = i; j < (int)FoldingCurve.size()-1; j++)FoldingCurve[j]->IsCalc = false;
            break;
        }
    }
    for(int i = mid - 1; i > 0; i--){
        if(calcCrossPt4Constraint(FoldingCurve[i], FoldingCurve[i+1]) != 0.0){
            for(int j = i; j < (int)FoldingCurve.size()-1; j++)FoldingCurve[j]->IsCalc = false;
            break;
        }
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
    AlignmentVertex4dDirection();

    parent->SortCurve();
    parent->AlignmentVertex4dDirection();
    //FoldingCurve.front()->first->p3 = (FoldingCurve.front()->second->p - FoldingCurve.front()->third->p).norm() * (FoldingCurve.front()->second->p3 - FoldingCurve.front()->third->p3).normalized() + FoldingCurve.front()->third->p3;
    //FoldingCurve.back()->first->p3 = (FoldingCurve.back()->second->p - FoldingCurve.back()->third->p).norm() * (FoldingCurve.back()->second->p3 - FoldingCurve.back()->third->p3).normalized() + FoldingCurve.back()->third->p3;
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

std::vector<std::vector<std::shared_ptr<Vertex>>> FoldLine::CalclateRegressionCurve(double a, const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool IsWriteCSV, bool IsStartEnd,
                                                                                    std::vector<std::vector<std::shared_ptr<Vertex>>>& Tri_fixside){
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

    int StartingIndex = (IsStartEnd)? 0: FoldingCurve.size()/2;
    std::vector<std::vector<std::shared_ptr<Vertex>>> Triangles;
    if(IsWriteCSV){
        std::ofstream ofs;
        std::filesystem::create_directory("./Optimization");
        std::string File = "./Optimization/RegressionCurve_";
        File += std::to_string(ValidFC.size() - 3) + ".csv" ;
        ofs.open(File, std::ios::out);
        ofs << "a(radian), a(degree), area, Ebend, Eruling"<<std::endl;
        for(double _a = 0.0; _a <= 2.0*std::numbers::pi; _a+= 1e-3){
            a = _a;
            cb_Folding(ValidFC, Poly_V, a, StartingIndex);
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
            ofs << _a << " , " << MathTool::rad2deg(_a) << " , " << area << " , " << Ebend << " , " << Eruling <<  std::endl;
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


double calclate_adjacent_angle(const std::shared_ptr<Vertex4d>&xbef, const std::shared_ptr<Vertex4d>&x, const std::shared_ptr<Vertex4d>& xnext, double a, double& phi1, Eigen::Vector3d& SrfN){
    Eigen::Vector3d e = (xbef->first->p3 - x->first->p3).normalized(), e2 = (xnext->first->p3 - x->first->p3).normalized() , Axis = (x->third->p3 - x->first->p3).normalized();
    double beta = std::acos(e.dot(e2)), k = 2.0 * std::numbers::pi - (std::acos(e2.dot(Axis)) + std::acos(e.dot(Axis)));
    double _phi1 = std::atan((std::cos(k) - std::cos(beta))/(std::sin(beta)*std::cos(a)- std::sin(k)));
    phi1 = (_phi1 < 0)? _phi1+std::numbers::pi: _phi1;
    double phi2 = k - phi1;
    double sin_a = sin(phi1)*sin(a)/sin(phi2);
    if(sin_a > 1)sin_a = 1; else if(sin_a < -1)sin_a = -1;
    double cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
    if(cos_a > 1)cos_a = 1; else if(cos_a < -1)cos_a = -1;
    double a2 = std::atan2(sin_a, cos_a);
    if(a2 < 0)a2 += 2*std::numbers::pi;//koko chuui
    SrfN = MathTool::ProjectionVector(e.cross(e2), e2, true);
    return a2;
}

//引数として渡されるFoldingCurveは間引かれたrulingがないものと仮定して実装する
double GetFlapAngle(int ind, const std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, double a_init, int StartingIndex){

    double a = a_init , phi1, a2;
    Eigen::Vector3d e, e2, befN, SrfN;

    a = (a_init < 0.0)? a_init + 2.0 * std::numbers::pi: (a_init > 2.0*std::numbers::pi)? a_init - 2.0*std::numbers::pi: a;
    if(ind == StartingIndex)return a;
    else if(ind > StartingIndex){
        for(int i = StartingIndex; i < (int)FoldingCurve.size() - 1; i++){
            e = (FoldingCurve[i-1]->first->p3 - FoldingCurve[i]->first->p3).normalized();
            e2 = (FoldingCurve[i+1]->first->p3 - FoldingCurve[i]->first->p3).normalized();
            if(i != StartingIndex){
                SrfN = MathTool::ProjectionVector(e.cross(e2), -e, true);
                a = update_flapangle(a2, befN, SrfN, e);
            }
            a2 = calclate_adjacent_angle(FoldingCurve[i-1], FoldingCurve[i], FoldingCurve[i+1], a, phi1, SrfN);
            if(i == (int)FoldingCurve.size() - 2 || i == ind)break;
            befN = SrfN;
        }
    }else{
        //右手系の回転
        Eigen::Vector3d e_cmn = (FoldingCurve[StartingIndex]->first->p3 - FoldingCurve[StartingIndex - 1]->first->p3).normalized();//共有する辺
        Eigen::Vector3d e_prev = (FoldingCurve[StartingIndex + 1]->first->p3 - FoldingCurve[StartingIndex]->first->p3).normalized();//左手系で使った辺
        Eigen::Vector3d e_next = (FoldingCurve[StartingIndex - 1]->first->p3 - FoldingCurve[StartingIndex - 2]->first->p3).normalized();//右手系で使う辺
        a = 2*std::numbers::pi - a;
        SrfN = MathTool::ProjectionVector(-e_next.cross(e_cmn), e_cmn, true);
        befN = MathTool::ProjectionVector(-e_cmn.cross(e_prev), -e_cmn, true);
        double dir = e_cmn.dot(befN.cross(SrfN)), tau = std::acos(SrfN.dot(befN));
        if(dir < 0)tau = 2*std::numbers::pi - tau;
        a = 2.0*std::numbers::pi  - (a + tau);
        if(a > 2*std::numbers::pi)a -= 2*std::numbers::pi; else if(a < 0)a += 2*std::numbers::pi;
        befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
        for(int i = StartingIndex - 1; i > 0; i--){
            e = (FoldingCurve[i + 1]->first->p3 - FoldingCurve[i]->first->p3).normalized();
            e2 = (FoldingCurve[i - 1]->first->p3 - FoldingCurve[i]->first->p3).normalized();
            SrfN = MathTool::ProjectionVector(e.cross(e2), -e, true);
            if(i != StartingIndex - 1){
                dir = e.dot(befN.cross(SrfN));
                tau = std::acos(SrfN.dot(befN));
                if(dir <= 0)tau = 2.0*std::numbers::pi - tau;
                if(a2 - tau <= 0)a= a2 - tau + 2.0 * std::numbers::pi;  else a = a2 - tau;
            }
            a2 = calclate_adjacent_angle(FoldingCurve[i + 1], FoldingCurve[i], FoldingCurve[i - 1], a, phi1, SrfN);
            if(i == 1 || i == ind)break;
            befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
        }
    }
    return a;
}

void FoldLine::applyAAAMethod(const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool IsStartEnd, double a){

    if(FoldingCurve.empty() && a_flap == -1)return;

    //SimplifyModel(validsize, isroot);
    //qDebug()<<"tempolary remove simplifyModel method from applyAAAMethod for test of optimization vertex.";
    //std::string msg = (!IsStartEnd)? "center": "end point";
    int StartingIndex = (IsStartEnd)? 0: FoldingCurve.size()/2;
    cb_Folding(FoldingCurve, Poly_V, a, StartingIndex);
    if(RulingsCrossed(FoldingCurve) < 1e-9){
        a_flap = a; //tol = _tol;
    }
    a_flap = a;
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

    double phi1;
    a2 = calclate_adjacent_angle(xbef, x, xnext, a, phi1, SrfN);
}

void CalcRuling_back(double a, std::shared_ptr<Vertex4d>& xbef, std::shared_ptr<Vertex4d>& x, std::shared_ptr<Vertex4d>& xnext, const std::vector<std::shared_ptr<Vertex>>& Poly_V,  double& a2, Eigen::Vector3d& SrfN, const Eigen::Vector3d &SpinAxis ){

    double phi1;
    a2 = calclate_adjacent_angle(xbef, x, xnext, a, phi1, SrfN);
    Eigen::Vector3d e = (xbef->first->p3 - x->first->p3).normalized(), e2 = (xnext->first->p3 - x->first->p3).normalized();
    Eigen::Vector3d N = (Eigen::AngleAxisd(a, e) * e.cross(e2)).normalized(); // 回転行列を作成
    Eigen::Vector3d r3d = (Eigen::AngleAxisd(phi1, N) * e).normalized();
    Eigen::Vector3d r2d = (Eigen::AngleAxisd(phi1, SpinAxis)* (xbef->first->p- x->first->p)).normalized();//展開図のruling方向
    Eigen::Vector3d crossPoint;
    for(int k = 0; k < (int)Poly_V.size(); k++){
        crossPoint = MathTool::calcCrossPoint_2Vector(x->first->p, 1000.0 * r2d + x->first->p, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p);
        if(MathTool::is_point_on_line(crossPoint, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) && r2d.dot(crossPoint - x->first->p) > 0)break;
    }
    x->second->p = crossPoint;
    x->second->p3 = (crossPoint - x->first->p).norm() * r3d + x->first->p3;
}

void _FoldingAAAMethod_right(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle, int StartingIndex){
   std::vector<std::shared_ptr<Vertex4d>> FC;
   for(auto&fc:FoldingCurve){if(fc->IsCalc)FC.push_back(fc);}
   Eigen::Vector3d e, e2, SrfN, befN;
   double tau, dir, a2;
   double a = (angle < 0.0)? angle + 2.0 * std::numbers::pi: (angle > 2.0*std::numbers::pi)? angle - 2.0*std::numbers::pi: angle;
   for(int ind = StartingIndex -1; ind > 0; ind--){
       e = (FC[ind + 1]->first->p3 - FC[ind]->first->p3).normalized();
       e2 = (FC[ind - 1]->first->p3 - FC[ind]->first->p3).normalized();
       SrfN = MathTool::ProjectionVector(e.cross(e2), -e, true);
       if(ind != StartingIndex-1){
           dir = e.dot(befN.cross(SrfN));
           tau = std::acos(SrfN.dot(befN));
           if(dir <= 0)tau = 2.0*std::numbers::pi - tau;
           if(a2 - tau <= 0)a= a2 - tau + 2.0 * std::numbers::pi;  else a = a2 - tau;
       }
       CalcRuling_back(a, FC[ind+1], FC[ind], FC[ind-1], Poly_V, a2, befN, Eigen::Vector3d::UnitZ());

       if(ind == 1)break;
       befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
   }
}

void _FoldingAAAMethod_center(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle, int StartingIndex){

   double a2 = 0;
    Eigen::Vector3d e, e2;
    Eigen::Vector3d N, r, befN, befedge;
    std::vector<std::shared_ptr<Vertex4d>> FC;
    for(auto&fc:FoldingCurve){if(fc->IsCalc)FC.push_back(fc);}

    double a = (angle < 0.0)? angle + 2.0 * std::numbers::pi: (angle > 2.0*std::numbers::pi)? angle - 2.0*std::numbers::pi: angle;
    double a_init = a;
    std::shared_ptr<Vertex4d> fc, fc_bef, fc_next;

    for(int ind = StartingIndex; ind < (int)FC.size() - 1; ind++){
        e = (FC[ind-1]->first->p3 - FC[ind]->first->p3).normalized();
        e2 = (FC[ind+1]->first->p3 - FC[ind]->first->p3).normalized();
        if(ind != StartingIndex){
            Eigen::Vector3d SrfN = MathTool::ProjectionVector(e.cross(e2), -e, true);
            a = update_flapangle(a2, befN, SrfN, e);
        }
        CalcRuling(a, FC[ind-1], FC[ind], FC[ind+1], Poly_V, a2, befN, -Eigen::Vector3d::UnitZ());
        if(ind == (int)FoldingCurve.size() -2)break;
        befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
    }
    if(StartingIndex == 0)return;
    //右手系の回転
    Eigen::Vector3d e_cmn = (FC[StartingIndex]->first->p3 - FC[StartingIndex-1]->first->p3).normalized();//共有する辺
    Eigen::Vector3d e_prev = (FC[StartingIndex+1]->first->p3 - FC[StartingIndex]->first->p3).normalized();//左手系で使った辺
    Eigen::Vector3d e_next = (FC[StartingIndex-1]->first->p3 - FC[StartingIndex-2]->first->p3).normalized();//右手系で使う辺
    a = (angle < 0.0)? angle + 2.0 * std::numbers::pi: (angle > 2.0*std::numbers::pi)? angle - 2.0*std::numbers::pi: angle;
    a = 2*std::numbers::pi - a;
    Eigen::Vector3d SrfN = MathTool::ProjectionVector(-e_next.cross(e_cmn), e_cmn, true);
    befN = MathTool::ProjectionVector(-e_cmn.cross(e_prev), -e_cmn, true);
    double dir = e_cmn.dot(befN.cross(SrfN));
    double tau = std::acos(SrfN.dot(befN));
    if(dir < 0)tau = 2*std::numbers::pi - tau;
    a = 2.0*std::numbers::pi  - (a + tau);
    if(a > 2*std::numbers::pi)a -= 2*std::numbers::pi;
    if(a < 0)a += 2*std::numbers::pi;
    //a = 2.0*std::numbers::pi -  a;
    befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
    for(int ind = StartingIndex -1; ind > 0; ind--){
        e = (FC[ind + 1]->first->p3 - FC[ind]->first->p3).normalized();
        e2 = (FC[ind - 1]->first->p3 - FC[ind]->first->p3).normalized();
        SrfN = MathTool::ProjectionVector(e.cross(e2), -e, true);
        if(ind != StartingIndex-1){
            dir = e.dot(befN.cross(SrfN));
            tau = std::acos(SrfN.dot(befN));
            if(dir <= 0)tau = 2.0*std::numbers::pi - tau;
            if(a2 - tau <= 0)a= a2 - tau + 2.0 * std::numbers::pi;  else a = a2 - tau;
        }
        CalcRuling_back(a, FC[ind+1], FC[ind], FC[ind-1], Poly_V, a2, befN, Eigen::Vector3d::UnitZ());

        if(ind == 1)break;
        befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
    }
}

void _FoldingAAAMethod_left(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle){

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
        std::shared_ptr<Vertex> v_clst = getClosestVertex(FoldingCurve[Vertices_Ind[0]]->second , FoldingCurve[Vertices_Ind[0]]->first, FoldingCurve, true);
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
        std::shared_ptr<Vertex> v_clst = getClosestVertex(FoldingCurve[Vertices_Ind.back()]->second , FoldingCurve[Vertices_Ind.back()]->first, FoldingCurve, true);
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

