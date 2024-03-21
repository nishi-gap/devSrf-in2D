#include "foldline.h"

std::ofstream ofs_obj, ofs_sbj;
bool hasContainSbjfunc = false;

const double eps = 1e-7;
inline Eigen::Vector3d _calcruling3d(const double& a, Eigen::Vector3d e, Eigen::Vector3d e2, Eigen::Vector3d axis, double& beta, std::vector<double>& Phi);
void CalcRuling(double a, std::shared_ptr<Vertex4d>& xbef, std::shared_ptr<Vertex4d>& x, std::shared_ptr<Vertex4d>& xnext, const std::vector<std::shared_ptr<Vertex>>& Poly_V, double& a2, Eigen::Vector3d& SrfN, const Eigen::Vector3d& SpinAixs);
void CalcRuling_back(double a, std::shared_ptr<Vertex4d>& xbef, std::shared_ptr<Vertex4d>& x, std::shared_ptr<Vertex4d>& xnext, const std::vector<std::shared_ptr<Vertex>>& Poly_V, double& a2, Eigen::Vector3d& SrfN, const Eigen::Vector3d& SpinAixs);
inline double update_flapangle(double a, const Eigen::Vector3d& befN, const Eigen::Vector3d& SrfN, const Eigen::Vector3d& e);
void _FoldingAAAMethod_center(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle, int StartingIndex);

void cb_Folding(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve, const std::vector<std::shared_ptr<Vertex>>& Poly_V, const double angle, int StartingPoint){
    _FoldingAAAMethod_center(FoldingCurve, Poly_V, angle, StartingPoint);
}

inline Eigen::Vector3d calcCrossPoint_2Vertex(const std::shared_ptr<Vertex>& p1, const std::shared_ptr<Vertex>& q1, const std::shared_ptr<Vertex>& p2, const std::shared_ptr<Vertex>& q2){
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
    return;
    tar->first->p = t * (src->line_parent->v->p - src->line_parent->o->p).normalized() + src->line_parent->o->p;
    tar->first->p3 = t * (src->line_parent->v->p3 - src->line_parent->o->p3).normalized() + src->line_parent->o->p3;
    return;

}

namespace RevisionVertices{
    using FoldLine3d = std::vector<std::shared_ptr<Vertex4d>>;

    struct OptimizeParam{
        FoldLine3d FC;
        std::vector<std::shared_ptr<Vertex>> Poly_V;
        FoldLine3d BasePt;
        int StartingIndex;
        double wsim, wp;
        void AddWeight(double _wsim, double _wp){ wsim =_wsim; wp = _wp;}      
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
    
    using ObjData = OptimizeParam;
    using ObjData_v = OptimizeParam_v;

    Eigen::Vector3d decideRulingDirectionOn3d(Eigen::Vector3d e, Eigen::Vector3d N, double a, double phi);

    double getK(const std::shared_ptr<Vertex>&xbef, const std::shared_ptr<Vertex4d>&x ,const std::shared_ptr<Vertex>&xnext){
        Eigen::Vector3d e = (xbef->p - x->first->p).normalized(), e2 = (xnext->p - x->first->p).normalized(), axis = (x->third->p - x->first->p).normalized();
        return 2*std::numbers::pi - std::acos(e.dot(axis)) - std::acos(e2.dot(axis));
    }
}

FoldLine::FoldLine(PaintTool _type)
{
    CtrlPts.clear();
    curveNum = 300;
    type = _type;
    CurvePts.clear();
    a_flap = -1;
    validsize = 0;
}

FoldLine::FoldLine(const FoldLine& _FL){
    CtrlPts = _FL.CtrlPts;
    curveNum = _FL.curveNum;
    type = _FL.type;
    CurvePts = _FL.CurvePts;
    a_flap = _FL.a_flap;
    validsize = _FL.validsize;
    FoldingCurve.clear();
    for(const auto& v4d: _FL.FoldingCurve)FoldingCurve.push_back(v4d);
}

std::shared_ptr<FoldLine> FoldLine::deepCopy(){
    std::shared_ptr<FoldLine> fl = std::make_shared<FoldLine>(type);
    fl->CtrlPts = CtrlPts;  fl->CurvePts = CurvePts; fl->a_flap = a_flap; fl->validsize = validsize;
    fl->FoldingCurve.resize(FoldingCurve.size());
    for(int i = 0; i < (int)FoldingCurve.size(); i++)fl->FoldingCurve[i] = FoldingCurve[i]->deepCopy();

    return fl;
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

void FoldLine::initialize_foldstate(const std::vector<std::shared_ptr<Vertex>>& Poly_V){
    if((int)FoldingCurve.size() < 3)return;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    int StartingPoint = FoldingCurve.size()/2;
    Eigen::Vector3d e, e2, SpinAxis, Nb, N4;
    int mid = ValidFC.size()/2;
    e = (ValidFC[mid-1]->first->p3 - ValidFC[mid]->first->p3).normalized(),e2 = (ValidFC[mid+1]->first->p3 - ValidFC[mid]->first->p3).normalized();
    SpinAxis = (ValidFC[mid]->third->p3 - ValidFC[mid]->first->p3).normalized();
    Nb = -(e.cross(e2)).normalized(), N4 = (e.cross(SpinAxis)).normalized();

    double a_ll = std::atan2(e.dot(Nb.cross(N4)), N4.dot(-Nb));
    if(a_ll < 0)a_ll += 2.0*std::numbers::pi;
    double a_con = a_ll + std::numbers::pi;
    if(a_con > 2.0*std::numbers::pi)a_con -= 2.0*std::numbers::pi;
    if(a_con < 0)a_con +=2.0*std::numbers::pi;
    qDebug() << "a_con = " << MathTool::rad2deg(a_con);
    cb_Folding(FoldingCurve, Poly_V, a_con, StartingPoint);
}

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
    setCurve(dim);
}

bool FoldLine::moveCtrlPt(Eigen::Vector3d& p, int movePtIndex, int dim){
    if((int)CtrlPts.size() < movePtIndex || movePtIndex < 0)return false;
    CtrlPts[movePtIndex] = p;
    setCurve(dim);
    if(FoldingCurve.empty())return false;
    return true;

}

std::shared_ptr<FoldLine> FoldLine::addCtrlPt(Eigen::Vector3d& p, int dim){
    //if((int)CtrlPts.size() <= dim)
    CtrlPts.push_back(p);
    bool hasCrv = setCurve(dim);
    return shared_from_this();
}

//type == 0 line
//type == 1 b-spline
bool FoldLine::setCurve(int dim){
    using namespace MathTool;
    CurvePts.clear();
    if(type == PaintTool::Crease){
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

double RulingsCrossed(std::vector<std::shared_ptr<Vertex4d>>& FoldingCurve){
    double f = 0.0;
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    int mid = ValidFC.size()/2;

    auto f_cross = [](const std::shared_ptr<Vertex4d>& Vmove, const std::shared_ptr<Vertex4d>& Vfix){
        Eigen::Matrix2d A; Eigen::Vector2d b;
        A(0,0) = Vmove->second->p.x() - Vmove->first->p.x(); A(0,1) = -(Vfix->second->p.x() - Vfix->first->p.x());
        A(1,0) = Vmove->second->p.y() - Vmove->first->p.y(); A(1,1) = -(Vfix->second->p.y() - Vfix->first->p.y());
        b(0) = Vfix->first->p.x() - Vmove->first->p.x(); b(1) = Vfix->first->p.y() - Vmove->first->p.y();
        Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
        return (ts(1) >= 1 || ts(1) <= 0)? 0.0: 1.0/std::abs(ts(1));
    };

    for(int i = mid; i > 0; i--) f += f_cross(ValidFC[i-1], ValidFC[i]);
    for(int i = mid; i < (int)ValidFC.size() -1; i++) f += f_cross(ValidFC[i+1], ValidFC[i]);
    return f;
}

Eigen::Vector3d RevisionVertices::decideRulingDirectionOn3d(Eigen::Vector3d e, Eigen::Vector3d N, double a, double phi){
    e = e.normalized();
    N = (Eigen::AngleAxisd(a, -e) * N).normalized(); // 回転行列を作成
    return (Eigen::AngleAxisd(phi, N) * e).normalized();
}

bool FoldLine::isbend(){
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    double fr = RulingsCrossed(ValidFC);
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

    //if(DebugMode::Singleton::getInstance().isdebug())qDebug() <<"constraint function = " << MathTool::rad2deg(a[0]) << "(" << a[0] << ")  , " << f ;
    return f;
 }

inline double calcCrossPt4Constraint(const std::shared_ptr<Vertex4d>& x, const std::shared_ptr<Vertex4d>& x2){
    Eigen::Matrix2d A; Eigen::Vector2d b;
    A(0,0) = x->first->p.x() - x->second->p.x(); A(0,1) = -(x2->first->p.x() - x2->second->p.x());
    A(1,0) = x->first->p.y() - x->second->p.y(); A(1,1) = -(x2->first->p.y() - x2->second->p.y());
    b(0) = x->first->p.x() - x2->first->p.x(); b(1) = x->first->p.y() - x2->first->p.y();
    Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
    return ((ts(0) >= 1 && ts(1) >= 1) || (ts(0) <= 0 && ts(1) <= 0))? 0.0: 1.0/std::min(std::abs(ts(0)), std::abs(ts(1)));
    return ((0 < ts(0) && ts(0) < 1) && ((0 < ts(1) && ts(1) < 1))) ? 1.0/std::min(ts(0), ts(1)): (!(0 < ts(0) && ts(0) < 1) && ((0 < ts(1) && ts(1) < 1)))? 1.0/ts(1): ((0 < ts(0) && ts(0) < 1) && !((0 < ts(1) && ts(1) < 1)))? 1.0/ts(0):  0;
}


double Fconv(const std::vector<std::shared_ptr<Vertex4d>>& FC, const std::vector<std::shared_ptr<Vertex4d>>& BasePt, int s, int e){
    double f = 0.0;
    for(int i = s; i <= e; i++){
        double kori = RevisionVertices::getK(BasePt[i-1]->first, BasePt[i], BasePt[i+1]->first);
        double k = RevisionVertices::getK(FC[i-1]->first, FC[i], FC[i+1]->first);
        double tmp = (k - std::numbers::pi)*(kori - std::numbers::pi);
        if(tmp < 0)f += abs(tmp);
    }
    return f;
}

double Fconst_mixed(const std::vector<double> &X, std::vector<double> &grad, void* f_data){
    RevisionVertices::ObjData_v *od = (RevisionVertices::ObjData_v *)f_data;

    int mid = od->FC.size()/2;
    auto func = [&mid](double x, std::vector<std::shared_ptr<Vertex4d>>& FC, std::vector<std::shared_ptr<Vertex4d>>& BasePt,
                                      std::vector<std::shared_ptr<Vertex>>& Poly_V, double a, int n, int s)->double{
        auto fruling = [](const std::shared_ptr<Vertex4d>& Vmove, const std::shared_ptr<Vertex4d>& Vfix){
            Eigen::Matrix2d A; Eigen::Vector2d b;
            A(0,0) = Vmove->second->p.x() - Vmove->first->p.x(); A(0,1) = -(Vfix->second->p.x() - Vfix->first->p.x());
            A(1,0) = Vmove->second->p.y() - Vmove->first->p.y(); A(1,1) = -(Vfix->second->p.y() - Vfix->first->p.y());
            b(0) = Vfix->first->p.x() - Vmove->first->p.x(); b(1) = Vfix->first->p.y() - Vmove->first->p.y();
            Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
            return (0.0 <= ts(1) && ts(1) < 1.0)? 1.0/ts(1): 0.0;
        };

        double fc = 0.0, fr = 0.0;
        update_Vertex(x, FC[n], BasePt[n]);
        cb_Folding(FC, Poly_V, a, s);
        if(n > mid){
            if(n <= (int)FC.size() - 1){
                fr += fruling(FC[n - 1], FC[n - 2]); fc += Fconv(FC, BasePt, s, n - 1);
            }
        }
        else{
            if(n >= 0){
                fr += fruling(FC[n + 1], FC[n + 2]); fc += Fconv(FC, BasePt, n + 1, s);
            }
        }
        return fr + fc;
    };

    double f = func(X[0], od->FC, od->BasePt, od->Poly_V, od->a, od->ind_l, od->StartingIndex);
    if(!grad.empty()){
        double fp = func(X[0] + eps, od->FC, od->BasePt, od->Poly_V, od->a, od->ind_l, od->StartingIndex);
        double fm = func(X[0] - eps, od->FC, od->BasePt, od->Poly_V, od->a, od->ind_l, od->StartingIndex);
        grad[0] = (fp - fm)/(2.0*eps);
    }
    ofs_sbj << X[0] << "," << f << "\n";
    return f;
}

double ObjFunc_Vertex(const std::vector<double> &X, std::vector<double> &grad, void* f_data){

    RevisionVertices::ObjData_v *od = (RevisionVertices::ObjData_v *)f_data;
    std::vector<std::shared_ptr<Vertex>> Poly_V = od->Poly_V;
    int mid = od->FC.size()/2;

    auto objfunc = [&mid](double x, std::vector<std::shared_ptr<Vertex4d>>& FC, std::vector<std::shared_ptr<Vertex4d>>& BasePt,
                          std::vector<std::shared_ptr<Vertex>>& Poly_V, double wsim, double wp, double a, int n, int s){
        update_Vertex(x, FC[n], BasePt[n]);
        cb_Folding(FC, Poly_V, a, s);
        double fa = 0.0, fsim = 0.0;
        auto fparalell = [](std::shared_ptr<Vertex4d>& x1, std::shared_ptr<Vertex4d>& x2){
            Eigen::Vector3d v1 = (x1->second->p - x1->first->p).normalized(), v2 = (x2->second->p - x2->first->p).normalized();
            return 1.0 - v1.dot(v2);
        };
        if(n > mid){
            if(n <= (int)FC.size() - 1)fa += wp * fparalell(FC[n - 2], FC[n - 1]);
        }
        else{
            if(n >= 0)fa += wp * fparalell(FC[n +1], FC[n + 2]);
        }
        fsim = wsim* (FC[n]->first->p - BasePt[n]->first->p).norm();

        return fa + fsim;
    };

    auto constfunc = [&mid](double x, std::vector<std::shared_ptr<Vertex4d>>& FC, std::vector<std::shared_ptr<Vertex4d>>& BasePt,
                       std::vector<std::shared_ptr<Vertex>>& Poly_V, double a, int n, int s)->double{
        auto fruling = [](const std::shared_ptr<Vertex4d>& Vmove, const std::shared_ptr<Vertex4d>& Vfix){
            Eigen::Matrix2d A; Eigen::Vector2d b;
            A(0,0) = Vmove->second->p.x() - Vmove->first->p.x(); A(0,1) = -(Vfix->second->p.x() - Vfix->first->p.x());
            A(1,0) = Vmove->second->p.y() - Vmove->first->p.y(); A(1,1) = -(Vfix->second->p.y() - Vfix->first->p.y());
            b(0) = Vfix->first->p.x() - Vmove->first->p.x(); b(1) = Vfix->first->p.y() - Vmove->first->p.y();
            Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
            return (0.0 <= ts(1) && ts(1) < 1.0)? 1.0/ts(1): 0.0;
        };

        double fc = 0.0, fr = 0.0;
        update_Vertex(x, FC[n], BasePt[n]);
        cb_Folding(FC, Poly_V, a, s);
        if(n > mid){
            if(n <= (int)FC.size() - 1){
                fr += fruling(FC[n - 1], FC[n - 2]); fc += Fconv(FC, BasePt, s, n - 1);
            }
        }
        else{
            if(n >= 0){
                fr += fruling(FC[n + 1], FC[n + 2]); fc += Fconv(FC, BasePt, n + 1, s);
            }
        }
        return fr + fc;
    };

    //geometory
    double wc = 1, wo = 1e-2;
    double fc = wc * constfunc(X[0], od->FC, od->BasePt, od->Poly_V, od->a, od->ind_l, od->StartingIndex);
    double fo = wo * objfunc(X[0], od->FC, od->BasePt, od->Poly_V, od->wsim, od->wp, od->a, od->ind_l, od->StartingIndex);
    if(!grad.empty()){
       double fop = wo * objfunc(X[0] + eps, od->FC, od->BasePt, od->Poly_V, od->wsim, od->wp, od->a, od->ind_l, od->StartingIndex);
       double fom = wo * objfunc(X[0] - eps, od->FC, od->BasePt, od->Poly_V, od->wsim, od->wp, od->a, od->ind_l, od->StartingIndex);
       grad[0] = (fop - fom)/(2.0*eps);
    }
    double f = (hasContainSbjfunc)? fo + fc: fo;
    ofs_obj << X[0] << "," << f << "\n";

    update_Vertex(X[0], od->FC[od->ind_r], od->BasePt[od->ind_r]);
    //qDebug() <<"objfunc = " << fo;
    return  f;
}

void SetOptimizationParameter(nlopt::opt& opt, const std::vector<double>&bnd_lower, const std::vector<double>&bnd_upper, double ftol, double xtol, double maxtime){
    opt.set_lower_bounds(bnd_lower);
    opt.set_upper_bounds(bnd_upper);
    opt.set_xtol_rel(xtol);
    //opt.set_ftol_rel(ftol);
    opt.set_maxtime(maxtime);//stop over this time
}

#include <chrono>
bool FoldLine::PropagateOptimization_Vertex(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double wp, double wsim){

    if(isbend())return true;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC, FC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    int bef_l = -1, bef_r = -1;
    int mid = ValidFC.size()/2;
    int ind_l = ValidFC.size() - 1, ind_r = 1;

    std::vector<Eigen::Vector3d> initX(ValidFC.size());
    for(int i = 0; i < (int)initX.size(); i++)initX[i] = ValidFC[i]->first->p;
    double fr;
    std::ofstream ofs;

    auto SetParam = [&](const std::shared_ptr<Vertex4d>& v4d, double&x, double&lower, double& upper){
        x = (v4d->first->p - v4d->third->p).norm();
        lower = std::min((v4d->first->p - v4d->third->p).norm(), 0.0);
        Eigen::Vector3d point;
        for(int k = 0; k < (int)Poly_V.size(); k++){
            point = MathTool::calcCrossPoint_2Vector(v4d->first->p, 1000.0 * (v4d->first->p - v4d->third->p).normalized() + v4d->first->p, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p);
            if(MathTool::is_point_on_line(point, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) && (v4d->first->p - v4d->third->p).normalized().dot(point - v4d->first->p) > 0)break;
        }
        upper = std::max((v4d->first->p - v4d->third->p).norm(), (point - v4d->third->p).norm());
    };

    auto apply_optimization = [&](nlopt::opt& opt, int ind)->double{
        if(ind < 0 || ind >= (int)ValidFC.size())return -1;
        int mid = ValidFC.size()/2;

        double ftol = 1e-5, xtol = 1e-5, maxtime = 2;

        std::vector<double> X(1), bnd_lower(1), bnd_upper(1);//ind_lに左右関係なくパラメータを定義→左右の交点位置の最適化は別々に行う
        double minf = 1e+10, minX = (ValidFC[ind]->first->p - ValidFC[ind]->third->p).norm();
        //qDebug()<<"init X = " << minX;
        SetParam(ValidFC[ind], X[0], bnd_lower[0], bnd_upper[0]);

        SetOptimizationParameter(opt, bnd_lower, bnd_upper, ftol, xtol, maxtime);
        double _minf;
        nlopt::result result;
        try {
            result = opt.optimize(X, _minf);
            qDebug() <<"Result is " << result;
        }catch (std::exception& e) {
            qDebug() << "nlopt failed: " << e.what();
        }

        update_Vertex(X[0], FC[ind], ValidFC[ind]);
        cb_Folding(FC, Poly_V, a_flap, mid);
        if(ind > mid){//左側
            double fr = calcCrossPt4Constraint(FC[ind - 2], FC[ind - 1]);
            if(fr == 0 && _minf <= minf){
                minX = X[0];
                minf = _minf;
            }
        }else {//右側
            double fr = calcCrossPt4Constraint(FC[ind + 2], FC[ind + 1]);
            if((fr == 0) && _minf <= minf){
                minX = X[0];
                minf = _minf;
            }
        }
        //qDebug()<<"ind " << ind << " , minf = " << minf << " , X = " << minX;
        update_Vertex(minX, FC[ind], ValidFC[ind]);
        cb_Folding(FC, Poly_V, a_flap, mid);
        return minX;
    };

    ind_l = std::min(mid+2, (int)ValidFC.size()-1); ind_r = std::max(mid-2, 0);
    FC.resize(ValidFC.size());

    //geometory
    nlopt::opt opt = nlopt::opt(nlopt::LN_AUGLAG, 1);
    RevisionVertices::ObjData_v od = {a_flap, FC, Poly_V, mid};// Specify the augmented Lagrangian algorithm
    od.AddBasePt(ValidFC);
    od.AddWeight(wsim, wp);

    if(!hasContainSbjfunc){
       nlopt::opt local_optimizer(nlopt::LN_NELDERMEAD, 1);// Specify the Nelder-Mead algorithm for solving each step
       opt.set_local_optimizer(local_optimizer);
       opt.add_equality_constraint(Fconst_mixed, &od);
    }else opt = nlopt::opt(nlopt::LN_NELDERMEAD, 1);
    opt.set_min_objective(ObjFunc_Vertex, &od);

    qDebug()<<"hasContainSbjfunc is " << hasContainSbjfunc;

    for(int i = 0; i < (int)FC.size(); i++){
       FC[i] = ValidFC[i]->deepCopy(); FC[i]->addline(ValidFC[i]->line_parent);
    }
    od.FC = FC;

    auto starttime = std::chrono::high_resolution_clock::now();
    do{

       od.SetCrossedIndex(ind_l, ind_l);
       std::string filepath_obj = "./Optimization/objfunc" + std::to_string(ind_l) + ".csv";
       std::string filepath_sbj = "./Optimization/sbjfunc" + std::to_string(ind_l) + ".csv";
       ofs_obj.open(filepath_obj, std::ios::out); ofs_sbj.open(filepath_sbj, std::ios::out);
       apply_optimization(opt, ind_l);
       ofs_obj.close(); ofs_sbj.close();

       od.SetCrossedIndex(ind_r, ind_r);
       filepath_obj = "./Optimization/objfunc" + std::to_string(ind_r) + ".csv";
       filepath_sbj = "./Optimization/sbjfunc" + std::to_string(ind_r) + ".csv";
       ofs_obj.open(filepath_obj, std::ios::out); ofs_sbj.open(filepath_sbj, std::ios::out);
       apply_optimization(opt, ind_r);
       ofs_obj.close(); ofs_sbj.close();

       qDebug()<<"ind_l = " << ind_l   << " , ind_r = " << ind_r;
       fr = RulingsCrossed(FC);

       if(((bef_l != -1 && bef_r != -1) && (bef_l == ind_l && bef_r == ind_r)))break;//同じ場所をloopしていたら修正を強制終了
       bef_l = ind_l; bef_r = ind_r;

       if(ind_l < (int)ValidFC.size())ind_l++;
       if(ind_r >= 0)ind_r--;
       if(ind_r < 0 && ind_l >= (int)ValidFC.size())break;

    }while(fr != 0.0);

    for(int i = 0; i < (int)initX.size(); i++){
       double sgn = (MathTool::is_point_on_line(initX[i], ValidFC[i]->first->p, ValidFC[i]->third->p))?1: -1;
       qDebug() << i << " : error  = " << sgn * (initX[i] - ValidFC[i]->first->p).norm();
       ofs << sgn * (initX[i] - ValidFC[i]->first->p).norm() << " ,";
    }
    auto endtime = std::chrono::high_resolution_clock::now();
    qDebug() << "Elapsed time: " <<  std::chrono::duration_cast<std::chrono::milliseconds>(endtime - starttime).count() << " millisecond\n";

    for(int i = 0; i < (int)FC.size(); i++){
       ValidFC[i]->first->p = FC[i]->first->p; ValidFC[i]->first->p3 = FC[i]->first->p3;
    }
    cb_Folding(ValidFC, Poly_V, a_flap, (int)ValidFC.size()/2);
    ofs.close();
    //hasContainSbjfunc = !hasContainSbjfunc;
    return (fr == 0.0) ? true: false;
}


bool FoldLine::Optimization_FlapAngle(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double wp, double wsim, int rank){

    auto apply_optimization = [&](nlopt::opt& opt, std::vector<double>& X, std::vector<std::shared_ptr<Vertex4d>>& ValidFC, const std::vector<std::shared_ptr<Vertex>>& Poly_V,
                            int StartingIndex,  double& minf, double& fruling){
        try {
            qDebug() <<"small a" ;
            nlopt::result result = opt.optimize(X, minf);

            cb_Folding(ValidFC, Poly_V, X[0], StartingIndex);
            qDebug() <<"result :  lower bound" <<result ;
            fruling = RulingsCrossed(ValidFC);
            qDebug() << "found minimum at f(" << MathTool::rad2deg(X[0]) << ") = " << minf << "  , fruling = " << fruling;
        }catch (std::exception& e) { qDebug() << "nlopt failed: " << e.what() ; }
    };

    if(FoldingCurve.empty())return false;

    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&fc: FoldingCurve){if(fc->IsCalc)ValidFC.push_back(fc);}
    double a = 0, a_min, a_max;
    double ftol = 1e-6, xtol = 1e-6, maxtime = 1;
    const double th_ruling = 1e-9;

    Eigen::Vector3d e, e2, SpinAxis, Nb, N4;
    int mid = ValidFC.size()/2;
    e = (ValidFC[mid-1]->first->p3 - ValidFC[mid]->first->p3).normalized(),e2 = (ValidFC[mid+1]->first->p3 - ValidFC[mid]->first->p3).normalized();
    SpinAxis = (ValidFC[mid]->third->p3 - ValidFC[mid]->first->p3).normalized();
    Nb = -(e.cross(e2)).normalized(), N4 = (e.cross(SpinAxis)).normalized();

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
    int StartingPoint = ValidFC.size()/2;

    //csvファイルへの書き出し
    if(DebugMode::Singleton::getInstance().isdebug()){
        std::ofstream ofs2;
        std::filesystem::create_directory("./Optimization");
        std::string AngleFile;
        AngleFile = "./Optimization/OptimArea_" + std::to_string(rank) + "_" + std::to_string(validsize) + ".csv";
        ofs2.open(AngleFile, std::ios::out);
        ofs2 << "a(radian),a(degree) ,Eruling ,Eruling(is in range) ,Earea(same weight), Earea(different weight)\n";
        for(double _a = 0; _a <= 2*std::numbers::pi; _a+= 1e-3){
            cb_Folding(ValidFC, Poly_V, _a, StartingPoint);
            double f = RulingsCrossed(ValidFC);
            std::string fr = (a_min <= _a && _a <= a_max)? std::to_string(RulingsCrossed(ValidFC)): "";
            ofs2 << _a << ", " << MathTool::rad2deg(_a) << " , " << f << ", "<< fr  << ", " << std::endl ;
        }ofs2.close();
    }

    RevisionVertices::ObjData od = {ValidFC, Poly_V, StartingPoint};

    od.AddWeight(wsim, wp);
    std::vector<double> X, minf, fruling;
    bnd_lower.push_back(a_min); bnd_upper.push_back(a_max);
    X.push_back(a_min + 0.5); minf.resize(2); fruling.resize(2);
    nlopt::opt opt = nlopt::opt(nlopt::LN_NELDERMEAD, 1);
    opt.set_min_objective(Fruling, &od);
    SetOptimizationParameter(opt, bnd_lower, bnd_upper, ftol, xtol, maxtime);
    apply_optimization(opt, X, ValidFC, Poly_V, StartingPoint, minf[0], fruling[0]);
    double minx = X[0];

    X[0] = a_max - 0.5;
    apply_optimization(opt, X, ValidFC, Poly_V, StartingPoint, minf[1], fruling[1]);

    if(fruling[0] >= th_ruling && fruling[1] > th_ruling){
        if(minf[0] < minf[1]) a_flap = minx;
        else a_flap = X[0];
        cb_Folding(ValidFC, Poly_V, a_flap, StartingPoint);
        qDebug()<<"could not avoid ruling cross";
        qDebug() << "result : smaller = " << MathTool::rad2deg(a_flap)  << "(" << a_flap << ") " << "  ,  f_ruling = " <<RulingsCrossed(FoldingCurve) ;
        return false;
    }
    if(fruling[0] < th_ruling && fruling[1] < th_ruling){
        if(minf[0] < minf[1]) a_flap = minx;
        else a_flap = X[0];
    }else if(fruling[0] < th_ruling && fruling[0] >= th_ruling) a_flap = minx;
    else if(fruling[0] >= th_ruling && fruling[1] < th_ruling)  a_flap = X[0];
    cb_Folding(ValidFC, Poly_V, a_flap, StartingPoint);
    qDebug() << "result : smaller = " << MathTool::rad2deg(a_flap)  << "(" << a_flap << "),  f_ruling = " <<RulingsCrossed(FoldingCurve) ;
    qDebug() << "finish";
    return true;
}

void FoldLine::CheckIsCrossedRulings(){
    if((int)FoldingCurve.size() < 3)return;
    int mid = FoldingCurve.size()/2;
    for(int i = mid + 1; i < (int)FoldingCurve.size()-1; i++){//左側
        if(calcCrossPt4Constraint(FoldingCurve[i], FoldingCurve[i-1]) != 0.0){
            for(int j = i; j < (int)FoldingCurve.size()-1; j++)FoldingCurve[j]->IsCalc = false;
            break;
        }
    }
    for(int i = mid - 1; i > 0; i--){//右側
        if(calcCrossPt4Constraint(FoldingCurve[i], FoldingCurve[i+1]) != 0.0){
            for(int j = i; j > 0; j--)FoldingCurve[j]->IsCalc = false;
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
                if(IS.type_mvk == 1 || IS.type_mvk == 2){
                    IS.v->line_parent->color = -IS.v->line_parent->color;
                }
            }else{
                if(IS.type_mvk == 0 || IS.type_mvk == 3){
                    IS.v->line_parent->color = -IS.v->line_parent->color;
                }
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

    if(parent == nullptr)return;
    if(parent->FoldingCurve.empty())return;
    int dim = 3;
    FoldingCurve.clear();
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

    if(FoldingCurve.empty()){
        for(auto&l: polyline_surface){
            std::shared_ptr<CrvPt_FL> p;
            if(type == PaintTool::Crease)p = getCrossPoint(CtrlPts, l->o, l->v, dim);
            if(p == nullptr)continue;
            FoldingCurve.push_back(std::make_shared<Vertex4d>(Vertex4d(p, l->v, l->o)));
        }
        SortCurve();
    }else{
        for(auto&l: polyline_surface){
            std::shared_ptr<CrvPt_FL> p;
            if(type == PaintTool::Crease)p = getCrossPoint(CtrlPts, l->o, l->v, dim);
            if(p == nullptr)continue;
            if(std::abs(p->s - FoldingCurve.front()->first->s) < std::abs(p->s - FoldingCurve.back()->first->s)){
                FoldingCurve.front()->third = l->o; FoldingCurve.front()->second = l->v; FoldingCurve.front()->first = p;
            }else{
                FoldingCurve.back()->third = l->o; FoldingCurve.back()->second = l->v; FoldingCurve.back()->first = p;
            }
        }
    }

    for(auto it = parent->FoldingCurve.begin() + 1; it != parent->FoldingCurve.end() - 1; it++){
        if(!(*it)->IsCalc)continue;
        std::shared_ptr<CrvPt_FL> p;
        if(type == PaintTool::Crease)p = getCrossPoint(CtrlPts, (*it)->first, (*it)->second, dim);
        if(p != nullptr){
            std::shared_ptr<Vertex4d> v4d = std::make_shared<Vertex4d>(p, (*it)->second, (*it)->first);
            v4d->addline((*it));
            FoldingCurve.push_back(v4d);
            (*it)->second = p;
        }
        if(type == PaintTool::Crease)p = getCrossPoint(CtrlPts, (*it)->first, (*it)->third, dim);
        if(p != nullptr){
            std::shared_ptr<Vertex4d> v4d = std::make_shared<Vertex4d>(p, (*it)->first, (*it)->third);
            v4d->addline((*it));
            FoldingCurve.push_back(v4d);
            (*it)->third = p;
        }
    }
    validsize = FoldingCurve.size();
    SortCurve();
    AlignmentVertex4dDirection();

    parent->SortCurve();
    parent->AlignmentVertex4dDirection();
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

void FoldLine::applyAAAMethod(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double a){

    if(FoldingCurve.empty() && a_flap == -1)return;
    int StartingIndex = FoldingCurve.size()/2;
    cb_Folding(FoldingCurve, Poly_V, a, StartingIndex);
    if(RulingsCrossed(FoldingCurve) < 1e-9){
        a_flap = a; //tol = _tol;
    }

    a_flap = a;
    return;
}

inline double update_flapangle(double a, const Eigen::Vector3d& befN, const Eigen::Vector3d& SrfN, const Eigen::Vector3d& e){
    double dir = (-e).dot(befN.cross(SrfN));
    double tau = (SrfN.dot(befN) < -1)? std::numbers::pi: (SrfN.dot(befN) > 1)? 0: std::acos(SrfN.dot(befN));
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

void CalcRuling_back(double a, std::shared_ptr<Vertex4d>& xbef, std::shared_ptr<Vertex4d>& x, std::shared_ptr<Vertex4d>& xnext, const std::vector<std::shared_ptr<Vertex>>& Poly_V,  double& a2, Eigen::Vector3d& SrfN, const Eigen::Vector3d &SpinAxis){
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
    double tau = (SrfN.dot(befN) < -1)? std::numbers::pi: (SrfN.dot(befN) > 1)? 0: std::acos(SrfN.dot(befN));
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
            tau = (SrfN.dot(befN) < -1)? std::numbers::pi: (SrfN.dot(befN) > 1)? 0: std::acos(SrfN.dot(befN));
            if(dir <= 0)tau = 2.0*std::numbers::pi - tau;
            if(a2 - tau <= 0)a= a2 - tau + 2.0 * std::numbers::pi;  else a = a2 - tau;
        }
        CalcRuling_back(a, FC[ind+1], FC[ind], FC[ind-1], Poly_V, a2, befN, Eigen::Vector3d::UnitZ());

        if(ind == 1)break;
        befN = MathTool::ProjectionVector(e.cross(e2), e2, true);
    }
}
