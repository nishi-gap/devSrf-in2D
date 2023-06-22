
#include "foldline.h"

const double eps = 1e-7;
using namespace MathTool;

std::string File_Ebend = "./Optimization/Ebend.csv";
std::string File_Eruling = "./Optimization/Eruling.csv";
std::ofstream ofs_Ebend, ofs_Eruling;

namespace RevisionVertices{
    using FoldLine3d = std::vector<Vertex4d>;
    struct OptimizeParam{
        FoldLine3d FC;
        std::vector<HalfEdge*> Edges;
        std::vector<Vertex*> Vertices, Poly_V;
        bool IsValidSigmoid;
        double phi02, phim1;
        std::vector<double*> res_Fbend, res_Fruling, res_a;
        int type;
        OptimizeParam(FoldLine3d& _FC, std::vector<HalfEdge*>& _Edges, std::vector<Vertex*>& _Vertices, std::vector<Vertex*>& _Poly_V, double _phi02, double _phim1, int _t, bool IVS):
            FC{_FC}, Edges{_Edges}, Vertices{_Vertices}, Poly_V{_Poly_V}, phi02{_phi02}, phim1{_phim1}, type{_t}, IsValidSigmoid{IVS}{}
        ~OptimizeParam(){}
    };
    struct OptimizeParam_v: public OptimizeParam{
        double a;
        OptimizeParam_v(double _a, FoldLine3d& _FC, std::vector<HalfEdge*>& _Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, double _phi02, double _phi01, int _t, bool IVS):
            a{_a}, OptimizeParam::OptimizeParam( _FC, _Edges, Vertices,  Poly_V, _phi02, _phi01, _t, IVS){}
        ~OptimizeParam_v(){}
    };

    using ObjData = OptimizeParam;
    using ObjData_v = OptimizeParam_v;

    inline double getK(const glm::f64vec3 o, const glm::f64vec3 x, const glm::f64vec3 x2);
    inline glm::f64vec3 getVertex(const glm::f64vec3& o, const glm::f64vec3& p, const double t);
    double F_vpos(FoldLine3d &FC, const std::vector<double>& T, std::vector<double>& grad);
    double F_k(FoldLine3d &FC, const std::vector<double>& T, std::vector<double>& grad);
    double F_ruling(const std::vector<double> &T, std::vector<double> &grad, void* f_data);
    double Const_Edev(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    double E_fair(const std::vector<glm::f64vec3>& X, FoldLine3d& FC);
    double E_sim(const std::vector<glm::f64vec3>& X, FoldLine3d& FC);
    double Minimize_SmoothSrf(const std::vector<double>& X, std::vector<double> &grad, void *f_data);
    double Minimize_Vpos(const std::vector<double> &T, std::vector<double> &grad, void* f_data);
}

void _FoldingAAAMethod(std::vector<Vertex4d>& FoldingCurve, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Poly_V,double a, double phi02, double phim1, int type);
Vertex* findAxisVertex(Vertex4d& R, Vertex *e1, Vertex *e2);
std::vector<Vertex4d> TrimPoints2(std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces, std::vector<Vertex*>& Vertices, std::vector<Vertex4d>& FoldingCurve, double tol);
bool IsRulingCrossed(glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint,  std::vector<Vertex*>& Poly_V);
std::vector<CrvPt_FL*> SetPointsOnCurve(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_v, std::vector<glm::f64vec3>& CtrlPts, int dim, int divSize);
inline glm::f64vec3 calcCrossPoint_2Vector(Vertex* p1, Vertex* q1, Vertex* p2, Vertex* q2);
inline glm::f64vec3 calcTargetDistanceOnPlane(Vertex *p, Vertex *o, Vertex *v1, Vertex *v2);
void Douglas_Peucker_algorithm(std::vector<Vertex4d>& FoldingCurve, std::vector<Vertex4d>& res, double tol = std::numbers::pi/9.0);
void Douglas_Peucker_algorithm2(std::vector<Vertex4d>& FoldingCurve, std::vector<Vertex4d>& res, int size);
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

double Fbend(std::vector<Vertex*>& Poly_V,std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices){
    auto Cot_jk = [](glm::f64vec3 ej, glm::f64vec3 ek){
        double a = glm::angle(ej, ek);
        return 1.0/std::tan(a);
    };

    double f = 0.0;
    //三角形分割
    //すべてのエッジに対して曲げエネルギーの計算
    //相和を返す
    /*
    std::vector<HalfEdge*> _edges = EdgeCopy(Edges, Vertices);
    std::vector<Face*> _faces; EdgeRecconection(Poly_V,_faces, _edges);

    for(const auto&e: _edges){
        if(e->pair == nullptr)continue;
        glm::f64vec3 N = e->face->getNormalVec();
        glm::f64vec3 Np = e->pair->face->getNormalVec();
        double phi = (abs(glm::dot(N,Np)) > 1)? std::numbers::pi: std::numbers::pi - std::acos(glm::dot(N,Np));
        f += 1.0/(phi*phi);
    }*/
    return f;
}
double Fbend2(std::vector<Vertex4d>& FoldingCurve, bool IsValidSigmoid){
    std::vector<int> Vertices_Ind;
    auto Sigmoid =[](double x, bool IsValidSigmoid){
        if(IsValidSigmoid)return 1.0/(1.0 + std::exp(-x));
        return x;
    };
    double f = 0.0;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}
    glm::f64vec3 Nt = glm::normalize(glm::cross(FoldingCurve[0].first->p3 - FoldingCurve[Vertices_Ind[1]].first->p3, FoldingCurve[Vertices_Ind[1]].second->p3 - FoldingCurve[Vertices_Ind[1]].first->p3));
    glm::f64vec3 Nb = glm::normalize(glm::cross(FoldingCurve[Vertices_Ind[1]].third->p3 - FoldingCurve[Vertices_Ind[1]].first->p3, FoldingCurve[0].first->p3 - FoldingCurve[Vertices_Ind[1]].first->p3));
    double phi = ((glm::dot(Nt,Nb)) > 1)? std::numbers::pi: ((glm::dot(Nt,Nb)) < -1)? 0:  std::numbers::pi - abs(std::acos(glm::dot(Nt,Nb)));
    f += Sigmoid(1.0/(phi*phi), IsValidSigmoid);
    for(int i = 1; i < (int)Vertices_Ind.size() - 1; i++){
        glm::f64vec3 Ntp = glm::normalize(glm::cross(FoldingCurve[Vertices_Ind[i]].first->p3 - FoldingCurve[Vertices_Ind[i+1]].first->p3,
                FoldingCurve[Vertices_Ind[i+1]].second->p3 - FoldingCurve[Vertices_Ind[i+1]].first->p3));
        glm::f64vec3 Nbp = glm::normalize(glm::cross(FoldingCurve[Vertices_Ind[i+1]].third->p3 - FoldingCurve[Vertices_Ind[i+1]].first->p3,
                FoldingCurve[Vertices_Ind[i]].first->p3 - FoldingCurve[Vertices_Ind[i+1]].first->p3));

        phi = ((glm::dot(Ntp,Nbp)) > 1)? std::numbers::pi: ((glm::dot(Ntp,Nbp)) < -1)? 0:  std::numbers::pi - abs(std::acos(glm::dot(Ntp,Nbp)));
        f += Sigmoid(1.0/(phi*phi), IsValidSigmoid);
        phi = ((glm::dot(Ntp,Nt)) > 1)? std::numbers::pi: ((glm::dot(Ntp,Nt)) < -1)? 0:  std::numbers::pi - std::acos(glm::dot(Ntp,Nt));
        f += Sigmoid(1.0/(phi*phi), IsValidSigmoid);
        Nt = Ntp; Nb = Nbp;
    }
    return f;
}

double Fruling(const std::vector<double> &a, std::vector<double> &grad, void* f_data)
{
    auto f_r =[](std::vector<Vertex4d>& FoldingCurve)->double{
        double f = 0.0;
        Eigen::Matrix2d A;
        Eigen::Vector2d b;
        std::vector<int> Vertices_Ind;
        for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}
        for(int i = 1; i < (int)Vertices_Ind.size() -1; i++){
            for(int j = -1; j <= 1; j +=2){
                int i2 = i + j;
                if(i2 <= 0 || i2 >= (int)Vertices_Ind.size() - 1)continue;
                A(0,0) = FoldingCurve[Vertices_Ind[i]].second->p.x - FoldingCurve[Vertices_Ind[i]].first->p.x;
                A(0,1) = -(FoldingCurve[Vertices_Ind[i2]].second->p.x - FoldingCurve[Vertices_Ind[i2]].first->p.x);
                A(1,0) = FoldingCurve[Vertices_Ind[i]].second->p.y - FoldingCurve[Vertices_Ind[i]].first->p.y;
                A(1,1) = -(FoldingCurve[Vertices_Ind[i2]].second->p.y - FoldingCurve[Vertices_Ind[i2]].first->p.y);
                b(0) = FoldingCurve[Vertices_Ind[i2]].first->p.x - FoldingCurve[Vertices_Ind[i]].first->p.x;
                b(1) = FoldingCurve[Vertices_Ind[i2]].first->p.y - FoldingCurve[Vertices_Ind[i]].first->p.y;
                Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
                if(0 < ts(0) && ts(0) < 1) f += 1.0/ts(0);
            }
        }
        return f;
    };
    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    double phi02 = od->phi02, phim1 = od->phim1;
    int type = od->type;
    RevisionVertices::FoldLine3d FoldingCurve = od->FC;
    double f = f_r(FoldingCurve);
    if(!grad.empty()){
        std::vector<HalfEdge*> Edges = od->Edges;
        std::vector<Vertex*> Poly_V = od->Poly_V;
        std::vector<Vertex> tmp;
        for(auto v: FoldingCurve)tmp.push_back(v.second);
        _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a[0] + eps, phi02, phim1, type);
       double fp = f_r(FoldingCurve);
        _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a[0] - eps, phi02, phim1, type);
       double fm = f_r(FoldingCurve);
       grad[0] = (fp - fm)/(2.0 * eps);
       //for(int j = 0; j < (int)FoldingCurve.size(); j++)FoldingCurve[j].second->p3 = tmp[j].p3;
    }
    if(ofs_Eruling.is_open()) ofs_Eruling << a[0] <<", " << glm::degrees(a[0]) << ", " <<  f <<  " ,  " << grad[0] << std::endl;
    return f;
 }

double ObjFunc(const std::vector<double> &a, std::vector<double> &grad, void* f_data){
    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    double phi02 = od->phi02, phim1 = od->phim1;
    int type = od->type;
    RevisionVertices::FoldLine3d FoldingCurve = od->FC;
    std::vector<HalfEdge*> Edges = od->Edges;
    std::vector<Vertex*>Poly_V = od->Poly_V;
    _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a[0], phi02, phim1, type);
    double f = Fbend2(FoldingCurve, od->IsValidSigmoid);
    if(!grad.empty()){
        _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a[0] + eps, phi02, phim1, type);
       double fp = Fbend2(FoldingCurve, od->IsValidSigmoid);
        _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a[0] - eps, phi02, phim1, type);
       double fm = Fbend2(FoldingCurve, od->IsValidSigmoid);
       grad[0] = (fp - fm)/(2.0 * eps);
    }  
    //std::cout <<"Fbend(" << glm::degrees(a[0]) << ")  =  " <<  f <<  " ,  " << grad[0] << std::endl;
    if(ofs_Ebend.is_open())ofs_Ebend << a[0] <<", " << glm::degrees(a[0]) << ", " <<  f <<  " ,  " << grad[0] << std::endl;
    return f;
}

double RevisionVertices::Const_Edev(const std::vector<double>& X, std::vector<double> &grad, void *f_data){

}
double RevisionVertices::E_fair(const std::vector<glm::f64vec3>& X, FoldLine3d& FC){

}
double RevisionVertices::E_sim(const std::vector<glm::f64vec3>& X, FoldLine3d& FC){
    int n = 0;
    double e = 0.0;
    for(auto&fc: FC){
        if(!fc.IsCalc){
            e += glm::distance(X[n], fc.first->p3_ori);
            n++;
        }
    }
    return e;
}

double RevisionVertices::Minimize_SmoothSrf(const std::vector<double>& X, std::vector<double> &grad, void *f_data){
    RevisionVertices::ObjData *od = (RevisionVertices::ObjData *)f_data;
    FoldLine3d FC = od->FC;
}

bool FoldLine::Optimization_SmooothSrf(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, int type){
    std::vector<double> X;
    for(auto&FC : FoldingCurve){
        if(!FC.IsCalc){X.push_back(FC.first->p3_ori.x); X.push_back(FC.first->p3_ori.y); X.push_back(FC.first->p3_ori.z);}
    }
    double phi02, phim1;
    RevisionVertices::ObjData od = {FoldingCurve, Edges, Vertices, Poly_V, phi02, phim1, type, false};
    nlopt::opt opt;
    opt = nlopt::opt(nlopt::LD_MMA, X.size());
    opt.set_min_objective(RevisionVertices::Minimize_SmoothSrf, &od);
    opt.add_inequality_constraint(Fruling, &od);
    opt.add_inequality_constraint(RevisionVertices::Const_Edev, &od);
    //opt.set_param("inner_maxeval", 100);
    opt.set_maxtime(2.0);//stop over this time
    opt.set_xtol_rel(1e-13);
}

bool FoldLine::Optimization_FlapAngle(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, int type, bool IsValidSigmoid){
    if(FoldingCurve.empty())return false;

    double phi02, phim1;
    std::vector<int> Vertices_Ind;for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}
    phi02 = glm::angle(glm::normalize(FoldingCurve[Vertices_Ind[0]].second->p - FoldingCurve[Vertices_Ind[0]].first->p),
            glm::normalize(FoldingCurve[Vertices_Ind[1]].first->p - FoldingCurve[Vertices_Ind[0]].first->p));
    phim1 = glm::angle(glm::normalize(FoldingCurve[Vertices_Ind.end()[-2]].first->p - FoldingCurve.back().first->p),
            glm::normalize(FoldingCurve.back().second->p - FoldingCurve.back().first->p));
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
         _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a, phi02, phim1, type);
         double f = 0.0;
         Eigen::Matrix2d A;
         Eigen::Vector2d b;
         for(int i = 1; i < (int)Vertices_Ind.size() -1; i++){
             for(int j = -1; j <= 1; j +=2){
                 int i2 = i + j;
                 if(i2 <= 0 || i2 >= (int)Vertices_Ind.size() - 1)continue;
                 A(0,0) = FoldingCurve[Vertices_Ind[i]].second->p.x - FoldingCurve[Vertices_Ind[i]].first->p.x;
                 A(0,1) = -(FoldingCurve[Vertices_Ind[i2]].second->p.x - FoldingCurve[Vertices_Ind[i2]].first->p.x);
                 A(1,0) = FoldingCurve[Vertices_Ind[i]].second->p.y - FoldingCurve[Vertices_Ind[i]].first->p.y;
                 A(1,1) = -(FoldingCurve[Vertices_Ind[i2]].second->p.y - FoldingCurve[Vertices_Ind[i2]].first->p.y);
                 b(0) = FoldingCurve[Vertices_Ind[i2]].first->p.x - FoldingCurve[Vertices_Ind[i]].first->p.x;
                 b(1) = FoldingCurve[Vertices_Ind[i2]].first->p.y - FoldingCurve[Vertices_Ind[i]].first->p.y;
                 Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
                 if(0 < ts(0) && ts(0) < 1) f += 1.0/ts(0);
             }
         }
         double fb = Fbend2(FoldingCurve, false);
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
    RevisionVertices::ObjData od = {FoldingCurve, Edges, Vertices, Poly_V, phi02, phim1, type, false};
    nlopt::opt opt;
    opt = nlopt::opt(nlopt::LD_MMA, 1);
    opt.set_min_objective(ObjFunc, &od);
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
         _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, _a, phi02, phim1, type);
         double f = 0.0;
         Eigen::Matrix2d A;
         Eigen::Vector2d b;
         for(int i = 1; i < (int)Vertices_Ind.size() -1; i++){
             for(int j = -1; j <= 1; j +=2){
                 int i2 = i + j;
                 if(i2 <= 0 || i2 >= (int)Vertices_Ind.size() - 1)continue;
                 A(0,0) = FoldingCurve[Vertices_Ind[i]].second->p.x - FoldingCurve[Vertices_Ind[i]].first->p.x;
                 A(0,1) = -(FoldingCurve[Vertices_Ind[i2]].second->p.x - FoldingCurve[Vertices_Ind[i2]].first->p.x);
                 A(1,0) = FoldingCurve[Vertices_Ind[i]].second->p.y - FoldingCurve[Vertices_Ind[i]].first->p.y;
                 A(1,1) = -(FoldingCurve[Vertices_Ind[i2]].second->p.y - FoldingCurve[Vertices_Ind[i2]].first->p.y);
                 b(0) = FoldingCurve[Vertices_Ind[i2]].first->p.x - FoldingCurve[Vertices_Ind[i]].first->p.x;
                 b(1) = FoldingCurve[Vertices_Ind[i2]].first->p.y - FoldingCurve[Vertices_Ind[i]].first->p.y;
                 Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);
                 if(0 < ts(0) && ts(0) < 1) f += 1.0/ts(0);
             }
         }
         double fb = Fbend2(FoldingCurve, false);
         double val = (f == 0)? fb: 1;
        ofs2 << _a << ", " << glm::degrees(_a) << " , " << f << ", " << fb << ", "<< val << std::endl;
    }ofs2.close();

    double minf_amin, minf_amax;
    double res_amin, res_amax;
    ofs_Ebend.open(File_Ebend, std::ios::out); ofs_Eruling.open(File_Eruling, std::ios::out);
      try {
        std::vector<double> _a{a_min + 0.1};
          nlopt::result result = opt.optimize(_a, minf_amin);

          std::cout <<"result :  lower bound" <<result << std::endl;
          std::cout << "found minimum at f(" << glm::degrees(_a[0]) << ") = " << std::setprecision(10) << minf_amin << std::endl;
          res_amin = _a[0];
      }
      catch (std::exception& e) {
          std::cout << "nlopt failed: " << e.what() << std::endl;
      }ofs_Ebend.close(); ofs_Eruling.close();
    try {
      std::vector<double> _a{a_max - 0.1};
        nlopt::result result = opt.optimize(_a, minf_amax);

        std::cout <<"result :  upper bound" <<result << std::endl;
        std::cout << "found minimum at f(" << glm::degrees(_a[0]) << ") = " << std::setprecision(10) << minf_amax << std::endl;
        res_amax = _a[0];
    }
    catch (std::exception& e) {
        std::cout << "nlopt failed: " << e.what() << std::endl;
    }
   if(minf_amax > minf_amin){std::cout << "result : smaller = " << glm::degrees(res_amin) << ",  f = " << minf_amin << std::endl; a2 = res_amin; }
   else{std::cout << "result : smaller = " << glm::degrees(res_amax) << ",  f = " << minf_amax << std::endl; a2 = res_amax; }
     _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, a2, phi02, phim1, type);
    std::cout << "finish"<<std::endl;
    return true;
}

void FoldLine::Optimization_Vertices(std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, std::vector<Vertex*>& Poly_V, int type, bool IsValidSigmoid){
    if(FoldingCurve.empty())return;
    glm::f64vec3 Nb = glm::normalize(glm::cross(FoldingCurve[2].first->p3 - FoldingCurve[1].first->p3, FoldingCurve[0].first->p3 - FoldingCurve[1].first->p3));
    glm::f64vec3 Nf = glm::normalize(glm::cross(FoldingCurve[1].first->p3 - FoldingCurve[0].first->p3, FoldingCurve[1].second->p3 - FoldingCurve[1].first->p3));
    glm::f64vec3 SpinAxis = glm::normalize(FoldingCurve[0].first->p3 - FoldingCurve[1].first->p3);
    double a = (glm::dot(Nb, Nf) < -1)? std::numbers::pi: (glm::dot(Nb, Nf) > 1)? 0: std::acos(glm::dot(Nb, Nf));
    if(glm::dot(glm::cross(Nb, Nf), SpinAxis) > 0)a = 2.0 * std::numbers::pi - a;

    double phi02 = glm::angle(glm::normalize(FoldingCurve[0].second->p - FoldingCurve[0].first->p),glm::normalize(FoldingCurve[1].first->p - FoldingCurve[0].first->p));
    double phim1 = glm::angle(glm::normalize(FoldingCurve.end()[-2].first->p - FoldingCurve.back().first->p),glm::normalize(FoldingCurve.back().second->p - FoldingCurve.back().first->p));

    glm::f64vec3 N4;
    for(auto&he: Edges){
        if(he->vertex->p3 == FoldingCurve[0].first->p3 && he->next->vertex->p3 == FoldingCurve[1].first->p3){N4 = he->face->getNormalVec();break;}
    }

    double a_ll = std::atan2(glm::dot(glm::cross(Nb,N4), glm::normalize(FoldingCurve[0].first->p3 - FoldingCurve[1].first->p3)),glm::dot(-Nb, N4));
    if(a_ll < 0)a_ll += 2.0*std::numbers::pi;
    double a_con = a_ll + std::numbers::pi;
    if(a_con > 2.0*std::numbers::pi)a_con -= 2.0*std::numbers::pi;
    if(a_con < 0)a_con +=2.0*std::numbers::pi;

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

    if(k < std::numbers::pi && IsMount){a_min = a_con + std::numbers::pi; a_max = 2.0 * std::numbers::pi;}
    if(k >= std::numbers::pi && IsMount){a_min = std::numbers::pi; a_max = std::numbers::pi + a_con;}
    if(k < std::numbers::pi && !IsMount){a_min = 0.0;a_max = a_con - std::numbers::pi;}
    if(k >= std::numbers::pi && !IsMount){a_min = a_con - std::numbers::pi; a_max = std::numbers::pi;}

    RevisionVertices::ObjData_v od_v = RevisionVertices::ObjData_v{a, FoldingCurve, Edges, Vertices, Poly_V, phi02, phim1, type, IsValidSigmoid};
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
         _FoldingAAAMethod(FoldingCurve, Edges, Poly_V, X.back(), phi02, phim1, type);
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
        }
        _FoldingAAAMethod(FC, data->Edges, data->Poly_V, a, data->phi02, data->phim1, data->type);

        f = fr(FC);
        for(int i = 0; i < (int)FC.size(); i++){
            FC[i].first->p = getVertex(FC[i].third->p2_ori, FC[i].first->p2_ori, X[i] + eps); FC[i].first->p3 = getVertex(FC[i].third->p3_ori, FC[i].first->p3_ori, X[i] + eps);
            _FoldingAAAMethod(FC, data->Edges, data->Poly_V, a, data->phi02, data->phim1, data->type);
            double fp = fr(FC);
            FC[i].first->p = getVertex(FC[i].third->p2_ori, FC[i].first->p2_ori, X[i] - eps); FC[i].first->p3 = getVertex(FC[i].third->p3_ori, FC[i].first->p3_ori, X[i] - eps);
            _FoldingAAAMethod(FC, data->Edges, data->Poly_V, a, data->phi02, data->phim1, data->type);
            double fm = fr(FC);
            grad[i] = (fp - fm)/(2.0*eps);
            FC[i].first->p = tmp[i][0]; FC[i].first->p3 = tmp[i][1];
        }
        if(grad.size() != FC.size()){
            _FoldingAAAMethod(FC, data->Edges, data->Poly_V, a + eps, data->phi02, data->phim1, data->type);
            double fp = fr(FC);
            _FoldingAAAMethod(FC, data->Edges, data->Poly_V, a - eps, data->phi02, data->phim1, data->type);
            double fm = fr(FC);
            grad.back() = (fp - fm)/(2.0*eps);
        }
    }
    std::cout <<"F_ruling(vertices) = " << f << std::endl;
    return f;
}



double RevisionVertices::Minimize_Vpos(const std::vector<double> &X, std::vector<double> &grad, void* f_data){
    ObjData_v *od = (ObjData_v *)f_data;
    double phi02 = od->phi02, phim1 = od->phim1;
    int type = od->type;
    FoldLine3d FC = od->FC;
    double fk, fvpos, fa;

    if(!grad.empty()){
        double a = (grad.size() != FC.size())? X.back(): od->a;
        std::vector<double> grad_vpos(grad.size()), grad_k(grad.size()), grad_a(grad.size(), 0);
        std::vector<HalfEdge*> Edges = od->Edges;
        std::vector<Vertex*> Poly_V = od->Poly_V;
        for(int i = 0; i < (int)FC.size(); i++ ){
            FC[i].first->p = getVertex(FC[i].third->p2_ori, FC[i].first->p2_ori, X[i]);
            FC[i].first->p3 = getVertex(FC[i].third->p3_ori, FC[i].first->p3_ori, X[i]);
        }
        _FoldingAAAMethod(FC, Edges, Poly_V, a, phi02, phim1, type);
        fk = F_k(FC, X, grad_k); fvpos = F_vpos(FC, X, grad_vpos);
        fa = abs(a - od->a); grad.back() = (abs(a + eps - od->a) - abs(a - eps - od->a)/(2.0*eps));
        for(int i = 0; i < (int)X.size(); i++){
            grad[i] = grad_k[i] + grad_vpos[i] + grad_a[i];
        }
    }
    std::cout <<"fk = " << fk <<", fvpos = " << fvpos <<", fa = " << fa <<  "  fk + fvpos = " << fk + fvpos + fa << std::endl;
    return fk + fvpos;
}

inline glm::f64vec3 calcCrossPoint_2Vector(Vertex* p1, Vertex* q1, Vertex* p2, Vertex* q2){
    Eigen::Matrix2d A; Eigen::Vector2d b;
    glm::f64vec3 v1 = q1->p - p1->p, v2 = q2->p - p2->p;
    b(0) = p2->p.x - p1->p.x; b(1) = p2->p.y - p1->p.y;
    A(0,0) = v1.x; A(0,1) = -v2.x;
    A(1,0) = v1.y; A(1,1) = -v2.y;
    Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
    return  x(0) * v1 + p1->p;
}
inline bool IsParallel(Vertex* p1, Vertex* q1, Vertex* p2, Vertex* q2){
    glm::f64vec3 v1 = glm::normalize(q1->p - p1->p), v2 = glm::normalize(q2->p - p2->p);
    if(abs(glm::dot(v1, v2)) >= 1.0 - 1e-5)return true;
    return false;
}

inline glm::f64vec3 calcTargetDistanceOnPlane(Vertex *p, Vertex *o, Vertex *v1, Vertex *v2){
    Eigen::Matrix2d A; Eigen::Vector2d b;
    b(0) = p->p.x - o->p.x; b(1) = p->p.y - o->p.y;
    A(0,0) = v1->p.x - o->p.x; A(0,1) = v2->p.x - o->p.x;
    A(1,0) = v1->p.y - o->p.y; A(1,1) = v2->p.y - o->p.y;
    Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
    return  x(0) * (v1->p3 - o->p3) + x(1) * (v2->p3 - o->p3) + o->p3;
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
            if(std::find_if(T_crs.begin(), T_crs.end(), [&e](CrvPt_FL* T){return T->p == e->vertex->p || T->p == e->next->vertex->p;}) != T_crs.end())continue;
            std::vector<double>arcT = BezierClipping(CtrlPts, e, dim);

            for(auto&t: arcT){
                if(t < 0 || 1 < t){std::cout<<"t is not correct value " << t << std::endl; continue;}
                if(std::find_if(T_crs.begin(), T_crs.end(), [&t](CrvPt_FL* T){return abs(T->s - t) < 1e-5;}) != T_crs.end())continue;
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
    std::sort(T_crs.begin(), T_crs.end(), [](CrvPt_FL* T, CrvPt_FL* T2){return T->s > T2->s;});

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
    for(auto&fc: FoldingCurve)fc.IsCalc = true;

    auto tmp = TrimPoints2(Edges, Faces, Vertices, FoldingCurve, tol);
    for(auto&V4d: FoldingCurve){
        if(std::find(tmp.begin(), tmp.end(), V4d) != tmp.end()){
        }else V4d.IsCalc = false;
    }
    std::vector<int>Vertices_Ind;
    for(int i = 0; i < FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}
    for(int i = 0; i < (int)Vertices_Ind.size() - 1; i++){
        for(int j = Vertices_Ind[i] + 1; j < Vertices_Ind[i+1]; j++){
            FoldingCurve[j].first->p =  calcCrossPoint_2Vector(FoldingCurve[j].first, FoldingCurve[j].third, FoldingCurve[Vertices_Ind[i]].first, FoldingCurve[Vertices_Ind[i+1]].first);
            FoldingCurve[j].first->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].first, FoldingCurve[Vertices_Ind[i]].first, FoldingCurve[Vertices_Ind[i]].third, FoldingCurve[Vertices_Ind[i+1]].third);
            FoldingCurve[j].third->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].third, FoldingCurve[Vertices_Ind[i]].first, FoldingCurve[Vertices_Ind[i]].third, FoldingCurve[Vertices_Ind[i+1]].third);
            FoldingCurve[j].second->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].second, FoldingCurve[Vertices_Ind[i]].first, FoldingCurve[Vertices_Ind[i]].second, FoldingCurve[Vertices_Ind[i+1]].second);
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
            FoldingCurve[j].first->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].first, FoldingCurve[Vertices_Ind[i]].first, FoldingCurve[Vertices_Ind[i]].third, FoldingCurve[Vertices_Ind[i+1]].third);
            FoldingCurve[j].third->p3 = calcTargetDistanceOnPlane(FoldingCurve[j].third, FoldingCurve[Vertices_Ind[i]].first, FoldingCurve[Vertices_Ind[i]].third, FoldingCurve[Vertices_Ind[i+1]].third);
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
    /*
   for(auto&fc: FoldingCurve){
        if(fc.IsCalc)continue;
        for(auto&e: Edges){
            if(e->edgetype == EdgeType::r && MathTool::is_point_on_line(fc.first->p, e->next->vertex->p, e->vertex->p)){
                e->r->Gradation = 0;
                break;
            }
        }
    }*/
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
}

void _FoldingAAAMethod(std::vector<Vertex4d>& FoldingCurve, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Poly_V,
                                 double a, double phi02, double phim1, int type){

    auto getClosestVertex = [&](Vertex *v, Vertex* o,  std::vector<Vertex4d> FoldingCurve){
       Vertex *V_max = nullptr;
       double t_max = -1;
       for(auto&fc: FoldingCurve){
           if(!fc.IsCalc)continue;
           if((glm::length(fc.second->p - v->p) < 1e-7|| glm::length(fc.second->p - o->p) < 1e-7))continue;
           if(MathTool::is_point_on_line(fc.second->p, v->p, o->p)){
               double t = glm::distance(fc.second->p, o->p)/glm::distance(v->p, o->p);
               if(t > t_max){
                   t_max = t; V_max = fc.second;
               }
           }
       }
       return V_max;
   };

    auto SetOnPlane = [&](Vertex4d& V, std::vector<Vertex*>& Poly_V, Vertex4d& p, Vertex4d& p2, Vertex *q, Vertex *q2, int IsEnd){
        //p, q: 左、p2, q2：右
        glm::f64vec3 vec;
        if(IsEnd == -1){vec = glm::normalize(q2->p - p2.first->p);}//左が端
        else if(IsEnd == 1){vec = glm::normalize(q->p - p.first->p);}//右が端
        else{
            if(!IsParallel(p.first, q, p2.first, q2))vec = glm::normalize(calcCrossPoint_2Vector(p.first, q, p2.first, q2) - V.first->p);
            else vec = glm::normalize(q->p - p.first->p);
        }
        if(glm::dot(glm::normalize(V.third->p - V.first->p), vec) > 0) vec = -1000. * vec + V.first->p;
        else vec = 1000. * vec + V.first->p;
        Vertex *Qv = new Vertex(vec);

        for(int k = 0; k < (int)Poly_V.size(); k++){
            glm::f64vec3 q;
            if(IsParallel(V.first, Qv, Poly_V[k], Poly_V[(k + 1) % (int)Poly_V.size()]))continue;
            q = calcCrossPoint_2Vector(V.first, Qv, Poly_V[k], Poly_V[(k + 1) % (int)Poly_V.size()]);
            if(MathTool::is_point_on_line(q, Poly_V[k]->p, Poly_V[(k + 1) % (int)Poly_V.size()]->p) &&
                    MathTool::is_point_on_line(q, V.first->p, Qv->p) ){
                V.second->p = q;
                V.second->p3 = calcTargetDistanceOnPlane(V.second, q2, p.first, p2.first);
            }
        }
    };

    double a2 = 0, l;
    glm::f64vec3 e, e2, x, e_bef;
    glm::f64vec3 N, r, n, befN;
    std::vector<int> Vertices_Ind;
    for(int i = 0; i < (int)FoldingCurve.size(); i++){if(FoldingCurve[i].IsCalc)Vertices_Ind.push_back(i);}
    double tau, dir, _tau;
    a = (a < 0.0)? a + 2.0 * std::numbers::pi: (a > 2.0*std::numbers::pi)? a - 2.0*std::numbers::pi: a;
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
            dir = glm::dot(-e, glm::cross(befN, srfN));
            if(type == 0)tau = std::acos(glm::dot(e_bef, e_next));
            if(type == 1)tau = std::acos(glm::dot(srfN, befN));
            _tau = std::acos(glm::dot(srfN, befN));
            if(dir <= 0)_tau = 2.0*std::numbers::pi - _tau;
            tau = _tau;
            a = (a2 - tau <= 0) ? a2 - tau + 2.0 * std::numbers::pi: a2 - tau;
        }
        Vertex *AxisV = findAxisVertex(fc, fc_bef.first, fc_next.first);
        if(AxisV == nullptr)continue;
        glm::f64vec3 Axis = (AxisV->p3 - fc.first->p3)/glm::length((AxisV->p3 - fc.first->p3));
        glm::f64vec3 Axis_2d = (AxisV->p - fc.first->p)/glm::length((AxisV->p - fc.first->p));
        double phi3 = std::acos(glm::dot((fc_next.first->p3 - fc.first->p3)/glm::length((fc_next.first->p3 - fc.first->p3)),Axis));
        double phi4 = std::acos(glm::dot((fc_bef.first->p3 - fc.first->p3)/glm::length((fc_bef.first->p3 - fc.first->p3)), Axis));

        double k = 2.0 * std::numbers::pi - phi3 - phi4;
        double _phi1 = std::atan((std::cos(k) - std::cos(beta))/(std::sin(beta)*std::cos(a)- std::sin(k)));
        double tmp2 = (std::sin(beta)*std::cos(a)- std::sin(k)), tmp3 = (std::cos(k) - std::cos(beta));
        double df_k = glm::degrees(std::numbers::pi - k);
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
        a2 = std::atan2(sin_a, cos_a);
        if(ind == (int)FoldingCurve.size() -2)break;
        e_bef = MathTool::ProjectionVector(e, e2, true);
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
                auto p = calcTargetDistanceOnPlane(FoldingCurve[Vertices_Ind[0]].second, v_clst,  FoldingCurve[Vertices_Ind[j]].first, FoldingCurve[Vertices_Ind[j+1]].first);
                FoldingCurve[Vertices_Ind[0]].second->p3 = p;
            }
            for(int i = Vertices_Ind[0] + 1; i < Vertices_Ind[1]; i++)
                SetOnPlane(FoldingCurve[i], Poly_V, FoldingCurve[Vertices_Ind[1]], FoldingCurve[Vertices_Ind[0]], FoldingCurve[Vertices_Ind[0]].second, v_clst, 1);
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
            HalfEdge *h_clst = nullptr;
            for(auto&h: v_clst->halfedge){if(h->edgetype != EdgeType::ol)h_clst = h;}
            if(h_clst != nullptr){
                auto p = calcTargetDistanceOnPlane(FoldingCurve[Vertices_Ind.back()].second, v_clst,  h_clst->next->vertex, h_clst->next->next->vertex);
                FoldingCurve[Vertices_Ind.back()].second->p3 = p;
            }else {
                std::cout <<"can't find correct edge left" << std::endl;
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
