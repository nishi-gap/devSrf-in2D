#include "mathtool.h"

namespace MathTool{
glm::f64vec3 bspline(std::vector<glm::f64vec3>&CtrlPts, double t, int dim, std::vector<double>Knot){
    glm::f64vec3 vec{0,0,0};
    for(int j = 0; j < (int)CtrlPts.size(); j++){
        double b = basis((int)CtrlPts.size(), j,dim, t + FLT_EPSILON,Knot);
        vec += CtrlPts[j] *  b;
    }
    return vec;
}

double distP2L(glm::f64vec3 la, glm::f64vec3 lb, glm::f64vec3& p, glm::f64vec3& q){
    double s = glm::dot((p - la), (lb - la))/glm::dot(lb - la, lb - la);
    if(0 <= s && s <= 1){
        q = la + s * (lb - la);
        return glm::distance(p,q);
    }
    return -1;
}



glm::f64vec3 getIntersectionPoint(glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3, glm::f64vec3& p4){
    double det = (p1.x - p2.x) * (p4.y - p3.y) - (p4.x - p3.x) * (p1.y - p2.y);
    double t = ((p4.y - p3.y) * (p4.x - p2.x) + (p3.x - p4.x) * (p4.y - p2.y)) / det;

    return glm::f64vec3{ t * p1.x + (1.0 - t) * p2.x, t * p1.y + (1.0 - t) * p2.y, 0};
}

bool IsIntersect(glm::f64vec3&p1, glm::f64vec3&p2, glm::f64vec3&p3, glm::f64vec3&p4, bool ConsiderEnd){
    auto set3pt = [](glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3) {
        return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
    };
    double t1 = set3pt(p1, p2, p3);
    double t2 = set3pt(p1, p2, p4);
    double t3 = set3pt(p3, p4, p1);
    double t4 = set3pt(p3, p4, p2);
    if(ConsiderEnd){
        if (t1 * t2 <= 0 && t3 * t4 <= 0) return true;//交点を持つ
    }
    else{
        if (t1 * t2 < 0 && t3 * t4 < 0) return true;//交点を持つ
    }
    return false;
}

void Triangulation(std::vector<glm::f64vec3>&input, std::vector<std::array<glm::f64vec3, 3>>&output){
    output.clear();
    std::vector<glm::f64vec3> Edges;
    std::copy(input.begin(), input.end(), back_inserter(Edges) );
    int n = Edges.size();
    while(Edges.size() >= 3){
        n = Edges.size();
        for(int i = 0; i < n; i++){
            int prev = (n + i - 1) % n;
            int next = (n + i + 1) % n;
            std::array<glm::f64vec3, 3> tri = {Edges[prev], Edges[i], Edges[next]};
            bool elimTriMesh = true;
            for(int j = 0; j < n - 3; j++){
                glm::f64vec3 p = Edges[(next + 1 + j) % n];
                bool check1 = hasPointInTriangle3D(p, tri);
                bool check2 = IsAngleLessThan180(tri[1], tri[0], tri[2]);
                if(check1 || !check2){
                    elimTriMesh = false;
                    break;
                }
            }
            if(elimTriMesh){
                Edges.erase(Edges.begin() + i);
                output.push_back(tri);
                break;
            }
        }

    }
}

bool IsAngleLessThan180(glm::f64vec3& o, glm::f64vec3& a, glm::f64vec3& b){
    glm::f64vec3 ao = glm::normalize(a - o);
    glm::f64vec3 bo = glm::normalize(b - o);

    glm::f64vec3 Normal = glm::cross(ao, bo);
    glm::f64vec3 BiNormal = glm::cross(ao, Normal);
    if(glm::dot(BiNormal, bo) < 0){return true;}
    return false;
}

bool hasPointInTriangle3D(glm::f64vec3 p, std::array<glm::f64vec3, 3>& V){
    double angle = 0.0;
    glm::f64vec3 v1, v2;
    int n = V.size();
    for(int i = 0; i < n; i++){
        v1 = glm::normalize(V[i] - p);
        v2 = glm::normalize(V[(i + 1) % n] - p);
        angle += std::acos(glm::dot(v1, v2));
    }
    if(angle >= 2 * std::numbers::pi - DBL_EPSILON)return true;
    return false;
}


std::vector<double> _bezierclipping(std::vector<glm::f64vec3>&CtrlPts_base, std::vector<glm::f64vec3>&CtrlPts_cur, std::array<glm::f64vec3, 2>& line, int dim){
    double t_min = 1, t_max = 0;
    glm::f64vec3 p = line[0], q = line[1];
    std::vector<glm::f64vec3> D = GrahamScan(CtrlPts_cur);
    int n = D.size();
    std::vector<double> T;
    for(int i = 0; i < n - 1; i++){
        glm::f64vec3 v = D[i], v2 = D[i + 1];
        double xp = ((p.y * q.x - p.x * q.y)*(v2.x - v.x) - (v.y * v2.x - v.x * v2.y)*(q.x - p.x)) / ((v2.y - v.y)*(q.x - p.x) - (v2.x - v.x)*(q.y - p.y));
        double yp = ((p.y * q.x - p.x * q.y)*(v2.y - v.y) - (v.y * v2.x - v.x * v2.y)*(q.y - p.y)) / ((v2.y - v.y)*(q.x - p.x) - (v2.x - v.x)*(q.y - p.y));
        glm::f64vec3 pt_new{xp, yp, 0};
        if(is_point_on_line(pt_new, p, q) && is_point_on_line(pt_new, v, v2)){ T.push_back(xp); }
    }
    if(T.size() == 0)return {};

    if(T.size() == 1){return{(p.x + q.x)/2}; }
    t_min = *std::min_element(T.begin(), T.end()); t_min = (t_min < 0) ? 0: (t_min > 1)? 1: t_min;
    t_max = *std::max_element(T.begin(), T.end()); t_max = (t_max < 0) ? 0: (t_max > 1) ? 1:  t_max;
    std::array<glm::f64vec3, 2> next_line = std::array{glm::f64vec3{t_min, 0,0}, glm::f64vec3{t_max, 0,0}};

    if(abs(t_max -  t_min) < DBL_EPSILON){
        return {t_max};
        //return {(t_max + t_min)/2.0};
    }
    std::pair<std::vector<glm::f64vec3>, std::vector<glm::f64vec3>> _bez = BezierSplit(CtrlPts_base, t_max);
    double bez_t = t_min / (t_max);
    _bez = BezierSplit(_bez.first, bez_t);
    if(abs(glm::distance(p,q) - abs(t_max - t_min)) < DBL_EPSILON){
        std::pair<std::vector<glm::f64vec3>, std::vector<glm::f64vec3>> bez_spl = BezierSplit(_bez.second, 0.5);
        std::vector<glm::f64vec3> b1 = bez_spl.first, b2 = bez_spl.second;
        std::array<glm::f64vec3, 2> next_line2; std::copy(next_line.begin(), next_line.end(), next_line2.begin());
        std::vector<double>t1 = _bezierclipping(CtrlPts_base, b1, next_line, dim), t2 = _bezierclipping(CtrlPts_base, b2, next_line2, dim);
        for(auto&t: t2)t1.push_back(t);
        return t1;
    }

    return _bezierclipping(CtrlPts_base, _bez.second, next_line, dim);
}

//http://www-ikn.ist.hokudai.ac.jp/~k-sekine/slides/convexhull.pdf
//https://kajindowsxp.com/graham-algo/
std::vector<glm::f64vec3> GrahamScan(std::vector<glm::f64vec3>& Q){
    std::vector<glm::f64vec3> S;
    if(Q.size() < 3)return Q;
    glm::f64vec3 p_ml = Q[0];
    for(auto&p: Q){
        if(p_ml.y > p.y)p_ml = p;
        else if(p_ml.y == p.y && p_ml.x > p.x) p_ml = p;
    }
    std::vector<std::pair<double, glm::f64vec3>>Args;
    for(auto&p: Q){
        if(p == p_ml)continue;
        double phi = atan2(p.y - p_ml.y, p.x - p_ml.x);
        bool hasSameAngle = false;
        for(auto& X: Args){
            if(X.first == phi){
                if(glm::distance(X.second, p_ml) < glm::distance(p, p_ml))X.second = p;
                hasSameAngle = true;
            }
        }
        if(!hasSameAngle)Args.push_back(std::make_pair(phi, p));
    }
    // compare only the first value
    std::sort(Args.begin(), Args.end(),[](auto const& x, auto const& y) {return x.first < y.first; });
    if(Args.size() >= 1)S.push_back(Args[Args.size() - 1].second);
    S.push_back(p_ml);
    glm::f64vec3 top, next;

    for(int i = 0; i < (int)Args.size(); i++){
        do{
            top = S.back();
            S.pop_back();
            next = S.back();
            if(SignedArea(next, Args[i].second,top) <= 0){
                S.push_back(top);
                S.push_back(Args[i].second);
                break;
            }
        }while(1);

    }
    return S;
}
//xy平面上に乗っていると仮定
double SignedArea(glm::f64vec3 a, glm::f64vec3 b, glm::f64vec3 p){
    glm::f64vec3 v = a-p, v2 = b-p;
    return v.x * v2.y - v.y*v2.x;
}

bool is_point_on_line(glm::f64vec3 p, glm::f64vec3 lp1, glm::f64vec3 lp2){
    double ac = glm::distance(p, lp1), bc = glm::distance(p, lp2), lp = glm::distance(lp1, lp2);
    if(ac < 1e-9 || bc < 1e-9)return true;
    if(abs(lp - ac - bc) < 1e-9) return true;
    return false;
}

//split at t for bezier curve
std::pair<std::vector<glm::f64vec3>, std::vector<glm::f64vec3>> BezierSplit(std::vector<glm::f64vec3> CtrlPts, double t){
    std::vector<glm::f64vec3>lp, rp;
    std::vector<std::vector<glm::f64vec3>> Tree = de_casteljau_algorithm(CtrlPts, t);
    Tree.insert(Tree.begin(), CtrlPts);
    for(auto& d: Tree){
        lp.push_back(d[0]);
        rp.push_back(d[d.size() - 1]);
    }
    return {lp, rp};

}

std::vector<std::vector<glm::f64vec3>> de_casteljau_algorithm(std::vector<glm::f64vec3> CtrlPts, double t){
    std::vector<glm::f64vec3> Q;
    std::vector<std::vector<glm::f64vec3>> _Q;
    glm::f64vec3 prev = CtrlPts[0];
    for(int i = 1; i < (int)CtrlPts.size(); i++){
        glm::f64vec3 new_p = (1. - t) * prev + t * CtrlPts[i];
        Q.push_back(new_p);
        prev = CtrlPts[i];
    }
    if(Q.size() == 0)return _Q;

    _Q.push_back(Q);
    if(Q.size() == 1)return _Q;

    std::vector<std::vector<glm::f64vec3>> tmp = de_casteljau_algorithm(Q,t);
    _Q.insert(_Q.end(), tmp.begin(), tmp.end());
    return _Q;
}


double basis(int n, int i, int p, double u, std::vector<double>& U){
    //if((U[0] <= u && u <= U[i]) && (U[i+p] <= u && u <= U[n+p]))return 0;
    if(p == 0){return (U[i] <= u && u < U[i+1]) ? 1: 0;}
    //double a = (U[i+p] != U[i])? (u - U[i])/(U[i+p] - U[i]): (u == U[i+p] && U[i] != 0) ? 1: 0;
    //double b = (U[i+p+1] != U[i+1])? (U[i+p+1] - u)/(U[i+p+1] - U[i+1])? (U[i+1] == u && U[i+p+1] != 0): 1: 0;
    double a = (U[i+p] != U[i])? (u - U[i])/(U[i+p] - U[i]): 0;
    double b = (U[i+p+1] != U[i+1])? (U[i+p+1] - u)/(U[i+p+1] - U[i+1]): 0;
    return a * basis(n,i,p-1,u,U) + b * basis(n,i+1,p-1,u,U);
}

double factorial(int n){
    return (n > 1)? factorial(n - 1) * n : 1;
}

double cmb(int n , int i){
    return factorial(n) / (factorial(i) * factorial(n - i));
}

double BernsteinBasisFunc(int n, int i, double t){
    return cmb(n, i) * std::pow(t, i) * std::pow(1.0 - t, n - i);
}

glm::f64vec3 bezier(std::vector<glm::f64vec3>& CtrlPts, double t, int dim){
    glm::f64vec3 v = {0.f,0.f, 0.f};
    for (int i = 0; i < dim+1; i++) v +=  BernsteinBasisFunc(dim, i, t)* CtrlPts[i];
    return v;
}

void diffBezier(std::vector<glm::f64vec3>& dP, std::vector<glm::f64vec3>& CtrlPts, double t, int dim){
    dP.resize(3);
    double eps = 1e-5;
    glm::f64vec3 v = bezier(CtrlPts, t, dim);
    glm::f64vec3 vm = bezier(CtrlPts, t - eps, dim);
    glm::f64vec3 vp = bezier(CtrlPts, t + eps, dim);
    glm::f64vec3 vm2 = bezier(CtrlPts, t - 2.0 * eps, dim);
    glm::f64vec3 vp2 = bezier(CtrlPts, t + 2.0 * eps, dim);
    dP[0] = (vp - vm)/(2.0 * eps);
    dP[1] = (vp - 2.0 * v + vm)/pow(eps, 2);
    dP[2] = (-vm2 + 2.0 * vm - 2.0 * vp + vp2)/(2.0 * pow(eps, 3));
}

std::vector<double> LSM_apply(std::vector<double>&_y, int dim){
    //type: 0->x, 1->y, 2->x*x, 3->xy
    auto sum = [](std::vector<double>&y, int type){
        double z = 0.0;
        for(int i = 0; i < (int)y.size(); i++)z += (type == 0) ? i: (type==1)?y[i]: (type==2)?i*i: i*y[i];
        return z;
    };
    std::string file = "result_lsm.csv";
    std::ofstream ofs(file);
    int n = _y.size();
    std::vector<double> Y(n);

    ofs << "row data(y), result(y)" << std::endl;
    if(dim == 1){
        double xy = sum(_y,3), xx = sum(_y,2), x = sum(_y,0), y = sum(_y,1);
        double a = ((double)n*xy - x*y)/((double)n * xx - x*x);
        double b = (xx*y - xy*x)/((double)n*xx - x*x);
        for(int i = 0; i < n; i++){
            Y[i] = a*(double)i + b;
            ofs << _y[i] << ", " << Y[i] << std::endl;
        }
    }

    return Y;
}

void PCA(std::vector<glm::f64vec3>& _X, std::vector<glm::f64vec3>& BasisVectors){
    using namespace Eigen;

    BasisVectors.resize(2);

    int n = _X.size();
    MatrixXd X(n, 3);
    for(int i = 0; i < n; i++){X(i,0) = _X[i].x; X(i,1) = _X[i].y; X(i,2) = _X[i].z;}

    MatrixXd Xstd = X.rowwise() - X.colwise().mean();
    for (int c = 0; c < 3; c++) {
      auto col = X.col(c).array();
      double std_dev = sqrt((col - col.mean()).square().sum() / col.size() - 1);
      std::cout <<"std_dev " << std_dev <<std::endl;
      if(!std::isfinite(std_dev)){
          if(c == 0){BasisVectors[0] = glm::f64vec3{0,1,0};BasisVectors[1] = glm::f64vec3{0,0,1};}
          else if(c == 1){BasisVectors[0] = glm::f64vec3{1,0,0};BasisVectors[1] = glm::f64vec3{0,0,1};}
          else{BasisVectors[0] = glm::f64vec3{1,0,0};BasisVectors[1] = glm::f64vec3{0,1,0};}
          return;
      }
      Xstd.col(c) /= std_dev;
    }

    MatrixXd S = (Xstd.adjoint() * Xstd) / double(Xstd.rows() - 1);
    std::cout<<"S" <<std::endl;
    std::cout<<S <<std::endl;

     SelfAdjointEigenSolver<MatrixXd> ES(S);
     VectorXd e = ES.eigenvalues();
     MatrixXd U = ES.eigenvectors();
     std::cout<<"eigen value : " << e << std::endl;
     std::cout<<"eigen vectors : " << U << std::endl;
     int minIndex, ind = 0;
     e.minCoeff(&minIndex);
     for(int i = 0; i < 3; i++){
         if(i == minIndex)continue;
         std::cout <<"basis vector " << i << " : " << U.col(i) << std::endl;
         BasisVectors[ind++] = glm::f64vec3{U(0,i), U(1,i), U(2,i)};
     }
}

std::vector<glm::f64vec3> getPlaneFromCurvePoints(std::vector<glm::f64vec3>& Points, std::vector<glm::f64vec3>& BasisVectors){
    std::vector<glm::f64vec3> PointsOnPlane;
    PCA(Points, BasisVectors);
    const glm::f64vec3 N = glm::normalize(glm::cross(BasisVectors[0], BasisVectors[1]));
    for(const auto&p: Points)PointsOnPlane.push_back(p - glm::dot(p, N)*N);

    return PointsOnPlane;
}

glm::f64vec3 ProjectionVector(glm::f64vec3 v, glm::f64vec3 n, bool Isnormalize){
    n = glm::normalize(n);
    glm::f64vec3 V = v - glm::dot(v,n)* n;
    return (Isnormalize)? V/glm::length(V): V;
}

double AngleIn2Edges(HalfEdge *p, HalfEdge *p2, bool Is3d){
  if(Is3d)return acos(glm::dot(glm::normalize(p->next->vertex->p3 - p->vertex->p3), glm::normalize(p2->next->vertex->p3 - p2->vertex->p3)));
  return acos(glm::dot(glm::normalize(p->next->vertex->p - p->vertex->p), glm::normalize(p2->next->vertex->p - p2->vertex->p)));
}

}

/*
BaseCRV::BaseCRV(PaintTool _type){
    isempty = true;
    type = _type;
    CtrlPts.clear();
    CurvePts.clear();
}

bool BaseCRV::addCtrlPt(glm::f64vec3 p, int dim){
    CtrlPts.push_back(p);
    bool res = setCurve(dim);
    return res;
}

bool BaseCRV::deleteCtrlPt(glm::f64vec3 p, int dim){
    double dst = 5.0;
    std::vector<glm::f64vec3>::iterator delV = CtrlPts.end();
    for(auto it = CtrlPts.begin(); it != CtrlPts.end(); ++it){
        if(glm::distance(p, *it) < dst){
            dst = glm::distance(p, *it); delV = it;
        }
    }
    if(delV == CtrlPts.end())return false;
    CtrlPts.erase(delV);
    updateCurve(dim);
    return true;
}

int BaseCRV::movePtIndex(glm::f64vec3 p){
    double dst = 5.0;
    int ind = -1;
    for(int i = 0; i < (int)CtrlPts.size(); i++){
        if(glm::distance(p, CtrlPts[i]) < dst){
            dst = glm::distance(p, CtrlPts[i]); ind = -1;
        }
    }
    return ind;
}

bool BaseCRV::moveCtrlPt(glm::f64vec3 p, int moveIndex, int dim){
    if(moveIndex == -1 || (int)CtrlPts.size() < moveIndex)return false;
    CtrlPts[moveIndex] = p;
    bool res = updateCurve(dim);
    return res;
}

bool BaseCRV::setCurve(int dim){
    if(((type == PaintTool::Bspline_r || type == PaintTool::FoldLine_bezier) && (int)CtrlPts.size() <= dim) ||
            ((type == PaintTool::Arc_r || type == PaintTool::FoldLine_arc) && (int)CtrlPts.size() < 3) ||
            ((type == PaintTool::Line_r || type == PaintTool::FoldLine_line) && (int)CtrlPts.size() < 2))return false;
    int CurveSize = 1000;
    CurvePts.resize(CurveSize);

    if(type == PaintTool::FoldLine_bezier){
        for(int i = 0; i < CurveSize; i++)CurvePts[i] +=  MathTool::bezier(CtrlPts, (double)i/(double)CurveSize, dim);
    }else if(type == PaintTool::Bspline_r){
        int knotSize = (int)CtrlPts.size() + dim + 1;
        std::vector<double>Knot(knotSize);
        for(int j = 0; j < knotSize; j++)Knot[j] = (double)j/(double)knotSize;
        for(int j = 0; j < dim + 1;j++){ Knot[j] = 0; Knot[knotSize - 1 - j] = 1;}
        for(int i = 0; i < CurveSize; i++)CurvePts[i] = MathTool::bspline(CtrlPts, (double)i/(double)CurveSize, dim, Knot);
    }else if(type == PaintTool::Arc_r || type == PaintTool::FoldLine_arc)Arc();
    else if(type == PaintTool::Line_r || type == PaintTool::FoldLine_line)Line();
    return true;
}

bool BaseCRV::updateCurve(int dim){
    int CurveSize = CurvePts.size();
    if(type == PaintTool::FoldLine_bezier){
        for(int i = 0; i < CurveSize; i++)CurvePts[i] +=  MathTool::bezier(CtrlPts, (double)i/(double)CurveSize, dim);
    }else if(type == PaintTool::Bspline_r){
        std::vector<double>Knot;
        for(int i = 0; i < CurveSize; i++)CurvePts[i] = MathTool::bspline(CtrlPts, (double)i/(double)CurveSize, dim, Knot);
    }else if(type == PaintTool::Arc_r || type == PaintTool::FoldLine_arc)Arc();
    else if(type == PaintTool::Line_r || type == PaintTool::FoldLine_line)Line();
    return true;
}

void BaseCRV::Arc(){
    if((int)CtrlPts.size() < 3)return;
    int CurveSize = CurvePts.size();
    for(int i = 0; i < CurveSize; i++){
        double l = glm::distance(CtrlPts[0], CtrlPts[1]);
        double l2 = glm::distance(CtrlPts[0], CtrlPts[2]);
        CtrlPts[2] = l/l2 * (CtrlPts[2] - CtrlPts[0]) + CtrlPts[0];
        double phi = acos(glm::dot(glm::normalize(CtrlPts[1] - CtrlPts[0]), glm::normalize(CtrlPts[2] - CtrlPts[0])));
        glm::f64vec3 axis = glm::cross(glm::normalize(CtrlPts[1] - CtrlPts[0]), glm::normalize(CtrlPts[2] - CtrlPts[0]));
        glm::f64mat4x4  T = glm::translate(CtrlPts[0]);
        glm::f64mat4x4 invT = glm::translate(-CtrlPts[0]);
        glm::f64mat4x4  R;
        R = glm::rotate(phi * (double)i/(double)CurveSize, axis);
        CurvePts[i] = T * R * invT * glm::f64vec4{CtrlPts[1],1};
    }
    return;
}

void BaseCRV::Line(){
    int CurveSize = CurvePts.size();
    if(CtrlPts.size() < 2)return;
    while(CtrlPts.size() != 2){
        CtrlPts.erase(CtrlPts.end() - 1);
        CtrlPts.shrink_to_fit();
    }
    for(int i = 0; i < CurveSize; i++)CurvePts[i] = (double)i/(double)(CurveSize - 1) * (CtrlPts[1] - CtrlPts[0]) + CtrlPts[0];
    return;
}

void SmoothCRV::InsertCtrlPt(glm::f64vec3& p){
    int ind;

    int InsertType;//-1：どこにものっかっていない　0：曲線上　1：制御点を結んだ線上
    if(CurvePts.empty())InsertType =  -1;
    else{
        double lc = 5, lp = 5;
        for(int i = 0; i < (int)CurvePts.size(); i++){
            if(glm::distance(p,CurvePts[i]) < lc){
                lc = glm::distance(p,CurvePts[i]); ind = i;
            }
        }
        glm::f64vec3 q;
        for(int i = 0; i < (int)CtrlPts.size() - 1; i++) {
            double l = MathTool::distP2L(CtrlPts[i], CtrlPts[i + 1], p, q);
            if(l != -1 && l < lp){
                lp = l; ind = i;
            }
        }
        InsertType = (lc == 5 && lp == 5)? -1: (lc < lp)? 0: 1;
    }

    double dist;
    std::vector<int>CurveIndexs((int)CtrlPts.size(), -1);
    for(int j = 0; j < (int)CtrlPts.size(); j++){
        dist = 100;
        for(int i = 0; i < (int)CurvePts.size(); i++){
            if(dist > glm::distance(CurvePts[i], CtrlPts[j])){
                dist = glm::distance(CurvePts[i], CtrlPts[j]); CurveIndexs[j] = i;
            }
        }
    }
    if(InsertType == -1){
        double minDist = 200.0;
        double l;
        glm::f64vec3 q;
        for(int i = 0; i < (int)CtrlPts.size() - 1; i++){
            l = MathTool::distP2L(CtrlPts[i], CtrlPts[i + 1], p, q);
            if(l != -1 && l < minDist){ minDist = l; InsertPoint = p; InsertPointSegment = i;}
        }
        if(minDist == 100.0){InsertPointSegment = -1;}
    }
    else if(InsertType == 0){//曲線上
        for(int i = 0; i < (int)CurveIndexs.size() - 1; i++){
            if(ind >=CurveIndexs[i] && ind < CurveIndexs[i + 1]){
                double s = (double)(ind - CurveIndexs[i])/(double)(CurveIndexs[i + 1] - CurveIndexs[i]);
                InsertPoint = s * (CtrlPts[i + 1] - CtrlPts[i]) + CtrlPts[i];
                InsertPointSegment = i;
                return;
            }
        }
    }else if(InsertType == 1){//直線上
        glm::f64vec3 q;
        double l = MathTool::distP2L(CtrlPts[ind], CtrlPts[ind + 1], p, q);
        if(l != -1){
            InsertPoint = q;
            InsertPointSegment = ind;
        }else{
            InsertPointSegment = -1;
        }
    }
    return;
}
void SmoothCRV::SetNewPoint(){
    if(InsertPointSegment == -1)return;
    CtrlPts.insert(CtrlPts.begin() + InsertPointSegment + 1, InsertPoint);
}

bool SmoothCRV::setPoint(std::vector<HalfEdge*>&outline, glm::f64vec3 N, glm::f64vec3& cp, std::vector<glm::f64vec3>& P){
    glm::f64vec3 N0 = N + cp, N1 = -N + cp;
    std::vector<glm::f64vec3> crossPoint;
    P.clear();
    bool IsIntersected = false;

    for(auto&edge: outline){
        glm::f64vec3 pt;
        bool hasCrossPoint = edge->hasCrossPoint(N0, N1, pt, false);
        if(hasCrossPoint){
            crossPoint.push_back(pt);
            IsIntersected = true;
        }
    }

    std::vector<double>dist;
    for(auto& p: crossPoint)dist.push_back(glm::distance(p, N1));
    for(int i = 0; i < (int)crossPoint.size(); i++){
        for(int j = i + 1; j < (int)crossPoint.size(); j++){
            if(dist[j] < dist[i]){
                std::swap(crossPoint[i], crossPoint[j]); std::swap(dist[i], dist[j]);
            }
        }
    }
    glm::f64vec3 minPt = N0;
    for(auto&p: crossPoint){
        minPt = (glm::distance(p, cp) < glm::distance(minPt, cp)) ? p : minPt;
    }
    std::vector<glm::f64vec3>::iterator itr = std::find(crossPoint.begin(), crossPoint.end(), minPt);
    if(itr == crossPoint.end())return IsIntersected;
    int minInd = std::distance(crossPoint.begin(), itr);
    minInd /= 2;
    for(int i = 0; i < (int)crossPoint.size(); i++){
        P.push_back(crossPoint[i]);
    }
    return IsIntersected;
}

bool SmoothCRV::setRulingVector(std::vector<HalfEdge*>& SurfaceEdge, int dim){
    if((type == PaintTool::Bspline_r && (int)CtrlPts.size() <= dim) || (type == PaintTool::Arc_r && (int)CtrlPts.size() < 3) || (type == PaintTool::Line_r && (int)CtrlPts.size() < 2)) return false;
    rulings.resize(DivSize - 1);
    Vertices.resize(2*(DivSize-1));
    double l = 1000;//適当に大きな値
    std::vector<glm::f64vec3> crossPoint;
    std::vector<std::vector<glm::f64vec3>> CrossPoints;
    int sind = -1, eind = (int)CurvePts.size() - 1;
    glm::f64vec3 V, N;

    if(type == PaintTool::Bspline_r){
            double t = 0.0;
            int knotSize = (int)CtrlPts.size() + dim + 1;
            std::vector<double>Knot(knotSize);
            for(int j = 0; j < knotSize; j++)Knot[j] = (double)j/(double)knotSize;
            for(int j = 0;j< dim + 1;j++){ Knot[j] = 0; Knot[knotSize - 1 - j] = 1;}
            int i = 0;
            int sind = -1, eind = (int)CurvePts.size();
            t = Knot[dim] + (Knot[(int)CtrlPts.size()] - Knot[dim]) * (1.0 / (double)CurveSize);
            while(i < (int)CurvePts.size()){
                glm::f64vec3 vec, vec2;
                vec = bspline(CtrlPts, t  + FLT_EPSILON,dim, Knot);
                vec2 = bspline(CtrlPts,t  - FLT_EPSILON,dim, Knot);
                V = (vec - vec2)/(2.0 * FLT_EPSILON);
                V = glm::normalize(V);
                N = l * glm::f64vec3{ -V.y, V.x, 0};
                if (glm::dot(glm::normalize(N), glm::f64vec3{0,1,0}) < 0) N *= -1;
                bool IsIntersected = setPoint(SurfaceEdge, N, CurvePts[i], crossPoint);
                CrossPoints.push_back(crossPoint);
                if(sind == -1 && IsIntersected)sind = i;
                if(!IsIntersected && sind != -1)eind = std::min(eind, i);
                i++;
                t += (Knot[(int)ControllPoints.size()] - Knot[curveDimention]) * (1.0/ (double)crvPtNum);
            }
    }else if(type == PaintTool::Arc_r){
        for(int i = 0; i < (int)CurvePts.size(); i++){
            if(i == 0)V = glm::normalize(CurvePts[1] - CurvePts[0]);
            else V = glm::normalize(CurvePts[i] - CurvePts[i - 1]);
            bool IsIntersected = setPoint(SurfaceEdge, N, CurvePts[i], crossPoint);
            CrossPoints.push_back(crossPoint);
            if(sind == -1 && IsIntersected)sind = i;
            if(!IsIntersected && sind != -1)eind = std::min(eind, i);
        }
    }else if(type == PaintTool::Line_r){
          V = glm::normalize(CtrlPts[1] - CtrlPts[0]);
          glm::f64vec3 N = l * glm::f64vec3{-V.y, V.x, 0};
          int sind = -1, eind = curveNum - 1;
          for(int i = 0; i < CurveSize; i++){
              bool IsIntersected = setPoint(SurfaceEdge, N, CurvePts[i], crossPoint);
              CrossPoints.push_back(crossPoint);
              if(sind == -1 && IsIntersected)sind = i;
              if(!IsIntersected && sind != -1)eind = std::min(eind, i);
          }
    }
    double n = double(sind);
    while(n < eind){
        Vertex *v1 = new Vertex(CrossPoints[floor(n)][0]), *v2 = new Vertex(CrossPoints[floor(n)][1]);
        HalfEdge *h1 = new HalfEdge(v1, EdgeType::r), *h2 = new HalfEdge(v2, EdgeType::r);
        h2->prev = h2->next = h2->pair =  h1; h1->prev = h1->next = h1->pair =  h2;
        Vertices.push_back(v1); Vertices.push_back(v2);
        rulings.push_back(h1); rulings.push_back(h2);
        n += (double)(eind - sind + 1)/(double)(DivSize - 1);
    }
    isempty = false;
}

bool SmoothCRV::updateRulingVector(const std::vector<HalfEdge*>& SurfaceEdge, int dim){

}
*/
