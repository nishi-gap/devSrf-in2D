#include "mathtool.h"

namespace MathTool{
glm::f64vec3 bspline(std::vector<glm::f64vec3>&CtrlPts, double t, int dim, std::vector<double>Knot){
    glm::f64vec3 vec{0,0,0};
    for(int j = 0; j < (int)CtrlPts.size(); j++){
        double b = basis(j,dim, t + FLT_EPSILON,Knot);
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

double set3pt(glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3) {
    return (p2.x - p1.x) * (p3.y - p1.y) - (p2.y - p1.y) * (p3.x - p1.x);
}

bool IsIntersect(glm::f64vec3&p1, glm::f64vec3&p2, glm::f64vec3&p3, glm::f64vec3&p4){
    double t1 = set3pt(p1, p2, p3);
    double t2 = set3pt(p1, p2, p4);
    double t3 = set3pt(p3, p4, p1);
    double t4 = set3pt(p3, p4, p2);

    if (t1 * t2 < 0 && t3 * t4 < 0) return true;//交点を持つ
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
    if(angle >= 2 * M_PI - FLT_EPSILON)return true;
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

    if(abs(t_max - t_min) < 1e-7){
        return {(t_max + t_min)/2};
    }
    std::pair<std::vector<glm::f64vec3>, std::vector<glm::f64vec3>> _bez = BezierSplit(CtrlPts_base, t_max, dim);
    double bez_t = t_min / (t_max);
    _bez = BezierSplit(_bez.first, bez_t, dim);
    if(abs(glm::distance(p,q) - abs(t_max - t_min)) < 1e-7){

        std::pair<std::vector<glm::f64vec3>, std::vector<glm::f64vec3>> bez_spl = BezierSplit(_bez.second, 0.5,  dim);
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
    if(Q.size() < 3)return S;
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

bool is_point_on_line(glm::f64vec3& p, glm::f64vec3& lp1, glm::f64vec3& lp2){
    double ac = glm::distance(p, lp1), bc = glm::distance(p, lp2), lp = glm::distance(lp1, lp2);
    if(abs(lp - ac - bc) < 1e-5) return true;
    return false;
}

//split at t for bezier curve
std::pair<std::vector<glm::f64vec3>, std::vector<glm::f64vec3>> BezierSplit(std::vector<glm::f64vec3> CtrlPts, double t, int dim){
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


double basis(int i, int p, double u, std::vector<double>& U){
    if(p == 0) return (U[i] <= u && u < U[i+1]) ? 1: 0;
    double a = (U[i+p] != U[i])? (u - U[i])/(U[i+p] - U[i]): 0;
    double b = (U[i+p+1] != U[i+1])? (U[i+p+1] - u)/(U[i+p+1] - U[i+1]): 0;
    return a * basis(i,p-1,u,U) + b * basis(i+1,p-1,u,U);
}

double factorial(int n){
    return (n > 1)? factorial(n - 1) * n : 1;
}

double cmb(int n , int i){
    return factorial(n) / (factorial(i) * factorial(n - i));
}

}


