#include "mathtool.h"

namespace DebugMode{
    Singleton& Singleton::getInstance(){
        static Singleton instance;
        return instance;
    }
    void Singleton::switchval(){ val = !val;}
    bool Singleton::isdebug() const{return val;}
}

namespace MathTool{
    double rad2deg(double a){return a * 180.0/std::numbers::pi;}
    double deg2rad(double a){return a * std::numbers::pi/180.0;}

    Eigen::Vector3d PCA(std::vector<Eigen::Vector3d>& X, Eigen::Vector3d& o) {


        o = Eigen::Vector3d(0, 0, 0);
        for (auto& x : X)o += x;
        o /= static_cast<double>(X.size());

        Eigen::MatrixXd M(3, static_cast<int>(X.size()));
        for (int i = 0; i < (int)X.size(); i++)M.col(i) = X[i] - o;
        Eigen::MatrixXd S = (M * (M.transpose()))/(double)X.size();//共分散行列
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ES(S);
        Eigen::Vector3d e,v0,v1,v2;
        //hoge
        //e = ES.eigenvalues(); Eigen::MatrixXd U = ES.eigenvectors();
        //Eigen::Vector3d v0 = Eigen::Vector3d(U.col(0).x(), U.col(0).y(), U.col(0).z()), v1 = Eigen::Vector3d(U.col(1).x(), U.col(1).y(), U.col(1).z()), v2 = Eigen::Vector3d(U.col(2).x(), U.col(2).y(), U.col(2).z());
        double minelem = std::min({e(0), e(1), e(2)});
        if(minelem == e(0)){
            return (v1.cross(v2)).normalized();
        }else if(minelem == e(1)){
            return (v1.cross(v2)).normalized();
        }else{
            return (v1.cross(v2)).normalized();
        }
    }
    Eigen::Vector3d CrossPointLineAndPlane(Eigen::Vector3d e, Eigen::Vector3d e2, Eigen::Vector3d o, Eigen::Vector3d p, Eigen::Vector3d v){
        e = e.normalized(); e2 = e2.normalized(); v = (v - p).normalized();
        //Eigen::Vector3d e = (ValidFC[i]->first->p3 - ValidFC[i-1]->first->p3).normalized(), e2 = (ValidFC[i-2]->first->p3 - ValidFC[i-1]->first->p3).normalized();
        //Eigen::Vector3d v = (ValidFC[i+1]->first->p3 - ValidFC[i+1]->third->p3).normalized();
        Eigen::Matrix3d A;
        A.col(0) = e; A.col(1) = e2; A.col(2) = -v;
        Eigen::Vector3d b = p - o;
        //Eigen::Vector3d b = ValidFC[i+1]->third->p3 - ValidFC[i-1]->first->p3;
        Eigen::Vector3d x = A.colPivHouseholderQr().solve(b);
        return x(2)*v + p;
    }

    Eigen::Vector3d bspline(std::vector<Eigen::Vector3d>&CtrlPts, double t, int dim, std::vector<double>Knot){
        Eigen::Vector3d vec(0,0,0);
        for(int j = 0; j < (int)CtrlPts.size(); j++){
            double b = basis((int)CtrlPts.size(), j,dim, t + 1e-9,Knot);
            vec += CtrlPts[j] *  b;
        }
        return vec;
    }

    bool IsIntersect(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, const Eigen::Vector3d& p4, Eigen::Vector3d& q, bool ConsiderEnd){
        auto set3pt = [](const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) {
            return (p2(0) - p1(0)) * (p3(1) - p1(1)) - (p2(1) - p1(1)) * (p3(0) - p1(0));
        };

        double det = (p1.x() - p2.x()) * (p4.y() - p3.y()) - (p4.x() - p3.x()) * (p1.y() - p2.y());
        double t = ((p4.y() - p3.y()) * (p4.x() - p2.x()) + (p3.x() - p4.x()) * (p4.y() - p2.y())) / det;
        q = t * p1 + (1.0 - t)*p2;

        double t1 = set3pt(p1, p2, p3), t2 = set3pt(p1, p2, p4), t3 = set3pt(p3, p4, p1), t4 = set3pt(p3, p4, p2);
        if(ConsiderEnd){
            if (t1 * t2 <= 0 && t3 * t4 <= 0) return true;//交点を持つ
        }
        else{
            if (t1 * t2 < 0 && t3 * t4 < 0) return true;//交点を持つ
        }

        return false;
    }

    bool IsAngleLessThan180(Eigen::Vector3d& o, Eigen::Vector3d& a, Eigen::Vector3d& b){
        Eigen::Vector3d ao = (a - o), bo = (b - o);
        ao = ao.normalized(); bo = bo.normalized();
        Eigen::Vector3d Normal = ao.cross(bo);
        Eigen::Vector3d BiNormal = ao.cross(Normal);
        return (BiNormal.dot(bo) < 0)? true: false;
    }

    bool hasPointInTriangle3D(const Eigen::Vector3d& p, std::array<Eigen::Vector3d, 3>& V){
        double angle = 0.0;
        Eigen::Vector3d v1, v2;
        int n = V.size();
        for(int i = 0; i < n; i++){
            v1 = V[i] - p; v1 = v1.normalized();
            v2 = V[(i + 1) % n] - p; v2 = v2.normalized();
            angle += std::acos(v1.dot(v2));
        }

        return (angle >= 2 * std::numbers::pi - 1e-9)? true: false;
    }


    std::vector<double> _bezierclipping(const std::vector<Eigen::Vector3d>&CtrlPts_base, std::vector<Eigen::Vector3d>&CtrlPts_cur, std::array<Eigen::Vector3d, 2>& line, int dim){
        double t_min = 1, t_max = 0;
        Eigen::Vector3d p = line[0], q = line[1];
        auto D = GrahamScan(CtrlPts_cur);
        int n = D.size();
        std::vector<double> T;
        for(int i = 0; i < n - 1; i++){
            double xp = ((p.y() * q.x() - p.x() * q.y())*(D[i + 1].x() - D[i].x()) - (D[i].y() * D[i + 1].x() - D[i].x() * D[i + 1].y())*(q.x() - p.x())) / ((D[i + 1].y() - D[i].y())*(q.x() - p.x()) - (D[i + 1].x() - D[i].x())*(q.y() - p.y()));
            double yp = ((p.y() * q.x() - p.x() * q.y())*(D[i + 1].y() - D[i].y()) - (D[i].y() * D[i + 1].x() - D[i].x() * D[i + 1].y())*(q.y() - p.y())) / ((D[i + 1].y() - D[i].y())*(q.x() - p.x()) - (D[i + 1].x() - D[i].x())*(q.y() - p.y()));
            Eigen::Vector3d pt_new(xp, yp, 0);
            if(is_point_on_line(pt_new, p, q) && is_point_on_line(pt_new, D[i], D[i+1])){ T.push_back(xp); }
        }
        if(T.size() == 0)return {};
        if(T.size() == 1)return{(p.x() + q.x())/2};
        t_min = *std::min_element(T.begin(), T.end()); t_min = (t_min < 0) ? 0.0: (t_min > 1)? 1.0: t_min;
        t_max = *std::max_element(T.begin(), T.end()); t_max = (t_max < 0) ? 0.0: (t_max > 1) ? 1.0:  t_max;
        std::array<Eigen::Vector3d, 2> next_line = std::array{Eigen::Vector3d(t_min, 0.0,0.0), Eigen::Vector3d(t_max, 0.0,0.0)};

        if(abs(t_max -  t_min) < 1e-9)return {t_max};
        
        std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> _bez = BezierSplit(CtrlPts_base, t_max);
        double bez_t = t_min / (t_max);
        _bez = BezierSplit(_bez.first, bez_t);
        if(abs((p - q).norm() - abs(t_max - t_min)) < 1e-9){
            std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> bez_spl = BezierSplit(_bez.second, 0.5);
            std::vector<Eigen::Vector3d> b1 = bez_spl.first, b2 = bez_spl.second;
            std::array<Eigen::Vector3d, 2> next_line2; std::copy(next_line.begin(), next_line.end(), next_line2.begin());
            std::vector<double>t1 = _bezierclipping(CtrlPts_base, b1, next_line, dim), t2 = _bezierclipping(CtrlPts_base, b2, next_line2, dim);
            for(auto&t: t2)t1.push_back(t);
            return t1;
        }

        return _bezierclipping(CtrlPts_base, _bez.second, next_line, dim);
    }

    //http://www-ikn.ist.hokudai.ac.jp/~k-sekine/slides/convexhull.pdf
    //https://kajindowsxp.com/graham-algo/
    std::vector<Eigen::Vector3d> GrahamScan(const std::vector<Eigen::Vector3d>& Q){
        std::vector<Eigen::Vector3d> S;
        if((int)Q.size() < 3)return Q;
        Eigen::Vector3d p_ml = Q[0];
        for(auto&p: Q){
            if(p_ml.y() > p.y())p_ml = p;
            else if(p_ml.y() == p.y() && p_ml.x() > p.x()) p_ml = p;
        }
        std::vector<std::pair<double, Eigen::Vector3d>>Args;
        for(auto&p: Q){
            if(p == p_ml)continue;
            double phi = atan2(p.y() - p_ml.y(), p.x() - p_ml.x());
            bool hasSameAngle = false;
            for(auto& X: Args){
                if(X.first == phi){
                    if((X.second - p_ml).norm() < (p - p_ml).norm())X.second = p;
                    hasSameAngle = true;
                }
            }
            if(!hasSameAngle)Args.push_back(std::make_pair(phi, p));
        }
        // compare only the first value
        std::sort(Args.begin(), Args.end(),[](auto const& x, auto const& y) {return x.first < y.first; });
        if((int)Args.size() >= 1)S.push_back(Args[Args.size() - 1].second);
        S.push_back(p_ml);
        Eigen::Vector3d top, next;
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

    Eigen::Vector3d calcCrossPoint_2Vector(Eigen::Vector3d p1, Eigen::Vector3d q1, Eigen::Vector3d p2, Eigen::Vector3d q2){
        double t = ((p2.x() - p1.x())*(p2.y() - q2.y()) - (p2.x() - q2.x())*(p2.y() - p1.y()))/((q1.x() - p1.x()) * (p2.y() - q2.y()) - (p2.x() - q2.x())*(q1.y() - p1.y()));
        return Eigen::Vector3d(t * (q1.x() - p1.x()) + p1.x(), t * (q1.y() - p1.y()) + p1.y(), 0);
    }

    //xy平面上に乗っていると仮定
    double SignedArea(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& p){
        Eigen::Vector3d v = a-p, v2 = b-p;
        return v.x() * v2.y() - v.y() * v2.x();
    }

    bool is_point_on_line(Eigen::Vector3d p, Eigen::Vector3d lp1, Eigen::Vector3d lp2){
        double ac = (p - lp1).norm(), bc = (p - lp2).norm(), lp = (lp1 - lp2).norm();
        return ((ac < 1e-7 || bc < 1e-7) || abs(lp - ac - bc) < 1e-7)? true: false;
    }

    //split at t for bezier curve
    std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> BezierSplit(const std::vector<Eigen::Vector3d>& CtrlPts, double t){
        std::vector<Eigen::Vector3d>lp, rp;
        std::vector<std::vector<Eigen::Vector3d>> Tree = de_casteljau_algorithm(CtrlPts, t);
        Tree.insert(Tree.begin(), CtrlPts);
        for(auto& d: Tree){
            lp.push_back(d[0]);
            rp.push_back(d[d.size() - 1]);
        }
        return {lp, rp};
    }

    std::vector<std::vector<Eigen::Vector3d>> de_casteljau_algorithm(const std::vector<Eigen::Vector3d>& CtrlPts, double t){
        std::vector<Eigen::Vector3d> Q;
        std::vector<std::vector<Eigen::Vector3d>> _Q;
        Eigen::Vector3d prev = CtrlPts[0];
        for(int i = 1; i < (int)CtrlPts.size(); i++){
            Eigen::Vector3d new_p = (1. - t) * prev + t * CtrlPts[i];
            Q.push_back(new_p);
            prev = CtrlPts[i];
        }
        if((int)Q.size() == 0)return _Q;
        if((int)Q.size() == 1)return {Q};
        _Q.push_back(Q);
        std::vector<std::vector<Eigen::Vector3d>> tmp = de_casteljau_algorithm(Q,t);
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

    double factorial(int n){ return (n > 1)? factorial(n - 1) * n : 1;}

    double cmb(int n , int i){ return factorial(n) / (factorial(i) * factorial(n - i));}

    double BernsteinBasisFunc(int n, int i, double t){return cmb(n, i) * std::pow(t, i) * std::pow(1.0 - t, n - i);}

    Eigen::Vector3d bezier(std::vector<Eigen::Vector3d>& CtrlPts, double t, int dim){
        Eigen::Vector3d v(0.,0., 0.);
        for (int i = 0; i < dim + 1; i++) v +=  BernsteinBasisFunc(dim, i, t)* CtrlPts[i];
        return v;
    }


    Eigen::Vector3d ProjectionVector(const Eigen::Vector3d& v, Eigen::Vector3d n, bool Isnormalized){
        n = n.normalized();
        Eigen::Vector3d V = v - v.dot(n)* n;
        return (Isnormalized)? V/V.norm(): V;
    }

}
