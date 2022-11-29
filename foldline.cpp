#include "foldline.h"

const double eps = 1e-5;

FoldLine::FoldLine(int crvNum, int rsize, int _type)
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
    he = he2 = nullptr;
    CtrlPts_res.clear();
    Curve_res.clear();
}

std::vector<glm::f64vec3> FoldLine::getCtrlPt(){return CtrlPts;}

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

bool FoldLine::addCtrlPt(glm::f64vec3& p, int dim, OUTLINE *outline, std::vector<Face*>& Faces, std::vector<HalfEdge*>&Edges, std::vector<Vertex*>& Vertices){
    //if((int)CtrlPts.size() <= dim)
    CtrlPts.push_back(p);
    bool hasCrv = setCurve(dim);
    /*
    else{

        if(he != nullptr){
            vx2 = new Vertex(p); he2 = new HalfEdge(vx2, EdgeType::fl);
            he2->prev = he; he->next = he2;
        }
        if(he == nullptr){
            vx = new Vertex(p); he = new HalfEdge(vx, EdgeType::fl);
        }
    }


    point.clear();
    if(hasCrv && he != nullptr && he2 != nullptr){
        std::vector<double> T = BezierClipping(CtrlPts, he, dim);
        for(auto&t: T){
            if(t != -1){
                glm::f64vec3 points = glm::f64vec3{0,0,0};
                for (int i = 0; i < int(CtrlPts.size()); i++) points += cmb(dim, i) * std::pow(t, i) * std::pow(1 - t, dim - i) * CtrlPts[i];
                point.push_back(points);
            }
        }

    }*/
    return hasCrv;
}

bool FoldLine::ChangeColor(OUTLINE *outline, int val, int dim){
    color += val;
    color = (color > 255) ? 255. : (color < -255)? -255. : color;
    //bool res = applyCurvedFolding(dim, outline);
    //return res;
}

//type == 0 line
//type == 1 b-spline
bool FoldLine::setCurve(int dim){
    int crvPtNum = CurvePts.size();
    if(type == 0){
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
    else if(type == 1){
        if((int)CtrlPts.size() <= dim) return false;
        int knotSize = (int)CtrlPts.size() + dim + 1;
        std::vector<double>T(knotSize);
        for(int j = 0; j < knotSize; j++)T[j] = (double)j/(double)knotSize;
        //端点をくっつける
        for(int j = 0;j< dim + 1;j++){ T[j] = 0; T[knotSize - 1 - j] = 1;}
        double t = T[dim];
        for(int i = 0; i < crvPtNum; i++){
            glm::f64vec3 vec = bspline(CtrlPts, t, dim, T);
            CurvePts[i] = vec;
            t += (T[(int)CtrlPts.size()] - T[dim])/(double)(crvPtNum);
        }
    }
    else if(type == 2){
        if((int)CtrlPts.size() <= dim)return false;
        double t = 0.0;
        glm::f64vec3 v;
        for(int n = 0; n < crvPtNum; n++){
            v = {0.f,0.f, 0.f};
            for (int i = 0; i < int(CtrlPts.size()); i++) {
                v += cmb(dim, i) * std::pow(t, i) * std::pow(1 - t, dim - i) * CtrlPts[i];
            }
            CurvePts[n] = v;
            t += 1/(double)crvPtNum;
        }
    }
    return true;
}

//数値微分
//http://www.slis.tsukuba.ac.jp/~fujisawa.makoto.fu/cgi-bin/wiki/index.php?4%B3%AC%C8%F9%CA%AC
//https://shuhoyo.hatenablog.com/entry/frenet-serre
//https://ja.wikipedia.org/wiki/%E3%83%95%E3%83%AC%E3%83%8D%E3%83%BB%E3%82%BB%E3%83%AC%E3%81%AE%E5%85%AC%E5%BC%8F
//https://sci.tea-nifty.com/blog/2010/08/3-8550.html
bool FoldLine::applyCurvedFolding(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, int dim){
    if(Faces.empty() || ((type == 1 || type == 2) && (int)CtrlPts.size() <= dim))return false;

    Rulings_2dL.clear(); Rulings_2dR.clear(); Rulings_3dL.clear(); Rulings_3dR.clear();
    std::vector<Vertex*> _CrvPts;
    double tau, k3d;
    int knotSize = (int)CtrlPts.size() + dim + 1;
    std::vector<double>Knot(knotSize);
    for(int j = 0; j < knotSize; j++)Knot[j] = (double)j/(double)knotSize;
    for(int j = 0;j< dim + 1;j++){ Knot[j] = 0; Knot[knotSize - 1 - j] = 1;}
    double t = (type == 1) ? Knot[dim] : 0;
    glm::f64vec3 d, dr, dr2, dr3;
    glm::f64vec3 dh_p,dh_2p, dh_m, dh_2m, dh_3p, dh_3m;
    glm::f64vec3 T, N, B;
    Vertex *vl, *vr;
    double angle = color/510. * M_PI;
    double da = 0.;//折角度は現在一定のため

    glm::f64vec3 crossPts;
    for(int i = 0; i < maxRsize; i++){
        if(type == 0){

        }else if(type == 1){
            if((int)CtrlPts.size() <= dim) return false;
            if(i == 0 || i == (int)CurvePts.size() - 1){
                t += (Knot[(int)CtrlPts.size()] - Knot[dim])/(double)(CurvePts.size());
                continue;
            }
            d = bspline(CtrlPts, t, dim, Knot);
            dh_p = bspline(CtrlPts, t + eps, dim, Knot);
            dh_2p = bspline(CtrlPts, t + 2 * eps, dim, Knot);
            dh_3p = bspline(CtrlPts, t + 3 * eps, dim, Knot);
            dh_m = bspline(CtrlPts, t - eps, dim, Knot);
            dh_2m = bspline(CtrlPts, t - 2 * eps, dim, Knot);
            dh_3m = bspline(CtrlPts, t - 3 * eps, dim, Knot);

            dr = (dh_p - dh_m)/ (2 * eps);
            dr2 = (dh_p - 2. * d + dh_m)/(eps * eps);
            dr3 = (-dh_3p + 8. * dh_2p - 13. * dh_p + 13. * dh_m - 8. * dh_2m + dh_3m)/(8. * std::pow(eps, 3));
            //glm::f64vec3 T = glm::normalize(dr), N = glm::normalize(glm::cross(dr, glm::cross(dr2, dr))), B = glm::cross(T, N);
            //k2d = glm::length(glm::cross(dr, dr2))/std::pow(glm::length(dr),3);
            tau = (glm::dot(dr,glm::cross(dr2, dr3)))/std::pow(glm::length(glm::cross(dr, dr2)), 2);

            t += (Knot[(int)CtrlPts.size()] - Knot[dim])/(double)(CurvePts.size());
        }else if(type == 2){
            int j = i * CurvePts.size() / (double)maxRsize;
            dr = -3. * (1 - t) * (1 - t) * CtrlPts[0] + (9 * t * t - 12 * t + 3) * CtrlPts[1] + (-9 * t * t + 6 * t) * CtrlPts[2] + 3 * t * t * CtrlPts[3];
            dr2 = 6. * (1 - t) * CtrlPts[0] + (18 * t - 12) * CtrlPts[1] + (-18 * t + 6) * CtrlPts[2] + 6 * t * CtrlPts[3];
            dr3 = -6. * (CtrlPts[0] - CtrlPts[3]) + 18. * (CtrlPts[1] - CtrlPts[2]);
            T = glm::normalize(dr);
            N = glm::normalize(glm::cross(dr, glm::cross(dr2, dr)));
            B = glm::normalize(glm::cross(dr, dr2));
            if(glm::dot(glm::f64vec3{0,0,-1}, B) < 0){B *= -1; N *= -1;}
            k3d = glm::length(glm::cross(dr, dr2))/std::pow(glm::length(dr), 3);
            tau = glm::dot(dr, glm::cross(dr2, dr3))/std::pow(glm::length(glm::cross(dr, dr2)), 2);
            //k2d = k3d * cos(angle);
            double phi_bl = (std::isfinite(k3d*sin(angle)/tau)) ? atan(tan(k3d*sin(angle)/tau)) : M_PI/2, phi_br = (std::isfinite(k3d*sin(angle)/tau)) ? atan(tan(k3d*sin(angle)/tau)) : M_PI/2;

            glm::f64vec3 rr3d = glm::normalize(cos(phi_br) * T + sin(phi_br) * cos(angle) * N + sin(phi_br) * sin(angle) * B);
            glm::f64vec3 rl3d = glm::normalize(cos(phi_bl) * T - sin(phi_bl) * cos(angle) * N + sin(phi_bl) * sin(angle) * B);
            glm::f64vec3 rl2d = glm::normalize(glm::rotate(phi_bl, glm::f64vec3{0,0,1}) * glm::f64vec4{T,1});
            glm::f64vec3 rr2d = glm::normalize(glm::rotate(-phi_br, glm::f64vec3{0,0,1}) * glm::f64vec4{T,1});

            bool res = false, res2 = false;
            res= setPoint(Faces, rl2d, CurvePts[j],crossPts);
            if(res){
                double l = glm::distance(crossPts, CurvePts[j]);
                vl = new Vertex(crossPts); vl->p3 = rl3d * l + CurvePts[j];
                Vertices.push_back(vl);
                splitOnPoint(vl, Edges);
                //Rulings_2dL.push_back(std::array{crossPts, CurvePts[j]});
                //Rulings_3dL.push_back(std::array{rl3d * l + CurvePts[j], CurvePts[j]});
                //std::cout << "left  " <<  glm::to_string(crossPts) << " " << glm::to_string(rl3d * l + CurvePts[j]) << "  " << glm::to_string(rl2d) << "  " << glm::to_string(rl3d) << std::endl;
            }
            res2 = setPoint(Faces, rr2d, CurvePts[j],crossPts);
            if(res2){
                double l = glm::distance(crossPts, CurvePts[j]);
                vr = new Vertex(crossPts); vr->p3 = rr3d * l + CurvePts[j];
                Vertices.push_back(vr);
                splitOnPoint(vr, Edges);
                //Rulings_2dR.push_back(std::array{crossPts, CurvePts[j]});
                //Rulings_3dR.push_back(std::array{rr3d * l + CurvePts[j], CurvePts[j]});
                //std::cout << "right  " << glm::to_string(crossPts) << " " << glm::to_string(rr3d * l + CurvePts[j]) <<  "  " <<glm::to_string(rr2d) << "  " << glm::to_string(rr3d)  << std::endl;
            }
            if(res && res2){
                devide(vl, vr, Faces, Edges, EdgeType::r);
                Vertex *v = new Vertex(CurvePts[j]);
                _CrvPts.push_back(v);
                Vertices.push_back(v);
                splitOnPoint(v, Edges);
            }
            t += 1/(double)maxRsize;
        }
    }
    int Csize = _CrvPts.size();
    for(int i = 1; i < Csize; i++){
        devide(_CrvPts[i - 1], _CrvPts[i], Faces, Edges, EdgeType::fl);
    }
    return true;
}

void FoldLine::deform(){

}

bool FoldLine::setPoint(std::vector<Face*>& Faces, glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint){
    double l = 1000;
    glm::f64vec3 N0 = l * N + cp;
    glm::f64vec3 v, v2;
    bool IsIntersected = false;
    double minDist = 1000;

    for(auto&f: Faces){
        if(!cn(f, cp))continue;
        HalfEdge *h = f->halfedge;
        do{
            v = h->vertex->p;
            v2 = h->next->vertex->p;
            if (IsIntersect(v, v2, N0, cp)) {
                glm::f64vec3 p = getIntersectionPoint(v, v2, N0, cp);
                double dist = glm::distance(p, cp);
                if(dist < minDist){
                    crossPoint = p;
                    minDist = dist;
                }
                IsIntersected = true;
            }
            h = h->next;
        }while(h != f->halfedge);
    }

    return IsIntersected;
}

#include <fstream>
#include <sstream>
bool FoldLine::modify2DRulings(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices,int dim){
    glm::f64vec3 d, dr, dr2, dr3, T, N, B;
    glm::f64vec3 dh_p, dh_m;
    double k3d, tau, k3d_bef;
    std::vector<glm::f64vec3> T2d;
    std::vector<double> K2D, K2D_bef;
    auto rad_2d = [](double k, double tau, double a, double da) { return atan2(k*sin(a),(da + tau)); };
    std::vector<HalfEdge*> SearchedEdge;
    CrossPts.clear();
    if(type == 0){

    }else if(type == 2){
        std::string file = "result.csv";
        std::ofstream ofs(file);
        double a, a_bef;
        int EdgeSize = Edges.size();
        for(int i = 0; i < EdgeSize; i++){
            if(Edges[i]->r == nullptr || std::find(SearchedEdge.begin(), SearchedEdge.end(), Edges[i]) != SearchedEdge.end() || std::find(SearchedEdge.begin(), SearchedEdge.end(), Edges[i]->pair) != SearchedEdge.end())continue;
            SearchedEdge.push_back(Edges[i]);
            std::vector<double>arcT = BezierClipping(CtrlPts, Edges[i], dim);
            for(auto&t: arcT){
                if(t < 0 || 1 < t){std::cout<<" t is not found " << std::endl; continue;}
                glm::f64vec3 pt_new{0,0,0};
                for (int i = 0; i < int(CtrlPts.size()); i++) pt_new += cmb(dim, i) * std::pow(t, i) * std::pow(1 - t, dim - i) * CtrlPts[i];               
                if(!is_point_on_line(pt_new, Edges[i]->vertex->p, Edges[i]->next->vertex->p)) continue;
                double sa = glm::distance(pt_new, Edges[i]->vertex->p), sc = glm::distance(Edges[i]->vertex->p, Edges[i]->next->vertex->p);
                Vertex *v = new Vertex(pt_new); v->p3= sa/sc * (Edges[i]->next->vertex->p3 - Edges[i]->vertex->p3) + Edges[i]->vertex->p3;
                CrossPts.push_back(v);

                dr = -3. * (1 - t) * (1 - t) * CtrlPts[0] + (9 * t * t - 12 * t + 3) * CtrlPts[1] + (-9 * t * t + 6 * t) * CtrlPts[2] + 3 * t * t * CtrlPts[3];
                dr2 = 6. * (1 - t) * CtrlPts[0] + (18 * t - 12) * CtrlPts[1] + (-18 * t + 6) * CtrlPts[2] + 6 * t * CtrlPts[3];
                K2D.push_back(glm::length(glm::cross(dr, dr2))/std::pow(glm::length(dr), 3));
                T2d.push_back(glm::normalize(dr));
                double t2 = t - eps;
                glm::f64vec3 dr_bef = -3. * (1 - t2) * (1 - t2) * CtrlPts[0] + (9 * t2 * t2 - 12 * t2 + 3) * CtrlPts[1] + (-9 * t2 * t2 + 6 * t2) * CtrlPts[2] + 3 * t2 * t2 * CtrlPts[3];
                glm::f64vec3 dr2_bef = 6. * (1 - t2) * CtrlPts[0] + (18 * t2 - 12) * CtrlPts[1] + (-18 * t2 + 6) * CtrlPts[2] + 6 * t2 * CtrlPts[3];
                K2D_bef.push_back(glm::length(glm::cross(dr_bef, dr2_bef))/std::pow(glm::length(dr_bef), 3));
            }
        }
        std::vector<double>Knot;
        Curve_res = GlobalSplineInterpolation(CrossPts, CtrlPts_res, Knot);
        for(int j = 0; j < (int)CrossPts.size(); j++){
            int ind = -1;
            double dist = 1;
            for(int i = 0; i < (int)Curve_res.size(); i++){
                if(glm::distance(Curve_res[i], CrossPts[j]->p3) < dist){
                    dist = glm::distance(Curve_res[i], CrossPts[j]->p3); ind = i;
                }
            }
            if(ind < 3) ind = 3;
            else if(ind + 4 > (int)Curve_res.size())ind = (int)Curve_res.size() - 4;
            double t = Knot[dim] + ind*(Knot[CtrlPts_res.size()] - Knot[dim])/Curve_res.size();
            d = bspline(CtrlPts_res, t, dim, Knot); dh_p = bspline(CtrlPts_res, t + eps, dim, Knot); dh_m = bspline(CtrlPts_res, t - eps, dim, Knot);
            //dh_2p = bspline(CtrlPts_res, t + 2 * eps, dim, Knot); dh_2m = bspline(CtrlPts_res, t - 2 * eps, dim, Knot);
            //dh_3p = bspline(CtrlPts_res, t + 3 * eps, dim, Knot); bspline(CtrlPts_res, t - 3 * eps, dim, Knot);
            dr = (dh_p - dh_m)/ (2 * eps);
            dr2 = (dh_p - 2. * d + dh_m)/(eps * eps);
            dr3 = 3. * (dh_p - 2. * eps * dr - dh_m)/std::pow(eps, 3);
            T = glm::normalize(dr);
            N = glm::normalize(glm::cross(dr, glm::cross(dr2, dr)));
            B = glm::normalize(glm::cross(dr, dr2));
            if(glm::dot(glm::f64vec3{0,0,-1}, B) < 0){B *= -1; N *= -1;}
            //std::cout << glm::to_string(T) << ", " << glm::to_string(N) << ", " << glm::to_string(B) << " , " << glm::dot(dr, glm::cross(dr2, dr3)) << " ,  "<<std::pow(glm::length(glm::cross(dr, dr2)), 2) <<std::endl;
            tau = glm::dot(dr, glm::cross(dr2, dr3))/std::pow(glm::length(glm::cross(dr, dr2)), 2);
            k3d = glm::length(glm::cross(dr, dr2))/std::pow(glm::length(dr), 3);
            k3d = (K2D[j] > k3d) ? K2D[j] : k3d;
            a = (k3d != 0) ? acos(K2D[j]/k3d): 0;

            t -= eps;
            d = bspline(CtrlPts_res, t, dim, Knot); dh_p = bspline(CtrlPts_res, t + eps, dim, Knot); dh_m = bspline(CtrlPts_res, t - eps, dim, Knot);
            glm::f64vec3 dr_bef = (dh_p - dh_m)/ (2 * eps);
            glm::f64vec3 dr2_bef = (dh_p - 2. * d + dh_m)/(eps * eps);
            k3d_bef = glm::length(glm::cross(dr_bef, dr2_bef))/std::pow(glm::length(dr_bef), 3);
            k3d_bef = (K2D_bef[j] > k3d_bef) ? K2D_bef[j]: k3d_bef;
            a_bef = (k3d_bef != 0) ? acos(K2D_bef[j]/k3d_bef): 0;
            double da = (a - a_bef)/eps;
            double phi_bl = rad_2d(k3d, tau, a, -da), phi_br = rad_2d(k3d, tau, a, da);

            glm::f64vec3 rr3d = glm::normalize(cos(phi_br) * T + sin(phi_br) * cos(a) * N + sin(phi_br) * sin(a) * B);
            glm::f64vec3 rl3d = glm::normalize(cos(phi_bl) * T - sin(phi_bl) * cos(a) * N + sin(phi_bl) * sin(a) * B);
            glm::f64vec3 rl2d = glm::normalize(glm::rotate(phi_bl, glm::f64vec3{0,0,1}) * glm::f64vec4{T2d[j],1});
            glm::f64vec3 rr2d = glm::normalize(glm::rotate(phi_br, glm::f64vec3{0,0,1}) * glm::f64vec4{T2d[j],1});
            //std::cout<<glm::to_string(rl2d) << " , " << glm::to_string(rr2d)<<std::endl;
            std::array<HalfEdge*, 2> he_new = splitOnPoint(CrossPts[j], Edges);
            /*
            glm::f64vec3 pt_new;
            bool res = setPoint(Faces, rl2d, CrossPts[j]->p,pt_new);
            if(res){
                double l = glm::distance(pt_new, CrossPts[j]->p3);

                if(glm::dot(rl2d, glm::normalize(he_new[0]->next->vertex->p - pt_new)) > 0){
                    he_new[0]->next->vertex->p3 = rl3d * l + CrossPts[j]->p3;
                    he_new[0]->next->vertex->p = pt_new;
                }else{
                    he_new[1]->next->vertex->p3 = rl3d * l + CrossPts[j]->p3;
                    he_new[1]->next->vertex->p = pt_new;
                }
            }
            bool res2 = setPoint(Faces, rr2d, CrossPts[j]->p,pt_new);
            if(res2){
                double l = glm::distance(CrossPts[j]->p3, pt_new);
                if(glm::dot(rr2d, glm::normalize(he_new[1]->next->vertex->p - pt_new)) > 0){
                    he_new[1]->next->vertex->p3 = rr3d * l + CrossPts[j]->p3;
                    he_new[1]->next->vertex->p = pt_new;
                }else{
                    he_new[0]->next->vertex->p3 = rr3d * l + CrossPts[j]->p3;
                    he_new[0]->next->vertex->p = pt_new;
                }
            }*/
            ofs << k3d << " , " << k3d_bef << ", "<<K2D[j] << " , " << tau << " , " << a << ", " << phi_bl << " , " << phi_br<<  std::endl;
        }

    }

    return false;
}

bool FoldLine::BezierCrvOn3dSrf(std::vector<glm::f64vec3>& CtrlPts, double t, int dim, std::vector<Face*>& Faces, glm::f64vec3& v_3d){
    glm::f64vec3 v_2d{0,0,0};
    for (int i = 0; i < int(CtrlPts.size()); i++) v_2d += cmb(dim, i) * std::pow(t, i) * std::pow(1 - t, dim - i) * CtrlPts[i];

    for(auto&f: Faces){
        if(cn(f,v_2d)){
            HalfEdge *h = f->halfedge;
            glm::f64vec3 v1 = h->next->vertex->p - h->vertex->p;
            glm::f64vec3 v2 = h->prev->vertex->p - h->vertex->p;
            glm::f64vec3 p = v_2d - h->vertex->p;
            Eigen::Matrix2d A; A(0,0) = v1.x; A(1,0) = v1.y; A(0,1) = v2.x; A(1,1) = v2.y;
            Eigen::Vector2d b; b(0) = p.x; b(1) = p.y;
            Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
            v_3d = x(0) * (h->next->vertex->p3 - h->vertex->p3) + x(1) * (h->prev->vertex->p3 - h->vertex->p3) + h->vertex->p3;
            //std::cout<< glm::length(h->next->vertex->p3 - h->next->vertex->p) << " , " << glm::length(h->prev->vertex->p3 - h->prev->vertex->p) << " , " <<  glm::length(h->vertex->p3 - h->vertex->p) << std::endl;
            return true;
        }
    }
    std::cout<<"not found in BezierCrvOn3dSrf"<<std::endl;
    return false;
}

std::array<HalfEdge*, 2> FoldLine::splitOnPoint(Vertex *v, std::vector<HalfEdge*>& Edges){
    std::array<HalfEdge*, 2> he_new{nullptr, nullptr};
    for(auto&h: Edges){
        if(is_point_on_line(v->p, h->vertex->p, h->next->vertex->p)){
            double t = glm::length(v->p - h->vertex->p)/glm::length(h->next->vertex->p - h->vertex->p);
            v->p3 = h->vertex->p3 + t * (h->next->vertex->p3 - h->vertex->p3);
            EdgeType _type = h->edgetype;
            HalfEdge *h_new = new HalfEdge(v, _type); h_new->r = h->r;
            if(h->pair != nullptr){
                HalfEdge *_pair = h->pair;
                HalfEdge *h2_new = new HalfEdge(v, _type); h2_new->r = h->r;
                h_new->face = h->face; h2_new->face = _pair->face;
                h_new->next = h->next; h_new->prev = h; h->next = h_new;
                h2_new->next = _pair->next; h2_new->prev = _pair; _pair->next = h2_new;
                h2_new->pair = h; h_new->pair = _pair; _pair->pair = h_new; h->pair = h2_new;
                Edges.push_back(h_new); Edges.push_back(h2_new);
                he_new[0] = h_new; he_new[1] = h2_new;
            }else{
                std::cout <<"pair not found in splitOnPoint"<<std::endl;
                h_new->face = h->face; h_new->next = h->next; h_new->prev = h; h->next = h_new;
                Edges.push_back(h_new);
            }
            return he_new;
        }
    }

    return he_new;
}

void FoldLine::devide2Faces(std::vector<HalfEdge*>& Inserted, std::vector<HalfEdge*>& Edges, std::vector<Face*>& Faces){
    std::vector<Face*> F_devided;
    for(int i = 0; i < (int)Inserted.size(); i++){
        if(std::find(F_devided.begin(),F_devided.end(), Inserted[i]->face) != F_devided.end())continue;
        HalfEdge *he = Inserted[i];
        HalfEdge *he_pair = nullptr;
        do{
            if(he != Inserted[i] && std::find(Inserted.begin(), Inserted.end(), he) != Inserted.end()){
                he_pair = he;
                break;
            }
            he = he->next;
        }while(he != Inserted[i]);
        if(he_pair == nullptr)continue;
        HalfEdge *h_new = new HalfEdge(Inserted[i]->vertex, EdgeType::fl);
        HalfEdge *h2_new = new HalfEdge(he_pair->vertex, EdgeType::fl);
        h_new->prev = Inserted[i]->prev; h_new->next = he_pair;
        h2_new->prev = he_pair->prev; h2_new->next = Inserted[i];
        Inserted[i]->prev = h2_new; he_pair->prev = h_new;
        h_new->pair = h2_new; h2_new->pair = h_new;
        HalfEdge *h = Inserted[i];
        do{
            h->face = he->face;
            h = h->next;
        }while(h != Inserted[i]);
        Face *f_new = new Face(he_pair);
        h = he_pair;
        do{
            h->face = f_new;
            h = h->next;
        }while(h != he_pair);
        Edges.push_back(h_new); Edges.push_back(h2_new); Faces.push_back(f_new);
        F_devided.push_back((Inserted[i]->face)); F_devided.push_back(f_new);
    }

}

void FoldLine::devide(Vertex *v1, Vertex *v2, std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, EdgeType _type){
    for(auto&f: Faces){
        HalfEdge *h = f->halfedge, *h1 = nullptr, *h2 = nullptr;
        do{
            if(h->vertex == v1)h1 = h;
            if(h->vertex == v2)h2 = h;
            h = h->next;
        }while(h != f->halfedge);
        if(h1 != nullptr && h2 != nullptr){
            HalfEdge *h1_new = new HalfEdge(h1->vertex, _type), *h2_new = new HalfEdge(h2->vertex, _type);
            h2_new->pair = h1_new; h1_new->pair = h2_new;

            h1_new->prev = h1->prev; h1_new->next = h2;
            h2_new->prev = h2->prev; h2_new->next = h1;
            h1->prev->next = h1_new; h2->prev->next = h2_new;
            h1->prev = h2_new; h2->prev = h1_new;

            HalfEdge *h = h1;
            f->halfedge = h1;
            do{
                h->face = f;
                h = h->next;
            }while(h != h1);
            h = h2;
            Face* face2 = new Face(h2);
            Faces.push_back(face2);
            face2->halfedge = h2;
            do{
                h2->face = face2;
                h = h->next;
            }while(h != h2);
            Edges.push_back(h1_new); Edges.push_back(h2_new);
            return;
        }
    }
}

