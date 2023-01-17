#include "foldline.h"

const double eps = 1e-5;
using namespace MathTool;

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
    he = he2 = nullptr;
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

bool FoldLine::addCtrlPt(glm::f64vec3& p, int dim){
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

//数値微分
//http://www.slis.tsukuba.ac.jp/~fujisawa.makoto.fu/cgi-bin/wiki/index.php?4%B3%AC%C8%F9%CA%AC
//https://shuhoyo.hatenablog.com/entry/frenet-serre
//https://ja.wikipedia.org/wiki/%E3%83%95%E3%83%AC%E3%83%8D%E3%83%BB%E3%82%BB%E3%83%AC%E3%81%AE%E5%85%AC%E5%BC%8F
//https://sci.tea-nifty.com/blog/2010/08/3-8550.html
bool FoldLine::applyCurvedFolding(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, int dim){
    using namespace MathTool;
    if(Faces.empty() || (type == PaintTool::FoldLine_bezier && (int)CtrlPts.size() <= dim))return false;

    Rulings_2dL.clear(); Rulings_2dR.clear(); Rulings_3dL.clear(); Rulings_3dR.clear();
    std::vector<Vertex*> _CrvPts;
    double tau, k3d;
    int knotSize = (int)CtrlPts.size() + dim + 1;
    std::vector<double>Knot(knotSize);
    for(int j = 0; j < knotSize; j++)Knot[j] = (double)j/(double)knotSize;
    for(int j = 0;j< dim + 1;j++){ Knot[j] = 0; Knot[knotSize - 1 - j] = 1;}
    double t =  0;
    glm::f64vec3 d, dr, dr2, dr3;
    glm::f64vec3 dh_p,dh_2p, dh_m, dh_2m, dh_3p, dh_3m;
    glm::f64vec3 T, N, B;
    Vertex *vl, *vr;
    double angle = color/510. * M_PI;
    //double da = 0.;//折角度は現在一定のため

    glm::f64vec3 crossPts;
    for(int i = 0; i < maxRsize; i++){
        if(type == PaintTool::FoldLine_line){

        }else if(type == PaintTool::FoldLine_arc){
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
        }else if(type == PaintTool::FoldLine_bezier){
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
            //res= setPoint(Faces, rl2d, CurvePts[j],crossPts);
            if(res){
                double l = glm::distance(crossPts, CurvePts[j]);
                vl = new Vertex(crossPts); vl->p3 = rl3d * l + CurvePts[j];
                Vertices.push_back(vl);
                std::vector<HalfEdge*> H_new;
                for(auto&h : Edges){
                    H_new = h->Split(vl,Edges);
                    if(!H_new.empty())break;
                }
                //Rulings_2dL.push_back(std::array{crossPts, CurvePts[j]});
                //Rulings_3dL.push_back(std::array{rl3d * l + CurvePts[j], CurvePts[j]});
                //std::cout << "left  " <<  glm::to_string(crossPts) << " " << glm::to_string(rl3d * l + CurvePts[j]) << "  " << glm::to_string(rl2d) << "  " << glm::to_string(rl3d) << std::endl;
            }
            //res2 = setPoint(Faces, rr2d, CurvePts[j],crossPts);
            if(res2){
                double l = glm::distance(crossPts, CurvePts[j]);
                vr = new Vertex(crossPts); vr->p3 = rr3d * l + CurvePts[j];
                Vertices.push_back(vr);
                std::vector<HalfEdge*> H_new;
                for(auto&h : Edges){
                    H_new = h->Split(vr,Edges);
                    if(!H_new.empty())break;
                }
                //Rulings_2dR.push_back(std::array{crossPts, CurvePts[j]});
                //Rulings_3dR.push_back(std::array{rr3d * l + CurvePts[j], CurvePts[j]});
                //std::cout << "right  " << glm::to_string(crossPts) << " " << glm::to_string(rr3d * l + CurvePts[j]) <<  "  " <<glm::to_string(rr2d) << "  " << glm::to_string(rr3d)  << std::endl;
            }
            if(res && res2){
                devide(vl, vr, Faces, Edges, EdgeType::r);
                Vertex *v = new Vertex(CurvePts[j]);
                _CrvPts.push_back(v);
                Vertices.push_back(v);
                std::vector<HalfEdge*> H_new;
                for(auto&h : Edges){
                    H_new = h->Split(v,Edges);
                    if(!H_new.empty())break;
                }

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

bool FoldLine::setPoint(const std::vector<glm::f64vec3>& edge_outline, glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint){
    double l = 1000;
    bool IsIntersected = false;
    double minDist = 1000;
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    int n = edge_outline.size();
    for(int i = 0; i < n; i++){
        glm::f64vec3 v = edge_outline[(i+1) % n] - edge_outline[i];
        A(0,0) = v.x; A(0,1) = -l*N.x; A(1,0) = v.y; A(1,1) = -l * N.y;
        b(0) = cp.x - edge_outline[i].x; b(1) = cp.y - edge_outline[i].y;
        if(abs(glm::dot(glm::normalize(v), N))>= 1-FLT_EPSILON)continue;
        Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
        if (0 <= x(1) && x(1) <= 1) {
            glm::f64vec3 p = x(1)*l*N + cp;
            double dist = glm::distance(p, cp);
            if(dist < minDist){
                crossPoint = p;
                minDist = dist;
            }
            IsIntersected = true;
        }
    }
    return IsIntersected;
}



inline double FoldLine::rad_2d(double k, double tau, double a, double da){
    double x = k*sin(a), y = da + tau;
    //if(abs(x) < FLT_EPSILON && abs(y) < FLT_EPSILON)return 0;
    return M_PI/2 - atan2(y, x);
}
bool FoldLine::modify2DRulings(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, const std::vector<HalfEdge*>& edge_outline, int dim){
    using namespace MathTool;
    //auto diff1 =[](std::vector<glm::f64vec3>& P, double t){return -3. * std::pow(1.0 - t,2) * P[0] + (9.0*t*t - 12*t + 3.) * P[1] + (-9.0* t*t + 6.*t) * P[2] + 3.0*std::pow(t,2)*P[3];};
    //auto diff2 =[](std::vector<glm::f64vec3>& P, double t){return 6. * (1. - t) * P[0] + (18. * t - 12.) * P[1] + (-18. * t + 6.) * P[2] + 6. * t * P[3];};
    //auto diff3 =[](std::vector<glm::f64vec3>& P){return -6.0*P[0]+ 18.0*P[1] -18.0*P[2] + 6.0*P[3];};
    //auto diffPts = [diff1, diff2, diff3](std::vector<glm::f64vec3>& P, double t){ return std::vector<glm::f64vec3>{diff1(P,t), diff2(P,t), diff3(P)};};
    auto setface = [](HalfEdge *he, Face *f){HalfEdge *h = he; do{h->face = f; h = h->next;}while(he != h);};
    auto VectorDigAlign = [](glm::f64vec3& v, int dig = 13){
        glm::f64vec3 vi = glm::f64vec3{(int)v.x, (int)v.y, (int)v.z};
        glm::f64vec3 vf = v - vi;
        int n = pow(10, dig); int x = vf.x * n, y = vf.y * n, z = vf.z * n;
        v = glm::f64vec3{vi.x + (double)x/(double)n, vi.y + (double)y/(double)n, vi.z + (double)z/(double)n}; };

    std::string file = "result.csv"; std::ofstream resultfile(file);
    std::string file2 = "curvepoint.csv"; std::ofstream curvepointfile(file2);
    curvepointfile << "t, kappa2d,  kappa, tau, CtrlPtsx, CtrlPtsy, CtrlPtsz, Px, Py, Pz,Tx, Ty, Tz, Nx, Ny, Nz, Bx, By, Bz"<<std::endl;
    resultfile << "k3d , k2d , tau , a, da,  phi_bl , phi_br, length(p - p3), k2d_bef, k3d_bef, k2d_next, k3d_next" << std::endl;

    std::vector<glm::f64vec3>edges_ol;
    std::vector<HalfEdge*> edges_he;
    for(auto&e: edge_outline){
        glm::f64vec3 p = e->vertex->p;
        edges_he.push_back(new HalfEdge(new Vertex(p), EdgeType::ol));
        edges_ol.push_back(p);
    }
    for(int i = 0; i < (int)edges_he.size(); i++)edges_he[i]->next = edges_he[(i + 1) %(int)edges_he.size()];
    Face *face_ol = new Face(edges_he[0]);

    std::vector<HalfEdge*> SearchedEdge;  
    if(type == PaintTool::FoldLine_line){
    }else if(type == PaintTool::FoldLine_arc){
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
            CrvPt_FL P = CrvPt_FL{v, v3, 0};
            P.k2d = P.k2d_m = P.k2d_p = 1/glm::distance(CtrlPts[0], CtrlPts[1]);
            T_crs.push_back(P);

        }
    }else if(type == PaintTool::FoldLine_bezier){
        for(auto& e: Edges){
            if((e->edgetype == EdgeType::r && (std::find(SearchedEdge.begin(), SearchedEdge.end(), e) != SearchedEdge.end() || std::find(SearchedEdge.begin(), SearchedEdge.end(), e->pair) != SearchedEdge.end()))
                    || e->edgetype != EdgeType::r)continue;
            SearchedEdge.push_back(e);
            std::vector<double>arcT = BezierClipping(CtrlPts, e, dim);
            for(auto&t: arcT){
                if(t < 0 || 1 < t){std::cout<<"t is not correct value " << t << std::endl; continue;}
                glm::f64vec3 v2{0,0,0};
                for (int i = 0; i < int(CtrlPts.size()); i++) v2 += MathTool::BernsteinBasisFunc(dim, i, t) * CtrlPts[i];
                if(!is_point_on_line(v2, e->vertex->p, e->next->vertex->p)) continue;
                double sa = glm::distance(v2, e->vertex->p), sc = glm::distance(e->vertex->p, e->next->vertex->p);
                glm::f64vec3 v3 = sa/sc * (e->next->vertex->p3 - e->vertex->p3) + e->vertex->p3;
                VectorDigAlign(v2); VectorDigAlign(v3);
                CrvPt_FL P(v2, v3, t);
                T_crs.push_back(P);
            }
        }
    }
    std::vector<double> Knot, Knot2d;
    std::vector<glm::f64vec3> dP(3);
    std::sort(T_crs.begin(), T_crs.end());

    double CurveLen3d = 0.0,  CurveLen2d = 0.0;
    Curve_res = GlobalSplineInterpolation(T_crs, CtrlPts_res, Knot, CurveLen3d, true, dim);
    GlobalSplineInterpolation(T_crs, CtrlPts_res2d, Knot2d, CurveLen2d, false, dim);
    for(auto&t: T_crs){
        int index = &t - &T_crs[0];
         diff(t.s, Knot, dP, CtrlPts_res, index, dim);
         t.T3d = glm::normalize(dP[0]);
         t.N3d = glm::normalize(glm::cross(dP[0], glm::cross(dP[1], dP[0])));
         t.B3d = glm::normalize(glm::cross(dP[0], dP[1]));
         t.k3d = glm::length(glm::cross(dP[0], dP[1]))/std::pow(glm::length(dP[0]), 3); //t.digAlign(t.k3d);
         t.tau = glm::dot(dP[0], glm::cross(dP[1], dP[2]))/std::pow(glm::length(glm::cross(dP[0], dP[1])), 2); // t.digAlign(t.tau);
         diff(t.s, Knot2d, dP, CtrlPts_res2d, index, dim);
         glm::f64vec3 T = glm::normalize(dP[0]), N = glm::f64vec3{-T.y, T.x, 0}, B = glm::cross(T, N);
         t.k2d = glm::length(glm::cross(dP[0], dP[1]))/std::pow(glm::length(dP[0]), 3); //t.digAlign(t.k2d);
         t.T2d = T;
         t.N2d = N;
         t.B2d = B;
         diff(t.s - eps, Knot2d, dP, CtrlPts_res2d, index, dim); t.k2d_m = glm::length(glm::cross(dP[0], dP[1]))/std::pow(glm::length(dP[0]), 3);// t.digAlign(t.k2d_m);
         diff(t.s + eps, Knot2d, dP, CtrlPts_res2d, index, dim); t.k2d_p = glm::length(glm::cross(dP[0], dP[1]))/std::pow(glm::length(dP[0]), 3); //t.digAlign(t.k2d_p);

        //弧長パラメータ上で中心差分法をおこなってみる -　今は０から１の範囲のｔ
         curvepointfile << t.s << ", " << t.k2d <<" , " << t.k3d << ", " << t.tau <<", " << t.p.x <<", " <<t.p.y <<", " <<t.p.z  <<", "<< t.p3.x << ", " <<t.p3.y <<", "<< t.p3.z <<", " << t.T3d.x << ", " << t.T3d.y <<", " << t.T3d.z <<", " << t.N3d.x << ", " << t.N3d.y<<", " << t.N3d.z << ", " << t.B3d.x << ", " << t.B3d.y << ", " << t.B3d.z << std::endl;
         if(index == 0){
             for(auto&edge: Edges){
                 if(is_point_on_line(t.p,edge->vertex->p, edge->next->vertex->p)){
                     glm::f64vec3 f_front = edge->face->getNormalVec();
                     if(edge->edgetype == EdgeType::r){
                         glm::f64vec3 f_front2 = edge->pair->face->getNormalVec();
                         f_front = glm::normalize((f_front + f_front2)/2.0);
                     }std::swap(t.B3d, t.N3d);
                     if(glm::dot(t.B3d, f_front) < 0)t.B3d *= -1;
                     if(abs(glm::dot(t.B3d, f_front)) < abs(glm::dot(t.N3d, f_front))){glm::f64vec3 tmp = t.B3d; t.B3d = t.N3d; t.N3d = tmp;}
                 }
             }
             if(glm::dot(glm::normalize(glm::cross(t.B3d, t.T3d)), t.N3d) < 0) t.N3d *= -1;
         }
         if(index > 0 && glm::dot(T_crs[index - 1].N3d, t.N3d) < 0){ t.N3d *= -1;}
         if(index > 0 && glm::dot(T_crs[index - 1].B3d, t.B3d) < 0){ t.B3d *= -1;}

         t.k3d = (t.k2d > t.k3d) ? t.k2d : t.k3d;
         t.a = (t.k3d != 0) ? acos((long double)t.k2d/(long double)t.k3d): 0;

         diff(t.s - eps, Knot, dP, CtrlPts_res, index, dim);
         t.k3d_m = glm::length(glm::cross(dP[0], dP[1]))/std::pow(glm::length(dP[0]), 3); t.k3d_m = (t.k2d_m > t.k3d_m) ? t.k2d_m : t.k3d_m;//t.digAlign(t.k3d_m); t.k3d_m = (t.k2d_m > t.k3d_m) ? t.k2d_m : t.k3d_m;
         diff(t.s + eps, Knot, dP, CtrlPts_res, index, dim);
         t.k3d_p = glm::length(glm::cross(dP[0], dP[1]))/std::pow(glm::length(dP[0]), 3); t.k3d_p = (t.k2d_p > t.k3d_p) ? t.k2d_p : t.k3d_p;//t.digAlign(t.k3d_p); t.k3d_p = (t.k2d_p > t.k3d_p) ? t.k2d_p : t.k3d_p;
         double a_m = (t.k3d_m != 0) ? acos((long double)t.k2d_m/(long double)t.k3d_m): 0;
         double a_p = (t.k3d_p != 0) ? acos((long double)t.k2d_p/(long double)t.k3d_p): 0;
         t.da = (a_p - a_m)/(2.0*eps);
    }
    //T_crs[0].da = T_crs[1].da; T_crs[T_crs.size()-1].da = T_crs[T_crs.size()-2].da;

    //LSM_apply_tau();
    for(auto&t: T_crs){
        double phi_bl = rad_2d(t.k3d, t.tau, t.a, -t.da), phi_br = rad_2d(t.k3d, t.tau, t.a, t.da);
         glm::f64vec3 rr3d = glm::normalize(cos(phi_br) * t.T3d + sin(phi_br) * cos(t.a) * t.N3d + sin(phi_br) * sin(t.a) * t.B3d);
         glm::f64vec3 rl3d = glm::normalize(cos(phi_bl) * t.T3d - sin(phi_bl) * cos(t.a) * t.N3d + sin(phi_bl) * sin(t.a) * t.B3d);
         glm::f64vec3 rl2d = (glm::rotate(phi_bl, t.B2d) * glm::f64vec4{t.T2d,1});
         glm::f64vec3 rr2d = (glm::rotate(-phi_br, t.B2d) * glm::f64vec4{t.T2d,1});

         resultfile <<std::setprecision(15) << t.k3d << ", "<<t.k2d << " , " << t.tau << " , " << t.a << " , " << t.da << ", " << phi_bl * 180/M_PI << " , " << phi_br * 180/M_PI  <<", " << glm::length(t.p - t.p3)<< ", " << t.k2d_m << ", " <<t.k3d_m << ", " << t.k2d_p << ", "<< t.k3d_p << std::endl;

         //Rulings_2dL.push_back({rl2d, t.p}); Rulings_2dR.push_back({rr2d, t.p}); Rulings_3dL.push_back({rl3d, t.p3}); Rulings_3dR.push_back({rr3d, t.p3});
         std::vector<HalfEdge*> H_new;

         for(auto&he: Edges){
             H_new = he->Split(&t,Edges);
             if(H_new.empty())continue;
             if(H_new.size() >0){
                 HalfEdge *h = H_new[0];
                 glm::f64vec3 crossptl, crossptr;
                 setPoint(edges_ol, rl2d, t.p, crossptl);setPoint(edges_ol, rr2d, t.p, crossptr);
                 double ll = glm::distance(crossptl, t.p), lr = glm::distance(crossptr, t.p);
                 if(glm::dot(rl2d, glm::normalize(h->next->vertex->p - h->vertex->p)) > 0){
                     h->next->vertex->p = crossptl;
                     h->prev->vertex->p = crossptr;

                     h->next->vertex->p3 = ll*rl3d + t.p3;
                     h->prev->vertex->p3 = lr*rr3d + t.p3;
                 }else{
                     h->next->vertex->p = crossptr;
                     h->prev->vertex->p = crossptl;

                     h->next->vertex->p3 = lr*rr3d + t.p3;
                     h->prev->vertex->p3 = ll*rl3d + t.p3;
                 }
             }
             break;
         }
    }

    std::vector<Face*> Faces_new;
    for(auto&f: Faces){
        HalfEdge *h = f->halfedge;
        std::vector<HalfEdge*> V;
        do{
            for(auto&t: T_crs)if(t.p == h->vertex->p){V.push_back(h);break;}
            h = h->next;
        }while(h != f->halfedge);
        if(V.size() == 2){
            HalfEdge *h1 = new HalfEdge(V[0]->vertex, EdgeType::fl);
            HalfEdge *h2 = new HalfEdge(V[1]->vertex, EdgeType::fl);
            h1->pair = h2; h2->pair = h1; h2->face = V[0]->face; h1->face = V[1]->face;
            h1->prev = V[0]->prev; h2->prev = V[1]->prev;
            h1->next = V[1]; h2->next = V[0];
            V[1]->prev->next = h2; V[0]->prev->next = h1;
            Edges.push_back(h1); Edges.push_back(h2);
            setface(h1, h1->face);
            Face *fn = new Face(h2);
            setface(h2, fn);
            Faces_new.push_back(fn);
            FoldingCurve.push_back(h1); FoldingCurve.push_back(h2);
        }
    }
    Faces.insert(Faces.end(), Faces_new.begin(), Faces_new.end()); // 連結

    for(auto&h : edges_he)delete h;
    edges_he.clear();
    edges_ol.clear();
    delete face_ol;
    return true;
}

bool FoldLine::BezierCrvOn3dSrf(std::vector<glm::f64vec3>& CtrlPts, double t, int dim, std::vector<Face*>& Faces, glm::f64vec3& v_3d){
    glm::f64vec3 v_2d{0,0,0};
    for (int i = 0; i < int(CtrlPts.size()); i++) v_2d += MathTool::BernsteinBasisFunc(dim, i, t) * CtrlPts[i];

    for(auto&f: Faces){
        if(f->IsPointInFace(v_2d)){
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

void FoldLine::ProjectBezierOn3d(int dim){
    using namespace Eigen;
    int t_size = T_crs.size();
    //if(t_size < n)return;
    dim = T_crs.size() - 1;
    MatrixXd A(t_size, t_size);
    for(int i = 0; i < t_size; i++){
        for(int j = 0; j < t_size; j++) A(i,j) = MathTool::BernsteinBasisFunc(dim, j, T_crs[i].s); 
    }

    MatrixXd b(t_size, 3);
    for(int i = 0; i < t_size; i++){ b(i, 0) = T_crs[i].p3.x; b(i,1) = T_crs[i].p3.y;  b(i,2) = T_crs[i].p3.z;}

    MatrixXd X = A.fullPivLu().solve(b);
    CtrlPts_res.resize(t_size);
    for(int i = 0; i < t_size; i++){ CtrlPts_res[i].x = X(i,0); CtrlPts_res[i].y = X(i,1); CtrlPts_res[i].z = X(i,2);}

    for(int i = 0; i < t_size; i++){
        Vector3d v{0,0,0};
        for(int j = 0; j < t_size; j++)v += A(i,j)*X.row(j);
        std::cout << i << " : " << (b.row(i) - v.transpose()).norm() << std::endl;
    }

    int csize = 500;
    double t = T_crs[0].s;
    glm::f64vec3 v{0,0,0};
    Curve_res.assign(csize, v);
    for(int i = 0; i < csize; i++){
        for(int j = 0; j < t_size; j++)Curve_res[i] += MathTool::BernsteinBasisFunc(dim, j, t) * CtrlPts_res[j];
        t += ((T_crs.end() - 1)->s - T_crs[0].s)/(double)(csize - 1);
    }
}

void FoldLine::setCoeff(std::vector<double>& a, double t, std::vector<double>& Knot, int j){
    a.assign(8, 0);
    double b[8];
    b[0] = (Knot[j+3] - Knot[j]) * (Knot[j+2] - Knot[j]) * (Knot[j+1] - Knot[j]);
    b[1] = (Knot[j+3] - Knot[j]) * (Knot[j+2] - Knot[j+1]) * (Knot[j+2] - Knot[j]);
    b[2] = (Knot[j+3] - Knot[j]) * (Knot[j+3] - Knot[j+1]) * (Knot[j+2] - Knot[j+1]);
    b[3] = (Knot[j+3] - Knot[j]) * (Knot[j+3] - Knot[j+1]) * (Knot[j+3] - Knot[j+2]);
    b[4] = (Knot[j+4] - Knot[j+1]) * (Knot[j+3] - Knot[j+1]) * (Knot[j+2] - Knot[j+1]);
    b[5] = (Knot[j+4] - Knot[j+1]) * (Knot[j+3] - Knot[j+1]) * (Knot[j+3] - Knot[j+2]);
    b[6] = (Knot[j+4] - Knot[j+1]) * (Knot[j+4] - Knot[j+2]) * (Knot[j+3] - Knot[j+2]);
    b[7] = (Knot[j+4] - Knot[j+1]) * (Knot[j+4] - Knot[j+2]) * (Knot[j+4] - Knot[j+3]);
    a[0] = (Knot[j] <= t && t < Knot[j+1] && b[0] != 0) ? 1.0/b[0] : 0;
    a[1] = (Knot[j+1] <= t && t < Knot[j+2] && b[1] != 0) ? -1.0/b[1] : 0;
    a[2] = (Knot[j+1] <= t && t < Knot[j+2] && b[2] != 0) ? -1.0/b[2] : 0;
    a[3] = (Knot[j+2] <= t && t < Knot[j+3] && b[3] != 0) ?  1.0/b[3] : 0;
    a[4] = (Knot[j+1] <= t && t < Knot[j+2] && b[4] != 0) ? -1.0/b[4] : 0;
    a[5] = (Knot[j+2] <= t && t < Knot[j+3] && b[5] != 0) ?  1.0/b[5] : 0;
    a[6] = (Knot[j+2] <= t && t < Knot[j+3] && b[6] != 0) ?  1.0/b[6] : 0;
    a[7] = (Knot[j+3] <= t && t < Knot[j+4] && b[7] != 0) ?  -1.0/b[7] : 0;
    return;
}

void FoldLine::diff(double t, std::vector<double>& Knot, std::vector<glm::f64vec3>& dP, std::vector<glm::f64vec3>& CtrlPts, int index, const int n_times){
    auto dr1 = [](double t, double a, double b, double c){return (t-a)*(t-b) + (t-b)*(t-c) + (t-a)*(t-c);};
    auto dr2 = [](double t, double a, double b, double c){return 2 * ((t-a) + (t-b) + (t-c));};
    auto dr3 = [](){return 6.0;};
    if(index == 0) t += eps;
    else if(index == (int)T_crs.size() - 1) t -=  3.0*eps;
    dP.assign(n_times, glm::f64vec3{0,0,0});
    if(n_times == 3){
        std::vector<double>a;
        for(int j = 0; j < (int)CtrlPts.size(); j++){

            setCoeff(a, t, Knot, j);
            dP[0] += (a[0]*dr1(t,Knot[j],Knot[j],Knot[j]) + a[1]*dr1(t,Knot[j],Knot[j],Knot[j+2]) + a[2]*dr1(t,Knot[j],Knot[j+1],Knot[j+3]) + a[3]*dr1(t,Knot[j],Knot[j+3],Knot[j+3]) +
                    a[4]*dr1(t,Knot[j+4],Knot[j+1],Knot[j+1]) + a[5]*dr1(t,Knot[j+1],Knot[j+3],Knot[j+4]) + a[6]*dr1(t,Knot[j+4],Knot[j+4],Knot[j+2]) + a[7]*dr1(t,Knot[j+4],Knot[j+4],Knot[j+4])) * CtrlPts[j];
            dP[1] += (a[0]*dr2(t,Knot[j],Knot[j],Knot[j]) + a[1]*dr2(t,Knot[j],Knot[j],Knot[j+2]) + a[2]*dr2(t,Knot[j],Knot[j+1],Knot[j+3]) + a[3]*dr2(t,Knot[j],Knot[j+3],Knot[j+3]) +
                    a[4]*dr2(t,Knot[j+4],Knot[j+1],Knot[j+1]) + a[5]*dr2(t,Knot[j+1],Knot[j+3],Knot[j+4]) + a[6]*dr2(t,Knot[j+4],Knot[j+4],Knot[j+2]) + a[7]*dr2(t,Knot[j+4],Knot[j+4],Knot[j+4])) * CtrlPts[j];
            double coef = 0.0;
            for(auto&_a: a)coef += _a*dr3();
            dP[2] += coef * CtrlPts[j];
        }
    }else{
        glm::f64vec3 v = bspline(CtrlPts, t, n_times, Knot);
        glm::f64vec3 vm = bspline(CtrlPts, t - eps, n_times, Knot);
        glm::f64vec3 vp = bspline(CtrlPts, t + eps, n_times, Knot);
        glm::f64vec3 vm2 = bspline(CtrlPts, t - 2.0 * eps, n_times, Knot);
        glm::f64vec3 vp2 = bspline(CtrlPts, t + 2.0 * eps, n_times, Knot);
        dP[0] = (vp - vm)/(2.0 * eps);
        dP[1] = (vp - 2.0 * v + vm)/pow(eps, 2);
        dP[2] = (-vm2 + 2.0 * vm - 2.0 * vp + vp2)/(2.0 * pow(eps, 3));
    }
    return;
}

void FoldLine::LSM_apply_tau(int dim){
    //type: 0->x, 1->y, 2->x*x, 3->xy
    auto sum = [](std::vector<double>&y, int type){
        double z = 0.0;
        for(int i = 0; i < (int)y.size(); i++)z += (type == 0) ? i: (type==1)?y[i]: (type==2)?i*i: i*y[i];
        return z;
    };
    std::string file = "result_lsm.csv";
    std::ofstream ofs(file);
    int n = T_crs.size();
    std::vector<double>_y(n);
    for(int i = 0; i < n; i++)_y[i] = T_crs[i].tau;

    ofs << "row data(y), result(y)" << std::endl;
    if(dim == 1){
        double xy = sum(_y,3), xx = sum(_y,2), x = sum(_y,0), y = sum(_y,1);
        double a = ((double)n*xy - x*y)/((double)n * xx - x*x);
        double b = (xx*y - xy*x)/((double)n*xx - x*x);
        for(int i = 0; i < n; i++){
            ofs << _y[i] << ", " << a*(double)i + b << std::endl;
            T_crs[i].tau = a*(double)i + b;
        }
    }

}


double optimizer::Fvert(std::vector<double>& X){
    double F = 0.0;
    std::vector<glm::f64vec3> m_new;
    int psize = 8;
    double px, py;
    glm::f64vec3 c, c2;
    int i = 0;
    for(auto&f: ptchDevSrf){
        HalfEdge *h = f->halfedge;
        glm::f64vec3 o = h->vertex->p3;
        glm::f64vec3 e1 = glm::normalize(h->next->vertex->p3 - o);
        glm::f64vec3 e2 = glm::normalize(h->prev->vertex->p3 - o);
        do{
            px = X[psize * i]; py = X[psize * i + 1];
            c = glm::f64vec3{X[psize * i + 2], X[psize * i + 3], X[psize * i + 4]};
            c2 = glm::f64vec3{X[psize * i + 5], X[psize * i + 6], X[psize * i + 7]};
            m_new.push_back(h->vertex->p3 + c2 + glm::cross(c, o) + px*glm::cross(c, e1) + py*glm::cross(c, e2));
            h = h->next;
            i++;
        }while(h != f->halfedge);
    }
    int k = m_new.size();
    for(int i = 0; i < k; i++){
        for(int j = i+1; j < k; j++)F += pow(glm::length(m_new[i] - m_new[j]), 2);
    }
    return F;
}

double optimizer::Ffit(){
    double fitval = 0.0;
    for(auto&f: Faces){
        HalfEdge *h = f->halfedge;

        do{
            h = h->next;
        }while(h != f->halfedge);
    }
    return fitval;
}

double optimizer::Ffair(){
    double f = 0.0;
    for(int i = 0; i < (int)FoldingCurve.size() - 1; i++)f += pow(glm::length(FoldingCurve[i-1]->p3 - 2.0 * FoldingCurve[i]->p3 + FoldingCurve[i+1]->p3), 2);
    return f;
}

double optimizer::Fconv(){
    auto sigTriArea = [](glm::f64vec3 a, glm::f64vec3 b, glm::f64vec3 c){};
    double f = 0.0;

}

void optimizer::alignedNewRulingDirection(){
    int psize = FoldingCurve.size();

    for(int i = 1; i < psize - 1; i++){
        HalfEdge *h = getEdge(FoldingCurve[i-1], FoldingCurve[i], ptchDevSrf), *h2 = getEdge(FoldingCurve[i], FoldingCurve[i+1], ptchDevSrf);
        glm::f64vec3 e = glm::normalize(h->vertex->p3 - h->next->vertex->p3), e2 = glm::normalize(h2->next->vertex->p3 - h2->vertex->p3);
        double a = acos(glm::dot(-e, glm::normalize(h->next->next->vertex->p3 - h->next->vertex->p3))), b = acos(glm::dot(e2, glm::normalize(h2->prev->vertex->p3 - h2->vertex->p3)));
        double c = cos(a + b);
        glm::f64vec3 x = h->pair->prev->vertex->p3;
    }
}

void optimizer::apply(double wfit, double wfair){

    double Px, Py, cx, cy, cz, c2x, c2y, c2z; //(px, py), c(cx, cy, cz), c-(c2x, c2y, c2z)
    std::vector<double> X;
    for(int i = 0; i < (int)Faces.size(); i++){
        glm::f64vec3 o = Faces[i]->halfedge->vertex->p3;

    }
    double F = Fvert(X) + wfit * Ffit() + wfair * Ffair();

}

HalfEdge* optimizer::getEdge(Vertex* start, Vertex* end, const std::vector<Face*>& searchedFace){
    for(const auto& f: searchedFace){
        HalfEdge *h = f->halfedge;
        do{
            if(h->vertex == start && h->next->vertex == end)return h;
            h = h->next;
        }while(h != f->halfedge);
    }
    return nullptr;
}
