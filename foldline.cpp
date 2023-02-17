
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

bool FoldLine::modify2DRulings(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, const std::vector<HalfEdge*>& edge_outline, int dim, int t_type){
    using namespace MathTool;
    /*
    auto VectorDigAlign = [](glm::f64vec3& v, int dig = 13){
        glm::f64vec3 vi = glm::f64vec3{(int)v.x, (int)v.y, (int)v.z};
        glm::f64vec3 vf = v - vi;
        int n = pow(10, dig); int x = vf.x * n, y = vf.y * n, z = vf.z * n;
        v = glm::f64vec3{vi.x + (double)x/(double)n, vi.y + (double)y/(double)n, vi.z + (double)z/(double)n}; };
    T_crs.clear();

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
        int divSize = 1;

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
                    )continue;
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
    }else if(type == PaintTool::FoldLine_test){
        int faceVerticesSize = edges_ol.size();
        for(int i = 0; i < faceVerticesSize; i++){
            bool hasPoint = MathTool::IsIntersect(edges_ol[i], edges_ol[(i+1)%faceVerticesSize], CtrlPts[0], CtrlPts[1]);
            if(hasPoint){
                glm::f64vec3 NewPoint = MathTool::getIntersectionPoint(edges_ol[i], edges_ol[(i+1)%faceVerticesSize], CtrlPts[0], CtrlPts[1]);
                CrvPt_FL P(NewPoint, i);
                T_crs.push_back((P));
            }
        }


    }

    std::sort(T_crs.begin(), T_crs.end());
#if 0
    std::vector<double> Knot, Knot2d;
    std::vector<glm::f64vec3> dP(3);


    double CurveLen3d = 0.0,  CurveLen2d = 0.0;
    Curve_res = GlobalSplineInterpolation(T_crs, CtrlPts_res, Knot, CurveLen3d, true, dim, t_type);
    GlobalSplineInterpolation(T_crs, CtrlPts_res2d, Knot2d, CurveLen2d, false, dim, t_type);
    std::cout << "error of length " << CurveLen3d - CurveLen2d << std::endl;

    for(auto&t: T_crs){
        int index = &t - &T_crs[0];
         diff(t.s, Knot, dP, CtrlPts_res, index, dim);

         t.T3d = glm::normalize(dP[0]);
         t.k3d = glm::length(glm::cross(dP[0], dP[1]))/std::pow(glm::length(dP[0]), 3);
         t.tau = glm::dot(dP[0], glm::cross(dP[1], dP[2]))/std::pow(glm::length(glm::cross(dP[0], dP[1])), 2);

         //diff(t.s - eps, Knot, dP, CtrlPts_res, index, dim); glm::f64vec3 Tm = glm::normalize(dP[0]);
         //diff(t.s + eps, Knot, dP, CtrlPts_res, index, dim); glm::f64vec3 Tp = glm::normalize(dP[0]);
         //t.N3d = glm::normalize((Tp - Tm)/(2.0*eps));
         //t.B3d = glm::normalize(glm::cross(t.T3d, t.N3d));
         t.N3d = glm::normalize(glm::cross(dP[0], glm::cross(dP[1], dP[0])));
         t.B3d = glm::normalize(glm::cross(dP[0], dP[1]));


         diff(t.s, Knot2d, dP, CtrlPts_res2d, index, dim);
         glm::f64vec3 T = glm::normalize(dP[0]), N = glm::f64vec3{-T.y, T.x, 0}, B = glm::cross(T, N);
         t.k2d = glm::length(glm::cross(dP[0], dP[1]))/std::pow(glm::length(dP[0]), 3);
         t.T2d = T;
         t.N2d = N;
         t.B2d = B;
         diff(t.s - eps, Knot2d, dP, CtrlPts_res2d, index, dim); t.k2d_m = glm::length(glm::cross(dP[0], dP[1]))/std::pow(glm::length(dP[0]), 3);
         diff(t.s + eps, Knot2d, dP, CtrlPts_res2d, index, dim); t.k2d_p = glm::length(glm::cross(dP[0], dP[1]))/std::pow(glm::length(dP[0]), 3);

        //弧長パラメータ上で中心差分法をおこなってみる -　今は０から１の範囲のｔ
         curvepointfile << t.s << ", " << t.k2d <<" , " << t.k3d << ", " << t.tau <<", " << t.p.x <<", " <<t.p.y <<", " <<t.p.z  <<", "<< t.p3.x << ", " <<t.p3.y <<", "<< t.p3.z <<", " << t.T3d.x << ", " << t.T3d.y <<", " << t.T3d.z <<", " << t.N3d.x << ", " << t.N3d.y<<", " << t.N3d.z << ", " << t.B3d.x << ", " << t.B3d.y << ", " << t.B3d.z << std::endl;
             for(auto&edge: Edges){
                 if(is_point_on_line(t.p,edge->vertex->p, edge->next->vertex->p)){
                     glm::f64vec3 f_front = edge->face->getNormalVec();
                     if(edge->edgetype == EdgeType::r){
                         glm::f64vec3 f_front2 = edge->pair->face->getNormalVec();
                         f_front = glm::normalize((f_front + f_front2)/2.0);
                     }//std::swap(t.B3d, t.N3d);
                     //if(glm::dot(t.B3d, f_front) < 0)t.B3d *= -1;
                     if(abs(glm::dot(t.B3d, f_front)) < abs(glm::dot(t.N3d, f_front))){
                         //std::cout<<"hello"<<std::endl;
                         //glm::f64vec3 tmp = t.B3d; t.B3d = t.N3d; t.N3d = tmp;
                     }
                 }
             }
             //if(glm::dot(glm::normalize(glm::cross(t.B3d, t.T3d)), t.N3d) < 0) t.N3d *= -1;

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
         double phi_bl = rad_2d(t.k3d, t.tau, t.a, -t.da), phi_br = rad_2d(t.k3d, t.tau, t.a, t.da);
         resultfile <<std::setprecision(15) << t.k3d << ", "<<t.k2d << " , " << t.tau << " , " << t.a << " , " << t.da << ", " << phi_bl * 180/M_PI << " , " << phi_br * 180/M_PI  <<", " << glm::length(t.p - t.p3)<< ", " << t.k2d_m << ", " <<t.k3d_m << ", " << t.k2d_p << ", "<< t.k3d_p << std::endl;
    }

    //return true;
    for(auto&t: T_crs){
        double phi_bl = rad_2d(t.k3d, t.tau, t.a, -t.da), phi_br = rad_2d(t.k3d, t.tau, t.a, t.da);
         glm::f64vec3 rr3d = glm::normalize(cos(phi_br) * t.T3d + sin(phi_br) * cos(t.a) * t.N3d + sin(phi_br) * sin(t.a) * t.B3d);
         glm::f64vec3 rl3d = glm::normalize(cos(phi_bl) * t.T3d - sin(phi_bl) * cos(t.a) * t.N3d + sin(phi_bl) * sin(t.a) * t.B3d);
         glm::f64vec3 rl2d = (glm::rotate(phi_bl, t.B2d) * glm::f64vec4{t.T2d,1});
         glm::f64vec3 rr2d = (glm::rotate(-phi_br, t.B2d) * glm::f64vec4{t.T2d,1});       

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
#endif
    for(auto&t: T_crs){
        std::vector<HalfEdge*> H_new;
        for(auto&he: Edges){
            H_new = he->Split(&t,Edges);
            if(H_new.size() > 0)break;
        }
    }
    std::vector<Face*> Faces_new;
    FoldingCurve.clear();
    double tbef = -1;
    int faceNum = Faces.size();

    for(int i = 0; i < faceNum; i++){
        Face *f = Faces[i];
        HalfEdge *h = f->halfedge;
        tbef = -1;
        std::vector<HalfEdge*> V;
        do{
            for(auto&t: T_crs)if(t.p == h->vertex->p){
                if(tbef == -1)V.push_back(h);
                else if(tbef != -1 && tbef < t.s)V.push_back(h);
                else if(tbef != -1 && tbef > t.s)V.insert(V.begin(), h);
                tbef = t.s;
                break;
            }
            h = h->next;
        }while(h != f->halfedge);

        if(V.size() == 2){

            HalfEdge *h1 = new HalfEdge(V[0]->vertex, EdgeType::fl);
            HalfEdge *h2 = new HalfEdge(V[1]->vertex, EdgeType::fl);
            if(type == PaintTool::FoldLine_test ||type == PaintTool::FoldLine_line){
                ruling *r = new ruling(V[0]->vertex, V[1]->vertex);
                h1->r = r; h2->r = r;
            }
            h1->pair = h2; h2->pair = h1;
            V[1]->prev->next = h2; V[0]->prev->next = h1;
            h1->prev = V[0]->prev; h2->prev = V[1]->prev;
            h1->next = V[1]; h2->next = V[0];

            Edges.push_back(h1); Edges.push_back(h2);

            Face *fn = new Face(h2);
            f->ReConnect(h1);
            fn->ReConnect(h2);
            Faces_new.push_back(fn);
            FoldingCurve.push_back(h1);
            if(f == Faces.back())FoldingCurve.push_back(h2);
        }

    }

    Faces.insert(Faces.end(), Faces_new.begin(), Faces_new.end()); // 連結
    return true;
    */
}

bool FoldLine::SplitFace4DebugAAAMethod(glm::f64vec3& NewPoint, std::vector<Face*> &faces, std::vector<HalfEdge*>& edges, std::vector<Vertex*>& vertices){
    CtrlPts.push_back(NewPoint);
    if(type != PaintTool::FoldLine_test || faces.empty())return false;

    glm::f64vec3 CrossPoint;
    int faceSize = faces.size();
    glm::f64vec3 p = NewPoint, q = glm::f64vec3{p.x, p.y + 1000, 0}; p.y -= 200;
    for(int i = 0; i < faceSize; i++){
        Face *f = faces[i];
        std::vector<HalfEdge*> _NewEdges;
        std::vector<glm::f64vec3> CrossPoints;

        for(auto&e: edges){
            bool _hasCrossPoint = e->hasCrossPoint2d(p,q, CrossPoint, true);
            if(_hasCrossPoint){
                if(std::find(CrossPoints.begin(), CrossPoints.end(), CrossPoint) == CrossPoints.end())CrossPoints.push_back(CrossPoint);
            }
       }

        while(!CrossPoints.empty()){
            CrossPoint = CrossPoints.back();
            CrossPoints.pop_back();

            HalfEdge *EndPoint = nullptr;
            for(auto&_e: edges){
                if(glm::distance(_e->vertex->p, CrossPoint) < 1e-5)EndPoint = _e;
            }
            if(EndPoint != nullptr){
                _NewEdges.push_back(EndPoint);
                continue;
            }
            for(auto&e: edges){
                auto res = e->Split(new Vertex(CrossPoint), edges);
                if(!res.empty()){
                    double t = glm::length(res[0]->prev->vertex->p - res[0]->vertex->p)/glm::length(res[0]->prev->vertex->p - res[0]->next->vertex->p);
                    res[0]->vertex->p3 = res[0]->prev->vertex->p3 + t * (res[0]->next->vertex->p3 - res[0]->prev->vertex->p3);
                    _NewEdges.push_back(res[0]);
                    break;
                }
            }
        }

        if(_NewEdges.size() == 2){
            HalfEdge *h1 = new HalfEdge(_NewEdges[0]->vertex, EdgeType::r);
            HalfEdge *h2 = new HalfEdge(_NewEdges[1]->vertex, EdgeType::r);
            ruling *r = new ruling(_NewEdges[0]->vertex, _NewEdges[1]->vertex);
            h1->r = r; h2->r = r;
            h1->pair = h2; h2->pair = h1;
            _NewEdges[1]->prev->next = h2; _NewEdges[0]->prev->next = h1;
            h1->prev = _NewEdges[0]->prev; h2->prev = _NewEdges[1]->prev;
            h1->next = _NewEdges[1]; h2->next = _NewEdges[0];

            Face *fn = new Face(h2);
            f->ReConnect(h1);
            fn->ReConnect(h2);
            faces.push_back(fn); edges.push_back(h1); edges.push_back(h2);
            vertices.push_back(_NewEdges[0]->vertex);vertices.push_back(_NewEdges[1]->vertex);
            return true;
        }
    }

    return false;
}

void FoldLine::applyAAAMethod(const std::vector<glm::f64vec3>& edges_ol, double a){

    double phi02 = AngleIn2Edges(FoldingCurve[0], FoldingCurve[0]->pair->next);
    double phim1 = AngleIn2Edges(FoldingCurve[0], FoldingCurve[0]->pair->next);
    _FoldingAAAMethod(a, phi02, phim1, edges_ol);
}


void FoldLine::_FoldingAAAMethod(double & a, double phi02, double phim1, const std::vector<glm::f64vec3>& edges_ol){
    auto IsSame = [](double c, double b){
        return (abs(c - b) <= 1e-5)? true: false;
    };
    double a2 = 0, l;
    glm::f64vec3 e, e2, x;
    glm::f64vec3 N, r, n;
    glm::f64vec3 e_bef;
    std::string file = "axis.csv"; std::ofstream AXIS(file);
    std::string file2 = "curvepoint.csv"; std::ofstream curvepointfile(file2);
    curvepointfile << a << ", " << FoldingCurve[0]->vertex->p3.x <<", " << FoldingCurve[0]->vertex->p3.y << ", " << FoldingCurve[0]->vertex->p3.z << std::endl;
    {
        e = glm::normalize(FoldingCurve[0]->vertex->p3 - FoldingCurve[1]->vertex->p3);
        N = glm::normalize(glm::cross(glm::f64vec3{glm::rotate (a, e)  * glm::f64vec4{FoldingCurve[0]->face->getNormalVec(),1}}, e));//phi02の回転軸
        r = glm::rotate(phi02, N) * glm::f64vec4{-e,1};r = glm::normalize(r);//新しいruling方向
        l = glm::length(FoldingCurve[0]->pair->next->next->vertex->p - FoldingCurve[0]->vertex->p);
        FoldingCurve[0]->pair->next->next->vertex->p3 = l * r + FoldingCurve[0]->vertex->p3;
    }

    int flsize = FoldingCurve.size();
    for(int i = 1; i < flsize - 1; i++){
        x = FoldingCurve[i]->vertex->p3;
        curvepointfile << a << ", " << x.x <<", " << x.y << ", " << x.z << std::endl;
        AXIS << x.x << ", " << x.y << ", " << x.z << ", " << FoldingCurve[i-1]->next->next->vertex->p3.x << ", " << FoldingCurve[i-1]->next->next->vertex->p3.y << ", " << FoldingCurve[i-1]->next->next->vertex->p3.z << std::endl;
        e = glm::normalize(FoldingCurve[i-1]->vertex->p3 - x);
        e2 = glm::normalize(FoldingCurve[i+1]->vertex->p3 - x);
        double beta = glm::angle(e,e2);
        //beta = M_PI - abs(FoldingCurve[i]->prev->r->Gradation/255.0)*M_PI;
        glm::f64vec3 ax2 = IsSame(beta, std::numbers::pi) ? FoldingCurve[i]->face->getNormalVec():ProjectionVector(e2,e);ax2 = glm::normalize(ax2);//折角度aの回転軸
        if(i != 1){
            glm::f64vec3 e_next = IsSame(beta, std::numbers::pi) ? FoldingCurve[i]->face->getNormalVec(): (ProjectionVector(e,e2)); e_next = glm::normalize(e_next);
            double tau = (glm::dot(e_bef, e_next) > 1) ? 0: (glm::dot(e_bef, e_next) < -1) ? std::numbers::pi: glm::angle(e_bef,e_next);
            a = (a2 - tau < 0) ? a2 - tau + 2.0 * std::numbers::pi: a2 - tau;
            std::cout << "---------------------------------" << std::endl;
            if(!std::isfinite(a))std::cout << "a is not finite " << i << std::endl;
            if(!std::isfinite(a2))std::cout << "a' is not finite " << i << std::endl;
            std::cout << i << " : tau = " << tau << " , a = " <<  glm::degrees(a) << " , a' = " << glm::degrees(a2) << std::endl;

        }

        double phi3 = AngleIn2Edges(FoldingCurve[i-1]->next, FoldingCurve[i]);
        double phi4 = AngleIn2Edges(FoldingCurve[i-1]->next, FoldingCurve[i-1]->pair);
        double k = 2.0 * std::numbers::pi - phi3 - phi4;
        if(!(beta <= phi3 + phi4 && phi3 + phi4 <= 2.0 * std::numbers::pi - beta)){std::cout <<"exceed boundary condition   " << beta << " , " << phi3+phi4<< std::endl; exit(0);}
        double phi1 = IsSame(glm::sin(beta)*glm::cos(a), glm::sin(k))? std::numbers::pi/2.0: atan((glm::cos(k) - glm::cos(beta))/(glm::sin(beta)*glm::cos(a)- glm::sin(k)));
        //phi1 *= (cos(k) - cos(beta) < 0)? -1: 1; phi1 = (phi1 < 0) ? phi1 + std::numbers::pi: phi1;//0 ~ piの範囲に収める
        //phi1 = glm::atan((cos(k) - cos(beta))/(sin(beta)*cos(a)- sin(k)));
        if(phi1 < 0){ phi1 += std::numbers::pi;std::cout <<"phi1 < 0" << std::endl;}
        double phi2 = k - phi1;
        std::cout << std::setprecision(10) <<  glm::degrees(phi1) << ", " << glm::degrees(phi2) << ", " << glm::degrees(phi3) << ", " << glm::degrees(phi4)<< " : " <<  glm::degrees(k) <<  std::endl;
        std::cout << (glm::cos(k) - glm::cos(beta)) << " , " << (glm::sin(beta)*glm::cos(a)- glm::sin(k)) << std::endl;
        if(phi1 < 0 || phi2 < 0){ std::cout <<"minus phi1 or phi2" << std::endl; exit(0);}
        N = glm::normalize(glm::cross(e,glm::f64vec3{glm::rotate (a, e)  * glm::f64vec4{ax2,1}}));//phi1の回転軸
        r = glm::rotate(phi1, N) * glm::f64vec4{e,1};r = glm::normalize(r);//新しいruling方向
        n = glm::rotate(phi1, glm::f64vec3{0,0,1})*glm::f64vec4{glm::normalize(FoldingCurve[i-1]->vertex->p- FoldingCurve[i]->vertex->p),1};n = glm::normalize(n);//展開図のruling方向
        //std::cout << "rotation  " << glm::to_string(glm::rotate (a, e)  * glm::f64vec4{ax2,1}) << std::endl;
        std::cout << i << "   ruling direction : r = " <<  glm::to_string(r) << " , n = " <<glm::to_string(n) << ", N = " << glm::to_string(N) << " , e = " << glm::to_string(e) << std::endl;
        //std::cout << glm::to_string(glm::normalize(FoldingCurve[i-1]->pair->prev->vertex->p - FoldingCurve[i]->vertex->p)) << " , " << glm::to_string(glm::normalize(FoldingCurve[i-1]->vertex->p - FoldingCurve[i]->vertex->p)) << std::endl;
        std::cout <<"r, e  " << glm::angle(e,r) - phi1 << " , " << glm::angle(r, e2) - phi2 << std::endl;
        glm::f64vec3 crossPoint;
        bool hasPointOnEdge = setPoint(edges_ol, n, FoldingCurve[i]->vertex->p, crossPoint);
        if(hasPointOnEdge){
            FoldingCurve[i-1]->pair->prev->vertex->p = crossPoint;
            l = glm::distance(crossPoint, FoldingCurve[i]->vertex->p);
            FoldingCurve[i-1]->pair->prev->vertex->p3 = l * r + FoldingCurve[i]->vertex->p3;
        }else std::cout << "cross point could not be founded in " << i << std::endl;
        double sin_a = (sin(phi1)*sin(a)/sin(phi2) > 1) ? 1: (sin(phi1)*sin(a)/sin(phi2) < -1)? -1: (sin(phi1)*sin(a)/sin(phi2));
        //double cos_a = ((cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta)) > 1) ? 1: ((cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta)) < -1)? -1:  (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));

        glm::f64vec3 r2 = ProjectionVector(r,e2); r2 = glm::normalize(r2);      
        glm::f64vec3 ax_beta = IsSame(beta, std::numbers::pi) ? FoldingCurve[i]->face->getNormalVec(): e;
        ax_beta = ProjectionVector(ax_beta, e2); ax_beta = glm::normalize(ax_beta);
        double cos_a = (glm::dot(r2, ax_beta) > 1) ? 1: (glm::dot(r2, ax_beta) < -1)? -1: glm::dot(r2, ax_beta);//値が期待しているのと違う
        cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
        //std::cout << "ruling " << glm::to_string(r) << " , " << glm::to_string(r2) << " , " << glm::to_string(ax_beta) << " , " << glm::to_string(e2) << std::endl;
        //double tan_a = atan2(sin(phi1)*sin(a)*sin(beta),(cos(phi1) - cos(phi2)*cos(beta))*cos(a));

        //tan_a = (tan_a < 0) ? tan_a + 2.0*std::numbers::pi: tan_a;
        //a2 = (sin(phi1)*sin(a)/sin(phi2) > 1) ? asin(1.0): (sin(phi1)*sin(a)/sin(phi2) < -1)? asin(-1.0): asin(sin(phi1)*sin(a)/sin(phi2));//a'
        //a2 = (sin_a < 0 && cos_a < 0)? 2.0*std::numbers::pi - acos(cos_a): (sin_a < 0 && cos_a >= 0)? 2.0*std::numbers::pi + asin(sin_a): acos(cos_a);
        a2 = (sin_a >= 0 && cos_a >= 0)? asin(sin_a): (sin_a >= 0 && cos_a < 0)?std::numbers::pi - asin(sin_a): (sin_a < 0 && cos_a < 0)? std::numbers::pi + asin(sin_a): 2.0*std::numbers::pi + asin(sin_a);
        std::cout <<"sin " << sin_a << " , cos " << cos_a << " , a' =  " << glm::degrees(a2) << std::endl;
        if(i == (int)FoldingCurve.size() -2)break;
        e_bef = IsSame(beta, std::numbers::pi) ? FoldingCurve[i]->face->getNormalVec(): (ProjectionVector(e,e2)); e_bef = glm::normalize(e_bef);
    }

    {
        e = glm::normalize(FoldingCurve[flsize-2]->vertex->p3 - FoldingCurve[flsize-1]->vertex->p3);
        N = glm::normalize(glm::cross(e,glm::f64vec3{glm::rotate(a, e)  * glm::f64vec4{FoldingCurve[flsize-2]->face->getNormalVec(),1}}));//phi02の回転軸
        r = glm::rotate(phim1, N) * glm::f64vec4{e,1};r = glm::normalize(r);//新しいruling方向
        l = glm::length(FoldingCurve[flsize-2]->pair->prev->vertex->p - FoldingCurve[flsize-1]->vertex->p);
        FoldingCurve[flsize-2]->pair->prev->vertex->p3 = l * r + FoldingCurve[flsize-1]->vertex->p3;
    }

        std::cout <<"||||||||||||||||||||||||||||||" << std::endl;
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
    //else if(index == (int)T_crs.size() - 1) t -=  3.0*eps;
    dP.assign(n_times, glm::f64vec3{0,0,0});
    if(false){
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


double FoldLine::AngleIn2Edges(HalfEdge *p, HalfEdge *p2, bool Is3d){
  if(Is3d)return glm::angle(glm::normalize(p->next->vertex->p3 - p->vertex->p3), glm::normalize(p2->next->vertex->p3 - p2->vertex->p3));
  return glm::angle(glm::normalize(p->next->vertex->p - p->vertex->p), glm::normalize(p2->next->vertex->p - p2->vertex->p));
}

