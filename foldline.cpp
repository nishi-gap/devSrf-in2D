
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

bool FoldLine::moveCtrlPt(glm::f64vec3& p, int movePtIndex){
    double d = 10.0;
    if(CtrlPts.size() < movePtIndex || movePtIndex < 0)return false;
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

bool FoldLine::ChangeColor(OUTLINE *outline, int val, int dim){
    //color += val;
    //color = (color > 255) ? 255. : (color < -255)? -255. : color;
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

bool FoldLine::setPoint(const std::vector<Vertex*>& Poly_v, glm::f64vec3 N, glm::f64vec3& cp, glm::f64vec3& crossPoint){
    double l = 1000;
    bool IsIntersected = false;
    double minDist = 1000;
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    N = glm::normalize(N);
    int n = Poly_v.size();
    for(int i = 0; i < n; i++){
        glm::f64vec3 v = Poly_v[(i+1) % n]->p - Poly_v[i]->p;
        A(0,0) = v.x; A(0,1) = -l*N.x; A(1,0) = v.y; A(1,1) = -l * N.y;
        b(0) = cp.x - Poly_v[i]->p.x; b(1) = cp.y - Poly_v[i]->p.y;
        if(abs(glm::dot(glm::normalize(v), N))>= 1-DBL_EPSILON)continue;
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

double FoldLine::Fbend(std::vector<glm::f64vec3>&C, std::vector<glm::f64vec3>&C2){
    double f = 0.0;

    std::vector<glm::f64vec3> S;
    auto getBisector = [](glm::f64vec3 p, glm::f64vec3 q, glm::f64vec3 p2, glm::f64vec3 q2){
        glm::f64vec3 v = glm::normalize(q - p), v2 = glm::normalize(q2 - p2);
        if(abs(glm::dot(v, v2)) >= 1 - DBL_EPSILON)return v;
        return glm::normalize(v + v2);
    };
    auto calcCrossPoint = [](glm::f64vec3 p, glm::f64vec3 q, Vertex *p2, Vertex *q2){
        glm::f64vec3 v2 = (q2->p - p2->p);
        Eigen::Matrix2d A;
        Eigen::Vector2d b;
        b(0) = p2->p.x - p.x; b(1) = p2->p.y - p.y;
        A(0,0) = q.x; A(0,1) = -v2.x; A(1,0) = q.y; A(1,1) = -v2.y;
        Eigen::Vector2d ts = A.colPivHouseholderQr().solve(b);

        return ts(1);
    };

    auto modifyS = [](std::vector<double>&S){
        double s_min = *std::min_element(S.begin(), S.end()), s_max = *std::max_element(S.begin(), S.end());
        for(auto&s: S){
            if(s_min < 0 && s_max <= 1)s += abs(s_min);
            if(0 <= s_min && s_max > 1 )s -= (abs(s_max) - 1);
        }
    };
    //i = 0
    S.push_back(getBisector(FoldingCurve[0]->vertex->p, FoldingCurve[0]->pair->next->next->vertex->p, FoldingCurve[1]->vertex->p, FoldingCurve[1]->prev->vertex->p));
    for(int i = 1; i < (int)FoldingCurve.size()-1; i++){
        //R_i-1とR_iのなす角度を二等分する直線を求める
        S.push_back(getBisector(FoldingCurve[i]->vertex->p, FoldingCurve[i]->prev->vertex->p, FoldingCurve[i+1]->vertex->p, FoldingCurve[i+1]->prev->vertex->p));
    }   
    glm::f64vec3 c_prev, c2_prev;
    std::vector<double> Cs, C2s;
    double s = 0.5; Cs.push_back(s);
    c_prev = s * (FoldingCurve[0]->pair->next->next->vertex->p - FoldingCurve[0]->vertex->p) + FoldingCurve[0]->vertex->p;
    //C.push_back(s * (FoldingCurve[0]->pair->next->next->vertex->p3 - FoldingCurve[0]->vertex->p3) + FoldingCurve[0]->vertex->p3);
    s = 0.75; C2s.push_back(s);
    c2_prev = s *(FoldingCurve[0]->pair->next->next->vertex->p - FoldingCurve[0]->vertex->p) + FoldingCurve[0]->vertex->p;
    //C2.push_back(0.75 * (FoldingCurve[0]->pair->next->next->vertex->p3 - FoldingCurve[0]->vertex->p3) + FoldingCurve[0]->vertex->p3);
    for(int i = 0; i < (int)S.size(); i++){
        glm::f64vec3 S_orth = glm::f64vec3{-S[i].y, S[i].x, 0};
        s = calcCrossPoint(c_prev, S_orth, FoldingCurve[i+1]->vertex, FoldingCurve[i+1]->prev->vertex);
        Cs.push_back(s);
       //C.push_back(s * (FoldingCurve[i+1]->prev->vertex->p3 - FoldingCurve[i+1]->vertex->p3) + FoldingCurve[i+1]->vertex->p3);
       c_prev = s * (FoldingCurve[i+1]->prev->vertex->p - FoldingCurve[i+1]->vertex->p) + FoldingCurve[i+1]->vertex->p;
       s = calcCrossPoint(c2_prev, S_orth, FoldingCurve[i+1]->vertex, FoldingCurve[i+1]->prev->vertex);
       C2s.push_back(s);
       //C2.push_back(s * (FoldingCurve[i+1]->prev->vertex->p3 - FoldingCurve[i+1]->vertex->p3) + FoldingCurve[i+1]->vertex->p3);
       c2_prev = s * (FoldingCurve[i+1]->prev->vertex->p - FoldingCurve[i+1]->vertex->p) + FoldingCurve[i+1]->vertex->p;
    }
    modifyS(Cs); modifyS(C2s);
    C.push_back(Cs[0] * (FoldingCurve[0]->pair->next->next->vertex->p3 - FoldingCurve[0]->vertex->p3) + FoldingCurve[0]->vertex->p3);
    C2.push_back(C2s[0] * (FoldingCurve[0]->pair->next->next->vertex->p3 - FoldingCurve[0]->vertex->p3) + FoldingCurve[0]->vertex->p3);
    for(int i = 1; i <= (int)S.size();i++){
        C.push_back(Cs[i] * (FoldingCurve[i]->prev->vertex->p3 - FoldingCurve[i]->vertex->p3) + FoldingCurve[i]->vertex->p3);
        C2.push_back(C2s[i] * (FoldingCurve[i]->prev->vertex->p3 - FoldingCurve[i]->vertex->p3) + FoldingCurve[i]->vertex->p3);
    }
    /*
    std::cout << glm::degrees(glm::angle(glm::normalize(C[1] - C[0]), glm::normalize(FoldingCurve[0]->pair->next->next->vertex->p3 - FoldingCurve[0]->vertex->p3))) << std::endl;
    std::cout << glm::degrees(glm::angle(glm::normalize(C[0] - C[1]), glm::normalize(FoldingCurve[1]->prev->vertex->p3 - FoldingCurve[1]->vertex->p3))) << std::endl;
    std::cout << "-----------"<<std::endl;
    for(int i = 1; i < (int)S.size(); i++){
        std::cout << glm::degrees(glm::angle(glm::normalize(C[i+1] - C[i]), glm::normalize(FoldingCurve[i]->prev->vertex->p3 - FoldingCurve[i]->vertex->p3))) << std::endl;
        std::cout << glm::degrees(glm::angle(glm::normalize(C[i] - C[i+1]), glm::normalize(FoldingCurve[i+1]->prev->vertex->p3 - FoldingCurve[i+1]->vertex->p3))) << std::endl;
        std::cout <<"-----------"<<std::endl;
    }*/
    //直線に対して直交するように伸ばす直線がR_i-1,R_iのなす角度が等しくなるようば場所にm_iを配置
    //式(2),(3)を実装する
    for(int i = 1; i < (int)C.size()-1;i++){
        double h = glm::length(C[i-1] - C2[i-1]);
        double L = glm::length((C[i] - C[i-1])/2.0) + glm::length((C[i] - C[i+1])/2.0);
        double L2 = glm::length((C2[i] - C2[i-1])/2.0) + glm::length((C2[i] - C2[i+1])/2.0);
        double w = h * (log(L2) - log(L))/(log(L2) - log(L));
        f += w * std::pow(glm::length(FoldingCurve[i]->face->getNormalVec() - FoldingCurve[i-1]->face->getNormalVec()), 2);
    }
    //Cではみ出している場所を導出
    //すべての点が入るように修正
    //C2も同様にする
    //CとC2が交差しないようにする(tsの値を保持しておけば比較は簡単にできるC2のts > Cのts)

    return f;
}

bool FoldLine::Optimization(std::vector<glm::f64vec3>&C, std::vector<glm::f64vec3>&C2){
    std::cout << Fbend(C, C2) << std::endl;
    return false;
}

bool FoldLine::modify2DRulings(std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges, std::vector<Vertex*>& Vertices, const std::vector<Vertex*>& Poly_v, int dim, int t_type){
    using namespace MathTool;

    auto VectorDigAlign = [](glm::f64vec3& v, int dig = 13){
        glm::f64vec3 vi = glm::f64vec3{(int)v.x, (int)v.y, (int)v.z};
        glm::f64vec3 vf = v - vi;
        int n = pow(10, dig); int x = vf.x * n, y = vf.y * n, z = vf.z * n;
        v = glm::f64vec3{vi.x + (double)x/(double)n, vi.y + (double)y/(double)n, vi.z + (double)z/(double)n}; };
    std::vector<CrvPt_FL*> T_crs;

    std::vector<glm::f64vec3>edges_ol;
    std::vector<HalfEdge*> edges_he;
    double t_max = -1, t_min = 1;
    for(auto&e: Poly_v){
        glm::f64vec3 p = e->p;
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
            CrvPt_FL *P = new CrvPt_FL{v, v3, 0};
            P->k2d = P->k2d_m = P->k2d_p = 1/glm::distance(CtrlPts[0], CtrlPts[1]);
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
                //VectorDigAlign(v2); VectorDigAlign(v3);
                t_max = std::max(t_max, t); t_min = std::min(t_min, t);
                CrvPt_FL *P = new CrvPt_FL(v2, v3, t);
                T_crs.push_back(P);
            }
        }
    }else if(type == PaintTool::FoldLine_test){
        int faceVerticesSize = edges_ol.size();
        for(int i = 0; i < faceVerticesSize; i++){
            bool hasPoint = MathTool::IsIntersect(edges_ol[i], edges_ol[(i+1)%faceVerticesSize], CtrlPts[0], CtrlPts[1], true);
            if(hasPoint){
                glm::f64vec3 NewPoint = MathTool::getIntersectionPoint(edges_ol[i], edges_ol[(i+1)%faceVerticesSize], CtrlPts[0], CtrlPts[1]);
                CrvPt_FL *P = new CrvPt_FL(NewPoint, i);
                T_crs.push_back(P);
            }
        }
    }

    for(int i = 0; i < (int)T_crs.size(); i++){
        for(int j = i + 1; j < (int)T_crs.size(); j++){
            if(T_crs[i]->s < T_crs[j]->s){
                auto tmp = T_crs[i];
                T_crs[i] = T_crs[j];
                T_crs[j] = tmp;
                i = 0;
            }
        }
    }

    for(auto&t: T_crs){
        Vertices.push_back(t);
    }
    for(auto&t: T_crs){
        std::vector<HalfEdge*> H_new;
        for(auto&he: Edges){
            H_new = he->Split(t, Edges);
            if(H_new.size() > 0)
            {
                break;
            }
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

            HalfEdge_FL *h1 = new HalfEdge_FL(InsertPoints[0], EdgeType::fl);
            HalfEdge_FL *h2 = new HalfEdge_FL(InsertPoints[1], EdgeType::fl);
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
            if(h2->v->s == t_max || h2->v->s == t_min)FoldingCurve.push_back(h2);
            if(h1->v->s == t_max || h1->v->s == t_min)FoldingCurve.push_back(h1);
            if(h1->v->s - h2->v->s > 0)FoldingCurve.push_back(h2);
            else FoldingCurve.push_back(h1);
            if(f == Faces.back())FoldingCurve.push_back(h1);
        }
    }
    for(int i = 0; i < (int)FoldingCurve.size(); i++){
        for(int j = i + 1; j < (int)FoldingCurve.size(); j++){
            if(FoldingCurve[i]->v->s < FoldingCurve[j]->v->s){
                auto tmp = FoldingCurve[i];
                FoldingCurve[i] = FoldingCurve[j];
                FoldingCurve[j] = tmp;
                i = 0;
            }
        }
    }
    FoldingCurve.erase(std::unique(FoldingCurve.begin(), FoldingCurve.end()), FoldingCurve.end());
    Faces.insert(Faces.end(), Faces_new.begin(), Faces_new.end()); // 連結

    return true;
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
                    vertices.push_back(res[0]->vertex);
                    for(auto&h: res){
                        if(h->face == f)_NewEdges.push_back(h);
                    }

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
            _NewEdges[0]->face->ReConnect(h1);
            fn->ReConnect(h2);
            faces.push_back(fn); edges.push_back(h1); edges.push_back(h2);
            return true;
        }
    }

    return false;
}

void FoldLine::drawRulingInAllAngles(std::vector<std::array<glm::f64vec3, 2>>& _Rulings){
    auto IsSame = [](double c, double b){
        return (abs(c - b) <= 1e-5)? true: false;
    };
    double a2 = 0;
    glm::f64vec3 e, e2, x;
    glm::f64vec3 N, r;
    glm::f64vec3 e_bef, befN;
    double veclen = 40;
    std::array<glm::f64vec3, 2> newVec;
    _Rulings.clear();
    double angle = 0.0;
    int flsize = FoldingCurve.size();
    while(angle <= 2*std::numbers::pi){
        double a = angle;
        for(int i = 1; i < flsize - 1; i++){

            x = FoldingCurve[i]->vertex->p3;
            e = glm::normalize(FoldingCurve[i-1]->vertex->p3 - x);
            e2 = glm::normalize(FoldingCurve[i+1]->vertex->p3 - x);
            double beta = glm::angle(e,e2);
            //beta = (1.0 - abs(FoldingCurve[i]->prev->r->Gradation/255.0)) * std::numbers::pi;
            //if(FoldingCurve[i]->next->r->Gradation < 0)IsVallyFoldLine = true;

            if(i != 1){
                glm::f64vec3 e_next = IsSame(beta, std::numbers::pi) ? FoldingCurve[i]->face->getNormalVec(): MathTool::ProjectionVector(e + e2, -e); e_next = glm::normalize(e_next);
                //double tau = glm::orientedAngle(e_bef,e_next,-e);
                double tau = std::acos(glm::dot(glm::normalize(glm::cross(e, e2)), befN));
                if(glm::dot(glm::cross(befN, glm::normalize(glm::cross(e, e2) ) ), e ) > 0)tau = 2.0 * std::numbers::pi - tau;
                if(tau < 0) tau += 2.0 * std::numbers::pi;
                if(tau > 2.0 * std::numbers::pi)tau -= 2.0*std::numbers::pi;
                a = (a2 - tau < 0) ? a2 - tau + 2.0 * std::numbers::pi: a2 - tau;
                if(!std::isfinite(a))std::cout << "a is not finite " << i << std::endl;
                if(!std::isfinite(a2))std::cout << "a' is not finite " << i << std::endl;
            }

            double phi3 = glm::angle(e2, glm::normalize(FoldingCurve[i]->pair->next->next->vertex->p3 - x));
            double phi4 = glm::angle(e, glm::normalize(FoldingCurve[i]->pair->next->next->vertex->p3 - x));
            double k = 2.0 * std::numbers::pi - phi3 - phi4;
            double phi1 = IsSame(glm::sin(beta)*glm::cos(a), glm::sin(k))? std::numbers::pi/2.0: atan2((glm::cos(k) - glm::cos(beta)),(glm::sin(beta)*glm::cos(a)- glm::sin(k)));
            if(phi1 < 0){phi1 = std::numbers::pi + phi1;}
            double phi2 = k - phi1;
            double th = 1e-4;
            if(abs((phi1 + phi3) - std::numbers::pi) <= th && abs((phi1 + phi3) - std::numbers::pi) <= th)std::cout << i << " : flat-foldable , " << glm::degrees(a) << std::endl;
            if(abs(phi1 - k/2.0) <= th && abs(phi2 - k/2.0) <= th)std::cout << i << " : equal design angles ,  " << glm::degrees(a)<<std::endl;
            if(abs(phi1 + phi4 - std::numbers::pi) <= th && abs(phi2 + phi3 - std::numbers::pi) <= th)std::cout << i << " : cotinuation , " << glm::degrees(a) << std::endl;
            N = glm::normalize(glm::rotate (a, -e)  * glm::f64vec4{glm::normalize(glm::cross(e2, e)),1});
            r = glm::rotate(phi1, -N) * glm::f64vec4{e,1};r = glm::normalize(r);//新しいruling方向
            newVec = std::array<glm::f64vec3, 2>{veclen * r + x, x};_Rulings.push_back(newVec);

            double sin_a = (sin(phi1)*sin(a)/sin(phi2) > 1) ? 1: (sin(phi1)*sin(a)/sin(phi2) < -1)? -1: (sin(phi1)*sin(a)/sin(phi2));
            glm::f64vec3 r2 = ProjectionVector(r,e2); r2 = glm::normalize(r2);
            glm::f64vec3 ax_beta = IsSame(beta, std::numbers::pi) ? FoldingCurve[i]->face->getNormalVec(): e;
            ax_beta = ProjectionVector(ax_beta, e2); ax_beta = glm::normalize(ax_beta);
            //double cos_a = (glm::dot(r2, ax_beta) > 1) ? 1: (glm::dot(r2, ax_beta) < -1)? -1: glm::dot(r2, ax_beta);//値が期待しているのと違う
            double cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
            a2 = (sin_a >= 0 && cos_a >= 0)? asin(sin_a): (sin_a >= 0 && cos_a < 0)?std::numbers::pi - asin(sin_a): (sin_a < 0 && cos_a < 0)? std::numbers::pi + abs(asin(sin_a)): 2.0*std::numbers::pi + asin(sin_a);
            if(i == (int)FoldingCurve.size() -2)break;
            //e_bef = IsSame(beta, std::numbers::pi) ? FoldingCurve[i]->face->getNormalVec(): (ProjectionVector(e,e2)); e_bef = glm::normalize(e_bef);
            befN = glm::normalize(glm::cross(e, e2));;
        }

        angle += 1e-4;
    }

}

void FoldLine::applyAAAMethod(std::vector<Vertex*>& Poly_v,  std::vector<Face*>& Faces, std::vector<HalfEdge*>& edges, double a){
    if(FoldingCurve.empty())return;
    double phi02 = glm::angle(glm::normalize(FoldingCurve[0]->pair->next->next->vertex->p - FoldingCurve.front()->vertex->p),
                              glm::normalize(FoldingCurve[0]->next->vertex->p - FoldingCurve.front()->vertex->p));

    double phim1 = glm::angle(glm::normalize(FoldingCurve.end()[-2]->vertex->p - FoldingCurve.back()->vertex->p),
                              glm::normalize(FoldingCurve.back()->prev->vertex->p - FoldingCurve.back()->vertex->p));
    //std::cout << "phi02 = " << glm::degrees(phi02) << " , phim1 = " << glm::degrees(phim1) << std::endl;
    _FoldingAAAMethod(a, phi02, phim1, Poly_v, Faces, edges);
}

//o原点, v修正したい頂点
Vertex *getClosestVertex(Vertex *v, Vertex* o,  std::vector<HalfEdge*> Edges){
    Vertex *V_max = nullptr;
    double t_max = -1;
    for(auto&e: Edges){
        if((e->vertex == o || e->vertex == v) || e->edgetype != EdgeType::ol)continue;
        if(abs(glm::distance(e->vertex->p, v->p) + glm::distance(e->vertex->p, o->p) - glm::distance(v->p, o->p)) < 1e-5){
            double t = glm::distance(e->vertex->p, o->p)/glm::distance(v->p, o->p);
            if(t > t_max){
                t_max = t;
                V_max = e->vertex;
            }
        }
    }
    return V_max;
}

glm::f64vec3 calcTargetDistanceOnPlane(Vertex *p, Vertex *o, Vertex *v1, Vertex *v2){
    Eigen::Matrix2d A;
    Eigen::Vector2d b;
    b(0) = p->p.x - o->p.x; b(1) = p->p.y - o->p.y;
    A(0,0) = v1->p.x - o->p.x; A(0,1) = v2->p.x - o->p.x;
    A(1,0) = v1->p.y - o->p.y; A(1,1) = v2->p.y - o->p.y;
    Eigen::Vector2d x = A.colPivHouseholderQr().solve(b);
    return  x(0) * (v1->p3 - o->p3) + x(1) * (v2->p3 - o->p3) + o->p3;
}

void FoldLine::_FoldingAAAMethod(double & a, double phi02, double phim1,  std::vector<Vertex*>& Poly_v,  std::vector<Face*>& Faces, std::vector<HalfEdge*>& edges){
    auto IsSame = [](double c, double b){
        return (abs(c - b) <= 1e-5)? true: false;
    };
    double a2 = 0, l;
    glm::f64vec3 e, e2, x, e_bef;
    glm::f64vec3 N, r, n, N0, befN;
    double veclen = 40;
    std::array<glm::f64vec3, 2> newVec;

    SingleRuling.clear();NewRuling2d.clear();
    int flsize = FoldingCurve.size();

    for(int i = 1; i < flsize - 1; i++){
        x = FoldingCurve[i]->vertex->p3;
        e = glm::normalize(FoldingCurve[i-1]->vertex->p3 - x);
        e2 = glm::normalize(FoldingCurve[i+1]->vertex->p3 - x);
        double beta = glm::angle(e,e2);

        if(i != 1){
            glm::f64vec3 e_next = IsSame(beta, std::numbers::pi) ? FoldingCurve[i]->face->getNormalVec(): MathTool::ProjectionVector(e + e2, -e); e_next = glm::normalize(e_next);
            double tau = std::acos(glm::dot(glm::normalize(glm::cross(e, e2)), befN));
            if(glm::dot(glm::cross(befN, glm::normalize(glm::cross(e, e2) ) ), e ) > 0)tau = 2.0 * std::numbers::pi - tau;
            //double tau = (glm::dot(e, glm::cross(e_bef, e_next)) > 0)? std::acos(glm::dot(e_bef, e_next)): 2.0*std::numbers::pi - abs(std::acos(glm::dot(e_bef, e_next)));
            //double tau = glm::orientedAngle(e_bef,e_next,-e);
            if(tau < 0) tau += 2.0 * std::numbers::pi;
            if(tau > 2.0 * std::numbers::pi)tau -= 2.0*std::numbers::pi;
            a = (a2 - tau < 0) ? a2 - tau + 2.0 * std::numbers::pi: a2 - tau;
            std::cout << "---------------------------------" << std::endl;
            if(!std::isfinite(a))std::cout << "a is not finite " << i << std::endl;
            if(!std::isfinite(a2))std::cout << "a' is not finite " << i << std::endl;
            std::cout << i << " : tau = " << glm::degrees(tau) << " , a = " <<  glm::degrees(a) << " , a' = " << glm::degrees(a2) << std::endl;

        }

        double phi3 = glm::angle(e2, glm::normalize(FoldingCurve[i]->pair->next->next->vertex->p3 - x));
        double phi4 = glm::angle(e, glm::normalize(FoldingCurve[i]->pair->next->next->vertex->p3 - x));
        double k = 2.0 * std::numbers::pi - phi3 - phi4;
        double _phi1 = IsSame(glm::sin(beta)*glm::cos(a), glm::sin(k))? std::numbers::pi/2.0: atan2((std::cos(k) - std::cos(beta)),(std::sin(beta)*std::cos(a)- std::sin(k)));

        //if(phi1 < 0){phi1 = std::numbers::pi + phi1;}
        double phi1 = (_phi1 < 0)? _phi1+std::numbers::pi: _phi1;
        double phi2 = k - phi1;
         std::cout <<"angle "<< i << " :  " << "phi1 = " <<  std::setprecision(10) <<  glm::degrees(phi1) << ", phi2 = " << glm::degrees(phi2) << ", phi3 = " <<
                 glm::degrees(phi3) << ", phi4 = " << glm::degrees(phi4)<< " : k = " <<  glm::degrees(k) << " , beta = " << glm::degrees(beta) << " , sum = " << phi1+phi2+phi3+phi4 - 2.0*std::numbers::pi<<  std::endl;
        N = glm::normalize(glm::rotate (a, -e)  * glm::f64vec4{glm::cross(e2, e),1.0});
        if(i == 1)N0 = N;
        r = glm::rotate(phi1, -N) * glm::f64vec4{e, 1.0};r = glm::normalize(r);//新しいruling方向
        n = glm::rotate(phi1, glm::f64vec3{0,0,-1.0})* glm::f64vec4{(glm::normalize(FoldingCurve[i-1]->vertex->p- FoldingCurve[i]->vertex->p)), 1.0};n = glm::normalize(n);//展開図のruling方向
        newVec = std::array<glm::f64vec3, 2>{veclen * e + x, x};SingleRuling.push_back(newVec);
        newVec = std::array<glm::f64vec3, 2>{veclen * N + x, x};SingleRuling.push_back(newVec);
        newVec = std::array<glm::f64vec3, 2>{veclen * r + x, x};SingleRuling.push_back(newVec);
        newVec = std::array<glm::f64vec3, 2>{veclen * (e + e2) + x, x};SingleRuling.push_back(newVec);

        glm::f64vec3 crossPoint;
        bool hasPointOnEdge = setPoint(Poly_v, n, FoldingCurve[i]->vertex->p, crossPoint);
        if(hasPointOnEdge){
            l = glm::distance(crossPoint, FoldingCurve[i]->vertex->p);
            FoldingCurve[i]->prev->vertex->p = crossPoint;
            FoldingCurve[i]->prev->vertex->p3 = l * r + FoldingCurve[i]->vertex->p3;
        }else std::cout << "cross point could not be founded in " << i << std::endl;
        double sin_a = (sin(phi1)*sin(a)/sin(phi2) > 1) ? 1: (sin(phi1)*sin(a)/sin(phi2) < -1)? -1: (sin(phi1)*sin(a)/sin(phi2));
        double cos_a = (cos(phi1) - cos(phi2)*cos(beta))/(sin(phi2)*sin(beta));
        a2 = (sin_a >= 0 && cos_a >= 0)? asin(sin_a): (sin_a >= 0 && cos_a < 0)?std::numbers::pi - asin(sin_a): (sin_a < 0 && cos_a < 0)? std::numbers::pi + abs(asin(sin_a)): 2.0*std::numbers::pi + asin(sin_a);
        if(i == (int)FoldingCurve.size() -2)break;
        e_bef = IsSame(beta, std::numbers::pi) ? FoldingCurve[i]->face->getNormalVec(): MathTool::ProjectionVector(e + e2, e2); e_bef = glm::normalize(e_bef);
        befN = glm::normalize(glm::cross(e, e2));
    }

    {//i = 0(端の面)
        Vertex *v_clst = getClosestVertex(FoldingCurve[0]->pair->next->next->vertex , FoldingCurve[0]->vertex, edges);
        if(v_clst == nullptr){
            x = FoldingCurve[0]->vertex->p3;
            e = glm::normalize(FoldingCurve[1]->vertex->p3 - x);       
            r = glm::rotate(phi02, N0) * glm::f64vec4{e,1};r = glm::normalize(r);//新しいruling方向
            l = glm::length(FoldingCurve[0]->pair->next->next->vertex->p - FoldingCurve[0]->vertex->p);
            FoldingCurve[0]->pair->next->next->vertex->p3 = l * r + FoldingCurve[0]->vertex->p3;
        }else{
            HalfEdge *h_clst = nullptr;
            for(auto&h: v_clst->halfedge){if(h->edgetype == EdgeType::ol)h_clst = h;}
            if(h_clst != nullptr){
                auto p = calcTargetDistanceOnPlane(FoldingCurve[0]->pair->next->next->vertex, h_clst->prev->vertex, v_clst, h_clst->prev->prev->vertex);
                FoldingCurve[0]->pair->next->next->vertex->p3 = p;
            }else {
                std::cout <<"can't find correct edge" << std::endl;
            }

        }
    }

    {
        Vertex *v_clst = getClosestVertex(FoldingCurve.back()->prev->vertex , FoldingCurve.back()->vertex, edges);
        if(v_clst == nullptr){
            x = FoldingCurve.back()->vertex->p3;
            e = glm::normalize(FoldingCurve.end()[-2]->vertex->p3 - x);
            auto ee = glm::normalize(FoldingCurve.end()[-2]->prev->vertex->p3- FoldingCurve.end()[-2]->vertex->p3);
            N = glm::normalize(glm::cross(ee, -e));
            r = glm::rotate(phim1, N) * glm::f64vec4{e,1};r = glm::normalize(r);//新しいruling方向
            l = glm::length(FoldingCurve.back()->prev->vertex->p - FoldingCurve.back()->vertex->p);
            FoldingCurve.back()->prev->vertex->p3 = l * r + x;
        }else{
            HalfEdge *h_clst = nullptr;
            for(auto&h: v_clst->halfedge){if(h->edgetype != EdgeType::ol)h_clst = h;}
            if(h_clst != nullptr){
                auto p = calcTargetDistanceOnPlane(FoldingCurve.back()->prev->vertex, v_clst, h_clst->next->vertex, h_clst->next->next->vertex);
               FoldingCurve.back()->prev->vertex->p3 = p;
            }else {
                std::cout <<"can't find correct edge" << std::endl;
            }
        }
    }

    std::cout <<"||||||||||||||||||||||||||||||" << std::endl;
}


void FoldLine::EdgeRecconection(std::vector<Vertex*>& Poly_V, std::vector<Face*>& Faces, std::vector<HalfEdge*>& Edges){
    int n = Poly_V.size();


    for(int i = 0; i < n; i++){
        Vertex *v1 = Poly_V[i], *v2 = Poly_V[(i + 1) % n];
        HalfEdge *h1 = nullptr, *h2 = nullptr;
        std::vector<PointOnLine> PoL;
        for(auto&e: Edges){
            if((e->vertex == v1) && e->edgetype == EdgeType::ol)h1 = e;
            if((e->vertex == v2) && e->edgetype == EdgeType::ol)h2 = e;

            if((e->vertex == v1 || e->vertex == v2) || e->edgetype != EdgeType::ol)continue;
            if(abs(glm::distance(e->vertex->p, v1->p) + glm::distance(e->vertex->p, v2->p) - glm::distance(v1->p, v2->p)) < 1e-5){
                double t = glm::distance(e->vertex->p, v2->p)/glm::distance(v1->p, v2->p);
                auto p = PointOnLine(t, e->vertex);
                PoL.push_back(p);
            }
        }
        std::sort(PoL.begin(), PoL.end());

        if(PoL.empty()){//polygon同士の頂点で一度つなぎなおす
            h1->prev = h2; h2->next = h1;
        }else{
            for(int j = 1; j < (int)PoL.size(); j++){
                HalfEdge *edge1 = nullptr, *edge2= nullptr;
                for(auto&h: PoL[j-1].v->halfedge){
                    if(h->edgetype == EdgeType::ol)edge1 = h;
                }
                for(auto&h: PoL[j].v->halfedge){
                    if(h->edgetype != EdgeType::ol)edge2 = h;
                }
                edge1->next = edge2; edge2->prev = edge1;
            }
            for(auto&e: PoL.front().v->halfedge){
                if(e->edgetype == EdgeType::ol)continue;
                e->prev = h2;
                h2->next = e;
            }
            for(auto&e: PoL.back().v->halfedge){
                if(e->edgetype != EdgeType::ol)continue;
                e->next = h1;
                h1->prev = e;
            }
        }
    }
    std::vector<bool> IsReconnected(Edges.size(), false);
    int faceIndex = 0;

    for(int i = 0; i < (int)IsReconnected.size(); i++){
        if(IsReconnected[i])continue;
        HalfEdge *h = Edges[i];
        do{
            int j = std::distance(Edges.begin(), std::find(Edges.begin(), Edges.end(), h));
            IsReconnected[j] = true;
            Edges[j]->face = Faces[faceIndex];
            h = h->next;
        }while(h != Edges[i]);
        faceIndex++;
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


double FoldLine::AngleIn2Edges(HalfEdge *p, HalfEdge *p2, bool Is3d){
  if(Is3d)return glm::angle(glm::normalize(p->next->vertex->p3 - p->vertex->p3), glm::normalize(p2->next->vertex->p3 - p2->vertex->p3));
  return glm::angle(glm::normalize(p->next->vertex->p - p->vertex->p), glm::normalize(p2->next->vertex->p - p2->vertex->p));
}

