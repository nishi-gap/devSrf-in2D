#ifndef MATHTOOL_H
#define MATHTOOL_H

#include <glm/gtx/transform.hpp>
#include <glm/gtx/string_cast.hpp>
#include <vector>
#include <iostream>
#include <tuple>
#include <Eigen/Dense>
#include <utility>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numbers>
#include <queue>

enum class EdgeType{
    none,
    ol,//outline
    r,//ruling(from curve line)
    cl,//curve line
    fl,//fold line
    r_fl,//ruling(fold line)
};

enum class CurveType{
    none,
    bezier3,
    bsp3,
    line,
    arc,
    DebugTest,
};

struct ColorPoint{
    double color, angle;
    ColorPoint(): color(200), angle(std::numbers::pi/2.0){}
    ColorPoint(double _c, double _a): color(_c), angle(_a){}
};

enum class PaintTool{
    None,
    Reset,
    deform,

    AddCurve,
    Bezier_r,
    Bspline_r,
    Line_r,
    Arc_r,

    SetColor,
    NewGradationMode,

    Rectangle_ol,
    Polygon_ol,
    Polyline_ol,
    EditVertex_ol,
    Move_ol,
    Const_ol,
    ConnectVertices_ol,

    MoveCtrlPt,
    InsertCtrlPt,
    DeleteCtrlPt,
    DeleteCurve,

    FoldLine_bezier,
    FoldLine_line,
    FoldLine_arc,
    FoldLine_test,
    FoldlineColor,
    FoldLine_move,

    CheckDevelopability,

};

namespace MathTool{

double distP2L(glm::f64vec3 la, glm::f64vec3 lb, glm::f64vec3& p, glm::f64vec3& q);//点と線分の距離, s:laからlbへの比率(垂線が内部にあれば0 ~ 1)

bool IsIntersect(glm::f64vec3&p1, glm::f64vec3&p2, glm::f64vec3&p3, glm::f64vec3&p4, bool ConsiderEnd = false);
glm::f64vec3 getIntersectionPoint(glm::f64vec3& p1, glm::f64vec3& p2, glm::f64vec3& p3, glm::f64vec3& p4);
glm::f64vec3 calcCrossPoint_2Vector(glm::f64vec3 p1, glm::f64vec3 q1, glm::f64vec3 p2, glm::f64vec3 q2);

void Triangulation(std::vector<glm::f64vec3>&input, std::vector<std::array<glm::f64vec3, 3>>&output);
bool hasPointInTriangle3D(glm::f64vec3 p, std::array<glm::f64vec3, 3>& V);
bool IsAngleLessThan180(glm::f64vec3& o, glm::f64vec3& a, glm::f64vec3& b);

std::vector<double> _bezierclipping(const std::vector<glm::f64vec3>&CtrlPts_base, std::vector<glm::f64vec3>&CtrlPts_cur, std::array<glm::f64vec3, 2>& line, int dim);//交点が一つのみの場合
bool is_point_on_line(glm::f64vec3 p, glm::f64vec3 lp1, glm::f64vec3 lp2);
std::pair<std::vector<glm::f64vec3>, std::vector<glm::f64vec3>> BezierSplit(const std::vector<glm::f64vec3>& CtrlPts, double t);
std::vector<std::vector<glm::f64vec3>> de_casteljau_algorithm(const std::vector<glm::f64vec3>& CtrlPts, double t);

std::vector<glm::f64vec3> GrahamScan(std::vector<glm::f64vec3>& Q);//凸包の計算
double SignedArea(glm::f64vec3 a, glm::f64vec3 b, glm::f64vec3 p);

glm::f64vec3 bspline(std::vector<glm::f64vec3>&CtrlPts, double t, int dim, std::vector<double>Knot);

double factorial(int n);
double cmb(int n, int i);
double BernsteinBasisFunc(int n, int i, double t);
glm::f64vec3 bezier(std::vector<glm::f64vec3>& CtrlPts, double t, int dim);
void diffBezier(std::vector<glm::f64vec3>& dP, std::vector<glm::f64vec3>& CtrlPts, double t, int dim);

double basis(int n, int i, int p, double u, std::vector<double>& U);

std::vector<double> LSM_apply(std::vector<double>&y, int dim = 1);

std::vector<glm::f64vec3> getPlaneFromCurvePoints(std::vector<glm::f64vec3>& Points, std::vector<glm::f64vec3>& BasisVectors);
glm::f64vec3 ProjectionVector(glm::f64vec3 v, glm::f64vec3 n, bool Isnormalize = false);

}

template <typename T>
class NTreeNode{
public:
  T data;
  std::vector<NTreeNode<T>*> children;
  NTreeNode(const T& val): data(val){}
  void releaseAll(){

  }

};

template <typename T>
class NTree {
private:
    NTreeNode<T>* root;

public:
    NTree(const T& val) { root = new NTreeNode<T>(val);}
    NTree(){root = nullptr;}

    bool empty(){return (root == nullptr)? true: false;}
    void insert(const T& parentVal, const T& val){
        NTreeNode<T>* newNode = new NTreeNode<T>(val);
        insertRecursive(root, parentVal, newNode);
    }
    NTreeNode<T>* GetRoot(){return (root != nullptr)? root: nullptr;}

    NTreeNode<T>* getParent(const T& val){
        if (root == nullptr)return nullptr;
        std::queue<NTreeNode<T>*> q;
        q.push(root);
        while (!q.empty()) {
            NTreeNode<T>* cur = q.front(); q.pop();
            if(cur->data == val)return cur;
            for (NTreeNode<T>* child : cur->children)  q.push(child);
        }
        return nullptr;
    }

    std::vector<NTreeNode<T>*> GetChildren(NTreeNode<T>* parent){
        return parent->children;
    }

    void insertRecursive(NTreeNode<T>* node, const T& parentVal, NTreeNode<T>* newNode){
        if (node == nullptr) return;
        if (node->data == parentVal) {
            node->children.push_back(newNode);
            return;
        }
        for (NTreeNode<T>* child : node->children) {
            insertRecursive(child, parentVal, newNode);
        }
    }

    void erase(const T& val){
        NTreeNode<T>* Tree = find(val);
        if(Tree == nullptr)return;
        NTreeNode<T>* par = getParent(val);
        for(const auto&child: Tree->children)par->children.push_back(child);
        delete Tree;//内部の変数のアドレスは解放されていない
    }

    void clear(){
        if (root == nullptr)return;
        std::queue<NTreeNode<T>*> q;
        q.push(root);
        while (!q.empty()) {
            NTreeNode<T>* cur = q.front(); q.pop();
            for (NTreeNode<T>* child : cur->children)q.push(child);
            delete cur;
        }
        root = nullptr;
    }

    void changeRoot(const T& val){
        NTreeNode<T>* newNode = new NTreeNode<T>(val);
        NTreeNode<T>* tmp = root;
        root = newNode;
        root->children.push_back(tmp);
    }
    NTreeNode<T>* find(const T& val){
        if (root == nullptr)return nullptr;
        std::queue<NTreeNode<T>*> q;
        q.push(root);
        while (!q.empty()) {
            NTreeNode<T>* cur = q.front(); q.pop();
            if(cur->data == val)return cur;
            for (NTreeNode<T>* child : cur->children){
                if(child != nullptr)q.push(child);
            }
        }
        return nullptr;
    }
    void printTree(NTreeNode<T>* node, int depth = 0){
        if (node == nullptr)return;
        for (int i = 0; i < depth; ++i)std::cout << "*";
        std::cout << node->data << std::endl;
        for (NTreeNode<T>* child : node->children)printTree(child, depth + 1);
    }
    void print(){printTree(root);}
};

#endif // MATHTOOL_H
