#ifndef MATHTOOL_H
#define MATHTOOL_H

#include <cmath>
#include <iomanip>
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
#include <memory>

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

namespace DebugMode{
    class Singleton{
    public:
        static Singleton& getInstance();
        bool isdebug() const;
        void switchval();
        int valuecanged(int n);
        void init();
    private:
        Singleton(){}
        Singleton(const Singleton&) = delete;
        void operator=(const Singleton&) = delete;
        bool val;
    };
}

namespace MathTool{
double rad2deg(double a);
double distP2L(const Eigen::Vector3d& la, const Eigen::Vector3d& lb, const Eigen::Vector3d& p, Eigen::Vector3d& q);//点と線分の距離, s:laからlbへの比率(垂線が内部にあれば0 ~ 1)

bool IsIntersect(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3, const Eigen::Vector3d& p4, bool ConsiderEnd = false);
Eigen::Vector3d getIntersectionPoint(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2,const Eigen::Vector3d p3, const Eigen::Vector3d& p4);
Eigen::Vector3d calcCrossPoint_2Vector(const Eigen::Vector3d& p1, const Eigen::Vector3d& q1, const Eigen::Vector3d& p2, const Eigen::Vector3d& q2);

bool hasPointInTriangle3D(const Eigen::Vector3d& p, std::array<Eigen::Vector3d, 3>& V);
bool IsAngleLessThan180(Eigen::Vector3d& o, Eigen::Vector3d& a, Eigen::Vector3d& b);

std::vector<double> _bezierclipping(const std::vector<Eigen::Vector3d>&CtrlPts_base, std::vector<Eigen::Vector3d>&CtrlPts_cur, std::array<Eigen::Vector3d, 2>& line, int dim);//交点が一つのみの場合
bool is_point_on_line(Eigen::Vector3d p, Eigen::Vector3d lp1, Eigen::Vector3d lp2);
std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>> BezierSplit(const std::vector<Eigen::Vector3d>& CtrlPts, double t);
std::vector<std::vector<Eigen::Vector3d>> de_casteljau_algorithm(const std::vector<Eigen::Vector3d>& CtrlPts, double t);

std::vector<Eigen::Vector3d> GrahamScan(const std::vector<Eigen::Vector3d>& Q);//凸包の計算
double SignedArea(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& p);

Eigen::Vector3d bspline(std::vector<Eigen::Vector3d>&CtrlPts, double t, int dim, std::vector<double>Knot);

double factorial(int n);
double cmb(int n, int i);
double BernsteinBasisFunc(int n, int i, double t);
Eigen::Vector3d bezier(std::vector<Eigen::Vector3d>& CtrlPts, double t, int dim);

double basis(int n, int i, int p, double u, std::vector<double>& U);

std::vector<double> LSM_apply(std::vector<double>&y, int dim = 1);
Eigen::Vector3d ProjectionVector(const Eigen::Vector3d& v, Eigen::Vector3d n, bool Isnormalize = false);

template<typename T>
void swap(T &a, T& b){ T c = a; a = b; b = c;}

}

template <typename T>
class NTreeNode{
public:
  T data;
  std::vector<std::shared_ptr<NTreeNode<T>>> children;
  NTreeNode(const T& val): data(val){}
};

template <typename T>
class NTree {
private:
    std::shared_ptr<NTreeNode<T>> root;

public:
    NTree(const T& val) { root = std::make_shared<NTreeNode<T>>(NTreeNode<T>(val));}
    NTree(){root = std::shared_ptr<NTreeNode<T>>(nullptr);}

    bool empty(){return (root == nullptr)? true: false;}
    void insert(const T& parentVal, const T& val){
        std::shared_ptr<NTreeNode<T>> newNode = std::make_shared<NTreeNode<T>>(NTreeNode<T>(val));
        insertRecursive(root, parentVal, newNode);
    }
    std::shared_ptr<NTreeNode<T>> GetRoot(){return (root != nullptr)? root: std::shared_ptr<NTreeNode<T>>(nullptr);}

    std::shared_ptr<NTreeNode<T>> getParent(const T& val){
        if (root == nullptr)return std::shared_ptr<NTreeNode<T>>(nullptr);
        std::queue<std::shared_ptr<NTreeNode<T>>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<NTreeNode<T>> cur = q.front(); q.pop();
            if(cur->data == val)return cur;
            for (std::shared_ptr<NTreeNode<T>> child : cur->children)  q.push(child);
        }
        return std::shared_ptr<NTreeNode<T>>(nullptr);
    }

    std::vector<std::shared_ptr<NTreeNode<T>>> GetChildren(const std::shared_ptr<NTreeNode<T>>& parent){
        return parent->children;
    }

    void insertRecursive(const std::shared_ptr<NTreeNode<T>>& node, const T& parentVal, const std::shared_ptr<NTreeNode<T>>& newNode){
        if (node == nullptr) return;
        if (node->data == parentVal) {
            node->children.push_back(newNode);
            return;
        }
        for (const std::shared_ptr<NTreeNode<T>>& child : node->children) {
            insertRecursive(child, parentVal, newNode);
        }
    }

    void erase(const T& val){
        std::shared_ptr<NTreeNode<T>> Tree = find(val);
        if(Tree == nullptr)return;
        std::shared_ptr<NTreeNode<T>> par = getParent(val);
        for(const auto&child: Tree->children)par->children.push_back(child);
        delete Tree;//内部の変数のアドレスは解放されていない
    }

    void clear(){
        if (root == nullptr)return;
        std::queue<std::shared_ptr<NTreeNode<T>>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<NTreeNode<T>> cur = q.front(); q.pop();
            for (const std::shared_ptr<NTreeNode<T>>& child : cur->children)q.push(child);
            cur.reset();
        }
        root = nullptr;
    }

    void changeRoot(const T& val){
        std::shared_ptr<NTreeNode<T>> newNode = std::make_shared<NTreeNode<T>>(NTreeNode<T>(val));
        std::shared_ptr<NTreeNode<T>> tmp = root;
        root = newNode;
        root->children.push_back(tmp);
    }
    std::shared_ptr<NTreeNode<T>> find(const T& val){
        if (root == nullptr)return nullptr;
        std::queue<std::shared_ptr<NTreeNode<T>>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<NTreeNode<T>> cur = q.front(); q.pop();
            if(cur->data == val)return cur;
            for (std::shared_ptr<NTreeNode<T>> child : cur->children){
                if(child != nullptr)q.push(child);
            }
        }
        return nullptr;
    }
    void printTree(const std::shared_ptr<NTreeNode<T>>& node, int depth = 0){
        if (node == nullptr)return;
        for (int i = 0; i < depth; ++i)std::cout << "*";
        std::cout << node->data << std::endl;
        for (const std::shared_ptr<NTreeNode<T>>& child : node->children)printTree(child, depth + 1);
    }
    void print(){printTree(root);}
};

#endif // MATHTOOL_H
