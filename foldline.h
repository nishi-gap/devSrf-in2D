#ifndef FOLDLINE_H
#define FOLDLINE_H

#include "setrulings.h"
#include <fstream>
#include <sstream>
#include <filesystem>
#include <nlopt.hpp>
#include <QDebug>
#include <thread>

template <typename T>
class PointOnEndEdge{
public:
    std::shared_ptr<T> v;
    double t;
    PointOnEndEdge(std::shared_ptr<T>&_v, double _t): v(_v), t(_t){}
};

class FoldLine : public std::enable_shared_from_this<FoldLine>
{
public:
    FoldLine(PaintTool _type);
    std::vector<Eigen::Vector3d>CtrlPts;
    std::vector<Eigen::Vector3d> point;
    std::vector<std::shared_ptr<Vertex4d>> FoldingCurve;
    std::vector<Eigen::Vector3d> CurvePts;
    std::vector<std::array<Eigen::Vector3d, 2>>  AllRulings;
    PaintTool type;
    double a_flap;
    double tol;
    int validsize;

    std::shared_ptr<FoldLine> deepCopy();
    bool isbend();
    bool addCtrlPt(Eigen::Vector3d& p, int dim);
    bool delCtrlPt(Eigen::Vector3d& p, int dim, std::shared_ptr<OUTLINE>& outline);
    bool moveCtrlPt(Eigen::Vector3d& p, int movePtIndex, int dim);
    bool setCurve(int dim);
    double getColor();
    bool RevisionCrosPtsPosition();

    bool Optimization_EndPoint(const std::vector<std::shared_ptr<Vertex>>& Poly_V);
    bool Optimization_FlapAngle(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double wb, double wp, int rank, int alg, bool IsStartEnd, int OptimizationAlgorithm, bool OptimizeAngleFor3Rulings);
    bool Optimization_Vertex(const std::vector<std::shared_ptr<Vertex>>& Poly_V);
    bool PropagateOptimization_Vertex(const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool IsStartEnd, int VertexMoveAlg, int OptimizationAlgorithm, double range, double warea, double wsim);
    void ReviseCenterFlapAngle(const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool IsStartEnd, int AlgType);
    void _movevertex(double t, const std::vector<std::shared_ptr<Vertex>>& Poly_V);

    std::vector<std::vector<Eigen::Vector3d>> Optimization_SmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v, bool IsConnectEndPoint);
    std::vector<std::vector<Eigen::Vector3d>> Optimization_PlanaritySrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void ReassignColor();
    void TrimLines(int size);
    //void SimplifyModel(double tol, bool isroot);
    void SimplifyModel(int iselim, bool isroot);
    bool SimpleSmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void SortCurve(bool ascending = false);
    void AlignmentVertex4dDirection();
    void CheckIsCrossedRulings();
    void initialize_foldstate(bool IsStartEnd, const std::vector<std::shared_ptr<Vertex>>& Poly_V);
    void reassignruling(std::shared_ptr<FoldLine>& parent, const std::vector<std::shared_ptr<Line>>& Surface, const std::vector<std::shared_ptr<Line>>& Rulings);

    void applyAAAMethod(const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool IsStartEnd, double a);
    void revisecrossedruling(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void drawRulingInAllAngles(std::vector<std::array<Eigen::Vector3d, 2>>& _Rulings);

    void test_rotate(double a);
    void FittingEndPoint_flattencurve(const Eigen::Vector3d& initRight, const Eigen::Vector3d& initLeft);
    std::vector<std::vector<std::shared_ptr<Vertex>>> CalclateRegressionCurve(double a, const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool IsWriteCSV, bool IsStartEnd, std::vector<std::vector<std::shared_ptr<Vertex>>>& Tri_fixside);
private:
    double color;

    int curveNum;

    std::vector<std::shared_ptr<CrvPt_FL>> Points_On_Curve;
    bool ReviseVertexPos(const std::vector<std::shared_ptr<Vertex>>& Poly_V, int EndIndex_left, int EndIndex_right, int AlgOptim, double range, double warea, double wsim);
    //Eigen::MatrixXd GlobalSplineInterpolation(std::vector<double>&Knot, bool is3d, int dim);
};


class NTreeNode{
public:
    std::shared_ptr<FoldLine> data;
    std::vector<std::shared_ptr<NTreeNode>> children;
    NTreeNode(const std::shared_ptr<FoldLine>& val): data(val){}
    bool IsChild(const std::shared_ptr<FoldLine>& val){
        std::queue<std::shared_ptr<NTreeNode>> q; q.push(std::shared_ptr<NTreeNode>(this));
        while (!q.empty()) {
            std::shared_ptr<NTreeNode> cur = q.front(); q.pop();
            for (std::shared_ptr<NTreeNode> child : cur->children){
                if(child->data == val)return true;
                q.push(child);
            }
        }
        return false;
    }
};


class NTree {
private:
    std::shared_ptr<NTreeNode> root;

public:
    NTree(const FoldLine& val) { root = std::make_shared<NTreeNode>(NTreeNode(std::make_shared<FoldLine>(val)));}
    NTree(const std::shared_ptr<FoldLine>& val) { root = std::make_shared<NTreeNode>(NTreeNode(val));}
    NTree(){root = std::make_shared<NTreeNode>(nullptr);}

    bool empty(){return (root == nullptr)? true: false;}
    void insert(const std::shared_ptr<FoldLine>& parentVal, const std::shared_ptr<FoldLine>& val){
        std::shared_ptr<NTreeNode> newNode = std::make_shared<NTreeNode>(val);
        insertRecursive(root, parentVal, newNode);
    }
    std::shared_ptr<NTreeNode> GetRoot(){return (root != nullptr)? root: std::shared_ptr<NTreeNode>(nullptr);}

    std::shared_ptr<NTreeNode> getParent(const std::shared_ptr<FoldLine>& val){
        if (root == nullptr)return std::shared_ptr<NTreeNode>(nullptr);
        std::queue<std::shared_ptr<NTreeNode>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<NTreeNode> cur = q.front(); q.pop();
            for (std::shared_ptr<NTreeNode> child : cur->children){
                if(child->data == val)return cur;
                q.push(child);
            }
        }
        return std::shared_ptr<NTreeNode>(nullptr);
    }

    std::vector<std::shared_ptr<NTreeNode>> GetChildren(const std::shared_ptr<NTreeNode>& parent){
        return parent->children;
    }

    void insertRecursive(const std::shared_ptr<NTreeNode>& node, const std::shared_ptr<FoldLine>& parentVal, const std::shared_ptr<NTreeNode>& newNode){
        if (node == nullptr) return;
        if (node->data == parentVal) {
            node->children.push_back(newNode);
            return;
        }
        for (const std::shared_ptr<NTreeNode>& child : node->children) {
            insertRecursive(child, parentVal, newNode);
        }
    }

    void erase(const std::shared_ptr<FoldLine>& val){
        std::shared_ptr<NTreeNode> Tree = find(val);
        if(Tree == nullptr)return;
        std::shared_ptr<NTreeNode> par = getParent(val);
        for(const auto&child: Tree->children)par->children.push_back(child);
        //delete Tree;//内部の変数のアドレスは解放されていない
    }

    void clear(){
        if (root == nullptr)return;
        std::queue<std::shared_ptr<NTreeNode>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<NTreeNode> cur = q.front(); q.pop();
            for (const std::shared_ptr<NTreeNode>& child : cur->children)q.push(child);
            cur.reset();
        }
        root = nullptr;
    }

    void changeRoot(const std::shared_ptr<FoldLine>& val){
        std::shared_ptr<NTreeNode> newNode = std::make_shared<NTreeNode>(NTreeNode(val));
        std::shared_ptr<NTreeNode> tmp = root;
        root = newNode;
        root->children.push_back(tmp);
    }
    std::shared_ptr<NTreeNode> find(const std::shared_ptr<FoldLine>& val){
        if (root == nullptr)return nullptr;
        std::queue<std::shared_ptr<NTreeNode>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<NTreeNode> cur = q.front(); q.pop();
            if(cur->data == val)return cur;
            for (std::shared_ptr<NTreeNode> child : cur->children){
                if(child != nullptr)q.push(child);
            }
        }
        return nullptr;
    }
    void printTree(const std::shared_ptr<NTreeNode>& node, int depth = 0){
        if (node == nullptr)return;
        for (int i = 0; i < depth; ++i)std::cout << "*";
        std::cout << node->data << std::endl;
        for (const std::shared_ptr<NTreeNode>& child : node->children)printTree(child, depth + 1);
    }
    void print(){printTree(root);}
    int getLayerNum(){
        int rank = 0;
        if (root == nullptr)return rank;
        std::queue<std::shared_ptr<NTreeNode>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<NTreeNode> cur = q.front(); q.pop();
            for (std::shared_ptr<NTreeNode> child : cur->children){
                if(child != nullptr)q.push(child);
            }
            rank++;
        }
        return rank;
    }
};

#endif // FOLDLINE_H
