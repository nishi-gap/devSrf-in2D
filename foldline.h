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
    FoldLine(){}
    FoldLine(PaintTool _type);
    FoldLine(const FoldLine& _FL);
    FoldLine& operator = (FoldLine& _FL){
        CtrlPts = _FL.CtrlPts;
        curveNum = _FL.curveNum;
        type = _FL.type;
        CurvePts = _FL.CurvePts;
        a_flap = _FL.a_flap;
        validsize = _FL.validsize;
        FoldingCurve.clear();
        for(const auto& v4d: _FL.FoldingCurve)FoldingCurve.push_back(v4d);

        return *this;
    }
    std::shared_ptr<FoldLine> deepCopy();

    std::vector<Eigen::Vector3d>CtrlPts;
    std::vector<std::shared_ptr<Vertex4d>> FoldingCurve;
    std::vector<Eigen::Vector3d> CurvePts;
    PaintTool type;
    double a_flap;
    int validsize;


    bool isbend();
    std::shared_ptr<FoldLine> addCtrlPt(Eigen::Vector3d& p, int dim);
    bool delCtrlPt(Eigen::Vector3d& p, int dim, std::shared_ptr<OUTLINE>& outline);
    bool moveCtrlPt(Eigen::Vector3d& p, int movePtIndex, int dim);
    bool setCurve(int dim);

    bool Optimization_FlapAngle(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double wp, double wsim, int rank);
    bool PropagateOptimization_Vertex(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double wp, double wsim);

    void ReassignColor();

    void SimplifyModel(int iselim, bool isroot);
    bool SimpleSmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void SortCurve(bool ascending = false);
    void AlignmentVertex4dDirection();
    void initialize_foldstate(const std::vector<std::shared_ptr<Vertex>>& Poly_V);
    void reassignruling(std::shared_ptr<FoldLine>& parent, const std::vector<std::shared_ptr<Line>>& Surface, const std::vector<std::shared_ptr<Line>>& Rulings);

    void applyAAAMethod(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double a);
private:
    int curveNum;
    std::vector<std::shared_ptr<CrvPt_FL>> Points_On_Curve;
};


template<typename T> class TreeNode : public std::enable_shared_from_this<TreeNode<T>> {
public:
    T data;
    std::vector<std::shared_ptr<TreeNode>> children;
    TreeNode(const T& val) :data(val) {}
    TreeNode() = default;
    bool IsChild(const T& val) {
        if(val == data)return false;
        std::queue<std::shared_ptr<TreeNode>> q;
        q.push(std::shared_ptr<TreeNode>(this));
        while (!q.empty()) {
            std::shared_ptr<TreeNode> cur = q.front(); q.pop();
            for (std::shared_ptr<TreeNode> child : cur->children) {
                if (child->data == val)return true;
                q.push(child);
            }
        }
        return false;
    }

    // 帰りがけ順の実装
    std::vector<T> postOrderTraversal() {
        std::vector<T> result;
        for (auto& child : children) {
            auto childResult = child->postOrderTraversal();
            result.insert(result.end(), childResult.begin(), childResult.end());
        }
        result.push_back(data);
        return result;
    }
};

template<typename T> class NTree : public std::enable_shared_from_this<NTree<T>> {
private:
    std::shared_ptr<TreeNode<T>> root;

    void insertRecursive(const std::shared_ptr<TreeNode<T>>& node, const T& parentVal, const std::shared_ptr<TreeNode<T>>& newNode) {
        if (node == nullptr) return;
        if (node->data == parentVal) {
            node->children.push_back(newNode);
            return;
        }
        for (const std::shared_ptr<TreeNode<T>>& child : node->children)  insertRecursive(child, parentVal, newNode);
    }

    std::shared_ptr<TreeNode<T>> copyNode(const std::shared_ptr<TreeNode<T>>& node) {
        if (node == nullptr)return nullptr;
        auto NewNode = std::make_shared<TreeNode<T>>(node->data);
        for (auto& child : node->children)child = copyNode(child);
        return NewNode;
    }

public:
    NTree(const T& val) { root = std::make_shared<TreeNode<T>>(val); }
    NTree() { root = std::make_shared<TreeNode<T>>(); }
    NTree(const NTree& _T){ root = _T.GetRoot();}
    NTree& operator=(const NTree& _T){root = _T.GetRoot(); return *this;}

    //rootは必ず要素が空とする、折り目はrootのchildren以降に入れる
    void AddLast(const T& val) {
        std::shared_ptr<TreeNode<T>> newNode = std::make_shared<TreeNode<T>>(val);
        if (root == nullptr)  root = std::make_shared<TreeNode<T>>();
        std::queue<std::shared_ptr<TreeNode<T>>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<TreeNode<T>> cur = q.front(); q.pop();
            if (cur->children.empty()) { cur->children.push_back(newNode); return; }
            for (std::shared_ptr<TreeNode<T>> child : cur->children) q.push(child);
        }
    }

    void Clear() {
        if (root == nullptr)return;
        std::queue<std::shared_ptr<TreeNode<T>>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<TreeNode<T>> cur = q.front(); q.pop();
            for (const std::shared_ptr<TreeNode<T>>& child : cur->children)q.push(child);
            cur.reset();
        }
        root = nullptr;
    }

    void ChangeRoot(const T& val) {
        std::shared_ptr<TreeNode<T>> newNode = std::make_shared<TreeNode<T>>(val);
        std::shared_ptr<TreeNode<T>> tmp = root;
        root = newNode;
        root->children.push_back(tmp);
    }

    void DeepCopy(const NTree& _T) { root = copyNode(_T.root); }

    void Erase(const T& val) {
        std::shared_ptr<TreeNode<T>> Node = find(val);
        if (Node == nullptr)return;
        std::shared_ptr<TreeNode<T>> par = GetParent(val);
        for (const auto& child : Node->children)par->children.push_back(child);
    }

    std::shared_ptr<TreeNode<T>> find(const T& val) {
        if (root == nullptr)return nullptr;
        std::queue<std::shared_ptr<TreeNode<T>>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<TreeNode<T>> cur = q.front(); q.pop();
            if (cur->data == val)return cur;
            for (std::shared_ptr<TreeNode<T>> child : cur->children) {
                if (child != nullptr)q.push(child);
            }
        }
        return nullptr;
    }

    std::vector<std::shared_ptr<TreeNode<T>>> GetChildren(const std::shared_ptr<TreeNode<T>>& parent) { return parent->children; }

    std::shared_ptr<TreeNode<T>> GetRoot() const{ return (root != nullptr) ? root : std::shared_ptr<TreeNode<T>>(nullptr); }

    std::shared_ptr<TreeNode<T>> GetParent(const T& val) {
        if (root == nullptr)return std::shared_ptr<TreeNode<T>>(nullptr);
        std::queue<std::shared_ptr<TreeNode<T>>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<TreeNode<T>> cur = q.front(); q.pop();
            for (std::shared_ptr<TreeNode<T>> child : cur->children) {
                if (child->data == val)return cur;
                q.push(child);
            }
        }
        return std::shared_ptr<TreeNode<T>>(nullptr);
    }

    std::shared_ptr<TreeNode<T>> GetParent(const std::shared_ptr<TreeNode<T>>& node) {
        if (root == nullptr)return nullptr;
        std::queue<std::shared_ptr<TreeNode<T>>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<TreeNode<T>> cur = q.front(); q.pop();
            for (std::shared_ptr<TreeNode<T>> child : cur->children) {
                if (child == node)return cur;
                q.push(child);
            }
        }
        return nullptr;
    }

    void Insert(const T& parentVal, const T& val) {
        if(find(val) != nullptr)return;//木の要素で重複が起きないようにする
        std::shared_ptr<TreeNode<T>> newNode = std::make_shared<TreeNode<T>>(val);
        insertRecursive(root, parentVal, newNode);
    }

    void Insert(const std::shared_ptr<TreeNode<T>>& node, const T& val){
        if(find(val) != nullptr)return;//木の要素で重複が起きないようにする
        std::shared_ptr<TreeNode<T>> newNode = std::make_shared<TreeNode<T>>(val);
        insertRecursive(root, node->data, newNode);
    }

    //parの一つ下の階層にchildがいる場合のみ使う(parとchildの間に要素を挿入して par -> _data -> childの関係にする)
    void InsertNode(const T& par, const T& child, const T& _data) {
        if(par == _data || child == _data)return;
        auto parNode = find(par);
        if (parNode == nullptr)return;
        auto itr = std::find_if(parNode->children.begin(), parNode->children.end(), [&](const std::shared_ptr<TreeNode<T>>&  t){ return t->data == child; });
        if (itr == parNode->children.end())return;
        insert(par, _data);
        Planting(_data, child);
    }

    std::vector<std::shared_ptr<TreeNode<T>>> NTree2Array() {
        std::vector<std::shared_ptr<TreeNode<T>>> List;
        if (root == nullptr) return List;
        std::queue<std::shared_ptr<TreeNode<T>>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<TreeNode<T>> cur = q.front(); q.pop();
            List.push_back(cur);
            for (std::shared_ptr<TreeNode<T>> child : cur->children)  q.push(child);

        }
        return List;
    }

    //par: 新たに親となる要素, tar: 挿し木を行う要素
    void Planting(const T&  par, const T& tar) {
        auto cutNode = GetParent(tar), childNode = find(tar);
        if (cutNode == nullptr || childNode == nullptr || (cutNode == childNode))return;

        auto it = std::find(cutNode->children.begin(), cutNode->children.end(), childNode);
        if (it != cutNode->children.end())  cutNode->children.erase(it);
        auto newparent = find(par);
        newparent->children.push_back(childNode);
    }

    void Planting(const T&  par, std::shared_ptr<TreeNode<T>>& tar) {
        auto cutNode = GetParent(tar), newparent = find(par);
        if (cutNode == nullptr || tar == nullptr)return;

        auto it = std::find(cutNode->children.begin(), cutNode->children.end(), tar);
        if (it != cutNode->children.end())  cutNode->children.erase(it);
        newparent->children.push_back(tar);
    }

    void SwapNode(const T& n1, const T& n2) {
        std::shared_ptr<TreeNode<T>> branch1 = find(n1), branch2 = find(n2);
        if (branch1 == nullptr || branch2 == nullptr)return;
        branch1->data = n2; branch2->data = n1;
    }

    void printTree(const std::shared_ptr<TreeNode<T>>& node, int depth = 0) {
        if (node == nullptr)return;
        for (int i = 0; i < depth; ++i)std::cout << "*";
        std::cout << node->data << std::endl;
        for (const std::shared_ptr<TreeNode<T>>& child : node->children)printTree(child, depth + 1);
    }

    // 帰りがけ順の実行
    std::vector<T> postOrderTraversal(const T& val) {
        auto node = find(val);
        if (node != nullptr)  return node->postOrderTraversal();
        return std::vector<T>();
    }

    void print_inverse() {
        std::vector<T> result = root->postOrderTraversal();
        for (auto& data : result)std::cout << data << std::endl;
        std::cout << "///////////////////" << std::endl;
    }

    void print() { printTree(root); std::cout << "///////////////////" << std::endl; }

    void print_bwf() {
        if (root == nullptr)return;
        std::queue<std::shared_ptr<TreeNode<T>>> q;
        q.push(root);
        while (!q.empty()) {
            std::shared_ptr<TreeNode<T>> cur = q.front(); q.pop();
            if(!cur->children.empty())std::cout << cur->data << " -> ";
            for (std::shared_ptr<TreeNode<T>> child : cur->children) {
                std::cout << child->data << " ";
                q.push(child);
            }
            if (!cur->children.empty())std::cout << std::endl;
        }
        std::cout << "//////////////////////" << std::endl;
    }

};

#endif // FOLDLINE_H
