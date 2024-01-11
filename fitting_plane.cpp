#include "fitting_plane.h"
//https://programming-surgeon.com/script/fit-plane/
Eigen::Vector3d find_plane(const std::vector<Eigen::Vector3d>& X, Eigen::Vector3d& com){
    Eigen::MatrixXd A = Eigen::MatrixXd(3, X.size());
    for(int i = 0; i < (int)X.size(); i++){A.col(i) = X[i];}
    com = A.rowwise().mean();
    for(int i = 0; i < (int)X.size(); i++)A.col(i) -= com;
    A = A * A.transpose();
    //Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ES(A);//固有値問題
    Eigen::VectorXd e = ES.eigenvalues();
    Eigen::MatrixXd U = ES.eigenvectors();
    Eigen::VectorXd::Index minId;
    double min = e.minCoeff(&minId);
    Eigen::Vector3d N = (U.col(minId)).transpose();
    return N;
}

void flatten_lsp(std::shared_ptr<FoldLine>& FldLine){
    std::vector<Eigen::Vector3d> data;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&v: FldLine->FoldingCurve){
        if(v->IsCalc){
            data.push_back(v->first->p3);
            ValidFC.push_back(v);
        }
    }
    Eigen::Vector3d com, o, e, e2;
    Eigen::Vector3d N = find_plane(data, com);
    std::shared_ptr<Vertex> v, p;
    int mid = data.size()/2;
    for(int i = mid-1; i > 0; i--){
        v = ValidFC[i-1]->first;
        p = ValidFC[i-1]->third;
        o = ValidFC[i+1]->first->p3; e = (ValidFC[i]->first->p3 - o), e2 = (ValidFC[i+2]->first->p3 - o);
        Eigen::Vector3d PointOnPlane = MathTool::CrossPointLineAndPlane(e, e2, o , p->p3, v->p3);
        ValidFC[i-1]->first->p3 = PointOnPlane;
        ValidFC[i-1]->first->p = (PointOnPlane - p->p3).norm()*(v->p - p->p).normalized() + p->p;
    }

    //左側
    for(int i = mid+1; i < (int)ValidFC.size()-1; i++){
        v = ValidFC[i+1]->first;
        p = ValidFC[i+1]->third;
        o = ValidFC[i-1]->first->p3; e = (ValidFC[i]->first->p3 - o).normalized(), e2 = (ValidFC[i-2]->first->p3 - o).normalized();
        Eigen::Vector3d PointOnPlane = MathTool::CrossPointLineAndPlane(e, e2, o , p->p3, v->p3);
        ValidFC[i+1]->first->p3 = PointOnPlane;
        ValidFC[i+1]->first->p = (PointOnPlane - p->p3).norm()*(v->p - p->p).normalized() + p->p;
    }
}
