#include "fitting_plane.h"
//https://programming-surgeon.com/script/fit-plane/
Eigen::Vector3d find_plane(const std::vector<Eigen::Vector3d>& X, Eigen::Vector3d& com){
    Eigen::MatrixXd A = Eigen::MatrixXd(3, X.size());
    for(int i = 0; i < (int)X.size(); i++){A.col(i) = X[i];}
    com = A.rowwise().mean();
    for(int i = 0; i < (int)X.size(); i++)A.col(i) -= com;
    A = A * A.transpose();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ES(A);//固有値問題
    Eigen::VectorXd e = ES.eigenvalues();
    Eigen::MatrixXd U = ES.eigenvectors();
    Eigen::VectorXd::Index minId;
    double min = e.minCoeff(&minId);
    Eigen::Vector3d N = (U.col(minId)).transpose();
    return N;
}

void _flatten_lsp(std::shared_ptr<FoldLine>& FldLine){
    std::vector<Eigen::Vector3d> data;
    std::vector<std::shared_ptr<Vertex4d>> ValidFC;
    for(auto&v: FldLine->FoldingCurve){
        if(v->IsCalc){
            data.push_back(v->first->p3);
            ValidFC.push_back(v);
        }
    }
    //N,o: 面の法線ベクトルと点、V,p:rulingのベクトルと通る点
    auto PointOnPlaneAndLine = [](Eigen::Vector3d N, const Eigen::Vector3d& o, Eigen::Vector3d V, const Eigen::Vector3d& p){
        N = N.normalized(); V = V.normalized();
        double t = N.dot(o - p)/N.dot(V);
        return t*V + p;
    };
    int mid = data.size()/2;
    Eigen::Vector3d com, o = ValidFC[mid]->first->p3;
    Eigen::Vector3d N = find_plane(data, com);
    Eigen::Vector3d V, p;
    for(int i = 0; i < (int)ValidFC.size(); i++){
        V = ValidFC[i]->first->p3 - ValidFC[i]->third->p3;
        p = ValidFC[i]->third->p3;
        Eigen::Vector3d PointOnPlane = PointOnPlaneAndLine(N, com, V, p);
        ValidFC[i]->first->p3 = PointOnPlane;
        ValidFC[i]->first->p = (PointOnPlane - p).norm()*(ValidFC[i]->first->p - ValidFC[i]->third->p).normalized() + ValidFC[i]->third->p;
    }
    return;
    for(int i = mid-1; i >= 0; i--){
        V = ValidFC[i]->first->p3 - ValidFC[i]->third->p3;
        p = ValidFC[i]->third->p3;
        Eigen::Vector3d PointOnPlane = PointOnPlaneAndLine(N, o, V, p);
        ValidFC[i]->first->p3 = PointOnPlane;
        ValidFC[i]->first->p = (PointOnPlane - p).norm()*(ValidFC[i]->first->p - ValidFC[i]->third->p).normalized() + ValidFC[i]->third->p;
    }

    //左側
    for(int i = mid+1; i < (int)ValidFC.size(); i++){
        V = ValidFC[i]->first->p3 - ValidFC[i]->third->p3;
        p = ValidFC[i]->third->p3;
        Eigen::Vector3d PointOnPlane = PointOnPlaneAndLine(N, o, V, p);
        ValidFC[i]->first->p3 = PointOnPlane;
        ValidFC[i]->first->p = (PointOnPlane - p).norm()*(ValidFC[i]->first->p - ValidFC[i]->third->p).normalized() + ValidFC[i]->third->p;
    }
}
