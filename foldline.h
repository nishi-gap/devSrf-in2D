#ifndef FOLDLINE_H
#define FOLDLINE_H

#include "setrulings.h"
#include <fstream>
#include <sstream>
#include <filesystem>
#include <nlopt.hpp>
#include <QDebug>

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
    double a_flap;
    double tol;
    int validsize;

    bool isbend();
    bool addCtrlPt(Eigen::Vector3d& p, int dim);
    bool delCtrlPt(Eigen::Vector3d& p, int dim, std::shared_ptr<OUTLINE>& outline);
    bool moveCtrlPt(Eigen::Vector3d& p, int movePtIndex, int dim);

    double getColor();
    bool RevisionCrosPtsPosition();
    bool Optimization_FlapAngle(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double wb, double wp, int rank, int alg, bool IsStartEnd);
    bool Optimization_Vertex(const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool IsStartEnd);
    std::vector<std::vector<Eigen::Vector3d>> Optimization_SmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v, bool IsConnectEndPoint);
    std::vector<std::vector<Eigen::Vector3d>> Optimization_PlanaritySrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void ReassignColor();
    void TrimLines(int size);
    //void SimplifyModel(double tol, bool isroot);
    void SimplifyModel(int iselim, bool isroot);
    bool SimpleSmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void SortCurve(bool ascending = false);
    void reassignruling(std::shared_ptr<FoldLine>& parent, const std::vector<std::shared_ptr<Line>>& Surface, const std::vector<std::shared_ptr<Line>>& Rulings);

    void applyAAAMethod(const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool IsStartEnd, double a, double _tol, bool isroot);
    void revisecrossedruling(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void drawRulingInAllAngles(std::vector<std::array<Eigen::Vector3d, 2>>& _Rulings);

    std::vector<std::vector<std::shared_ptr<Vertex>>> CalclateRegressionCurve(double a, const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool IsWriteCSV, bool IsStartEnd, std::vector<std::vector<std::shared_ptr<Vertex>>>& Tri_fixside);
private:
    double color;
    bool setCurve(int dim);
    int curveNum;
    PaintTool type;
    std::vector<std::shared_ptr<CrvPt_FL>> Points_On_Curve;
    //Eigen::MatrixXd GlobalSplineInterpolation(std::vector<double>&Knot, bool is3d, int dim);
};

#endif // FOLDLINE_H
