#ifndef FOLDLINE_H
#define FOLDLINE_H

#include <cmath>
#include <numbers>
#include <setrulings.h>
#include <utility>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <nlopt.hpp>

class FoldLine
{
public:
    FoldLine(PaintTool _type);
    std::vector<Eigen::Vector3d>CtrlPts;
    bool addCtrlPt(Eigen::Vector3d& p, int dim);
    bool delCtrlPt(Eigen::Vector3d& p, int dim, std::shared_ptr<OUTLINE>& outline);
    bool moveCtrlPt(Eigen::Vector3d& p, int movePtIndex);
    std::vector<Eigen::Vector3d> CurvePts;
    double getColor();
    bool RevisionCrosPtsPosition();
    bool Optimization_FlapAngle(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double wb, double wp, bool ConstFunc = true);
    std::vector<std::vector<Eigen::Vector3d>> Optimization_SmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v, bool IsConnectEndPoint);
    std::vector<std::vector<Eigen::Vector3d>> Optimization_PlanaritySrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void ReassignColor(std::vector<std::shared_ptr<Line>>& Rulings, ColorPoint& CP);
    void SimplifyModel(double tol);
    bool SimpleSmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void SortCurve(bool ascending = false);
    void reassinruling(std::shared_ptr<FoldLine>& parent);
    std::vector<Eigen::Vector3d> point;
    std::vector<Vertex4d> FoldingCurve;
    void applyAAAMethod(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double a, bool begincenter);
    std::vector<std::array<Eigen::Vector3d, 2>>  AllRulings;
    void drawRulingInAllAngles(std::vector<std::array<Eigen::Vector3d, 2>>& _Rulings);

private:
    double color;
    bool setCurve(int dim);
    PaintTool type;
    std::vector<std::shared_ptr<CrvPt_FL>> Points_On_Curve;
};

#endif // FOLDLINE_H
