#ifndef FOLDLINE_H
#define FOLDLINE_H

#include <cmath>
#include <numbers>
#include <setrulings.hpp>
#include <utility>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <nlopt.hpp>

class FoldLine : public std::enable_shared_from_this<FoldLine>
{
public:
    FoldLine(PaintTool _type);
    std::vector<Eigen::Vector3d>CtrlPts;
    std::vector<Eigen::Vector3d> point;
    std::vector<Vertex4d> FoldingCurve;
    std::vector<Eigen::Vector3d> CurvePts;
    std::vector<std::array<Eigen::Vector3d, 2>>  AllRulings;
    double a_flap;

    bool isbend();
    bool addCtrlPt(Eigen::Vector3d& p, int dim);
    bool delCtrlPt(Eigen::Vector3d& p, int dim, std::shared_ptr<OUTLINE>& outline);
    bool moveCtrlPt(Eigen::Vector3d& p, int movePtIndex, int dim);

    double getColor();
    bool RevisionCrosPtsPosition();
    bool Optimization_FlapAngle(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double wb, double wp, bool ConstFunc = true);
    std::vector<std::vector<Eigen::Vector3d>> Optimization_SmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v, bool IsConnectEndPoint);
    std::vector<std::vector<Eigen::Vector3d>> Optimization_PlanaritySrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void ReassignColor(std::vector<std::shared_ptr<Line>>& Rulings, ColorPoint& CP);
    void SimplifyModel(double tol);
    bool SimpleSmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v, const std::vector<std::shared_ptr<FoldLine>>& FL);
    void SortCurve(bool ascending = false);
    void reassignruling(std::shared_ptr<FoldLine>& parent);

    void applyAAAMethod(const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool begincenter, double a = -1); 
    void drawRulingInAllAngles(std::vector<std::array<Eigen::Vector3d, 2>>& _Rulings);

private:
    double color;
    bool setCurve(int dim);
    int curveNum;
    PaintTool type;
    std::vector<std::shared_ptr<CrvPt_FL>> Points_On_Curve;

};


#endif // FOLDLINE_H
