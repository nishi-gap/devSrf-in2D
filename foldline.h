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
#include <QDebug>
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
    bool Optimization_FlapAngle(const std::vector<std::shared_ptr<Vertex>>& Poly_V, double wb, double wp, bool ConstFunc = true);
    std::vector<std::vector<Eigen::Vector3d>> Optimization_SmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v, bool IsConnectEndPoint);
    std::vector<std::vector<Eigen::Vector3d>> Optimization_PlanaritySrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void ReassignColor();
    void TrimLines(int size);
    void SimplifyModel(double tol, bool isroot);
    void SimplifyModel(int iselim, bool isroot);
    bool SimpleSmooothSrf(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void SortCurve(bool ascending = false);
    void reassignruling(std::shared_ptr<FoldLine>& parent);

    void applyAAAMethod(const std::vector<std::shared_ptr<Vertex>>& Poly_V, bool begincenter, double a, double _tol, bool isroot);
    void revisecrossedruling(const std::vector<std::shared_ptr<Vertex>>& Poly_v);
    void drawRulingInAllAngles(std::vector<std::array<Eigen::Vector3d, 2>>& _Rulings);


private:
    double color;
    bool setCurve(int dim);
    int curveNum;
    PaintTool type;
    std::vector<std::shared_ptr<CrvPt_FL>> Points_On_Curve;


};


#endif // FOLDLINE_H
