#ifndef FITTING_PLANE_H
#define FITTING_PLANE_H

#include <iostream>
#include <Eigen/Dense>
#include "foldline.h"

Eigen::Vector3d find_plane(const std::vector<Eigen::Vector3d>& X, Eigen::Vector3d& com);
void _flatten_lsp(std::shared_ptr<FoldLine>& FldLine);

#endif // FITTING_PLANE_H
