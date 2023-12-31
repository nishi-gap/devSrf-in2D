#ifndef TRANSFORM_H
#define TRANSFORM_H
#include<QPointF>
#include "foldline.h"

//折り目
void AffinTrans(std::shared_ptr<FoldLine>& Crease, const QPointF& befPos, const QPointF& curPos);
void AffinScale(std::shared_ptr<FoldLine>& Crease,const QPointF& basePos, const QPointF& befPos, const QPointF& curPos);
void AffinRotate(std::shared_ptr<FoldLine>& Crease, const QPointF& basePos, const QPointF& befPos, const QPointF& curPos);


#endif // TRANSFORM_H
