#include "module.h"

double Net::calcHPWL() const
{
    if (_termList.empty()) return 0.0;

    double minX = _termList[0]->getCenterX();
    double maxX = minX;
    double minY = _termList[0]->getCenterY();
    double maxY = minY;

    for (size_t i = 1; i < _termList.size(); ++i) {
        double cx = _termList[i]->getCenterX();
        double cy = _termList[i]->getCenterY();
        if (cx < minX) minX = cx;
        if (cx > maxX) maxX = cx;
        if (cy < minY) minY = cy;
        if (cy > maxY) maxY = cy;
    }
    return (maxX - minX) + (maxY - minY);
}
