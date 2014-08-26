#include "../volumeElimination/volumeElimination.h"
#include <vil/vil_image_view.h>

#ifndef _GRIDSV_H_
#define _GRIDSV_H_

template<class T>
class gridSV
{
public:
	gridSV(vil_image_view<T> &im, int voxelsize);
	int getLabel(int i, int j, int k);
private:
	int _width;
	int _height;
	int _layers;
	volumeElimination::vector3i _label;
	double gridSV<T>::cost(double v1, double v2);
};

#include "gridSV.cpp"

#endif