#ifndef __PREPARE_H
#define __PREPARE_H

#include "field.h"
#include "vtk.h"

namespace Prepare {
    int convertVTK(std::vector<std::string>&,Int);
    int probe(std::vector<std::string>&,Int);
}

#endif
