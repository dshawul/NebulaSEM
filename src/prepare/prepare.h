#ifndef __PREPARE_H
#define __PREPARE_H

#include "field.h"
#include "vtk.h"

namespace Prepare {
    void convertVTK(std::vector<std::string>&,Int,Int);
    void probe(std::vector<std::string>&,Int,Int);
}

#endif
