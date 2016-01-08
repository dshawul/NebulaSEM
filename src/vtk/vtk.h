#ifndef __VTK_H
#define __VTK_H

#include "field.h"

namespace Vtk {
    void write_vtk(Int);
    extern bool write_polyhedral;
    extern bool write_cell_value;
}

#endif
