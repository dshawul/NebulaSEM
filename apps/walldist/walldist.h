#ifndef __WALLDIST_H
#define __WALLDIST_H

#include <istream>
#include "my_types.h"

namespace Mesh {
    void calc_walldist(Int,Int = 1);
}
void walldist(std::istream&);

#endif
