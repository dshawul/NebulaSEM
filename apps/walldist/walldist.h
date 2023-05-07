#ifndef __WALLDIST_H
#define __WALLDIST_H

#include <istream>
#include "types.h"

namespace Mesh {
    void calc_walldist(Int,Int = 1);
}
void walldist(std::istream&);

#endif
