#ifndef __WRAPPER_
#define __WRAPPER_

#include <fstream>

/**
  Solver wrappers
 */
namespace Solver {

    extern std::ifstream input;

    void Initialize(int argc, char* argv[]);
    void Finalize();
}

#endif
