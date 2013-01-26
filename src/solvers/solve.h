#ifndef __SOLVE_H
#define __SOLVE_H

#include "field.h"

void Solve(const MeshMatrix<Scalar>&);
void Solve(const MeshMatrix<Vector>&); 
void Solve(const MeshMatrix<STensor>&); 
void Solve(const MeshMatrix<Tensor>&); 

#endif
