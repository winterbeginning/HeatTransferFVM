#ifndef _fvMatrix_
#define _fvMatrix_

#include "SpaceMatrix.hpp"
#include "Field.hpp"
#include "Solver.hpp"

template <typename ValueType>
class fvMatrix : public SpaceMatrix<ValueType>
{
public:
    Field<ValueType>& psi;     // 关联的场
    Solver<ValueType>* solver; // 绑定的求解器

    fvMatrix(Field<ValueType>& psi)
        : SpaceMatrix<ValueType>(psi.mesh.numCells), psi(psi), solver(nullptr)
    {
    }

    void setSolver(Solver<ValueType>& s)
    {
        this->solver = &s;
    }

    void solve()
    {
        if (solver)
        {
            psi.internalField = solver->solve(*this, psi.internalField);
        }
        else
        {
            std::cerr << "Error: No solver bound to matrix!" << std::endl;
        }
    }

    // 支持直接调用 solve(solver)
    void solve(Solver<ValueType>& s)
    {
        psi.internalField = s.solve(*this, psi.internalField);
    }
};

#endif
