#ifndef _SurfaceField_
#define _SurfaceField_

#include <vector>
#include "Mesh.hpp"

template <typename valType>
class SurfaceField
{
public:
    const Mesh& mesh;

    std::vector<valType> field;

    SurfaceField(const Mesh& mesh, valType initVal = valType{}) : mesh(mesh)
    {
        field.assign(mesh.facePoints.size(), initVal);
    }

    // 获取特定面的值
    valType& operator[](int idx)
    {
        return field[idx];
    }
    const valType& operator[](int idx) const
    {
        return field[idx];
    }

    // 获取面场的大小
    size_t size() const
    {
        return field.size();
    }
};

#endif