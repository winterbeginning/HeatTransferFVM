#ifndef _SurfaceField_
#define _SurfaceField_

#include <vector>
#include <string>
#include "Mesh.hpp"

template <typename valType>
class SurfaceField
{
public:
    std::string name;
    const Mesh& mesh;

    std::vector<valType> internalFaceField;

    SurfaceField(const Mesh& mesh, valType initVal = valType{}) : mesh(mesh)
    {
        internalFaceField.assign(mesh.nInternalFace, initVal);
    }

    // 获取特定单元的值 (辅助函数)
    valType& operator[](int idx)
    {
        return internalFaceField[idx];
    }
    const valType& operator[](int idx) const
    {
        return internalFaceField[idx];
    }

    // 获取场的大小
    size_t size() const
    {
        return internalFaceField.size();
    }
};

#endif