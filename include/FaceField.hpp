#ifndef _FaceField_
#define _FaceField_

#include <vector>
#include "Mesh.hpp"

template <typename ValueType>
class FaceField
{
public:
    const Mesh& mesh;

    std::vector<ValueType> field;

    FaceField(const Mesh& mesh, ValueType initVal = ValueType{}) : mesh(mesh)
    {
        field.assign(mesh.facePoints.size(), initVal);
    }

    // 获取特定面的值
    ValueType& operator[](int idx)
    {
        return field[idx];
    }
    const ValueType& operator[](int idx) const
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