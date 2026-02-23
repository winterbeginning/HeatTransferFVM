#ifndef _SpaceMatrix_
#define _SpaceMatrix_

#include <vector>
#include <map>

template <typename ValueType>
class SpaceMatrix
{
public:
    int n;
    // 每次组装用的临时存储，结束后会被压缩为 CSR 格式
    std::vector<std::map<int, ValueType>> assembly_A;
    std::vector<ValueType> b;

    // CSR 数据结构
    std::vector<ValueType> csr_val;
    std::vector<int> csr_col;
    std::vector<int> csr_ptr;

    SpaceMatrix(int n) : n(n), assembly_A(n), b(n, ValueType{})
    {
    }

    // 核心组装接口
    void addToA(int i, int j, ValueType value)
    {
        if (i < 0 || i >= n || j < 0 || j >= n)
            return;
        assembly_A[i][j] += value;
    }

    void addTob(int i, ValueType value)
    {
        if (i < 0 || i >= n)
            return;
        b[i] += value;
    }

    // 将 map 形式压缩为 CSR 形式
    void compress()
    {
        csr_val.clear();
        csr_col.clear();
        csr_ptr.assign(n + 1, 0);

        int current_nnz = 0;
        for (int i = 0; i < n; ++i)
        {
            csr_ptr[i] = current_nnz;
            for (auto const& [j, value] : assembly_A[i])
            {
                csr_val.push_back(value);
                csr_col.push_back(j);
                current_nnz++;
            }
        }
        csr_ptr[n] = current_nnz;
    }

    // 清理数据，但不改变 n
    void clear()
    {
        for (auto& row : assembly_A)
            row.clear();
        std::fill(b.begin(), b.end(), 0.0);
        csr_val.clear();
        csr_col.clear();
        csr_ptr.clear();
    }

    // CSR 矩阵-向量乘法
    std::vector<ValueType> multiply(const std::vector<ValueType>& x) const
    {
        std::vector<ValueType> res(n, ValueType{});
        for (int i = 0; i < n; ++i)
        {
            for (int k = csr_ptr[i]; k < csr_ptr[i + 1]; ++k)
            {
                res[i] += csr_val[k] * x[csr_col[k]];
            }
        }
        return res;
    }

    // 获取对角线元素（基于 CSR 查找）
    ValueType diag(int i) const
    {
        for (int k = csr_ptr[i]; k < csr_ptr[i + 1]; ++k)
        {
            if (csr_col[k] == i)
                return csr_val[k];
        }
        return ValueType{};
    }

    int size() const
    {
        return n;
    }
};

#endif
