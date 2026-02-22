#pragma once

#include <vector>
#include <map>

class SpaceMatrix
{
private:
    int n;
    // 使用 std::vector<std::map<int, double>> 进行装配，兼顾稀疏性和查找效率
    // 在装配阶段，map 可以自动处理索引排序和重复项加和
    std::vector<std::map<int, double>> A;

    std::vector<double> b;

public:
    SpaceMatrix(int n) : n(n), A(n), b(n, 0.0)
    {
    }

    // 设置/增加值
    void addToA(int i, int j, double val)
    {
        if (i < 0 || i >= n || j < 0 || j >= n)
            return;
        A[i][j] += val;
    }

    void setA(int i, int j, double val)
    {
        if (i < 0 || i >= n || j < 0 || j >= n)
            return;
        A[i][j] = val;
    }

    double getA(int i, int j) const
    {
        auto it = A[i].find(j);
        if (it != A[i].end())
            return it->second;
        return 0.0;
    }

    std::map<int, double> getA(int i) const
    {
        return A.at(i);
    }

    std::vector<std::map<int, double>>& getA()
    {
        return A;
    }

    void addTob(int i, double val)
    {
        if (i < 0 || i >= n)
            return;
        b[i] += val;
    }

    void setb(int i, double val)
    {
        if (i < 0 || i >= n)
            return;
        b[i] = val;
    }

    double getb(int i) const
    {
        if (i < 0 || i >= n)
            return 0;
        return b.at(i);
    }

    std::vector<double>& getb()
    {
        return b;
    }

    std::vector<double> getb() const
    {
        return b;
    }

    // 矩阵-向量乘法
    std::vector<double> multiply(const std::vector<double>& x) const
    {
        std::vector<double> res(n, 0.0);
        for (int i = 0; i < n; ++i)
        {
            for (auto const& [j, val] : A[i])
            {
                res[i] += val * x[j];
            }
        }
        return res;
    }

    // 获取对角线元素
    double diag(int i) const
    {
        auto it = A[i].find(i);
        if (it != A[i].end())
            return it->second;
        return 0.0;
    }

    void clear()
    {
        for (auto& row : A)
            row.clear();
    }

    size_t size() const
    {
        return n;
    }
};
