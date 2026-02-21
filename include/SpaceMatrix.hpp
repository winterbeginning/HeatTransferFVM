#pragma once

#include <vector>
#include <map>

class SpaceMatrix
{
public:
    int n;
    // 使用 std::vector<std::map<int, double>> 进行装配，兼顾稀疏性和查找效率
    // 在装配阶段，map 可以自动处理索引排序和重复项加和
    std::vector<std::map<int, double>> data;

    SpaceMatrix(int n) : n(n), data(n)
    {
    }

    // 设置/增加值
    void addValue(int i, int j, double val)
    {
        if (i < 0 || i >= n || j < 0 || j >= n)
            return;
        data[i][j] += val;
    }

    void setValue(int i, int j, double val)
    {
        if (i < 0 || i >= n || j < 0 || j >= n)
            return;
        data[i][j] = val;
    }

    double getValue(int i, int j) const
    {
        auto it = data[i].find(j);
        if (it != data[i].end())
            return it->second;
        return 0.0;
    }

    // 矩阵-向量乘法
    std::vector<double> multiply(const std::vector<double>& x) const
    {
        std::vector<double> res(n, 0.0);
        for (int i = 0; i < n; ++i)
        {
            for (auto const& [j, val] : data[i])
            {
                res[i] += val * x[j];
            }
        }
        return res;
    }

    // 获取对角线元素
    double diag(int i) const
    {
        auto it = data[i].find(i);
        if (it != data[i].end())
            return it->second;
        return 0.0;
    }

    void clear()
    {
        for (auto& row : data)
            row.clear();
    }
};
