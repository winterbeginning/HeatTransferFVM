#!/bin/bash

echo "========================================="
echo "非正交修正对比测试"
echo "========================================="

cd /home/winter/HeatTransferFVM

# 创建结果存储数组
declare -A results

# 函数：运行测试并提取结果
run_test() {
    local mesh_type=$1
    local non_orth=$2
    local test_name=$3
    
    echo ""
    echo "========================================="
    echo "测试: $test_name"
    echo "========================================="
    
    # 创建main.cpp
    cat > src/main.cpp << EOF
#include <iostream>
#include "Mesh.hpp"
#include "FiniteVolume.hpp"
#include "Solver.hpp"
#include <algorithm>

int main()
{
    std::cout << "--- Steady State Heat Conduction ---" << std::endl;
    Mesh mesh;
EOF

    if [ "$mesh_type" == "poly" ]; then
        cat >> src/main.cpp << 'EOF'
    mesh.ReadFromFile("/home/winter/HeatTransferFVM/polyMesh");
EOF
    else
        cat >> src/main.cpp << 'EOF'
    mesh.createSquareMesh(50, 50, 1.0, 1.0);
EOF
    fi

    cat >> src/main.cpp << EOF
    std::cout << "--- Mesh Cells " << mesh.numCells << " ---" << std::endl;

    FiniteVolume fvm(mesh);
    fvm.setSolveOption(false, true, false);  // 只开启扩散
    fvm.NonOrthogonalCorrection = $non_orth;  // 设置非正交修正开关
    fvm.properties.setProperties(1.0, 1.0, 1.0, 1.0);

    Field<double>& T = fvm.T;
    T.fill(0.0);
    T.setBoundary("left", 0.0, 0.0, 1.0);
    T.setBoundary("down", 0.0, 0.0, 0.0);
    T.setBoundary("right", 0.0, 0.0, 0.0);
    T.setBoundary("top", 100.0, 0.0, 1.0);
EOF

    if [ "$mesh_type" == "poly" ]; then
        cat >> src/main.cpp << 'EOF'
    // 三维楔形网格的Z方向边界
    T.setBoundary("Base", 0.0, 0.0, 0.0);  // Neumann: 绝热
    T.setBoundary("Top", 0.0, 0.0, 0.0);   // Neumann: 绝热
EOF
    fi

    cat >> src/main.cpp << 'EOF'

    Solver solver(5000, 1e-6, false, SolverType::GAUSS_SEIDEL);  // 关闭输出
    fvm.solve(TimeScheme::STEADY, solver);

    // 输出温度统计
    double maxT = *std::max_element(T.internalField.begin(), T.internalField.end());
    double minT = *std::min_element(T.internalField.begin(), T.internalField.end());
    double avgT = 0.0;
    for(auto t : T.internalField) avgT += t;
    avgT /= T.internalField.size();
    
    std::cout << "最大温度: " << maxT << " °C" << std::endl;
    std::cout << "最小温度: " << minT << " °C" << std::endl;
    std::cout << "平均温度: " << avgT << " °C" << std::endl;
    
    return 0;
}
EOF

    # 编译并运行
    cmake --build build > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "编译失败!"
        return 1
    fi
    
    local output=$(./build/HeatTransferFVM 2>&1)
    echo "$output" | tail -5
    
    # 提取结果
    local max_temp=$(echo "$output" | grep "最大温度" | awk '{print $2}')
    local avg_temp=$(echo "$output" | grep "平均温度" | awk '{print $2}')
    
    results["${test_name}_max"]=$max_temp
    results["${test_name}_avg"]=$avg_temp
}

# 测试1: 非结构网格 + 非正交修正开启
run_test "poly" "true" "非结构网格+非正交修正ON"

# 测试2: 非结构网格 + 非正交修正关闭
run_test "poly" "false" "非结构网格+非正交修正OFF"

# 测试3: 正交网格 + 非正交修正开启
run_test "square" "true" "正交网格+非正交修正ON"

# 测试4: 正交网格 + 非正交修正关闭
run_test "square" "false" "正交网格+非正交修正OFF"

# 输出对比表格
echo ""
echo "========================================="
echo "对比结果汇总"
echo "========================================="
printf "%-30s | %12s | %12s\n" "测试配置" "最大温度(°C)" "平均温度(°C)"
echo "--------------------------------------------------------------------"
printf "%-30s | %12s | %12s\n" "非结构网格+非正交修正ON" "${results[非结构网格+非正交修正ON_max]}" "${results[非结构网格+非正交修正ON_avg]}"
printf "%-30s | %12s | %12s\n" "非结构网格+非正交修正OFF" "${results[非结构网格+非正交修正OFF_max]}" "${results[非结构网格+非正交修正OFF_avg]}"
printf "%-30s | %12s | %12s\n" "正交网格+非正交修正ON" "${results[正交网格+非正交修正ON_max]}" "${results[正交网格+非正交修正ON_avg]}"
printf "%-30s | %12s | %12s\n" "正交网格+非正交修正OFF" "${results[正交网格+非正交修正OFF_max]}" "${results[正交网格+非正交修正OFF_avg]}"
echo "--------------------------------------------------------------------"

# 计算差异
echo ""
echo "关键观察:"
echo "1. 非结构网格: 非正交修正的影响"
if [ -n "${results[非结构网格+非正交修正ON_avg]}" ] && [ -n "${results[非结构网格+非正交修正OFF_avg]}" ]; then
    diff=$(echo "${results[非结构网格+非正交修正ON_avg]} - ${results[非结构网格+非正交修正OFF_avg]}" | bc)
    echo "   平均温度差异: $diff °C"
fi

echo ""
echo "2. 正交网格: 非正交修正应该无影响(结果应该完全相同)"
if [ -n "${results[正交网格+非正交修正ON_avg]}" ] && [ -n "${results[正交网格+非正交修正OFF_avg]}" ]; then
    diff=$(echo "${results[正交网格+非正交修正ON_avg]} - ${results[正交网格+非正交修正OFF_avg]}" | bc 2>/dev/null || echo "0")
    echo "   平均温度差异: $diff °C (应该≈0)"
fi

echo ""
echo "========================================="
echo "测试完成!"
echo "========================================="
