#!/bin/bash
cd /home/winter/HeatTransferFVM

echo "=========================================="
echo "NON-ORTHOGONAL MESH (OpenFOAM polyMesh)"
echo "=========================================="

for mode in 0 1 2; do
    echo ""
    echo "--- Testing MODE $mode ---"
    
    # 设置模式和网格
    sed -i "s/int testMode = [0-9];/int testMode = $mode;/" src/main.cpp
    sed -i 's|// mesh.createSquareMesh|mesh.createSquareMesh|' src/main.cpp
    sed -i 's|mesh.ReadFromFile|// mesh.ReadFromFile|' src/main.cpp
    sed -i 's|mesh.createSquareMesh|// mesh.createSquareMesh|' src/main.cpp
    sed -i 's|// mesh.ReadFromFile|mesh.ReadFromFile|' src/main.cpp
    
    cd build && make -j$(nproc) > /dev/null 2>&1
    ./HeatTransferFVM 2>&1 | grep -E "(MODE|Method|Maximum|converged in|Temperature|Min T|Max T|Avg T)" | head -20
    cd ..
done

echo ""
echo ""
echo "=========================================="
echo "ORTHOGONAL MESH (createSquareMesh)"
echo "=========================================="

for mode in 0 1 2; do
    echo ""
    echo "--- Testing MODE $mode ---"
    
    # 设置模式和网格
    sed -i "s/int testMode = [0-9];/int testMode = $mode;/" src/main.cpp
    sed -i 's|mesh.ReadFromFile|// mesh.ReadFromFile|' src/main.cpp
    sed -i 's|// mesh.createSquareMesh|mesh.createSquareMesh|' src/main.cpp
    
    cd build && make -j$(nproc) > /dev/null 2>&1
    ./HeatTransferFVM 2>&1 | grep -E "(MODE|Method|Maximum|converged in|Temperature|Min T|Max T|Avg T)" | head -20
    cd ..
done
