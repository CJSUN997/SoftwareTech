#ifndef PARAMETERIZATION_H
#define PARAMETERIZATION_H
#include <vector>
#include "Eigen/Sparse"
#include <String>
#include "Mesh.h"
#include <iostream>
class Mesh;
class Parameterization 
{


public:
	Parameterization(Mesh & mesh) :obj(mesh){};

	~Parameterization();
	void calculate();

private:
	// ========== function =============
	void boundaryVerticesParameterize();
	void innerVerticesParameterize();
	
	void ComputeShapPreservingWeight(int currIndex, std::vector<float> &vertexWeight);
	bool sortNeighborVertex(std::vector<int> &neighborVertexIndices, std::vector<int> &localNeighborVertexIndices, int index);

	// ========= values ===============

	Mesh & obj;
};

#endif // !PARAMETERIZATION_H