#include "Parameterization.h"
#include <fstream>      // std::ifstream


Parameterization::~Parameterization()
{
}


void Parameterization::calculate()
{
	this->boundaryVerticesParameterize();
	this->innerVerticesParameterize();
}


void Parameterization::boundaryVerticesParameterize()
{

	/*************************/
	/* insert your code here */
	/*************************/
}


void Parameterization::innerVerticesParameterize()
{

	/*************************/
	/* insert your code here */
	/*************************/
}





void Parameterization::ComputeShapPreservingWeight(int index, std::vector<float> &vertexWeight)
{

}


bool Parameterization::sortNeighborVertex(std::vector<int> &neighborVertexIndices, std::vector<int> &localNeighborVertexIndices, int index)
{
	return false;
}

