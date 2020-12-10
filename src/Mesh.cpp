#include <list>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Mesh.h"
#include <stdio.h>
#include <queue>


inline Eigen::Vector3d MinVector3d(Eigen::Vector3d v1, Eigen::Vector3d v2) {
    return Eigen::Vector3d(std::min(v1(0), v2(0)),
                           std::min(v1(1), v2(1)),
                           std::min(v1(2), v2(2)));
}

inline Eigen::Vector3d MaxVector3d(Eigen::Vector3d v1, Eigen::Vector3d v2) {
    return Eigen::Vector3d(std::max(v1(0), v2(0)),
                           std::max(v1(1), v2(1)),
                           std::max(v1(2), v2(2)));
}

OneRingHEdge::OneRingHEdge(const Vertex* v) {
    if (v == NULL) start = next = NULL;
    else start = next = v->HalfEdge();
}

HEdge* OneRingHEdge::NextHEdge() {
    HEdge* ret = next;
    if (next && next->Prev()->Twin() != start)
        next = next->Prev()->Twin();
    else
        next = NULL;
    return ret;
}

Mesh::~Mesh() {
    Clear();
}

const HEdgeList& Mesh::Edges() const {
    return heList;
}

const HEdgeList& Mesh::BoundaryEdges() const {
    return bheList;
}

const VertexList& Mesh::Vertices() const {
    return vList;
}

const FaceList& Mesh::Faces() const {
    return fList;
}

// load a .obj mesh definition file
bool Mesh::LoadObjFile(const char* filename) {
    if (filename == NULL || strlen(filename) == 0) return false;
    std::ifstream ifs(filename);
    if (ifs.fail()) return false;

    Clear();

    std::string line;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;

        std::istringstream iss(line);
        std::string type;
        iss >> type;
        // vertex
        if (type.compare("v") == 0) {
            double x, y, z;
            iss >> x >> y >> z;
            AddVertex(new Vertex(x, y, z));
        }
        // face
        else if (type.compare("f") == 0) {
            int index[3];
            iss >> index[0] >> index[1] >> index[2];
            AddFace(index[0] - 1, index[1] - 1, index[2] - 1);
        }
    }
    ifs.close();

    size_t i;
    Eigen::Vector3d box = this->MaxCoord() - this->MinCoord();
    for (i = 0; i < vList.size(); i++) vList[i]->SetPosition(vList[i]->Position() / box(0));

    Eigen::Vector3d tot = Eigen::Vector3d::Zero();
    for (i = 0; i < vList.size(); i++) tot += vList[i]->Position();
    Eigen::Vector3d avg = tot / vList.size();
    for (i = 0; i < vList.size(); i++) vList[i]->SetPosition(vList[i]->Position() - avg);

    HEdgeList list;
    for (i = 0; i < bheList.size(); i++)
        if (bheList[i]->Start()) list.push_back(bheList[i]);
    bheList = list;

    for (i = 0; i < vList.size(); i++) {
        vList[i]->adjHEdges.clear();
        vList[i]->SetIndex((int)i);
        vList[i]->SetFlag(0);
    }

    return true;
}

void Mesh::AddVertex(Vertex* v) {
    vList.push_back(v);
}

void Mesh::AddFace(int v1, int v2, int v3) {
    int i;
    HEdge *he[3], *bhe[3];
    Vertex* v[3];
    Face* f;

    // obtain objects
    for (i = 0; i < 3; i++) he[i] = new HEdge();
    for (i = 0; i < 3; i++) bhe[i] = new HEdge(true);
    v[0] = vList[v1];
    v[1] = vList[v2];
    v[2] = vList[v3];
    f = new Face();

    // connect prev-next pointers
    SetPrevNext(he[0], he[1]);
    SetPrevNext(he[1], he[2]);
    SetPrevNext(he[2], he[0]);
    SetPrevNext(bhe[0], bhe[1]);
    SetPrevNext(bhe[1], bhe[2]);
    SetPrevNext(bhe[2], bhe[0]);

    // connect twin pointers
    SetTwin(he[0], bhe[0]);
    SetTwin(he[1], bhe[2]);
    SetTwin(he[2], bhe[1]);

    // connect start pointers for bhe
    bhe[0]->SetStart(v[1]);
    bhe[1]->SetStart(v[0]);
    bhe[2]->SetStart(v[2]);
    for (i = 0; i < 3; i++) he[i]->SetStart(v[i]);

    // connect start pointers
    // connect face-hedge pointers
    for (i = 0; i < 3; i++) {
        v[i]->SetHalfEdge(he[i]);
        v[i]->adjHEdges.push_back(he[i]);
        SetFace(f, he[i]);
    }
    v[0]->adjHEdges.push_back(bhe[1]);
    v[1]->adjHEdges.push_back(bhe[0]);
    v[2]->adjHEdges.push_back(bhe[2]);

    // mearge boundary if in need
    for (i = 0; i < 3; i++) {
        Vertex* start = bhe[i]->Start();
        Vertex* end = bhe[i]->End();
        for (size_t j = 0; j < end->adjHEdges.size(); j++) {
            HEdge* curr = end->adjHEdges[j];
            if (curr->IsBoundary() && curr->End() == start) {
                SetPrevNext(bhe[i]->Prev(), curr->Next());
                SetPrevNext(curr->Prev(), bhe[i]->Next());
                SetTwin(bhe[i]->Twin(), curr->Twin());
                bhe[i]->SetStart(NULL); // mark as unused
                curr->SetStart(NULL); // mark as unused
                break;
            }
        }
    }

    // finally add hedges and faces to list
    for (i = 0; i < 3; i++) heList.push_back(he[i]);
    for (i = 0; i < 3; i++) bheList.push_back(bhe[i]);
    fList.push_back(f);
}

void Mesh::Clear() {
    size_t i;
    for (i = 0; i < heList.size(); i++) delete heList[i];
    for (i = 0; i < bheList.size(); i++) delete bheList[i];
    for (i = 0; i < vList.size(); i++) delete vList[i];
    for (i = 0; i < fList.size(); i++) delete fList[i];
    for (i = 0; i < BoundaryLoop.size(); i++) delete BoundaryLoop[i];
    heList.clear();
    bheList.clear();
    vList.clear();
    fList.clear();
    BoundaryLoop.clear();
}

Eigen::Vector3d Mesh::MinCoord() const {
    Eigen::Vector3d minCoord = Eigen::Vector3d::Zero();
    for (size_t i = 0; i < vList.size(); i++)
        minCoord = MinVector3d((vList[i])->Position(), minCoord);
    return minCoord;
}

Eigen::Vector3d Mesh::MaxCoord() const {
    Eigen::Vector3d maxCoord = Eigen::Vector3d::Zero();
    for (size_t i = 0; i < vList.size(); i++)
        maxCoord = MaxVector3d((vList[i])->Position(), maxCoord);
    return maxCoord;
}

void Mesh::DisplayMeshInfo() {

	int Num_vertices = vList.size();
	int Num_hedges = heList.size();
	int Num_bhedges = bheList.size();
	int Num_edges;
	int Num_fcaes = fList.size();

	int boundaries;
	int component;
	int genus;
	int Euler;

	Num_edges = (Num_hedges + Num_bhedges) / 2; 
	boundaries = CountBoundaryLoops();
	component = CountConnectedComponents();
    //欧拉公式
	Euler = Num_vertices - Num_edges + Num_fcaes; 
	genus = (2 - Euler - boundaries) / 2; 

	std::cout << "Number of Vertices: " << Num_vertices << std::endl;
	std::cout << "Number of Edges: " << Num_edges << std::endl;
	std::cout << "Number of Faces: " << Num_fcaes << std::endl;
	std::cout << "Number of Boundaries: " << boundaries << std::endl;
	std::cout << "Number of Genus: " << genus << std::endl;
	std::cout << "Number of Components: " << component << std::endl;
	std::cout << "Euler characteristic: " << Euler << std::endl;

    std::cout << "*********Keyboard Function description:********** \n"
        << "1. u(U): explicit Umbrellasmooth \n"
        << "2. s(S): implicit Umbrellasmooth \n"
        << "3  k(K): Skeleton Extraction \n"
        << "4. 1(1): Viewing mode \n"
        << "5. 2(2): Selection mode \n"
        << "6.  ESC: exit" << std::endl;

}

// compute the normal of each vertex
void Mesh::ComputeVertexNormals() {

	const double PI = 3.14159265;
	Eigen::Vector3d t1, t2, nor, Ver1;

	int i, k, j;
	for (i = 0; i < vList.size(); i++) {

		//OneRingVertex ring(vList[i]);  // I didn't use the one ring function. it's a little bit confused especially for new comer of C++;

		k = vList[i]->Valence(); //get the degree of this vertex;
		HEdge* nextHedge;
		Vertex* nextVertex;

		t1 = Eigen::Vector3d(0, 0, 0); //define two vectors;
		t2 = Eigen::Vector3d(0, 0, 0);

		if (vList[i]->IsBoundary()) {

			nextHedge = vList[i]->HalfEdge()->Twin();  //get the first half edge and its twin, one ring vertex;
			nextVertex = nextHedge->Start();

			//make sure the starting half edge is boundary edge, so the next vertex is also on the boundary;
			//while (!(nextVertex->IsBoundary())){
			//nextHedge = nextHedge->Next()->Twin();
			//nextVertex = nextHedge->Start();
			//}

			Ver1 = nextVertex->Position(); //now, the starting vertex is denoted as Ver1;

			if (k == 2) {
				nextVertex = nextHedge->Next()->Twin()->Start();
				t1 = Ver1 - nextVertex->Position();
				t2 = Ver1 + nextVertex->Position() - vList[i]->Position();
			}
			else if (k == 3) {
				nextHedge = nextHedge->Next()->Twin();
				nextVertex = nextHedge->Start();
				t2 = nextVertex->Position() - vList[i]->Position();
				nextHedge = nextHedge->Next()->Twin();
				nextVertex = nextHedge->Start();
				t1 = Ver1 - nextVertex->Position();
			}
			else {
				for (j = 1; j < k - 1; j++) {
					nextHedge = nextHedge->Next()->Twin();
					nextVertex = nextHedge->Start();
					t2[0] += (2 * cos(PI / (k - 1)) - 2)* sin(j*PI / (k - 1))* nextVertex->Position()[0];
					t2[1] += (2 * cos(PI / (k - 1)) - 2)* sin(j*PI / (k - 1))* nextVertex->Position()[1];
					t2[2] += (2 * cos(PI / (k - 1)) - 2)* sin(j*PI / (k - 1))* nextVertex->Position()[2];

				}
				nextHedge = nextHedge->Next()->Twin();
				nextVertex = nextHedge->Start();
				t1 = Ver1 - nextVertex->Position();
				t2 = t2 + sin(PI / (k - 1))* t1;
			}

		}
		else {
			nextHedge = vList[i]->HalfEdge()->Twin();  //get the first half edge that points to this vertex;
			nextVertex = nextHedge->Start();

			for (j = 0; j < k; j++) {

				t1 += cos(2 * PI*j / k)*(nextVertex->Position() - vList[i]->Position());
				t2 += sin(2 * PI*j / k)*(nextVertex->Position() - vList[i]->Position());

				nextHedge = nextHedge->Next()->Twin(); //get next half edge that points to this vertex;
				nextVertex = nextHedge->Start();

			}

		}

		//set the normalied product of two vertors as vertex normal
		nor = -t1.cross(t2); // why it nagates itself??
		nor.normalize();
		vList[i]->SetNormal(nor);
	}


}

// compute the vertex curvature of the graph
void Mesh::ComputeVertexCurvatures() {
    const double PI = 3.14159265;
    Eigen::Vector3d p1,p2,p3 = Eigen::Vector3d(0, 0, 0);
    Eigen::Vector3d Vtemp = Eigen::Vector3d(0, 0, 0); //to calculate the sum of weighted (Vj-Vi)
       
    std::vector<double> Curvature;
    double MeanCur, GaussCur, MaxCur, MinCur; 
    double Cur_max=0; // use for curvature normalization
    double cotAlpha, cotBeta,Areas, dot_res, v_angle = 0; // use for curvature calculation
    double Cur_count[1000] = { 0 };
    double Cur_Cum[1000] = { 0 };//column diagram & Cumulative histogram
    
    for (int i = 0; i < vList.size(); i++) {
        HEdgeList HList ;
        Vtemp = Eigen::Vector3d(0, 0, 0);
        Areas = 0;
        if (vList[i]->IsBoundary()) {
            MeanCur =0;
            GaussCur = 0;
        } 
        else {
            // calculate mean curvature
            OneRingHEdge ring(vList[i]);
            HEdge* curr = NULL;
            //迭代计算，通过while（curr = ring.NextHEdge())可以遍历以vList[i]为起点的所有边
            while (curr = ring.NextHEdge()) {
                HList.push_back(curr);
            }
            //calculate cotangent weighted vector and the areas of triangles
            for (int j=0; j < HList.size(); j++) {
                // cotangent weighted vector
                cotAlpha = Cot(HList[j]->Start()->Position(), 
                               HList[j]->Twin()->Prev()->Start()->Position(), 
                               HList[j]->Twin()->Start()->Position());
                cotBeta = Cot(HList[j]->Start()->Position(), 
                              HList[j]->Next()->End()->Position(), 
                              HList[j]->End()->Position());
                Vtemp += (cotAlpha+cotBeta)*(HList[j]->End()->Position() - HList[j]->Start()->Position());
                // areas of triangles
                Areas += TriArea(HList[j]->Start()->Position(),
                                HList[j]->Twin()->Prev()->Start()->Position(),
                                HList[j]->Twin()->Start()->Position());
            }
            MeanCur = sqrt(Vtemp.dot(Vtemp)) /(2*Areas);
            //calculate gauss curvature
            v_angle = 0;
            for (int j = 0; j < HList.size(); j++) {
                p1 = HList[j]->Start()->Position();
                p2 = HList[j]->Next()->End()->Position();
                p3 = HList[j]->End()->Position();
                dot_res = (p1 - p2).normalized().dot((p1 - p3).normalized());
                dot_res = dot_res < -1 ? -1 : dot_res;
                dot_res = dot_res > 1 ? 1 : dot_res;
                v_angle += std::acos(dot_res); // total angle
            }          
            GaussCur = 2 * (2* PI - v_angle) / Areas;
        }
        //calculate max curvature and min curvature
        MaxCur = MeanCur + sqrt(sqrt(MeanCur * MeanCur - GaussCur));
        MinCur = MeanCur - sqrt(sqrt(MeanCur * MeanCur - GaussCur));

        // set the Color of verteics
        // choose different types of curvature to show on the model
        char chooseCurvature ='H';
        switch (chooseCurvature){
            case 'H'://Mean Curvature
                Curvature.push_back(MeanCur);
                vList[i]->SetColor(Eigen::Vector3d(MeanCur,0, 0));
                Cur_max = MeanCur > Cur_max ? MeanCur : Cur_max;
                break;
            case 'G'://Gauss Curvature                
                GaussCur = abs(GaussCur);
                //std::cout << GaussCur << std::endl;
                Curvature.push_back(GaussCur);
                vList[i]->SetColor(Eigen::Vector3d(GaussCur, 0, 0));
                Cur_max = GaussCur > Cur_max ? GaussCur : Cur_max;
                break;
            case 'I':// Minimum curvature
                Curvature.push_back(MinCur);
                vList[i]->SetColor(Eigen::Vector3d(MinCur, 0, 0));
                Cur_max = MinCur > Cur_max ? MinCur : Cur_max;
                break;
            case 'A'://Maximum curvature
                Curvature.push_back(MaxCur);
                vList[i]->SetColor(Eigen::Vector3d(MaxCur, 0, 0));
                Cur_max = MaxCur > Cur_max ? MaxCur : Cur_max;
                break;
            default:
                vList[i]->SetColor(Eigen::Vector3d(1, 1, 1));
                break;
        }
    }
    //std::cout << "*************************" << std::endl;
    // now we record the wanted curvature in the vector Curvature then we need to normalize the curvature
    //and do the histogram equalization to enhance the contrast of curvature 
    for (int i = 0; i < Curvature.size(); i++) {
        Curvature[i] = Curvature[i] / Cur_max;
        if (int(Curvature[i] * 1000) == 1000) Cur_count[1000 - 1] += 1;
        else  Cur_count[int(Curvature[i] * 1000)] +=1;
        if (i == Curvature.size() - 1) {
            for (int j = 0; j < 1000; j++) {
                Cur_Cum[j] = (j==0)? Cur_count[j]:Cur_count[j] + Cur_Cum[j - 1];
            }
            for (int j = 0; j < Curvature.size(); j++) {
                //int num = Curvature[j] * Curvature.size() ;
                int num = Curvature[j] * 1000;
                Curvature[j] = Cur_Cum[num] / Curvature.size();
                //std::cout <<  Curvature[j] << std::endl;
            }
        }
    }
    //std::cout << "*************************" << std::endl;
    //set color
    for (int i = 0; i < vList.size(); i++) {
        double cur_temp = Curvature[i];
        if (cur_temp < 0.5) {
            vList[i]->SetColor(Eigen::Vector3d(0, 1, 1-cur_temp));//(0,1,1)->(0,1,0.5)
        }
        else if (cur_temp <0.7) {
            vList[i]->SetColor(Eigen::Vector3d((cur_temp-0.5)*5, 1, 0.5-(cur_temp-0.5)*2.5));//(0,1,0.5)->(1,1,0)
            //vList[i]->SetColor(Eigen::Vector3d((cur_temp - 0.5) * 5, 1- (cur_temp - 0.5) * 5, 0.5));//(0,1,0.5)->(1,0,0.5)
        }
        else if (cur_temp < 0.8) {
            vList[i]->SetColor(Eigen::Vector3d(1, 1-(cur_temp-0.7)*5, 0));//(1,1,0)->(1,0.5,0)
            //vList[i]->SetColor(Eigen::Vector3d(1, (cur_temp - 0.7) * 5, 0.5- (cur_temp - 0.7) * 5));//(1,0,0.5)->(1,0.5,0)
        }
        else if (cur_temp < 0.9) {
            vList[i]->SetColor(Eigen::Vector3d(1,0.5-(cur_temp -0.8)*5 ,0));//(1,0.5,0)->(1,0,0)
        }
        else vList[i]->SetColor(Eigen::Vector3d(1- (cur_temp - 0.9) *5 ,0,(cur_temp-0.9)*10 ));  //(1,0,0)->(0.5,0,1)
    }
}

// umbrella smoothing
// uniformWeights: true for uniform-weight Laplacian, false for cotangent-weight Laplacian

void Mesh::UmbrellaSmooth(bool uniformWeights) {
    /*************************/
    /**insert your code here**/
    /*************************/

}

// implicit umbrella smoothing
// uniformWeights: true for uniform-weight Laplacian, false for cotangent-weight Laplacian
void Mesh::ImplicitUmbrellaSmooth(bool uniformWeights) {
    /*************************/
    /* insert your code here */
    /*************************/		

}

//Skeleton Extraction
//by Mesh Contraction
void Mesh::SkeletonParamInit()
{
    VerticesAreaCompute(&originalArea);
    areaRatio = Eigen::VectorXd::Ones(vList.size());
    for (int i = 0; i < originalArea.size(); i++) {
        m_ml += originalArea[i];
    }
    m_ml =150* sqrt(m_ml /vList.size());
    m_mh = 3;
    m_ml = 10;
}

bool Mesh::Skeleton()
{
    /*************************/
   /***date：2020-11-27******/
   /*editor:SUN CHANGJIANG***/
   /*************************/

   //in this hw we need to calculate the Laplacian Martex first
   //and I use cotangent weight there

   //compute the L matrix using cotangent weight
    SparseMatrixBuilder* Lapl = new SparseMatrixBuilder();
    for (int i = 0; i < vList.size(); i++) {

        //compute the summation of weight around a vertex;
        int k = vList[i]->Valence();
        HEdge* nextHedge = vList[i]->HalfEdge()->Twin();
        Eigen::Vector3d preVertex, currVertex, nexVertex;
        double cotAphla, cotBeta, sumWeight = 0.0;

        // compute weight for each one-ring vertex
        for (int r = 0; r < k; r++) {
            currVertex = nextHedge->Start()->Position();
            preVertex = nextHedge->Twin()->Prev()->Start()->Position();
            nexVertex = nextHedge->Next()->Twin()->Start()->Position();
            cotAphla = Mesh::Cot(vList[i]->Position(), nexVertex, currVertex);
            cotBeta = Mesh::Cot(vList[i]->Position(), preVertex, currVertex);
            Lapl->AddEntry(i, nextHedge->Start()->Index(), m_ml * (cotAphla + cotBeta));
            sumWeight += cotAphla + cotBeta;
            nextHedge = nextHedge->Next()->Twin();
        }
        Lapl->AddEntry(i, i, -sumWeight * m_ml); //the value Lii = -1 because the weight of other vertex have divided sumWeight        
    }
    m_ml = m_ml * SL;

    // add the matrix WH in the end of Lapl
    VerticesAreaCompute(&currentArea);
    for (int i = 0; i < vList.size(); i++) {
        areaRatio[i] = m_mh * sqrt(originalArea[i] / currentArea[i]);
        //std::cout << "current areaRatio of the vertex" << i << "is :" << areaRatio[i] << std::endl;
        Lapl->AddEntry(vList.size() + i, i, areaRatio[i]);
    }

    Eigen::SparseMatrix<double> A(2 * vList.size(), vList.size());
    Eigen::SparseMatrix<double> M(vList.size(), vList.size());

    A = Lapl->ToSparseMatrix(2 * vList.size(), vList.size());
    M = A.transpose() * A;
    solver = new SparseLinearSystemSolver(M);

    // structure the Column vector b
    //get delta vectors of x, y, z coord for each vertex before transformation
    vx = Eigen::VectorXd::Zero(vList.size());
    vy = Eigen::VectorXd::Zero(vList.size());
    vz = Eigen::VectorXd::Zero(vList.size());

    for (int i = 0; i < vList.size(); i++) {
        vx[i] = vList[i]->Position()(0);
        vy[i] = vList[i]->Position()(1);
        vz[i] = vList[i]->Position()(2);
    }
    Eigen::VectorXd bx = Eigen::VectorXd::Zero(2 * vList.size());
    Eigen::VectorXd by = Eigen::VectorXd::Zero(2 * vList.size());
    Eigen::VectorXd bz = Eigen::VectorXd::Zero(2 * vList.size());
    for (int i = 0; i < vList.size(); i++) {
        bx[vList.size() + i] = areaRatio[i] * vx[i];
        by[vList.size() + i] = areaRatio[i] * vy[i];
        bz[vList.size() + i] = areaRatio[i] * vz[i];
    }

    Eigen::VectorXd bx2 = Eigen::VectorXd::Zero(vList.size());
    Eigen::VectorXd by2 = Eigen::VectorXd::Zero(vList.size());
    Eigen::VectorXd bz2 = Eigen::VectorXd::Zero(vList.size());
    bx2 = A.transpose() * bx;
    by2 = A.transpose() * by;
    bz2 = A.transpose() * bz;
    bx = solver->Solve(bx2);
    by = solver->Solve(by2);
    bz = solver->Solve(bz2);

    // mesh Contraction
    for (int i = 0; i < vList.size(); i++) {
        Eigen::Vector3d newPosition = Eigen::Vector3d(bx[i], by[i], bz[i]);
        vList[i]->SetPosition(newPosition);
    }
    //
    meshcontractioncount++;
    std::cout << "Mesh contraction step: " << meshcontractioncount << " finished" << std::endl;

    if (currentArea.back() / originalArea.back() <= 0.1) {
        std::cout << "Mesh contract to the threshold(0.1) " << std::endl;
        return true;
    }
    else
        return false;
}
//use the bhelist directly
int Mesh::CountBoundaryLoopsV1() {
    
    /*
    start from bhelist[0] as boundary1, iterate through this vector and find the next boundary half edge 
    which the start vertex is the end vertex of boundary1 and call it boundary2. then repeat this step, until we 
    find a boundary half edge which the end vertex is the start vertex of boundary1.so we get a boundary loop successfully.
    then we can remember this loop, delect it from bhelist, and find another loop from bhelist. finally, we can get 
    all the boundary loops. then we can calculate the amount of loops  and draw it in the window. 
    */
    //int BoundaryLoopNum = 0;
    HEdgeList BdEdges = bheList;
    HEdgeList* Boundaryloop = new HEdgeList;
    HEdge* bhetemp = BdEdges[0];
    int bhe_hand = 0;
    while (!BdEdges.empty()) {
        for (size_t i = 0; i < BdEdges.size(); i++) {
            if (bhetemp->End() == BdEdges[i]->Start()) {
                Boundaryloop->push_back(bhetemp); // add boundary half edge in the boundaryloop
                bhetemp = BdEdges[i];                                            
                BdEdges.erase(BdEdges.begin() + bhe_hand);
                if (bhe_hand < i) bhe_hand = i - 1;
                else bhe_hand = i;
                break;
            }
            else if (i == (BdEdges.size()-1)) {
                Boundaryloop->push_back(bhetemp);
                BoundaryLoop.push_back(Boundaryloop);
                Boundaryloop = new HEdgeList;
                BdEdges.erase(BdEdges.begin() + bhe_hand);
                bhetemp = BdEdges[0];
                bhe_hand = 0;           
            }
        }
    }  
	return BoundaryLoop.size();
}
//get the bhelist by BFS
int Mesh::CountBoundaryLoops()
{
    HEdgeList Bhelist;
    Vertex* v1;

    std::queue<Vertex*> VertexQueue;
    VertexQueue.push(vList[0]);
    vList[0]->SetValid(false);
    //loop through the vList by BFS, and get the boundary half edge list named Bhelist
    while (!VertexQueue.empty()) {
        v1 = VertexQueue.front();
        OneRingHEdge ring(v1);
        HEdge* curr = NULL;
        //迭代计算，通过while（curr = ring.NextHEdge())可以遍历以v1为起点的所有边
        while (curr = ring.NextHEdge()) {
            if (curr->IsBoundary()) {
                Bhelist.push_back(curr);
            }
            if (curr->End()->IsValid()) {
                VertexQueue.push(curr->End());
                curr->End()->SetValid(false);
            }
        }
        VertexQueue.pop();
    }
    /************calculate the num of boundaryloop ************/
    if (Bhelist.empty()) return 0;
    HEdgeList BdEdges = Bhelist;
    HEdgeList* Boundaryloop = new HEdgeList;
    HEdge* bhetemp = BdEdges[0];
    int bhe_hand = 0;
    while (!BdEdges.empty()) {
        for (size_t i = 0; i < BdEdges.size(); i++) {
            if (bhetemp->End() == BdEdges[i]->Start()) {
                Boundaryloop->push_back(bhetemp); // add boundary half edge in the boundaryloop
                bhetemp = BdEdges[i];
                BdEdges.erase(BdEdges.begin() + bhe_hand);
                if (bhe_hand < i) bhe_hand = i - 1;
                else bhe_hand = i;
                break;
            }
            else if (i == (BdEdges.size() - 1)) {
                Boundaryloop->push_back(bhetemp);
                BoundaryLoop.push_back(Boundaryloop);
                Boundaryloop = new HEdgeList;
                BdEdges.erase(BdEdges.begin() + bhe_hand);
                bhetemp = BdEdges[0];
                bhe_hand = 0;
            }
        }
    }
    return BoundaryLoop.size();
}

void Mesh::TraverseComponents(HEdge* he)
{
	return;
}


int Mesh::CountConnectedComponents() 
{
    /*************************/
    /* insert your code here */
    /*************************/
	return 0;
}

void Mesh::GroupingVertexFlags() {
    // set vertex flag to be 255 initially
    for (size_t i = 0; i < vList.size(); i++)
        if (vList[i]->Flag() != 0)
            vList[i]->SetFlag(255);

    int id = 0;
    VertexList tmpList;
    for (int i = 0; i < vList.size(); i++)
        if (vList[i]->Flag() == 255) {
            id++;
            vList[i]->SetFlag(id);
            tmpList.push_back(vList[i]);
            while (! tmpList.empty()) {
                Vertex* v = tmpList.back();
                tmpList.pop_back();
                OneRingVertex ring = OneRingVertex(v);
                while (Vertex* v2 = ring.NextVertex()) {
                    if (v2->Flag() == 255) {
                        v2->SetFlag(id);
                        tmpList.push_back(v2);
                    }
                }
            }
        }
}

void Mesh::VerticesAreaCompute(std::vector<double> *vec)
{
    double totalAreas = 0;
    for (int i = 0; i < vList.size(); i++) {
        HEdgeList HList;
        OneRingHEdge ring(vList[i]);
        HEdge* curr = NULL;
        
        double Areas = 0;
        //迭代计算，通过while（curr = ring.NextHEdge())可以遍历以vList[i]为起点的所有边
        while (curr = ring.NextHEdge()) {
            HList.push_back(curr);
        }
        //calculate cotangent weighted vector and the areas of triangles
        for (int j = 0; j < HList.size(); j++) {
            // areas of triangles
            Areas += TriArea(HList[j]->Start()->Position(),
                HList[j]->Twin()->Prev()->Start()->Position(),
                HList[j]->Twin()->Start()->Position());
        }
        totalAreas += Areas / 3;
        vec->push_back(Areas/3);
        //std::cout << "current area of the vertex" << i << "is :" << Areas << std::endl;
    }
    std::cout << "the total areas is: " << totalAreas << std::endl;
    vec->push_back(totalAreas);
}

void Mesh::SetPrevNext(HEdge* e1, HEdge* e2) {
    e1->SetNext(e2);
    e2->SetPrev(e1);
}

void Mesh::SetTwin(HEdge* e1, HEdge* e2) {
    e1->SetTwin(e2);
    e2->SetTwin(e1);
}

void Mesh::SetFace(Face* f, HEdge* e) {
    f->SetHalfEdge(e);
    e->SetFace(f);
}

double Mesh::Cot(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) {
    Eigen::Vector3d v1 = p1 - p2;
    Eigen::Vector3d v2 = p3 - p2;

    double _dot_res = v1.normalized().dot(v2.normalized());
    if (_dot_res < -1.0) {
        _dot_res = -1.0;
    }
    else if (_dot_res >  1.0) {
        _dot_res = 1.0;
    }
    return 1.0 / std::tan(std::acos(_dot_res));
}

double Mesh::TriArea(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) {
    Eigen::Vector3d v1 = p2 - p1;
    Eigen::Vector3d v2 = p3 - p1;
    return v1.cross(v2).norm() / 2.0;
}

void Vertex::SetParameterPosition(float x, float y)
{
	parameterposition[0] = x;
	parameterposition[1] = y;
}
