#include <iostream>
#include <TPZGenGrid3D.h>
#include <TPZVTKGeoMesh.h>
#include <TPZHDivApproxCreator.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <TPZLinearAnalysis.h>
#include "TPZSSpStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSimpleTimer.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZVTKGenerator.h"
#include "TPZStructMatrixOMPorTBB.h"
#include <fstream>
#include <string>
#include <vector>

using namespace std;



enum EMatid {EDomain, EBC};
#define USE_STD
//#define USE_TBB
//#define USE_OMP

int main (int argc, char * const argv[]) {

fstream textfile;
ofstream resultfile;
int threadnum =0;
textfile.open("../PTests.txt",ios::in);
if (textfile.is_open()){ //checking whether the file is open
  string tp;
  while(getline(textfile, tp)){ //read data from file object and put it into string.
     threadnum =  stoi(tp);
  }
  textfile.close(); //close the file object.
  }  
  const bool isSolveAndPostProc = false;

  // Create geometric mesh
  const TPZVec<REAL> minX = {-1.,-1.,-1.};
  const TPZVec<REAL> maxX = {1.,1.,1.};
  const int ndivperside = 15;
  const TPZVec<int> nelDiv = {ndivperside,ndivperside,ndivperside};
  const MMeshType elType = MMeshType::EHexahedral;
  TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
  TPZGeoMesh *gmesh = gen3d.BuildVolumetricElements(EDomain);
  gmesh = gen3d.BuildBoundaryElements(EBC, EBC, EBC, EBC, EBC, EBC);
  
  const bool isPrintGmesh = true;
  if (isPrintGmesh){
    ofstream out("geomesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
  }
  
  // Create computational mesh
  const int pOrder = 3;
  TPZHDivApproxCreator hdivCreator(gmesh);
  hdivCreator.SetDefaultOrder(pOrder);
  hdivCreator.ProbType() = ProblemType::EDarcy;
  
  // Insert volumetric material
  TPZMixedDarcyFlow* matdarcy = new TPZMixedDarcyFlow(EDomain,gmesh->Dimension());
  matdarcy->SetConstantPermeability(1.);
  hdivCreator.InsertMaterialObject(matdarcy);
  
  // Insert boundary conditions
  const REAL appliedPressure = 1.;
  const int dirType = 0;
  TPZFMatrix<STATE> val1(1,1,0.);
  TPZManVector<STATE> val2(1,appliedPressure);
  TPZBndCondT<STATE> *BCond1 = matdarcy->CreateBC(matdarcy, EBC, dirType, val1, val2);
  hdivCreator.InsertMaterialObject(BCond1);

  // Create computational
  TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace();
  
  // Create analysis
  const int nThreads = threadnum;
  TPZLinearAnalysis an(cmesh,true);
#if defined(USE_STD)
  TPZSSpStructMatrix<STATE,TPZStructMatrixOR<STATE>> matsp(cmesh);
#elif defined(USE_TBB)
  TPZSSpStructMatrix<STATE,TPZStructMatrixOMPorTBB<STATE>> matsp(cmesh);
  matsp.SetTBBorOMP(true);
#elif defined(USE_OMP)
  TPZSSpStructMatrix<STATE,TPZStructMatrixOMPorTBB<STATE>> matsp(cmesh);
  matsp.SetTBBorOMP(false);
#else
  DebugStop();
#endif
  matsp.SetNumThreads(nThreads);
  an.SetStructuralMatrix(matsp);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELDLt);
  an.SetSolver(step);

  // Setup timer and assemble
  cout << "\n--------- Starting assemble -----------" << endl;
  TPZSimpleTimer time_assemble("Assemble timer");
  an.Assemble();
  cout << "Total time for assemble = " << time_assemble.ReturnTimeDouble() / 1000. << " seconds" << endl;
  resultfile.open("../Results.txt",ios_base::app);
  resultfile <<  "Total time for assemble = " << time_assemble.ReturnTimeDouble() / 1000. << " seconds"  << endl;
  
  // Solve problem (Not needed for timing!)
  if(isSolveAndPostProc) an.Solve();
  
  // Post process if needed
  if(isSolveAndPostProc){
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), cmesh);
    const std::string plotfile = "PostProcess"; //sem o .vtk no final
    constexpr int vtkRes{0};
    
    TPZManVector<std::string,2> fields = {"Flux","Pressure"};
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    
    vtk.Do();
  }
  
  return 0;
}
