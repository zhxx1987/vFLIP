#ifndef _coarse_fine_solver
#define _coarse_fine_solver
#include "Smoke_solver3D.h"


using namespace  std;
namespace ivock{

	class CoarseFineSmokeSolver
	{
	public:
		CoarseFineSmokeSolver();
		~CoarseFineSmokeSolver();
		void Finalize();
		uint _nx,_ny,_nz,_N;
		float _h;
		SmokeSolver3D fine_solver;
		SmokeSolver3D coarse_solver;
		void initialize(uint nx, uint ny, uint nz, float L);

		void FineProjection();
		void SetFaceBoundary(buffer3Df &fine, buffer3Df &coarse, int offx, int offy, int offz);
		void FineToCoarseOperator(buffer3Df &fine, buffer3Df &coarse, int offx, int offy, int offz);
		void FineToCoarse();
		void CoarseToFineOperator(buffer3Df &fine, buffer3Df &coarse);
		void ComputeSmallDifference(buffer3Df &fine, buffer3Df &coarse, buffer3Df &diff);
		void LargeScaleAdvection(float dt);
		void SmallScaleAdvection(float dt);
		void TimeStep(float dt);
		Vec3f getCoarseVel(Vec3f &pos);
		Vec3f getFineVel(Vec3f &pos);
		Vec3f traceCoarse(float dt);
		Vec3f traceFine(float dt);
	};

}



#endif