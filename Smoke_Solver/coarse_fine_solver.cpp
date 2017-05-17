#include "coarse_fine_solver.h"




void ivock::CoarseFineSmokeSolver::TimeStep( float dt )
{
	//given u_fine, u_coarse
	//fine_solver.advect(rho,T);
	//U~ = coarse_solver.advect(U);
	//compute u* = u - U~;
	//stretch-advect curl(u*) with U;
	//apply small scale vorticity
	//FineToCoarse(rho, T, u)
	//CoarseSolver: apply force, projection
	//FineSolver:Flux cell boundary condition, apply f*, projection
}

void ivock::CoarseFineSmokeSolver::SmallScaleAdvection( float dt )
{
	ComputeSmallDifference(fine_solver._un,coarse_solver._un,fine_solver._un);
	ComputeSmallDifference(fine_solver._vn,coarse_solver._vn,fine_solver._vn);
	ComputeSmallDifference(fine_solver._wn,coarse_solver._wn,fine_solver._wn);


	fine_solver.compute_curl();
	CoarseToFineOperator(fine_solver._un, coarse_solver._un);
	CoarseToFineOperator(fine_solver._vn, coarse_solver._vn);
	CoarseToFineOperator(fine_solver._wn, coarse_solver._wn);

	fine_solver.stretch(dt);
	fine_solver.advect_field(dt,fine_solver._wxn,fine_solver._wxnp1);
	fine_solver.advect_field(dt,fine_solver._wyn,fine_solver._wynp1);
	fine_solver.advect_field(dt,fine_solver._wzn,fine_solver._wznp1);
}

void ivock::CoarseFineSmokeSolver::LargeScaleAdvection( float dt )
{
	coarse_solver.advect_field(dt,coarse_solver._un,coarse_solver._unp1);
	coarse_solver.advect_field(dt,coarse_solver._vn,coarse_solver._vnp1);
	coarse_solver.advect_field(dt,coarse_solver._wn,coarse_solver._wnp1);
	CoarseToFineOperator(fine_solver._unp1, coarse_solver._unp1);
	CoarseToFineOperator(fine_solver._vnp1, coarse_solver._vnp1);
	CoarseToFineOperator(fine_solver._wnp1, coarse_solver._wnp1);
}

void ivock::CoarseFineSmokeSolver::ComputeSmallDifference(buffer3Df &fine, buffer3Df &coarse, buffer3Df &diff)
{
	diff.setZero();
	int slice = fine._blockx * fine._blocky;
	tbb::parallel_for((size_t)0,
		(size_t)(fine._blockx*fine._blocky*fine._blockz),
		(size_t)1,
		[&](size_t thread_idx)
	{
		int bk = thread_idx/(slice);
		int bj = (thread_idx%slice)/fine._blockx;
		int bi = thread_idx%(fine._blockx);
		for (int kk=0;kk<8;kk++)for(int jj=0;jj<8;jj++)for(int ii=0;ii<8;ii++)
		{
			int i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;

			if(i>=0&&i<fine._nx
				&&j>=0&&j<fine._ny
				&&k>=0&&k<fine._nz)
			{
				Vec3f fpos = (Vec3f(i,j,k) - Vec3f(fine._ox,fine._oy,fine._oz))*_h;
				diff(i,j,k) = fine(i,j,k) - coarse.sample_linear(fpos[0],fpos[1],fpos[2]);
			}
		}
	});
}
void ivock::CoarseFineSmokeSolver::FineToCoarseOperator(buffer3Df &fine, buffer3Df &coarse, int offx, int offy, int offz)
{
	int slice = coarse._blockx * coarse._blocky;
	tbb::parallel_for((size_t)0,
		(size_t)(coarse._blockx*coarse._blocky*coarse._blockz),
		(size_t)1,
		[&](size_t thread_idx)
	{
		int bk = thread_idx/(slice);
		int bj = (thread_idx%slice)/coarse._blockx;
		int bi = thread_idx%(coarse._blockx);
		for (int kk=0;kk<8;kk++)for(int jj=0;jj<8;jj++)for(int ii=0;ii<8;ii++)
		{
			int i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;

			if(i>=0&&i<coarse._nx
				&&j>=0&&j<coarse._ny
				&&k>=0&&k<coarse._nz)
			{
				int count = 0;
				float q = 0;
				for(int iii = max(0,i*_N - offx); iii <min(fine._nx,i*_N + _N - offx);iii++)
				{
					for(int jjj = max(0, j*_N - offy); jjj<min(fine._ny, j*_N + _N - offy);jjj++)
					{
						for (int kkk = max(0, k*_N - offz); kkk<min(fine._nz, k*_N + _N - offz);kkk++)
						{
							count++;
							q += fine(iii,jjj,kkk);
						}
					}
				}
				coarse(i,j,k) = (count>0)?q/(float)count:0;
			}
		}
	});
}
void ivock::CoarseFineSmokeSolver::FineToCoarse()
{
	FineToCoarseOperator(fine_solver._un, coarse_solver._un,_N/2,0,0);
	FineToCoarseOperator(fine_solver._vn, coarse_solver._vn,0,_N/2,0);
	FineToCoarseOperator(fine_solver._wn, coarse_solver._wn,0,0,_N/2);
	FineToCoarseOperator(fine_solver._rho,coarse_solver._rho,0,0,0);
	FineToCoarseOperator(fine_solver._Tbf,coarse_solver._Tbf,0,0,0);

}

void ivock::CoarseFineSmokeSolver::initialize( uint nx, uint ny, uint nz, float L )
{
	_nx = nx;
	_ny = ny;
	_nz = nz;
	_h  = L/(float)_nx;
	coarse_solver.init(nx/_N,ny/_N,nz/_N,L);
	coarse_solver.init(nx, ny, nz, L);
}

void ivock::CoarseFineSmokeSolver::Finalize()
{
	fine_solver.Finalize();
	coarse_solver.Finalize();
}

ivock::CoarseFineSmokeSolver::~CoarseFineSmokeSolver()
{
	Finalize();
}

ivock::CoarseFineSmokeSolver::CoarseFineSmokeSolver()
{
	_N=8;
}

void ivock::CoarseFineSmokeSolver::CoarseToFineOperator( buffer3Df &fine, buffer3Df &coarse)
{
	int slice = fine._blockx * fine._blocky;
	tbb::parallel_for((size_t)0,
		(size_t)(fine._blockx*fine._blocky*fine._blockz),
		(size_t)1,
		[&](size_t thread_idx)
	{
		int bk = thread_idx/(slice);
		int bj = (thread_idx%slice)/fine._blockx;
		int bi = thread_idx%(fine._blockx);
		for (int kk=0;kk<8;kk++)for(int jj=0;jj<8;jj++)for(int ii=0;ii<8;ii++)
		{
			int i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;

			if(i>=0&&i<fine._nx
				&&j>=0&&j<fine._ny
				&&k>=0&&k<fine._nz)
			{
				Vec3f fpos = (Vec3f(i,j,k) - Vec3f(fine._ox,fine._oy,fine._oz))*_h;
				fine(i,j,k) = coarse.sample_linear(fpos[0],fpos[1],fpos[2]);
			}
		}
	});
}

void ivock::CoarseFineSmokeSolver::FineProjection()
{
	
	//for each coarse cell face velocity
	//get its corresponding velocity face in fine grid

	//fine cell set boundary velocity

}

void ivock::CoarseFineSmokeSolver::SetFaceBoundary( buffer3Df &fine, buffer3Df &coarse, int offx, int offy, int offz )
{
	int slice = coarse._blockx * coarse._blocky;
	tbb::parallel_for((size_t)0,
		(size_t)(coarse._blockx*coarse._blocky*coarse._blockz),
		(size_t)1,
		[&](size_t thread_idx)
	{
		int bk = thread_idx/(slice);
		int bj = (thread_idx%slice)/coarse._blockx;
		int bi = thread_idx%(coarse._blockx);
		for (int kk=0;kk<8;kk++)for(int jj=0;jj<8;jj++)for(int ii=0;ii<8;ii++)
		{
			int i=bi*8+ii, j=bj*8+jj, k=bk*8+kk;

			if(i>=0&&i<coarse._nx
				&&j>=0&&j<coarse._ny
				&&k>=0&&k<coarse._nz)
			{
				
				float q = 0;
				for(int iii = max(0,i*_N); iii <min(fine._nx,i*_N + offx);iii++)
				{
					for(int jjj = max(0, j*_N); jjj<min(fine._ny, j*_N +  offy);jjj++)
					{
						for (int kkk = max(0, k*_N); kkk<min(fine._nz, k*_N + offz);kkk++)
						{
							
							q += fine(iii,jjj,kkk);
						}
					}
				}
				float coarse_flux = coarse(i,j,k);
				float diff_flux = coarse_flux - q/(float)(_N*_N);
				for(int iii = max(0,i*_N); iii <min(fine._nx,i*_N + offx);iii++)
				{
					for(int jjj = max(0, j*_N); jjj<min(fine._ny, j*_N +  offy);jjj++)
					{
						for (int kkk = max(0, k*_N); kkk<min(fine._nz, k*_N + offz);kkk++)
						{

							fine(iii,jjj,kkk) += diff_flux;
						}
					}
				}
			}
		}
	});
}
