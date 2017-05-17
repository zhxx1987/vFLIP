#ifndef __buffer_h__
#define __buffer_h__

#include "array.h"
#include <vector>
#include "vec.h"
using namespace std;
namespace ivock{

	template<class T>
	class Buffer3D{
	public:
		Buffer3D(){ _nx=_ny=_nz=_n=0; _hx=_hy=_hz=0.0; _data = new array1D<T>; }
		~Buffer3D(){}
		array1D<T> *_data;
		vector<vector<vector<T*>>> index;
		uint _blockN;
		Vec3f origin;


		// f should be between 0 and 1
		inline void cubic_interp_weights(T f, T& wneg1, T& w0, T& w1, T& w2)
		{
			T f2(f*f), f3(f2*f);
			wneg1=-T(1./3)*f+T(1./2)*f2-T(1./6)*f3;
			w0=1-f2+T(1./2)*(f3-f);
			w1=f+T(1./2)*(f2-f3);
			w2=T(1./6)*(f3-f);
		}

		template<class S>
		inline S cubic_interp(const S& value_neg1, const S& value0, 
			const S& value1, const S& value2, T f)
		{
			T wneg1, w0, w1, w2;
			cubic_interp_weights(f, wneg1, w0, w1, w2);
			return wneg1*value_neg1 + w0*value0 + w1*value1 + w2*value2;
		}




		uint _nx, _ny, _nz, _n;
		double _hx, _hy, _hz;
		double _ox, _oy, _oz;
		uint _physical_zstride, _physical_ystride, _physical_nx, _physical_ny, _physical_nz, _physical_n,_blockx,_blocky,_blockz;
		void setOrigin(Vec3f & c)
		{
			origin = c;
		}
		void init(Buffer3D<T> &f)
		{
			origin = f.origin;
			_nx=f._nx;_ny=f._ny;_nz=f._nz;
			_n=_nx*_ny*_nz;
			_hx=_hy=_hz=f._hx;
			_ox=f._ox; _oy=f._oy; _oz=f._oz;
			_blockN = 8;

			_blockx = (!(_nx%_blockN))?_nx/_blockN:(_nx/_blockN + 1);
			_blocky = (!(_ny%_blockN))?_ny/_blockN:(_ny/_blockN + 1);
			_blockz = (!(_nz%_blockN))?_nz/_blockN:(_nz/_blockN + 1);
			_physical_nx = _blockx*_blockN;
			_physical_ny = _blocky*_blockN;
			_physical_nz = _blockz*_blockN;
			_physical_n  = _physical_nx*_physical_ny*_physical_nz;
			_data->alloc(_physical_n);
			_data->setZero();


			//index.clear();
			//index.resize(_blockz);
			//for (int k=0;k<_blockz;k++)
			//{
			//	index[k].resize(_blocky);
			//}
			//for (int k=0;k<_blockz;k++)for(int j=0;j<_blocky;j++)
			//{
			//	index[k][j].resize(_blockx);
			//}
			//for (int k=0;k<_blockz;k++)for(int j=0;j<_blocky;j++) for (int i=0;i<_blockx;i++)
			//{
			//	index[k][j][i] = &(_data->getPtr()[((k*_blocky+j)*_blockx+i)*_blockN*_blockN*_blockN]);
			//}
		}
		void init(uint nx, uint ny, uint nz, double h, double ox, double oy, double oz)
		{
			origin = Vec3f(0);
			_nx=nx;_ny=ny;_nz=nz;
			_n=nx*ny*nz;
			_hx=_hy=_hz=h;
			_ox=ox; _oy=oy; _oz=oz;
			_blockN = 8;

			//data are now going to be stored by blocks
			_blockx = (!(_nx%_blockN))?_nx/_blockN:(_nx/_blockN + 1);
			_blocky = (!(_ny%_blockN))?_ny/_blockN:(_ny/_blockN + 1);
			_blockz = (!(_nz%_blockN))?_nz/_blockN:(_nz/_blockN + 1);
			_physical_nx = _blockx*_blockN;
			_physical_ny = _blocky*_blockN;
			_physical_nz = _blockz*_blockN;
			_physical_n  = _physical_nx*_physical_ny*_physical_nz;
			_data->alloc(_physical_n);
			_data->setZero();


			//index.clear();
			//index.resize(_blockz);
			//for (int k=0;k<_blockz;k++)
			//{
			//	index[k].resize(_blocky);
			//}
			//for (int k=0;k<_blockz;k++)for(int j=0;j<_blocky;j++)
			//{
			//	index[k][j].resize(_blockx);
			//}
			//for (int k=0;k<_blockz;k++)for(int j=0;j<_blocky;j++) for (int i=0;i<_blockx;i++)
			//{
			//	index[k][j][i] = &(_data->getPtr()[((k*_blocky+j)*_blockx+i)*_blockN*_blockN*_blockN]);
			//}



		}

		void init(uint nx, uint ny, uint nz)
		{
			origin = Vec3f(0);
			_nx=nx;_ny=ny;_nz=nz;
			_n=nx*ny*nz;
			_hx=_hy=_hz=0;
			_ox=0; _oy=0; _oz=0;
			_blockN = 8;


			_blockx = (!(_nx%_blockN))?_nx/_blockN:(_nx/_blockN + 1);
			_blocky = (!(_ny%_blockN))?_ny/_blockN:(_ny/_blockN + 1);
			_blockz = (!(_nz%_blockN))?_nz/_blockN:(_nz/_blockN + 1);
			_physical_nx = _blockx*_blockN;
			_physical_ny = _blocky*_blockN;
			_physical_nz = _blockz*_blockN;
			_physical_n  = _physical_nx*_physical_ny*_physical_nz;

			_data->alloc(_physical_n);
			_data->setZero();

			//index.clear();
			//index.resize(_blockz);
			//for (int k=0;k<_blockz;k++)
			//{
			//	index[k].resize(_blocky);
			//}
			//for (int k=0;k<_blockz;k++)for(int j=0;j<_blocky;j++)
			//{
			//	index[k][j].resize(_blockx);
			//}
			//for (int k=0;k<_blockz;k++)for(int j=0;j<_blocky;j++) for (int i=0;i<_blockx;i++)
			//{
			//	index[k][j][i] = &(_data->getPtr()[((k*_blocky+j)*_blockx+i)*_blockN*_blockN*_blockN]);
			//}



		}

		array1D<T> *getArray(){return _data;}
		uint getSize() {return _n;}
		void setZero() {_data->setZero();}
		void free() {
			_data->free();
			//for (int k=0;k<_blockz;k++)for(int j=0;j<_blocky;j++)
			//{
			//	index[k][j].resize(0);
			//}
			//for (int k=0;k<_blockz;k++)
			//{
			//	index[k].resize(0);
			//}
			//index.resize(0); 
		}
		void copy(Buffer3D<T> &b)
		{
			memcpy(_data->getPtr(), (b.getArray())->getPtr(),sizeof(T)*_physical_n);

		}
		const T& operator()(uint i, uint j, uint k) const
		{
			uint I = i>>3, J=j>>3, K=k>>3;
			uint ii = i&7, jj = j&7, kk = k&7;
			uint idx = ((K*_blockx*_blocky + J*_blockx +I)<<9) + (kk<<6) + (jj<<3) +ii;
			return (_data->getPtr()[idx]);
			//return (index[K][J][I])[(((kk<<3)+jj)<<3)+ii];
		}
		T& operator()(uint i, uint j, uint k)
		{

			uint I = i>>3, J=j>>3, K=k>>3;
			uint ii = i&7, jj = j&7, kk = k&7;
			uint idx = ((K*_blockx*_blocky + J*_blockx +I)<<9) + (kk<<6) + (jj<<3) +ii;
			return (_data->getPtr()[idx]);
			//return (index[K][J][I])[(((kk<<3)+jj)<<3)+ii];
		}
		//indexing

		T at(int i, int j, int k)
		{
			uint ti=std::min(std::max(i,(int)0),(int)_nx-1), tj=std::min(std::max(j,(int)0),(int)_ny-1),tk=std::min(std::max(k,(int)0),(int)_nz-1);
			//if(i>=0 && i<_nx && j>=0 && j<_ny && k>=0 && k<_nz )
			{
				uint I = ti>>3, J=tj>>3, K=tk>>3;
				uint ii = ti&7, jj = tj&7, kk = tk&7;
				uint idx = ((K*_blockx*_blocky + J*_blockx +I)<<9) + (kk<<6) + (jj<<3) +ii;
				return (_data->getPtr()[idx]);
				//return (index[K][J][I])[(((kk<<3)+jj)<<3)+ii];
			}
			//else
			//{
			//	return (T)0;
			//}
		}
		inline double lerp(double a, double b, double c)
		{
			return (1.0-c)*(double)a + c*(double)b;
		}
		//sampling
		T sample_linear(double world_x, double world_y, double world_z)
		{
			double grid_x = (world_x-origin[0])/_hx + _ox;
			double grid_y = (world_y-origin[1])/_hy + _oy;
			double grid_z = (world_z-origin[2])/_hz + _oz;

			int grid_i = (int)floor(grid_x);
			int grid_j = (int)floor(grid_y);
			int grid_k = (int)floor(grid_z);

			double cx = grid_x - (double)grid_i;
			double cy = grid_y - (double)grid_j;
			double cz = grid_z - (double)grid_k;

			double v1 = lerp(lerp((double)at(grid_i,grid_j,grid_k), (double)at(grid_i+1,grid_j,grid_k),cx),
				lerp((double)at(grid_i,grid_j+1,grid_k), (double)at(grid_i+1, grid_j+1, grid_k),cx),
				cy);
			double v2 = lerp(lerp((double)at(grid_i,grid_j,grid_k+1), (double)at(grid_i+1,grid_j,grid_k+1),cx),
				lerp((double)at(grid_i,grid_j+1,grid_k+1), (double)at(grid_i+1, grid_j+1, grid_k+1),cx),
				cy);
			return (T)(lerp(v1, v2, cz));

		}
		Vec3f sample_grad(double world_x, double world_y, double world_z)
		{
			T right = sample_linear(world_x+_hx,world_y,world_z);
			T left  = sample_linear(world_x-_hx,world_y,world_z);
			T up    = sample_linear(world_x,world_y+_hx,world_z);
			T bottom= sample_linear(world_x,world_y-_hx,world_z);
			T front = sample_linear(world_x,world_y,world_z-_hx);
			T back  = sample_linear(world_x,world_y,world_z+_hx);
			return 0.5f*Vec3f((right-left)/_hx, (up-bottom)/_hx, (back-front)/_hx);
		}
		T sample_cubic(double world_x, double world_y, double world_z)
		{
			double gx =(world_x-origin[0])/_hx + _ox;
			double gy =(world_y-origin[1])/_hy + _oy;
			double gz =(world_z-origin[2])/_hz + _oz;

			int i = (int)floor(gx);
			int j = (int)floor(gy);
			int k = (int)floor(gz);

			double x = gx - (double)i;
			double y = gy - (double)j;
			double z = gz - (double)k;

			double fn=cubic_interp((double)(at(i-1,j-1,k-1)), (double)(at(i,j-1,k-1)), 
				(double)(at(i+1,j-1,k-1)), (double)(at(i+2,j-1,k-1)), x);
			double f0=cubic_interp((double)(at(i-1,j,k-1)), (double)(at(i,j,k-1)), 
				(double)(at(i+1,j,k-1)), (double)(at(i+2,j,k-1)), x);
			double f1=cubic_interp((double)(at(i-1,j+1,k-1)), (double)(at(i,j+1,k-1)), 
				(double)(at(i+1,j+1,k-1)), (double)(at(i+2,j+1,k-1)), x);
			double f2=cubic_interp((double)(at(i-1,j+2,k-1)), (double)(at(i,j+2,k-1)), 
				(double)(at(i+1,j+2,k-1)), (double)(at(i+2,j+2,k-1)), x);

			double fzn = cubic_interp(fn, f0, f1, f2, y);

			fn=cubic_interp((double)(at(i-1,j-1,k)), (double)(at(i,j-1,k)), 
				(double)(at(i+1,j-1,k)), (double)(at(i+2,j-1,k)), x);
			f0=cubic_interp((double)(at(i-1,j,k)), (double)(at(i,j,k)), 
				(double)(at(i+1,j,k)), (double)(at(i+2,j,k)), x);
			f1=cubic_interp((double)(at(i-1,j+1,k)), (double)(at(i,j+1,k)), 
				(double)(at(i+1,j+1,k)), (double)(at(i+2,j+1,k)), x);
			f2=cubic_interp((double)(at(i-1,j+2,k)), (double)(at(i,j+2,k)), 
				(double)(at(i+1,j+2,k)), (double)(at(i+2,j+2,k)), x);

			double fz0 = cubic_interp(fn, f0, f1, f2, y);


			fn=cubic_interp((double)(at(i-1,j-1,k+1)), (double)(at(i,j-1,k+1)), 
				(double)(at(i+1,j-1,k+1)), (double)(at(i+2,j-1,k+1)), x);
			f0=cubic_interp((double)(at(i-1,j,k+1)), (double)(at(i,j,k+1)), 
				(double)(at(i+1,j,k+1)), (double)(at(i+2,j,k+1)), x);
			f1=cubic_interp((double)(at(i-1,j+1,k+1)), (double)(at(i,j+1,k+1)), 
				(double)(at(i+1,j+1,k+1)), (double)(at(i+2,j+1,k+1)), x);
			f2=cubic_interp((double)(at(i-1,j+2,k+1)), (double)(at(i,j+2,k+1)), 
				(double)(at(i+1,j+2,k+1)), (double)(at(i+2,j+2,k+1)), x);

			double fz1 = cubic_interp(fn, f0, f1, f2, y);


			fn=cubic_interp((double)(at(i-1,j-1,k+2)), (double)(at(i,j-1,k+2)), 
				(double)(at(i+1,j-1,k+2)), (double)(at(i+2,j-1,k+2)), x);
			f0=cubic_interp((double)(at(i-1,j,k+2)), (double)(at(i,j,k+2)), 
				(double)(at(i+1,j,k+2)), (double)(at(i+2,j,k+2)), x);
			f1=cubic_interp((double)(at(i-1,j+1,k+2)), (double)(at(i,j+1,k+2)), 
				(double)(at(i+1,j+1,k+2)), (double)(at(i+2,j+1,k+2)), x);
			f2=cubic_interp((double)(at(i-1,j+2,k+2)), (double)(at(i,j+2,k+2)), 
				(double)(at(i+1,j+2,k+2)), (double)(at(i+2,j+2,k+2)), x);

			double fz2 = cubic_interp(fn, f0, f1, f2, y);


			double res = cubic_interp(fzn, fz0, fz1, fz2, z);


			
			return res;

		}


		T sample_cube_lerp(double world_x, double world_y, double world_z, float &v0, float &v1,float &v2,float &v3,float &v4,float &v5,float &v6,float &v7)
		{
			double grid_x = (world_x-origin[0])/_hx + _ox;
			double grid_y = (world_y-origin[1])/_hy + _oy;
			double grid_z = (world_z-origin[2])/_hz + _oz;

			int grid_i = (int)floor(grid_x);
			int grid_j = (int)floor(grid_y);
			int grid_k = (int)floor(grid_z);

			double cx = grid_x - (double)grid_i;
			double cy = grid_y - (double)grid_j;
			double cz = grid_z - (double)grid_k;

			v0 = (double)at(grid_i,grid_j,grid_k); v1 = (double)at(grid_i+1,grid_j,grid_k);
			v2 = (double)at(grid_i,grid_j+1,grid_k); v3 = (double)at(grid_i+1, grid_j+1, grid_k);
			v4 = (double)at(grid_i,grid_j,grid_k+1); v5 = (double)at(grid_i+1,grid_j,grid_k+1);
			v6 = (double)at(grid_i,grid_j+1,grid_k+1); v7 = (double)at(grid_i+1, grid_j+1, grid_k+1);

			double iv1 = lerp(lerp((double)at(grid_i,grid_j,grid_k), (double)at(grid_i+1,grid_j,grid_k),cx),
				lerp((double)at(grid_i,grid_j+1,grid_k), (double)at(grid_i+1, grid_j+1, grid_k),cx),
				cy);
			double iv2 = lerp(lerp((double)at(grid_i,grid_j,grid_k+1), (double)at(grid_i+1,grid_j,grid_k+1),cx),
				lerp((double)at(grid_i,grid_j+1,grid_k+1), (double)at(grid_i+1, grid_j+1, grid_k+1),cx),
				cy);
			return (T)(lerp(iv1, iv2, cz));

		}
		void sample_cube(double world_x, double world_y, double world_z, float &v0, float &v1,float &v2,float &v3,float &v4,float &v5,float &v6,float &v7)
		{
			double grid_x = (world_x-origin[0])/_hx + _ox;
			double grid_y = (world_y-origin[1])/_hy + _oy;
			double grid_z = (world_z-origin[2])/_hz + _oz;

			int grid_i = (int)floor(grid_x);
			int grid_j = (int)floor(grid_y);
			int grid_k = (int)floor(grid_z);

			double cx = grid_x - (double)grid_i;
			double cy = grid_y - (double)grid_j;
			double cz = grid_z - (double)grid_k;

			v0 = (double)at(grid_i,grid_j,grid_k); v1 = (double)at(grid_i+1,grid_j,grid_k);
			v2 = (double)at(grid_i,grid_j+1,grid_k); v3 = (double)at(grid_i+1, grid_j+1, grid_k);
			v4 = (double)at(grid_i,grid_j,grid_k+1); v5 = (double)at(grid_i+1,grid_j,grid_k+1);
			v6 = (double)at(grid_i,grid_j+1,grid_k+1); v7 = (double)at(grid_i+1, grid_j+1, grid_k+1);

		}
	};

	typedef Buffer3D<float> buffer3Df;
	typedef Buffer3D<double> buffer3Dd;
	typedef Buffer3D<char> buffer3Dc;
	typedef Buffer3D<int> buffer3Di;
	typedef Buffer3D<uint> buffer3Dui;
}

#endif