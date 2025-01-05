/*!
 @file cuda_utils.h

 @brief CUDA用の関数・定義など
  - cppファイルからもインクルードするので，CUDAカーネル関数やデバイス(__device__)関数はここには書かない
  - 共通のデバイス関数はcuda_utils.cuの方に書く

 @author Makoto Fujisawa
 @date 2023-02
*/


#ifndef _RX_CUDA_UTILS_H_
#define _RX_CUDA_UTILS_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

#include <helper_cuda.h>
#include <helper_cuda_gl.h>
//#include <helper_math.h>

#include "vector_types.h"
#include "vector_functions.h"


//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
typedef unsigned int uint;
typedef unsigned char uchar;


// 1ブロックあたりのスレッド数
static int THREAD_NUM = 256;

// 各種定数
#define PI_F      3.14159265359f
#define DEG2RAD_F 0.0174533f
#define DIM       3
#define QUAT      4           //海老沢追加 四元数のサイズ
#define VEL_EPSILON 1.0e-4    //海老沢追加 粒子を静止させるかどうかを判定
#define ANGVEL_EPSILON 1.0e-6 //海老沢追加 回転速度を静止させるかどうかを判定
#define MU 0.3                //海老沢追加 静止摩擦係数 0.3

const float DEG_TO_RAD = 0.0174532925199432;	//! degree -> radian の変換係数(pi/180.0)
const float RAD_TO_DEG = 57.295779513082320;	//! radian -> degree の変換係数(180.0/pi)


#define MAX_BOX_NUM 10
struct CuBox
{
	float3 cen, ext;
	float3 rot[3], inv_rot[3];
	int flg;
};

#define MAX_SPHERE_NUM 10
struct CuSphere
{
	float3 cen;
	float rad;
	unsigned int flg;
};



//-----------------------------------------------------------------------------
// 近傍探索用データ
//-----------------------------------------------------------------------------
struct Cell
{
	uint3 GridSize;
	float3 WorldOrigin;
	float3 WorldMax;
	float3 CellWidth;

	float3* dSortedPos;			//!< ソート済みパーティクル座標
	float3* dSortedVel;			//!< ソート済みパーティクル速度

	uint* dSortedIndex;			//!< ソート済みパーティクルインデックス
	uint* dGridParticleHash;	//!< 各パーティクルのグリッドハッシュ値(ソート用キー)
	uint* dCellStart;			//!< ソートリスト内の各セルのスタートインデックス
	uint* dCellEnd;				//!< ソートリスト内の各セルのエンドインデックス
	uint  uNumParticles;		//!< 総パーティクル数
	uint  uNumCells;			//!< 総セル数
};


//-----------------------------------------------------------------------------
// 各種パラメータ格納用構造体
//  - GPU-CPU間でパラメータをまとめてやり取りするために使用
//-----------------------------------------------------------------------------
struct SceneParameter 
{
	// ユーザ入力パラメータ
	int max_particles;			//!< 最大粒子数

	float3 boundary_cen;		//!< シミュレーション空間の中心
	float3 boundary_ext;		//!< シミュレーション空間の大きさ(各辺の長さの1/2)
	float3 boundary_min;		//!< シミュレーション空間の最小座標
	float3 boundary_max;		//!< シミュレーション空間の最大座標
	float dens;					//!< 初期密度
	float mass;					//!< 粒子の質量
	uint kernel_particles;		//!< 有効半径h以内の粒子数

	float dt;					//!< 時間ステップ幅
	int step;					//!< 現ステップ数
	float t;					//!< 現時刻(=step*dt)

	float viscosity;			//!< 動粘性係数(or XSPH人工粘性用の係数)
	float vorticity;			//!< vorticity confinement用係数
	float gas_k;				//!< ガス定数
	float gamma;				//!< Tait方程式の係数([Becker2007]ではγ=7)
	float B;					//!< Tait方程式の係数(B=ρ0*cs^2/gamma, [Becker2007]ではcs≈88.5[m/s], B≈1119[kPa])

	float res;					//!< 境界での反発係数
	float3 gravity;				//!< 重力加速度(ベクトル)
	int use_inlet;				//!< 流入境界条件の有無

	int mesh_n;					//!< 表面メッシュ生成:セル分割数(最大)
	float mesh_thr;				//!< 表面メッシュ生成:閾値(表面を表す関数値)

	// ユーザ入力パラメータから計算されるパラメータ
	float effective_radius;		//!< カーネル関数W計算時に用いる有効半径h
	float h;					//!< 有効半径
	float kernel_radius;		//!< 近傍探索半径(Wの範囲が0～2hという場合に対応するために必要)
	float rest_dens;			//!< 有効半径h決定後に改めて計算された初期密度(実際にはこちらを用いる)
	float particle_radius;		//!< 粒子半径(主に描画用)

	float aw;					//!< カーネルの定数係数
	float ag;					//!< 勾配計算用カーネルの定数係数
	float al;					//!< ラプラシアン計算用カーネルの定数係数

	// 一定速度場，渦速度場による移動用
	float3 d;					//!< 1方向へ移動するときの移動ベクトル
	float vscl;					//!< 渦速度場のスケール
	float tmax;					//!< 渦速度場の時間スケール

	Cell cell;					//!< 近傍探索用グリッドセル情報

#if MAX_BOX_NUM
	CuBox box[MAX_BOX_NUM];
#endif
	uint num_box;

#if MAX_SPHERE_NUM
	CuSphere sphere[MAX_SPHERE_NUM];
#endif
	uint num_sphere;

	// 視点
	int view;					//!< 視点設定がある場合は1, 四元数による視点設定があるときは2
	float3 view_trans;
	float3 view_rot;
	float4 view_quat;
	float3 bgcolor;
};
inline void InitSceneParameter(SceneParameter &p)
{
	// シミュレーションパラメータの初期化
	p.max_particles = 50000;
	p.boundary_cen = make_float3(0.0, 0.0, 0.0);
	p.boundary_ext = make_float3(1.0, 1.0, 1.0);
	p.boundary_min = make_float3(-1.0, -1.0, -1.0);
	p.boundary_max = make_float3( 1.0,  1.0,  1.0);
	p.dens = p.rest_dens = (float)998.29;
	p.mass = (float)0.04;
	p.kernel_particles = (float)20.0;
	p.dt = 0.001;
	p.viscosity = 1.0e-3;
	p.vorticity = 1.0e-3;

	p.gas_k = 3.0;
	p.gamma = 7.0;
	p.B = 1000.0;

	p.res = 0.0;
	p.gravity = make_float3(0.0, -9.82, 0.0);
	p.use_inlet = 0;

	p.mesh_n = 64;
	p.mesh_thr = 100.0f;

	p.t = 0.0f;
	p.step = 0;

	p.view = 0;
	p.view_trans = make_float3(0.0, 0.0, -5.0);
	p.view_rot   = make_float3(0.0, 0.0, 0.0);
	p.view_quat  = make_float4(1.0, 0.0, 0.0, 0.0);
	p.bgcolor    = make_float3(1.0, 1.0, 1.0);

	p.cell.GridSize = make_uint3(1, 1, 1);
	p.cell.WorldOrigin = make_float3(0, 0, 0);
	p.cell.WorldMax = make_float3(1, 1, 1);
	p.cell.CellWidth = make_float3(1, 1, 1);

	p.cell.dSortedPos = 0;			// ソート済みパーティクル座標
	p.cell.dSortedVel = 0;			// ソート済みパーティクル速度

	p.cell.dSortedIndex = 0;		// ソート済みパーティクルインデックス
	p.cell.dGridParticleHash = 0;	// 各パーティクルのグリッドハッシュ値(ソート用キー)
	p.cell.dCellStart = 0;			// ソートリスト内の各セルのスタートインデックス
	p.cell.dCellEnd = 0;			// ソートリスト内の各セルのエンドインデックス
	p.cell.uNumParticles = 0;		// 総パーティクル数
	p.cell.uNumCells = 1;			// 総セル数

	p.num_box = 0;
	p.num_sphere = 0;
}


//-----------------------------------------------------------------------------
// CUDA変数用出力オペレータの定義
//-----------------------------------------------------------------------------
//! float3の出力オペレータ
inline std::ostream& operator<<(std::ostream& out, const float3& a)
{
	return out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
}

//! uint3の出力オペレータ
inline std::ostream& operator<<(std::ostream& out, const uint3& a)
{
	return out << "(" << a.x << ", " << a.y << ", " << a.z << ")";
}


//-----------------------------------------------------------------------------
// CUDAエラーチェック用関数
//-----------------------------------------------------------------------------
#define CUCHECK(val) CudaCheckFunc(val, #val, __FILE__, __LINE__)

template<typename T>
inline bool CudaCheckFunc(T rtrn_val, const char* func, const char* file, const int line)
{
	if(rtrn_val){
		fprintf(stderr, "CUDA error at %s line %d : %s (error code = %d), function=%s\n",
			file, line, cudaGetErrorString(rtrn_val), (int)(rtrn_val), func);
		return true;
	} else
	{
		return false;
	}
}

#define CUERROR(msg) CudaLastError(msg, __FILE__, __LINE__)

inline bool CudaLastError(const char* msg, const char* file, const int line)
{
	cudaError_t err = cudaGetLastError();
	if(err != cudaSuccess){
		fprintf(stderr, "CUDA last error at %s line %d : %s (error code = %d)\n",
			file, line, cudaGetErrorString(err), (int)(err));
		return true;
	} else
	{
		return false;
	}
}


//-----------------------------------------------------------------------------
// 行列
//  - 大きさ3のfloat3配列で3x3行列を表現
//-----------------------------------------------------------------------------
//! 3x3行列の初期化(単位行列)
inline void Identity(float3 m[3])
{
	m[0] = make_float3(1.0f, 0.0f, 0.0f);
	m[1] = make_float3(0.0f, 1.0f, 0.0f);
	m[2] = make_float3(0.0f, 0.0f, 1.0f);
}

//! matrix3x3の逆行列
inline void Inverse(float3 m[3], float3 inv_m[3])
{
	float d = m[0].x*m[1].y*m[2].z- 
			  m[0].x*m[2].y*m[1].z+ 
			  m[1].x*m[2].y*m[0].z- 
			  m[1].x*m[0].y*m[2].z+ 
			  m[2].x*m[0].y*m[1].z- 
			  m[2].x*m[1].y*m[0].z;

	if(d == 0) d = 1;

	inv_m[0] = make_float3( (m[1].y*m[2].z-m[1].z*m[2].y)/d, -(m[0].y*m[2].z-m[0].z*m[2].y)/d, (m[0].y*m[1].z-m[0].z*m[1].y)/d);
	inv_m[1] = make_float3(-(m[1].x*m[2].z-m[1].z*m[2].x)/d,  (m[0].x*m[2].z-m[0].z*m[2].x)/d, -(m[0].x*m[1].z-m[0].z*m[1].x)/d);
	inv_m[2] = make_float3( (m[1].x*m[2].y-m[1].y*m[2].x)/d, -(m[0].x*m[2].y-m[0].y*m[2].x)/d, (m[0].x*m[1].y-m[0].y*m[1].x)/d);
}


//! オイラー角から回転行列を生成
inline void EulerToMatrix(float3 m[3], float pitch, float yaw, float roll)
{
	yaw   = DEG_TO_RAD*(yaw);
	pitch = DEG_TO_RAD*(pitch);
	roll  = DEG_TO_RAD*(roll);

	float cy = cos(yaw); 
	float sy = sin(yaw); 
	float cp = cos(pitch); 
	float sp = sin(pitch); 
	float cr = cos(roll);
	float sr = sin(roll);

	float cc = cy*cr; 
	float cs = cy*sr; 
	float sc = sy*cr; 
	float ss = sy*sr;

	m[0] = make_float3(cc+sp*ss, cs-sp*sc, -sy*cp);
	m[1] = make_float3(-cp*sr,   cp*cr,    -sp);
	m[2] = make_float3(sc-sp*cs, ss+sp*cc, cy*cp);
}



#endif // #ifndef _RX_CUDA_UTILS_H_