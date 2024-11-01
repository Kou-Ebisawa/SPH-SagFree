/*! 
  @file sph.cuh
	
  @brief CUDA関数の宣言
		 - CUDAを呼び出すC++のコードは基本的にこのファイルをインクルードする
 
  @author Makoto Fujisawa
  @date 2023-02
*/

#ifndef _SPH_CUH_
#define _SPH_CUH_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "cuda_utils.h"

//-----------------------------------------------------------------------------
// CUDA関数
//-----------------------------------------------------------------------------
extern "C"
{

//-----------------------------------------------------------------------------
// 粒子処理
void CuSphDensity(float* drestdens,float* ddens, float* dvol, int n);
void CuSphPressure(float* drestdens,float* dpres, float* ddens, int n);

void CuSphVorticity(float* dvort, float* dvel, float* ddens, float* dvol, int* datt, int n);
void CuSphForces(float* drestdens,float* dacc, float* dvel, float* ddens, float* dpres, float* dvort, float* dvol, int* datt,float3 power, int n);
void CuSphViscosityForces(float* drestdens,float* dacc, float* dvel, float* ddens, float* dvol, int* datt, int n);
void CuSphXSPHViscosity(float* dvel, float* ddens, float* dvol, int* datt, int n);
void CuSphIntegrate(float* dpos, float* dvel, float* dacc, int* datt, int* fix,int n);//海老沢変更 fixを追加
void CuSphIntegrateV(float* dvel, float* dacc, int* datt, int n);
void CuSphIntegrateP(float* dpos, float* dvel, int* datt, int* dfix, int n);

//海老沢追加--------------------------------------------------------------------------------
//XPBDの制約
void CuXPBDConstraint(float* dpos, float* dmas, float* dlen, float* dkss, float* dkbt,float* dquat,float* domega, float* dlamb_ss,float* dlamb_bt,int* dfix, float dt, int n,int iter);
//衝突制約
void CuCollisionConstraint(float* dpos, float* dvel, int* dfix, float3 center, float rad, float dt, int n);
//時間積分
void CuIntegrate(float* dpos, float* dcurpos, float* dvel, float dt, int n);
//外力計算
void CuCalExternalForces(float* dpos, float* dmass, float* dvel, int* dfix,float3 gravity, float3 wind, float dt, int n);
//PBDの位置ベース法
void CuPBDStretchingConstraint(float* dpos, float* dmas, float* dlen, float* dkss, float* dquat, int* dfix, int n, int iter);
//現在のある配列の出力(デバック用)
void CuPrint3Dfloat(float* dpos,float* dvel,float* dacc,int n);
//接線の更新
void CuTangUpdate(float* dpos, float* dtang, int* dfix, int n);
//重力のみの計算
void CuOnlyGravity(float* dvel, float* dmass, int* dfix,float dt, int n);
//いくつかのパラメータに0を代入
void CuSetParametersZero(float* dangvel,float* dfss,float*dpbf_lambda, int n);
//角加速度の更新
void CuAngVelUpdate(float* dangvel, float* dquat,int* dfix,float dt, int n);
//各加速度の時間積分
void CuAngVelIntegrate(float* dangvel,float* dcurquat, float* dquat,int* dfix,float dt, int n);
//四元数の設定
void CuQuatSet(float* dcurquat, float* dquat, int* dfix, int n);
//基準となる姿勢の変更
void CuRestDensSet(float* dpos, float* dRestDens, float* dvol, int n);
//一律の基準密度の設定
void CuRestTotalDens(float* drestdens,float dens, int n);
//SagFree処理-------
//グローバルフォースステップ
void CuGlobalForceStep(float* dfss,float* dmass,int* dlast_index,float3 gravity,int num_elastic);
//ローカルフォースステップ
void CuLocalForceStep(float* dpos, float* dlen, float* dquat,float* dcurquat, float* dkss, float* dfss, int* dfix, int n);
//グローバルトルクステップ
void CuGlobalTorqueStep(float* dpos, float* dquat, float* domega, float* dlen, float* dkss, float* dkbt, int* dfix, int* dlast_index, int num_elastic);
//グローバルトルクステップ(Videoを参考にした方)
void CuVideoGlobalTorqueStep(float* dpos, float* dquat, float* domega, float* dlen, float* dkss, float* dkbt, int* dfix, int* dlast_index, int num_elastic);
//ローカルトルクステップ
void CuLocalTorqueStep(float* dquat, float* domega, float* dlen, float* dkbt, int* dfix, int n);
//-----------------
//密度制約
void CuPbfConstraint(float* dpos, float* ddens, float* drestdens, float* dpbf_lambda, float* dvol, int n);
//PBFで解く場合の外力計算
void CuPbfExternalForces(float* dacc, int* datt, float3 power, int n);
//------------------------------------------------------------------------------------------

// 粒子体積計算
void CuSphCalVolume(float* dvol, int *datt, int n, float v);

// 粒子密度場生成(表面メッシュ生成用)
void CuSphDensityInGrid(float* dF, float* dvol, int* datt, int n, int3 gnum, float3 gmin, float3 glen);

// 粒子描画色計算
void CuColorScalar(float* dcol, int* datt, float* dval, int n, float3 c1, float3 c2, float2 range);
void CuColorVector(float* dcol, int* datt, float* dval, int n, float3 c1, float3 c2, float2 range);
void CuColorConstant(float* dcol, int* datt, float3 col, int n);

// 近傍粒子探索用
void CuCalcHash(uint* dhash, uint* dindex, float* dpos, int n);
void CuSort(unsigned int* dhash, uint* dindex, uint n);
void CuReorderDataAndFindCellStart(Cell cell, float* dpos, float* dvel, uint n);

// パラメータをGPUメモリに転送
void CuSetParameters(const SceneParameter* hparams);

// デバッグ用
float CuCalAverage(float* data, int n);		// スカラー値が入った配列の値の平均値を求めて返す
float CuCalAverageV(float* data, int n);	// ベクトル値が入った配列のベクトルの長さの平均値を求めて返す
void CuScan(float* dScanData, float* dData, int num);

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// GPU処理補助
void CuInit();
void CuSetDevice(int id);
void CuDeviceProp(void);

void CuAllocateArray(void** devPtr, size_t size);
void CuSetArrayValue(void* devPtr, int val, size_t size);
void CuFreeArray(void* devPtr);

void CuCopyArrayD2D(void* dDst, void* dSrc, int size);
void CuCopyArrayFromDevice(void* host, void* device, cudaGraphicsResource** resource, int offset, int size);
void CuCopyArrayToDevice(void* device, const void* host, int offset, int size);

void CuThreadSync(void);

void CuRegisterGLBufferObject(unsigned int vbo, cudaGraphicsResource** resource);
void CuUnregisterGLBufferObject(cudaGraphicsResource* resource);
void* CuMapGLBufferObject(cudaGraphicsResource** resource);
void CuUnmapGLBufferObject(cudaGraphicsResource* resource);
//-----------------------------------------------------------------------------

} // extern "C"


#endif // #ifdef _SPH_CUH_