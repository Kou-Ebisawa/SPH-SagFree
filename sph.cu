/*! 
  @file sph.cu
	
  @brief CUDA : SPH法

  @author Makoto Fujisawa
  @date 2023-02
*/


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <cstdio>
#include <GL/glew.h>
#if __APPLE__
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <GL/gl.h>
	#include <GL/glu.h>
#endif

#include "sph_kernel.cu"

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

//-----------------------------------------------------------------------------
// CUDA関数
//-----------------------------------------------------------------------------
extern "C"
{
/*!
 * パラメータをGPUへ転送
 * @param[in] hparams ホスト(CPU)メモリに格納されたパラメータ
 */
void CuSetParameters(const SceneParameter* hparams)
{
	CUCHECK(cudaMemcpyToSymbol(params, hparams, sizeof(SceneParameter), 0, cudaMemcpyHostToDevice));
}

/*!
 * スレッド数からブロック/グリッド数の計算(必要スレッド数n以上になるように設定)
 * @param[in] n 必要スレッド数
 * @param[out] block,grid ブロック数/グリッド数
 */
void CuCalGridN(int n, dim3& block, dim3& grid)
{
	// スレッド数の設定(n以上になるように設定)
	block = dim3(THREAD_NUM, 1, 1); // 1ブロック当りスレッド数
	grid = dim3((n+block.x-1)/block.x, 1, 1); // 1グリッドあたりのブロック数
}

/*!
 * CUDAカーネルを使って並列計算:粒子密度の計算
 *  - 粒子位置はcell.dSortedPosから取得するので引数として渡す必要なし
 * @param[out] ddens 粒子密度(デバイスメモリ)
 * @param[in] dvol 粒子体積(デバイスメモリ)
 * @param[in] n 粒子数
 */
void CuSphDensity(float* drestdens,float* ddens, float* dvol,float* dmas, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxSphDensity<<<grid, block>>>(drestdens,ddens, dvol,dmas, n);	// カーネル実行
	cudaThreadSynchronize();
}

/*!
 * CUDAカーネルを使って並列計算:粒子圧力値を密度から計算
 * @param[out] dpres 粒子圧力(デバイスメモリ)
 * @param[in] ddens 粒子密度(デバイスメモリ)
 * @param[in] n 粒子数
 */
void CuSphPressure(float* drestdens,float* dpres, float* ddens, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxSphPressure<<<grid, block>>>(drestdens,dpres, ddens, n);	// カーネル実行
	cudaThreadSynchronize();
}

/*!
* CUDAカーネルを使って並列計算:各粒子の渦度計算[Macklin2013]
*  - M. Macklin & M. M\"{u}ller, "Position Based Fluids", ACM ToG, 32(4), pp.104:1-104:12, 2013.
*  - グリッド法向けの元の手法は 
*    R. Fedkiw; J. Stam & H. Jensen, "Visual simulation of smoke", Proc. SIGGRAPH 2001, pp.15-22, 2001.
* @param[out] dvort 各粒子の渦度ベクトル(デバイスメモリ)
* @param[in] dvel 粒子速度配列(デバイスメモリ)
* @param[in] ddens 粒子密度(デバイスメモリ)
* @param[in] dvol 粒子体積(デバイスメモリ)
* @param[in] datt 粒子属性(0で流体,1で境界)(デバイスメモリ)
* @param[in] n 粒子数
*/
void CuSphVorticity(float* dvort, float* dvel, float* ddens, float* dvol, int* datt, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxSphVorticity<<<grid, block>>>(dvort,dvel, ddens, dvol, datt, n);	// カーネル実行
	cudaThreadSynchronize();
}

/*!
 * CUDAカーネルを使って並列計算:粒子に働く力(圧力項&外力項)の計算
 * @param[out] dacc 各粒子に働く力(加速度dv/dt)(デバイスメモリ)
 * @param[in] dvel 粒子速度配列(デバイスメモリ)
 * @param[in] ddens 粒子密度(デバイスメモリ)
 * @param[in] dpres 粒子圧力(デバイスメモリ)
 * @param[in] dvort 各粒子の渦度ベクトル(デバイスメモリ)
 * @param[in] dvol  粒子体積(デバイスメモリ)
 * @param[in] datt 粒子属性(0で流体,1で境界)(デバイスメモリ)
 * @param[in] n 粒子数
 */
void CuSphForces(float* drestdens,float* dacc, float* dvel, float* ddens, float* dpres, float* dvort, float* dvol,float* dmas, int* datt,float3 power,float* dfss, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxSphForces<<<grid, block>>>(drestdens,dacc, dvel, ddens, dpres, dvort, dvol,dmas, datt,power,dfss, n);	// カーネル実行
	cudaThreadSynchronize();
}
/*!
* CUDAカーネルを使って並列計算:粒子に働く力(粘性項)の計算[Becker2007]
*  - M. Becker & M. Teschner, "Weakly Compressible SPH for Free Surface Flows", Proc. SCA2007, pp.209-217, 2007.
* @param[inout] dacc 各粒子に働く力(加速度dv/dt)(デバイスメモリ)
* @param[in] dvel 粒子速度配列(デバイスメモリ)
* @param[in] ddens 粒子密度(デバイスメモリ)
* @param[in] dvol 粒子体積(デバイスメモリ)
* @param[in] datt 粒子属性(0で流体粒子，それ以外で境界粒子)(デバイスメモリ)
* @param[in] n 粒子数
*/
void CuSphViscosityForces(float* drestdens,float* dacc, float* dvel, float* ddens, float* dvol,float* dmas, int* datt, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxSphViscosity<<<grid, block>>>(drestdens,dacc, dvel, ddens, dvol,dmas, datt, n);	// カーネル実行
	cudaThreadSynchronize();
}

/*!
* CUDAカーネルを使って並列計算:XSPH人工粘性の計算[Schechter2012]
*  - H. Schechter & R. Bridson, "Ghost SPH for animating water", ACM ToG, 31(4), pp.61:1-61:8, 2012.
*  - 力(加速度)として粘性項を加えるのではなく，速度を直接更新する
*  - 流体の物質的性質としての粘性ではなく計算安定性のための計算といった方がよさそう
* @param[inout] dvel 粒子速度配列(デバイスメモリ)
* @param[in] ddens 粒子密度(デバイスメモリ)
* @param[in] dvol 粒子体積(デバイスメモリ)
* @param[in] datt 粒子属性(0で流体粒子，それ以外で境界粒子)(デバイスメモリ)
* @param[in] n 粒子数
*/
void CuSphXSPHViscosity(float* dvel, float* ddens, float* dvol,float* dmas, int* datt, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxSphXSPHViscosity<<<grid, block>>>(dvel, ddens, dvol,dmas, datt, n);	// カーネル実行
	cudaThreadSynchronize();
}


/*!
 * CUDAカーネルを使って並列計算:加速度に従って位置と速度を更新
 * @param[inout] dpos 粒子位置配列(デバイスメモリ)
 * @param[inout] dvel 粒子速度配列(デバイスメモリ)
 * @param[in] dacc 各粒子に働く力(加速度)を格納した配列(デバイスメモリ)
 * @param[in] dvol 粒子体積(デバイスメモリ)
 * @param[in] datt 粒子属性(0で流体粒子，それ以外で境界粒子)(デバイスメモリ)
 * fix:海老沢追加 固定点を表す
 * @param[in] n 粒子数
 */
void CuSphIntegrate(float* dpos, float* dvel, float* dacc, int* datt, int* dfix,int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxSphIntegrate<<<grid, block>>>(dpos, dvel, dacc, datt,dfix, n);	// カーネル実行
	cudaThreadSynchronize();
}


/*!
* CUDAカーネルを使って並列計算:加速度に従って速度のみを更新
*  - XSPH用
* @param[inout] dvel 粒子速度配列(デバイスメモリ)
* @param[in] dacc 各粒子に働く力(加速度)を格納した配列(デバイスメモリ)
* @param[in] dvol 粒子体積(デバイスメモリ)
* @param[in] datt 粒子属性(0で流体粒子，それ以外で境界粒子)(デバイスメモリ)
* @param[in] n 粒子数
*/
void CuSphIntegrateV(float* dvel, float* dacc, int* datt, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxSphIntegrateVelocity<<<grid, block>>>(dvel, dacc, datt, n);	// カーネル実行
	cudaThreadSynchronize();
}
/*!
* CUDAカーネルを使って並列計算:速度を従って位置を更新
*  - XSPH用
* @param[inout] dpos 粒子位置配列(デバイスメモリ)
* @param[in] dvel 粒子速度配列(デバイスメモリ)
* @param[in] dvol 粒子体積(デバイスメモリ)
* @param[in] datt 粒子属性(0で流体粒子，それ以外で境界粒子)(デバイスメモリ)
* @param[in] n 粒子数
*/
void CuSphIntegrateP(float* dpos, float* dvel, int* datt,int*dfix, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxSphIntegratePosition<<<grid, block>>>(dpos, dvel, datt,dfix, n);	// カーネル実行
	cudaThreadSynchronize();
}



/*!
* CUDAカーネルを使って並列計算:粒子体積の計算
*  - 粒子位置はcell.dSortedPosから取得するので引数として渡す必要なし
* @param[out] dvol 粒子体積(デバイスメモリ)
* @param[in] datt 粒子属性(0で流体粒子，それ以外で境界粒子)(デバイスメモリ)
* @param[in] n 処理する粒子数(offsetからの相対的な位置)
* @param[in] v 流体粒子の場合の粒子体積値
*/
void CuSphCalVolume(float* dvol, int *datt, int n, float v)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxSphCalVolume<<<grid, block>>>(dvol, datt, n, v);	// カーネル実行
	cudaThreadSynchronize();
}

/*!
* CUDAカーネルを使って並列計算:粒子密度の計算
*  - 粒子位置はcell.dSortedPosから取得するので引数として渡す必要なし
* @param[out] dF 密度値を格納するグリッドセル配列(デバイスメモリ)
* @param[in] dvol 粒子体積(デバイスメモリ)
* @param[in] datt 粒子属性(0で流体,1で境界)(デバイスメモリ)
* @param[in] n 粒子数
* @param[in] gnum グリッド数
* @param[in] gmin グリッド最小座標
* @param[in] glen グリッド幅
*/
void CuSphDensityInGrid(float* dF, float* dvol, int* datt, int n, int3 gnum, float3 gmin, float3 glen)
{
	// 総グリッドセル数=スレッド数としてブロック/グリッドサイズを計算
	int numcell = gnum.x*gnum.y*gnum.z;
	dim3 block, grid;
	CuCalGridN(numcell, block, grid);	
	CxSphDensityAtCell<<<grid, block>>>(dF, dvol, datt, n, gnum, gmin, glen);	// カーネル実行
	cudaThreadSynchronize();
}


/*!
* CUDAカーネルを使って並列計算:粒子の描画色を密度から計算
* @param[out] dcol 粒子色配列(デバイスメモリ)
* @param[in] dval  粒子物理量配列(デバイスメモリ)
* @param[in] n 粒子数
* @param[in] c1,c2 物理量が最小,最大のときの色(間の色は線形補間で求められる)
* @param[in] range x要素に最小値，y要素に最大値
*/
void CuColorScalar(float* dcol, int* datt, float* dval, int n, float3 c1, float3 c2, float2 range)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxColorScalar<<<grid, block>>>(dcol, datt, dval, n, c1, c2, range);	// カーネル実行
	cudaThreadSynchronize();
}
void CuColorVector(float* dcol, int* datt, float* dval, int n, float3 c1, float3 c2, float2 range)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxColorVector<<<grid, block>>>(dcol, datt, dval, n, c1, c2, range);	// カーネル実行
	cudaThreadSynchronize();
}
/*!
* CUDAカーネルを使って並列計算:粒子の描画色設定 - 一定の色
* @param[out] dcol 粒子色配列(デバイスメモリ)
* @param[in] col 描画色
* @param[in] n 粒子数
*/
void CuColorConstant(float* dcol, int* datt, float3 col, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxColorConstant<<<grid, block>>>(dcol, datt, col, n);	// カーネル実行
	cudaThreadSynchronize();
}


/*!
 * 各粒子のグリッドハッシュ値計算(近傍探索用)
 * @param[out] dhash 各粒子のグリッドハッシュ値を格納した配列
 * @param[out] dsortedidx 各粒子のインデックスを格納した配列(後からハッシュ値でソートされる -> 現時点ではまだソード済みではない)
 * @param[in] dpos 粒子位置を格納した配列
 * @param[in] n 粒子数
 */
void CuCalcHash(uint* dhash, uint* dindex, float* dpos, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算
	CxCalcHash<<<grid, block>>>(dhash, dindex, dpos, n);	// カーネル実行
	cudaThreadSynchronize();
}

/*!
 * thrust::sort_by_keyによるハッシュ値に基づくソート
 * @param[in] dhash ハッシュ値
 * @param[in] dindex インデックス(パーティクル，ポリゴンなど)
 * @param[in] n データ数
 */
void CuSort(unsigned int* dhash, uint* dindex, uint n)
{
	thrust::sort_by_key(thrust::device_ptr<unsigned int>(dhash),
					    thrust::device_ptr<unsigned int>(dhash+n),
					    thrust::device_ptr<unsigned int>(dindex));
	cudaThreadSynchronize();
}

/*!
 * パーティクル配列をソートされた順番に並び替え，各セルの始まりと終わりのインデックスを検索
 * @param[in] cell 近傍探索用グリッドデータ
 * @param[in] dpos 粒子位置
 * @param[in] dvel 粒子速度
 */
void CuReorderDataAndFindCellStart(Cell cell, float* dpos, float* dvel, uint n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// 粒子数=スレッド数としてブロック/グリッドサイズを計算

	// セルスタート位置配列の初期化
	CUCHECK(cudaMemset(cell.dCellStart, 0xffffffff, cell.uNumCells*sizeof(uint)));

	// シェアードメモリのサイズ
	uint smemSize = sizeof(uint)*(THREAD_NUM+1);

	// カーネル実行
	CxReorderDataAndFindCellStartD<<<grid, block, smemSize>>>(cell, dpos, dvel, n);
	cudaThreadSynchronize();
}

//海老沢追加-----------------------------------------------------------------------------------------------------
//XPBDの伸び・せん断，曲げ・ねじれ制約の処理
//dpos:位置
//dmas:質量
//dlen:基準長
//dkss:伸び剛性
//dkbt:曲げ剛性
//dquat:姿勢(四元数)
//domega:基準ダルボーベクトル
//dlamb_ss:XPBDの伸び・せん断制約に用いるλ
//dlamb_bt:XPBDの曲げ・ねじれ制約に用いるλ
//dfix:固定点(毛髪の開始点)を示す配列(1なら固定点,0ならそれ以外)
//dt:タイムステップ
//n:粒子数
//iter:反復回数
//example_flag:形状によって，処理を一部変える
void CuXPBDConstraint(float* dpos,float* dcurpos,float* dmas, float* dlen, float* dkss,float* dkbt, float* dquat,float* dcurquat, float* domega, float* dlamb_ss,float* dlamb_bt,int* dfix, float dt,int n,int iter,bool example_flag) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	//XPBDの処理のために，λを0にする
	CxSetLambdaZero << <grid, block >> > (dlamb_ss,dlamb_bt, n);
	cudaThreadSynchronize();
	//伸び制約の反復
	for (int i = 0; i < iter; i++) {
		//全ての制約を同時に実行すると，衝突が発生するため，奇数と偶数に分けて実行する
		//偶数番目のidを実行
		CxStretchingShearConstraint << <grid, block >> > (dpos, dcurpos, dmas, dlen, dkss, dquat, dcurquat, dlamb_ss, dfix, dt, n, 0, i, example_flag);
		//奇数番目のidを実行
		CxStretchingShearConstraint << <grid, block >> > (dpos, dcurpos, dmas, dlen, dkss, dquat, dcurquat, dlamb_ss, dfix, dt, n, 1, i, example_flag);
		
		CxBendTwistConstraint << <grid, block >> > (dmas, dquat, dcurquat, domega, dkbt, dlamb_bt, dlen, dfix, dt, n, 0, i, example_flag);
		CxBendTwistConstraint << <grid, block >> > (dmas, dquat, dcurquat, domega, dkbt, dlamb_bt, dlen, dfix, dt, n, 1, i, example_flag);
	}
}

//衝突制約
//ここでは，一番実装が容易な球との衝突だけを扱う
//dpos:位置
//dvel:速度
//dfix:固定点(毛髪の開始点)を示す配列(1なら固定点,0ならそれ以外)
//center:毛髪との衝突を扱いたい球の中心
//rad:毛髪との衝突を扱いたい球の半径
//dt:タイムステップ
//n:粒子数
void CuCollisionConstraint(float* dpos, float* dvel, int* dfix, float3 center, float rad, float dt, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxCollisionConstraint << <grid, block >> > (dpos, dvel, dfix, center, rad, dt, n);
	cudaThreadSynchronize();
}

//海老沢追加
//時間積分
//位置ベース法に従い，現在の位置と位置修正後の位置から速度を求め，位置を更新
//dpos:位置(位置修正後)
//dcurpos:前ステップの位置(位置修正前)
//dvel:速度
//dt:タイムステップ
//n:粒子数
//vel_control:一定以下の速度の場合に切り捨てを行うかどうかを指定
void CuIntegrate(float* dpos,float* dcurpos,float* dvel,float dt,int n,bool vel_control) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxIntegrate << <grid, block >> > (dpos, dcurpos, dvel, dt, n, vel_control);
	cudaThreadSynchronize();
}

//海老沢追加
//風や重力をイメージした外力計算
//デバック用
void CuCalExternalForces(float* dpos,float*dvel,float* dmass,int* dfix,float3 gravity, float3 wind, float dt, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxCalExternalForces << <grid, block >> > (dpos, dvel, dmass, dfix, gravity, wind, dt, n);
	cudaThreadSynchronize();
}

//海老沢追加
//位置ベース法(拡張位置ベース法でなく，伸び・せん断制約のみ)
//デバック用
void CuPBDStretchingConstraint(float* dpos, float* dmas, float* dlen, float* dkss, float* dquat, int* dfix, int n, int iter) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	for (int i = 0; i < iter; i++) {
		CxStretchingConstraint<<<grid,block>>>(dpos, dmas, dlen, dkss, dquat, dfix, n, 0);
		cudaThreadSynchronize();
		CxStretchingConstraint<<<grid,block>>>(dpos, dmas, dlen, dkss, dquat, dfix, n, 1);
	}
}

//デバック用の配列特定位置の出力
void CuPrint3Dfloat(float* dpos,float* dvel,float* dacc,int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxPrint3Dfloat << <grid, block >> > (dpos, dvel, dacc, n);
	cudaThreadSynchronize();
}

//接線の更新
//kajiya-kayモデルでのレンダリングに利用
//dpos:位置
//dtang:エッジごとの接線
//dfix:固定点(毛髪の開始点)を示す配列(1なら固定点,0ならそれ以外)
//n:粒子数
void CuTangUpdate(float* dpos, float* dtang, int* dfix, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxTangUpdate << <grid, block >> > (dpos, dtang, dfix, n);
	cudaThreadSynchronize();
}

//角速度など初期からデバイスメモリに設定をしているものの初期値を0に設定
//dangvel:角速度
//dfss:エッジごとにかかる力(GlobalForceStepで求める)
//dpbf_lambda:密度制約の計算過程に必要なλをメモリ確保
//n:粒子数
void CuSetParametersZero(float* dangvel, float* dfss, float* dpbf_lambda, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxSetParametersZero << <grid, block >> > (dangvel,dfss,dpbf_lambda, n);
	cudaThreadSynchronize();
}

//角加速度の更新
//dangvel:角速度
//dquat:姿勢(四元数)
//dfix:固定点(毛髪の開始点)を示す配列(1なら固定点,0ならそれ以外)
//dt:タイムステップ
//n:粒子数
void CuAngVelUpdate(float* dangvel, float* dquat,int* dfix,float dt, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxAngVelUpdate << <grid, block >> > (dangvel, dquat,dfix, dt, n);
	cudaThreadSynchronize();
}

//各加速度の時間積分
//dangvel:角速度
//dcurquat:前ステップの姿勢(位置修正前)
//dquat:現在の姿勢(位置修正後)
//dfix:固定点(毛髪の開始点)を示す配列(1なら固定点,0ならそれ以外)
//dt:タイムステップ
//n:粒子数
//vel_control:角速度が一定以下なら切り捨てを行うかどうかを判定
void CuAngVelIntegrate(float* dangvel,float* dcurquat, float* dquat,int* dfix,float dt, int n,bool vel_control) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxAngVelIntegrate << <grid, block >> > (dangvel, dcurquat, dquat, dfix, dt, n, vel_control);
	cudaThreadSynchronize();
}

//基準となる密度の設定
//dpos:位置
//dRestDens:粒子ごとに設定する基準密度
//dvol:体積
//n:粒子数
void CuRestDensSet(float* dpos,float* dRestDens, float* dvol,float* dmas, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxRestDensSet << <grid, block >> > (dpos, dRestDens, dvol, dmas, n);
	cudaThreadSynchronize();
}

//一律の基準となる密度の設定
//デバック用
void CuRestTotalDens(float* drestdens,float dens, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxRestTotalDens << <grid, block >> > (drestdens, dens, n);
	cudaThreadSynchronize();
}

//グローバルフォースステップ
//重力などによりエッジにかかる力を求める
//dfss:エッジごとにかかる力
//dmass:質量
//last_index:毛髪ごとの最後の粒子のインデックスを格納
//gravity:重力
//num_elastic:ここでは，毛髪ごとに並列計算するため，毛髪の数を渡す
void CuGlobalForceStep(float* dpos,float* dfss,float* dmass, int* last_index, float3 gravity,float* ddens,float* drestdens,float* dvol, int num_elastic) {
	dim3 block, grid;
	CuCalGridN(num_elastic, block, grid);
	CxGlobalForceStep << <grid, block >> > (dpos,dfss, dmass, last_index, gravity, ddens, drestdens, dvol, num_elastic);
	cudaThreadSynchronize();
}

//ローカルフォースステップ
//グローバルフォースステップで求めたエッジごとの力から，変形を防ぐための基準長や姿勢を求める
//dpos:位置
//dlen:基準長
//dquat:姿勢
//dcurquat:前ステップの姿勢(シミュレーション開始前のcurquatはquatと一致するため，更新後の値を代入)
//dkss:伸び剛性
//dfix:固定点(毛髪の開始点)を示す配列(1なら固定点,0ならそれ以外)
//n:粒子数(エッジごとに並列計算)
void CuLocalForceStep(float* dpos, float* dlen, float* dquat,float* dcurquat, float* dkss, float* dfss, int* dfix, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxLocalForceStep << <grid, block >> > (dpos, dlen, dquat, dcurquat, dkss, dfss, dfix, n);
	cudaThreadSynchronize();
}

//グローバルトルクステップ
//毛髪ごとにフォースステップで生じたトルクを打ち消す基準ダルボーベクトルを求める
//dpos:位置
//dquat:姿勢
//domega:基準ダルボーベクトル
//dlen:基準長
//dkss:伸び剛性
//dkbt:曲げ剛性
//dfix:固定点(毛髪の開始点)を示す配列(1なら固定点,0ならそれ以外)
//last_index:毛髪ごとの最後の粒子のインデックスを格納
//num_elastic:ここでは，毛髪ごとに並列計算するため，毛髪の数を渡す
void CuGlobalTorqueStep(float* dpos, float* dquat, float* domega, float* dlen, float* dkss, float* dkbt, int* dfix, int* dlast_index, int num_elastic) {
	dim3 block, grid;
	CuCalGridN(num_elastic, block, grid);
	//下から求める
	CxGlobalTorqueStep << <grid, block >> > (dpos, dquat, domega, dlen, dkss, dkbt, dfix, dlast_index, num_elastic);
	//上から求める
	//CxGlobalTorqueStep_Upstair << <grid, block >> > (dpos, dquat, domega, dlen, dkss, dkbt, dfix, dlast_index, num_elastic);
	cudaThreadSynchronize();
}

//ローカルトルクステップ
//基準ダルボーベクトルを適切な形で正規化する
//dquat:姿勢
//domega:基準ダルボーベクトル
//deln:基準長
//dkbt:曲げ剛性
//dfix:固定点(毛髪の開始点)を示す配列(1なら固定点,0ならそれ以外)
//n:粒子数(基準ダルボーベクトルごとに並列計算)
void CuLocalTorqueStep(float* dquat,float* domega, float* dlen, float* dkbt, int* dfix, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxLocalTorqueStep << <grid, block >> > (dquat, domega, dlen, dkbt, dfix, n);
	cudaThreadSynchronize();
}

//密度制約の計算
//pos:位置
//ddens:現在の密度
//drestdens:基準密度
//dpbf_lambda:制約に用いるλ
//dvol:体積
//n:粒子数
void CuPbfConstraint(float* dpos,float* ddens,float* drestdens,float*dpbf_lambda,float*dvol,float* dmas,int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxSphDensity << <grid, block >> > (drestdens, ddens, dvol,dmas, n);//密度計算
	CxPbfLambda << <grid, block >> > (ddens,drestdens, dpbf_lambda, dvol,dmas, n);//制約に用いるλを求める
	cudaThreadSynchronize();
	CxPbfConstraint << <grid, block >> > (dpos, drestdens, dpbf_lambda, dvol,dmas, n);//制約処理
	cudaThreadSynchronize();
}

//PBFで解く場合の外力項の計算
//dacc:加速度
//datt:粒子属性(0で流体,1で境界)
//power:風などの力
//n:粒子数
void CuPbfExternalForces(float* dacc, int* datt, float3 power,bool wind_flag, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxPbfExternalForces << <grid, block >> > (dacc, datt, power, wind_flag, n);
	cudaThreadSynchronize();
}

//摩擦制約
void CuFrictionConstraint(float* dpos, float* dcurpos, float* drestdens, float* dvol, float* ddens, int* dfix, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	//各粒子との摩擦力を合計した後，静止摩擦かどうかを判定
	//CxFrictionConstraint << <grid, block >> > (dpos, dcurpos, drestdens, dvol, ddens, dfix, n);
	//各粒子と静止摩擦かを判定した後，合計
	CxFrictionAllParticlesConstraint << <grid, block >> > (dpos, dcurpos, drestdens, dvol, ddens, dfix, n);
	cudaThreadSynchronize();
}

//摩擦制約の後，姿勢を修正する
void CuFrictionConstraint_withQuat(float* dpos, float* dcurpos, float* drestdens, float* dvol, float* ddens, float* dquat, float* dlen, int* dfix, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxFrictionConstraint_withQuat << <grid, block >> > (dpos, dcurpos, drestdens, dvol, ddens, dquat, dlen, dfix, n);
	cudaThreadSynchronize();
}

//2頂点から姿勢を設定
void CuQuatSet(float* dpos, float* dquat, int* dfix, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxQuatSet << <grid, block >> > (dpos, dquat, dfix, n);
	cudaThreadSynchronize();
}

//トルクを計算したい
void CuCalcTorque(float* dpos,float* dmas, float* dquat, float* dfss, float* dlength,float* dkss, int* dfix, float3 gravity, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxCalcTorque << <grid, block >> > (dpos, dmas, dquat, dfss, dlength, dkss, dfix, gravity, n);
	cudaThreadSynchronize();
}

//--------------------------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// GPU処理補助関数
//-----------------------------------------------------------------------------
/*!
 * CUDAデバイスの設定 - idを直接指定
 * @param[in] id デバイスID
 */
void CuSetDevice(int id)
{
	int device_count = 0;
	cudaGetDeviceCount(&device_count);
	if(id < 0 || id >= device_count){
		id = 0;
	}
	cudaSetDevice(id);

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, id);

	std::cout << " ---- GPU Info ----" << std::endl;
	std::cout << " Device        : " << prop.name << std::endl;
	std::cout << " Global mem    : " << prop.totalGlobalMem << " Byte" << std::endl;
	std::cout << " Constant mem  : " << prop.totalConstMem  << " Byte" << std::endl;
	std::cout << " Thresds/Block : " << prop.maxThreadsPerBlock << std::endl;
	std::cout << std::endl;
	THREAD_NUM = prop.maxThreadsPerBlock; // 1ブロックあたりのスレッド数最大値
}

/*!
 * CUDAデバイスの設定
 *  - コマンドライン引数に基づきCUDAデバイスを設定((例)-device 0)
 * @param[in] argc コマンドライン引数の数
 * @param[in] argv コマンドライン引数リスト(argv[0]は実行ファイル名)
 */
void CuInit()
{
	CuSetDevice(0);
}


/*!
 * デバイスメモリの確保
 * @param[out] dPtr デバイスメモリへのポインタ
 * @param[in] size 確保サイズ(メモリ上のサイズ)
 */
void CuAllocateArray(void** dPtr, size_t size)
{
	CUCHECK(cudaMalloc(dPtr, size));
}

/*!
 * デバイスメモリの解放
 * @param[in] devPtr デバイスメモリへのポインタ
 */
void CuFreeArray(void* dPtr)
{
	CUCHECK(cudaFree(dPtr));
}

/*!
 * デバイスメモリ領域の初期化
 * @param[in] dPtr デバイスメモリへのポインタ
 * @param[in] val 初期値
 * @param[in] size 初期化する領域のサイズ(メモリ上のサイズ)
 */
void CuSetArrayValue(void* dPtr, int val, size_t size)
{
	CUCHECK(cudaMemset(dPtr, val, size));
}

/*!
 * デバイスメモリ間コピー
 * @param[in] dDst コピー先
 * @param[in] dSrc コピー元
 * @param[in] size コピーサイズ(メモリ上のサイズ)
 */
void CuCopyArrayD2D(void* dDst, void* dSrc, int size)
{
	CUCHECK(cudaMemcpy(dDst, dSrc, size, cudaMemcpyDeviceToDevice));
}


/*!
 * VBOをマッピング
 * @param[in] vbo VBO,PBO名
 */
void* CuMapGLBufferObject(cudaGraphicsResource** resource)
{
	void* ptr;
	CUCHECK(cudaGraphicsMapResources(1, resource, 0));
	size_t num_bytes;
	CUCHECK(cudaGraphicsResourceGetMappedPointer((void**)&ptr, &num_bytes, *resource));
	return ptr;
}

/*!
 * VBOをアンマップ
 * @param[in] vbo VBO,PBO名
 */
void CuUnmapGLBufferObject(cudaGraphicsResource* resource)
{
	CUCHECK(cudaGraphicsUnmapResources(1, &resource, 0));
}

/*!
 * PBO,VBOバッファをCUDAに登録
 * @param[in] vbo VBO,PBO名
 */
void CuRegisterGLBufferObject(unsigned int vbo, cudaGraphicsResource** resource)
{
	CUCHECK(cudaGraphicsGLRegisterBuffer(resource, vbo, cudaGraphicsMapFlagsNone));
}

/*!
 * PBO,VBOバッファをCUDAから削除
 * @param[in] vbo VBO,PBO名
 */
void CuUnregisterGLBufferObject(cudaGraphicsResource* resource)
{
	CUCHECK(cudaGraphicsUnregisterResource(resource));
}

/*!
 * デバイスからホストメモリへのコピー
 * @param[in] hDst コピー先ホストメモリ(最低size分確保されていること)
 * @param[in] dSrc コピー元デバイスメモリ
 * @param[in] vbo dSrcがVBOの場合，VBOのID．そうでない場合は0を指定
 * @param[in] size コピーサイズ(メモリ上のサイズ)
 */
void CuCopyArrayFromDevice(void* hDst, void* dSrc, cudaGraphicsResource** resource, int offset, int size)
{
	if(resource) dSrc = CuMapGLBufferObject(resource);

	CUCHECK(cudaMemcpy(hDst, (char*)dSrc+offset, size, cudaMemcpyDeviceToHost));

	if(resource) CuUnmapGLBufferObject(*resource);
}

/*!
 * ホストからデバイスメモリへのコピー
 * @param[in] dDst コピー先デバイスメモリ(最低size分確保されていること)
 * @param[in] hSrc コピー元ホストメモリ
 * @param[in] offset コピー先オフセット
 * @param[in] size コピーサイズ(メモリ上のサイズ)
 */
void CuCopyArrayToDevice(void* dDst, const void* hSrc, int offset, int size)
{
	CUCHECK(cudaMemcpy((char*)dDst+offset, hSrc, size, cudaMemcpyHostToDevice));
}

/*!
 * スレッド同期
 */
void CuThreadSync(void)
{
	CUCHECK(cudaThreadSynchronize());
}

/*!
 * デバイスプロパティの表示
 */
void CuDeviceProp(void)
{
	int n;	//デバイス数
	CUCHECK(cudaGetDeviceCount(&n));

	for(int i = 0; i < n; ++i){
		cudaDeviceProp dev;

		// デバイスプロパティ取得
		CUCHECK(cudaGetDeviceProperties(&dev, i));

		printf("device %d\n", i);
		printf(" device name : %s\n", dev.name);
		printf(" total global memory : %d (MB)\n", (int)dev.totalGlobalMem/1024/1024);
		printf(" shared memory / block : %d (KB)\n", (int)dev.sharedMemPerBlock/1024);
		printf(" register / block : %d\n", dev.regsPerBlock);
		printf(" warp size : %d\n", dev.warpSize);
		printf(" max pitch : %d (B)\n", (int)dev.memPitch);
		printf(" max threads / block : %d\n", dev.maxThreadsPerBlock);
		printf(" max size of each dim. of block : (%d, %d, %d)\n", dev.maxThreadsDim[0], dev.maxThreadsDim[1], dev.maxThreadsDim[2]);
		printf(" max size of each dim. of grid  : (%d, %d, %d)\n", dev.maxGridSize[0], dev.maxGridSize[1], dev.maxGridSize[2]);
		printf(" clock rate : %d (MHz)\n", dev.clockRate/1000);
		printf(" total constant memory : %d (KB)\n", (int)dev.totalConstMem/1024);
		printf(" compute capability : %d.%d\n", dev.major, dev.minor);
		printf(" alignment requirement for texture : %d\n", (int)dev.textureAlignment);
		printf(" device overlap : %s\n", (dev.deviceOverlap ? "ok" : "not"));
		printf(" num. of multiprocessors : %d\n", dev.multiProcessorCount);
		printf(" kernel execution timeout : %s\n", (dev.kernelExecTimeoutEnabled ? "on" : "off"));
		printf(" integrated : %s\n", (dev.integrated ? "on" : "off"));
		printf(" host memory mapping : %s\n", (dev.canMapHostMemory ? "on" : "off"));

		printf(" compute mode : ");
		if(dev.computeMode == cudaComputeModeDefault) printf("default mode (multiple threads can use) \n");
		else if(dev.computeMode == cudaComputeModeExclusive) printf("exclusive mode (only one thread will be able to use)\n");
		else if(dev.computeMode == cudaComputeModeProhibited) printf("prohibited mode (no threads can use)\n");

	}

	//printf("Device with Maximum GFLOPS : %d\n", gpuGetMaxGflopsDeviceId());
}

/*!
 * thrust::exclusive_scanの呼び出し
 * @param[out] dScanData scan後のデータ
 * @param[in] dData 元データ
 * @param[in] num データ数
 */
void CuScan(float* dScanData, float* dData, int num)
{
	thrust::exclusive_scan(thrust::device_ptr<float>(dData),
		thrust::device_ptr<float>(dData + num),
		thrust::device_ptr<float>(dScanData));
}

/*!
 * デバッグ用 : スカラー値が入った配列の値の平均値を求めて返す
 * @param[out] dcol 粒子色配列(デバイスメモリ)
 * @param[in] ddens 粒子密度配列(デバイスメモリ)
 * @param[in] n 粒子数
 * @param[in] userparam ユーザーパラメータ(任意)
 */
float CuCalAverage(float* data, int n)
{
	if (n == 0) return 0;
	float avg = 0.0f;

	float* data_scan = 0;
	CuAllocateArray((void**)&data_scan, n * sizeof(float));

	// 合計値を求めるためにscan(prefix sum)を計算
	CuScan(data_scan, data, n);

	// Exclusive scan (最後の要素が0番目からn-2番目までの合計になっている)なので，
	// Scan前配列の最後(n-1番目)の要素と合計することでポリゴン数を計算
	float lval, lsval;
	CUCHECK(cudaMemcpy((void*)&lval, (void*)(data + n - 1), sizeof(float), cudaMemcpyDeviceToHost));
	CUCHECK(cudaMemcpy((void*)&lsval, (void*)(data_scan + n - 1), sizeof(float), cudaMemcpyDeviceToHost));
	float total = lval + lsval;
	avg = total / n;

	if (data_scan != 0) CuFreeArray(data_scan);

	return avg;
}
float CuCalAverageV(float* vdata, int n)
{
	if (n == 0) return 0;
	float avg = 0.0f;

	float* data = 0;
	float* data_scan = 0;
	CuAllocateArray((void**)&data, n * sizeof(float));
	CuAllocateArray((void**)&data_scan, n * sizeof(float));

	dim3 block, grid;
	CuCalGridN(n, block, grid);	// データ数=スレッド数としてブロック/グリッドサイズを計算
	CxVectorToScalar<<<grid, block>>>(vdata, data, n);	// カーネル実行
	cudaThreadSynchronize();

	// 合計値を求めるためにscan(prefix sum)を計算
	CuScan(data_scan, data, n);

	// Exclusive scan (最後の要素が0番目からn-2番目までの合計になっている)なので，
	// Scan前配列の最後(n-1番目)の要素と合計することでポリゴン数を計算
	float lval, lsval;
	CUCHECK(cudaMemcpy((void*)&lval, (void*)(data + n - 1), sizeof(float), cudaMemcpyDeviceToHost));
	CUCHECK(cudaMemcpy((void*)&lsval, (void*)(data_scan + n - 1), sizeof(float), cudaMemcpyDeviceToHost));
	float total = lval + lsval;
	avg = total / n;

	if (data_scan != 0) CuFreeArray(data_scan);

	return avg;
}

}   // extern "C"
