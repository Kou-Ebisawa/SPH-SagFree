/*!
  @file sph_kernel.cu

  @brief CUDA : SPH法
         - CUDAカーネル及びデバイス関数を記述
         - ホスト関数が書かれた*.cuファイルからのみインクルード(決してcppからインクルードしないように)

  @author Makoto Fujisawa
  @date 2023-02
 */


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "helper_math.h"
#include <math_constants.h>

#include "cuda_utils.h"
#include "cuda_utils.cu"


//-----------------------------------------------------------------------------
// 定数(デバイスメモリ)
//-----------------------------------------------------------------------------
__device__ __constant__ SceneParameter params;	// 各種パラメータ


//-----------------------------------------------------------------------------
// device関数 - デバイス(GPU)で実行・デバイス関数からのみ呼び出し可
//-----------------------------------------------------------------------------
/*!
* グリッド位置計算
* @param[in] p 座標
* @return グリッド座標
*/
__device__
inline int3 calcGridPos(float3 p)
{
    int3 grid;
    grid.x = floor((p.x-params.cell.WorldOrigin.x)/params.cell.CellWidth.x);
    grid.y = floor((p.y-params.cell.WorldOrigin.y)/params.cell.CellWidth.y);
    grid.z = floor((p.z-params.cell.WorldOrigin.z)/params.cell.CellWidth.z);

    grid.x = min(max(grid.x, 0), params.cell.GridSize.x-1);
    grid.y = min(max(grid.y, 0), params.cell.GridSize.y-1);
    grid.z = min(max(grid.z, 0), params.cell.GridSize.z-1);

    return grid;
}

/*!
* グリッド座標から1次元配列中での位置を計算
* @param[in] gridPos グリッド座標
* @return アドレス
*/
__device__
inline uint calcGridHash(int3 gridPos)
{
    return __umul24(__umul24(gridPos.z, params.cell.GridSize.y), params.cell.GridSize.x)+__umul24(gridPos.y, params.cell.GridSize.x)+gridPos.x;
}

/*!
* 粒子の衝突処理洋館数
* @param[inout] p,v 粒子位置,速度
* @param[in] dt タイムステップ幅
*/
__device__
void collision(float3 &p, float3 &v, float dt)
{
    float d;
    float3 nrm, cp;
    float res = params.res;

    // ボックス形状のオブジェクトとの衝突
#if MAX_BOX_NUM
    for(int i = 0; i < params.num_box; ++i){
        if(params.box[i].flg == 0) continue;
        collisionPointBox(p, params.box[i].cen, params.box[i].ext+make_float3(params.particle_radius), params.box[i].rot, params.box[i].inv_rot, cp, d, nrm);
        if(d < 0.0){
            res = (res > 0) ? (res*fabs(d)/(dt*length(v))) : 0.0f;
            v -= (1+res)*nrm*dot(nrm, v);
            p = cp;
        }
    }
#endif

    // 球形状のオブジェクトとの衝突
#if MAX_SPHERE_NUM
    for(int i = 0; i < params.num_sphere; ++i){
        if(params.sphere[i].flg == 0) continue;
        collisionPointSphere(p, params.sphere[i].cen, params.sphere[i].rad+params.particle_radius, cp, d, nrm);
        if(d < 0.0){
            res = (res > 0) ? (res*fabs(d)/(dt*length(v))) : 0.0f;
            v -= (1+res)*nrm*dot(nrm, v);
            p = cp;
        }
    }
#endif

    // シミュレーション空間の境界(AABB)との衝突
    float3 l0 = params.boundary_min;
    float3 l1 = params.boundary_max;
    if(distPointAABB(p, 0.5*(l1+l0), 0.5*(l1-l0), cp, d, nrm)){
        res = (res > 0) ? (res*fabs(d)/(dt*length(v))) : 0.0f;
        v -= (1+res)*nrm*dot(nrm, v);
        p = cp;
    }
}


//-----------------------------------------------------------------------------
// global関数 - デバイス(GPU)で実行・ホスト関数からのみ呼び出し可
//-----------------------------------------------------------------------------
/*!
 * SPH法による粒子密度の計算(Poly6カーネル)
 * @param[out] ddens 粒子密度
 * @param[in] dvol 粒子体積
 * @param[in] n 粒子数
 */
__global__ 
void CxSphDensity(float*drestdens,float* ddens, float* dvol, int n)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    float3 pos0 = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    //float m = params.mass;
    float a = params.aw;
    float rest_dens = params.rest_dens;
    //インデックスの計算
    uint sid = params.cell.dSortedIndex[id];
    //海老沢SPH追加
    rest_dens = drestdens[sid];

    // 粒子を中心として半径h内に含まれるグリッド(caclGridPos内で境界処理あり)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0-make_float3(h));
    grid_pos1 = calcGridPos(pos0+make_float3(h));

    // 周囲のグリッドセルも含めて近傍探索して密度計算
    float dens = 0.0f;
    for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
        for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
            for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // グリッドハッシュ値

                // セル内のパーティクルのスタートインデックス
                uint startIndex = params.cell.dCellStart[ghash];
                if(startIndex != 0xffffffff){	// セルが空でないかのチェック
                    // セル内のパーティクルで反復
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for(uint j = startIndex; j < endIndex; ++j){
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 rij = pos0-pos1;
                        float r = length(rij);
                        if(r <= h){
                            // Poly6カーネルで密度を計算 (rho = Σ m Wij)
                            float q = h*h-r*r;
                            // 流体粒子はrest_dens*dvol[sj] = mとなるように設定してある
                            // 境界粒子は想体積と初期密度から複数層境界粒子があった場合の仮想質量Φ=ρ0*Vbを求めて使う 
                            float m = rest_dens*dvol[sj];

                            dens += m*a*q*q*q;
                        }
                    }
                }
            }
        }
    }
    ddens[sid] = dens;
}

/*!
 * 粒子に働く力の計算
 *  - 重力，浮力など
 * @param[out] dpres 粒子圧力
 * @param[in] ddens 粒子密度
 * @param[in] n 粒子数
 */
__global__ 
void CxSphPressure(float* drestdens,float* dpres, float* ddens, int n)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    float p = 0.0f;

    // 理想気体の状態方程式に基づく計算[Muller2003]
    //p = params.gas_k * (ddens[id] - params.rest_dens);

    // Tait方程式に基づく計算:WCSPH[Becker2007]
    float rdens = ddens[id]/params.rest_dens;
    //海老沢SPH追加
    rdens = ddens[id] / drestdens[id];

    //p = params.B * (powf(rdens, params.gamma)-1.0f);
    // gamma=7固定の場合
    p = params.B * (rdens*rdens*rdens*rdens*rdens*rdens*rdens-1.0f);

    // 負圧の場合は0にする
    p = clamp(p, 0.0, 1.0e6);

    // 粒子圧力の更新(グローバルメモリを書き換える)
    dpres[id] = p;
}



/*!
* 粒子渦度の計算
*  - vorticity confinementによる乱流表現のための渦度計算
* @param[out] dvort 各粒子の渦度ベクトル(デバイスメモリ)
* @param[in] dvel 粒子速度
* @param[in] ddens 粒子密度
* @param[in] dvol 粒子体積
* @param[in] datt 粒子属性(0で流体,1で境界)
* @param[in] n 粒子数
*/
__global__ 
void CxSphVorticity(float* dvort, float* dvel, float* ddens, float* dvol, int* datt, int n)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    uint sid = params.cell.dSortedIndex[id];
    if(datt[sid] != 0){  // 境界粒子の場合は渦度=0
        float3 v0 = make_float3(0.0f);
        dvort[DIM*sid+0] = v0.x;  dvort[DIM*sid+1] = v0.y; dvort[DIM*sid+2] = v0.z;
        return;
    }

    // 粒子iの変数値
    float3 pos0 = params.cell.dSortedPos[id];
    float3 vel0 = make_float3(dvel[DIM*sid], dvel[DIM*sid+1], dvel[DIM*sid+2]);

    float h = params.h;
    float m = params.mass;
    float a = params.ag;

    // パーティクル周囲のグリッド(caclGridPos内で境界処理あり)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0-make_float3(h));
    grid_pos1 = calcGridPos(pos0+make_float3(h));

    // 周囲のグリッドセルも含めて近傍探索して圧力による力
    float3 f = make_float3(0.0f);
    for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
        for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
            for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);

                // セル内のパーティクルのスタートインデックス
                uint startIndex = params.cell.dCellStart[ghash];
                if(startIndex != 0xffffffff){	// セルが空でないかのチェック
                    // セル内のパーティクルで反復
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for(uint j = startIndex; j < endIndex; ++j){
                        uint sj = params.cell.dSortedIndex[j];
                        if(sj == sid) continue;

                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 vel1 = make_float3(dvel[DIM*sj], dvel[DIM*sj+1], dvel[DIM*sj+2]);

                        float3 rij = pos0-pos1;
                        float3 vji = vel1-vel0;
                        float dens1 = ddens[sj];

                        float r = length(rij);
                        if(r <= h && r > 0.0001){
                            float q = h-r;
                            float3 gw = a*q*q*rij/r;
                            f += m/dens1*cross(vji, gw);
                        }
                    }
                }
            }
        }
    }

    // 粒子の渦度ベクトルの更新(グローバルメモリを書き換える)
    dvort[DIM*sid+0] = f.x;
    dvort[DIM*sid+1] = f.y;
    dvort[DIM*sid+2] = f.z;
}


/*!
 * 粒子に働く力の計算
 *  - 密度を一定にするような圧力
 *  - 重力
 *  - vorticity confinement
 * @param[out] dacc 粒子に働く力(加速度)
 * @param[in] dvel 粒子速度配列
 * @param[in] ddens 粒子密度
 * @param[in] dpres 粒子圧力
 * @param[in] dvort 各粒子の渦度ベクトル
 * @param[in] dvol  粒子体積
 * @param[in] datt 粒子属性(0で流体,1で境界)
 * @param[in] n 粒子数
 */
//海老沢power追加
__global__ 
void CxSphForces(float* drestdens,float* dacc, float* dvel, float* ddens, float* dpres, float* dvort, float* dvol, int* datt,float3 power, int n)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    uint sid = params.cell.dSortedIndex[id];
    if(datt[sid] != 0){  // 境界粒子の場合は粒子にかかる力=0
        float3 v0 = make_float3(0.0f);
        dacc[DIM*sid+0] = v0.x;  dacc[DIM*sid+1] = v0.y; dacc[DIM*sid+2] = v0.z;
        return;
    }

    // 粒子iの変数値
    float3 pos0 = params.cell.dSortedPos[id];
    float3 omega0 = make_float3(dvort[DIM*sid], dvort[DIM*sid+1], dvort[DIM*sid+2]);
    float dens0 = ddens[sid];
    float pres0 = dpres[sid];
    //int3 grid = calcGridPos(pos0);
    float prsi = pres0/(dens0*dens0);

    float h = params.h;
    //float m = params.mass;
    float a = params.ag;
    float rest_dens = params.rest_dens;
    //海老沢SPH追加
    rest_dens = drestdens[sid];

    // パーティクル周囲のグリッド(caclGridPos内で境界処理あり)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0-make_float3(h));
    grid_pos1 = calcGridPos(pos0+make_float3(h));

    // 周囲のグリッドセルも含めて近傍探索して圧力による力
    float3 f = make_float3(0.0f);
    float3 eta = make_float3(0.0f);
    for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
        for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
            for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);

                // セル内のパーティクルのスタートインデックス
                uint startIndex = params.cell.dCellStart[ghash];
                if(startIndex != 0xffffffff){	// セルが空でないかのチェック
                    // セル内のパーティクルで反復
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for(uint j = startIndex; j < endIndex; ++j){
                        uint sj = params.cell.dSortedIndex[j];
                        if(sj == sid) continue;

                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 omega1 = make_float3(dvort[DIM*sj], dvort[DIM*sj+1], dvort[DIM*sj+2]);

                        float3 rij = pos0-pos1;

                        float dens1 = ddens[sj];
                        float pres1 = dpres[sj];
                        float prsj = pres1/(dens1*dens1);

                        float r = length(rij);
                        if(r <= h && r > 0.0001){
                            float q = h-r;
                            float3 gw = a*q*q*rij/r;

                            // 流体粒子はrest_dens*dvol[sj] = mとなるように設定してある
                            // 境界粒子は想体積と初期密度から複数層境界粒子があった場合の仮想質量Φ=ρ0*Vbを求めて使う 
                            float m = rest_dens*dvol[sj];
                            //海老沢SPH追加
                           // m = 0.01;

                            f += -m*(prsi+prsj)*gw; // 圧力項の計算
                            eta += (m/dens1)*length(omega1)*gw;
                        }
                    }
                }
            }
        }
    }

    // 重力
    f += params.gravity+power;
    //f += power;

    // Vorticity Confinement
    if(length(eta) > 1e-3){
        f += params.vorticity*dens0*cross(normalize(eta), omega0);
    }

    // 粒子にかかる外力(加速度)の更新(グローバルメモリを書き換える)
    dacc[DIM*sid+0] = f.x; 
    dacc[DIM*sid+1] = f.y;
    dacc[DIM*sid+2] = f.z;
}

/*!
* 粒子に働く粘性力の計算
*  - XSPHでなく粘性項を力(加速度)として計算する方法[Becker2007]
* @param[inout] dacc 粒子に働く力(加速度)
* @param[in] dvel 粒子速度
* @param[in] ddens 粒子密度
* @param[in] dvol 粒子体積
* @param[in] datt 粒子属性(0で流体粒子，それ以外で境界粒子)
* @param[in] n 粒子数
*/
__global__ 
void CxSphViscosity(float* drestdens,float* dacc, float* dvel, float* ddens, float* dvol, int* datt, int n)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    uint sid = params.cell.dSortedIndex[id];
    if(datt[sid] != 0) return;

    // 粒子iの変数値
    float3 pos0 = params.cell.dSortedPos[id];
    float3 vel0 = make_float3(dvel[DIM*sid], dvel[DIM*sid+1], dvel[DIM*sid+2]);
    float dens0 = ddens[sid];

    float h = params.h;
    float m = params.mass;
    float a = params.ag;
    float alpha = params.viscosity;  // 粘性定数．論文だと[0.08, 0.5]と行っているがこれだと大きすぎる...
    float cs = 88.5;
    float eps = 0.001*h*h;
    float rest_dens = params.rest_dens;
    //海老沢SPH追加
    rest_dens = drestdens[sid];

    // パーティクル周囲のグリッド(caclGridPos内で境界処理あり)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0-make_float3(h));
    grid_pos1 = calcGridPos(pos0+make_float3(h));

    // 周囲のグリッドセルも含めて近傍探索して圧力による力
    float3 f = make_float3(0.0f);
    for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
        for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
            for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);

                // セル内のパーティクルのスタートインデックス
                uint startIndex = params.cell.dCellStart[ghash];
                if(startIndex != 0xffffffff){	// セルが空でないかのチェック
                    // セル内のパーティクルで反復
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for(uint j = startIndex; j < endIndex; ++j){
                        uint sj = params.cell.dSortedIndex[j];
                        if(sj == sid) continue;

                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 vel1 = make_float3(dvel[DIM*sj], dvel[DIM*sj+1], dvel[DIM*sj+2]);
                        float3 rij = pos0-pos1;

                        float dens1 = ddens[sj];
                        float vx = dot(vel0-vel1, rij);

                        float r = length(rij);
                        if(r <= h && r > 0.0001 && vx < 0.0f){
                            float nu = (2.0f*alpha*h*cs)/(dens0+dens1);
                            float visc = -nu*(vx/(r*r+eps));
                            float q = h-r;
                            //float m = rest_dens*dvol[sj];
                            f += -m*visc*a*q*q*rij/r;
                        }
                    }
                }
            }
        }
    }

    // 粒子にかかる外力(加速度)の更新(グローバルメモリを書き換える)
    dacc[DIM*sid+0] += f.x; 
    dacc[DIM*sid+1] += f.y;
    dacc[DIM*sid+2] += f.z;
}

/*!
* 粒子に働く力の計算
*  - XSPH Artificial Viscosityによる速度更新[Schechter2012]
* @param[inout] dvel 粒子速度
* @param[in] ddens 粒子密度
* @param[in] dvol 粒子体積
* @param[in] datt 粒子属性(0で流体粒子，それ以外で境界粒子)
* @param[in] n 粒子数
*/
__global__ 
void CxSphXSPHViscosity(float* dvel, float* ddens, float* dvol, int* datt, int n)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    uint sid = params.cell.dSortedIndex[id];
    if(datt[sid] != 0) return;

    // 粒子iの変数値
    float3 pos0 = params.cell.dSortedPos[id];
    float3 vel0 = make_float3(dvel[DIM*sid], dvel[DIM*sid+1], dvel[DIM*sid+2]);

    float h = params.h;
    float m = params.mass;
    float a = params.aw;
    float eps = params.viscosity;
    float rest_dens = params.rest_dens;

    // パーティクル周囲のグリッド(caclGridPos内で境界処理あり)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0-make_float3(h));
    grid_pos1 = calcGridPos(pos0+make_float3(h));

    // 周囲のグリッドセルも含めて近傍探索して圧力による力
    float3 dv = make_float3(0.0f);
    for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
        for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
            for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);

                // セル内のパーティクルのスタートインデックス
                uint startIndex = params.cell.dCellStart[ghash];
                if(startIndex != 0xffffffff){	// セルが空でないかのチェック
                    // セル内のパーティクルで反復
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for(uint j = startIndex; j < endIndex; ++j){
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 vel1 = make_float3(dvel[DIM*sj], dvel[DIM*sj+1], dvel[DIM*sj+2]);
                        float3 rij = pos0-pos1;

                        float dens1 = ddens[sj];

                        float r = length(rij);
                        if(r <= h && r > 0.0001){
                            float q = h*h-r*r;
                            //float m = rest_dens*dvol[sj];
                            dv += (m/dens1)*(vel1-vel0)*a*q*q*q;
                        }
                    }
                }
            }
        }
    }

    dv *= eps;

    // 粒子速度の更新(グローバルメモリを書き換える)
    dvel[DIM*sid+0] += dv.x; 
    dvel[DIM*sid+1] += dv.y;
    dvel[DIM*sid+2] += dv.z;
}


/*!
 * 粒子を前進オイラー法で移動させる
 *  - 位置の速度による積分
 *  - 境界処理を含む
 * @param[inout] dpos 粒子位置
 * @param[inout] dvel 粒子速度
 * @param[in] dacc 粒子に働く力(加速度)
 * @param[in] datt 粒子属性(0で流体粒子，それ以外で境界粒子)
 * fix:海老沢追加 1ならば固定点
* @param[in] n 粒子数
 */
__global__ 
void CxSphIntegrate(float* dpos, float* dvel, float* dacc, int* datt,int* dfix, int n)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)
    if(datt[id] != 0){  // 境界粒子の場合は速度を0にして位置は変えない
        float3 v0 = make_float3(0.0f);
        dvel[DIM*id+0] = v0.x;  dvel[DIM*id+1] = v0.y; dvel[DIM*id+2] = v0.z;
        return;
    }
    //海老沢追加
    //固定点ならば、速度を0にして処理をスキップ
    if (dfix[id] == 1||id==0||dfix[id-1]==1) {
        float3 v0 = make_float3(0.0f);
        dvel[DIM * id + 0] = v0.x;  dvel[DIM * id + 1] = v0.y; dvel[DIM * id + 2] = v0.z;
        return;
    }
    // 粒子位置,速度,力
    float3 p = make_float3(dpos[DIM*id+0], dpos[DIM*id+1], dpos[DIM*id+2]);
    float3 v = make_float3(dvel[DIM*id+0], dvel[DIM*id+1], dvel[DIM*id+2]);
    float3 a = make_float3(dacc[DIM*id+0], dacc[DIM*id+1], dacc[DIM*id+2]);
    float dt = params.dt;

    // 更新位置，速度の更新
    v += dt*a;
    p += dt*v;

    // 周囲境界との衝突処理
    collision(p, v, dt);

    // 粒子位置・速度の更新(グローバルメモリを書き換える)
    dpos[DIM*id+0] = p.x;  dpos[DIM*id+1] = p.y; dpos[DIM*id+2] = p.z;
    dvel[DIM*id+0] = v.x;  dvel[DIM*id+1] = v.y; dvel[DIM*id+2] = v.z;
}

//海老沢PBF追加
//SPHをPBFで解く場合の圧力計算を除いた場合
__global__
void CxPbfExternalForces(float* dacc, int* datt,float3 power, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return;

    uint sid = params.cell.dSortedIndex[id];
    if (datt[sid] != 0) {  // 境界粒子の場合は粒子にかかる力=0
        float3 v0 = make_float3(0.0f);
        dacc[DIM * sid + 0] = v0.x;  dacc[DIM * sid + 1] = v0.y; dacc[DIM * sid + 2] = v0.z;
        return;
    }

    float3 acc = make_float3(dacc[DIM * id], dacc[DIM * id + 1], dacc[DIM * id + 2]);
    acc = params.gravity + power;

    dacc[DIM * id] = acc.x;
    dacc[DIM * id + 1] = acc.y;
    dacc[DIM * id + 2] = acc.z;
}

/*!
* 粒子を前進オイラー法で移動させる
*  - 速度の時間積分のみ, XSPHのみ，境界処理なし
* @param[inout] dvel 粒子速度
* @param[in] dacc 粒子に働く力(加速度)
* @param[in] datt 粒子属性(0で流体粒子，それ以外で境界粒子)
* @param[in] n 粒子数
*/
__global__ 
void CxSphIntegrateVelocity(float* dvel, float* dacc, int* datt, int n)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)
    if(datt[id] != 0){  // 境界粒子の場合は速度を0にして位置は変えない
        float3 v0 = make_float3(0.0f);
        dvel[DIM*id+0] = v0.x;  dvel[DIM*id+1] = v0.y; dvel[DIM*id+2] = v0.z;
        return;
    }

    // 粒子位置,速度,力
    float3 v = make_float3(dvel[DIM*id+0], dvel[DIM*id+1], dvel[DIM*id+2]);
    float3 a = make_float3(dacc[DIM*id+0], dacc[DIM*id+1], dacc[DIM*id+2]);
    float dt = params.dt;

    // 更新位置，速度の更新
    v += dt*a;

    // 粒子位置・速度の更新(グローバルメモリを書き換える)
    dvel[DIM*id+0] = v.x;  dvel[DIM*id+1] = v.y; dvel[DIM*id+2] = v.z;
}
/*!
* 粒子を前進オイラー法で移動させる
*  - 位置の速度による積分のみ，XSPH用，境界処理を含む
* @param[inout] dpos 粒子位置
* @param[inout] dvel 粒子速度
* @param[in] datt 粒子属性(0で流体粒子，それ以外で境界粒子)
* @param[in] n 粒子数
*/
//海老沢 dfix追加
__global__ 
void CxSphIntegratePosition(float* dpos, float* dvel, int* datt,int* dfix, int n)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n || datt[id] != 0) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    //海老沢追加
    //固定点ならば、速度を0にして処理をスキップ
    if (dfix[id] == 1 || id == 0 || dfix[id - 1] == 1) {
        float3 v0 = make_float3(0.0f);
        dvel[DIM * id + 0] = v0.x;  dvel[DIM * id + 1] = v0.y; dvel[DIM * id + 2] = v0.z;
        return;
    }

    // 粒子位置,速度,力
    float3 p = make_float3(dpos[DIM*id+0], dpos[DIM*id+1], dpos[DIM*id+2]);
    float3 v = make_float3(dvel[DIM*id+0], dvel[DIM*id+1], dvel[DIM*id+2]);
    float dt = params.dt;

    // 更新位置，速度の更新
    p += dt*v;

    // 周囲境界との衝突処理
    collision(p, v, dt);

    // 粒子位置・速度の更新(グローバルメモリを書き換える)
    dpos[DIM*id+0] = p.x;  dpos[DIM*id+1] = p.y; dpos[DIM*id+2] = p.z;
    dvel[DIM*id+0] = v.x;  dvel[DIM*id+1] = v.y; dvel[DIM*id+2] = v.z;
}


/*!
* 境界粒子処理のための粒子体積計算
*  - 境界粒子は "Versatile Rigid-Fluid Coupling for Incompressible SPH", 2.2 式(3)の上のV_bi で計算
*  - 流体粒子は V=mass/rest_dens
* @param[out] dvol 粒子体積
* @param[in] datt 粒子属性(0で流体粒子，それ以外で境界粒子)
* @param[in] n 処理する粒子数(offsetからの相対的な位置)
* @param[in] v 流体粒子の場合の粒子体積値
*/
__global__ 
void CxSphCalVolume(float* dvol, int *datt, int n, float v)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    uint sid = params.cell.dSortedIndex[id];
    int att = datt[sid];
    if(att == 0){   // 流体粒子の場合
        dvol[sid] = v;
        return;
    }

    float3 pos0 = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    float m = params.mass;
    float a = params.aw;

    // 粒子を中心として半径h内に含まれるグリッド(caclGridPos内で境界処理あり)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0-make_float3(h));
    grid_pos1 = calcGridPos(pos0+make_float3(h));

    // 周囲のグリッドセルも含めて近傍探索して体積計算
    float mw = 0.0f;    // ΣmW
    for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
        for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
            for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // グリッドハッシュ値
                uint startIndex = params.cell.dCellStart[ghash];                // セル内のパーティクルのスタートインデックス
                if(startIndex != 0xffffffff){	// セルが空でないかのチェック
                    // セル内のパーティクルで反復
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for(uint j = startIndex; j < endIndex; ++j){
                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 rij = pos0-pos1;
                        float r = length(rij);
                        if(r <= h){
                            float q = h*h-r*r;
                            mw += m*a*q*q*q;
                        }
                    }
                }
            }
        }
    }

    // 計算した体積をグローバルメモリに書き込み
    dvol[sid] = m/mw;
}


/*!
* グリッド上での密度を計算(表面メッシュ生成用)
* @param[out] dF 密度値を格納するグリッドセル配列
* @param[in] dvol 粒子体積
* @param[in] datt 粒子属性(0で流体,1で境界)
* @param[in] n 粒子数
* @param[in] gnum グリッド数
* @param[in] gmin グリッド最小座標
* @param[in] glen グリッド幅
*/
__global__
void CxSphDensityAtCell(float* dF, float* dvol, int* datt, int n, 
                        int3 gnum, float3 gmin, float3 glen)
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    int3 gridPos = calcGridPos(id, gnum);

    if(gridPos.x < gnum.x && gridPos.y < gnum.y && gridPos.z < gnum.z){
        float3 pos0;    // グリッドセル中心座標
        pos0.x = gmin.x+(gridPos.x)*glen.x;
        pos0.y = gmin.y+(gridPos.y)*glen.y;
        pos0.z = gmin.z+(gridPos.z)*glen.z;

        float h = params.effective_radius;
        float m = params.mass;
        float a = params.aw;

        int3 grid_pos0, grid_pos1;
        grid_pos0 = calcGridPos(pos0-make_float3(h));
        grid_pos1 = calcGridPos(pos0+make_float3(h));

        float dens = 0.0f;
        for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
            for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
                for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                    int3 ngrid = make_int3(x, y, z);
                    uint ghash = calcGridHash(ngrid);   // グリッドハッシュ値

                    // セル内のパーティクルのスタートインデックス
                    uint startIndex = params.cell.dCellStart[ghash];
                    if(startIndex != 0xffffffff){	// セルが空でないかのチェック
                        // セル内のパーティクルで反復
                        uint endIndex = params.cell.dCellEnd[ghash];
                        for(uint j = startIndex; j < endIndex; ++j){
                            uint sj = params.cell.dSortedIndex[j];
                            if(datt[sj] != 0) continue; // 境界粒子は表面メッシュ生成には使わない
                            float3 pos1 = params.cell.dSortedPos[j];
                            float3 rij = pos0-pos1;
                            float r = length(rij);
                            if(r <= h){
                                // Poly6カーネルで密度を計算 (rho = Σ m Wij)
                                float q = h*h-r*r;
                                dens += m*a*q*q*q;
                            }
                        }
                    }
                }
            }
        }

        dF[gridPos.x+gridPos.y*gnum.x+gridPos.z*gnum.x*gnum.y] = dens;
    }

}


//-----------------------------------------------------------------------------
// 粒子の描画色の計算
//-----------------------------------------------------------------------------
/*!
* 粒子の描画色を粒子の持つ物理量から計算
* @param[out] dcol 粒子色配列(デバイスメモリ)
* @param[in] dval 粒子物理量配列(デバイスメモリ)
* @param[in] n 粒子数
*/
__global__ 
void CxColorScalar(float* dcol, int* datt, float* dval, int n, float3 c1, float3 c2, float2 range)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    // 物理量と属性
    float val = dval[id];
    int att = datt[id];

    // 粒子色の計算
    float l = range.y-range.x;
    float t = clamp((val-range.x)/l, 0.0f, 1.0f);
    float3 col = lerp(c1, c2, t);

    // 粒子色のグローバルメモリへの格納
    dcol[4*id+0] = col.x;  dcol[4*id+1] = col.y; dcol[4*id+2] = col.z;
    dcol[4*id+3] = (1.0f-(float)(att));
}
__global__ 
void CxColorVector(float* dcol, int* datt, float* dval, int n, float3 c1, float3 c2, float2 range)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    // 物理量と属性
    float val = length(make_float3(dval[DIM*id], dval[DIM*id+1], dval[DIM*id+2]));
    int att = datt[id];

    // 粒子色の計算
    float l = range.y-range.x;
    float t = clamp((val-range.x)/l, 0.0f, 1.0f);
    float3 col = lerp(c1, c2, t);

    // 粒子色のグローバルメモリへの格納
    dcol[4*id+0] = col.x;  dcol[4*id+1] = col.y; dcol[4*id+2] = col.z;
    dcol[4*id+3] = (1.0f-(float)(att));
}

/*!
* 粒子の描画色設定 : すべて同じ色
* @param[out] dcol 粒子色配列(デバイスメモリ)
* @param[in] ddens 粒子密度配列(デバイスメモリ)
* @param[in] n 粒子数
*/
__global__ 
void CxColorConstant(float* dcol, int* datt, float3 col, int n)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)
    
    // 粒子属性
    int att = datt[id];

    // 粒子色のグローバルメモリへの格納
    dcol[4*id+0] = col.x;  dcol[4*id+1] = col.y; dcol[4*id+2] = col.z;
    dcol[4*id+3] = (1.0f-(float)(att));
}

/*!
* デバッグ用 : ベクトル配列からベクトルの大きさ(スカラー値)の配列を計算
* @param[in] dvdata ベクトル値配列(デバイスメモリ)
* @param[out] dsdata スカラー値配列(デバイスメモリ)
* @param[in] n 粒子数
*/
__global__ 
void CxVectorToScalar(float* dvdata, float* dsdata, int n)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    // ベクトルの大きさを求める
    float val = length(make_float3(dvdata[DIM*id], dvdata[DIM*id+1], dvdata[DIM*id+2]));

    // 結果のグローバルメモリへの格納
    dsdata[id] = val;
}


//-----------------------------------------------------------------------------
// ハッシュ
//-----------------------------------------------------------------------------
/*!
 * 各粒子のグリッドハッシュ値計算
 * @param[out] dhash 各粒子のグリッドハッシュ値を格納した配列
 * @param[out] dsortedidx 各粒子のインデックスを格納した配列(後からハッシュ値でソートされる -> 現時点ではまだソード済みではない)
 * @param[in] dpos 粒子位置を格納した配列
 * @param[in] n 粒子数
 */
__global__
void CxCalcHash(uint* dhash, uint* dsortedidx, float* dpos, uint n)
{
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    // 粒子位置
    float3 p = make_float3(dpos[DIM*id+0], dpos[DIM*id+1], dpos[DIM*id+2]);

    // 粒子位置から含まれるグリッドセルのハッシュ値を計算
    int3 grid = calcGridPos(p);
    uint hash = calcGridHash(grid);

    dhash[id] = hash;
    dsortedidx[id] = id;
}

/*!
 * パーティクルデータをソートして，ハッシュ内の各セルの最初のアドレスを検索
 * @param[in] cell パーティクルグリッドデータ
 * @param[in] dpos 粒子位置配列
 * @param[in] dvel 粒子速度配列
 */
__global__
void CxReorderDataAndFindCellStartD(Cell cell, float* dpos, float* dvel, uint n)
{
    // シェアードメモリ
    extern __shared__ uint sharedHash[];	// サイズ : blockSize+1

    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    uint hash;
    if(id < n){ // シェアードメモリ使用のためにスレッド同期を行うので(id >= n)の時もreturnはしない
        hash = cell.dGridParticleHash[id];	// ハッシュ値
        sharedHash[threadIdx.x+1] = hash;	// ハッシュ値をシェアードメモリに格納

        if(id > 0 && threadIdx.x == 0){
            // 各シェアードメモリの最初は隣のグリッドのパーティクルのハッシュ値を格納
            sharedHash[0] = cell.dGridParticleHash[id-1];
        }
    }

    __syncthreads();    // スレッド同期(他のスレッドがシェアードメモリ格納処理を終えるまで待つ)

    if(id < n){
        // インデックス0である，もしくは，一つ前のパーティクルのグリッドハッシュ値が異なる場合，
        // パーティクルは分割領域の最初
        if(id == 0 || hash != sharedHash[threadIdx.x]){
            cell.dCellStart[hash] = id;
            if(id > 0){
                // 一つ前のパーティクルは，一つ前の分割領域の最後
                cell.dCellEnd[sharedHash[threadIdx.x]] = id;
            }
        }

        // インデックスが最後ならば，分割領域の最後
        if(id == cell.uNumParticles-1){
            cell.dCellEnd[hash] = id+1;
        }

        // 位置と速度のデータを並び替え
        // ソートしたインデックスで参照も可能だが探索時のグローバルメモリアクセスを極力抑えるためにデータそのものを並び替える
        uint sid = cell.dSortedIndex[id];
        float3 pos = make_float3(dpos[DIM*sid+0], dpos[DIM*sid+1], dpos[DIM*sid+2]);
        float3 vel = make_float3(dvel[DIM*sid+0], dvel[DIM*sid+1], dvel[DIM*sid+2]);

        cell.dSortedPos[id] = pos;
        cell.dSortedVel[id] = vel;
    }
}

//---------------------------------------------
//以下、海老沢追加 四元数はfloat4(x,y,z,w)として扱う(出村さん参考)

//出村さんのコードから
//二つの四元数の積を取る
__host__ __device__ float4 quatProduct(float4 a, float4 b) {
    return make_float4(
        a.x * b.w + a.w * b.x - a.z * b.y + a.y * b.z,
        a.y * b.w + a.z * b.x + a.w * b.y - a.x * b.z,
        a.z * b.w - a.y * b.x + a.x * b.y + a.w * b.z,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z
    );
}

//出村さんのコードから
//四元数の共役を取る
__host__ __device__ float4 quatConjugate(float4 quat) {
    return make_float4(-quat.x, -quat.y, -quat.z, quat.w);
}

//出村さんのコードから
//3次元ベクトルと四元数の積を取る
__host__ __device__ float3 rotVecByQuat(float3 vec, float4 quat) {
    float4 vecq = make_float4(vec, 0.0f);
    float4 vecq_dash = quatProduct(quatProduct(quat, vecq), quatConjugate(quat));
    return make_float3(vecq_dash);
};

//3次元ベクトルの長さを求める
__host__ __device__ float Length(float3 vec) {
    return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

//4次元ベクトルの長さを求める
__host__ __device__ float Length(float4 vec) {
    return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z + vec.w * vec.w);
}

//XPBDに用いるλを全て0にする
//dlamb_ss:伸び制約のλ
//dlamb_bt:曲げ制約のλ
__global__
void CxSetLambdaZero(float* dlamb_ss,float* dlamb_bt,int n) {
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    dlamb_ss[DIM * id] = dlamb_ss[DIM * id + 1] = dlamb_ss[DIM * id + 2] = 0.f;
    dlamb_bt[QUAT * id] = dlamb_bt[QUAT * id + 1] = dlamb_bt[QUAT * id + 2] = dlamb_bt[QUAT * id + 3] = 0.f;
}

//XPBDの伸び・せん断制約
//dpos:位置
//dmas:質量
//dlen:基準長
//dkss:曲げ剛性
//dquat:姿勢
//dlamb_ss:XPBDのλ
//fix:固定点かどうか(1ならば固定点)
//dt:タイムステップ
//n:粒子数
//odd_even:偶数のスレッドIDか奇数のスレッドIDかを判別
__global__
void CxStrechingShearConstraint(float* dpos, float* dmas, float* dlen, float* dkss, float* dquat, float* dlamb_ss, int* dfix, float dt, int n,int odd_even,int iter) {
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    id = id * 2 + odd_even;
    
    if (id >= n-1) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない) 最後の粒子はスキップ
    if (dfix[id] == 1 || dfix[id + 1] == 1) return;//弾性体毎に辺の数だけ行う
    //printf("%d\n", id);
    //if(iter==0) printf("Number %d\n", id);

    float3 pos0 = make_float3(dpos[DIM * id], dpos[DIM * id + 1], dpos[DIM * id + 2]);
    float3 pos1 = make_float3(dpos[DIM * (id + 1)], dpos[DIM * (id + 1) + 1], dpos[DIM * (id + 1) + 2]);
    //if(id==9) printf("id %d pos0 x:%f,y:%f,z:%f,pos1 x:%f,y:%f,z:%f\n",id, pos0.x, pos0.y, pos0.z, pos1.x, pos1.y, pos1.z);

    float mass = dmas[id];//現在は粒子の重さは全て等しいとしている
    float length = dlen[id];
    float kss = dkss[id];
    float4 quat = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);//四元数は(x,y,z,w)になるように受け取る
    float3 lambda_ss = make_float3(dlamb_ss[DIM * id], dlamb_ss[DIM * id + 1], dlamb_ss[DIM * id + 2]);
    //printf("Length %f", length);
    //quat = make_float4(0.f, 0.f, 0.f, 1.f);
    //if (id == 9) printf("quat x:%f,y:%f,z:%f,w:%f \n", quat.x, quat.y, quat.z, quat.w);

    float w1 = 1 / mass;//重み
    float w2 = 1 / mass;//重み
    float wq = 1/length;//辺の重み
    //wq = 1.0e2;//(角速度ありで気持ち悪い動き)
    //wq = 7.5e3;
    wq = 1.0;
    //wq = 2.5e2;
    //wq = length * 5;

    if (dfix[id-1] == 1) {
        w1 = 0;
        wq = 0;
    }

    //w=1.0;
    //wq=1.0e-4;

    float alpha = 1 / (kss * dt * dt);//α
    //alpha = 0;
    float3 ds3 = rotVecByQuat(make_float3(0.f, 0.f, 1.f), quat);

    float3 dist = pos1 - pos0;
    float distLen2 = dot(dist,dist);//分母のlはこれを意味してるかも

    float weight = w1 + w2 + length * length * (4 * wq + alpha);//分母
    float3 mole = length * (pos0 - pos1 + length * ds3 - alpha * length * lambda_ss);//分子
    //mole = length * (pos1 - pos0 + length * ds3 - alpha * length * lambda_ss);
    float3 lambda = mole / weight;//Δλ
    //if (id == 0) printf("lambda x:%f,y:%f,z:%f\n", lambda.x, lambda.y, lambda.z);

    //重み追加
    float3 delta_pos0 = -w1*lambda / length;//Δx0(論文と符号逆)
    float3 delta_pos1 = w2*lambda / length;//Δx1(論文と符号逆
    float4 q_e3_bar = quatProduct(quat, make_float4(0.f, 0.f, -1.f, 0.f));
    //float4 delta_quat = -wq*2.f * quatProduct(quatProduct(make_float4(lambda, 0.f), quat), make_float4(0.f, 0.f, -1.f, 0.f));//Δqs
    float4 inter_quat = quatProduct(make_float4(lambda, 0.f), q_e3_bar);
    float4 delta_quat = -wq * 2.f * inter_quat;
    float4 new_quat = quat + delta_quat;//更新後のqs
    new_quat = normalize(new_quat);

    if (dfix[id-1] == 0) {
        //printf("dpos %f", Length(delta_pos0));
        dpos[id * DIM] += delta_pos0.x;
        dpos[id * DIM + 1] += delta_pos0.y;
        dpos[id * DIM + 2] += delta_pos0.z;
    }

    dpos[(id + 1) * DIM] += delta_pos1.x;
    dpos[(id + 1) * DIM + 1] += delta_pos1.y;
    dpos[(id + 1) * DIM + 2] += delta_pos1.z;
    
    if (dfix[id] == 0) {
        dquat[id * QUAT] = new_quat.x;
        dquat[id * QUAT + 1] = new_quat.y;
        dquat[id * QUAT + 2] = new_quat.z;
        dquat[id * QUAT + 3] = new_quat.w;
    }

    dlamb_ss[id * DIM] += lambda.x;
    dlamb_ss[id * DIM + 1] += lambda.y;
    dlamb_ss[id * DIM + 2] += lambda.z;
}

//曲げねじれ制約の方向を決定
//delta_omega=cur_omega-rest_omega
//delta_omega_plus=cur_omega+rest_omega
__host__ __device__ 
int deltaOmegaSign(float4 delta_omega, float4 delta_omega_plus) {
    if (dot(delta_omega, delta_omega) > dot(delta_omega_plus, delta_omega_plus)) return -1;

    return 1;
}


//曲げねじれ制約の追加
//dquat:辺の姿勢(四元数)
//domega:基準ダルボーベクトル
//dkbt:曲げ剛性
//dlamb_bt:曲げ制約に用いるλ
//dfix:固定点かどうかを表す(1ならば固定点)
//dt:タイムステップ
//n:粒子数
//odd_even:偶数のスレッドIDか奇数のスレッドIDかを判別
__global__ 
void CxBendTwistConstraint(float* dquat, float* domega, float* dkbt, float* dlamb_bt, float* dlength, int* dfix, float dt, int n, int odd_even, int iter) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    id = id * 2 + odd_even;
    if (id >= n - 2) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない) 最後の2粒子はスキップ
    if (dfix[id + 1] == 1 || dfix[id + 2] == 1) return;//末端の粒子とそのひとつ前の粒子では、エッジが足りない

    float kbt = dkbt[id];
    float4 quat1 = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);
    float4 quat2 = make_float4(dquat[QUAT * (id + 1)], dquat[QUAT * (id + 1) + 1], dquat[QUAT * (id + 1) + 2], dquat[QUAT * (id + 1) + 3]);
    float4 rest_omega = make_float4(domega[QUAT * id], domega[QUAT * id + 1], domega[QUAT * id + 2], domega[QUAT * id + 3]);
    float4 lambda_bt = make_float4(dlamb_bt[QUAT * id], dlamb_bt[QUAT * id + 1], dlamb_bt[QUAT * id + 2], dlamb_bt[QUAT * id + 3]);

    float dlen1 = dlength[id];
    float dlen2 = dlength[id + 1];

    float wq1 = 1 / dlen1;//1/(1.0e-3)*dlen1
    float wq2 = 1 / dlen2;
    float wq = 1.0e-6;
    float alpha = 1 / (kbt * dt * dt);
    
    wq1 = wq2 = 1.0;

    if (dfix[id] == 1) {
        wq1 = 0;
    }

    float weight = wq1 + wq2 + alpha;
    float4 cur_omega = quatProduct(quatConjugate(quat1), quat2);
    int s = deltaOmegaSign(cur_omega - rest_omega, cur_omega + rest_omega);

    float4 delta_omega = cur_omega - (s * rest_omega);
    delta_omega.w = 0.f;//omega.w=0としておく
    float4 lambda = (-delta_omega - alpha * lambda_bt) / weight;
    //lambda.w = 0.f;

    //重み追加
    float4 delta_quat1 = wq1 * quatProduct(quat2, quatConjugate(lambda));
    //float4 delta_quat1 = wq1*quatProduct(quat2, lambda);
    float4 delta_quat2 = wq2 * quatProduct(quat1, lambda);
    float4 new_quat1 = normalize(quat1 + delta_quat1);
    float4 new_quat2 = normalize(quat2 + delta_quat2);

    //quat1の更新
    if (dfix[id] == 0) {
        dquat[QUAT * id] = new_quat1.x;
        dquat[QUAT * id + 1] = new_quat1.y;
        dquat[QUAT * id + 2] = new_quat1.z;
        dquat[QUAT * id + 3] = new_quat1.w;
    }
    //quat2の更新
    dquat[QUAT * (id + 1)] = new_quat2.x;
    dquat[QUAT * (id + 1) + 1] = new_quat2.y;
    dquat[QUAT * (id + 1) + 2] = new_quat2.z;
    dquat[QUAT * (id + 1) + 3] = new_quat2.w;

    //lambdaの更新
    dlamb_bt[QUAT * id] += lambda.x;
    dlamb_bt[QUAT * id + 1] += lambda.y;
    dlamb_bt[QUAT * id + 2] += lambda.z;
    dlamb_bt[QUAT * id + 3] += lambda.w;
}

//衝突制約の実装
__global__
void CxCollisionConstraint(float* dpos, float* dvel, int* dfix, float3 center, float rad, float dt, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n)return;
    if (dfix[id] == 1) return;

    float3 pos = make_float3(dpos[DIM * id], dpos[DIM * id + 1], dpos[DIM * id + 2]);
    float3 vel = make_float3(dvel[DIM * id], dvel[DIM * id + 1], dvel[DIM * id + 2]);

    float3 d = pos - center;
    float l = rad - length(d);
    if (l <= 0) return;

    float3 norm = normalize(d);

    dpos[DIM * id] += l * norm.x;
    dpos[DIM * id + 1] += l * norm.y;
    dpos[DIM * id + 2] += l * norm.z;

    float3 addVel = -dot(norm, vel) * norm;

    dvel[DIM * id] += addVel.x;
    dvel[DIM * id + 1] += addVel.y;
    dvel[DIM * id + 2] += addVel.z;
}

//時間積分
//dpos:更新後の位置
//dcurpos:更新前の位置
//dvel:速度(更新した速度を代入)
//dt:タイムステップ
//n:粒子数
__global__
void CxIntegrate(float* dpos, float* dcurpos,float* dvel,float dt,int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    float3 pos = make_float3(dpos[DIM * id], dpos[DIM * id + 1], dpos[DIM * id + 2]);
    float3 cur_pos = make_float3(dcurpos[DIM * id], dcurpos[DIM * id + 1], dcurpos[DIM * id + 2]);
    float3 vel = (pos-cur_pos)/dt;
    
    //時間積分の際に，一定速度以下の粒子については変化なしとして，固定してしまう
    if (length(vel) < VEL_EPSILON) {
        //速度を0に固定
        vel = make_float3(0.f);
        //位置を更新前の値に戻す
        pos = cur_pos;
        dpos[DIM * id] = pos.x;
        dpos[DIM * id + 1] = pos.y;
        dpos[DIM * id + 2] = pos.z;
    }

    dvel[DIM * id] = vel.x;
    dvel[DIM * id + 1] = vel.y;
    dvel[DIM * id + 2] = vel.z;

    //if (Length(make_float3(dvel[DIM * id], dvel[DIM * id + 1], dvel[DIM * id + 2])) > 2.5f) printf("Not Stable Number : %d\n", id);

    dcurpos[DIM * id] = pos.x;
    dcurpos[DIM * id + 1] = pos.y;
    dcurpos[DIM * id + 2] = pos.z;
}

//外力計算
//風をイメージした加速度を与える
//全ての粒子に同じ加速度を与える
//dpos:位置
__global__
void CxCalExternalForces(float* dpos,float* dvel,float* dmass,int* dfix,float3 gravity, float3 wind,float dt, int n){
int id = blockDim.x * blockIdx.x + threadIdx.x;
if (id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない){
if (dfix[id] == 1) return;
if (dfix[id - 1] == 1)return;

float mass = dmass[id];

dvel[DIM * id] += (gravity.x+wind.x) * dt;
dvel[DIM * id + 1] += (gravity.y+wind.y) * dt;
dvel[DIM * id + 2] += (gravity.z+wind.z) * dt;

dpos[DIM * id] += dvel[DIM * id] * dt;
dpos[DIM * id + 1] += dvel[DIM * id + 1] * dt;
dpos[DIM * id + 2] += dvel[DIM * id + 2] * dt;
}

//位置ベース法の伸び制約
//XPBDの伸び制約に問題があるようなので、確認用
__global__
void CxStretchingConstraint(float* dpos, float* dmas, float* dlen, float* dkss, float* dquat, int* dfix, int n,int odd_even) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    id = id * 2 + odd_even;

    if (id >= n - 1) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない) 最後の粒子はスキップ
    if (dfix[id + 1] == 1) return;//弾性体毎に辺の数だけ行う

    float3 pos0 = make_float3(dpos[DIM * id], dpos[DIM * id + 1], dpos[DIM * id + 2]);
    float3 pos1 = make_float3(dpos[DIM * (id + 1)], dpos[DIM * (id + 1) + 1], dpos[DIM * (id + 1) + 2]);
    //if(id==1) printf("pos0 x:%f,y:%f,z:%f \n", pos0.x, pos0.y, pos0.z);

    float mass = dmas[id];//現在は粒子の重さは全て等しいとしている
    float length = dlen[id];
    float kss = dkss[id];
    float4 quat = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);//四元数は(x,y,z,w)になるように受け取る
    //if(id==1) printf("mass:%f,length:%f,quat x:%f,y:%f,z:%f,w:%f \n", mass, length, quat.x, quat.y, quat.z, quat.w);

    float3 e3 = make_float3(0, 0, 1);
    float3 d3 = rotVecByQuat(e3, quat);

    float w = 1.0 / mass;//直線用
    float wq = 1.0e5;//直線用
    wq = 1 / length;

    //Δposを計算
    float3 gamma = (pos1 - pos0) / length - d3;
    gamma /= (2*w) / length + 4.0f * wq * length;
    //if (id == 9) printf("before pos.x:%f,pos.y:%f,pos.z:%f\n", pos0.x, pos0.y, pos0.z);
    float ks = 1.0;
    gamma *= ks;

    float3 delta_pos0 = w * gamma;
    float3 delta_pos1 = -w * gamma;
    //if (id == 9) printf("updated pos.x:%f,pos.y:%f,pos.z:%f\n", pos0.x+delta_pos0.x, pos0.y+delta_pos0.y, pos0.z+delta_pos0.z);

    //if (id == 1) printf("delta_pos0 x:%f,y:%f,z:%f \n", delta_pos0.x, delta_pos0.y, delta_pos0.z);

    // calc delta_q
    float4 e3q = make_float4(e3, 0.0f);
    float4 q_e3_bar = quatProduct(quat, quatConjugate(e3q)); // calc q*e3_bar
    float4 gammaq = make_float4(gamma, 0.0f);
    //if (id == 1) printf("gamma_quat x:%f,y:%f,z:%f,w:%f \n", gammaq.x, gammaq.y, gammaq.z, gammaq.w);
    float4 inter_quat = quatProduct(gammaq, q_e3_bar);
    float4 delta_quat = wq * length * quatProduct(gammaq, q_e3_bar);//2.0追加
    //if (id == 1) printf("inter quat  x:%f,y:%f,z:%f,w:%f \n", inter_quat.x, inter_quat.y, inter_quat.z, inter_quat.w);
    //if (id == 1) printf("delta_quat x:%f,y:%f,z:%f,w:%f \n", delta_quat.x, delta_quat.y, delta_quat.z, delta_quat.w);

    float4 new_quat = quat + delta_quat;//更新後のqs
    new_quat = normalize(new_quat);

    if (dfix[id] == 0) {
        //printf("dpos %f", Length(delta_pos0));
        dpos[id * DIM] += delta_pos0.x;
        dpos[id * DIM + 1] += delta_pos0.y;
        dpos[id * DIM + 2] += delta_pos0.z;
    }

    dpos[(id + 1) * DIM] += delta_pos1.x;
    dpos[(id + 1) * DIM + 1] += delta_pos1.y;
    dpos[(id + 1) * DIM + 2] += delta_pos1.z;

    dquat[id * QUAT] = new_quat.x;
    dquat[id * QUAT + 1] = new_quat.y;
    dquat[id * QUAT + 2] = new_quat.z;
    dquat[id * QUAT + 3] = new_quat.w;
}

//位置、速度、加速度を出力
__global__
void CxPrint3Dfloat(float* dpos, float* dvel, float* dacc, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return;

    printf("id %d\n", id);
    printf("pos x:%f,y:%f,z:%f\n", dpos[id * DIM], dpos[id * DIM + 1], dpos[id * DIM + 2]);
    printf("vel x:%f,y:%f,z:%f\n", dvel[id * DIM], dvel[id * DIM + 1], dvel[id * DIM + 2]);
    printf("acc x:%f,y:%f,z:%f\n", dacc[id * DIM], dacc[id * DIM + 1], dacc[id * DIM + 2]);
}

//接線の更新
__global__
void CxTangUpdate(float* dpos, float* dtang, int* dfix, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return;
    if (dfix[id] == 1)return;

    dtang[DIM * id] = dpos[DIM * id] - dpos[DIM * (id - 1)];
    dtang[DIM * id + 1] = dpos[DIM * id + 1] - dpos[DIM * (id - 1) + 1];
    dtang[DIM * id + 2] = dpos[DIM * id + 2] - dpos[DIM * (id - 1) + 2];
}

//重力のみをかける
__global__
void CxOnlyGravity(float* dvel, float* dmass, int* dfix,float dt, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return;
    if (dfix[id] == 1)return;
    float gravity = -9.81;

    dvel[DIM * id + 1] += dmass[DIM * id] * gravity * dt;
}

__global__
void CxSetParametersZero(float* dangvel,float* dfss,float* dpbf_lambda, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return;
    dangvel[DIM * id] = dangvel[DIM * id + 1] = dangvel[DIM * id + 2] = 0;
    dfss[DIM * id] = dfss[DIM * id + 1] = dfss[DIM * id + 2] = 0;
    dpbf_lambda[id] = 0.f;
}

//角加速度の更新
__global__
void CxAngVelUpdate(float* dangvel, float* dquat,int* dfix,float dt, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n-1) return;
    if (dfix[id + 1] == 1) return;

    float4 quat = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);
    float3 avel = make_float3(dangvel[DIM * id], dangvel[DIM * id + 1], dangvel[DIM * id + 2]);

    float4 avelq = make_float4(avel, 0.0f);
    quat += 0.5f * quatProduct(quat, avelq) * dt;
    quat = normalize(quat);

    dquat[QUAT * id] = quat.x;
    dquat[QUAT * id + 1] = quat.y;
    dquat[QUAT * id + 2] = quat.z;
    dquat[QUAT * id + 3] = quat.w;
}

//各加速度の時間積分
__global__
void CxAngVelIntegrate(float* dangvel,float* dcurquat, float* dquat,int* dfix,float dt, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n-1) return;
    if (dfix[id + 1] == 1) return;

    float4 quat = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);
    float4 cur_quat = make_float4(dcurquat[QUAT * id], dcurquat[QUAT * id + 1], dcurquat[QUAT * id + 2], dcurquat[QUAT * id + 3]);

    float4 delta_rot = quatProduct(quatConjugate(cur_quat), quat);
    float3 new_AngVel = 2.0f * make_float3(delta_rot.x, delta_rot.y, delta_rot.z) / dt;

    //printf("newAngVel x:%f,y:%f,z:%f,w:%f\n", delta_rot.x, delta_rot.y, delta_rot.z,delta_rot.w);

    //角速度が一定以下なら，動いていないとして固定する
    if (length(new_AngVel) < ANGVEL_EPSILON) {
        new_AngVel = make_float3(0.f);
        quat = cur_quat;
        dquat[QUAT * id] = quat.x;
        dquat[QUAT * id + 1] = quat.y;
        dquat[QUAT * id + 2] = quat.z;
        dquat[QUAT * id + 3] = quat.w;
    }

    dangvel[DIM * id] = new_AngVel.x;
    dangvel[DIM * id + 1] = new_AngVel.y;
    dangvel[DIM * id + 2] = new_AngVel.z;

    dcurquat[QUAT * id] = quat.x;
    dcurquat[QUAT * id + 1] = quat.y;
    dcurquat[QUAT * id + 2] = quat.z;
    dcurquat[QUAT * id + 3] = quat.w;
}

//四元数の設定
__global__
void CxQuatSet(float* dcurquat, float* dquat, int* dfix, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n - 1) return;
    if (dfix[id + 1] == 1) return;

    float4 quat = make_float4(1, 0, 0, 1);
    quat = normalize(quat);

    dquat[QUAT * id] = dcurquat[QUAT * id] = quat.x;
    dquat[QUAT * id + 1] = dcurquat[QUAT * id + 1] = quat.y;
    dquat[QUAT * id + 2] = dcurquat[QUAT * id + 2] = quat.z;
    dquat[QUAT * id + 3] = dcurquat[QUAT * id + 3] = quat.w;
}

//基準密度を最初に計算した密度に設定
__global__
void CxRestDensSet(float* dpos,float* dRestDens,float* dvol,int n) {
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    float3 pos0 = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    //float m = params.mass;
    float a = params.aw;
    float rest_dens = params.rest_dens;

    // 粒子を中心として半径h内に含まれるグリッド(caclGridPos内で境界処理あり)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0 - make_float3(h));
    grid_pos1 = calcGridPos(pos0 + make_float3(h));

    // 周囲のグリッドセルも含めて近傍探索して密度計算
    float dens = 0.0f;
    for (int z = grid_pos0.z; z <= grid_pos1.z; ++z) {
        for (int y = grid_pos0.y; y <= grid_pos1.y; ++y) {
            for (int x = grid_pos0.x; x <= grid_pos1.x; ++x) {
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // グリッドハッシュ値

                // セル内のパーティクルのスタートインデックス
                uint startIndex = params.cell.dCellStart[ghash];
                if (startIndex != 0xffffffff) {	// セルが空でないかのチェック
                    // セル内のパーティクルで反復
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for (uint j = startIndex; j < endIndex; ++j) {
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 rij = pos0 - pos1;
                        float r = length(rij);
                        if (r <= h) {
                            // Poly6カーネルで密度を計算 (rho = Σ m Wij)
                            float q = h * h - r * r;
                            // 流体粒子はrest_dens*dvol[sj] = mとなるように設定してある
                            // 境界粒子は想体積と初期密度から複数層境界粒子があった場合の仮想質量Φ=ρ0*Vbを求めて使う 
                            float m = rest_dens * dvol[sj];
                            //海老沢SPH追加
                            //m = 0.01;
                            dens += m * a * q * q * q;
                        }
                    }
                }
            }
        }
    }


    // 計算した密度をグローバルメモリに書き込み
    uint sid = params.cell.dSortedIndex[id];
    //最低限の密度を設定
    //dRestDens[sid] = fmaxf(dens,rest_dens);
    dRestDens[sid] = dens;
}

__global__
void CxRestTotalDens(float* drestdens, float dens, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    drestdens[id] = dens;
}



//以下SagFreeの処理を移植--------------------------------------------------------------------------------------------------------------------------------
//LU分解を用いて、連立1次方程式を解く
//1なら成功,0なら失敗
__device__ __host__
int LUDecomp(float A[][4], int n) {
    if (n <= 0) return 0;

    for (int i = 0; i < n; i++) {
        //l_ijの計算(i>=j)
        for (int j = 0; j <= i; ++j) {
            float lu = A[i][j];
            for (int k = 0; k < j; k++) {
                lu -= A[i][k] * A[k][j];//l_ik*u_kj
            }
            A[i][j] = lu;
        }

        //u_ijの計算(i<j)
        for (int j = i + 1; j < n; ++j) {
            float lu = A[i][j];
            for (int k = 0; k < i; ++k) {
                lu -= A[i][k] * A[k][j];
            }
            A[i][j] = lu / A[i][i];
        }
    }

    return 1;
}

//A:LU分解された行列
//b:右辺ベクトル
//x:結果ベクトル
//n:行列の大きさ
__device__ __host__
int LUSolver(const float A[][4], const float b[], float x[], int n) {
    if (n <= 0) return 0;

    //前進代入
    //LY=bからYを計算
    for (int i = 0; i < n; ++i) {
        float bly = b[i];
        for (int j = 0; j < i; ++j) {
            bly -= A[i][j] * x[j];
        }
        //if (A[i][i] < 1.0e-6) printf("A trace error!\n");
        x[i] = bly / A[i][i];
    }

    //後退代入
    //UX=YからXを計算
    for (int i = n - 1; i >= 0; --i) {
        float yux = x[i];
        for (int j = i + 1; j < n; ++j) {
            yux -= A[i][j] * x[j];
        }
        x[i] = yux;
    }

    return 1;
}

//3次元ベクトルから四元数をもとめる
__device__ __host__
float4 quatFromDirector(float3 d3) {
    d3 = normalize(d3);
    float3 e3 = make_float3(0.f, 0.f, 1.0f);//z軸正を前提
    float3 w = cross(e3, d3);
    float4 q = make_float4(w, dot(e3, d3));
    q.w += length(q);
    return normalize(q);
}

//二つのベクトルから四元数を求める
__device__ __host__
float4 quatFromTwoVectors(float3 a, float3 b) {
    float3 v0 = normalize(a);
    float3 v1 = normalize(b);
    float c = dot(v0, v1);

    float3 axis = cross(v0, v1);
    float s = sqrt((1 + c) * 2);
    float3 vec = axis / s;
    float w = s * 0.5;

    float4 quat = make_float4(vec.x,vec.y,vec.z,w);
    return normalize(quat);
}

//グローバルフォースステップ
//ひとまず摩擦項を考えずに、単純に重力から求めることとする
__global__
void CxGlobalForceStep(float* dfss,float* dmass,int* dlast_ind,float3 gravity,int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    if (id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    int min;
    if (id == 0) min = 1;
    else min = dlast_ind[id - 1] + 2;
    int max = dlast_ind[id] - 1;
    for (int i = max; i > min - 1; i--) {//i=0の時は力0
        float mass = dmass[i];

        float3 prev_fss;
        if (i == max) prev_fss = make_float3(0.f);
        else prev_fss = make_float3(dfss[(i + 1) * DIM], dfss[(i + 1) * DIM + 1], dfss[(i + 1) * DIM + 2]);

        dfss[i * DIM] = -(mass * gravity.x - prev_fss.x);
        dfss[i * DIM + 1] = -(mass * gravity.y - prev_fss.y);
        dfss[i * DIM + 2] = -(mass * gravity.z - prev_fss.z);
    }
    
}

//ローカルフォースステップ
__global__
void CxLocalForceStep(float* dpos, float* dlen, float* dquat,float* dcurquat, float* dkss, float* dfss, int* dfix, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    if (id >= n - 1) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない) 最後の粒子はスキップ
    if (dfix[id] == 1 || dfix[id + 1] == 1) return;//弾性体毎に辺の数だけ行う、最初のエッジは固定して扱う

    float3 pos1 = make_float3(dpos[DIM * id], dpos[DIM * id + 1], dpos[DIM * id + 2]);
    float3 pos2 = make_float3(dpos[DIM * (id + 1)], dpos[DIM * (id + 1) + 1], dpos[DIM * (id + 1) + 2]);

    float3 fss = make_float3(dfss[DIM * id], dfss[DIM * id + 1], dfss[DIM * id + 2]);
    float fs_Len = Length(fss);

    float l0 = dlen[id];
    float tmp_ks = dkss[id];
    float4 quat = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);//四元数は(x,y,z,w)になるように受け取る

    float3 e3 = make_float3(0, 0, 1);
    float3 d3 = normalize(rotVecByQuat(e3, quat));

    //判別式を満たすかどうかを反復で行う----------------
    float delta = 10.f;
    float B = dot((pos1 - pos2), fss) / (tmp_ks)+1;
    float AC = dot(pos1 - pos2, pos1 - pos2) * fs_Len * fs_Len / (tmp_ks * tmp_ks);//Length2をdotで表現
    float discrim = B * B - 4 * AC;
    /*while (1) {
        if (discrim > 1.0e-2) break;
        tmp_ks += delta;
    }*/

    //---------------------------------------------------
    //判別式を満たさない場合、定式的にkssを求める--------
    if (discrim < 0) {
        //2次式より求めており、本来+-sqrt(...)であるが、kss>0より、+のみを考える
        tmp_ks = -dot(fss, pos1 - pos2) + 2 * fs_Len * Length(pos1 - pos2);
    }
    //---------------------------------------------------
    
    //伸び剛性を更新
    dkss[id] = tmp_ks;

    //二つの候補となる長さを求める(float型のままだと精度に影響があるため、double型に変更)
    double a = fs_Len * fs_Len;
    double b = -tmp_ks * (2.f * dot(fss, pos1 - pos2) + tmp_ks);
    double c = tmp_ks * tmp_ks * dot(pos1 - pos2, pos1 - pos2);//Length2をdotで表現

    double l1 = sqrt((-b + sqrt(abs(b * b - 4.f * a * c))) / (2.f * a));
    double l2 = sqrt((-b - sqrt(abs(b * b - 4.f * a * c))) / (2.f * a));

    //2つの候補のうち、より現在の長さに近いものを選ぶ
    float rest_length;
    if (abs(l1 - l0) > abs(l2 - l0) && abs(l2) > 1.0e-10) rest_length = l2;
    else if (abs(l1) > 1.0e-10) rest_length = l1;
    else {
        printf("Error Occured! in LocalTorqueStep!");
    }

    //姿勢の更新
    float3 new_ds = (rest_length * fss) / tmp_ks - (pos1 - pos2) / rest_length;

    //四元数をベクトルから求める(二つの手法があるが、ほとんど変わらないと推測)
    float4 new_qs = quatFromDirector(new_ds);
    //float4 new_qs = quatFromTwoVectors(e3, new_ds);

    //姿勢ベクトルの確認
    float3 d3_from_qs = rotVecByQuat(e3, new_qs);

    //力の確認用
    float3 new_Fss = (tmp_ks / rest_length) * ((pos1 - pos2) / rest_length + d3_from_qs);

    //基準長の更新
    dlen[id] = rest_length;
    //姿勢の更新
    dquat[id * QUAT] = dcurquat[id * QUAT] = new_qs.x;
    dquat[id * QUAT + 1] = dcurquat[id * QUAT + 1] = new_qs.y;
    dquat[id * QUAT + 2] = dcurquat[id * QUAT + 2] = new_qs.z;
    dquat[id * QUAT + 3] = dcurquat[id * QUAT + 3] = new_qs.w;
}

//伸び・せん断制約のトルクを求める
__device__ __host__
float4 StretchingShearTorque(float4 qs, float3 pos1, float3 pos2, float len, float kss) {
    float3 V = (pos1 - pos2) / len;

    float4 torqueSS;
    torqueSS.x = 4 * (qs.x * qs.x * qs.x + qs.x * qs.y * qs.y + qs.x * qs.z * qs.z + qs.x * qs.w * qs.w + V.z * qs.x - V.x * qs.z + V.y * qs.w);//x
    torqueSS.y = 4 * (qs.y * qs.y * qs.y + qs.y * qs.x * qs.x + qs.y * qs.z * qs.z + qs.y * qs.w * qs.w + V.z * qs.y - V.x * qs.w - V.y * qs.z);//y
    torqueSS.z = 4 * (qs.z * qs.z * qs.z + qs.z * qs.x * qs.x + qs.z * qs.y * qs.y + qs.z * qs.w * qs.w - V.z * qs.z - V.x * qs.x - V.y * qs.y);//z
    torqueSS.w = 4 * (qs.w * qs.w * qs.w + qs.w * qs.x * qs.x + qs.w * qs.y * qs.y + qs.w * qs.z * qs.z - V.z * qs.w - V.x * qs.y + V.y * qs.x);//w
    torqueSS = 1.f / 2.f * kss * torqueSS;

    return torqueSS;
}

//曲げ・ねじれ制約のトルクを求める(現在利用せず)
__device__ __host__
float4 BendTwistTorque(float4 q1, float4 q2, float4 Darboux, float kbt) {
    float4 torqueBT;
    float sum = q2.x * q2.x + q2.y * q2.y + q2.z * q2.z + q2.w * q2.w;
    torqueBT.x = sum * q1.x - q2.x * Darboux.w + q2.y * Darboux.z - q2.z * Darboux.y + q2.w * Darboux.x;
    torqueBT.y = sum * q1.y - q2.x * Darboux.z - q2.y * Darboux.w + q2.z * Darboux.x + q2.w * Darboux.y;
    torqueBT.z = sum * q1.z + q2.x * Darboux.y - q2.y * Darboux.x - q2.z * Darboux.w + q2.w * Darboux.z;
    torqueBT.w = sum * q1.w - q2.x * Darboux.x - q2.y * Darboux.y - q2.z * Darboux.z - q2.w * Darboux.w;

    return kbt * torqueBT;
}

//トルクから求める場合に利用(現在利用せず)
__device__ __host__
float4 SolveTorqueSolver(float4 q1, float4 q2, float4 torqueSS, float4 torqueBT, float kbt) {
    float4 rightForm = -(torqueSS + torqueBT);
    rightForm = rightForm / kbt;

    float4 New_Darboux;

    float sum = q2.x * q2.x + q2.y * q2.y + q2.z * q2.z + q2.w * q2.w;
    float b[4];
    b[0] = rightForm.x - sum * q1.x;
    b[1] = rightForm.y - sum * q1.y;
    b[2] = rightForm.z - sum * q1.z;
    b[3] = rightForm.w - sum * q1.w;

    float A[4][4];
    A[0][0] = q2.w;
    A[0][1] = -q2.z;
    A[0][2] = q2.y;
    A[0][3] = -q2.x;

    A[1][0] = q2.z;
    A[1][1] = q2.w;
    A[1][2] = -q2.x;
    A[1][3] = -q2.y;

    A[2][0] = -q2.y;
    A[2][1] = q2.x;
    A[2][2] = q2.w;
    A[2][3] = -q2.z;

    A[3][0] = -q2.x;
    A[3][1] = -q2.y;
    A[3][2] = -q2.z;
    A[3][3] = -q2.w;

    float x[4];
    int n1=LUDecomp(A, 4);
    //if (n1 == 0) printf("LUDecomp failure!!\n");
    int n2=LUSolver(A, b, x, 4);
    //if (n2 == 0)printf("LUSolver failure!!\n");

    New_Darboux = make_float4(x[0], x[1], x[2], x[3]);
    //printf("%d:New_Darboux x:%f,y:%f,z:%f,w:%\nf", New_Darboux.x, New_Darboux.y, New_Darboux.z, New_Darboux.w);
    return New_Darboux;
}

//四元数の逆元を返す
__device__ __host__
float4 QuatInverse(float4 quat) {
    return quatConjugate(quat) / dot(quat, quat);
}

//動画の方を参考に線形システム的に解く
//実際には、一番下のエッジから順番に解いていく
//一つ目のエッジは完全に固定するので、二つ目のエッジから考える
//隣接エッジの更新後の値を使う必要があるので、毛髪単位での並列化になる
__global__
void CxVideoGlobalTorqueStep(float* dpos, float* dquat, float* domega, float* dlen, float* dkss, float* dkbt, int* dfix, int* dlast_index, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    if (id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない) 最後の粒子はスキップ
    int last_edge = dlast_index[id] - 1;//最後の粒子-1が最後のエッジ
    float min;
    if (id == 0) min = 0;
    else min = dlast_index[id - 1]+1;//最初の粒子

    //一番下のエッジ(i==last_edge)の処理は他と違い、一つのみの基準ダルボーベクトルが関わる。一度きりなので、ループを使わずに書く---
    float3 init_pos1 = make_float3(dpos[last_edge * DIM], dpos[last_edge * DIM + 1], dpos[last_edge * DIM + 2]);
    float3 init_pos2 = make_float3(dpos[(last_edge + 1) * DIM], dpos[(last_edge + 1) * DIM + 1], dpos[(last_edge + 1) * DIM + 2]);

    float4 init_quat0= make_float4(dquat[(last_edge - 1) * QUAT], dquat[(last_edge - 1) * QUAT + 1], dquat[(last_edge - 1) * QUAT + 2], dquat[(last_edge - 1) * QUAT + 3]);
    float4 init_quat1= make_float4(dquat[last_edge * QUAT], dquat[last_edge * QUAT + 1], dquat[last_edge * QUAT + 2], dquat[last_edge * QUAT + 3]);

    float init_l0 = dlen[last_edge];
    float init_kss = dkss[last_edge];
    float init_kbt = dkbt[last_edge];

    float init_K = 2 * init_kbt / init_l0;
    float4 init_kq_inv = QuatInverse(init_K * init_quat1);

    float4 init_torqueSS = StretchingShearTorque(init_quat1, init_pos1, init_pos2, init_l0, init_kss);
    init_torqueSS = quatProduct(init_torqueSS, init_kq_inv);

    float4 init_Cur_Omega_Prev = quatProduct(quatConjugate(init_quat0), init_quat1);
    float4 init_Rest_Omega_Prev = init_Cur_Omega_Prev - init_torqueSS;
    
    //init_Rest_Omega_Prev = quatConjugate(init_Rest_Omega_Prev);
    //ダルボーベクトルの方向を曲げねじれ制約のように求める必要があるが、ほぼ1であると推測できるため、1で扱う(s*Omegaをそのままダルボーベクトルの配列に入れる)
    domega[(last_edge - 1) * QUAT] = init_Rest_Omega_Prev.x;
    domega[(last_edge - 1) * QUAT + 1] = init_Rest_Omega_Prev.y;
    domega[(last_edge - 1) * QUAT + 2] = init_Rest_Omega_Prev.z;
    domega[(last_edge - 1) * QUAT + 3] = init_Rest_Omega_Prev.w;
    //-----------------------------------------------------------------------------------------------------------------------------

    //他のエッジの処理
    for (int i = last_edge-1; i > min; i--) {//iのひとつ前の基準ダルボーベクトルを求める
        float3 pos1 = make_float3(dpos[i * DIM], dpos[i * DIM + 1], dpos[i * DIM + 2]);
        float3 pos2 = make_float3(dpos[(i + 1) * DIM], dpos[(i + 1) * DIM + 1], dpos[(i + 1) * DIM + 2]);

        float4 quat0 = make_float4(dquat[(i - 1) * QUAT], dquat[(i - 1) * QUAT + 1], dquat[(i - 1) * QUAT + 2], dquat[(i - 1) * QUAT + 3]);
        float4 quat1 = make_float4(dquat[i * QUAT], dquat[i * QUAT + 1], dquat[i * QUAT + 2], dquat[i * QUAT + 3]);
        float4 quat2 = make_float4(dquat[(i + 1) * QUAT], dquat[(i + 1) * QUAT + 1], dquat[(i + 1) * QUAT + 2], dquat[(i + 1) * QUAT + 3]);

        float l0 = dlen[i];
        float kss = dkss[i];
        float kbt = dkbt[i];

        //伸び・せん断制約のトルクをこれで割る
        float K = 2 * kbt / l0;
        float4 kq_inv = QuatInverse(K * quat1);

        //伸び・せん断制約のトルクを求める
        float4 torqueSS = StretchingShearTorque(quat1, pos1, pos2, l0, kss);
        torqueSS = quatProduct(torqueSS, kq_inv);

        float4 Cur_Omega_Prev = quatProduct(quatConjugate(quat0), quat1);//現在のエッジとひとつ前のエッジのダルボーベクトル
        float4 Cur_Omega_Next = quatConjugate(quatProduct(quatConjugate(quat1), quat2));//現在のエッジと一つ際のエッジのダルボーベクトル
        
        float4 Rest_Omega_Next = make_float4(domega[i * QUAT], domega[i * QUAT + 1], domega[i * QUAT + 2], domega[i * QUAT + 3]);
        //最終的に求める基準ダルボーベクトル
        float4 Rest_Omega_Prev = Cur_Omega_Next + Cur_Omega_Prev - Rest_Omega_Next - torqueSS;//ひとつ前のエッジとの間の基準ダルボーベクトル
        Rest_Omega_Prev = quatConjugate(Rest_Omega_Prev);

        //ダルボーベクトルの方向を曲げねじれ制約のように求める必要があるが、ほぼ1であると推測できるため、1で扱う(s*Omegaをそのままダルボーベクトルの配列に入れる)
        domega[(i - 1) * QUAT] = Rest_Omega_Prev.x;
        domega[(i - 1) * QUAT + 1] = Rest_Omega_Prev.y;
        domega[(i - 1) * QUAT + 2] = Rest_Omega_Prev.z;
        domega[(i - 1) * QUAT + 3] = Rest_Omega_Prev.w;
    }
}

//LocalTorqueStepを解くのに使う
__device__ __host__
float4 solveInverseRot(float4 quat, float4 tar, float& bendK) {
    const float SAFETY_FACTOR = min(quat.w, 0.002f);//minが扱えるか、0.00002f,0.002f

    tar -= dot(tar, quat) * quat;
    float4 Omega = -tar / bendK;

    if (dot(Omega,Omega) > SAFETY_FACTOR * SAFETY_FACTOR) {//Length2をdotに置き換え
        bendK = Length(tar) / SAFETY_FACTOR;
        //printf("bendK(Updated):%f\n", bendK);
        Omega = -tar / bendK;
    }

    float4 d = sqrt(1 - dot(Omega, Omega)) * quat;//Length2をdotに置き換え
    Omega += d;

    //printf("Last_Omega x:%f,y:%f,z:%f,w:%f\n", Omega.x, Omega.y, Omega.z, Omega.w);

    return Omega;
}

//ローカルトルクステップ
__global__
void CxLocalTorqueStep(float* dquat, float* domega, float* dlen,float* dkbt, int* dfix, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    if (id >= n - 2) return; // 粒子数を超えるスレッドIDのチェック
    if (dfix[id + 1] == 1 || dfix[id + 2] == 1) return;

    float4 quat1 = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);
    float4 quat2 = make_float4(dquat[QUAT * (id + 1)], dquat[QUAT * (id + 1) + 1], dquat[QUAT * (id + 1) + 2], dquat[QUAT * (id + 1) + 3]);
    float4 qs = quatProduct(quatConjugate(quat1), quat2);

    float4 omega = make_float4(domega[QUAT * id], domega[QUAT * id + 1], domega[QUAT * id + 2], domega[QUAT * id + 3]);

    float length = dlen[id];
    float tmp_kbt = dkbt[id];

    float4 last_omega = solveInverseRot(qs, omega, tmp_kbt);

    //曲げねじれ制約の剛性の更新
    dkbt[id] = tmp_kbt;
    //基準ダルボーベクトルの更新
    domega[QUAT * id] = last_omega.x;
    domega[QUAT * id + 1] = last_omega.y;
    domega[QUAT * id + 2] = last_omega.z;
    domega[QUAT * id + 3] = last_omega.w;
}

//密度制約
//lambdaの計算
__global__
void CxPbfLambda(float* ddens,float* drestdens,float* dpbf_lambda,float* dvol,int n) {//dpdf_lambdaは未定義
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    float3 pos_i = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    //float m = params.mass;
    float a = params.aw;
    float rest_dens = params.rest_dens;
    //インデックスの計算
    uint sid = params.cell.dSortedIndex[id];
    //海老沢SPH追加
    rest_dens = drestdens[sid];
    float dens = ddens[sid];

    // 粒子を中心として半径h内に含まれるグリッド(caclGridPos内で境界処理あり)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos_i - make_float3(h));
    grid_pos1 = calcGridPos(pos_i + make_float3(h));
    //---------------------------------------------

    //密度から制約の計算に用いるλを求める
    float C = dens / rest_dens - 1.0f;
    if (C > 0.f) {
        float lambda_denom_i = 0.f;
        float3 grad_i_C = make_float3(0.f);
        for (int z = grid_pos0.z; z <= grid_pos1.z; ++z) {
            for (int y = grid_pos0.y; y <= grid_pos1.y; ++y) {
                for (int x = grid_pos0.x; x <= grid_pos1.x; ++x) {
                    int3 ngrid = make_int3(x, y, z);
                    uint ghash = calcGridHash(ngrid);   // グリッドハッシュ値

                    // セル内のパーティクルのスタートインデックス
                    uint startIndex = params.cell.dCellStart[ghash];
                    if (startIndex != 0xffffffff) {	// セルが空でないかのチェック
                        // セル内のパーティクルで反復
                        uint endIndex = params.cell.dCellEnd[ghash];
                        for (uint j = startIndex; j < endIndex; ++j) {
                            uint sj = params.cell.dSortedIndex[j];
                            float3 pos_j = params.cell.dSortedPos[j];
                            float rest_dens_j = drestdens[sj];

                            //有効半径内に存在するかを確認
                            float3 r_ij = pos_i - pos_j;
                            float r = length(r_ij);

                            if (r <= 1.0e-3f) continue;
                            if (r < h) {
                                float q = h - r;
                                float m = rest_dens * dvol[sj];//初期密度と体積から計算
                                //Grad_J_Cを求める
                                float3 grad_j_C = -m / rest_dens_j * (-params.ag * q * q * r_ij / r);//params.agにカーネルの勾配定数(spikyカーネル)が格納されている.
                                lambda_denom_i += dot(grad_j_C, grad_j_C);
                                grad_i_C += grad_j_C;
                            }
                        }
                    }
                }
            }
        }
        lambda_denom_i += dot(grad_i_C, grad_i_C);
        //λに値を格納(float型)
        dpbf_lambda[sid] = -C / (lambda_denom_i + 1.0e-6f);//εを1.0e-6fで定義
    }

    else {
        dpbf_lambda[sid] = 0.f;
    }
}

__global__
void CxPbfConstraint(float*dpos,float* drestdens,float* dpbf_lambda,float* dvol,int n) {
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)

    float3 pos_i = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    //インデックスの計算
    uint sid = params.cell.dSortedIndex[id];
    //pbfの計算に必要なλを持ってくる
    float pbf_lambda_i = dpbf_lambda[sid];

    // 粒子を中心として半径h内に含まれるグリッド(caclGridPos内で境界処理あり)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos_i - make_float3(h));
    grid_pos1 = calcGridPos(pos_i + make_float3(h));

    // 周囲のグリッドセルも含めて近傍探索して密度計算
    float3 delta_pos_i = make_float3(0.f);
    for (int z = grid_pos0.z; z <= grid_pos1.z; ++z) {
        for (int y = grid_pos0.y; y <= grid_pos1.y; ++y) {
            for (int x = grid_pos0.x; x <= grid_pos1.x; ++x) {
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // グリッドハッシュ値
                // セル内のパーティクルのスタートインデックス
                uint startIndex = params.cell.dCellStart[ghash];
                if (startIndex != 0xffffffff) {	// セルが空でないかのチェック
                    // セル内のパーティクルで反復
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for (uint j = startIndex; j < endIndex; ++j) {
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos_j = params.cell.dSortedPos[j];
                        float restdens_j = drestdens[sj];
                        float pbf_lambda_j = dpbf_lambda[sj];

                        float3 r_ij = pos_i - pos_j;
                        float r = length(r_ij);
                        if (r <= 1.0e-3) continue;
                        if (r < h) {
                            float q;
                            q = h - r;
                            float m = restdens_j * dvol[sj];//初期密度と体積から計算
                            delta_pos_i += m / restdens_j * (pbf_lambda_i + pbf_lambda_j) * (params.ag * q * q * r_ij / r);
                        }
                    }
                }
            }
        }
    }

    //結果の格納
    dpos[DIM * sid] += delta_pos_i.x;
    dpos[DIM * sid + 1] += delta_pos_i.y;
    dpos[DIM * sid + 2] += delta_pos_i.z;
}

//摩擦制約の実装
//摩擦制約はとりあえず，速度更新の際に，一度のみ行うように設定する．
//反復する場合にはXPBDにしないと，反復に依存して修正量が変化すると考えられる
//周辺粒子が大事なので，Sortされた位置を利用する
__global__
void CxFrictionConstraint(float* dpos, float* dcurpos, int* dfix, int n) {
    // グリッド,ブロック内のスレッド位置を粒子インデックスとする
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // 粒子数を超えるスレッドIDのチェック(余りが出ないようにブロック数などが設定できるなら必要ない)
    if (dfix[id] == 1) return;//固定点であれば，スキップ

    float3 pos_i = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    //インデックスの計算
    uint sid = params.cell.dSortedIndex[id];

    // 粒子を中心として半径h内に含まれるグリッド(caclGridPos内で境界処理あり)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos_i - make_float3(h));
    grid_pos1 = calcGridPos(pos_i + make_float3(h));

    for (int z = grid_pos0.z; z <= grid_pos1.z; ++z) {
        for (int y = grid_pos0.y; y <= grid_pos1.y; ++y) {
            for (int x = grid_pos0.x; x <= grid_pos1.x; ++x) {
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // グリッドハッシュ値
                // セル内のパーティクルのスタートインデックス
                uint startIndex = params.cell.dCellStart[ghash];
                if (startIndex != 0xffffffff) {	// セルが空でないかのチェック
                    // セル内のパーティクルで反復
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for (uint j = startIndex; j < endIndex; ++j) {
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos_j = params.cell.dSortedPos[j];
                        float3 r_ij = pos_i - pos_j;
                        float r = length(r_ij);
                        if (r <= 1.0e-3) continue;
                        if (r < h) {
                            //ここに摩擦制約を書く

                        }
                    }
                }
            }
        }
    }

    //位置修正による更新
    dpos[DIM * sid];
    dpos[DIM * sid + 1];
    dpos[DIM * sid + 2];
}
