/*!
  @file sph.h

  @brief SPHによる流体シミュレーション

  @author Makoto Fujisawa
  @date 2023-02
*/

#ifndef _SPH_H_
#define _SPH_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "utils.h"

#include "nnsearch.h"	// グリッド分割による近傍探索
#include "pdist.h"		// 粒子配置生成用
#include "kernel.h"		// カーネル関数の定義
#include "solid.h"		// 固体境界の定義

#include "cuda_utils.h"

#include <deque>


//-----------------------------------------------------------------------------
// MARK:SPHクラスの宣言
//  - Miles Macklin and Matthias Muller, "Position Based Fluids", Proc. SIGGRAPH 2013, 2013. 
//  - http://blog.mmacklin.com/publications/
//-----------------------------------------------------------------------------
class SPH
{
public:
	enum ParticleArray
	{
		P_POSITION = 0,
		P_VELOCITY,
		P_ACC,//海老沢追加
		P_FORCE,
		P_VORTICITY,
		P_DENSITY,
		P_PRESSURE,
		P_VOLUME, 
		P_ATTRIBUTE,
		//海老沢追加-----------------------
		P_MASS,
		E_LENGTH,
		E_KSS,
		E_KBT,
		E_QUAT,
		E_DARBOUX,
		E_LAMBDA_SS,
		E_LAMBDA_BT,
		P_FIX,
		P_CURPOS,
		E_TANG,
		E_CURQUAT,
		E_ANGVEL,
		P_RESTDENS,
		E_FSS,
		P_LAST_IND,
		P_PBF_LAMBDA,
		//--------------------------------
		P_END,
	};
	enum ParticleColorArray
	{
		C_CONSTANT = 0,
		C_DENSITY,
		C_PRESSURE,
		C_VORTICITY,
		C_RAMP,
		C_END,
	};

	static const int MAX_STEPS = 50000;	//!< シミュレーション最大ステップ数

protected:
	bool m_initialized;				//!< 初期化済みかどうかのフラグ

	uint m_num_particles;			//!< 現在の粒子数
	uint m_max_particles;			//!< 最大粒子数

	glm::vec3 m_envmin;				//!< 環境のAABB最小座標
	glm::vec3 m_envmax;				//!< 環境のAABB最大座標

	int m_colortype;				//!< 粒子色分け描画設定(カラーVBO用)

	vector<rxPdist> m_pdist;		//!< シーンに追加する粒子情報

	// 境界・固体
	vector<rxSolid*> m_solids;		//!< 固体物体

	// デバイスメモリ(VBO)
	uint m_vbo_pos;					//!< 粒子座標VBO
	uint m_vbo_col;					//!< カラーVBO
	//海老沢追加
	uint m_vbo_tang;//法線VBO

	//! OpenGL(のVBO,PBO)-CUDA間のマッピングを扱うためのハンドル
	cudaGraphicsResource* m_cgr_pos;
	cudaGraphicsResource* m_cgr_col;
	//海老沢追加
	cudaGraphicsResource* m_cgr_tang;

	// デバイスメモリ
	float* m_d_vel;					//!< 粒子速度
	int* m_d_attr;					//!< 粒子属性(流体粒子:0，境界粒子:1)
	float* m_d_vol;					//!< 粒子体積([Akinci et al.,SIG2012]の境界粒子計算用(Φ))

	float* m_d_acc;					//!< 粒子にかかる力(加速度)(ベクトル)
	float* m_d_dens;				//!< 粒子密度(スカラー)
	float* m_d_pres;				//!< 粒子圧力(スカラー)
	float* m_d_vort;				//!< 粒子の渦度(ベクトル)

	//海老沢追加
	float* m_d_mass;//粒子質量
	float* m_d_rest_length;//基準長
	float* m_d_kss;//エッジの伸び剛性
	float* m_d_kbt;//エッジの曲げ剛性
	float* m_d_quat;//四元数
	float* m_d_rest_darboux;//エッジ間の基準ダルボーベクトル
	float* m_d_lambda_ss;//XPBDのSSのλ
	float* m_d_lambda_bt;//XPBDのBTのλ
	int* m_d_fix;//固定点かどうか(1ならば固定点)
	float* m_d_curpos;//時間積分に用いる更新前の位置
	float* m_d_curquat;//時間積分に用いる更新前の位置
	float* m_d_angvel;//回転速度を格納
	float* m_d_rest_density;//基準となる密度を格納

	float* m_d_fss;//エッジにかかる力を格納
	int* m_d_last_index;//毛髪ごとの最後の点のインデックスを格納

	float* m_d_pbf_lambda;//密度制約に用いるλを格納

public:
	SceneParameter m_params;		//!< シミュレーションパラメータ(GPUとのやりとりにも使う)

	bool m_use_bparticles;			//!< 境界粒子使用のON/OFFフラグ
	int m_method_visc;				//!< 粘性項計算方法(0:粘性項なし,1:WCSPHの論文の方法,2:XSPH人工粘性)

	uint m_num_bparticles;			//!< 境界粒子の数
	uint m_offset_bparticles;		//!< 境界粒子のスタート地点

	//海老沢追加---------------------------------------------
	int m_numElastic;  //弾性体の数
	// 頂点配列オブジェクト
	uint m_vao_pos;//計算点配列VAO

	//SPHの際に与える風の強さ
	float3 m_wind_power;
	//風の力による変形を適用するかどうか
	bool m_wind_flag;
	//衝突を行う球の設定
	float3 m_center;
	float m_rad;
	//重みに関するフラグを付ける
	bool m_example_flag;
	//-------------------------------------------------------
public:
	//! コンストラクタ
	SPH();

	//! デストラクタ
	virtual ~SPH();

	// シミュレーション空間情報の取得
	glm::vec3 GetMax(void) const { return m_envmax; }
	glm::vec3 GetMin(void) const { return m_envmin; }
	glm::vec3 GetDim(void) const { return m_envmax - m_envmin; }
	glm::vec3 GetCen(void) const { return 0.5f * (m_envmax + m_envmin); }
	glm::vec3 GetDimB(void) const { return m_solids[0]->GetMax() - m_solids[0]->GetMin(); }
	glm::vec3 GetCenB(void) const { return 0.5f * (m_solids[0]->GetMax() + m_solids[0]->GetMin()); }
	//NNGrid* GetNNGrid(void){ return m_nn; }

	//海老沢追加
	void ChangeNumParticles(int n);

	// 粒子情報
	int	GetNumParticles() const { return m_num_particles; }
	int	GetMaxParticles() const { return m_max_particles; }
	int GetNumBoundaryParticles() const { return m_num_bparticles; }
	int	GetNumFluidParticles() const { return m_num_particles - m_num_bparticles; }
	float GetSpacing() { return m_params.particle_radius * 2.0; }	//! 粒子配置間隔

	// シミュレーション設定
	void SetGravity(float x) { m_params.gravity = make_float3(0.0, x, 0.0); }	//!< 重力

	// 粒子VBO
	unsigned int GetCurrentReadBuffer(void) const { return m_vbo_pos; }
	unsigned int GetColorBuffer(void) const { return m_vbo_col; }
	//海老沢追加 接線を保存
	unsigned int GetTangBuffer(void) const { return m_vbo_tang; }

	// 描画用カラー設定
	void SetColorType(int type) { m_colortype = type; }
	int  GetColorType(void) const { return m_colortype; }

public:
	// シミュレーションステップ
	bool Update(float dt, int step = 0);

	// シーンの設定
	void SetBoxObstacle(glm::vec3 cen, glm::vec3 ext, glm::vec3 ang, glm::vec3 vel, int flg);
	void SetSphereObstacle(glm::vec3 cen, float rad, glm::vec3 vel, int flg);

	// ホスト<->VBO間転送
	void GetArrayFromDevice(int type, void* hdata, int num = -1);
	void SetArrayToDevice(int type, const void* data, int start, int count);
	void SetColorVBO(int type = -1, int picked = -1);

	//海老沢追加 XPBDに必要なパラメータをGPUに渡す
	void SetXPBD_Params(int max_particles, const void* mass, const void* length, const void* kss, const void* kbt, const void* quat, const void* darboux,const void* fix,const void* last_index);
	//void SetXPBD_Params(int max_particles, const void* mass, const void* quat);
	//海老沢追加 SagFree処理を行う
	void SagFree(glm::vec3 gravity);

	// 陰関数値計算
	float GetImplicit(float x, float y, float z);
	static float GetImplicit_s(float x, float y, float z, void* p) { return (float)(((SPH*)p)->GetImplicit(x, y, z)); }
	void CalImplicitField(int n[3], glm::vec3 minp, glm::vec3 d, float *hF);

	// SPH情報出力
	void OutputSetting(string fn);

	// 描画関数
	void DrawCells(glm::vec3 col);
	void DrawObstacles(int drw);

public:
	// シミュレーション初期化
	void Initialize(const SceneParameter& env,int num_elastic);
	void Allocate(int max_particles);
	void Finalize(void);

	// 分割セルに粒子を格納
	void SetParticlesToCell(float* dpos, float* dvel, int n, float h);

	// 境界粒子配置
	int SetBoundaryParticles(bool set = true);

	// 粒子設定
	void Reset(void);
	bool Set(const vector<glm::vec3>& ppos, const vector<glm::vec3>& pvel, const vector<glm::vec3>& etang, int particles);
	int  Add(rxPdist pd);

	// 境界の大きさ変更
	void ChangeBoundary(const glm::vec3& sl)
	{
		if (!m_solids.empty()) {
			static glm::vec3 minp = m_solids[0]->GetMin(), maxp = m_solids[0]->GetMax();
			m_solids[0]->SetMin(maxp - (maxp - minp) * sl);
			cout << "bnd : " << glm::to_string(0.5f * (maxp - minp) * sl) << endl;
		}
	}
	//海老沢追加
	//球を動かす
	void MoveSphere(float3 vel, float dt);


	int OutputParticles(string fn);
	int InputParticles(string fn);


protected:
	// rest densityの計算
	float calRestDensity(float h);

	// 近傍探索用セルの設定
	void setupNNCell(Cell& c, glm::vec3 vMin, glm::vec3 vMax, double h, int n);

	// 粒子追加
	int add(rxPdist pd, int start = -1, bool h2d = true);

	//-----------------------------------------------------------------------------
	// VBO,PBO関係
	//-----------------------------------------------------------------------------
	/*!
	 * VBO,PBOの作成
	 * @param[out] name 名前(VBO,PBOのID)
	 * @param[in] size バッファサイズ(バイト)
	 */
	void createVBO(GLuint& name, int size)
	{
		// VBO作成とバインド
		glGenBuffers(1, &name);
		glBindBuffer(GL_ARRAY_BUFFER, name);

		// バッファの作成と初期化
		glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_DRAW);

		// アンバインド
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}

	/*!
	 * VBO,PBOの削除
	 * @param[inout] name 名前(VBO,PBOのID)
	 */
	void deleteVBO(GLuint& name)
	{
		glBindBuffer(GL_ARRAY_BUFFER, name);
		glDeleteBuffers(1, &name);
		name = 0;
	}

	//海老沢追加
	void createVAO(GLuint& name) {
		//VAO作成とバインド
		glGenVertexArrays(1, &name);
		glBindVertexArray(name);
	}

	void BindAttrib(void) {
		glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
		glEnableVertexAttribArray(0);
	}

};






#endif	// _SPH_H_