/*!
  @file sph.cpp
	
  @brief SPH法の実装
 
  @author Makoto Fujisawa
  @date   2023-02
*/

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "sph.h"

//#include "rx_pcube.h"
#include "gldraw.h"

#include "rx_timer.h"	// 計算時間計測

#include "sph.cuh"	// CUDAによるSPH計算


//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
#ifndef P_USE_TIMER
#define P_USE_TIMER
#endif


/*!
 * デバッグ用:デバイスメモリの中身をファイルに出力
 * @param[in] name ファイル名
 * @param[in] ddata デバイスメモリの配列ポインタ
 */
template<typename T = float>
static inline void OutputDeviceMemory(string name, T* ddata, int n)
{
	vector<T> hdata(n);
	// CPU→GPUデータ転送
	CuCopyArrayFromDevice(&hdata[0], ddata, 0, 0, n * sizeof(T));
	// ファイル出力
	ofstream fo;
	fo.open(name.c_str());
	for (int i = 0; i < n; ++i) fo << hdata[i] << ", ";
	fo << endl;
	fo.close();


}


//-----------------------------------------------------------------------------
// SPHクラスの実装
//-----------------------------------------------------------------------------
/*!
 * コンストラクタ
 * @param[in] use_opengl VBO使用フラグ
 */
SPH::SPH() :
	m_initialized(false),
	m_d_vel(0),
	m_d_attr(0),
	m_d_vol(0),
	m_d_acc(0),
	m_d_dens(0),
	m_d_pres(0),
	m_d_vort(0),
	//海老沢追加
	m_d_mass(0),
	m_d_rest_length(0),
	m_d_kss(0),
	m_d_kbt(0),
	m_d_quat(0),
	m_d_rest_darboux(0),
	m_d_lambda_ss(0),
	m_d_lambda_bt(0),
	m_d_fix(0),
	m_d_curpos(0),
	m_d_curquat(0),
	m_d_angvel(0),
	m_d_rest_density(0),
	m_d_fss(0),
	m_d_last_index(0),
	m_d_pbf_lambda(0)
{
	m_num_particles = 0;
	m_num_bparticles = 0;

	m_use_bparticles = true;
	m_method_visc = 2;

	m_colortype = C_CONSTANT;
}

/*!
 * デストラクタ
 */
SPH::~SPH()
{
	Finalize();
}


/*!
 * シミュレーションの初期化
 * @param[in] max_particles 最大粒子数
 * @param[in] boundary_ext 境界の大きさ(各辺の長さの1/2)
 * @param[in] dens 初期密度
 * @param[in] mass 粒子の質量
 * @param[in] kernel_particle 有効半径h以内の粒子数
 */
void SPH::Initialize(const SceneParameter &env,int num_elastic)
{
	cout << "[SPH::Initialize]" << endl;
	CuInit();

	m_params = env;
	m_params.rest_dens = env.dens;
	//海老沢追加--------------------------------------------------------------------------
	m_params.rest_dens = 8.0e4;//9.8e3,8.0e3
	m_params.mass = 0.1;
	cout << "max particles " << m_params.max_particles << endl;

	//風の設定(GUIで一つの値を適用したり、解除したりできる)
	m_wind_power = make_float3(0.f, 0.f, 0.f);
	//衝突をする球との判定
	m_center = make_float3(-5.0, 0.0, 0.0);
	m_rad = 0.35;
	//ExampleRodの場合，重みや処理の一部を変更する必要があるため，フラグで管理する
	m_example_flag = false;
	//------------------------------------------------------------------------------------ 
	// 有効半径などの計算
	float volume = m_params.max_particles*m_params.mass/ m_params.rest_dens;
	//float volume = 1.0e-4;//海老沢追加

	cout << "max_particles:" << m_params.max_particles << " kernel_particles:" << m_params.kernel_particles << " volume:" << volume << endl;

	m_params.effective_radius = pow(((3.0*m_params.kernel_particles*volume)/(4.0* m_params.max_particles*RX_PI)), 1.0/3.0);
	m_params.particle_radius = pow((RX_PI/(6.0* m_params.kernel_particles)), 1.0/3.0)* m_params.effective_radius;
	m_params.kernel_radius = m_params.effective_radius;
	
	float h = m_params.effective_radius;
	float r = m_params.particle_radius;
	m_params.h = m_params.effective_radius;

	cout << "effective_radius:" << h << endl;

	// カーネル関数の定数
	m_params.aw = KernelCoefPoly6(h, 3, 1);
	m_params.ag = KernelCoefSpiky(h, 3, 2);
	m_params.al = KernelCoefVisc(h, 3, 3);

	// 初期密度の再計算
	m_params.rest_dens = calRestDensity(h);

	cout << "rest_dens " << m_params.rest_dens << endl;

	//m_paramsのパラメータは中から変更可能

	cout << "particle : " << endl;
	cout << " n_max = " << m_params.max_particles << endl;
	cout << " h = " << m_params.effective_radius << endl;
	cout << " r = " << m_params.particle_radius << endl;
	cout << " dens = " << m_params.rest_dens << endl;
	cout << " mass = " << m_params.mass << endl;
	cout << " kernel_particles = " << m_params.kernel_particles << endl;
	cout << " volume = " << volume << endl << endl;
	cout << " viscosity = " << m_params.viscosity << endl;
	cout << " coefs for kernel = " << m_params.aw << ", " << m_params.ag << ", " << m_params.al << endl;

	// シミュレーション環境の大きさ
	glm::vec3 bcen = glm::vec3(env.boundary_cen.x, env.boundary_cen.y, env.boundary_cen.z);
	glm::vec3 bext = glm::vec3(env.boundary_ext.x, env.boundary_ext.y, env.boundary_ext.z);
	m_envmin = bcen-bext;
	m_envmax = bcen+bext;

	cout << "simlation range : " << glm::to_string(m_envmin) << " - " << glm::to_string(m_envmax) << endl;

	// 粒子の衝突判定に用いる境界の大きさ．描画のために半径分小さくしておく
	m_params.boundary_min = make_float3(m_envmin[0]+r, m_envmin[1]+r, m_envmin[2]+r);
	m_params.boundary_max = make_float3(m_envmax[0]-r, m_envmax[1]-r, m_envmax[2]-r);
	
	// 境界設定
	m_solids.push_back(new rxSolidBox(m_envmin, m_envmax, -1));

	m_num_particles = m_params.max_particles;

	m_numElastic = num_elastic;//弾性体の数
	Allocate(m_params.max_particles);
}

/*!
 * メモリの確保
 *  - 最大粒子数で確保
 * @param[in] max_particles 最大粒子数
 */
void SPH::Allocate(int max_particles)
{
	//m_num_particles = max_particles;
	m_max_particles = max_particles;

	unsigned int vec_size = m_max_particles*DIM;
	unsigned int float_size = m_max_particles;
	unsigned int mem_vec_size = sizeof(float)*vec_size;
	unsigned int mem_float_size = sizeof(float)*float_size;

	//海老沢追加
	unsigned int quat_size = m_max_particles * QUAT;
	unsigned int mem_quat_size = sizeof(float) * quat_size;

	// 境界粒子設定
	//SetBoundaryParticles(true);

	// GPUメモリ確保
	createVBO(m_vbo_pos, mem_vec_size);
	createVBO(m_vbo_col, sizeof(float)*m_max_particles*4);	// 色については表示非表示を切り替えるためにαを使うので4チャンネルにする
	//海老沢追加
	createVBO(m_vbo_tang, mem_vec_size);//接線を法線として保存する

	// VBOバッファをCUDAに登録
	CuRegisterGLBufferObject(m_vbo_pos, &m_cgr_pos);
	CuRegisterGLBufferObject(m_vbo_col, &m_cgr_col);
	//海老沢追加
	CuRegisterGLBufferObject(m_vbo_tang, &m_cgr_tang);

	// デバイスメモリ確保
	CuAllocateArray((void**)&m_d_vel, mem_vec_size);
	CuAllocateArray((void**)&m_d_acc, mem_vec_size);
	CuAllocateArray((void**)&m_d_vort, mem_vec_size);
	CuAllocateArray((void**)&m_d_dens, m_max_particles * sizeof(float));
	CuAllocateArray((void**)&m_d_pres, m_max_particles * sizeof(float));
	CuAllocateArray((void**)&m_d_attr, m_max_particles * sizeof(int));
	CuAllocateArray((void**)&m_d_vol, m_max_particles * sizeof(float));

	//海老沢追加
	CuAllocateArray((void**)&m_d_mass, mem_float_size);						//質量
	CuAllocateArray((void**)&m_d_rest_length, mem_float_size);				//基準長
	CuAllocateArray((void**)&m_d_kss, mem_float_size);						//伸び剛性
	CuAllocateArray((void**)&m_d_kbt, mem_float_size);						//曲げ剛性
	CuAllocateArray((void**)&m_d_quat, mem_quat_size);						//四元数
	CuAllocateArray((void**)&m_d_rest_darboux, mem_quat_size);				//基準ダルボーベクトル
	CuAllocateArray((void**)&m_d_lambda_ss, mem_vec_size);					//XPBDのSSのλ
	CuAllocateArray((void**)&m_d_lambda_bt, mem_quat_size);					//XPBDのBTのλ
	CuAllocateArray((void**)&m_d_fix, m_max_particles * sizeof(int));		//固定点かどうか
	CuAllocateArray((void**)&m_d_curpos, mem_vec_size);						//前ステップの位置
	CuAllocateArray((void**)&m_d_curquat, mem_quat_size);					//前ステップの位置
	CuAllocateArray((void**)&m_d_angvel, mem_vec_size);						//回転速度
	CuAllocateArray((void**)&m_d_rest_density, mem_float_size);				//個々に初期密度を設定
	CuAllocateArray((void**)&m_d_fss, mem_vec_size);						//エッジにかかる力
	CuAllocateArray((void**)&m_d_last_index, m_numElastic * sizeof(int));	//毛髪の最後のインデックスを格納
	CuAllocateArray((void**)&m_d_pbf_lambda, mem_float_size);


	// 近傍粒子探索用セルのセットアップ
	glm::vec3 nnmin = m_envmin - glm::vec3(4.0 * m_params.particle_radius);
	glm::vec3 nnmax = m_envmax + glm::vec3(4.0 * m_params.particle_radius);
	setupNNCell(m_params.cell, nnmin, nnmax, m_params.effective_radius, m_max_particles);

	CuSetParameters(&m_params);

	SetColorVBO(SPH::C_CONSTANT);

	m_initialized = true;

	//一律の初期密度
	//CuRestTotalDens(m_d_rest_density,m_params.rest_dens, m_num_particles);
}

//海老沢追加
void SPH::SagFree(void) {
	float* pos = (float*)CuMapGLBufferObject(&m_cgr_pos);

	CuGlobalForceStep(pos,m_d_fss, m_d_mass, m_d_last_index, make_float3(0, -9.81, 0), m_d_dens, m_d_rest_density, m_d_vol, m_numElastic);

	CuLocalForceStep(pos, m_d_rest_length, m_d_quat,m_d_curquat, m_d_kss, m_d_fss, m_d_fix, m_num_particles);
	//トルクの計算
	//cout << "before Global Torque Step------------------------------------------------------" << endl;
	//CuCalcTorque(pos, m_d_mass, m_d_quat, m_d_fss, m_d_rest_length, m_d_kss, m_d_fix, make_float3(0, -9.81, 0), m_num_particles);

	CuGlobalTorqueStep(pos, m_d_quat, m_d_rest_darboux, m_d_rest_length, m_d_kss, m_d_kbt, m_d_fix, m_d_last_index, m_numElastic);

	CuLocalTorqueStep(m_d_quat, m_d_rest_darboux, m_d_rest_length, m_d_kbt, m_d_fix, m_num_particles);

	// GPUメモリ領域をアンマップ
	CuUnmapGLBufferObject(m_cgr_pos);
}

/*!
 * 空間分割法の準備
 * @param[out] cell 分割グリッドデータ
 * @param[out] gridsize 各軸のグリッド数
 * @param[out] cell_width セル幅
 * @param[in] vMin 環境の最小座標
 * @param[in] vMax 環境の最大座標
 * @param[in] h 有効半径
 */
void SPH::setupNNCell(Cell& c, glm::vec3 vMin, glm::vec3 vMax, double h, int n)
{
	if (h < 1e-6) return;

	glm::vec3 world_size = vMax - vMin;
	glm::vec3 world_origin = vMin;

	c.CellWidth = make_float3(h, h, h);
	c.GridSize.x = (int)(world_size[0] / c.CellWidth.x) + 1;
	c.GridSize.y = (int)(world_size[1] / c.CellWidth.y) + 1;
	c.GridSize.z = (int)(world_size[2] / c.CellWidth.z) + 1;

	c.uNumCells = c.GridSize.x * c.GridSize.y * c.GridSize.z;
	c.uNumParticles = n;

	c.WorldOrigin = make_float3(vMin[0], vMin[1], vMin[2]);
	c.WorldMax = make_float3(vMax[0], vMax[1], vMax[2]);


	if(m_params.cell.dSortedPos) CuFreeArray(m_params.cell.dSortedPos);
	if(m_params.cell.dSortedVel) CuFreeArray(m_params.cell.dSortedVel);
	if(m_params.cell.dGridParticleHash) CuFreeArray(m_params.cell.dGridParticleHash);
	if(m_params.cell.dSortedIndex) CuFreeArray(m_params.cell.dSortedIndex);
	if(m_params.cell.dCellStart) CuFreeArray(m_params.cell.dCellStart);
	if(m_params.cell.dCellEnd) CuFreeArray(m_params.cell.dCellEnd);

	// ソート済み位置・速度，ハッシュ値，ソート済みインデックス，セル始点・終点位置を格納する配列のメモリ確保
	int mem_size = sizeof(float) * n * DIM;
	CuAllocateArray((void**)&m_params.cell.dSortedPos, mem_size);
	CuAllocateArray((void**)&m_params.cell.dSortedVel, mem_size);
	CuAllocateArray((void**)&m_params.cell.dGridParticleHash, n * sizeof(uint));
	CuAllocateArray((void**)&m_params.cell.dSortedIndex, n * sizeof(uint));
	CuAllocateArray((void**)&m_params.cell.dCellStart, m_params.cell.uNumCells * sizeof(uint));
	CuAllocateArray((void**)&m_params.cell.dCellEnd, m_params.cell.uNumCells * sizeof(uint));

	cout << "[grid for nn search] " << endl;
	cout << "  size   : " << c.GridSize << endl;
	cout << "  num    : " << c.uNumCells << endl;
	cout << "  origin : " << c.WorldOrigin << endl;
	cout << "  width  : " << c.CellWidth << endl;
}

/*!
 * 確保したメモリの解放
 */
void SPH::Finalize(void)
{
	if(!m_initialized) return;

	// メモリ解放
	int num_solid = (int)m_solids.size();
	for(int i = 0; i < num_solid; ++i){
		delete m_solids[i];
	}
	m_solids.clear();	

	//if (m_vao_pos != 0) glDeleteVertexArrays(1, &m_vao_pos);//海老沢追加

	if(m_vbo_pos != 0) { CuUnregisterGLBufferObject(m_cgr_pos); deleteVBO(m_vbo_pos); }
	if(m_vbo_col != 0) { CuUnregisterGLBufferObject(m_cgr_col); deleteVBO(m_vbo_col); }
	if (m_vbo_tang != 0) { CuUnregisterGLBufferObject(m_cgr_tang); deleteVBO(m_vbo_tang); }
	if(m_d_vel != 0) CuFreeArray(m_d_vel);
	if(m_d_acc != 0) CuFreeArray(m_d_acc);
	if(m_d_vort != 0) CuFreeArray(m_d_vort);
	if(m_d_dens != 0) CuFreeArray(m_d_dens);
	if(m_d_pres != 0) CuFreeArray(m_d_pres);
	//海老沢追加----------------------------------------------
	if (m_d_mass != 0) CuFreeArray(m_d_mass);
	if (m_d_rest_length != 0)CuFreeArray(m_d_rest_length);
	if (m_d_kss != 0) CuFreeArray(m_d_kss);
	if (m_d_kbt != 0) CuFreeArray(m_d_kbt);
	if (m_d_quat != 0) CuFreeArray(m_d_quat);
	if (m_d_rest_darboux != 0) CuFreeArray(m_d_rest_darboux);
	if (m_d_lambda_ss != 0) CuFreeArray(m_d_lambda_ss);
	if (m_d_lambda_bt != 0) CuFreeArray(m_d_lambda_bt);
	if (m_d_fix != 0)CuFreeArray(m_d_fix);
	if (m_d_curpos != 0)CuFreeArray(m_d_curpos);
	if (m_d_curquat != 0)CuFreeArray(m_d_curquat);
	if (m_d_angvel != 0)CuFreeArray(m_d_angvel);
	if (m_d_rest_density != 0)CuFreeArray(m_d_rest_density);
	if (m_d_fss != 0)CuFreeArray(m_d_fss);
	if (m_d_last_index != 0)CuFreeArray(m_d_last_index);
	if (m_d_pbf_lambda != 0)CuFreeArray(m_d_pbf_lambda);
	//--------------------------------------------------------

	CuFreeArray(m_params.cell.dSortedPos);
	CuFreeArray(m_params.cell.dSortedVel);
	CuFreeArray(m_params.cell.dGridParticleHash);
	CuFreeArray(m_params.cell.dSortedIndex);
	CuFreeArray(m_params.cell.dCellStart);
	CuFreeArray(m_params.cell.dCellEnd);


	m_num_particles = 0;
	m_max_particles = 0;
}


/*!
 * SPHを1ステップ進める
 * @param[in] dt 時間ステップ幅
 * @retval ture  計算完了
 * @retval false 最大ステップ数を超えています
 */
bool SPH::Update(float dt, int step)
{
	// 流体オブジェクトの追加
	//if(!m_pdist.empty()){
	//	vector<rxPdist>::iterator itr = m_pdist.begin();
	//	for(; itr != m_pdist.end(); ++itr){
	//		if(step == itr->steps[0]){
	//			int num = add(*itr, -1, true);
	//			itr->steps.pop_front();
	//		}
	//	}
	//}

	if (m_num_particles == 0) return false;

	assert(m_initialized);
	float h = m_params.effective_radius;
	static bool init = true;

	CuSetParameters(&m_params);

	// VBOにGPUメモリ領域をマッピングしてそのポインタを取得
	float* pos = (float*)CuMapGLBufferObject(&m_cgr_pos);
	float* col = (float*)CuMapGLBufferObject(&m_cgr_col);
	float* tang = (float*)CuMapGLBufferObject(&m_cgr_tang);

	// 近傍探索用グリッドセルに粒子情報を格納(実際にはハッシュを計算しているだけで格納はしていない)
	SetParticlesToCell(pos, m_d_vel, m_num_particles, h);
	
	//海老沢追加
	//シンプルなsphかpbfかをここで決める
	bool sph = false;
	bool pbf = true;

	//海老沢追加
	//速度，加速度が一定(VEL_EPSILON,ANGVEL_EPSILON)以下の場合，速度を切り捨てる
	bool vel_control = false;
	//-------------------------------------------------------------------------
	
	// CUDAカーネルによる粒子処理
	CuSphDensity(m_d_rest_density, m_d_dens, m_d_vol, m_d_mass, m_num_particles);	// 密度計算
	if (sph) CuSphPressure(m_d_rest_density, m_d_pres, m_d_dens, m_num_particles);	// 圧力計算
	CuSphVorticity(m_d_vort, m_d_vel, m_d_dens, m_d_vol, m_d_attr, m_num_particles);	// 渦度計算
	if (sph) CuSphForces(m_d_rest_density, m_d_acc, m_d_vel, m_d_dens, m_d_pres, m_d_vort, m_d_vol, m_d_mass, m_d_attr, m_wind_power, m_d_fss, m_num_particles);	// 粒子にかかる力(圧力項&外力項)の計算
	//海老沢追加
	if (pbf) CuPbfExternalForces(m_d_acc, m_d_attr, m_wind_power, m_num_particles);
	m_method_visc = 2;
	if(m_method_visc == 1){	// WCSPHの論文での粘性項の計算．粘性項を粒子にかかる力として計算する
		CuSphViscosityForces(m_d_rest_density, m_d_acc, m_d_vel, m_d_dens, m_d_vol, m_d_mass, m_d_attr, m_num_particles);	// 粘性項の計算
		CuSphIntegrate(pos, m_d_vel, m_d_acc, m_d_attr, m_d_fix,m_num_particles);	// 速度＆位置の更新(海老沢追加 fix)
	}
	else if(m_method_visc == 2){	// XSPH人工粘性の計算．速度に対して粘性を付加するので速度更新>粘性計算>位置更新の順番となる
		CuSphIntegrateV(m_d_vel, m_d_acc, m_d_attr, m_num_particles);	// 速度のみを更新
		CuSphXSPHViscosity(m_d_vel, m_d_dens, m_d_vol, m_d_mass, m_d_attr, m_num_particles);	// 粘性項の計算(XSPH)
		CuSphIntegrateP(pos, m_d_vel, m_d_attr,m_d_fix, m_num_particles);	// 位置を速度から更新
	}
	else{	// 粘性項なし(時間積分のみ)
		CuSphIntegrate(pos, m_d_vel, m_d_acc, m_d_attr, m_d_fix,m_num_particles);	// 速度＆位置の更新(海老沢追加 fix)
	}

	//RXTIMER("force calculation");
	//海老沢追加---------------------------------------------------------------------------------------------------------------------	
	//角加速度の更新
	CuAngVelUpdate(m_d_angvel, m_d_quat,m_d_fix, dt, m_num_particles);

	//風や重力をイメージした外力の計算(加速度を設定)
	float3 gravity = make_float3(0.f, -9.81, 0.f);
	//gravity = make_float3(0.f, 0.f, 0.f);
	float3 wind = m_wind_power;
	//CuCalExternalForces(pos, m_d_vel,m_d_mass, m_d_fix, gravity,wind, dt, m_num_particles);

	//密度制約の利用
	if (pbf) {
		for (int i = 0; i < 10; i++) {//反復回数元は2回
			//密度制約のためにソートしなおす
			SetParticlesToCell(pos, m_d_vel, m_num_particles, h);
			//密度制約(仮想体積とした場合，不安定になる．)
			CuPbfConstraint(pos, m_d_dens, m_d_rest_density, m_d_pbf_lambda, m_d_vol,m_d_mass, m_num_particles);
		}
	}

	int iter = 64;
	//XPBDの制約の処理(伸び・せん断制約，曲げねじれ制約(中で反復)
	CuXPBDConstraint(pos, m_d_curpos, m_d_mass, m_d_rest_length, m_d_kss, m_d_kbt, m_d_quat, m_d_curquat, m_d_rest_darboux, m_d_lambda_ss, m_d_lambda_bt, m_d_fix, dt, m_num_particles, iter, m_example_flag);

	//摩擦制約の実装--------------------------------------
	SetParticlesToCell(pos, m_d_vel, m_num_particles, h);

	CuFrictionConstraint(pos, m_d_curpos, m_d_rest_density, m_d_vol, m_d_dens, m_d_fix, m_num_particles);
	//摩擦制約の後，姿勢を修正する
	//CuFrictionConstraint_withQuat(pos, m_d_curpos, m_d_rest_density, m_d_vol, m_d_dens, m_d_quat, m_d_rest_length, m_d_fix, m_num_particles);
	//----------------------------------------------------

	//時間積分---------------------------------------------------
	//vel_controlで一定以下の速度を切り捨てるか判定
	CuIntegrate(pos, m_d_curpos, m_d_vel, dt, m_num_particles, vel_control);

	//角加速度の積分
	CuAngVelIntegrate(m_d_angvel, m_d_curquat, m_d_quat, m_d_fix, dt, m_num_particles, vel_control);
	//-----------------------------------------------------------

	//衝突制約
	//球を動かす
	float3 SphereVel = make_float3(3.0, 0.0, 0.0);
	//MoveSphere(SphereVel, dt);

	CuCollisionConstraint(pos, m_d_vel, m_d_fix, m_center, m_rad, dt, m_num_particles);

	//CuPrint3Dfloat(pos, m_d_vel, m_d_acc, m_num_particles);
	
	// 接線の更新
	CuTangUpdate(pos, tang, m_d_fix, m_num_particles);
	//OutputParticles("debug_update.txt");
	//--------------------------------------------------------------------------------------------------------------------------------

	// GPUメモリ領域をアンマップ
	CuUnmapGLBufferObject(m_cgr_pos);
	CuUnmapGLBufferObject(m_cgr_col);
	CuUnmapGLBufferObject(m_cgr_tang);

	// 粒子の描画色の計算
	SetColorVBO(-1);

	//RXTIMER("color(vbo)");
	//RXTIMER_PRINT;

	init = false;
	return true;
}

/*!
 * rest densityの計算
 *  - 近傍に粒子が敷き詰められているとして密度を計算する
 * @param[in] h 有効半径
 * @return rest density
 */
float SPH::calRestDensity(float h)
{
	float r0 = 0.0;
	float l = 2 * m_params.particle_radius;
	int n = (int)ceil(m_params.kernel_radius / l) + 1;
	for (int x = -n; x <= n; ++x) {
		for (int y = -n; y <= n; ++y) {
			for (int z = -n; z <= n; ++z) {
				glm::vec3 rij = glm::vec3(x * l, y * l, z * l);
				r0 += m_params.mass * KernelPoly6(glm::length(rij), h, m_params.aw);
			}
		}
	}
	return r0;
}


//-----------------------------------------------------------------------------
// 近傍探索
//-----------------------------------------------------------------------------
/*!
 * 全粒子を分割セルに格納
 */
void SPH::SetParticlesToCell(float *ppos, float *pvel, int n, float h)
{
	// 近傍探索高速化用グリッドデータの作成
	// 分割セルのハッシュを計算
	CuCalcHash(m_params.cell.dGridParticleHash, m_params.cell.dSortedIndex, ppos, n);

	// ハッシュに基づきパーティクルをソート
	CuSort(m_params.cell.dGridParticleHash, m_params.cell.dSortedIndex, n);

	// パーティクル配列をソートされた順番に並び替え，
	// 各セルの始まりと終わりのインデックスを検索
	CuReorderDataAndFindCellStart(m_params.cell, ppos, pvel, n);

}



/*!
 * 探索用グリッドの描画
 * @param[in] col 粒子が含まれるセルの色
 */
void SPH::DrawCells(glm::vec3 col)
{
	//glDisable(GL_LIGHTING);
	//glLineWidth(1.0);
	//if(m_nn) m_nn->DrawCells(col);
}

/*!
 * 固体障害物の描画
 */
void SPH::DrawObstacles(int drw)
{
	vector<rxSolid*>::iterator i = m_solids.begin(); i++;	// 最初の1つはシーン全体の境界なので描画からは飛ばす
	for(; i != m_solids.end(); ++i){
		(*i)->Draw(drw);
	}
}





/*!
 * VBOからホストメモリへデータを転送，取得
 *  - CPU版ではデータ転送はなく単純に各データを返す
 * @param[in] type データの種類
 * @return ホストメモリ上のデータ
 */
void SPH::GetArrayFromDevice(int type, void* hdata, int num)
{
	assert(m_initialized);
	if(hdata == 0) return;

	if(num == -1) num = m_num_particles;
	void* ddata = 0;

	cudaGraphicsResource **graphics_resource = 0;
	int stride = DIM*sizeof(float);

	// 取得するデータ毎の設定
	switch(type){
	default:
	case P_POSITION: graphics_resource = &m_cgr_pos; break;
	case P_VELOCITY: ddata = m_d_vel; ;  break;
	case P_ACC: ddata = m_d_acc; ; break;//海老沢追加
	case P_FORCE: ddata = m_d_acc; break;
	case P_VORTICITY: ddata = m_d_vort;  break;
	case P_DENSITY: ddata = m_d_dens; stride = sizeof(float);  break;
	case P_PRESSURE: ddata = m_d_pres; stride = sizeof(float);  break;
	case P_ATTRIBUTE: ddata = m_d_attr; stride = sizeof(int); break;
	case P_VOLUME: ddata = m_d_vol; stride = sizeof(float); break;
		//以下、海老沢追加------------------------------------------------------------------------
	case P_MASS: ddata = m_d_mass; stride = sizeof(float); break;
	case E_LENGTH:ddata = m_d_rest_length; stride = sizeof(float); break;
	case E_KSS: ddata = m_d_kss; stride = sizeof(float); break;
	case E_KBT: ddata = m_d_kbt; stride = sizeof(float); break;
	case E_QUAT: ddata = m_d_quat; stride = QUAT*sizeof(float); break;
	case E_DARBOUX: ddata = m_d_rest_darboux; stride = QUAT * sizeof(float); break;
	case E_LAMBDA_SS: ddata = m_d_lambda_ss; stride = DIM * sizeof(float); break;
	case E_LAMBDA_BT: ddata = m_d_lambda_bt; stride = QUAT * sizeof(float); break;
	case P_FIX:ddata = m_d_fix; stride = sizeof(int); break;
	case E_TANG:graphics_resource = &m_cgr_tang; break;
	case E_CURQUAT:ddata = m_d_curquat; stride = QUAT * sizeof(float); break;
	case E_ANGVEL:ddata = m_d_angvel; stride = DIM * sizeof(float); break;
	case P_RESTDENS:ddata = m_d_rest_density; stride = sizeof(float); break;
	case E_FSS:ddata = m_d_fss; stride = DIM*sizeof(float); break;
	case P_LAST_IND:ddata = m_d_last_index; stride = sizeof(int); break;
	case P_PBF_LAMBDA:ddata = m_d_pbf_lambda; stride = sizeof(float); break;
		//----------
	}
	if(ddata || graphics_resource){
		// GPUメモリからCPUメモリへのデータ転送
		// 粒子位置(P_POSITION)はVBOメモリ(OpenGL管理)領域に格納されているのでcudaGraphicsResourceを使って取得
		CuCopyArrayFromDevice(hdata, ddata, graphics_resource, 0, num*stride);
	}
}

/*!
 * ホストメモリからVBO/GPUメモリへデータを転送
 * @param[in] type データの種類
 * @param[in] data ホストメモリ上のデータ
 * @param[in] start データの開始インデックス
 * @param[in] count 追加数
 */
void SPH::SetArrayToDevice(int type, const void* data, int start, int count)
{
	assert(m_initialized);
 
	switch(type){
	default:
	case P_POSITION:
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_pos);
		glBufferSubData(GL_ARRAY_BUFFER, start*DIM*sizeof(float), count*DIM*sizeof(float), data);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		break;

	case P_VELOCITY:
		CuCopyArrayToDevice(m_d_vel, data, start*DIM*sizeof(float), count*DIM*sizeof(float));
		break;

	case P_ATTRIBUTE:
		CuCopyArrayToDevice(m_d_attr, data, start*sizeof(int), count*sizeof(int));
		break;

	case P_VOLUME:
		CuCopyArrayToDevice(m_d_vol, data, start*sizeof(float), count*sizeof(float));
		break;

	//以下、海老沢追加----------------------------------------------------------------------------------------
	case P_MASS:
		CuCopyArrayToDevice(m_d_mass, data, start * sizeof(float), count * sizeof(float));
		break;

	case E_LENGTH:
		CuCopyArrayToDevice(m_d_rest_length, data, start * sizeof(float), count * sizeof(float));//Eから始まるものはエッジの数=粒子数-1(Darbouxは-2)
		break;

	case E_KSS:
		CuCopyArrayToDevice(m_d_kss, data, start * sizeof(float), count * sizeof(float));
		break;

	case E_KBT:
		CuCopyArrayToDevice(m_d_kbt, data, start * sizeof(float), count * sizeof(float));
		break;

	case E_QUAT:
		CuCopyArrayToDevice(m_d_quat, data, start * QUAT * sizeof(float), count * QUAT * sizeof(float));
		CuCopyArrayToDevice(m_d_curquat, data, start * QUAT * sizeof(float), count * QUAT * sizeof(float));//前ステップの姿勢
		break;

	case E_DARBOUX:
		CuCopyArrayToDevice(m_d_rest_darboux, data, start * QUAT * sizeof(float), count * QUAT * sizeof(float));//粒子数-2
		break;

	case E_LAMBDA_SS:
		CuCopyArrayToDevice(m_d_lambda_ss, data, start * DIM * sizeof(float), count * DIM * sizeof(float));
		break;

	case E_LAMBDA_BT:
		CuCopyArrayToDevice(m_d_lambda_bt, data, start * DIM * sizeof(float), count * DIM * sizeof(float));
		break;
	case P_FIX:
		CuCopyArrayToDevice(m_d_fix, data, start * sizeof(int), count * sizeof(int));
		break;
	case P_CURPOS:
		CuCopyArrayToDevice(m_d_curpos, data, start * DIM * sizeof(float), count * DIM * sizeof(float));
		break;
	case E_TANG:
		glBindBuffer(GL_ARRAY_BUFFER, m_vbo_tang);
		glBufferSubData(GL_ARRAY_BUFFER, start * DIM * sizeof(float), count * DIM * sizeof(float), data);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
		break;
	case E_FSS://いらない
		CuCopyArrayToDevice(m_d_fss, data, start * DIM * sizeof(float), count * DIM * sizeof(float));
		break;
	case P_LAST_IND:
		CuCopyArrayToDevice(m_d_last_index, data, start * sizeof(int), count * sizeof(int));
		break;
	case P_PBF_LAMBDA:
		CuCopyArrayToDevice(m_d_pbf_lambda, data, start * sizeof(float), count * sizeof(float));
		break;
	}
}

//海老沢追加
void SPH::SetXPBD_Params(int max_particles, const void* mass, const void* length,const void* kss, const void* kbt, const void* quat, const void* darboux,const void* fix,const void* last_index) {
	//直接全部のパラメータを指定して設定
	SetArrayToDevice(P_MASS, mass, 0, max_particles);
	SetArrayToDevice(E_LENGTH, length, 0, max_particles);
	SetArrayToDevice(E_KSS, kss, 0, max_particles);
	SetArrayToDevice(E_KBT, kbt, 0, max_particles);
	SetArrayToDevice(E_QUAT, quat, 0, max_particles);
	SetArrayToDevice(E_DARBOUX, darboux, 0, max_particles);
	SetArrayToDevice(P_FIX, fix, 0, max_particles);
	SetArrayToDevice(P_LAST_IND, last_index, 0, m_numElastic);//別の配列の長さを指定,個々の毛髪の最後の粒子を保存
	return;
}

/*!
 * カラー値用VBOの編集
 * @param[in] type 色のベースとする物性値
 */
void SPH::SetColorVBO(int type, int picked)
{
	float* col = (float*)CuMapGLBufferObject(&m_cgr_col);

	// 粒子の描画色の計算
	if(type == -1) type = m_colortype;
	switch(type){
	default: 
	case C_CONSTANT: CuColorConstant(col, m_d_attr, make_float3(0.4, 0.4, 1.0), m_num_particles); break;
	case C_DENSITY: CuColorScalar(col, m_d_attr, m_d_dens, m_num_particles, make_float3(0.2, 0.2, 0.2), make_float3(0.6, 0.6, 1.0), make_float2(300.0, 900.0)); break;
	case C_PRESSURE: CuColorScalar(col, m_d_attr, m_d_pres, m_num_particles, make_float3(0.2, 0.2, 0.2), make_float3(1.0, 1.0, 1.0), make_float2(0.0, 1000.0)); break;
	case C_VORTICITY: CuColorVector(col, m_d_attr, m_d_vort, m_num_particles, make_float3(0.0, 0.2, 0.0), make_float3(0.8, 1.0, 0.8), make_float2(0.0, 3.0)); break;
	}

	CuUnmapGLBufferObject(m_cgr_col);
}




//-----------------------------------------------------------------------------
// MARK:陰関数値
//-----------------------------------------------------------------------------
float SPH::GetImplicit(float x, float y, float z)
{
	return 0.0;
}

/*!
* パーティクルからグリッドの陰関数値を計算
* @param[in] pnx,pny,pnz グリッド数の指数 nx=2^pnx
* @param[in] minp グリッドの最小座標
* @param[in] d グリッド幅
* @param[out] hF 陰関数値(nx×ny×nzの配列)
*/
void SPH::CalImplicitField(int n[3], glm::vec3 minp, glm::vec3 d, float *hF)
{
	unsigned int memSize = sizeof(float)*n[0]*n[1]*n[2];

	float *dF = 0;
	CuAllocateArray((void**)&dF, memSize);

	// VBOにGPUメモリ領域をマッピングしてそのポインタを取得
	float* pos = (float*)CuMapGLBufferObject(&m_cgr_pos);

	// 近傍探索用グリッドセルに粒子情報を格納(実際にはハッシュを計算しているだけで格納はしていない)
	SetParticlesToCell(pos, m_d_vel, m_num_particles, m_params.effective_radius);

	int3   gnum = make_int3(n[0], n[1], n[2]);
	float3 gmin = make_float3(minp[0], minp[1], minp[2]);
	float3 glen = make_float3(d[0], d[1], d[2]);

	// 粒子密度を用いたボリュームデータ
	CuSphDensityInGrid(dF, m_d_vol, m_d_attr, m_num_particles, gnum, gmin, glen);

	// GPUメモリ領域をアンマップ
	CuUnmapGLBufferObject(m_cgr_pos);

	CuCopyArrayFromDevice(hF, dF, 0, 0, memSize);

	if(dF) CuFreeArray(dF);
}




//-----------------------------------------------------------------------------
// MARK:シミュデータの出力
//-----------------------------------------------------------------------------
/*!
 * シミュレーション設定(粒子数，範囲，密度，質量など)
 * @param[in] fn 出力ファイル名
 */
void SPH::OutputSetting(string fn)
{
	ofstream fout;
	fout.open(fn.c_str());
	if(!fout){
		cout << fn << " couldn't open." << endl;
		return;
	}

	fout << m_num_particles << endl;
	fout << m_envmin[0] << " " << m_envmin[1] << " " << m_envmin[2] << endl;
	fout << m_envmax[0] << " " << m_envmax[1] << " " << m_envmax[2] << endl;
	fout << m_params.rest_dens << endl;
	fout << m_params.mass << endl;
	fout << m_params.kernel_particles << endl;

	fout.close();
}





/*!
 * シーンのリセット
 */
void SPH::Reset(void)
{
	if(m_num_particles){
		// VBOにGPUメモリ領域をマッピングしてそのポインタを取得
		float* pos = (float*)CuMapGLBufferObject(&m_cgr_pos);
		float* col = (float*)CuMapGLBufferObject(&m_cgr_col);

		// 近傍探索用グリッドセルに粒子情報を格納(実際にはハッシュを計算しているだけで格納はしていない)
		SetParticlesToCell(pos, m_d_vel, m_num_particles, m_params.effective_radius);

		// 粒子体積(流体粒子は単純に粒子質量/静止密度)
		CuSphCalVolume(m_d_vol, m_d_attr, m_num_particles, m_params.mass/m_params.rest_dens);

		// 密度計算(描画色決定用)
		CuSphDensity(m_d_rest_density, m_d_dens, m_d_vol, m_d_mass, m_num_particles);

		// チェック用
		//float avg_val = CuCalAverage(m_d_vol, m_num_particles);
		//cout << "avg. value = " << avg_val << endl;

		// GPUメモリ領域をアンマップ
		CuUnmapGLBufferObject(m_cgr_pos);
		CuUnmapGLBufferObject(m_cgr_col);

		SetColorVBO(-1, -1);
	}
}



/*!
 * 粒子データのセット
 * @param[in] ppos 粒子座標
 * @param[in] pvel 粒子速度
 */
bool SPH::Set(const vector<glm::vec3>& ppos, const vector<glm::vec3>& pvel,const vector<glm::vec3>& ptang,int particles)//海老沢変更中
{
	if (ppos.empty() || (int)ppos.size() != particles) {//m_numparticlesからparticlesに変更中
		return false;
	}

	vector<float> vpos, vvel,etang;
	vector<int> vatt;

	int p = 0;
	for (uint i = 0; i < particles; ++i) {//m_num_particlesからparticlesに変更
		glm::vec3 p0 = ppos[i];
		glm::vec3 v0 = pvel[i];
		glm::vec3 t0 = ptang[i];
		for (uint j = 0; j < 3; ++j) {
			vpos.push_back(p0[j]);
			vvel.push_back(v0[j]);
			etang.push_back(t0[j]);
		}
		if (DIM == 4) {
			vpos.push_back(0.0f);
			vvel.push_back(0.0f);
			etang.push_back(0.0f);
		}
		vatt.push_back(0);
	}

	// GPUメモリに粒子データを転送
	SetArrayToDevice(P_POSITION, &vpos[0], 0, m_num_particles);
	SetArrayToDevice(P_VELOCITY, &vvel[0], 0, m_num_particles);
	SetArrayToDevice(P_ATTRIBUTE, &vatt[0], 0, m_num_particles);

	//海老沢追加 時間積分用にひとつ前のステップを記憶
	SetArrayToDevice(P_CURPOS, &vpos[0], 0, m_num_particles);
	//接線を転送
	SetArrayToDevice(E_TANG, &etang[0], 0, m_num_particles);

	// 粒子描画色設定など
	Reset();

	//海老沢追加----------------------------------------------------------------------------------------------
	//一部パラメータの初期化
	CuSetParametersZero(m_d_angvel, m_d_fss, m_d_pbf_lambda, m_num_particles);//角速度とエッジにかかる力を0に

	float h = m_params.effective_radius;
	float* pos = (float*)CuMapGLBufferObject(&m_cgr_pos);
	//個々に初期の密度を設定
	SetParticlesToCell(pos, m_d_vel, m_num_particles, h);
	CuSphCalVolume(m_d_vol, m_d_attr, m_num_particles, m_params.mass / m_params.rest_dens);
	//初期密度の設定
	CuRestDensSet(pos, m_d_rest_density, m_d_vol, m_d_mass, m_num_particles);
	//密度計算
	CuSphDensity(m_d_rest_density, m_d_dens, m_d_vol, m_d_mass, m_num_particles);
	CuUnmapGLBufferObject(m_cgr_pos);
	//-------------------------------------------------------------------------------------------------------

	return true;
}

/*!
 * 粒子配置の追加(step==0ならそのままシーンに追加)
 * @param[in] pd 粒子配置
 */
int SPH::Add(rxPdist pd)
{
	if(pd.steps.empty() || pd.steps[0] == 0){
		add(pd, -1);
		if(pd.steps.size() >= 2){
			pd.steps.pop_front();
			m_pdist.push_back(pd);
		}
	} else{
		m_pdist.push_back(pd);
	}
	return pd.pdist.size();
}



/*!
 * 粒子配置からシミュレーション空間への粒子の追加
 * @param[in] pd 粒子配置
 * @param[in] start 粒子データ内の追加開始位置(-1で今のデータの後ろに追加)
 * @param[in] h2d ホスト(CPU)からデバイス(GPU)メモリにデータ転送するかどうかのフラグ(初期値true)
 */
int SPH::add(rxPdist pd, int start, bool h2d)
{
	// 粒子配列への追加位置の調整
	uint index;
	if(start < 0){
		index = m_num_particles;
		start = m_num_particles;
	} else{
		index = start;
	}

	int count = 0;
	bool over = false;

	vector<float> vpos, vvel;
	vector<int> vatt;

	// 粒子配置から粒子配列へデータを転送
	vector<glm::vec3>::iterator i = pd.pdist.begin();
	for(; i != pd.pdist.end() ; ++i){
		// 追加位置が配列サイズを超えていたら最初(0)に戻ってデータを上書きする
		if(index >= m_max_particles){
			index = 0;
			over = true;
		}

		// 粒子位置・速度・属性の設定
		for(uint j = 0; j < 3; ++j){
			vpos.push_back((*i)[j]);
			vvel.push_back(pd.vel[j]);
		}
		if(DIM == 4){
			vpos.push_back(0.0f);
			vvel.push_back(0.0f);
		}
		vatt.push_back(pd.attr);

		index++;
		count++;
	}

	cout << "num : " << count << " (" << start << " -> " << m_num_particles+count << ")" << endl;

	// 現在の粒子数情報の更新
	if(over){
		m_num_particles = m_max_particles;
	} else{
		m_num_particles += count;
	}

	// GPUメモリに粒子データを転送
	SetArrayToDevice(P_POSITION, &vpos[0], start, count);
	SetArrayToDevice(P_VELOCITY, &vvel[0], start, count);
	SetArrayToDevice(P_ATTRIBUTE, &vatt[0], start, count);

	// 粒子描画色設定など
	Reset();

	return count;
}


/*!
* ボックス型障害物
* @param[in] cen ボックス中心座標
* @param[in] ext ボックスの大きさ(辺の長さの1/2)
* @param[in] ang ボックスの角度(オイラー角)
* @param[in] vel 初期速度
* @param[in] flg 有効/無効フラグ
*/
void SPH::SetBoxObstacle(glm::vec3 cen, glm::vec3 ext, glm::vec3 ang, glm::vec3 vel, int flg)
{
	// 描画用にm_solidsに登録
	rxSolidBox *box = new rxSolidBox(cen-ext, cen+ext, 1);
	glm::quat q(glm::radians(ang));
	box->SetRotation(q);
	m_solids.push_back(box);

	// CUDA側に登録
	if(m_params.num_box < MAX_BOX_NUM){
		CuBox box;
		box.cen = make_float3(cen[0], cen[1], cen[2]);
		box.ext = make_float3(ext[0], ext[1], ext[2]);
		EulerToMatrix(box.rot, ang[0], ang[1], ang[2]);
		Inverse(box.rot, box.inv_rot);
		box.flg = flg;
		m_params.box[m_params.num_box] = box;
		m_params.num_box++;
	}
}

/*!
* 球型障害物
* @param[in] cen 球体中心座標
* @param[in] rad 球体の半径
* @param[in] vel 初期速度
* @param[in] flg 有効/無効フラグ
*/
void SPH::SetSphereObstacle(glm::vec3 cen, float rad, glm::vec3 vel, int flg)
{
	// 描画用にm_solidsに登録
	rxSolidSphere *sphere = new rxSolidSphere(cen, rad, 1);
	m_solids.push_back(sphere);

	// CUDA側に登録
	if(m_params.num_sphere < MAX_SPHERE_NUM){
		CuSphere sphere;
		sphere.cen = make_float3(cen[0], cen[1], cen[2]);
		sphere.rad = rad;
		sphere.flg = flg;
		m_params.sphere[m_params.num_sphere] = sphere;
		m_params.num_sphere++;
	}
}

/*!
 * 境界粒子を追加/削除
 * @param[out] set 追加/削除のフラグ
 */
int SPH::SetBoundaryParticles(bool set)
{
	if(set){
		// 境界粒子設定
		int nb = 0;
		for(rxSolid *sld : m_solids){
			sld->SetOffset(m_params.particle_radius);
			vector<float> bpos;
			nb = GenerateParticlesOnSurf(rxSolid::GetImplicitG_s, sld, sld->GetMin(), sld->GetMax(), m_params.particle_radius, bpos, DIM);

			if(nb > 0){
				// 粒子配列上での境界粒子のスタート地点
				m_offset_bparticles = m_num_particles;

				// 境界粒子情報の追加
				if(m_num_particles+nb > m_max_particles) nb = m_max_particles-m_num_particles;
				m_num_particles += nb;
				m_num_bparticles += nb;

				// 境界粒子速度と属性(速度は0,属性は1)
				vector<float> bvel(m_num_bparticles*DIM, 0.0f);
				vector<int> batt(m_num_bparticles, 1);

				// GPUメモリに粒子データを転送
				SetArrayToDevice(P_POSITION, &bpos[0], m_offset_bparticles, m_num_bparticles);
				SetArrayToDevice(P_VELOCITY, &bvel[0], m_offset_bparticles, m_num_bparticles);
				SetArrayToDevice(P_ATTRIBUTE, &batt[0], m_offset_bparticles, m_num_bparticles);

				// 境界粒子に合わせて探索領域を拡張(流体粒子の動く範囲外まで拡張する必要はなさそう)
				//glm::vec3 minp = m_solids.back()->GetMin()-glm::vec3(4.0*m_params.particle_radius);
				//glm::vec3 maxp = m_solids.back()->GetMax()+glm::vec3(4.0*m_params.particle_radius);
				//// m_envmax,m_envminに合わせて拡張
				//for(int d = 0; d < 3; ++d){
				//	if(maxp[d] < m_envmax[d]) maxp[d] = m_envmax[d];
				//	if(minp[d] > m_envmin[d]) minp[d] = m_envmin[d];
				//}
				//setupNNCell(m_params.cell, minp, maxp, m_params.effective_radius, m_max_particles);

				Reset();

				////Dump<float>("_volb.txt", m_vol, m_num_bparticles, 1);
			}
		}
		cout << "boundary particles : " << nb << endl;
		return nb;
	}
	else{
		//for(uint i = m_offset_bparticles; i < m_num_particles; ++i){
		//	m_pos[DIM*i+0] = m_pos[DIM*i+1] = m_pos[DIM*i+2] = 0.0;
		//	m_vel[DIM*i+0] = m_vel[DIM*i+1] = m_vel[DIM*i+2] = 0.0;
		//	m_attr[i] = 0;
		//}
		//float vol = m_params.mass/m_params.rest_dens;
		//for(uint i = 0; i < m_max_particles; ++i) m_vol[i] = vol;

		//m_num_particles = m_offset_bparticles;
		//m_num_bparticles = 0;

		//// 分割セルに粒子を登録
		//m_nn->SetObjectToCell(m_pos, m_num_particles, DIM*sizeof(float));

		//return 0;
	}

	return 0;
}



/*!
 * 粒子情報の出力
 * @param[in] fn 出力ファイル名
 */
int SPH::OutputParticles(string fn)
{
	ofstream fout;
	fout.open(fn.c_str(), ios::out|ios::binary);
	if(!fout){
		cout << fn << " couldn't open." << endl;
		return 0;
	}
	else {
		cout << fn << " has opened." << endl;
	}
	vector<float> pos(m_num_particles*DIM), vel(m_num_particles*DIM),tang(m_num_particles*DIM);
	//GetArrayFromDevice(P_POSITION, &pos[0], m_num_particles);
	GetArrayFromDevice(P_VELOCITY, &vel[0], m_num_particles);
	//GetArrayFromDevice(E_TANG, &tang[0], m_num_particles);

	//海老沢追加
	vector<float> mass(m_num_particles), length(m_num_particles), kss(m_num_particles), kbt(m_num_particles), quat(m_num_particles * QUAT), darboux(m_num_particles * QUAT), curquat(m_num_particles * QUAT), angvel(m_num_particles * DIM),restdens(m_num_particles),fss(m_num_particles*DIM),pbf_lambda(m_num_particles);
	vector<int> fix(m_num_particles),last_index(m_numElastic);
	GetArrayFromDevice(P_MASS, &mass[0], m_num_particles);
	GetArrayFromDevice(E_LENGTH, &length[0], m_num_particles);
	GetArrayFromDevice(E_KSS, &kss[0], m_num_particles);
	GetArrayFromDevice(E_KBT, &kbt[0], m_num_particles);
	GetArrayFromDevice(E_QUAT, &quat[0], m_num_particles);
	GetArrayFromDevice(E_DARBOUX, &darboux[0], m_num_particles);
	GetArrayFromDevice(P_FIX, &fix[0], m_num_particles);
	GetArrayFromDevice(E_CURQUAT, &curquat[0], m_num_particles);
	GetArrayFromDevice(E_ANGVEL, &angvel[0], m_num_particles);
	GetArrayFromDevice(P_RESTDENS, &restdens[0], m_num_particles);
	GetArrayFromDevice(E_FSS, &fss[0], m_num_particles);
	GetArrayFromDevice(P_LAST_IND, &last_index[0], m_numElastic);
	GetArrayFromDevice(P_PBF_LAMBDA, &pbf_lambda[0], m_num_particles);

	fout.write((char*)&m_num_particles, sizeof(uint));
	//fout << m_num_particles << endl;
	//for(uint i = 0; i < m_num_particles; ++i){
	//	fout << "pos" << i << ": " << endl;
	//	for(int j = 0; j < 3; ++j){
	//		//fout.write((char*)&pos[DIM*i+j], sizeof(float));
	//		fout << j << ": ";
	//		fout << pos[DIM * i + j] << endl;
	//	}
	//}

	/*for(uint i = 0; i < m_num_particles; ++i) {
		fout << "tang" << i << ": " << endl;
		for (int j = 0; j < 3; ++j) {
			fout << j << ": ";
			fout << tang[DIM * i + j] << endl;
		}
	}*/

	/*for(uint i = 0; i < m_num_particles; ++i){
		for(int j = 0; j < 3; ++j){
			fout.write((char*)&vel[DIM*i+j], sizeof(float));
		}
	}*/
	//海老沢追加
	for (uint i = 0; i < m_num_particles; ++i) {
		fout << "mass" << i << ": ";
		fout << mass[i] << endl;
	}
	for (uint i = 0; i < m_num_particles; ++i) {
		fout << "rest_length" << i << ": ";
		fout << length[i] << endl;
	}
	for (uint i = 0; i < m_num_particles; ++i) {
		fout << "kss" << i << ": ";
		fout << kss[i] << endl;
	}
	for (uint i = 0; i < m_num_particles; ++i) {
		fout << "kbt" << i << ": ";
		fout << kbt[i] << endl;
	}
	for (uint i = 0; i < m_num_particles; ++i) {
		fout << "quat" << i << ": " << endl;
		for (int j = 0; j < 4; ++j) {
			fout << j << ": " << quat[QUAT * i + j] << " " << endl;
		}
	}
	for (uint i = 0; i < m_num_particles; ++i) {
		fout << "darboux" << i << ": " << endl;
		for (int j = 0; j < 4; ++j) {
			fout << j << ": " << darboux[QUAT * i + j] << " " << endl;
		}
	}
	for (uint i = 0; i < m_num_particles; ++i) {
		fout << "fix" << i << ": ";
		fout << fix[i] << endl;
	}
	for (uint i = 0; i < m_num_particles; ++i) {
		fout << "curquat" << i << ": " << endl;
		for (int j = 0; j < 4; ++j) {
			fout << j << ": " << curquat[QUAT * i + j] << " " << endl;
		}
	}
	for (uint i = 0; i < m_num_particles; ++i) {
		fout << "angvel" << i << ": " << endl;
		for (int j = 0; j < 3; ++j) {
			//fout.write((char*)&pos[DIM*i+j], sizeof(float));
			fout << j << ": ";
			fout << angvel[DIM * i + j] << endl;
		}
	}
	for (uint i = 0; i < m_num_particles; ++i) {
		fout << "restdens" << i << ": " << endl;
		fout << restdens[i] << endl;
	}
	for (uint i = 0; i < m_numElastic; i++) {
		fout << "last_index" << i << ": " << endl;
		fout << last_index[i] << endl;
	}
	for (uint i = 0; i < m_num_particles; ++i) {
		fout << "fss" << i << ": " << endl;
		for (int j = 0; j < 3; ++j) {
			fout << j << ": ";
			fout << fss[DIM * i + j] << endl;
		}
	}
	for (uint i = 0; i < m_num_particles; ++i) {
		fout << "pbf_lambda" << i << endl;
		fout << pbf_lambda[i] << endl;
	}
	fout.close();

	pos.clear(); vel.clear(); mass.clear(); quat.clear();

	return 1;
}


/*!
 * ファイルから粒子情報を読み込む
 * @param[in] stp ステップ数(sph_ステップ数.datのファイル名から読み込む)
 * @param[out] ppos 粒子座標
 * @param[out] pvel 粒子速度
 */
int SPH::InputParticles(string fn)
{
	ifstream fin;
	fin.open(fn.c_str(), ios::in|ios::binary);
	if(!fin){
		cout << fn << " couldn't find." << endl;
		return 0;
	}

	uint n;
	fin.read((char*)&n, sizeof(uint));
	m_num_particles = n;
	vector<float> pos(m_num_particles*DIM), vel(m_num_particles*DIM);

	for(uint i = 0; i < n; ++i){
		for(int j = 0; j < 3; ++j){
			fin.read((char*)&pos[DIM*i+j], sizeof(float));
		}
	}

	if(!fin.eof()){
		for(uint i = 0; i < n; ++i){
			for(int j = 0; j < 3; ++j){
				fin.read((char*)&vel[DIM*i+j], sizeof(float));
			}
		}
	}

	fin.close();

	SetArrayToDevice(P_POSITION, &pos[0], 0, m_num_particles);
	SetArrayToDevice(P_VELOCITY, &vel[0], 0, m_num_particles);
	SetColorVBO(m_colortype);

	return 1;
}

//海老沢追加
void SPH::ChangeNumParticles(int n) {
	m_num_particles = n;
	return;
}

//海老沢追加
//球を動かす
void SPH::MoveSphere(float3 vel,float dt) {
	m_center.x += dt * vel.x;
	m_center.y += dt * vel.y;
	m_center.z += dt * vel.z;
}

