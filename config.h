/*!
  @file config.h
	
  @brief シーン設定をファイルから読み込む
 
  @author Makoto Fujisawa
  @date 2012-08,2022-07
*/

#ifndef _RX_CONFIG_H_
#define _RX_CONFIG_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
// STL
#include <vector>
#include <string>

// ユーティリティ
#include "utils.h"

// シミュレーション
#include "sph.h"
#include "pdist.h"

#include "rx_mesh.h"

// 設定ファイル
#include "rx_atom_ini.h"

using namespace std;


//-----------------------------------------------------------------------------
// 文字列処理関数
//-----------------------------------------------------------------------------

/*!
* "(x, y, z)"の形式の文字列からglm::vec3型へ変換
*  - (x)となっていたら(x, x, x)とする．
* @param[in] s 文字列
* @param[out] v 値
* @return 要素記述数
*/
inline int StringToVec3(const string &s, glm::vec3 &v)
{
	int vcount = 0;
	size_t pos;
	v = glm::vec3(0.0);
	if((pos = s.find('(')) != string::npos){
		while(pos != string::npos && vcount < 3){
			size_t pos1 = pos;
			if((pos1 = s.find(',', pos+1)) != string::npos){
				v[vcount] = atof(s.substr(pos+1, (pos1-(pos+1))).c_str());
				vcount++;
				pos = pos1;
			}
			else if((pos1 = s.find(')', pos+1)) != string::npos){
				v[vcount] = atof(s.substr(pos+1, (pos1-(pos+1))).c_str());
				vcount++;
				break;
			}
			else{
				break;
			}
		}
	}
	if(vcount < 3){
		for(int i = vcount; i < 3; ++i){
			v[i] = v[vcount-1];
		}
	}

	return vcount;
}

/*!
* "(x, y, z, w)"の形式の文字列からVec4型へ変換
*  - (x)となっていたら(x, x, x, x)とする．
* @param[in] s 文字列
* @param[out] v 値
* @return 要素記述数
*/
inline int StringToVec4(const string &s, glm::vec4 &v)
{
	int vcount = 0;
	size_t pos;
	v = glm::vec4(0.0);
	if((pos = s.find('(')) != string::npos){
		while(pos != string::npos && vcount < 4){
			size_t pos1 = pos;
			if((pos1 = s.find(',', pos+1)) != string::npos){
				v[vcount] = atof(s.substr(pos+1, (pos1-(pos+1))).c_str());
				vcount++;
				pos = pos1;
			} else if((pos1 = s.find(')', pos+1)) != string::npos){
				v[vcount] = atof(s.substr(pos+1, (pos1-(pos+1))).c_str());
				vcount++;
				break;
			} else{
				break;
			}
		}
	}
	if(vcount < 4){
		for(int i = vcount; i < 4; ++i){
			v[i] = v[vcount-1];
		}
	}

	return vcount;
}

/*!
 * 文字列から値(glm::vec3,Vec4)を取得
 * @param[out] val 値
 * @param[in] str 文字列
 * @param[in] rel trueでシミュレーション空間の大きさに対する係数として計算
 * @param[in] cen シミュレーション空間の中心座標
 * @param[in] ext シミュレーション空間の大きさ(各辺の長さの1/2)
 */
inline void GetValueFromString(glm::vec3 &val, const string &str, bool rel = false, glm::vec3 cen = glm::vec3(0.0), glm::vec3 ext = glm::vec3(1.0))
{
	int n = StringToVec3(str, val);
	if(rel){
		val = cen+(n == 1 ? glm::vec3(glm::min(ext[0], ext[1], ext[2])) : ext)*val;
	}
}
inline void GetValueFromString(glm::vec4 &val, const string &str)
{
	int n = StringToVec4(str, val);
}

/*!
* 文字列から値(glm::float3,float4)を取得
* @param[out] val 値
* @param[in] str 文字列
* @param[in] rel trueでシミュレーション空間の大きさに対する係数として計算
* @param[in] cen シミュレーション空間の中心座標
* @param[in] ext シミュレーション空間の大きさ(各辺の長さの1/2)
*/
inline void GetValueFromString(float3 &val, const string &str, bool rel = false, glm::vec3 cen = glm::vec3(0.0), glm::vec3 ext = glm::vec3(1.0))
{
	glm::vec3 gval;
	int n = StringToVec3(str, gval);
	if(rel){
		gval = cen+(n == 1 ? glm::vec3(glm::min(ext[0], ext[1], ext[2])) : ext)*gval;
	}
	val.x = gval[0]; val.y = gval[1]; val.z = gval[2];
}
inline void GetValueFromString(float4 &val, const string &str)
{
	glm::vec4 gval;
	int n = StringToVec4(str, gval);
	val.x = gval[0]; val.y = gval[1]; val.z = gval[2]; val.w = gval[3];
}

/*!
 * 文字列から値(double)を取得
 * @param[out] val 値
 * @param[in] str 文字列
 * @param[in] rel trueでシミュレーション空間の大きさに対する係数として計算
 * @param[in] cen シミュレーション空間の中心座標
 * @param[in] ext シミュレーション空間の大きさ(各辺の長さの1/2)
 */
inline void GetValueFromString(double &val, const string &str, bool rel = false, glm::vec3 cen = glm::vec3(0.0), glm::vec3 ext = glm::vec3(1.0))
{
	val = atof(str.c_str());
	if(rel){
		val = glm::min(ext[0], ext[1], ext[2])*val;
	}
}



//-----------------------------------------------------------------------------
//! SceneConfigクラス - SPHのシーン設定をファイルから読み込む
//-----------------------------------------------------------------------------
class SceneConfig
{
	// クラスのメンバをコールバック関数にする場合の関数宣言と関数定義
	#define RXSETFUNC_DECL(fnc) static void fnc(string *ns, string *vl, int n, string hd, void* x); inline void fnc(string *ns, string *vl, int n, string hd);
	#define RXSETFUNC(cls, fnc) static void cls::fnc(string *ns, string *vl, int n, string hd, void* x){ ((cls*)x)->fnc(ns, vl, n, hd); } 

protected:
	SPH *m_pSolver;						//!< ソルバ
	SceneParameter m_Env;				//!< シーン環境設定
	
	string m_strCurrentScene;			//!< 現在のシーンの名前
	vector<string> m_vSceneFiles;		//!< シーンファイルリスト
	int m_iSceneFileNum;				//!< シーンファイルの数

	vector<string> m_vSceneTitles;		//!< シーンファイルリスト
	int m_iCurrentSceneIdx;				//!< 現在のシーンファイル

	vector<rxPolygons> m_vSolidPoly;	//!< 固体メッシュ

	const int MAX_SCENE_NUM = 16;

public:
	//! デフォルトコンストラクタ
	SceneConfig() : m_pSolver(0)
	{
		m_vSceneFiles.resize(12, "");	// シーンファイルリスト
		m_vSceneTitles.resize(12, "");	// シーンタイトルリスト
		m_iCurrentSceneIdx = -1;		// 現在のシーンファイル
		m_strCurrentScene = "null";		// 現在のシーンの名前
		m_iSceneFileNum = 0;
		InitSceneParameter(m_Env);

		Clear();
	}

	//! デストラクタ
	~SceneConfig(){}
	
	//! 設定初期化
	void Clear(void)
	{
		m_pSolver = 0;
		m_vSolidPoly.clear();
	}


	//! シーンタイトルリスト
	vector<string> GetSceneTitles(void) const { return m_vSceneTitles; }

	//! 現在のシーン
	int GetCurrentSceneIdx(void) const { return m_iCurrentSceneIdx; }

	//! シミュレーション環境
	SceneParameter GetEnv(void) const { return m_Env; }

	//! シミュレーションクラス
	void Set(SPH *solver){ m_pSolver = solver; }

	//! 固体ポリゴン
	int GetSolidPolyNum(void) const { return (int)m_vSolidPoly.size(); }
	vector<rxPolygons>& GetSolidPolys(void){ return m_vSolidPoly; }

public:
	/*!
	 * 設定からパラメータの読み込み
	 *  - 追加オブジェクト(Box型の液体とか)以外の全体的なパラメータ設定
	 * @param[in] names 項目名リスト
	 * @param[in] values 値リスト
	 * @param[in] n リストのサイズ
	 * @param[in] header ヘッダ名
	 */
	RXSETFUNC(SceneConfig, SetSpace)
	void SetSpace(string *names, string *values, int n, string header)
	{
		m_Env.use_inlet = 0;
		int idr = 0;
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")      GetValueFromString(m_Env.boundary_cen, values[i], false);
			else if(names[i] == "ext") GetValueFromString(m_Env.boundary_ext, values[i], false);
			else if(names[i] == "max_particle_num")  m_Env.max_particles = atoi(values[i].c_str());
			else if(names[i] == "density")			 m_Env.dens = atof(values[i].c_str());
			else if(names[i] == "mass")				 m_Env.mass = atof(values[i].c_str());
			else if(names[i] == "kernel_particles")  m_Env.kernel_particles = atof(values[i].c_str());
			else if(names[i] == "inlet_boundary")	 m_Env.use_inlet = atoi(values[i].c_str());
			else if(names[i] == "dt")				 m_Env.dt = atof(values[i].c_str());
			else if(names[i] == "viscosity")		 m_Env.viscosity = atof(values[i].c_str());
			else if(names[i] == "vorticity")		 m_Env.vorticity = atof(values[i].c_str());
			else if(names[i] == "gas_stiffness")	 m_Env.gas_k = atof(values[i].c_str());
			// 表面メッシュ
			else if(names[i] == "mesh_res_max")		m_Env.mesh_n = atof(values[i].c_str());
			else if(names[i] == "mesh_thr")			m_Env.mesh_thr = atoi(values[i].c_str());
			// 視点
			else if(names[i] == "view_trans"){		GetValueFromString(m_Env.view_trans, values[i], false); m_Env.view = 1; }
			else if(names[i] == "view_rot"){		GetValueFromString(m_Env.view_rot, values[i], false); m_Env.view = 1; }
			else if(names[i] == "view_quat"){		GetValueFromString(m_Env.view_quat, values[i]); m_Env.view = 2; }
			else if(names[i] == "background"){		GetValueFromString(m_Env.bgcolor, values[i], false); m_Env.view = 1; }
		}
	}

	/*!
	 * 粒子や固体オブジェクトの設定読み込み
	 * - 新しいオブジェクトを設定ファイルに追加する場合は，読み込み関数(SetLiquid*やSetSolid*)を追加した後，
	 *   この関数内に設定ファイル内の用語との対応関係を記述
	 */
	bool LoadSceneFromFile(void)
	{
		bool ok = true;
		rxINI* cfg = new rxINI();
		cfg->SetHeaderFunc("liquid box", &SceneConfig::SetLiquidBox, this);
		cfg->SetHeaderFunc("liquid box (r)", &SceneConfig::SetLiquidBox, this);
		cfg->SetHeaderFunc("liquid hollow box", &SceneConfig::SetLiquidHollowBox, this);
		cfg->SetHeaderFunc("liquid hollow box (r)", &SceneConfig::SetLiquidHollowBox, this);
		cfg->SetHeaderFunc("liquid sphere", &SceneConfig::SetLiquidSphere, this);
		cfg->SetHeaderFunc("liquid sphere (r)", &SceneConfig::SetLiquidSphere, this);
		cfg->SetHeaderFunc("solid box", &SceneConfig::SetSolidBox, this);
		cfg->SetHeaderFunc("solid box (r)", &SceneConfig::SetSolidBox, this);
		cfg->SetHeaderFunc("solid sphere", &SceneConfig::SetSolidSphere, this);
		cfg->SetHeaderFunc("solid sphere (r)", &SceneConfig::SetSolidSphere, this);
		cfg->SetHeaderFunc("solid polygon", &SceneConfig::SetSolidPolygon, this);
		cfg->SetHeaderFunc("solid polygon (r)", &SceneConfig::SetSolidPolygon, this);
		cfg->SetHeaderFunc("inlet line", &SceneConfig::SetInletLine, this);
		cfg->SetHeaderFunc("inlet line (r)", &SceneConfig::SetInletLine, this);
		if(!(cfg->Load(m_strCurrentScene))){
			cout << "Failed to open the " << m_strCurrentScene << " file!" << endl;
			ok = false;
		}
		delete cfg;
		return ok;
	}


	//! 液体 : 箱形
	RXSETFUNC(SceneConfig, SetLiquidBox)
	void SetLiquidBox(string* names, string* values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		glm::vec3 cen(0.0), ext(0.0), vel(0.0);
		int step = 0;
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")       GetValueFromString(cen, values[i], rel, m_pSolver->GetCen(), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "ext")  GetValueFromString(ext, values[i], rel, glm::vec3(0.0), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "vel")  GetValueFromString(vel, values[i], false);
			else if(names[i] == "step")	step = atoi(values[i].c_str());
		}
		cout << "set liquid box : " << glm::to_string(cen) << ", " << glm::to_string(ext) << ", " << glm::to_string(vel) << endl;
		m_pSolver->Add(MakeBox(cen, ext, vel, m_pSolver->GetSpacing(), 0, step));
	}

	//! 液体 : 箱形(中空)
	RXSETFUNC(SceneConfig, SetLiquidHollowBox)
	void SetLiquidHollowBox(string* names, string* values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		glm::vec3 cen(0.0), ext0(0.0), ext1(0.0), vel(0.0);
		int step = 0;
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")       GetValueFromString(cen, values[i], rel, m_pSolver->GetCen(), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "ext0") GetValueFromString(ext0, values[i], rel, glm::vec3(0.0), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "ext1") GetValueFromString(ext1, values[i], rel, glm::vec3(0.0), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "vel")  GetValueFromString(vel, values[i], false);
			else if(names[i] == "step")	step = atoi(values[i].c_str());
		}
		cout << "set liquid box(hollow) : " << glm::to_string(cen) << ", " << glm::to_string(ext0) << "-" << glm::to_string(ext1) << ", " << glm::to_string(vel) << endl;
		m_pSolver->Add(MakeHollowBox(cen, ext0, ext1, vel, m_pSolver->GetSpacing(), 0, step));
	}

	//! 液体 : 球
	RXSETFUNC(SceneConfig, SetLiquidSphere)
	void SetLiquidSphere(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		glm::vec3 cen(0.0), vel(0.0);
		double rad = 0.0;
		int step = 0;
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")       GetValueFromString(cen, values[i], rel, m_pSolver->GetCen(), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "rad")  GetValueFromString(rad, values[i], rel, glm::vec3(0.0), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "vel")  GetValueFromString(vel, values[i], false);
			else if(names[i] == "step")	step = atoi(values[i].c_str());
		}
		cout << "set liquid sphere : " << glm::to_string(cen) << ", " << rad << endl;
		m_pSolver->Add(MakeSphere(cen, rad, vel, m_pSolver->GetSpacing(), 0, step));
	}

	//! 液体流入 : 線分
	RXSETFUNC(SceneConfig, SetInletLine)
	void SetInletLine(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		glm::vec3 pos1(0.0), pos2(0.0), vel(0.0), up(0.0, 1.0, 0.0);
		int  span = -1, accum = 1, steps = -1;
		double spacing = 1.0;
		for(int i = 0; i < n; ++i){
			if(names[i] == "pos1")         GetValueFromString(pos1, values[i], rel, m_pSolver->GetCen(), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "pos2")    GetValueFromString(pos2, values[i], rel, m_pSolver->GetCen(), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "vel")     GetValueFromString(vel,  values[i], false);
			else if(names[i] == "up")      GetValueFromString(up,   values[i], false);
			else if(names[i] == "span")    span = atoi(values[i].c_str());
			else if(names[i] == "accum")   accum = atoi(values[i].c_str());
			else if(names[i] == "spacing") spacing = atof(values[i].c_str());
			else if(names[i] == "steps")   steps = atof(values[i].c_str());
		}

		cout << "set inlet boundary : " << glm::to_string(pos1) << "-" << glm::to_string(pos2) << ", " << glm::to_string(vel) << endl;
		cout << "                     span=" << span << ", up=" << glm::to_string(up) << ", accum=" << accum << ", spacing=" << spacing << endl;
		int num_of_inlets = m_pSolver->Add(MakeInletLine(pos1, pos2, vel, up, accum, span, spacing+m_pSolver->GetSpacing(), 0, steps));
	}

	//! 固体 : 箱形
	RXSETFUNC(SceneConfig, SetSolidBox)
	void SetSolidBox(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		glm::vec3 cen(0.0), ext(0.0), ang(0.0), vel(0.0);
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")      GetValueFromString(cen, values[i], rel, m_pSolver->GetCen(), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "ext") GetValueFromString(ext, values[i], rel, glm::vec3(0.0), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "ang") GetValueFromString(ang, values[i], false);
			else if(names[i] == "vel") GetValueFromString(vel, values[i], false);
		}
		m_pSolver->SetBoxObstacle(cen, ext, ang, vel, 1);
		cout << "set solid box : " << glm::to_string(cen) << ", " << glm::to_string(ext) << ", " << glm::to_string(ang) << endl;
	}

	//! 固体 : 球
	RXSETFUNC(SceneConfig, SetSolidSphere)
	void SetSolidSphere(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		glm::vec3 cen(0.0), move_pos1(0.0), move_pos2(0.0), vel(0.0);
		int  move = 0, move_start = -1;
		double rad = 0.0, move_max_vel = 0.0, lap = 1.0;
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")      GetValueFromString(cen, values[i], rel, m_pSolver->GetCen(), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "rad") GetValueFromString(rad, values[i], rel, glm::vec3(0.0), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "vel")  GetValueFromString(vel, values[i], false);
			else if(names[i] == "move_pos1") GetValueFromString(move_pos1, values[i], rel, m_pSolver->GetCen(), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "move_pos2") GetValueFromString(move_pos2, values[i], rel, m_pSolver->GetCen(), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "move") move = atoi(values[i].c_str());
			else if(names[i] == "move_start") move_start = atoi(values[i].c_str());
			else if(names[i] == "move_max_vel") move_max_vel = atof(values[i].c_str());
			else if(names[i] == "lap") lap = atof(values[i].c_str());
		}
		m_pSolver->SetSphereObstacle(cen, rad, vel, 1);
		cout << "set solid sphere : " << glm::to_string(cen) << ", " << rad << endl;
	}

	//! 固体 : ポリゴン
	RXSETFUNC(SceneConfig, SetSolidPolygon)
	void SetSolidPolygon(string *names, string *values, int n, string header)
	{
		bool rel = (header.find("(r)") != string::npos);
		string fn_obj;
		glm::vec3 cen(0.0), ext(0.0), ang(0.0), vel(0.0);
		for(int i = 0; i < n; ++i){
			if(names[i] == "cen")       GetValueFromString(cen, values[i], rel, m_pSolver->GetCen(), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "ext")  GetValueFromString(ext, values[i], rel, glm::vec3(0.0), 0.5f*m_pSolver->GetDim());
			else if(names[i] == "ang")  GetValueFromString(ang, values[i], false);
			else if(names[i] == "vel")  GetValueFromString(vel, values[i], false);
			else if(names[i] == "file") fn_obj = values[i];
		}
		if(!fn_obj.empty()){
			//m_pSolver->SetPolygonObstacle(fn_obj, cen, ext, ang, vel);
			cout << "set solid polygon : " << fn_obj << endl;
			cout << "                  : " << glm::to_string(cen) << ", " << glm::to_string(ext) << ", " << glm::to_string(ang) << endl;
		}
	}


	/*!
	 * シミュレーション空間の設定読み込み
	 */
	bool LoadSpaceFromFile(void)
	{
		bool ok = true;
		rxINI *cfg = new rxINI();
		cfg->SetHeaderFunc("space", &SceneConfig::SetSpace, this);
		if(!(cfg->Load(m_strCurrentScene))){
			cout << "Failed to open the " << m_strCurrentScene << " file!" << endl;
			ok = false;
		}
		delete cfg;
		return ok;
	}
	
	void LoadSpaceFromFile(const string input)
	{
		// SPH設定をファイルから読み込み
		ifstream fsin;
		fsin.open(input.c_str());

		glm::vec3 bmin, bmax, bcen, bext;
		fsin >> m_Env.max_particles;
		fsin >> bmin[0] >> bmin[1] >> bmin[2];
		fsin >> bmax[0] >> bmax[1] >> bmax[2];
		fsin >> m_Env.dens;
		fsin >> m_Env.mass;
		fsin >> m_Env.kernel_particles;

		bcen = 0.5f*(bmax+bmin);
		bext = 0.5f*(bmax-bmin);
		m_Env.boundary_ext = make_float3(bcen[0], bcen[1], bcen[2]);
		m_Env.boundary_cen = make_float3(bext[0], bext[1], bext[2]);

		fsin.close();

		cout << "[SPH - " << input << "]" << endl;
		cout << " num. of particles : " << m_Env.max_particles << endl;
		cout << " boundary min      : " << glm::to_string(bmin) << endl;
		cout << " boundary max      : " << glm::to_string(bmax) << endl;
		cout << " boundary cen      : " << glm::to_string(m_Env.boundary_cen) << endl;
		cout << " boundary ext      : " << glm::to_string(m_Env.boundary_ext) << endl;
		cout << " density           : " << m_Env.dens << endl;
		cout << " mass              : " << m_Env.mass << endl;
		cout << " kernel particles  : " << m_Env.kernel_particles << endl;
	}

	/*!
	 * 指定したフォルダにある設定ファイルの数とシーンタイトルを読み取る
	 * @param[in] dir 設定ファイルがあるフォルダ(何も指定しなければ実行フォルダ)
	 */
	void ReadSceneFiles(string dir = "")
	{
		m_vSceneFiles.resize(MAX_SCENE_NUM, "");	// シーンファイルリスト
		m_vSceneTitles.resize(MAX_SCENE_NUM, "");	// シーンタイトルリスト

		ifstream scene_ifs;
		string scene_fn = "null";
		for(int i = 1; i <= MAX_SCENE_NUM; ++i){
			if(ExistFile((scene_fn = CreateFileName(dir+"sph_scene_", ".cfg", i, 1)))){
				cout << "scene " << i << " : " << scene_fn << endl;
				m_vSceneFiles[i-1] = scene_fn;
				m_vSceneTitles[i-1] = scene_fn.substr(0, 11);

				// シーンタイトルの読み取り
				scene_ifs.open(scene_fn.c_str(), ios::in);
				string title_buf;
				getline(scene_ifs, title_buf);
				if(!title_buf.empty() && title_buf[0] == '#'){
					m_vSceneTitles[i-1] = title_buf.substr(2, title_buf.size()-2);
				}
				scene_ifs.close();

				m_iSceneFileNum++;
			}
		}

		if(m_iSceneFileNum){
			SetCurrentScene(0);
		}
	}

	/*!
	 * カレントのシーン設定
	 * @param[in] シーンインデックス
	 */
	bool SetCurrentScene(int idx)
	{
		if(idx < 0 || idx >= MAX_SCENE_NUM || m_vSceneFiles[idx].empty()){
			cout << "There is no scene files!" << endl;
			return false;
		}
		m_iCurrentSceneIdx = idx;
		m_strCurrentScene = m_vSceneFiles[m_iCurrentSceneIdx];
		return true;
	}

	/*!
	 * タイトルからカレントのシーン設定
	 * @param[in] label シーンタイトル
	 */
	bool SetCurrentSceneFromTitle(const string label)
	{
		int scene = 0;
		while(label.find(m_vSceneTitles[scene]) == string::npos || fabs((int)(m_vSceneTitles[scene].size()-label.size())) > 3) scene++;

		if(m_iCurrentSceneIdx != -1 && m_vSceneFiles[scene] != ""){
			cout << "scene " << scene+1 << " : " << m_vSceneFiles[scene] << endl;
			m_iCurrentSceneIdx = scene;
			m_strCurrentScene = m_vSceneFiles[scene];
			return true;
		}

		return false;
	}


	vector<string> GetSceneTitleList(void){ return m_vSceneTitles; }
};


#endif // #ifndef _RX_SPH_CONFIG_H_