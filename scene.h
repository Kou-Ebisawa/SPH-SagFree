/*!
  @file scene.h

  @brief GLFWによるOpenGL描画

  @author Makoto Fujisawa
  @date   2021-11
*/

#ifndef _RX_CONTROLLER_H_
#define _RX_CONTROLLER_H_



//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "utils.h"
#include "imgui.h"

// 設定ファイル
#include "rx_atom_ini.h"

// トラックボール＆テクスチャ
#include "rx_trackball.h"
#include "rx_texture.h"

// メッシュファイル読み込み関連
#include "rx_obj.h"

// ファイルによるシミュレーションシーン設定
#include "config.h"

//SagFreeとの統合(変更)---------------------------------------------------
#include "pbd.h"
//-------------------------------------------------------------------------

using namespace std;

//-----------------------------------------------------------------------------
// 定義/定数
//-----------------------------------------------------------------------------
#define RX_OUTPUT_TIME

const string RX_PROGRAM_NAME = "glfw_pbf";

//SagFree追加------------------------------------------------------------------
const int MAX_STEPS = 1000000;
const float GOUND_HEIGHT = -1.0;

const float DT = 0.01;

//PBDの描画フラグ
const int RXD_VERTEX = 0x0001;
const int RXD_EDGE = 0x0002;
const int RXD_FACE = 0x0004;
const int RXD_TETRA = 0x0008;	//!< 四面体(ワイヤフレーム)
const int RXD_TETCS = 0x0010;	//!< 四面体(断面)
const int RXD_BBOX = 0x0020;	//!< AABB(シミュレーション空間)
const int RXD_AXIS = 0x0040;   //!< 軸
const int RXD_FLOOR = 0x0080;	//!< 床

//-----------------------------------------------------------------------------
//以下SagFree追加部分--------------------------------------------------------------------------------
class ScenePBD
{
public:
	// 形状フラグ(海老沢追加)
	enum ElasticFigure
	{
		STRAIGHT = 0,
		CENTER_SPIRAL,
		NATURAL_SPIRAL,
		EXAMPLE,
		ON_SPHERE
	};
protected:
	static int m_winw;						//!< 描画ウィンドウの幅
	static int m_winh;						//!< 描画ウィンドウの高さ
	static rxTrackball m_view;				//!< トラックボール

	static float m_bgcolor[3];				//!< 背景色
	static bool m_animation_on;				//!< アニメーションON/OFF

	static int m_draw;						//!< 描画フラグ
	static int m_currentstep;				//!< 現在のステップ数
	static int m_simg_spacing;				//!< 画像保存間隔(=-1なら保存しない)

	// PBDによる弾性変形
	static float m_dt;
	static ElasticPBD* m_elasticbody;
	//弾性体の配列を作成
	static vector<ElasticPBD*> m_elasticbodies;
	
	static int m_num_elasticbodies;//弾性体の数
	static int m_allnum;//全ての粒子数の合計

	//海老沢追加 複数の弾性体の速度、位置を配列にまとめる(SPHに渡すとき用)
	static vector<glm::vec3> m_NewVelocity;
	static vector<glm::vec3> m_NewPos;

	//------------------------------------------------------------------------
	//! マウスピック
	static int m_picked;
	static float m_pickdist;	//!< ピックされた点までの距離

	//海老沢追加
	static SceneConfig sph_scene;
	static int sph_scene_idx;

public:
	//! コンストラクタ
	ScenePBD() {}

	//! デストラクタ
	~ScenePBD() {}

public:
	// コールバック関数
	static void Init(int argc, char* argv[]);
	static void Draw();
	static void Timer();
	static void Cursor(GLFWwindow* window, double xpos, double ypos);
	static void Mouse(GLFWwindow* window, int button, int action, int mods);
	static void Keyboard(GLFWwindow* window, int key, int mods);
	static void Resize(GLFWwindow* window, int w, int h);
	static void ImGui(GLFWwindow* window);
	static void Destroy();

	//海老沢追加------------------------------------------------
	static void initSPH(int max_particles,int num_elastic,float mass, bool use_bp = true);
	//XBPD->SPHでのデータのやり取りをする配列
	static void initArray(int max_particles);
	static void clearArray(void);

	//配列で弾性体を扱う
	static void initElasticbodies(int num_elasticbodies, int num_particles, float kss, float kbt,int type);
	static void XPBDtoSPH(void);
	static void SPHtoXPBD(void);
	static void XPBDParamsToDevice(int num_elasticbodies);
	static void DestroyElasticArray(void);
	
	//髪型のobjファイルを読み込む
	static bool readObjFile(const char* filename, vector<glm::vec3>& PosArray, vector<glm::ivec2>& IndexArray, vector<int>& FixArray);
	//髪型用のobjファイルに書き換え
	static bool MakeHairObjFile(const char* In_filename,const char* Out_filename);

	static void initElasticFromObj(vector<glm::vec3>PosArray, vector<glm::ivec2>IndexArray, vector<int>FixArray, float ks, float kbt,float mass,int& num_elastic,int &all_particles);

	static void drawHairObject(unsigned int vbo, int n, unsigned int color_vbo,unsigned int normal_vbo, float* data, int offset, double pscale, double prad, double czf = 1000, double czb = -1000);
	//風の強さを変える
	static void ChangeWindPower(float3 wind);

	//OBJファイルを(0,1)でクランプする
	static void FitVertices(vector<glm::vec3> &vertices);
	//-----------------------------------------------------------

private:
	// アニメーション切り替え
	static bool switchanimation(int on);

	// ファイル入出力
	static void savedisplay(const int& stp);

	// 視点
	static void resetview(void);

	// PBD
	static void initStraightRod(void);
	static void initCenterSpiralRod(void);
	static void initNaturalSpiralRod(void);
	static void initExampleRod(void);
	static void initMoreRod(bool sag_free_flag);

	// マウスピック
	static void clearPick(void);
};


#endif // #ifdef _RX_CONTROLLER_H_
