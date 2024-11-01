/*!
  @file pbd.h
	
  @brief PBDによる弾性体シミュレーション
 
  @author Makoto Fujisawa
  @date 2021
*/

#ifndef _RX_PBD_STRAND_H_
#define _RX_PBD_STRAND_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "utils.h"

#include "rx_mesh.h"

#pragma once 
#include <Eigen/Dense>

using namespace std;

//! 衝突判定用関数
typedef void (*CollisionFunc)(glm::vec3&, glm::vec3&, glm::vec3&, int);

//-----------------------------------------------------------------------------
// ElasticPBDクラス
//  - PBDによる弾性体シミュレーション
//-----------------------------------------------------------------------------
class ElasticPBD
{	
protected:
	// 形状データ
	int m_iNumVertices;					//!< 頂点数
	vector<glm::vec3> m_vCurPos;		//!< 現在の頂点位置
	vector<glm::vec3> m_vNewPos;		//!< 次のステップの頂点位置

	vector<glm::vec3> m_vVel;			//!< 頂点速度
	vector<bool> m_vFix;				//!< 頂点固定フラグ
	vector<float> m_vMass;				//!< 頂点質量(変形時の重み)

	rxPolygonsE m_poly;					//!< エッジ情報を含む3角形ポリゴン(元の頂点位置含む)
	int m_iNumEdge;						//!< エッジ数
	int m_iNumTris;						//!< ポリゴン数

	vector<float> m_vLengths;			//!< 頂点間の初期長さ(stretch constraint用)

	vector<int> m_vInEdge;				//!< 内部エッジならtrue

	// シミュレーションパラメータ
	glm::vec3 m_v3Min, m_v3Max;			//!< シミュレーション空間の大きさ
	glm::vec3 m_v3Gravity;				//!< 重力加速度ベクトル

	//海老沢追加------------------------------------------------------------------------
	int PBDstep;

	//衝突する球(1つのみ)
	glm::vec3 m_center;
	float m_rad;

	vector<glm::quat> m_vQuat;//現在の四元数を記憶
	vector<glm::quat> m_eOrgOmega;//初期のダルボーベクトルを記憶

	vector<glm::vec3> m_Fss;//辺ごとにかかる力

	vector<glm::vec3> m_LamdaSS;//Strech&ShearConstraintのラグランジュ乗数
	vector<glm::quat> m_LamdaBT;//Bend&TwistConstraintのラグランジュ乗数

	vector<glm::quat> m_tmpOrgOmega;

	vector<float> m_NewKss;//更新後のkssを保存
	vector<float> m_NewKbt;//更新後のkbtを保存

	vector<glm::vec3> m_vAcc;//SPH法から加速度を受け取るための配列
	vector<glm::vec3> m_vTang;//法線を保持する
	//-----------------------------------------------------------------------------------
	
	// 衝突処理関数の関数ポインタ
	CollisionFunc m_fpCollision;

public:
	int m_iNmax;						//!< 最大反復回数
	float m_fK;							//!< stiffnessパラメータ[0,1] 
	float m_fWind;						//!< 風の強さ

	bool m_bUseInEdge;					//!< 物体内部のエッジを制約計算に使うかどうかのフラグ
	int m_iObjectNo;					//!< オブジェクト番号

	float PBD_ks;                           //海老沢追加,初期設定の伸び剛性
	float PBD_kbt;                          //海老沢追加,初期設定の曲げ剛性

public:
	//! コンストラクタとデストラクタ
	ElasticPBD(int obj);
	~ElasticPBD();

	void Clear();

	void AddVertex(const glm::vec3 &pos, float mass);
	void AddEdge(int v0, int v1);
	void AddDarbouxVector(int i);

	// ストランド・メッシュ・四面体生成
	void GenerateStrand(glm::vec3 c1, glm::vec3 c2, int n);
	void GenerateCenterSpiral(glm::vec3 c1, glm::vec3 c2, int n);
	void GenerateNaturalSpiral(glm::vec3 c1, glm::vec3 c2, int n);
	void GenerateExampleRod(glm::vec3 c1, glm::vec3 c2, int n);

	//伸び剛性と曲げ剛性の設定
	void SetCoefficient(float ks, float kbt);;

	void Update(float dt);

	// OpenGL描画
	void Draw(int drw);

	// アクセスメソッド
	void SetSimulationSpace(glm::vec3 minp, glm::vec3 maxp){ m_v3Min = minp; m_v3Max = maxp; }
	void SetStiffness(float k){ m_fK = k; }
	void SetCollisionFunc(CollisionFunc func){ m_fpCollision = func; }

	glm::vec3 GetVertexPos(int i){ return ((i >= 0 && i < m_iNumVertices) ? m_vCurPos[i] : glm::vec3(0.0)); }
	//海老沢追加-----------------------------------------------
	int GetNumParticles(void) { return m_iNumVertices; }
	vector<glm::vec3> GetVertexPosPointer(void) { return m_vNewPos; }
	vector<glm::vec3> GetVertexVelPointer(void) { return m_vVel; }//初期速度は0なので実際には不要と考えられる
	vector<glm::vec3> GetTangPointer(void) { return m_vTang; }

	vector<float> GetVertexMass(void) { return m_vMass; }
	vector<float> GetInitial_Length(void) { return m_vLengths; }
	vector<float> GetInitial_Kss(void) { return m_NewKss; }
	vector<float> GetInitial_Kbt(void) { return m_NewKbt; }
	vector<glm::quat>GetQuat(void) { return m_vQuat; }
	vector<glm::quat>GetRestDarboux(void) { return m_eOrgOmega; }
	vector<bool>GetVertexFix(void) { return m_vFix; }

	void SetVertexAcc(vector<float> data, int start, int num);
	void SetVertexPos(vector<float> data, int start, int num);
	void SetVertexCurPos(vector<float> data, int start, int num);
	//---------------------------------------------------------

	// 固定点設定
	void FixVertex(int i, const glm::vec3 &pos);
	void FixVertex(int i);
	void UnFixVertex(int i);
	void UnFixAllVertex(void);
	bool IsFixed(int i) { return m_vFix[i]; }

	// 頂点選択(レイと頂点(球)の交差判定)
	int IntersectRay(const glm::vec3 &ray_origin, const glm::vec3 &ray_dir, float &t, float rad = 0.05);

	//海老沢追加 衝突する球を設定(とりあえず一つのみ)
	void SetCollisionSphere(glm::vec3 center, float rad);

	//海老沢追加 弾性体の全ての頂点を出力(デバック用)
	void PrintAllVertices(void);

	//SagFreeの処理をまとめる
	void SagFree(float ks, float kbt);

	//固定点と根元のy座標の判定
	bool judgePos(void);

protected:
	// PBD計算
	void calExternalForces(float dt);
	void genCollConstraints(float dt);
	//XPBDの制約
	void projectStretchingConstraint(float dt);
	void projectBendingConstraint(float dt);
	void projectCollisionConstraint(float dt);

	//PBDの制約
	void PBDConstraint(float ks);
	//伸びのみ
	void onlyStretch(float ks);

	//SagFreeの4段階計算
	void GlobalForceStep(void);
	void LocalForceStep(float ks);
	void GlobalTorqueStep(float kbt);
	void LocalTorqueStep(float kbt);

	//線形システムでGLobalTorqueStepを解く
	void VideoGlobalTorqueStep(float kbt);

	//GPUから持ってきたGlobalTorqueStep
	void GPU_GlobalTorqueStep(float kbt);

	void integrate(float dt);


	void clamp(glm::vec3 &pos) const
	{
		if(pos[0] < m_v3Min[0]) pos[0] = m_v3Min[0];
		if(pos[0] > m_v3Max[0]) pos[0] = m_v3Max[0];
		if(pos[1] < m_v3Min[1]) pos[1] = m_v3Min[1];
		if(pos[1] > m_v3Max[1]) pos[1] = m_v3Max[1];
		if(pos[2] < m_v3Min[2]) pos[2] = m_v3Min[2];
		if(pos[2] > m_v3Max[2]) pos[2] = m_v3Max[2];
	}

};


//! 四面体の体積計算
static inline float calVolume(const glm::vec3 &p0, const glm::vec3 &p1, const glm::vec3 &p2, const glm::vec3 &p3)
{
	return abs((1.0/6.0)*glm::dot(glm::cross(p1-p0, p2-p0), p3-p0));
}


#endif // _RX_PBD_STRAND_H_
