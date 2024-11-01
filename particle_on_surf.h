/*!
@file particle_on_surf.h

@brief 陰関数表面にパーティクルを配置

@author Makoto Fujisawa
@date   2013-06
*/


#ifndef _PARTICLE_ON_SURFACE_H_
#define _PARTICLE_ON_SURFACE_H_

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "utils.h"
#include "nnsearch.h"

#include "mc.h"

//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
typedef unsigned int uint;



//-----------------------------------------------------------------------------
// 陰関数等値面に粒子を配置するクラス
//  -参考: A. Witkin and P. Heckbert, "Using particles to sample and control implicit surfaces", SIGGRAPH1994.
//-----------------------------------------------------------------------------
class ParticleOnSurf
{
	struct SurfParticle
	{
		glm::vec3 Vel;
		glm::vec3 Ep;
		float Sigma;
		int Flag;

		SurfParticle() : Vel(glm::vec3(0.0)), Ep(glm::vec3(0.0)), Sigma(0.0), Flag(0) {}
	};

protected:
	uint m_num_particles;			//!< 粒子数

	float m_fParticleRadius;		//!< 粒子半径
	float m_fEffectiveRadius;		//!< 有効半径

	vector<float> m_vPos;			//!< 近傍探索ルーチンに渡すために位置だけ別管理
	vector<SurfParticle> m_vPs;	//!< 粒子情報
	__int64 m_iDim;					//!< m_vPosの1粒子が占める要素数(3 or 4) : C26451のwarningを避けるために__int64にしてある

									// Repulsionパラメータ
	float m_coef_alpha;	//!< repulsion amplitude
	float m_fEh;		//!< desired energy
	float m_fRho;		//!< 粒子のenergyをm_fEhに保つための係数
	float m_fPhi;		//!< 粒子位置を曲面上に保つための係数
	float m_fBeta;		//!< zero-divide防止用
	float m_fSh;		//!< desired repulsion radius (ユーザ指定)
	float m_fSmax;		//!< maximum repulsion radius

	float m_fGamma;	//!< equilibrium speed (σに対する倍数)
	float m_fNu;		//!< 粒子分裂のための係数
	float m_fDelta;	//!< 粒子削除のための係数

	void *m_pFuncPtr;
	glm::vec4 (*m_fpFunc)(float, float, float, void*);

	NNGrid *m_nn;			//!< 分割グリッドによる近傍探索
	vector< vector<int> > m_pnn;	//!< 近傍粒子

public:
	//! コンストラクタ
	ParticleOnSurf(int dim = 3)
	{
		m_num_particles = 0;
		m_iDim = static_cast<unsigned __int64>(dim);

		m_fParticleRadius = 0.1;
		m_fEffectiveRadius = 0.3;
		m_fSh = m_fParticleRadius;
		m_fSmax = 1.5*m_fSh;

		// 近傍探索セル
		m_nn = new NNGrid();

		m_coef_alpha = 6.0;
		m_fEh = 0.8*m_coef_alpha;
		m_fRho = 15.0;
		m_fPhi = 15.0;
		m_fBeta = 10.0;

		m_fGamma = 4.0;
		m_fNu = 0.2;
		m_fDelta = 0.7;

		m_fpFunc = 0;
		m_pFuncPtr = 0;
	}

	//! デストラクタ
	~ParticleOnSurf(){}

	/*!
	* 粒子初期化
	* @param[in] vrts 初期粒子配置
	* @param[in] rad  粒子半径
	* @param[in] minp,maxp 計算空間範囲
	* @param[in] func 表面を表す陰関数場(関数ポインタ)
	* @param[in] func_ptr funcの第一引数に渡す変数値(クラスオブジェクトなど)
	*/
	void Initialize(const vector<glm::vec3>& vrts, float rad, glm::vec3 minp, glm::vec3 maxp,
		glm::vec4(*func)(float, float, float, void*), void* func_ptr)
	{
		m_fParticleRadius = rad;
		m_fEffectiveRadius = 3.0*m_fParticleRadius;
		m_fSh = m_fParticleRadius;
		m_fSmax = 1.5*m_fSh;

		m_num_particles = (uint)vrts.size();

		m_fpFunc = func;
		m_pFuncPtr = func_ptr;

		// メモリ確保
		m_vPos.resize(m_num_particles*m_iDim, 0.0);
		m_vPs.resize(m_num_particles);


		// 初期粒子位置
		for(uint i = 0; i < m_num_particles; ++i){
			for(int j = 0; j < 3; ++j){
				m_vPos[m_iDim*i+j] = vrts[i][j];
			}
			m_vPs[i].Sigma = m_fSh;
		}

		// 分割セル設定
		m_nn->Setup(minp, maxp, m_fEffectiveRadius, m_num_particles);
		m_pnn.resize(m_num_particles);
	}
	/*!
	* 確保したメモリの解放
	*/
	void Finalize(void)
	{
		m_vPos.clear();
		m_vPs.clear();
	}


public:

	//! 粒子数
	int	GetNumParticles() const { return m_num_particles; }

	//! 粒子半径
	float GetParticleRadius(void){ return m_fParticleRadius; }

	//! 粒子データの取得
	float* GetPositionArray(void){ return &m_vPos[0]; }

	/*!
	* 粒子位置更新
	* @param[in] dt 更新用時間ステップ幅
	* @param[inout] num_iter 最大反復回数+実際の反復回数
	* @param[inout] eps 修飾判定用移動許容量+誤差
	*/
	inline void Update(float dt, int& num_iter, float& eps)
	{
		float v_avg = 0.0;
		int k;
		for(k = 0; k < num_iter; ++k){
			// 反発による速度の計算
			applyRepulsion2(dt);

			// 粒子位置の更新
			applyFloating(dt, v_avg);

			// 粒子追加/削除
			testFissionDeath();

			v_avg = sqrt(v_avg);
			if(v_avg < eps) break;
		}

		eps = v_avg;
		num_iter = k;
	}

protected:
	// 近傍取得
	/*!
	* 近傍粒子探索
	* @param[in] idx 探索中心粒子インデックス
	* @param[in] prts 粒子位置
	* @param[out] neighs 探索結果格納する近傍情報コンテナ
	* @param[in] h 有効半径
	*/
	void getNearestNeighbors(int idx, float* prts, vector<int>& neighs, float h = -1)
	{
		if(idx < 0 || idx >= (int)m_num_particles) return;
		glm::vec3 pos(prts[m_iDim*idx+0], prts[m_iDim*idx+1], prts[m_iDim*idx+2]);
		if(h < 0.0) h = m_fEffectiveRadius;
		m_nn->GetNN(pos, prts, m_num_particles, m_iDim*sizeof(float),neighs, h);
	}

	/*!
	* 近傍粒子探索
	* @param[in] idx 探索中心粒子インデックス
	* @param[out] neighs 探索結果格納する近傍情報コンテナ
	* @param[in] h 有効半径
	*/
	void getNearestNeighbors(glm::vec3 pos, vector<int>& neighs, float h)
	{
		if(h < 0.0) h = m_fEffectiveRadius;
		m_nn->GetNN(pos, &m_vPos[0], m_num_particles, m_iDim*sizeof(float), neighs, h);
	}

	// 分割セルに粒子を格納
	void setParticlesToCell(float* prts, int n, float h)
	{
		// 分割セルに粒子を登録
		m_nn->SetObjectToCell(prts, n, m_iDim*sizeof(float));

		// 近傍粒子探索
		if(h < 0.0) h = m_fEffectiveRadius;
		for(int i = 0; i < (int)m_num_particles; i++){
			m_pnn[i].clear();
			getNearestNeighbors(i, prts, m_pnn[i], h);
		}
	}
	void setParticlesToCell(void)
	{
		setParticlesToCell(&m_vPos[0], m_num_particles, m_fEffectiveRadius);
	}

protected:

	/*!
	* repulsion radius σ の更新
	*  - "Using particles to sample and control implicit surfaces"の式(10)〜(13)
	* @param[in] dt 更新用時間ステップ幅
	*/
	void updateSigma(float dt)
	{
		float h = m_fEffectiveRadius;

		for(uint i = 0; i < m_num_particles; ++i){
			glm::vec3 pos0(m_vPos[m_iDim*i+0], m_vPos[m_iDim*i+1], m_vPos[m_iDim*i+2]);
			float si = m_vPs[i].Sigma;

			float D = 0.0, Ds = 0.0;
			for(vector<int>::iterator itr = m_pnn[i].begin(); itr != m_pnn[i].end(); ++itr){
				int j = *itr;
				if(j < 0 || i == j) continue;

				glm::vec3 pos1(m_vPos[m_iDim*j+0], m_vPos[m_iDim*j+1], m_vPos[m_iDim*j+2]);

				glm::vec3 rij = pos1-pos0;

				float r = glm::length(rij);

				if(r <= h){
					float Eij = m_coef_alpha*exp(-r*r/(2.0*si*si));
					D += Eij;			// 現在の反発エネルギ(式(10)の上)
					Ds += r*r*Eij;		// エネルギDのσによる微分(式(13))
				}
			}
			Ds /= si*si*si;

			float Dv = -m_fRho*(D-m_fEh);// ターゲットエネルギに近づけるための線形フィードバック(式(10))
			float sv = Dv/(Ds+m_fBeta);	// σの変化量(式(12))

											// σの更新
			m_vPs[i].Sigma += sv*dt;
		}
	}

	/*!
	* 反発による速度の計算
	*  - 適応的な粒子の追加/削除のためのσの更新を含むバージョン
	*  - "Using particles to sample and control implicit surfaces"の4.3節，式(9)
	* @param[in] dt 更新用時間ステップ幅
	*/
	void applyRepulsion2(float dt)
	{
		float h = m_fEffectiveRadius;

		// 近傍探索セルへ粒子を格納 & 近傍探索
		setParticlesToCell();

		// repulsion radius の更新
		updateSigma(dt);

		for(uint i = 0; i < m_num_particles; ++i){
			glm::vec3 pos0(m_vPos[m_iDim*i+0], m_vPos[m_iDim*i+1], m_vPos[m_iDim*i+2]);
			float si = m_vPs[i].Sigma;

			glm::vec3 Ep(0.0);
			for(vector<int>::iterator itr = m_pnn[i].begin(); itr != m_pnn[i].end(); ++itr){
				int j = *itr;
				if(j < 0 || i == j) continue;

				glm::vec3 pos1(m_vPos[m_iDim*j+0], m_vPos[m_iDim*j+1], m_vPos[m_iDim*j+2]);

				glm::vec3 rij = pos1-pos0;

				float r = glm::length(rij);

				if(r <= h){
					float Eij = m_coef_alpha*exp(-r*r/(2.0*si*si));

					float sj = m_vPs[j].Sigma;
					float Eji = m_coef_alpha*exp(-r*r/(2.0*sj*sj));
					glm::vec3 rji = pos0-pos1;

					Ep += (rij/(si*si))*Eij-(rji/(sj*sj))*Eji;	// 式(9)
				}
			}
			Ep *= si*si;

			m_vPs[i].Ep = Ep;
		}
	}

	/*!
	* 反発による速度の計算
	*  - σの変更を含まないシンプルなバージョン
	*  - "Using particles to sample and control implicit surfaces"の4.1節
	* @param[in] dt 更新用時間ステップ幅
	*/
	void applyRepulsion(float dt)
	{
		float h = m_fEffectiveRadius;

		// 近傍探索セルへ粒子を格納 & 近傍探索
		setParticlesToCell();

		for(uint i = 0; i < m_num_particles; ++i){
			glm::vec3 pos0(m_vPos[m_iDim*i+0], m_vPos[m_iDim*i+1], m_vPos[m_iDim*i+2]);
			float si = m_vPs[i].Sigma;

			glm::vec3 Ep(0.0);
			for(vector<int>::iterator itr = m_pnn[i].begin(); itr != m_pnn[i].end(); ++itr){
				int j = *itr;
				if(j < 0 || i == j) continue;

				glm::vec3 pos1(m_vPos[m_iDim*j+0], m_vPos[m_iDim*j+1], m_vPos[m_iDim*j+2]);

				glm::vec3 rij = pos1-pos0;

				float r = glm::length(rij);

				if(r <= h){
					float Eij = m_coef_alpha*exp(-r*r/(2.0*si*si));	// 反発エネルギ
					Ep += rij*Eij;
				}
			}

			m_vPs[i].Ep = Ep;
		}
	}

	/*!
	* 粒子位置更新
	*  - 粒子が陰関数曲面上に載るように反発による速度を修正して，位置を更新
	* @param[in] dt 更新用時間ステップ幅
	* @param[out] v_avg 移動量の2乗平均値
	*/
	void applyFloating(float dt, float& v_avg)
	{
		v_avg = 0.0;
		if(!m_num_particles) return;

		float h = m_fEffectiveRadius;
		for(uint i = 0; i < m_num_particles; ++i){
			glm::vec3 p(m_vPos[m_iDim*i+0], m_vPos[m_iDim*i+1], m_vPos[m_iDim*i+2]);	// 現在の粒子座標
			glm::vec3 ve = m_vPs[i].Ep;	// 反発による速度

			glm::vec4 fv = m_fpFunc(p[0], p[1], p[2], m_pFuncPtr);	// 粒子座標における陰関数値とその勾配を取得
			glm::vec3 fx(fv[1], fv[2], fv[3]);	// 勾配
			float f = fv[0];				// 陰関数値

											// 粒子移流速度
			glm::vec3 v = ve-((glm::dot(fx, ve)+m_fPhi*f)/(glm::dot(fx, fx)))*fx;

			m_vPos[m_iDim*i+0] -= v[0]*dt;
			m_vPos[m_iDim*i+1] -= v[1]*dt;
			m_vPos[m_iDim*i+2] -= v[2]*dt;

			m_vPs[i].Vel = v;

			v_avg += glm::length2(v*dt);
		}
		v_avg /= m_num_particles;
	}

	/*!
	* 粒子の分裂と削除判定
	* @param[in] dt 更新用時間ステップ幅
	*/
	void testFissionDeath(void)
	{
		if(!m_num_particles) return;

		// 近傍探索セルへ粒子を格納 & 近傍探索
		setParticlesToCell();

		int num_remove = 0, num_fission = 0;
		float h = m_fEffectiveRadius;
		for(uint i = 0; i < m_num_particles; ++i){
			glm::vec3 p(m_vPos[m_iDim*i+0], m_vPos[m_iDim*i+1], m_vPos[m_iDim*i+2]);	// 粒子座標
			glm::vec3 v = m_vPs[i].Vel;	// 粒子速度
			float si = m_vPs[i].Sigma;	// 反発半径σ

			float vp = glm::length(v);

			// 粒子が平衡状態に近いかどうか
			if(vp < m_fGamma*si){
				// 反発エネルギDの計算
				float D = 0.0;
				for(vector<int>::iterator itr = m_pnn[i].begin(); itr != m_pnn[i].end(); ++itr){
					int j = *itr;
					if(j < 0 || i == j) continue;

					glm::vec3 pos1(m_vPos[m_iDim*j+0], m_vPos[m_iDim*j+1], m_vPos[m_iDim*j+2]);
					glm::vec3 rij = pos1-p;
					float r = glm::length(rij);
					if(r <= h){
						float Eij = m_coef_alpha*exp(-r*r/(2.0*si*si));
						D += Eij;			// 反発エネルギ
					}
				}
				float R = RX_FRAND();	// [0,1]の乱数

											// σが大きすぎる or エネルギーが適切でσが規定値以上 → 分裂(粒子追加)
				if((si > m_fSmax) || (D > m_fNu*m_fEh && si > m_fSh)){
					m_vPs[i].Flag = 1;	// 追加フラグをON
					num_fission++;
				}
				// σが小さすぎる and 乱数を使ったテストを通った → 削除
				else if((si < m_fDelta*m_fSh) && (R > si/(m_fDelta*m_fSh))){
					m_vPs[i].Flag = 2;	// 削除フラグをON
					num_remove++;
				}

			}

			// 表面から離れすぎた粒子も削除
			glm::vec4 fv = m_fpFunc(p[0], p[1], p[2], m_pFuncPtr);	// 粒子座標における陰関数値とその勾配を取得
			float f = fv[0];				// 陰関数値
			if(fabs(f) > 2.0*m_fSmax){
				m_vPs[i].Flag = 2;	// 削除フラグをON
				num_remove++;
			}
		}

		// 粒子削除
		if(num_remove){
			int cnt = 0;
			vector<SurfParticle>::iterator itr = m_vPs.begin();
			vector<float>::iterator jtr = m_vPos.begin();
			while(itr != m_vPs.end()){
				if(itr->Flag == 2){
					itr = m_vPs.erase(itr);
					jtr = m_vPos.erase(jtr, jtr+m_iDim);
					cnt++;
				} else{
					++itr;
					jtr += m_iDim;
				}
			}
			m_num_particles = (int)m_vPs.size();
			//cout << cnt << " particles are removed." << endl;
		}
	}


};




#endif // _PARTICLE_ON_SURFACE_H_

