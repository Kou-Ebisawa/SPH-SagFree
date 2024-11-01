/*!
  @file pdist.h

  @brief 粒子配置生成用の関数群

  @author Makoto Fujisawa
  @date 2020-06,2022-07
*/

#ifndef _PARTICLE_DISTRIBUTION_H_
#define _PARTICLE_DISTRIBUTION_H_

#include <vector>
#include <deque>

#include "utils.h"

#include "mc.h"
#include "particle_on_surf.h"


//! 粒子配置用構造体
struct rxPdist
{
	std::vector<glm::vec3> pdist;	//!< 追加粒子配置
	std::deque<int> steps;	//!< 追加ステップ(複数あり)
	glm::vec3 vel;			//!< 追加粒子の初期速度
	int attr;			//!< 追加粒子の属性
};

/*!
 * ボックス形状内に並べられた粒子配置の作成
 * @param[in] cen ボックスの中心座標
 * @param[in] ext ボックスの大きさ(各辺の長さの1/2)
 * @param[in] vel 初期速度
 * @param[in] spacing 粒子配置間隔
 * @param[in] attr 粒子属性(初期値0)
 * @param[in] step 追加ステップ数(初期値0)
 */
inline static rxPdist MakeBox(glm::vec3 cen, glm::vec3 ext, glm::vec3 vel, float spacing, int attr = 0, int step = 0)
{
	rxPdist pd;

	// 配置初期位置(原点)とランダムな揺らぎ量の設定
	int sx = (int)(ext[0]/spacing)-1;
	int sy = (int)(ext[1]/spacing)-1;
	int sz = (int)(ext[2]/spacing)-1;
	srand((unsigned)time(NULL));
	float jitter = spacing*0.001f;

	// 3次元格子状に粒子を配置
	int count = 0;
	for(int z = -sz; z <= sz; ++z){
		for(int y = -sy; y <= sy; ++y){
			for(int x = -sx; x <= sx; ++x){
				glm::vec3 p = cen+glm::vec3(x, y, z)*spacing;
				for(int l = 0; l < 3; ++l) p[l] += (RX_FRAND()*2.0f-1.0f)*jitter;
				pd.pdist.push_back(p);
				count++;
			}
		}
	}

	// 粒子速度・属性の設定
	pd.vel = vel;
	pd.attr = attr;

	// 追加タイミング(ステップ数)
	pd.steps.push_back(step);

	return pd;
}


/*!
 * ボックス形状内に並べられた粒子配置の作成
 * @param[in] cen ボックスの中心座標
 * @param[in] ext ボックスの大きさ(各辺の長さの1/2)
 * @param[in] vel 初期速度
 * @param[in] spacing 粒子配置間隔
 * @param[in] attr 粒子属性(初期値0)
 * @param[in] step 追加ステップ数(初期値0)
 */
inline static rxPdist MakeOne(glm::vec3 cen, glm::vec3 vel, int attr = 0, int step = 0)
{
	rxPdist pd;
	pd.pdist.push_back(cen);

	// 粒子速度・属性の設定
	pd.vel = vel;
	pd.attr = attr;

	// 追加タイミング(ステップ数)
	pd.steps.push_back(step);

	return pd;
}


/*!
 * 中空のボックス形状内に並べられた粒子配置の作成
 * @param[in] cen ボックスの中心座標
 * @param[in] ext0 ボックスの内部の大きさ(各辺の長さの1/2)
 * @param[in] ext1 ボックスの外部の大きさ(各辺の長さの1/2)
 * @param[in] vel 初期速度
 * @param[in] spacing 粒子配置間隔
 * @param[in] attr 粒子属性(初期値0)
 * @param[in] step 追加ステップ数(初期値0)
 */
inline static rxPdist MakeHollowBox(glm::vec3 cen, glm::vec3 ext0, glm::vec3 ext1, glm::vec3 vel, float spacing, int attr = 0, int step = 0)
{
	rxPdist pd;

	// 配置初期位置(原点)とランダムな揺らぎ量の設定
	int sx = (int)(ext1[0]/spacing)-1;
	int sy = (int)(ext1[1]/spacing)-1;
	int sz = (int)(ext1[2]/spacing)-1;
	srand((unsigned)time(NULL));
	float jitter = spacing*0.001f;

	// 3次元格子状に粒子を配置
	int count = 0;
	for(int z = -sz; z <= sz; ++z){
		for(int y = -sy; y <= sy; ++y){
			for(int x = -sx; x <= sx; ++x){
				glm::vec3 p = glm::vec3(x, y, z)*spacing;
				if(fabs(p[0]) < ext0[0] && fabs(p[1]) < ext0[1] && fabs(p[2]) < ext0[2]) continue;

				p += cen;
				for(int l = 0; l < 3; ++l) p[l] += (RX_FRAND()*2.0f-1.0f)*jitter;
				pd.pdist.push_back(p);
				count++;
			}
		}
	}

	// 粒子速度・属性の設定
	pd.vel = vel;
	pd.attr = attr;

	// 追加タイミング(ステップ数)
	pd.steps.push_back(step);

	return pd;
}


/*!
 * 球形状内に並べられた粒子配置の作成
 * @param[in] cen 球の中心座標
 * @param[in] ext 球の半径
 * @param[in] vel 初期速度
 * @param[in] spacing 粒子配置間隔
 * @param[in] attr 粒子属性(初期値0)
 * @param[in] step 追加ステップ数(初期値0)
 */
inline static rxPdist MakeSphere(glm::vec3 cen, float r, glm::vec3 vel, float spacing, int attr, int step)
{
	rxPdist pd;

	// ランダムな揺らぎ量の設定
	srand((unsigned)time(NULL));
	float jitter = spacing*0.01f;

	// グリッド数基準の半径
	int sr = (int)(r/spacing)+1;

	// 3次元格子状(球内)に粒子を配置
	int count = 0;
	for(int z = -sr; z <= sr; ++z){
		for(int y = -sr; y <= sr; ++y){
			for(int x = -sr; x <= sr; ++x){
				float dx[3];
				dx[0] = x*spacing;
				dx[1] = y*spacing;
				dx[2] = z*spacing;
				float l = sqrtf(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
				if((l <= r)) {
					glm::vec3 p(0.0);
					for(int j = 0; j < 3; ++j){
						p[j] = cen[j]+dx[j]+(RX_FRAND()*2.0f-1.0f)*jitter;
					}
					pd.pdist.push_back(p);
					count++;
				}
			}
		}
	}


	// 粒子速度・属性の設定
	pd.vel = vel;
	pd.attr = attr;

	// 追加タイミング(ステップ数)
	pd.steps.push_back(step);

	return pd;
}

/*!
 * 液体粒子の流入ラインの設定
 * @param[in] cen 球の中心座標
 * @param[in] ext 球の半径
 * @param[in] vel 初期速度
 * @param[in] spacing 粒子配置間隔
 * @param[in] attr 粒子属性(初期値0)
 * @param[in] step 追加ステップ数(初期値0)
 */
inline static rxPdist MakeInletLine(glm::vec3 pos1, glm::vec3 pos2, glm::vec3 vel, glm::vec3 up, int accum, int span, float spacing, int attr, int steps)
{
	rxPdist pd;

	glm::vec3 rel;	// 端点間の相対位置ベクトル
	float l = 0.0;		// 端点間の距離
	for(int i = 0; i < 3; ++i){
		rel[i] = pos2[i]-pos1[i];
		l += rel[i]*rel[i];
	}
	l = sqrt(l);

	int n = l/(spacing);	// 並べる粒子の数
	if(!n) return pd;

	// ランダムな揺らぎ量の設定
	srand((unsigned)time(NULL));
	float jitter = spacing*0.01f;

	// 3次元格子状(球内)に粒子を配置
	int count = 0;
	for(int j = 0; j < accum; ++j){
		for(int i = 0; i < n; ++i){
			glm::vec3 p;
			for(int k = 0; k < 3; ++k){
				p[k] = pos1[k]+rel[k]*(i+0.5)/n+(RX_FRAND()*2.0f-1.0f)*jitter;
			}
			pd.pdist.push_back(p);
			count++;
		}
		pos1 += up*spacing;
	}


	// 粒子速度・属性の設定
	pd.vel = vel;
	pd.attr = attr;

	// 追加タイミング(ステップ数)
	for(int i = 0; i < steps; i += span){
		pd.steps.push_back(i);
	}

	return pd;
}



/*!
 * 表面粒子生成
 *  - 陰関数場を入力として，MC法で初期粒子を生成，粒子間にバネを配置して位置を調整する
 *  - -参考: A. Witkin and P. Heckbert, "Using particles to sample and control implicit surfaces", SIGGRAPH1994.
 * @param[in] func 陰関数値取得関数ポインタ(第1引数にユーザー変数,第2-4引数に3次元座標を渡し，法線と陰関数値(Vec4)を返す関数)
 * @param[in] func_ptr funcの第1引数に渡すユーザ変数
 * @param[in] minp,maxp 計算領域
 * @param[in] rad 粒子半径
 * @param[out] ppos 生成した粒子を格納する配列(中身はすべて削除されて新しく情報が格納されるので注意)
 * @param[in] dim pposへ位置を格納するときの次元数
 *
 */
static int GenerateParticlesOnSurf(glm::vec4 (*func)(float, float, float, void*), void* func_ptr, glm::vec3 minp, glm::vec3 maxp,
								   float rad, vector<float> &ppos, int dim = 3)
{
	int num_particles = 0;

	vector<glm::vec3> vrts, nrms;
	vector<int> tris;

	minp -= glm::vec3(6.0*rad);
	maxp += glm::vec3(6.0*rad);

	// 初期配置粒子計算:メッシュ生成格子の決定
	float h = 2.0f*rad;
	float thr = 0.0f;
	int n[3];
	for(int i = 0; i < 3; ++i){
		n[i] = (int)((maxp[i]-minp[i])/h)+1;
	}

	// 初期配置粒子計算:Marching Cubesで等値面メッシュ生成
	//  -> 必要なのは頂点座標だけなのでメッシュインデックス作成を省けば少しは早くなりそう
	MCMesh mc;
	mc.CreateMesh(func, func_ptr, minp, h, n, thr, vrts, nrms, tris);

	// 頂点数
	num_particles = (int)vrts.size();

	if(num_particles){
		cout << num_particles << " particles are generated for the boundary." << endl;

		ParticleOnSurf sp;

		// 粒子初期配置
		sp.Initialize(vrts, rad, minp, maxp, func, func_ptr);

		// 粒子位置修正時の最大反復計算回数と許容誤差
		int iter = 50;
		float eps = 1.0e-4;

		// 粒子位置修正
		sp.Update(0.01, iter, eps);

		cout << iter << " iterations and the error value is " << eps << endl;

		// 粒子数の取得(位置修正時に必要に応じて粒子が削除されることがあるので必ず再取得しておくこと)
		num_particles = sp.GetNumParticles();

		// 固体粒子情報を格納するメモリ領域の確保(nmaxが-1の場合は新規確保，そうでなければ確保済み配列に追加)
		if(!ppos.empty()) ppos.clear();
		ppos.resize(dim*num_particles, 0.0f);

		// 固体粒子情報の取り出し
		float* spos = sp.GetPositionArray();
		for(int i = 0; i < num_particles; ++i){
			ppos[dim*(i)+0] = spos[dim*i+0];
			ppos[dim*(i)+1] = spos[dim*i+1];
			ppos[dim*(i)+2] = spos[dim*i+2];
		}
	}

	return num_particles;
}



#endif // #ifdef _RX_PARTICLE_DISTRIBUTION_H_
