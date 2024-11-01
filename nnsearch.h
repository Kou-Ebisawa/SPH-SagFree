/*!
@file nnsearch.h

@brief 矩形グリッド分割による近傍探索

@author Makoto Fujisawa
@date 2012-08,2022-03
*/

#ifndef _NNSEARCH_H_
#define _NNSEARCH_H_

#define NN_CELL_DRAW

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "utils.h"

#include <algorithm>
#include <vector>
#include <set>

#ifdef NN_CELL_DRAW
#include "gldraw.h"
#endif

using namespace std;

typedef unsigned int uint;


//-----------------------------------------------------------------------------
//! ハッシュ値によるソート用の構造体
//-----------------------------------------------------------------------------
struct HashSort
{
	uint hash;
	uint value;
};

/*!
 * ハッシュ値の比較関数
 * @param[in] left,right 比較する値
 * @return left < right
 */
inline bool LessHash(const HashSort &left, const HashSort &right)
{
	return left.hash < right.hash;
}


//-----------------------------------------------------------------------------
//! NNGridクラス - グリッド分割法による近傍探索
//-----------------------------------------------------------------------------
class NNGrid
{
	//! 探索用空間分割グリッドの各セル
	struct NNCell
	{
		uint* sorted_index;		//!< ハッシュ値でソートしたパーティクルインデックス
		uint* grid_hash;		//!< 各パーティクルのグリッドハッシュ値
		uint* starts;			//!< ソートリスト内の各セルのスタートインデックス
		uint* ends;				//!< ソートリスト内の各セルのエンドインデックス
		uint  num_cells;		//!< 総セル数
		NNCell() : sorted_index(0), grid_hash(0), starts(0), ends(0), num_cells(0) {}
	};


protected:
	// 空間分割格子
	NNCell m_cell;				//!< 探索用空間分割グリッド
	glm::ivec3 m_res;			//!< 格子の数
	glm::vec3 m_width;			//!< 格子一片の長さ

	glm::vec3 m_minp;			//!< 環境最小座標

	int m_sorted;					//!< 一度でもソートされたかどうかのフラグ

public:
	//! デフォルトコンストラクタ
	NNGrid()
	{
		m_sorted = 0;
	}

	//! デストラクタ
	~NNGrid()
	{
		if(m_cell.sorted_index) delete [] m_cell.sorted_index;
		if(m_cell.grid_hash) delete [] m_cell.grid_hash;
		if(m_cell.starts) delete [] m_cell.starts;
		if(m_cell.ends) delete [] m_cell.ends;
	}

public:
	/*!
	* 分割セルの初期設定
	* @param[in] vMin 環境の最小座標
	* @param[in] vMax 環境の最大座標
	* @param[in] h 影響半径
	*/
	void Setup(glm::vec3 vMin, glm::vec3 vMax, float h, int n)
	{
		if(h <= 0.0) return;

		glm::vec3 world_size = vMax-vMin;
		glm::vec3 world_origin = vMin;

		float max_axis = fabs(glm::max(world_size.x, world_size.y, world_size.z));

		float cell_width = h*1.01;	// 丸め誤差等を考慮して少しhより大きくする
		if(cell_width <= 0.0) return;

		for(int i = 0; i < 3; ++i){
			if(world_size[i] < h){
				m_res[i] = 1;
			}
			else{
				m_res[i] = (int)(world_size[i]/cell_width+0.5);
			}
		}

		// セル数を2の累乗にする場合
		//int d = (int)(log(max_axis/h)/log(2.0)+0.5);
		//int m = (int)(pow(2.0f, (float)d)+0.5);
		//float cell_width = max_axis/m;

		//for(int i = 0; i < 3; ++i){
		//	if(world_size[i] < h){
		//		m_res[i] = 1;
		//	}
		//	else{
		//		d = (int)(log(world_size[i]/cell_width)/log(2.0)+0.5);
		//		m_res[i] = (int)(pow(2.0, (float)d)+0.5);
		//	}
		//}

		m_width[0] = m_width[1] = m_width[2] = cell_width;
		m_minp = world_origin;

		m_cell.num_cells = m_res[0]*m_res[1]*m_res[2];

		//cout << "grid for nn search : " << endl;
		//cout << "  size   : " << m_res[0] << "x" << m_res[1] << "x" << m_res[2] << endl;
		//cout << "  num	: " << m_cell.num_cells << endl;
		//cout << "  origin : " << glm::to_string(m_minp) << endl;
		//cout << "  width  : " << glm::to_string(m_width) << endl;

		// 分割グリッド構造体の配列確保
		m_cell.sorted_index = new uint[n];
		m_cell.grid_hash = new uint[n];
		m_cell.starts = new uint[m_cell.num_cells];
		m_cell.ends = new uint[m_cell.num_cells];
	}

	/*!
	* オブジェクトを分割セルに格納
	*  - オブジェクトの属するグリッドハッシュを計算して格納する
	* @param[in] p 格納したい全オブジェクトの情報を記述した配列
	* @param[in] n パーティクル数
	* @param[in] stride オブジェクト情報配列の1要素の大きさ(byte数)
	*/
	void SetObjectToCell(void *p, uint n, int stride)
	{
		int mem_size1 = n*sizeof(uint);
		int mem_size2 = m_cell.num_cells*sizeof(uint);
		memset(m_cell.sorted_index, 0, mem_size1);
		memset(m_cell.grid_hash, 0, mem_size1);
		memset(m_cell.starts, 0xffffffff, mem_size2);
		memset(m_cell.ends, 0xffffffff, mem_size2);

		if(n == 0) return;

		// 各パーティクルのグリッドハッシュの計算
		float* pf = reinterpret_cast<float*>(p);
		for(uint i = 0; i < n; ++i){
			glm::vec3 pos(*pf, *(pf+1), *(pf+2));
			pf += stride/sizeof(float);

			// ハッシュ値計算
			uint hash = calGridHash(pos);

			m_cell.sorted_index[i] = i;
			m_cell.grid_hash[i] = hash;
		}

		// グリッドハッシュでソート
		vector<HashSort> hash_and_value;
		hash_and_value.resize(n);
		for(uint i = 0; i < n; ++i){
			hash_and_value[i].hash = m_cell.grid_hash[i];
			hash_and_value[i].value  = m_cell.sorted_index[i];
		}
		std::sort(hash_and_value.begin(), hash_and_value.end(), LessHash);
		for(uint i = 0; i < n; ++i){
			m_cell.sorted_index[i] = hash_and_value[i].value;
			m_cell.grid_hash[i] = hash_and_value[i].hash;
		}

		// パーティクル配列をソートされた順番に並び替え，
		// 各セルの始まりと終わりのインデックスを検索
		for(uint i = 0; i < n; i++){
			int hash = m_cell.grid_hash[i];

			if(i == 0){
				m_cell.starts[hash] = i;
				m_cell.ends[hash] = i;
			}
			else{
				int prev_hash = m_cell.grid_hash[i-1];

				if(i == 0 || hash != prev_hash){
					m_cell.starts[hash] = i;
					if(i > 0){
						m_cell.ends[prev_hash] = i;
					}
				}

				if(i == n-1){
					m_cell.ends[hash] = i+1;
				}
			}
		}
	}


	/*!
	* 近傍粒子探索
	* @param[in] pos 探索中心
	* @param[in] p パーティクル位置
	* @param[out] neighs 探索結果格納する近傍情報コンテナ
	* @param[in] h 有効半径
	*/
	void GetNN(glm::vec3 pos, void *p, uint n, int stride, vector<int> &neighs, float h)
	{
		// 分割セルインデックスの算出
		int x = (pos[0]-m_minp[0])/m_width[0];
		int y = (pos[1]-m_minp[1])/m_width[1];
		int z = (pos[2]-m_minp[2])/m_width[2];
		x = (x < 0 ? 0 : (x >= m_res[0] ? m_res[0]-1 : x));
		y = (y < 0 ? 0 : (y >= m_res[1] ? m_res[1]-1 : y));
		z = (z < 0 ? 0 : (z >= m_res[2] ? m_res[2]-1 : z));

		int numArdGrid = (int)(h/m_width[0])+1;
		for(int k = -numArdGrid; k <= numArdGrid; ++k){
			for(int j = -numArdGrid; j <= numArdGrid; ++j){
				for(int i = -numArdGrid; i <= numArdGrid; ++i){
					int i1 = x+i, j1 = y+j, k1 = z+k;
					if(i1 < 0 || i1 >= m_res[0] || j1 < 0 || j1 >= m_res[1] || k1 < 0 || k1 >= m_res[2]){
						continue;
					}
					getNeighborsInCell(pos, p, stride, i1, j1, k1, neighs, h);
				}
			}
		}
	}
	/*!
	* 近傍粒子探索(総当たり,デバッグ用)
	* @param[in] pos 探索中心
	* @param[in] p パーティクル位置
	* @param[out] neighs 探索結果格納する近傍情報コンテナ
	* @param[in] h 有効半径
	*/
	void GetNN_Direct(glm::vec3 pos0, void *p, uint n, int stride, vector<int> &neighs, float h)
	{
		float h2 = h*h;

		float* pf = reinterpret_cast<float*>(p);
		for(uint i = 0; i < n; i++){
			glm::vec3 pos1(*pf, *(pf+1), *(pf+2));
			pf += stride/sizeof(float);

			float d2 = glm::length2(pos0-pos1);
			if(d2 <= h2){
				neighs.push_back(i);
			}
		}
	}

	// OpenGL描画用セル情報の取得(デバッグ用)
	glm::ivec3 GetCellNum(void) const { return m_res; }
	int GetCell(int i, int j, int k, glm::vec3 *vrts, int *edgs, int voffset) const 
	{
		glm::vec3 minp = m_minp+glm::vec3(i*m_width[0], j*m_width[1], k*m_width[2]);
		*vrts++ = minp;
		*vrts++ = minp+glm::vec3(m_width[0], 0, 0);
		*vrts++ = minp+glm::vec3(m_width[0], m_width[1], 0);
		*vrts++ = minp+glm::vec3(0, m_width[1], 0);
		*vrts++ = minp+glm::vec3(0, 0, m_width[2]);
		*vrts++ = minp+glm::vec3(m_width[0], 0, m_width[2]);
		*vrts++ = minp+glm::vec3(m_width[0], m_width[1], m_width[2]);
		*vrts++ = minp+glm::vec3(0, m_width[1], m_width[2]);
		//const int faces[6][4] = { {0, 4, 7, 3}, {1, 2, 6, 5}, {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 7, 6, 2} };
		const int idxs[12][2] = { {0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6}, {6,7}, {7,4}, {0,4}, {1,5}, {2,6}, {3,7} };
		for(int e = 0; e < 12; ++e)
			for(int v = 0; v < 2; ++v)
				*edgs++ = idxs[e][v]+voffset;

		uint ghash = calGridHash(i, j, k);
		return (m_cell.starts[ghash] != 0xffffffff) ? (m_cell.ends[ghash]-m_cell.starts[ghash]) : 0;
	}

protected:
	/*!
	* 分割セル内のオブジェクトから近傍を検出
	* @param[in] pos 探索中心
	* @param[in] p パーティクル位置
	* @param[in] gi,gj,gk 対象分割セル
	* @param[out] neighs 探索結果格納する近傍情報コンテナ
	* @param[in] h 有効半径
	*/
	void getNeighborsInCell(glm::vec3 pos, void *p, int stride, int gi, int gj, int gk, vector<int> &neighs, float h)
	{
		float h2 = h*h;
		float* pf = reinterpret_cast<float*>(p);

		uint grid_hash = calGridHash(gi, gj, gk);

		uint start_index = m_cell.starts[grid_hash];
		if(start_index != 0xffffffff){	// セルが空でないかのチェック
			uint end_index = m_cell.ends[grid_hash];
			for(uint j = start_index; j < end_index; ++j){
				uint idx = m_cell.sorted_index[j];

				float* pfx = pf+idx*stride/sizeof(float);
				glm::vec3 xij = pos-glm::vec3(*pfx, *(pfx+1), *(pfx+2));

				float d2 = glm::length2(xij);
				if(d2 <= h2){
					neighs.push_back(idx);
				}
			}
		}
	}

	// グリッドハッシュの計算
	uint calGridHash(int x, int y, int z) const 
	{
		x = (x < 0 ? 0 : (x >= m_res[0] ? m_res[0]-1 : x));
		y = (y < 0 ? 0 : (y >= m_res[1] ? m_res[1]-1 : y));
		z = (z < 0 ? 0 : (z >= m_res[2] ? m_res[2]-1 : z));
		return z*m_res[0]*m_res[1]+y*m_res[0]+x;
	}
	uint calGridHash(glm::vec3 pos) const 
	{
		pos -= m_minp;

		// 分割セルインデックスの算出
		int x = pos[0]/m_width[0];
		int y = pos[1]/m_width[1];
		int z = pos[2]/m_width[2];
		return calGridHash(x, y, z);
	}

#ifdef NN_CELL_DRAW
public:
	/*!
	 * 探索用セルの描画
	 * @param[in] i,j,k グリッド上のインデックス
	 */
	inline void DrawCell(int i, int j, int k)
	{
		glPushMatrix();
		glTranslated(m_minp[0], m_minp[1], m_minp[2]);
		glTranslatef((i+0.5)*m_width[0], (j+0.5)*m_width[1], (k+0.5)*m_width[2]);
		glScalef(m_width[0], m_width[1], m_width[2]);
		DrawCubeVBO();
		glPopMatrix();
	}

	/*!
	* 探索用グリッドの描画
	* @param[in] col 粒子が含まれるセルの色
	*/
	inline void DrawCells(glm::vec3 col)
	{
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		// 粒子 or ポリゴンを含む全セルの描画
		glPushMatrix();
		for(int i = 0; i < m_res[0]; ++i){
			for(int j = 0; j < m_res[1]; ++j){
				for(int k = 0; k < m_res[2]; ++k){
					bool disp = false;
					uint grid_hash = calGridHash(i, j, k);
					if(m_cell.starts[grid_hash] != 0xffffffff){
						glColor3fv(glm::value_ptr(col));
						DrawCell(i, j, k);
					}
				}
			}
		}
		glPopMatrix();

		// 探索グリッド全体の大きさ
		glPushMatrix();
		glm::vec3 len(m_width[0]*m_res[0], m_width[1]*m_res[1], m_width[2]*m_res[2]);
		glm::vec3 cen = m_minp+0.5f*len;
		glTranslatef(cen[0], cen[1], cen[2]);
		glScalef(len[0], len[1], len[2]);
		glColor3d(0.5, 1.0, 0.5);
		DrawCubeVBO();
		glPopMatrix();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
#endif
};



#endif // #ifndef _NNSEARCH_H_