﻿/*!
  @file mc.h
	
  @brief 陰関数表面からのポリゴン生成(MC法)
	
	http://local.wasp.uwa.edu.au/~pbourke/geometry/polygonise/
 
  @author Raghavendra Chandrashekara (basesd on source code
			provided by Paul Bourke and Cory Gene Bloyd)
  @date   2010-03,2022-07
*/


#ifndef _MC_MESH_H_
#define _MC_MESH_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
// C標準
#include <cstdlib>

// STL
#include <map>
#include <vector>
#include <string>

#include <iostream>

#include "utils.h"
#include "rx_mesh.h"

#include "mc_tables.h"



//-----------------------------------------------------------------------------
// 名前空間
//-----------------------------------------------------------------------------
using namespace std;


//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
typedef unsigned int uint;


//-----------------------------------------------------------------------------
// MCMeshクラス
//-----------------------------------------------------------------------------
class MCMesh
{
	struct RxVertexID
	{
		uint newID;
		float x, y, z;
	};

	typedef std::map<uint, RxVertexID> ID2VertexID;

	struct RxTriangle
	{
		uint pointID[3];
	};

	typedef std::vector<RxTriangle> RxTriangleVector;

	struct RxScalarField
	{
		uint iNum[3];
		glm::vec3 fWidth;
		glm::vec3 fMin;
	};

protected:
	// メンバ変数
	uint m_nVertices;	//!< 等値面メッシュの頂点数
	uint m_nNormals;	//!< 等値面メッシュの頂点法線数(作成されていれば 法線数=頂点数)
	uint m_nTriangles;	//!< 等値面メッシュの三角形ポリゴン数

	ID2VertexID m_i2pt3idVertices;			//!< 等値面を形成する頂点のリスト
	RxTriangleVector m_trivecTriangles;		//!< 三角形ポリゴンを形成する頂点のリスト

	RxScalarField m_Grid;					//!< 分割グリッド情報

	// 陰関数値(スカラー値)取得用変数(どちらかのみ用いる)
	const float* m_ptScalarField;				//!< スカラー値を保持するサンプルボリューム
	float (*m_fpScalarFunc)(float, float, float, void*);	//!< スカラー値を返す関数ポインタ
	void *m_pScalarFuncPtr;

	float m_tIsoLevel;							//!< 閾値

	bool m_bValidSurface;					//!< メッシュ生成成功の可否


public:
	// コンストラクタ
	MCMesh()
	{
		// MARK:コンストラクタ
		m_Grid.fMin = glm::vec3(0.0);
		m_Grid.fWidth = glm::vec3(0.0);
		m_Grid.iNum[0] = 0;
		m_Grid.iNum[1] = 0;
		m_Grid.iNum[2] = 0;

		m_nTriangles = 0;
		m_nNormals = 0;
		m_nVertices = 0;

		m_ptScalarField = NULL;
		m_fpScalarFunc = 0;
		m_pScalarFuncPtr = 0;
		m_tIsoLevel = 0;
		m_bValidSurface = false;
	}


	// デストラクタ
	~MCMesh()
	{
		DeleteSurface();
	}

	
	/*!
	* 陰関数場から三角形メッシュを生成
	* @param[in] func 陰関数値取得用関数ポインタ
	* @param[in] func_ptr 陰関数値取得用関数の第四引数に渡すポインタ(funcがグローバル関数なら0,メンバ関数ならクラスのインスタンス)
	* @param[in] minp グリッドの最小座標
	* @param[in] h グリッドの幅
	* @param[in] n[3] グリッド数(x,y,z)
	* @param[in] threshold しきい値(陰関数値がこの値のところをメッシュ化)
	* @param[out] vrts 頂点座標
	* @param[out] nrms 頂点法線
	* @param[out] tris メッシュ
	* @retval true  メッシュ生成成功
	* @retval false メッシュ生成失敗
	*/
	bool CreateMesh(float (*func)(float, float, float, void*), void* func_ptr, glm::vec3 minp, float h, int n[3], float threshold, 
					vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, vector<int> &tris)
	{
		if(func == NULL) return false;

		int nx, ny, nz;
		nx = n[0]+1; ny = n[1]+1; nz = n[2]+1;

		// 陰関数場 → MC法のグリッドセル端点毎の陰関数値を格納した配列に変換
		vector<float> field(nx*ny*nz);
		for(int k = 0; k < nz; ++k){
			for(int j = 0; j < ny; ++j){
				for(int i = 0; i < nx; ++i){
					int idx = k*nx*ny+j*nx+i;
					glm::vec3 pos = minp+glm::vec3(i, j, k)*h;

					float val = func(pos[0], pos[1], pos[2], func_ptr);
					field[idx] = val;
				}
			}
		}

		// 表面メッシュ生成
		generateSurface(&field[0], minp, h, n, threshold, vrts, nrms, tris);

		return IsSurfaceValid();
	}
	/*!
	* 陰関数場から三角形メッシュを生成
	* - 値と勾配を含むglm::vec4を返す関数を用いるバージョン
	* - MC法は勾配は使わないので使われるのは0番目の要素のみ
	*/
	bool CreateMesh(glm::vec4 (*func)(float, float, float, void*), void* func_ptr, glm::vec3 minp, float h, int n[3], float threshold, 
					vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, vector<int> &tris)
	{
		if(func == NULL) return false;

		int nx, ny, nz;
		nx = n[0]+1; ny = n[1]+1; nz = n[2]+1;

		// 陰関数場 → MC法のグリッドセル端点毎の陰関数値を格納した配列に変換
		vector<float> field(nx*ny*nz);
		for(int k = 0; k < nz; ++k){
			for(int j = 0; j < ny; ++j){
				for(int i = 0; i < nx; ++i){
					int idx = k*nx*ny+j*nx+i;
					glm::vec3 pos = minp+glm::vec3(i, j, k)*h;

					float val = func(pos[0], pos[1], pos[2], func_ptr)[0];
					field[idx] = val;
				}
			}
		}

		// 表面メッシュ生成
		generateSurface(&field[0], minp, h, n, threshold, vrts, nrms, tris);

		return IsSurfaceValid();
	}

	/*!
	* 陰関数場から三角形メッシュを生成
	*  - min_p,h,n[3]はメッシュ生成用グリッドのパラメータであるとともに，入力されるフィールド値のパラメータでもある
	*    (つまり，fieldの大きさはn[0]*n[1]*n[2])
	* @param[in] field サンプルボリューム
	* @param[in] minp グリッドの最小座標
	* @param[in] h グリッドの幅
	* @param[in] n[3] グリッド数(x,y,z)
	* @param[in] threshold しきい値(陰関数値がこの値のところをメッシュ化)
	* @param[out] vrts 頂点座標
	* @param[out] nrms 頂点法線
	* @param[out] tris メッシュ
	* @retval true  メッシュ生成成功
	* @retval false メッシュ生成失敗
	*/
	bool CreateMeshFromField(float *field, glm::vec3 minp, double h, int n[3], float threshold, 
							 vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, vector<int> &tris)
	{
		if(field){
			// 表面メッシュ生成
			n[0] -= 1; n[1] -= 1; n[2] -= 1;
			generateSurface(field, minp, h, n, threshold, vrts, nrms, tris);
			return IsSurfaceValid();
		}
		return false;
	}

	//! 等値面作成が成功したらtrueを返す
	bool IsSurfaceValid() const { return m_bValidSurface; }

	//! 作成した等値面メッシュの破棄
	void DeleteSurface()
	{
		m_Grid.fWidth[0] = 0;
		m_Grid.fWidth[1] = 0;
		m_Grid.fWidth[2] = 0;
		m_Grid.iNum[0] = 0;
		m_Grid.iNum[1] = 0;
		m_Grid.iNum[2] = 0;

		m_nTriangles = 0;
		m_nNormals = 0;
		m_nVertices = 0;

		//m_vVertices.clear();
		//m_normals.clear();
		//m_vTriangles.clear();

		m_ptScalarField = NULL;
		m_tIsoLevel = 0;
		m_bValidSurface = false;
	}

	/*!
	 * メッシュ化に用いたグリッドの大きさ
	 * @param[out] fVolLength* グリッドの大きさ
	 * @return メッシュ生成されていれば1, そうでなければ-1
	 */
	int GetVolumeLengths(double& fVolLengthX, double& fVolLengthY, double& fVolLengthZ)
	{
		if(IsSurfaceValid()){
			fVolLengthX = m_Grid.fWidth[0]*m_Grid.iNum[0];
			fVolLengthY = m_Grid.fWidth[1]*m_Grid.iNum[1];
			fVolLengthZ = m_Grid.fWidth[2]*m_Grid.iNum[2];
			return 1;
		}
		else
			return -1;
	}
	// 作成したメッシュの情報
	uint GetNumVertices(void) const { return m_nVertices; }
	uint GetNumTriangles(void) const { return m_nTriangles; }
	uint GetNumNormals(void) const { return m_nNormals; }


protected:


	/*!
	* メッシュ生成
	* @param[in] sf 分割グリッド情報
	* @param[in] field サンプルボリューム
	* @param[in] threshold 閾値
	* @param[out] vrts メッシュ頂点
	* @param[out] nrms メッシュ頂点法線
	* @param[out] tris メッシュ幾何情報(頂点接続情報)
	*/
	void generateSurface(float *field, glm::vec3 minp, float h, int n[3], float threshold, 
						 vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, vector<int> &tris)
	{
		// MARK:GenerateSurface
		if(m_bValidSurface){
			DeleteSurface();
		}

		m_ptScalarField = field;
		m_tIsoLevel = threshold;
		m_Grid.iNum[0] = n[0];
		m_Grid.iNum[1] = n[1];
		m_Grid.iNum[2] = n[2];
		m_Grid.fWidth = glm::vec3(h);
		m_Grid.fMin = minp;

		uint slice0 = (m_Grid.iNum[0] + 1);
		uint slice1 = slice0*(m_Grid.iNum[1] + 1);

		// 等値面の生成
		for(uint z = 0; z < m_Grid.iNum[2]; ++z){
			for(uint y = 0; y < m_Grid.iNum[1]; ++y){
				for(uint x = 0; x < m_Grid.iNum[0]; ++x){
					// グリッド内の頂点配置情報テーブル参照用インデックスの計算
					uint tableIndex = 0;
					if(m_ptScalarField[z*slice1 + y*slice0 + x] < m_tIsoLevel)
						tableIndex |= 1;
					if(m_ptScalarField[z*slice1 + (y+1)*slice0 + x] < m_tIsoLevel)
						tableIndex |= 2;
					if(m_ptScalarField[z*slice1 + (y+1)*slice0 + (x+1)] < m_tIsoLevel)
						tableIndex |= 4;
					if(m_ptScalarField[z*slice1 + y*slice0 + (x+1)] < m_tIsoLevel)
						tableIndex |= 8;
					if(m_ptScalarField[(z+1)*slice1 + y*slice0 + x] < m_tIsoLevel)
						tableIndex |= 16;
					if(m_ptScalarField[(z+1)*slice1 + (y+1)*slice0 + x] < m_tIsoLevel)
						tableIndex |= 32;
					if(m_ptScalarField[(z+1)*slice1 + (y+1)*slice0 + (x+1)] < m_tIsoLevel)
						tableIndex |= 64;
					if(m_ptScalarField[(z+1)*slice1 + y*slice0 + (x+1)] < m_tIsoLevel)
						tableIndex |= 128;

					if(edgeTable[tableIndex] != 0){
						// エッジ上の頂点算出
						if(edgeTable[tableIndex] & 8){
							RxVertexID pt = calIntersection(x, y, z, 3);
							uint id = getEdgeID(x, y, z, 3);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
						if(edgeTable[tableIndex] & 1){
							RxVertexID pt = calIntersection(x, y, z, 0);
							uint id = getEdgeID(x, y, z, 0);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}
						if(edgeTable[tableIndex] & 256){
							RxVertexID pt = calIntersection(x, y, z, 8);
							uint id = getEdgeID(x, y, z, 8);
							m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
						}

						if(x == m_Grid.iNum[0] - 1){
							if(edgeTable[tableIndex] & 4){
								RxVertexID pt = calIntersection(x, y, z, 2);
								uint id = getEdgeID(x, y, z, 2);
								m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
							}
							if(edgeTable[tableIndex] & 2048){
								RxVertexID pt = calIntersection(x, y, z, 11);
								uint id = getEdgeID(x, y, z, 11);
								m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
							}
						}
						if(y == m_Grid.iNum[1] - 1){
							if(edgeTable[tableIndex] & 2){
								RxVertexID pt = calIntersection(x, y, z, 1);
								uint id = getEdgeID(x, y, z, 1);
								m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
							}
							if(edgeTable[tableIndex] & 512){
								RxVertexID pt = calIntersection(x, y, z, 9);
								uint id = getEdgeID(x, y, z, 9);
								m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
							}
						}
						if(z == m_Grid.iNum[2] - 1){
							if(edgeTable[tableIndex] & 16){
								RxVertexID pt = calIntersection(x, y, z, 4);
								uint id = getEdgeID(x, y, z, 4);
								m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
							}
							if(edgeTable[tableIndex] & 128){
								RxVertexID pt = calIntersection(x, y, z, 7);
								uint id = getEdgeID(x, y, z, 7);
								m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
							}
						}
						if((x==m_Grid.iNum[0] - 1) && (y==m_Grid.iNum[1] - 1))
							if(edgeTable[tableIndex] & 1024){
								RxVertexID pt = calIntersection(x, y, z, 10);
								uint id = getEdgeID(x, y, z, 10);
								m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
							}
						if((x==m_Grid.iNum[0] - 1) && (z==m_Grid.iNum[2] - 1))
							if(edgeTable[tableIndex] & 64){
								RxVertexID pt = calIntersection(x, y, z, 6);
								uint id = getEdgeID(x, y, z, 6);
								m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
							}
						if((y==m_Grid.iNum[1] - 1) && (z==m_Grid.iNum[2] - 1))
							if(edgeTable[tableIndex] & 32){
								RxVertexID pt = calIntersection(x, y, z, 5);
								uint id = getEdgeID(x, y, z, 5);
								m_i2pt3idVertices.insert(ID2VertexID::value_type(id, pt));
							}

						// ポリゴン生成
						for(uint i = 0; triTable[tableIndex][i] != 255; i += 3){
							RxTriangle triangle;
							uint pointID0, pointID1, pointID2;
							pointID0 = getEdgeID(x, y, z, triTable[tableIndex][i]);
							pointID1 = getEdgeID(x, y, z, triTable[tableIndex][i+1]);
							pointID2 = getEdgeID(x, y, z, triTable[tableIndex][i+2]);
							triangle.pointID[0] = pointID0;
							triangle.pointID[1] = pointID1;
							triangle.pointID[2] = pointID2;
							m_trivecTriangles.push_back(triangle);
						}
					}
				}
			}
		}


		renameVrtsAndTris(vrts, m_nVertices, tris, m_nTriangles);
		calNrms(vrts, m_nVertices, tris, m_nTriangles, nrms, m_nNormals);

		m_bValidSurface = true;
	}


	/*!
	* エッジIDの取得
	* @param[in] nX,nY,nZ グリッド位置
	* @param[in] nEdgeNo エッジ番号
	* @return エッジID
	*/
	uint getEdgeID(uint nX, uint nY, uint nZ, uint nEdgeNo)
	{
		switch(nEdgeNo){
		case 0:
			return getVertexID(nX, nY, nZ) + 1;
		case 1:
			return getVertexID(nX, nY + 1, nZ);
		case 2:
			return getVertexID(nX + 1, nY, nZ) + 1;
		case 3:
			return getVertexID(nX, nY, nZ);
		case 4:
			return getVertexID(nX, nY, nZ + 1) + 1;
		case 5:
			return getVertexID(nX, nY + 1, nZ + 1);
		case 6:
			return getVertexID(nX + 1, nY, nZ + 1) + 1;
		case 7:
			return getVertexID(nX, nY, nZ + 1);
		case 8:
			return getVertexID(nX, nY, nZ) + 2;
		case 9:
			return getVertexID(nX, nY + 1, nZ) + 2;
		case 10:
			return getVertexID(nX + 1, nY + 1, nZ) + 2;
		case 11:
			return getVertexID(nX + 1, nY, nZ) + 2;
		default:
			// Invalid edge no.
			return -1;
		}
	}

	/*!
	* 頂点IDの取得
	* @param[in] nX,nY,nZ グリッド位置
	* @return 頂点ID
	*/
	uint getVertexID(uint nX, uint nY, uint nZ)
	{
		return 3*(nZ*(m_Grid.iNum[1] + 1)*(m_Grid.iNum[0] + 1) + nY*(m_Grid.iNum[0] + 1) + nX);
	}


	/*!
	* 補間によりエッジ上の等値点を計算(サンプルボリュームより)
	* @param[in] nX,nY,nZ グリッド位置
	* @param[in] nEdgeNo エッジ番号
	* @return メッシュ頂点情報
	*/
	RxVertexID calIntersection(uint nX, uint nY, uint nZ, uint nEdgeNo)
	{
		double x1, y1, z1, x2, y2, z2;
		uint v1x = nX, v1y = nY, v1z = nZ;
		uint v2x = nX, v2y = nY, v2z = nZ;

		switch(nEdgeNo){
		case 0:
			v2y += 1;
			break;
		case 1:
			v1y += 1;
			v2x += 1;
			v2y += 1;
			break;
		case 2:
			v1x += 1;
			v1y += 1;
			v2x += 1;
			break;
		case 3:
			v1x += 1;
			break;
		case 4:
			v1z += 1;
			v2y += 1;
			v2z += 1;
			break;
		case 5:
			v1y += 1;
			v1z += 1;
			v2x += 1;
			v2y += 1;
			v2z += 1;
			break;
		case 6:
			v1x += 1;
			v1y += 1;
			v1z += 1;
			v2x += 1;
			v2z += 1;
			break;
		case 7:
			v1x += 1;
			v1z += 1;
			v2z += 1;
			break;
		case 8:
			v2z += 1;
			break;
		case 9:
			v1y += 1;
			v2y += 1;
			v2z += 1;
			break;
		case 10:
			v1x += 1;
			v1y += 1;
			v2x += 1;
			v2y += 1;
			v2z += 1;
			break;
		case 11:
			v1x += 1;
			v2x += 1;
			v2z += 1;
			break;
		}

		x1 = m_Grid.fMin[0]+v1x*m_Grid.fWidth[0];
		y1 = m_Grid.fMin[1]+v1y*m_Grid.fWidth[1];
		z1 = m_Grid.fMin[2]+v1z*m_Grid.fWidth[2];
		x2 = m_Grid.fMin[0]+v2x*m_Grid.fWidth[0];
		y2 = m_Grid.fMin[1]+v2y*m_Grid.fWidth[1];
		z2 = m_Grid.fMin[2]+v2z*m_Grid.fWidth[2];

		uint slice0 = (m_Grid.iNum[0] + 1);
		uint slice1 = slice0*(m_Grid.iNum[1] + 1);
		double val1 = m_ptScalarField[v1z*slice1 + v1y*slice0 + v1x];
		double val2 = m_ptScalarField[v2z*slice1 + v2y*slice0 + v2x];
		RxVertexID intersection = interpolate(x1, y1, z1, x2, y2, z2, val1, val2);

		return intersection;
	}


	/*!
	* グリッドエッジ両端の陰関数値から線型補間で等値点を計算
	* @param[in] fX1,fY1,fZ1 端点座標1
	* @param[in] fX2,fY2,fZ2 端点座標2
	* @param[in] tVal1 端点座標1でのスカラー値
	* @param[in] tVal2 端点座標2でのスカラー値
	* @return 頂点情報
	*/
	RxVertexID interpolate(double fX1, double fY1, double fZ1, double fX2, double fY2, double fZ2, double tVal1, double tVal2)
	{
		RxVertexID interpolation;
		double mu;

		mu = double((m_tIsoLevel - tVal1))/(tVal2 - tVal1);
		interpolation.x = fX1 + mu*(fX2 - fX1);
		interpolation.y = fY1 + mu*(fY2 - fY1);
		interpolation.z = fZ1 + mu*(fZ2 - fZ1);

		return interpolation;
	}


	/*!
	* メッシュ頂点，幾何情報を出力形式で格納
	* @param[out] vrts 頂点座標
	* @param[out] nvrts 頂点数
	* @param[out] tris 三角形ポリゴン幾何情報
	* @param[out] ntris 三角形ポリゴン数
	*/
	void renameVrtsAndTris(vector<glm::vec3> &vrts, uint &nvrts, vector<int> &tris, uint &ntris)
	{
		uint nextID = 0;
		ID2VertexID::iterator mapIterator = m_i2pt3idVertices.begin();
		RxTriangleVector::iterator vecIterator = m_trivecTriangles.begin();

		// Rename vertices.
		while(mapIterator != m_i2pt3idVertices.end()){
			(*mapIterator).second.newID = nextID;
			nextID++;
			mapIterator++;
		}

		// Now rename triangles.
		while(vecIterator != m_trivecTriangles.end()){
			for(uint i = 0; i < 3; i++){
				uint newID = m_i2pt3idVertices[(*vecIterator).pointID[i]].newID;
				(*vecIterator).pointID[i] = newID;
			}
			vecIterator++;
		}

		// Copy all the vertices and triangles into two arrays so that they
		// can be efficiently accessed.
		// Copy vertices.
		mapIterator = m_i2pt3idVertices.begin();
		nvrts = (int)m_i2pt3idVertices.size();
		vrts.resize(nvrts);
		for(uint i = 0; i < nvrts; i++, mapIterator++){
			vrts[i][0] = (*mapIterator).second.x;
			vrts[i][1] = (*mapIterator).second.y;
			vrts[i][2] = (*mapIterator).second.z;
		}
		// Copy vertex indices which make triangles.
		vecIterator = m_trivecTriangles.begin();
		ntris = (int)m_trivecTriangles.size();
		tris.resize(ntris*3);
		for(uint i = 0; i < ntris; i++, vecIterator++){
			tris[3*i+0] = (*vecIterator).pointID[0];
			tris[3*i+1] = (*vecIterator).pointID[1];
			tris[3*i+2] = (*vecIterator).pointID[2];
		}

		m_i2pt3idVertices.clear();
		m_trivecTriangles.clear();
	}

	/*!
	* 頂点法線計算
	* @param[in] vrts 頂点座標
	* @param[in] nvrts 頂点数
	* @param[in] tris 三角形ポリゴン幾何情報
	* @param[in] ntris 三角形ポリゴン数
	* @param[out] nrms 法線
	* @param[out] nnrms 法線数(=頂点数)
	*/
	void calNrms(const vector<glm::vec3> &vrts, uint nvrts, const vector<int> &tris, uint ntris, vector<glm::vec3> &nrms, uint &nnrms)
	{
		nnrms = nvrts;
		nrms.resize(nnrms);

		// Set all normals to 0.
		for(uint i = 0; i < nnrms; i++){
			nrms[i][0] = 0;
			nrms[i][1] = 0;
			nrms[i][2] = 0;
		}

		// Calculate normals.
		for(uint i = 0; i < ntris; i++){
			glm::vec3 vec1, vec2, normal;
			uint id0, id1, id2;
			id0 = tris[3*i+0];
			id1 = tris[3*i+1];
			id2 = tris[3*i+2];

			vec1 = vrts[id1]-vrts[id0];
			vec2 = vrts[id2]-vrts[id0];
			normal = glm::cross(vec1, vec2);

			nrms[id0] += normal;
			nrms[id1] += normal;
			nrms[id2] += normal;
		}

		// Normalize normals.
		for(uint i = 0; i < nnrms; i++){
			nrms[i] = glm::normalize(nrms[i]);
		}
	}

};




#endif // _RX_MC_MESH_H_

