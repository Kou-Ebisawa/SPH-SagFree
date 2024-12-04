/*! 
  @file rx_mesh.h

  @brief メッシュ構造の定義
 
  @author Makoto Fujisawa
  @date 2010
*/
// FILE --rx_mesh.h--

#ifndef _RX_MESH_H_
#define _RX_MESH_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <sstream>

#include <string>
#include <vector>
#include <map>
#include <set>


#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/ext.hpp"	// for glm::to_string()

using namespace std;

//-----------------------------------------------------------------------------
// 材質クラス
//-----------------------------------------------------------------------------
class rxMaterialOBJ
{
public:
	string name;		//!< 材質の名前

	glm::vec4 diffuse;		//!< 拡散色(GL_DIFFUSE)
	glm::vec4 specular;		//!< 鏡面反射色(GL_SPECULAR)
	glm::vec4 ambient;		//!< 環境光色(GL_AMBIENT)

	glm::vec4 color;			//!< 色(glColor)
	glm::vec4 emission;		//!< 放射色(GL_EMISSION)

	double shininess;	//!< 鏡面反射指数(GL_SHININESS)

	int illum;
	string tex_file;	//!< テクスチャファイル名
	unsigned int tex_name;	//!< テクスチャオブジェクト
};


typedef map<string, rxMaterialOBJ> rxMTL;


//-----------------------------------------------------------------------------
// メッシュクラス
//-----------------------------------------------------------------------------
// ポリゴン
class rxFace
{
public:
	vector<int> vert_idx;	//!< 頂点インデックス
	string material_name;	//!< 材質名
	vector<glm::vec2> texcoords;	//!< テクスチャ座標
	int attribute;			//!< 属性

public:
	rxFace() : attribute(0) {}

public:
	// オペレータによる頂点アクセス
	inline int& operator[](int i){ return vert_idx[i]; }
	inline int  operator[](int i) const { return vert_idx[i]; }

	// 関数による頂点アクセス
	inline int& at(int i){ return vert_idx.at(i); }
	inline int  at(int i) const { return vert_idx.at(i); }

	//! ポリゴン頂点数の変更
	void resize(int size)
	{
		vert_idx.resize(size);
	}

	//! ポリゴン頂点数の参照
	int size(void) const
	{
		return (int)vert_idx.size();
	}

	//! 初期化
	void clear(void)
	{
		vert_idx.clear();
		material_name = "";
		texcoords.clear();
	}
};

// 三角形ポリゴン
class rxTriangle : public rxFace
{
public:
	rxTriangle()
	{
		vert_idx.resize(3);
	}
};

// ポリゴンオブジェクト
class rxPolygons
{
public:
	vector<glm::vec3> vertices;	//!< 頂点座標
	vector<glm::vec3> normals;	//!< 頂点法線
	vector<rxFace> faces;	//!< ポリゴン
	rxMTL materials;		//!< 材質
	int open;				//!< ファイルオープンフラグ

	// Skinning関連
	vector< map<int, double> > weights;	//!< 各頂点のSkinning weight (対応するボーン番号(int型)と重み(double型)のマップ)


public:
	//! コンストラクタ
	rxPolygons() : open(0) {}
	//! デストラクタ
	~rxPolygons(){}

	//! 描画
	void Draw(int draw = 0x04, float dn = 0.02f, bool col = true);

protected:
	//! double版の材質設定
	void glMaterialdv(GLenum face, GLenum pname, const GLdouble *params);
};


//! エッジ
struct rxEdge
{
	int v[2];		//!< エッジ両端の頂点インデックス
	set<int> f;		//!< エッジにつながるポリゴン情報
	double len;		//!< エッジ長さ
	int attribute;	//!< 属性
};

//! 拡張ポリゴンオブジェクト(エッジ情報追加)
class rxPolygonsE : public rxPolygons
{
public:
	vector<rxEdge> edges;		//!< エッジ
	vector< set<int> > vedges;	//!< 頂点につながるエッジ情報
	vector< set<int> > vfaces;	//!< 頂点につながるポリゴン情報
	vector< set<int> > fedges;	//!< ポリゴンに含まれるエッジ情報

	vector<glm::vec3> vcol;			//!< 頂点色

public:
	//! コンストラクタ
	rxPolygonsE() : rxPolygons() {}
	//! デストラクタ
	virtual ~rxPolygonsE(){}

	//! 初期化
	void Clear(void)
	{
		vertices.clear();
		normals.clear();
		faces.clear();
		materials.clear();
		edges.clear();
		vedges.clear();
		vfaces.clear();
		fedges.clear();
	}

};



//海老沢追加
//境界の頂点とエッジの情報のみを持ったポリゴンを複数扱う
class Edge_Vert {
public:
	vector<glm::vec3> vertices;//頂点座標
	vector<rxEdge> edges;//エッジの情報

	void Edge_Vert::vertAdd(glm::vec3 vert) {
		//エッジが存在しない場合
		if (edges.size() == 0) {
			vertices.push_back(vert);
		}
		else {
			vertices.push_back(vert);
			//エッジの連結を変える
			edges[edges.size() - 1].v[1] = vertices.size() - 1;
		}
	}

	//頂点が既に一つ以上追加されていることを前提とする
	void Edge_Vert::edgeAdd(rxEdge edge) {
		//エッジのインデックスを書き換える
		edge.v[0] = vertices.size() - 1;
		//追加されなかった場合のために，最初のインデックスをさしておく
		edge.v[1] = 0;

		edges.push_back(edge);
	}

	//配列をクリア
	void Edge_Vert::Clear() {
		vertices.clear();
		edges.clear();
	}
};


//-----------------------------------------------------------------------------
// MARK:ポリゴンの描画
//-----------------------------------------------------------------------------
inline void rxPolygons::glMaterialdv(GLenum face, GLenum pname, const GLdouble *params)
{
	GLfloat col[4];
	col[0] = (GLfloat)params[0];
	col[1] = (GLfloat)params[1];
	col[2] = (GLfloat)params[2];
	col[3] = (GLfloat)params[3];
	glMaterialfv(face, pname, col);
}

/*!
 * ポリゴンの描画
 * @param[in] polys ポリゴンデータ
 * @param[in] draw 描画フラグ(下位ビットから頂点,エッジ,面,法線 - 1,2,4,8)
 */
inline void rxPolygons::Draw(int draw, float dn, bool col)
{
	// 頂点数とポリゴン数
	int vn = (int)vertices.size();
	int pn = (int)faces.size();
		
	if(draw & 0x02){
		// エッジ描画における"stitching"をなくすためのオフセットの設定
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
	}

	if(draw & 0x04){
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glEnable(GL_LIGHTING);
		if(col) glColor3d(0.0, 0.0, 1.0);
		bool use_mat = true;
		if(materials.empty()){
			// 材質が設定されていない，もしくは，エッジのみ描画の場合，GL_COLOR_MATERIALを用いる
			use_mat = false;

			glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
			glEnable(GL_COLOR_MATERIAL);
		}
		else{
			glDisable(GL_COLOR_MATERIAL);
		}


		bool tc = false;
		rxMaterialOBJ *cur_mat = NULL, *pre_mat = NULL;

		// すべてのポリゴンを描画
		for(int i = 0; i < pn; ++i){
			rxFace *face = &(faces[i]);
			int n = (int)face->vert_idx.size();

			if(use_mat){
				// 材質の設定
				cur_mat = &materials[face->material_name];

				// 現在の材質が前のと異なる場合のみ，OpenGLの材質情報を更新
				if(cur_mat != pre_mat){
					glColor4fv(glm::value_ptr(cur_mat->color));
					glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, glm::value_ptr(cur_mat->diffuse));
					glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, glm::value_ptr(cur_mat->specular));
					glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,  glm::value_ptr(cur_mat->ambient));
					glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, glm::value_ptr(cur_mat->emission));
					glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, (GLfloat)cur_mat->shininess);

					// テクスチャがあればバインド
					if(cur_mat->tex_name){
						glEnable(GL_TEXTURE_2D);
						glBindTexture(GL_TEXTURE_2D, cur_mat->tex_name);
						glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
					}
					else{
						glBindTexture(GL_TEXTURE_2D, 0);
						glDisable(GL_TEXTURE);
					}

					pre_mat = cur_mat;
				}
			}

			glm::vec3 cols[3] = { glm::vec3(1.0, 0.0, 0.0), glm::vec3(0.0, 0.0, 1.0), glm::vec3(0.0, 1.0, 0.0) };

			// ポリゴン描画
			glBegin(GL_POLYGON);
			for(int j = 0; j < n; ++j){
				int idx = face->vert_idx[j];
				if(idx >= 0 && idx < vn){
					glm::vec3 col(0.0);
					if(!weights.empty() && !weights[idx].empty()){
						int ncol = 0;
						map<int, double>::iterator iter;
						for(iter = weights[idx].begin(); iter != weights[idx].end(); ++iter){
							int joint_idx = iter->first;
							col += float(iter->second)*cols[joint_idx%3];
							ncol++;
						}
						//if(ncol != 0) col /= (double)ncol;
					}
					glColor3fv(glm::value_ptr(col));

					glNormal3fv(glm::value_ptr(normals[idx]));
					if(!face->texcoords.empty()) glTexCoord2fv(glm::value_ptr(face->texcoords[j]));
					glVertex3fv(glm::value_ptr(vertices[idx]));
				}
			}
			glEnd();
		}

		glBindTexture(GL_TEXTURE_2D, 0);
		glDisable(GL_TEXTURE_2D);
	}


	// 頂点描画
	if(draw & 0x01){
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);

		glDisable(GL_LIGHTING);
		if(col) glColor3d(1.0, 0.0, 0.0);
		glPointSize(5.0);
		glBegin(GL_POINTS);
		for(int i = 0; i < vn; ++i){
			glVertex3fv(glm::value_ptr(vertices[i]));
		}
		glEnd();
	}

	// エッジ描画
	if(draw & 0x02){
		glDisable(GL_LIGHTING);
		if(col) glColor3d(0.3, 1.0, 0.0);
		glLineWidth(1.0);
		for(int i = 0; i < pn; ++i){
			rxFace *face = &faces[i];
			int n = (int)face->vert_idx.size();
			
			glBegin(GL_LINE_LOOP);
			for(int j = 0; j < n; ++j){
				int idx = face->vert_idx[j];
				if(idx >= 0 && idx < vn){
					glVertex3fv(glm::value_ptr(vertices[idx]));
				}
			}
			glEnd();
		}
	}
	//glDisable(GL_POLYGON_OFFSET_FILL);

	// 法線描画
	if(draw & 0x08){
		glDisable(GL_LIGHTING);
		glColor3d(0.0, 0.9, 0.0);
		glLineWidth(1.0);
		for(int i = 0; i < vn; ++i){
			glBegin(GL_LINES);
			glVertex3fv(glm::value_ptr(vertices[i]));
			glVertex3fv(glm::value_ptr(vertices[i]+dn*normals[i]));
			glEnd();
		}
	}
}




//-----------------------------------------------------------------------------
// メッシュデータの読み込み，保存に便利な関数
//-----------------------------------------------------------------------------
/*!
 * 座標の回転
 * @param[in] pos0 回転したい座標値
 * @param[in] rot0 x,y,z軸周りの回転角(deg)
 * @return 開店後の座標値
 */
inline glm::vec3 Rotate3(const glm::vec3 &pos0, const glm::vec3 &rot)
{
	glm::vec3 S, C;
	S = glm::vec3(sin(rot[0]), sin(rot[1]), sin(rot[2]));
	C = glm::vec3(cos(rot[0]), cos(rot[1]), cos(rot[2]));

	glm::vec3 pos1;
	pos1[0] = C[1]*C[2]*pos0[0]-C[1]*S[2]*pos0[1]+S[1]*pos0[2];
	pos1[1] = (S[0]*S[1]*C[2]+C[0]*S[2])*pos0[0]+(-S[0]*S[1]*S[2]+C[0]*C[2])*pos0[1]-S[0]*C[1]*pos0[2];
	pos1[2] = (-C[0]*S[1]*C[2]+S[0]*S[2])*pos0[0]+(C[0]*S[1]*S[2]+S[0]*C[2])*pos0[1]+C[0]*C[1]*pos0[2];
	return pos1;
}

/*!
 * "文字列 数値"から数値部分の文字列のみを取り出す
 * @param[in] buf 元の文字列
 * @param[in] head 文字部分
 * @param[out] sub 数値部分の文字列
 */
inline bool StringToString(const string &buf, const string &head, string &sub)
{
	size_t pos = 0;
	if((pos = buf.find(head)) == string::npos) return false;
	pos += head.size();

	if((pos = buf.find_first_not_of(" 　\t", pos)) == string::npos) return false;
	
	sub = buf.substr(pos);
	return true;
}

/*!
 * "文字列 数値"から数値部分のみを取り出す
 * @param[in] buf 元の文字列
 * @param[in] head 文字部分
 * @param[out] val 数値
 */
inline bool StringToDouble(const string &buf, const string &head, double &val)
{
	string sub;
	if(StringToString(buf, head, sub)){
		sscanf(&sub[0], "%lf", &val);
		return true;
	}
	else{
		return false;
	}
}

/*!
 * "文字列 数値"から数値部分のみを取り出す
 * @param[in] buf 元の文字列
 * @param[in] head 文字部分
 * @param[out] val 数値
 */
inline bool StringToInt(const string &buf, const string &head, int &val)
{
	string sub;
	if(StringToString(buf, head, sub)){
		sscanf(&sub[0], "%d", &val);
		return true;
	}
	else{
		return false;
	}
}


/*!
 * "x y z"の形式の文字列からglm::vec2型へ変換
 * @param[in] s 文字列
 * @param[out] v 値
 * @return 要素記述数
 */
inline int StringToVec2s(const string &buf, const string &head, glm::vec2 &v)
{
	string sub;
	StringToString(buf, head, sub);

	glm::vec2 tmp;

	if(sscanf(&sub[0], "%f %f", &tmp[0], &tmp[1]) != 2){
		return 0;
	}
	v = tmp;

	return 2;
}

/*!
 * "x y z"の形式の文字列からglm::vec3型へ変換
 * @param[in] s 文字列
 * @param[out] v 値
 * @return 要素記述数
 */
inline int StringToVec3s(const string &buf, const string &head, glm::vec3 &v)
{
	string sub;
	StringToString(buf, head, sub);

	glm::vec3 tmp;
	float tmp_w;

	if(sscanf(&sub[0], "%f %f %f %f", &tmp[0], &tmp[1], &tmp[2], &tmp_w) != 4){
		if(sscanf(&sub[0], "%f %f %f", &tmp[0], &tmp[1], &tmp[2]) != 3){
			return 0;
		}
	}
	v = tmp;

	return 3;
}

inline string Vec3ToString(const glm::vec3 v)
{
	stringstream ss;
	ss << v[0] << " " << v[1] << " " << v[2];
	return ss.str();
	//string s;
	//sscanf(&s[0], "%f %f %f", v[0], v[1], v[2]);
	//return s;
}


/*!
 * 先頭の空白(スペース，タブ)を削除
 * @param[in] buf 元の文字列
 * @return 空白削除後の文字列
 */
inline string GetDeleteSpace(const string &buf)
{
	string buf1 = buf;

	size_t pos;
	while((pos = buf1.find_first_of(" 　\t")) == 0){
		buf1.erase(buf1.begin());
		if(buf1.empty()) break;
	}

	return buf1;
}

/*!
 * 先頭の空白(スペース，タブ)を削除
 * @param[inout] buf 処理文字列
 */
inline void DeleteHeadSpace(string &buf)
{
	size_t pos;
	while((pos = buf.find_first_of(" 　\t")) == 0){
		buf.erase(buf.begin());
		if(buf.empty()) break;
	}
}

/*!
 * 空白(スペース，タブ)を削除
 * @param[inout] buf 処理文字列
 */
inline void DeleteSpace(string &buf)
{
	size_t pos;
	while((pos = buf.find_first_of(" 　\t")) != string::npos){
		buf.erase(pos, 1);
	}
}

/*!
 * 文字列に含まれる指定された文字列の数を返す
 * @param[in] s 元の文字列
 * @param[in] c 検索文字列
 * @return 含まれる数
 */
inline int CountString(string &s, int offset, string c)
{
	int count = 0;
	size_t pos0 = offset, pos = 0;
	int n = (int)c.size();

	while((pos = s.find(c, pos0)) != string::npos){
		if(pos != pos0){
			count++;
		}
		else{
			s.erase(s.begin()+pos);
		}
		pos0 = pos+n;
	}

	// 最後の文字列除去
	if(s.rfind(c) == s.size()-n){
		count--;
	}

	return count;
}

/*!
 * 文字列からベクトル値を抽出
 *  - たとえば，"0, 1, 2"だと，"がend_str, ","がbrk_str
 * @param[in] data 元の文字列
 * @param[in] end_str ベクトルの区切り文字
 * @param[in] brk_str ベクトル要素間の区切り文字
 * @param[in] min_elem 最小要素数
 * @param[out] vecs ベクトル値
 * @return 見つかったベクトルの数
 */
template<class T> 
inline int ExtractVector(string data, const string &end_str, const string &brk_str, int min_elems, 
						 vector< vector<T> > &vecs)
{
	data += end_str;	// 後の処理のために区切り文字列を最後に足しておく
	int n = 0;
	size_t cpos[2] = {0, 0};
	while((cpos[1] = data.find(end_str, cpos[0])) != string::npos){
		// 区切り文字が見つかったら，前回の区切り文字位置との間の文字列を抜き出す
		string sub = data.substr(cpos[0], cpos[1]-cpos[0]);
		if(sub.empty()){
			cpos[0] = cpos[1]+end_str.size();
			break;
		}

		// 抜き出した文字列を各ベクトル要素に分割
		sub += brk_str;
		vector<T> val;
		size_t spos[2] = {0, 0};
		while((spos[1] = sub.find(brk_str, spos[0])) != string::npos){
			string val_str = sub.substr(spos[0], spos[1]-spos[0]);
			DeleteSpace(val_str);
			if(val_str.empty()){
				spos[0] = spos[1]+brk_str.size();
				continue;
			}

			val.push_back((T)atof(val_str.c_str()));
			spos[0] = spos[1]+brk_str.size();
		}
		if((int)val.size() >= min_elems){
			vecs.push_back(val);
			n++;
		}
		cpos[0] = cpos[1]+end_str.size();
	}

	return n;
}

/*!
 * 文字列から数値列を抽出
 * @param[in] data 元の文字列
 * @param[in] brk_str 要素間の区切り文字
 * @param[out] strs 各数値(string)
 * @return 見つかった要素の数
 */
inline int ExtractSubStr(string data, const string &brk_str, vector<string> &strs)
{
	data += brk_str;	// 後の処理のために区切り文字列を最後に足しておく
	int n = 0;
	size_t cpos[2] = {0, 0};
	while((cpos[1] = data.find(brk_str, cpos[0])) != string::npos){
		// 区切り文字が見つかったら，前回の区切り文字位置(cpos[0])との間の文字列を抜き出す
		string sub = data.substr(cpos[0], cpos[1]-cpos[0]);
		if(!sub.empty()){
			strs.push_back(sub);
			n++;
		}
		
		cpos[0] = cpos[1]+brk_str.size();
	}

	return n;
}


/*!
 * ポリゴンを三角形に分解
 * @param[out] ply 三角形に分解されたポリゴン列
 * @param[in] vidx ポリゴン(頂点数4以上)の頂点インデックス
 * @param[in] mat_name 材質名
 * @return 分解された三角形数
 */
inline int PolyToTri(vector<rxFace> &plys, 
					 const vector<int> &vidxs, 
					 string mat_name)
{
	int n = (int)vidxs.size();

	if(n <= 3) return 0;

	rxFace face;
	face.resize(3);
	face.material_name = mat_name;

	int num_tri = n-2;
	for(int i = 0; i < num_tri; ++i){
		face[0] = vidxs[0];
		face[1] = vidxs[i+1];
		face[2] = vidxs[i+2];

		if(face[1] != face[2] && face[0] != face[1] && face[2] != face[0]){
			plys.push_back(face);
		}
	}

	return num_tri;
}

/*!
 * ポリゴンを三角形に分解(with テクスチャ座標)
 * @param[out] ply 三角形に分解されたポリゴン列
 * @param[in] vidx ポリゴン(頂点数4以上)の頂点インデックス
 * @param[in] tidx ポリゴン(頂点数4以上)のテクスチャ座標インデックス
 * @param[in] vtc  テクスチャ座標ベクトル
 * @param[in] mat_name 材質名
 * @return 分解された三角形数
 */
inline int PolyToTri(vector<rxFace> &plys, 
					 const vector<int> &vidxs, 
					 const vector<int> &tidxs, 
					 const vector<glm::vec2> &vtc, 
					 string mat_name)
{
	int n = (int)vidxs.size();

	if(n <= 3) return 0;

	rxFace face;
	face.vert_idx.resize(3);
	face.material_name = mat_name;
	face.texcoords.resize(3);

	bool tc = !vtc.empty();

	int num_tri = n-2;
	for(int i = 0; i < num_tri; ++i){
		face[0] = vidxs[0];
		face[1] = vidxs[i+1];
		face[2] = vidxs[i+2];

		if(tc){
			face.texcoords[0] = vtc[tidxs[0]];
			face.texcoords[1] = vtc[tidxs[i+1]];
			face.texcoords[2] = vtc[tidxs[i+2]];
		}
		else{
			face.texcoords[0] = glm::vec2(0.0);
			face.texcoords[1] = glm::vec2(0.0);
			face.texcoords[2] = glm::vec2(0.0);
		}

		plys.push_back(face);
	}

	return num_tri;
}


/*!
 * ファイル名からフォルダパスのみを取り出す
 * @param[in] fn ファイル名(フルパス or 相対パス)
 * @return フォルダパス
 */
inline string ExtractDirPath(const string &fn)
{
	string::size_type pos;
	if((pos = fn.find_last_of("/")) == string::npos){
		if((pos = fn.find_last_of("\\")) == string::npos){
			return "";
		}
	}

	return fn.substr(0, pos);
}

/*!
 * ファイル名から拡張子を削除
 * @param[in] fn ファイル名(フルパス or 相対パス)
 * @return フォルダパス
 */
inline string ExtractPathWithoutExt(const string &fn)
{
	string::size_type pos;
	if((pos = fn.find_last_of(".")) == string::npos){
		return fn;
	}

	return fn.substr(0, pos);
}


/*!
 * 文字列小文字化
 * @param[inout] str 文字列
 */
inline void StringToLower(string &str)
{
	string::size_type i, size;

	size = str.size();

	for(i = 0; i < size; i++){
		if(str[i] >= 'A' && str[i] <= 'Z') str[i] += 32;
	}

	return;
}

/*!
 * 文字列が整数値を表しているかを調べる
 * @param[inout] str 文字列
 * @return 整数値ならtrue
 */
inline bool IsInteger(const string &str)
{
	if(str.find_first_not_of("-0123456789 \t") != string::npos) {
		return false;
	}

	return true;
}

/*!
 * 文字列が実数値を表しているかを調べる
 * @param[inout] str 文字列
 * @return 実数値ならtrue
 */
inline bool IsNumeric(const string &str)
{
	if(str.find_first_not_of("-0123456789. Ee\t") != string::npos) {
		return false;
	}

	return true;
}


/*!
 * バイナリファイルから型指定でデータを読み込む
 * @param[inout] file ファイル入力ストリーム
 * @return 読み込んだ数値
 */
template<class T> 
inline T ReadBinary(ifstream &file)
{
	T data;
	file.read((char*)&data, sizeof(T));
	return data;
}

/*!
 * バイナリファイルから型・サイズ指定でデータを読み込む
 * @param[inout] file ファイル入力ストリーム
 * @param[inout] byte 読み込むバイト数
 * @return 読み込んだ数値
 */
template<class T> 
inline T ReadBinary(ifstream &file, int byte)
{
	T data = 0;
	file.read((char*)&data, byte);
	return data;
}

/*!
 * バイナリファイルから2バイト分読み込んでint型に格納
 * @param[inout] file ファイル入力ストリーム
 * @return 読み込んだ数値
 */
inline int ReadInt(ifstream &file)
{
	int data = 0;
	file.read((char*)&data, 2);
	return data;
}

/*!
 * バイナリファイルから文字列を読み取る
 * @param[inout] file ファイル入力ストリーム
 * @param[out] name 文字列(何もなければ空)
 * @param[in] max_size 読み込む文字数(-1なら\0まで読み込む)
 * @return 文字列がない(ストリームの最初が\0)ならfalseを返す
 */
inline bool ReadString(ifstream &file, string &name, int max_size)
{
	char c = ReadBinary<char>(file);
	if((int)c == 0) return false;

	name.clear();
	name.push_back(c);
	do{
		c = ReadBinary<char>(file);
		name.push_back(c);
	}while(((int)c != 0) && (max_size == -1 || (int)name.size() < max_size));
	name.push_back((char)0);

	return true;
}


//-----------------------------------------------------------------------------
// Geometry processing
//-----------------------------------------------------------------------------
/*!
 * 法線反転
 * @param[in] vs 頂点リスト
 * @param[inout] ps メッシュリスト
 * @param[inout] ns 法線リスト
 */
static void ReverseNormals(const vector<glm::vec3> &vs, vector<int> &ps, vector<glm::vec3> &ns)
{
	if(ps.empty()) return;

	int nvp = 3;
	int pn = (int)ps.size()/nvp;

	vector<int> idx(nvp);
	for(int i = 0; i < pn; ++i){
		for(int j = 0; j < nvp; ++j){
			idx[j] = ps[nvp*i+j];
		}

		for(int j = 0; j < nvp; ++j){
			ps[nvp*i+j] = idx[nvp-j-1];
		}
	}

	if(ns.empty() || (int)ns.size() != pn){
		// 法線再計算
		ns.resize(pn);
		for(int i = 0; i < pn; ++i){
			ns[i] = glm::normalize(glm::cross(vs[ps[nvp*i+nvp-1]]-vs[ps[nvp*i+0]], vs[ps[nvp*i+1]]-vs[ps[nvp*i+0]]));
		}
	}
	else{
		// 法線反転
		for(int i = 0; i < pn; ++i){
			ns[i] = -ns[i];
		}
	}
}

/*!
 * メッシュ法線計算
 * @param[in] vs 頂点リスト
 * @param[inout] ps メッシュリスト
 * @param[inout] ns 法線リスト
 */
static void CalNormals(const vector<glm::vec3> &vs, vector<int> &ps, vector<glm::vec3> &ns)
{
	if(ps.empty()) return;

	int nvp = 3;
	int pn = (int)ps.size()/nvp;
	if(ns.empty() || (int)ns.size() != pn){
		// 法線再計算
		ns.resize(pn);
	}

	for(int i = 0; i < pn; ++i){
		int *idx = &ps[nvp*i];
		ns[i] = glm::normalize(glm::cross(vs[idx[1]]-vs[idx[0]], vs[idx[nvp-1]]-vs[idx[0]]));
		//ns[i] *= -1;
	}
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
static void CalVertexNormals(const vector<glm::vec3> &vrts, int nvrts, vector<rxTriangle> &tris, int ntris, 
							 vector<glm::vec3> &vnrms)
{
	int nnrms = nvrts;
	vnrms.resize(nnrms);
	
	// 頂点法線の初期化
	for(int i = 0; i < nnrms; i++){
		vnrms[i][0] = 0;
		vnrms[i][1] = 0;
		vnrms[i][2] = 0;
	}

	// 法線計算
	for(int i = 0; i < ntris; i++){
		glm::vec3 edge_vec1, edge_vec2, face_normal;

		// 面法線を計算
		edge_vec1 = vrts[tris[i][1]]-vrts[tris[i][0]];
		edge_vec2 = vrts[tris[i][2]]-vrts[tris[i][0]];
		face_normal = glm::normalize(glm::cross(edge_vec1, edge_vec2));

		// ポリゴンに所属する頂点の法線に積算
		vnrms[tris[i][0]] += face_normal;
		vnrms[tris[i][1]] += face_normal;
		vnrms[tris[i][2]] += face_normal;
	}

	// 頂点法線を正規化
	for(int i = 0; i < nnrms; i++){
		vnrms[i] = glm::normalize(vnrms[i]);
	}

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
static void CalVertexNormals(const vector<glm::vec3> &vrts, int nvrts, vector<rxFace> &tris, int ntris, 
							 vector<glm::vec3> &vnrms)
{
	int nnrms = nvrts;
	vnrms.resize(nnrms);
	
	// 頂点法線の初期化
	for(int i = 0; i < nnrms; i++){
		vnrms[i][0] = 0;
		vnrms[i][1] = 0;
		vnrms[i][2] = 0;
	}

	// 法線計算
	for(int i = 0; i < ntris; i++){
		glm::vec3 edge_vec1, edge_vec2, face_normal;
		int n = (int)tris[i].size();

		// 面法線を計算
		edge_vec1 = vrts[tris[i][1]]-vrts[tris[i][0]];
		edge_vec2 = vrts[tris[i][n-1]]-vrts[tris[i][0]];
		face_normal = glm::normalize(glm::cross(edge_vec1, edge_vec2));

		// ポリゴンに所属する頂点の法線に積算
		for(int j = 0; j < n; ++j){
			vnrms[tris[i][j]] += face_normal;
		}
	}

	// 頂点法線を正規化
	for(int i = 0; i < nnrms; i++){
		vnrms[i] = glm::normalize(vnrms[i]);
	}
}
/*!
 * 頂点法線計算
 * @param[in] polys ポリゴン
 */
static void CalVertexNormals(rxPolygons &polys)
{
	//(const vector<glm::vec3> &vrts, uint nvrts, vector<rxTriangle> &tris, uint ntris, vector<glm::vec3> &nrms)
	int pn = (int)polys.faces.size();
	int vn = (int)polys.vertices.size();

	vector<glm::vec3> fnrms;
	fnrms.resize(pn);

	int max_n = 0;

	// 面法線の計算
	for(int i = 0; i < pn; ++i){
		int n = (int)polys.faces[i].vert_idx.size()-1;
		fnrms[i] = glm::normalize(glm::cross(polys.vertices[polys.faces[i][1]]-polys.vertices[polys.faces[i][0]], 
											 polys.vertices[polys.faces[i][n]]-polys.vertices[polys.faces[i][0]]));

		n = n+1;
		if(n > max_n) max_n = n;
	}

	polys.normals.clear();
	polys.normals.resize(vn);
	for(int i = 0; i < vn; ++i){
		polys.normals[i] = glm::vec3(0.0);
	}

	// 頂点法線の計算
	vector<int> id;
	id.resize(max_n);
	for(int i = 0; i < pn; ++i){
		int n = (int)polys.faces[i].vert_idx.size();
		for(int j = 0; j < n; ++j){
			polys.normals[polys.faces[i][j]] += fnrms[i];
		}
	}

	// 法線正規化
	for(int i = 0; i < vn; i++){
		polys.normals[i] = glm::normalize(polys.normals[i]);
	}
}


/*!
 * 幾何情報を持たないポリゴン頂点列とポリゴン法線から頂点法線を計算
 * @param[in] vrts 幾何情報無しのポリゴン頂点列
 * @param[inout] nrms ポリゴン法線がそれぞれの頂点法線として格納されている
 */
static void CalVertexNormalWithoutGeometry(const vector<glm::vec3> &vrts, vector<glm::vec3> &nrms)
{
	int n = (int)vrts.size();

	int same_vrts[256];

	vector<int> find_vrts;
	find_vrts.resize(n);
	for(int i = 0; i < n; ++i){
		find_vrts[i] = 0;
	}

	int nvrts = 0;
	int avg_cnt = 0;

	double eps2 = 1e-4;
	for(int i = 0; i < n; ++i){
		if(find_vrts[i]) continue;
		
		glm::vec3 pos0 = vrts[i];
		glm::vec3 nrm0 = nrms[i];
		same_vrts[0] = i;
		int cnt = 1;
		find_vrts[i] = 1;
		for(int j = 0; j < n; ++j){
			if(find_vrts[j]) continue;

			glm::vec3 pos1 = vrts[j];

			if(glm::length2(pos0-pos1) < eps2){
				find_vrts[j] = 1;
				nrm0 += nrms[j];
				same_vrts[cnt] = j;
				cnt++;
			}
		}

		nrm0 /= (double)cnt;

		for(int j = 0; j < cnt; ++j){
			nrms[same_vrts[j]] = nrm0;
		}

		avg_cnt += cnt;
		nvrts++;
	}

	avg_cnt /= nvrts;
	//cout << avg_cnt << ", " << nvrts << endl;
}


/*!
 * 幾何情報を持たないポリゴン頂点列から幾何情報を生成
 * @param[inout] vrts 幾何情報無しのポリゴン頂点列
 * @param[out] idxs 幾何情報
 */
static void CalVertexGeometry(vector<glm::vec3> &vrts, vector< vector<int> > &idxs)
{
	int nv = (int)vrts.size();
	int np = nv/3;	// 三角形ポリゴン数

	// 幾何情報
	idxs.resize(np);
	for(int i = 0; i < np; ++i) idxs[i].resize(3);

	int same_vrts[256];	// 重複頂点のインデックスを格納

	vector<glm::vec3> compacted_vrts;	// 重複なし頂点列を格納するコンテナ

	// 重複除去処理済みフラグ格納コンテナの確保と初期化
	vector<int> find_vrts;
	find_vrts.resize(nv);
	for(int i = 0; i < nv; ++i){
		find_vrts[i] = 0;
	}

	int nvrts = 0;
	int avg_cnt = 0;

	double eps2 = 1e-4;
	for(int i = 0; i < nv; ++i){
		if(find_vrts[i]) continue;
		
		glm::vec3 pos0 = vrts[i];

		// まだ重複検索されていない頂点を格納
		compacted_vrts.push_back(pos0);
		idxs[i/3][i%3] = nvrts;

		same_vrts[0] = i;
		int cnt = 1;
		find_vrts[i] = 1;
		for(int j = 0; j < nv; ++j){
			if(find_vrts[j]) continue;

			glm::vec3 pos1 = vrts[j];

			// 距離が近い点を重複頂点とする
			if(glm::length2(pos0-pos1) < eps2){
				find_vrts[j] = 1;
				same_vrts[cnt] = j;

				idxs[j/3][j%3] = nvrts;	// compacted_vrts内のインデックス

				cnt++;	// 重複数をカウント
			}
		}

		avg_cnt += cnt;
		nvrts++;
	}

	avg_cnt /= nvrts;

	vrts = compacted_vrts;
}

/*!
 * オイラー角から回転行列に変換
 * @param[in] ang オイラー角
 * @param[out] mat 回転行列
 */
inline void EulerToMatrix(glm::vec3 ang, float mat[9])
{
	ang[1] = glm::radians(ang[1]);
	ang[0] = glm::radians(ang[0]);
	ang[2] = glm::radians(ang[2]);

	float cy = cos(ang[1]); 
	float sy = sin(ang[1]); 
	float cp = cos(ang[0]); 
	float sp = sin(ang[0]); 
	float cr = cos(ang[2]);
	float sr = sin(ang[2]);

	float cc = cy*cr; 
	float cs = cy*sr; 
	float sc = sy*cr; 
	float ss = sy*sr;

	mat[0] = cc+sp*ss;
	mat[1] = cs-sp*sc;
	mat[2] = -sy*cp;

	mat[3] = -cp*sr;
	mat[4] = cp*cr;
	mat[5] = -sp;

	mat[6] = sc-sp*cs;
	mat[7] = ss+sp*cc;
	mat[8] = cy*cp;
}

/*!
 * オイラー角から回転行列に変換
 * @param[in] ang オイラー角
 * @param[out] mat 回転行列
 */
inline void EulerToMatrix16(glm::vec3 ang, float mat[16])
{
	ang[1] = glm::radians(ang[1]);
	ang[0] = glm::radians(ang[0]);
	ang[2] = glm::radians(ang[2]);

	float cy = cos(ang[1]); 
	float sy = sin(ang[1]); 
	float cp = cos(ang[0]); 
	float sp = sin(ang[0]); 
	float cr = cos(ang[2]);
	float sr = sin(ang[2]);

	float cc = cy*cr; 
	float cs = cy*sr; 
	float sc = sy*cr; 
	float ss = sy*sr;

	mat[0]  = cc+sp*ss;
	mat[1]  = cs-sp*sc;
	mat[2]  = -sy*cp;

	mat[4]  = -cp*sr;
	mat[5]  = cp*cr;
	mat[6]  = -sp;

	mat[8]  = sc-sp*cs;
	mat[9]  = ss+sp*cc;
	mat[10] = cy*cp;
}

/*!
 * 頂点列をAABBに合うようにFitさせる(回転有り，アスペクト比無視)
 * @param[inout] poly ポリゴン
 * @param[in] cen AABB中心座標
 * @param[in] ext AABBの辺の長さ(1/2)
 * @param[in] ang 回転ベクトル
 */
static bool AffineVertices(rxPolygons &polys, glm::vec3 cen, glm::vec3 ext, glm::vec3 ang)
{
	int vn = (int)polys.vertices.size();
	if(vn <= 1) return false;

	// 現在のBBoxの大きさを調べる
	glm::vec3 minp, maxp;
	minp = maxp = polys.vertices[0];
	for(int i = 1; i < vn; ++i){
		glm::vec3 pos = polys.vertices[i];
		for(int i = 0; i < 3; ++i){
			if(pos[i] > maxp[i]) maxp[i] = pos[i];
			if(pos[i] < minp[i]) minp[i] = pos[i];
		}
	}
	
	glm::vec3 scale = (maxp-minp);
	glm::vec3 trans = (maxp+minp)/2.0f;

	for(int i = 0; i < 3; ++i){
		if(fabs(scale[i]) < RX_FEQ_EPS){
			scale[i] = 1.0;
		}
	}


	float mat[9];
	EulerToMatrix(-ang, mat);
	for(int i = 0; i < vn; ++i){
		glm::vec3 pos = polys.vertices[i];
		glm::vec3 nrm = polys.normals[i];

		glm::vec3 pos1 = ((pos-trans)/scale)*2.0f*ext;
		glm::vec3 nrm1 = nrm;

		pos[0] = mat[0]*pos1[0]+mat[1]*pos1[1]+mat[2]*pos1[2];
		pos[1] = mat[3]*pos1[0]+mat[4]*pos1[1]+mat[5]*pos1[2];
		pos[2] = mat[6]*pos1[0]+mat[7]*pos1[1]+mat[8]*pos1[2];
		nrm[0] = mat[0]*nrm1[0]+mat[1]*nrm1[1]+mat[2]*nrm1[2];
		nrm[1] = mat[3]*nrm1[0]+mat[4]*nrm1[1]+mat[5]*nrm1[2];
		nrm[2] = mat[6]*nrm1[0]+mat[7]*nrm1[1]+mat[8]*nrm1[2];

		pos += cen;

		polys.vertices[i] = pos;
		polys.normals[i]  = nrm;
	}

	return true;
}


/*!
 * 頂点列からのAABBの検索
 * @param[out] minp,maxp AABBの最大座標，最小座標
 * @param[in] vec_set 頂点列
 * @param[in] start_index,end_index 頂点列の検索範囲
 * @return 検索できたらtrue
 */
inline bool FindBBox(glm::vec3 &minp, glm::vec3 &maxp, 
					 const vector<glm::vec3> &vec_set, 
					 const int start_index, const int end_index)
{
	if((int)vec_set.size() == 0) return false;

	maxp = vec_set[start_index];
	minp = vec_set[start_index];

	for(int i = start_index+1; i < end_index; ++i){
		if(vec_set[i][0] > maxp[0]) maxp[0] = vec_set[i][0];
		if(vec_set[i][1] > maxp[1]) maxp[1] = vec_set[i][1];
		if(vec_set[i][2] > maxp[2]) maxp[2] = vec_set[i][2];
		if(vec_set[i][0] < minp[0]) minp[0] = vec_set[i][0];
		if(vec_set[i][1] < minp[1]) minp[1] = vec_set[i][1];
		if(vec_set[i][2] < minp[2]) minp[2] = vec_set[i][2];
	}

	return true;
}
/*!
 * 頂点列からのAABBの検索
 * @param[out] minp,maxp AABBの最大座標，最小座標
 * @param[in] vec_set 頂点列
 * @return 検索できたらtrue
 */
inline bool FindBBox(glm::vec3 &minp, glm::vec3 &maxp, const vector<glm::vec3> &vec_set)
{
	return FindBBox(minp, maxp, vec_set, 0, (int)vec_set.size());
}

/*!
 * 頂点列をAABBに合うようにFitさせる
 * @param[in] ctr AABB中心座標
 * @param[in] sl  AABBの辺の長さ(1/2)
 * @param[in] vec_set 頂点列
 * @param[in] start_index,end_index 頂点列の検索範囲
 */
inline bool FitVertices(const glm::vec3 &ctr, const glm::vec3 &sl, 
						vector<glm::vec3> &vec_set, 
						const int start_index, const int end_index)
{
	glm::vec3 ctr0, sl0, maxp, minp;

	// 現在のBBoxの大きさを調べる
	FindBBox(minp, maxp, vec_set, start_index, end_index);
			
	sl0  = (maxp-minp)/2.0f;
	ctr0 = (maxp+minp)/2.0f;

	int max_axis = ( ( (sl0[0] > sl0[1]) && (sl0[0] > sl0[2]) ) ? 0 : ( (sl0[1] > sl0[2]) ? 1 : 2 ) );
	int min_axis = ( ( (sl0[0] < sl0[1]) && (sl0[0] < sl0[2]) ) ? 0 : ( (sl0[1] < sl0[2]) ? 1 : 2 ) );
	float size_conv = sl[max_axis]/sl0[max_axis];

	// 全ての頂点をbboxにあわせて変換
	for(int i = start_index; i < end_index; ++i){
		vec_set[i] = (vec_set[i]-ctr0)*size_conv+ctr;
	}

	return true;
}

/*!
 * 頂点列をAABBに合うようにFitさせる
 * @param[in] ctr AABB中心座標
 * @param[in] sl  AABBの辺の長さ(1/2)
 * @param[in] vec_set 頂点列
 */
inline bool FitVertices(const glm::vec3 &ctr, const glm::vec3 &sl, vector<glm::vec3> &vec_set)
{
	FitVertices(ctr, sl, vec_set, 0, (int)vec_set.size());
	return true;
}

/*!
 * シミュレーション空間を覆うグリッドの算出
 *  - 最大グリッド数nmaxを指定して，最大辺がそのグリッド数になるように他の辺を分割する．
 *  - グリッドは立方体セルで構成される．
 * @param[in] minp,maxp シミュレーション空間の範囲
 * @param[in] nmax 最大グリッド数
 * @param[out] h 立方体セルの幅
 * @param[out] n 各軸のセル数
 * @param[in] extend 空間を完全に覆うための拡張幅[0,1]
 * @return 成功ならば0
 */
inline int CalMeshDiv(glm::vec3 &minp, const glm::vec3 maxp, const int nmax, double &h, int n[3], float extend)
{
	glm::vec3 l = maxp-minp;
	minp -= 0.5f*extend*l;
	l *= 1.0+extend;

	float max_l = 0;
	int max_axis = 0;
	for(int i = 0; i < 3; ++i){
		if(l[i] > max_l){
			max_l = l[i];
			max_axis = i;
		}
	}

	h = max_l/nmax;
	for(int i = 0; i < 3; ++i){
		n[i] = (int)(l[i]/h)+1;
	}

	return 0;
}


/*!
 * シミュレーション空間を覆うグリッドの算出
 *  - グリッド幅hを指定
 *  - グリッドは立方体セルで構成される．
 * @param[in] minp,maxp シミュレーション空間の範囲
 * @param[in] nmax 最大グリッド数
 * @param[out] h 立方体セルの幅
 * @param[out] n 各軸のセル数
 * @param[in] extend 空間を完全に覆うための拡張幅[0,1]
 * @return 成功ならば0
 */
inline int CalMeshDiv(glm::vec3 &minp, glm::vec3 maxp, double h, int n[3], float extend)
{
	glm::vec3 l = maxp-minp;
	minp -= 0.5f*extend*l;
	l *= 1.0+extend;

	for(int i = 0; i < 3; ++i){
		n[i] = (int)(l[i]/h)+1;
	}

	return 0;
}

/*!
 * メッシュ描画用のVBOを確保
 * @param[in] max_verts 最大頂点数
 * @param[in] uVrtVBO 頂点FBO
 * @param[in] uNrmVBO 法線FBO
 * @param[in] uTriVBO メッシュFBO
 */
inline bool AssignArrayBuffers(int max_verts, int dim, GLuint &uVrtVBO, GLuint &uNrmVBO, GLuint &uTriVBO)
{
	// 頂点VBO
	if(!uVrtVBO) glGenBuffers(1, &uVrtVBO);
	glBindBuffer(GL_ARRAY_BUFFER, uVrtVBO);
	glBufferData(GL_ARRAY_BUFFER, max_verts*dim*sizeof(float), 0, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// 法線VBO
	if(!uNrmVBO) glGenBuffers(1, &uNrmVBO);
	glBindBuffer(GL_ARRAY_BUFFER, uNrmVBO);
	glBufferData(GL_ARRAY_BUFFER, max_verts*dim*sizeof(float), 0, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// メッシュVBO
	if(!uTriVBO) glGenBuffers(1, &uTriVBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, uTriVBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, max_verts*3*3*sizeof(unsigned int), 0, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	return true;
}


/*!
 * FBOにデータを設定
 * @param[in] uVrtVBO 頂点FBO
 * @param[in] uNrmVBO 法線FBO
 * @param[in] uTriVBO メッシュFBO
 */
inline bool SetFBOFromArray(GLuint uVrtVBO, GLuint uNrmVBO, GLuint uTriVBO, 
							vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, vector<rxFace> &face)
{
	int nv = (int)vrts.size();
	int nm = (int)face.size();
	int nn = nv;

	// 頂点アレイに格納
	glBindBuffer(GL_ARRAY_BUFFER, uVrtVBO);
	float *vrt_ptr = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

	for(int i = 0; i < nv; ++i){
		for(int j = 0; j < 3; ++j){
			vrt_ptr[3*i+j] = (float)vrts[i][j];
		}
	}

	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// 法線情報の取得
	glBindBuffer(GL_ARRAY_BUFFER, uNrmVBO);
	float *nrm_ptr = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);

	for(int i = 0; i < nn; ++i){
		for(int j = 0; j < 3; ++j){
			nrm_ptr[3*i+j] = (float)nrms[i][j];
		}
	}

	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// 接続情報の取得
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, uTriVBO);
	unsigned int *tri_ptr = (unsigned int*)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_WRITE_ONLY);

	for(int i = 0; i < nm; ++i){
		for(int j = 0; j < 3; ++j){
			tri_ptr[3*i+j] = face[i][j];
		}
	}

	glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	return true;
}

/*!
 * VBOからホスト側配列にデータを転送
 * @param[in] uVrtVBO 頂点データVBO
 * @param[in] uNrmVBO 法線データVBO
 * @param[in] uTriVBO メッシュデータVBO
 * @param[out] vrts 頂点座標
 * @param[out] nrms 頂点法線
 * @param[out] face メッシュ
 * @param[in] nvrts 頂点数
 * @param[in] ntris 三角形ポリゴン数
 * @param[in] dim 頂点/法線の次元(3 or 4)
 */
inline bool SetArrayFromFBO(GLuint uVrtVBO, GLuint uNrmVBO, GLuint uTriVBO, 
							vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, vector<rxFace> &face, int nvrts, int ntris, int dim)
{
	// 頂点アレイに格納
	glBindBuffer(GL_ARRAY_BUFFER, uVrtVBO);
	float *vrt_ptr = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);

	vrts.resize(nvrts);
	for(int i = 0; i < nvrts; ++i){
		vrts[i][0] = vrt_ptr[dim*i];
		vrts[i][1] = vrt_ptr[dim*i+1];
		vrts[i][2] = vrt_ptr[dim*i+2];
	}

	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// 法線情報の取得
	glBindBuffer(GL_ARRAY_BUFFER, uNrmVBO);
	float *nrm_ptr = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);

	nrms.resize(nvrts);
	for(int i = 0; i < nvrts; ++i){
		nrms[i][0] = nrm_ptr[dim*i];
		nrms[i][1] = nrm_ptr[dim*i+1];
		nrms[i][2] = nrm_ptr[dim*i+2];
	}

	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	// 接続情報の取得
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, uTriVBO);
	unsigned  *tri_ptr = (unsigned int*)glMapBuffer(GL_ELEMENT_ARRAY_BUFFER, GL_READ_ONLY);

	face.resize(ntris);
	for(int i = 0; i < ntris; ++i){
		face[i].vert_idx.resize(3);
		face[i][0] = tri_ptr[3*i];
		face[i][1] = tri_ptr[3*i+1];
		face[i][2] = tri_ptr[3*i+2];
	}

	glUnmapBuffer(GL_ELEMENT_ARRAY_BUFFER);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	return true;
}



//-----------------------------------------------------------------------------
// エッジ関係
//-----------------------------------------------------------------------------

/*!
* 削除フラグ付きのポリゴンを削除
* @param[inout] src 変更するポリゴン
* @param[in] del_attr 削除する面に登録された属性
*/
static inline int DeletePolygon(rxPolygonsE *src, int del_attr)
{
	vector<rxFace>::iterator iter = src->faces.begin();
	int d = 0;
	while(iter != src->faces.end()){
		if(iter->attribute == del_attr){
			iter = src->faces.erase(iter);
			d++;
		}
		else{
			++iter;
		}
	}
	return d;
}


/*!
* アンブレラオペレータで頂点をFairing
* @param[inout] obj 拡張ポリゴンオブジェクト
* @return 成功の場合1
*/
static inline int VertexFairingByUmbrella(rxPolygonsE &obj)
{
	int vnum = (int)obj.vertices.size();
	if(!vnum || obj.vedges.empty()) return 0;

	for(int i = 0; i < vnum; ++i){
		glm::vec3 newpos = obj.vertices[i];

		// 頂点につながるエッジを走査
		int en = 0;
		set<int>::iterator eiter = obj.vedges[i].begin();
		do{
			int v1, v2;
			v1 = obj.edges[*eiter].v[0];
			v2 = obj.edges[*eiter].v[1];

			newpos += obj.vertices[(v1 != i ? v1 : v2)];
			en++;
		}while(++eiter != obj.vedges[i].end());

		newpos /= (double)(en+1.0);

		obj.vertices[i] = newpos;
	}

	return 1;
}

/*!
* 頂点に接続するポリゴンを探索
* @param[inout] obj 拡張ポリゴンオブジェクト
* @return 成功の場合1
*/
static inline int SearchVertexFace(rxPolygonsE &obj)
{
	int vnum = (int)obj.vertices.size();
	if(!vnum) return 0;

	obj.vfaces.clear();
	obj.vfaces.resize(vnum);

	rxFace* face;
	int pnum = (int)obj.faces.size();

	// 全ポリゴンを探索
	for(int i = 0; i < pnum; ++i){
		face = &obj.faces[i];
		for(int j = 0; j < face->size(); ++j){
			int v = face->at(j);
			obj.vfaces[v].insert(i);
		}
	}

	return 1;
}

/*!
* エッジデータを作成
* @param[inout] obj 拡張ポリゴンオブジェクト
* @return エッジ数
*/
static inline int SearchEdge(rxPolygonsE &obj)
{
	int edge_count = 0;
	int vnum = (int)obj.vertices.size();
	int pnum = (int)obj.faces.size();

	obj.vedges.clear();
	obj.vedges.resize(vnum);
	obj.fedges.clear();
	obj.fedges.resize(pnum);
	obj.edges.clear();

	rxFace* face;
	int vert_idx[2];


	// 全ポリゴンを探索
	for(int i = 0; i < pnum; ++i){
		face = &obj.faces[i];
		vnum = face->size();

		for(int j = 0; j < vnum; ++j){
			// エッジ頂点1
			vert_idx[0] = face->at(j);
			// エッジ頂点2
			vert_idx[1] = (j == vnum-1) ? face->at(0) : face->at(j+1);

			// 重複する稜線のチェック
			bool overlap = false;
			set<int>::iterator viter;

			// エッジ頂点1のチェック
			if((int)obj.vedges[vert_idx[0]].size() > 0){
				// 頂点につながっているエッジ(既に登録されたもの)を調べる
				viter = obj.vedges[vert_idx[0]].begin();
				do{
					int v1, v2;
					v1 = obj.edges[*viter].v[0];
					v2 = obj.edges[*viter].v[1];

					if( (vert_idx[0] == v1 && vert_idx[1] == v2) || (vert_idx[1] == v1 && vert_idx[0] == v2) ){
						overlap = true;
						obj.edges[*viter].f.insert(i);	// 面をエッジに加える
						obj.fedges[i].insert(*viter);	// エッジを面に加える
					}
				}while(++viter != obj.vedges[vert_idx[0]].end());
			}

			// エッジ頂点2のチェック
			if((int)obj.vedges[vert_idx[1]].size() > 0){
				// 頂点につながっているエッジ(既に登録されたもの)を調べる
				viter = obj.vedges[vert_idx[1]].begin();
				do{
					int v1, v2;
					v1 = obj.edges[*viter].v[0];
					v2 = obj.edges[*viter].v[1];

					if( (vert_idx[0] == v1 && vert_idx[1] == v2) || (vert_idx[1] == v1 && vert_idx[0] == v2) ){
						overlap = true;
						obj.edges[*viter].f.insert(i);	// 面をエッジに加える
						obj.fedges[i].insert(*viter);	// エッジを面に加える
					}
				}while(++viter != obj.vedges[vert_idx[1]].end());
			}

			// 重複する稜線なしの場合，新規に追加
			if(!overlap){
				obj.edges.push_back(rxEdge());
				obj.edges[edge_count].v[0] = vert_idx[0];
				obj.edges[edge_count].v[1] = vert_idx[1];
				obj.edges[edge_count].f.insert(i);

				obj.vedges[vert_idx[0]].insert(edge_count);
				obj.vedges[vert_idx[1]].insert(edge_count);
				obj.fedges[i].insert(edge_count);

				edge_count++;
			}

		}
	}

	return edge_count;
}


/*!
* エッジ長さを計算
* @param[inout] obj 拡張ポリゴンオブジェクト
* @return エッジ数
*/
static inline int CalEdgeLength(rxPolygonsE &obj)
{
	int en = (int)obj.edges.size();

	if(en == 0) return 0;

	int vn = (int)obj.vertices.size();

	rxEdge *edge;
	for(int i = 0; i < en; ++i){
		edge = &obj.edges[i];

		int v0 = edge->v[0];
		int v1 = edge->v[1];

		edge->len = glm::length(obj.vertices[v0]-obj.vertices[v1]);
	}

	return 1;
}




/*!
* アンブレラオペレータでSkinning weightをFairing
* @param[inout] obj 拡張ポリゴンオブジェクト
* @return 成功の場合1
*/
static inline int WeightFairingByUmbrella(rxPolygonsE &obj)
{
	int vnum = (int)obj.vertices.size();
	if(!vnum || obj.vedges.empty()) return 0;

	for(int i = 0; i < vnum; ++i){
		int n = static_cast<int>(obj.weights[i].size());	// 重みが設定されたボーンの数
		if(n <= 1) continue;

		// 	vector< map<int, double> > weights;	//!< 各頂点のSkinning weight (対応するボーン番号(int型)と重み(double型)のマップ)

		map<int, double>::iterator itr = obj.weights[i].begin(), itr2;
		for(; itr != obj.weights[i].end(); ++itr){
			int bone = itr->first;
			int nw = 0;
			double total_w = 0.0;

			// 頂点につながるエッジを走査
			int en = 0;
			set<int>::iterator eiter = obj.vedges[i].begin();
			do{
				int v1, v2;
				v1 = obj.edges[*eiter].v[0];
				v2 = obj.edges[*eiter].v[1];
				int v = (v1 != i ? v1 : v2);
				if((itr2 = obj.weights[v].find(bone)) != obj.weights[v].end()){
					total_w += itr2->second;
					nw++;
				}

				en++;
			}while(++eiter != obj.vedges[i].end());

			itr->second = (total_w+itr->second)/(en+1.0);
		}
	}

	return 1;
}


/*!
* Skinning weightの正規化
* @param[in] obj 拡張ポリゴンオブジェクト
* @return 成功の場合1
*/
static inline int NormalizeWeights(rxPolygonsE &obj)
{
	int vnum = (int)obj.vertices.size();
	if(!vnum || obj.vedges.empty()) return 0;

	for(int i = 0; i < vnum; ++i){
		double total_w = 0.0;
		map<int, double>::iterator itr;
		for(itr = obj.weights[i].begin(); itr != obj.weights[i].end(); ++itr){
			total_w += itr->second;
		}
		if(fabs(total_w) < 1.0e-6) total_w = 1.0;
		for(itr = obj.weights[i].begin(); itr != obj.weights[i].end(); ++itr){
			itr->second /= total_w;
		}
	}

	return 1;
}


/*!
* Skinning weightのテキスト出力(デバッグ用)
* @param[in] obj 拡張ポリゴンオブジェクト
* @param[in] filename 出力ファイル名
* @return 成功の場合1
*/
static inline int DebugWeights(const rxPolygonsE &obj, string filename)
{
	int vnum = (int)obj.vertices.size();
	if(!vnum || obj.vedges.empty()) return 0;

	// ファイル出力
	ofstream fo;
	fo.open(filename.c_str());

	fo << "[weights]" << endl;

	for(int i = 0; i < vnum; ++i){
		fo << "v" << i << " : ";
		map<int, double>::const_iterator itr = obj.weights[i].begin();
		for(; itr != obj.weights[i].end(); ++itr){
			int bone = itr->first;
			double weight = itr->second;
			fo << weight << "(" << bone << "), ";
		}
		fo << endl;
	}

	fo.close();


	return 1;
}






//海老沢追加
//ソートに用いる
inline bool edge_idx_compare(glm::ivec2 e1, glm::ivec2  e2) { return e1.x < e2.x; }

//エッジと頂点座標のみを出力
static inline bool SaveInterObj(string file_name,vector<Edge_Vert> inter_objs) {
	vector<glm::vec3> vertices;
	vector<rxEdge> edges;
	for (int i = 0; i < inter_objs.size(); i++) {
		Edge_Vert inter_obj = inter_objs[i];
		int prev_vert_size = vertices.size();

		for (int j = 0; j < inter_obj.vertices.size(); j++) {
			vertices.push_back(inter_obj.vertices[j]);
		}
		for (int j = 0; j < inter_obj.edges.size(); j++) {
			rxEdge edge = inter_obj.edges[j];
			edge.v[0] += prev_vert_size;
			edge.v[1] += prev_vert_size;
			edges.push_back(edge);
		}
	}

	ofstream file;

	file.open(file_name.c_str());
	if (!file || !file.is_open() || file.bad() || file.fail()) {
		cout << "rxOBJ::Save : Invalid file specified" << endl;
		return false;
	}

	// 頂点(v)
	int nv = (int)vertices.size();
	for (int i = 0; i < nv; ++i) {
		file << "v " << Vec3ToString(glm::vec3(vertices[i])) << endl;
	}

	//エッジ
	vector<glm::ivec2> edge_idx;//エッジを構成する二つのエッジを保持
	//エッジを構成する粒子のインデックスをivec(small,large)の形で保持
	int ne = (int)edges.size();
	for (int i = 0; i < ne; i++) {
		//元の値がインデックスなので+1しておく
		int v1 = edges[i].v[0] + 1;
		int v2 = edges[i].v[1] + 1;
		if (v1 < v2) {
			edge_idx.push_back(glm::ivec2(v1, v2));
		}
		else {
			edge_idx.push_back(glm::ivec2(v2, v1));
		}
	}

	//ソート処理が必要
	std::sort(edge_idx.begin(), edge_idx.end(), edge_idx_compare);

	//ファイルへの書き出し
	for (int i = 0; i < edge_idx.size(); i++) {
		file << "l " << edge_idx[i].x << " " << edge_idx[i].y << endl;
	}
}

//一つのメッシュのみをobjファイルに出力
static inline bool SaveOneMesh(string file_name, Edge_Vert Mesh) {
	ofstream file;

	file.open(file_name.c_str());
	if (!file || !file.is_open() || file.bad() || file.fail()) {
		cout << "rxOBJ::Save : Invalid file specified" << endl;
		return false;
	}

	// 頂点(v)
	int nv = (int)Mesh.vertices.size();
	for (int i = 0; i < nv; ++i) {
		file << "v " << Vec3ToString(glm::vec3(Mesh.vertices[i])) << endl;
	}

	//エッジ
	vector<glm::ivec2> edge_idx;//エッジを構成する二つのエッジを保持
	//エッジを構成する粒子のインデックスをivec(small,large)の形で保持
	int ne = (int)Mesh.edges.size();
	for (int i = 0; i < ne; i++) {
		//元の値がインデックスなので+1しておく
		int v1 = Mesh.edges[i].v[0] + 1;
		int v2 = Mesh.edges[i].v[1] + 1;
		if (v1 < v2) {
			edge_idx.push_back(glm::ivec2(v1, v2));
		}
		else {
			edge_idx.push_back(glm::ivec2(v2, v1));
		}
	}

	//ソート処理が必要
	std::sort(edge_idx.begin(), edge_idx.end(), edge_idx_compare);

	//ファイルへの書き出し
	for (int i = 0; i < edge_idx.size(); i++) {
		file << "l " << edge_idx[i].x << " " << edge_idx[i].y << endl;
	}
}





//内積配列のソートに用いる
inline bool dot_compare(std::pair<float,int> e1,std::pair<float,int> e2) { return e1.first < e2.first; }

//境界エッジの出力
static inline int SearchBoundaryEdge(rxPolygonsE& obj) {
	int num_e = (int)obj.edges.size();
	vector<int> boundary_edges;//境界エッジのインデックスを1次元配列で管理
	for (int i = 0; i < num_e; i++) {
		rxEdge edge = obj.edges[i];
		if (edge.f.size() == 1) {//エッジに含まれるポリゴンが一つだけなら，このエッジを追加
			/*for (auto j = edge.f.begin(); j != edge.f.end(); j++) {
				cout << "idx " << i << " f" << *j << endl;
			}*/
			boundary_edges.push_back(i);
		}
	}
	cout << "boundary edge count" << boundary_edges.size() << endl;

	int num_v = (int)obj.vertices.size();
	//境界エッジのみを含んだ，各頂点が持つエッジのインデックス(境界エッジを持たない頂点ではnullを持つ)
	vector<set<int>> new_vedges(num_v);

	//頂点から境界エッジに含まれるエッジのインデックスを追加
	for (int i = 0; i < num_v; i++) {
		set<int> edge_idxs= obj.vedges[i];
		//*itrがエッジのインデックス
		for (set<int>::iterator itr = edge_idxs.begin(); itr!=edge_idxs.end(); itr++) {
			//現在の粒子が持つエッジのインデックスが境界粒子の配列に含まれている場合
			if (find(boundary_edges.begin(), boundary_edges.end(), *itr) != boundary_edges.end()) {
				//cout << "add new vedges " << *itr << endl;
				new_vedges[i].insert(*itr);
			}
		}
	}
	cout << "edge_idx restored" << endl;

	vector<Edge_Vert> Meshes;
	int  x= 0;

	//上記でインデックスは取り出したが，実体はobjから取り出す
	while (!boundary_edges.empty()) {
		//ここにデータを記憶(一つの閉形式になるはず)
		Edge_Vert Mesh;

		int edge_idx = boundary_edges[0];
		//v_initを始点とする
		int v_init = obj.edges[edge_idx].v[0];

		//頂点を追加
		Mesh.vertAdd(obj.vertices[v_init]);
		//cout << "vert_id " << v_init << " add to Mesh" << endl;

		int v1, v2;
		//常にv2が次の頂点になるようにする
		v2 = obj.edges[edge_idx].v[1];
		//始点と一致するものが見つかるまでループする
		while (v_init != v2) {
			//-------------------------------------
			//edgeを記憶(実態はobj.edgesからとるが，v[0],v[1]のインデックスは変更される
			Mesh.edgeAdd(obj.edges[edge_idx]);
			//cout << "edge_id " << edge_idx << " add to Mesh" << endl;

			//edge_idxに一致する境界配列は削除する
			for (int i = 0; i < boundary_edges.size(); i++) {
				if (boundary_edges[i] == edge_idx) {
					boundary_edges.erase(boundary_edges.begin() + i);
					//cout << *(boundary_edges.begin() + i) << " erased!!" << endl;
					break;
				}
			}

			//次の頂点を記憶
			Mesh.vertAdd(obj.vertices[v2]);
			//cout << "vert_id " << v2 << " add to Mesh" << endl;
			//------------------------------------

			//頂点ごとに2組の境界エッジを持つことを前提としている
			set<int> edges = new_vedges[v2];
			auto itr = edges.begin();
			int first = *itr;
			itr++;
			int second = *itr;
			//cout << "first: " << first << " second: " << second << endl;
			//2組のうち，現在と異なるエッジ番号を選択する
			if (edge_idx == first) edge_idx = second;
			else edge_idx = first;

			//現在のv2と一致しない方を方を次のv2とする
			if (v2 == obj.edges[edge_idx].v[0]) {
				v1 = obj.edges[edge_idx].v[0];
				v2 = obj.edges[edge_idx].v[1];
			}
			else if(v2==obj.edges[edge_idx].v[1]) {
				v1 = obj.edges[edge_idx].v[1];
				v2 = obj.edges[edge_idx].v[0];
			}
			else {
				std::cerr << "Edge_Vert Error!! No Vertices can Add! Vertices are not included now Edge!!------------------------------------------------------------" << endl;
			}
		}
		rxEdge e;
		//最後のエッジを追加
		Mesh.edgeAdd(e);
		Mesh.vertAdd(obj.vertices[v2]);
		//最後のエッジを削除
		for (int i = 0; i < boundary_edges.size(); i++) {
			if (boundary_edges[i] == edge_idx) {
				boundary_edges.erase(boundary_edges.begin() + i);
				//cout << *(boundary_edges.begin() + i) << " erased!!" << endl;
				break;
			}
			cout << "Not Erased" << endl;
		}

		//メッシュを保存
		x++;
		Meshes.push_back(Mesh);
		cout << "Mesh Added! vertices size " << Mesh.vertices.size() << " edge size" << Mesh.edges.size() << " id: " << x << endl;
	}
	//edge_idxに一致する境界配列は削除する
	cout << "Mesh Add Finished!!" << endl;
	//メモリの開放
	boundary_edges.clear();
	new_vedges.clear();
	obj.Clear();

	//閉形式のメッシュを出力---------------------------------------------------------------
	for (int j = 0; j < Meshes.size(); j++) {
		string filename = "AssetsNotUsed/DebugMeshes/Mesh" + to_string(j) + ".obj";
		SaveOneMesh(filename, Meshes[j]);
	}
	//--------------------------------------------------------------------------------------


	//最終的な髪の毛の集合
	vector<Edge_Vert> Strands;

	//境界のみを持った閉形式のエッジと頂点の集合から縦方向のみのエッジを抽出
	for (int i = 0; i < Meshes.size(); i++) {
		vector<rxEdge> edges = Meshes[i].edges;
		vector<glm::vec3> vertices = Meshes[i].vertices;
		if (edges.empty() || vertices.empty()) {
			std::cerr << "Error:vertices or edges are empty!!------------------------------------------------------------------------------------------" << endl;
		}
		cout << "id " << i << "start" << endl;
		//内積の配列を作る(2つ目の要素にはインデックスを格納)
		vector<std::pair <float, int>> dot_array;

		for (int j = 0; j < vertices.size(); j++) {
			//一つ前の頂点との方向ベクトルだが，エッジは繰り返すので，最初と最後はオーバーフローをしないように
			glm::vec3 prev_dir, next_dir;
			//前の頂点との方向ベクトル
			if (j == 0) prev_dir = glm::normalize(vertices[vertices.size()-1] - vertices[0]);//最後のインデックスと最初のインデックス
			else prev_dir = glm::normalize(vertices[j - 1] - vertices[j]);

			//次の頂点との方向ベクトル
			if (j == vertices.size() - 1) next_dir = (vertices[0] - vertices[vertices.size() - 1]);
			else next_dir = glm::normalize(vertices[j + 1] - vertices[j]);
			//絶対値のみを比較したい
			dot_array.push_back(make_pair(abs(glm::dot(prev_dir, next_dir)), j));
		}
		cout << "get dot! id " << i << endl;

		//ソートして，最小から4つのインデックスを取る
		sort(dot_array.begin(), dot_array.end(), dot_compare);



		//距離が短いエッジを消す場合--------------------------------------------------------------------------------------------------------------
		//vector<int> min_idx;
		//for (int j = 0; j < 4; j++) {
		//	//最小の値を持つ4つのインデックス
		//	min_idx.push_back(dot_array[j].second);
		//	if (dot_array[j].first > 0.1) cout << "dot too high! id " << i << " value " << dot_array[j].first << endl;
		//}

		////切り捨てるエッジの条件を非常に簡略化し，インデックス間の数が少ない部分を切り捨てる
		////インデックスを小さい順に並べ替える(インデックスであり，値の大小は並び順)
		//sort(min_idx.begin(), min_idx.end());

		//int size = vertices.size();
		////インデックス間の距離(4つのインデックスの間にあるインデックスの数の差)
		//vector<int> dist{ min_idx[1] - min_idx[0], min_idx[2] - min_idx[1], min_idx[3] - min_idx[2], size + min_idx[0] - min_idx[3] };

		////最小距離とその対角にあるインデックス間のエッジを消す
		//auto iterator = std::min_element(dist.begin(), dist.end());
		//int id = std::distance(dist.begin(), iterator);

		////エッジを消す(11-12,21-22間のエッジを残す)
		//int min_idx11, min_idx12, min_idx21, min_idx22;
		//if (id == 0 || id == 2) {
		//	//1-2,3-0のエッジを残す
		//	min_idx11 = min_idx[1];
		//	min_idx12 = min_idx[2];

		//	min_idx21 = min_idx[3];
		//	min_idx22 = min_idx[0] + vertices.size();
		//}
		//else {
		//	//0-1,2-3のエッジを残す
		//	min_idx11 = min_idx[0];
		//	min_idx12 = min_idx[1];

		//	min_idx21 = min_idx[2];
		//	min_idx22 = min_idx[3];
		//}

		//Edge_Vert Strand1;
		//for (int k = min_idx11; k < min_idx12; k++) {
		//	Strand1.vertAdd(vertices[k]);
		//	Strand1.edgeAdd(edges[k]);
		//}
		//Strand1.vertAdd(vertices[min_idx12]);

		//Edge_Vert Strand2;
		//for (int k = min_idx21; k < min_idx22; k++) {
		//	int tmp = k;
		//	if (tmp > vertices.size()-1) tmp -= vertices.size();
		//	Strand2.vertAdd(vertices[tmp]);
		//	//cout << "did" << endl;
		//	//頂点の配列の最大値を超える場合，エッジの配列からオーバーフローする可能性がある．
		//	if (tmp > edges.size()-1) {
		//		//この場合，適当なエッジを追加しておく
		//		rxEdge e;
		//		Strand2.edgeAdd(e);
		//	}
		//	else {
		//		Strand2.edgeAdd(edges[tmp]);
		//	}
		//}
		////配列からオーバーフローしないようにする．
		//if (min_idx22 > vertices.size() - 1) min_idx22 -= vertices.size();
		//Strand2.vertAdd(vertices[min_idx22]);

		//Strands.push_back(Strand1);
		//Strands.push_back(Strand2);
		//cout << "finish id " << i << endl;
		//------------------------------------------------------------------------------------------------------------------------------------------------

		

		//y座標で判定
		//4つの頂点のそれぞれのy座標とインデックスを保持しておく-------------------------------------------------------------------------------------------
		vector<std::pair<float, int>> y_coord;
		for (int j = 0; j < 4; j++) {
			y_coord.push_back(make_pair(vertices[dot_array[j].second].y, dot_array[j].second));
		}
		//小さい順に並び変える
		sort(y_coord.begin(), y_coord.end(), dot_compare);

		//最大の2つが固定点となる始点として，それぞれのインデックスを保持
		int start1_ind = y_coord[2].second;
		int start2_ind = y_coord[3].second;

		//その他の二つを終点として，それぞれのインデックスを保持
		int last1_ind = y_coord[0].second;
		int last2_ind = y_coord[1].second;

		Edge_Vert Strand1;
		Strand1.vertAdd(vertices[start1_ind]);
		while (start1_ind != last1_ind && start1_ind != last2_ind) {
			//cout << "Strand1 start ind" << start1_ind << endl;

			//+方向に進んだ場合
			if (start1_ind == vertices.size() - 1) start1_ind = 0;
			else start1_ind++;

			if (start1_ind == start2_ind) {
				//もう一つのy最大についた場合には，初期化して-方向に進む
				start1_ind = y_coord[2].second;
				Strand1.Clear();
				break;
			}

			rxEdge e;
			Strand1.edgeAdd(e);
			Strand1.vertAdd(vertices[start1_ind]);
		}
		//breakし，Strand1が初期化された場合，-方向に進む
		if (Strand1.vertices.empty()) {
			//改めて初期点を追加
			Strand1.vertAdd(vertices[start1_ind]);
			while (start1_ind != last1_ind && start1_ind != last2_ind) {
				//-方向に進んだ場合
				if (start1_ind == 0) start1_ind = vertices.size() - 1;
				else start1_ind--;

				rxEdge e;
				Strand1.edgeAdd(e);
				Strand1.vertAdd(vertices[start1_ind]);
			}
		}
		
		cout << "finish Strand1 id " << i << endl;

		Edge_Vert Strand2;
		Strand2.vertAdd(vertices[start2_ind]);
		while (start2_ind != last1_ind && start2_ind != last2_ind) {
			//+方向に進んだ場合
			if (start2_ind == vertices.size() - 1) start2_ind = 0;
			else start2_ind++;

			if (start2_ind == start1_ind) {
				//もう一つのy最大についた場合には，初期化して-方向に進む
				start2_ind = y_coord[3].second;
				Strand2.Clear();
				break;
			}

			rxEdge e;
			Strand2.edgeAdd(e);
			Strand2.vertAdd(vertices[start2_ind]);
		}
		//breakし，Strand2が初期化された場合，-方向に進む
		if (Strand2.vertices.empty()) {
			//改めて初期点を追加
			Strand2.vertAdd(vertices[start2_ind]);
			while (start2_ind != last1_ind && start2_ind != last2_ind) {
				//-方向に進んだ場合
				if (start2_ind == 0) start2_ind = vertices.size() - 1;
				else start2_ind--;

				rxEdge e;
				Strand2.edgeAdd(e);
				Strand2.vertAdd(vertices[start2_ind]);
			}
		}
		cout << "finish Strand2 id " << i << endl;

		Strands.push_back(Strand1);
		Strands.push_back(Strand2);
	}

	cout << "start Save" << endl;
	//出力
	SaveInterObj("AssetsNotUsed/Inter_Test.obj", Strands);

	for (int j = 0; j < Strands.size(); j++) {
		string filename = "AssetsNotUsed/AllStrands/Strand" + to_string(j) + ".obj";
		SaveOneMesh(filename, Strands[j]);
	}
}



//頂点だけを取って，閾値を設定して，これを超えたら，毛髪が分かれると仮定する．
static void MakeHairObj(string filename,rxPolygonsE& obj,float T) {
	vector <Edge_Vert> g_hair;
	Edge_Vert Strand;
	for (int i = 0; i < obj.vertices.size();i++) {
		glm::vec3 v = obj.vertices[i];
		if (Strand.vertices.empty()) {
			//最初のみの処理
			Strand.vertAdd(v);
		}
		else {
			glm::vec3 last_vert = Strand.vertices[Strand.vertices.size() - 1];
			if (glm::length(v - last_vert) > T) {//直前のエッジとの距離が閾値を超えたら，一つの毛髪の列が終わったこととする．
				g_hair.push_back(Strand);
				Strand.Clear();
			}
			else {//毛髪が連続しているならエッジを追加する
				rxEdge e;
				Strand.edgeAdd(e);
			}
			Strand.vertAdd(v);
		}
	}
	SaveInterObj("AssetsNotUsed/" + filename + "_hair.obj", g_hair);

	//ランダムに1024本を選択して，そのObjファイルを作る
	vector<Edge_Vert> limit_hairs;
	srand((unsigned int)time(NULL));
	for (int i = 0; i < 1024; i++) {
		int num = g_hair.size();
		limit_hairs.push_back(g_hair[rand() % num]);
	}

	SaveInterObj("AssetsNotUsed/" + filename + "_hair_1024.obj", limit_hairs);
}

#endif // #ifndef _RX_MESH_H_
