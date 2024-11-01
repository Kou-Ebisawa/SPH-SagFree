/*!
  @file rx_utils.h

  @brief 各種便利な関数や定数など

  @author Makoto Fujisawa
  @date 2020-07
*/
// FILE --rx_utils.h--

#ifndef _RX_UTILS_H_
#define _RX_UTILS_H_

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
// C標準
#include <ctime>
#include <cmath>
//#include <cctype>
#include <cstdio>
#include <cassert>

#ifdef WIN32
#include <direct.h>
#endif

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

// STL
#include <vector>
#include <string>
#include <map>
#include <bitset>
#include <algorithm>

// OpenGL
#include <GL/glew.h>
#include <GLFW/glfw3.h>

// glm
#include "glm/glm.hpp"
#include "glm/gtc/matrix_transform.hpp"
#include "glm/gtc/type_ptr.hpp"
#include "glm/ext.hpp"	// for glm::to_string()

// 画像読み書き
#include "rx_bitmap.h"

// OpenGLのシェーダ
#include "rx_shaders.h"

using namespace std;


//-----------------------------------------------------------------------------
// 定義
//-----------------------------------------------------------------------------
const double RX_PI = 3.14159265358979323846;	// 円周率
const double RX_FEQ_EPS = 1.0e-8;				// 許容誤差
const double RX_FEQ_INF = 1.0e10;				// 大きい数の初期値設定用
const double RX_DEGREES_TO_RADIANS = 0.0174532925199432957692369076848;	//! degree -> radian の変換係数(pi/180.0)
const double RX_RADIANS_TO_DEGREES = 57.295779513082320876798154814114;	//! radian -> degree の変換係数(180.0/pi)

//! 許容誤差を含めた等値判定
template<class T>
inline bool RX_FEQ(const T &a, const T &b){ return (fabs(a-b) < RX_FEQ_EPS); }

//! 値のクランプ(クランプした値を返す)
template<class T>
inline T RX_CLAMP(const T &x, const T &a, const T &b){ return ((x < a) ? a : (x > b) ? b : x); }

//! 1次元線型補間
template<class T>
inline T RX_LERP(const T &a, const T &b, const T &t){ return a + t*(b-a); }

//! 乱数
inline double RX_RAND(const double &_min, const double &_max){ return (_max-_min)*(double(rand())/(1.0+RAND_MAX))+_min; }

//! [0,1]の浮動小数点数乱数
inline double RX_FRAND(){ return rand()/(double)RAND_MAX; }


//! スワップ
template<class T>
inline void RX_SWAP(T &a, T &b){ T c; c = a; a = b; b = c; }


//-----------------------------------------------------------------------------
// glm関係の関数
//-----------------------------------------------------------------------------
inline int IDX4(int row, int col){ return (row | (col<<2)); }

/*!
 * 1次元配列に格納された4x4行列とglm::vec3ベクトルの掛け算
 *  | m[0]  m[1]  m[2]  m[3]  |
 *  | m[4]  m[5]  m[6]  m[7]  |
 *  | m[8]  m[9]  m[10] m[11] |
 *  | m[12] m[13] m[14] m[15] |
 * @param[in] m 元の行列
 * @param[in] v ベクトル
 * @return 積の結果のベクトル(m*v)
 */
inline glm::vec3 MulMat4Vec3(GLfloat* m, const glm::vec3& v)
{
	glm::vec3 d(0.0);
	d[0] = (v[0]*m[IDX4(0, 0)]+v[1]*m[IDX4(0, 1)]+v[2]*m[IDX4(0, 2)]);
	d[1] = (v[0]*m[IDX4(1, 0)]+v[1]*m[IDX4(1, 1)]+v[2]*m[IDX4(1, 2)]);
	d[2] = (v[0]*m[IDX4(2, 0)]+v[1]*m[IDX4(2, 1)]+v[2]*m[IDX4(2, 2)]);
	return d;
}

//! オイラー角から回転行列への変換
inline void EulerToMatrix(float *m, float pitch, float yaw, float roll)
{
	yaw   = RX_DEGREES_TO_RADIANS*yaw;
	pitch = RX_DEGREES_TO_RADIANS*pitch;
	roll  = RX_DEGREES_TO_RADIANS*roll;

	float cy = cos(yaw); 
	float sy = sin(yaw); 
	float cp = cos(pitch); 
	float sp = sin(pitch); 
	float cr = cos(roll);
	float sr = sin(roll);

	float cc = cy*cr; 
	float cs = cy*sr; 
	float sc = sy*cr; 
	float ss = sy*sr;

	m[0]  = cc+sp*ss;
	m[1]  = cs-sp*sc;
	m[2]  = -sy*cp;
	m[3]  = 0.0;

	m[4]  = -cp*sr;
	m[5]  = cp*cr; 
	m[6]  = -sp;
	m[7]  = 0.0;

	m[8]  = sc-sp*cs;
	m[9]  = ss+sp*cc;
	m[10] = cy*cp;
	m[11] = 0.0;

	m[12] = 0.0;
	m[13] = 0.0;
	m[14] = 0.0;
	m[15] = 1.0;
}




//-----------------------------------------------------------------------------
// 文字列処理
//-----------------------------------------------------------------------------
/*!
 * 様々な型のstringへの変換(stringstreamを使用)
 * @param[in] x 入力
 * @return string型に変換したもの
 */
template<typename T>
inline std::string RX_TO_STRING(const T &x)
{
	std::stringstream ss;
	ss << x;
	return ss.str();
}

//! string型に<<オペレータを設定
template<typename T>
inline std::string &operator<<(std::string &cb, const T &a)
{
	cb += RX_TO_STRING(a);
	return cb;
}

/*!
 * "\n"が含まれるstringを複数のstringに分解する
 * @param[in] org 元の文字列
 * @param[in] div 分解後の文字列配列
 */
static inline void divideString(const string& org, vector<string>& div)
{
	size_t pos0 = 0, pos1 = 0;
	while(pos1 != string::npos){
		pos1 = org.find("\n", pos0);

		div.push_back(org.substr(pos0, pos1-pos0));

		pos0 = pos1+1;
	}
}


/*!
 * 整数値の下一桁を返す
 * @param[in] x 整数値
 * @return xの下一桁
 */
inline int LowOneDigit(const int &x)
{
	int x1 = (x < 0) ? -x : x;
	float a = 10;

	//INT_MAX = 2147483647
	for(int i = 10; i >= 1; i--){
		a = pow(10.0, (float)i);
		while(x1 > a){
			x1 -= (int)a;
		}
	}

	return x1;
}

/*!
 * 0付きの数字を生成
 * @param[in] n 数字
 * @param[in] d 桁数
 * @return 0付きの数字(string)
 */
inline string GenZeroNo(int n, const int &d)
{
	string zero_no = "";
	int dn = d-1;
	if(n > 0){
		dn = (int)(log10((float)n))+1;
	}
	else if(n == 0){
		dn = 1;
	}
	else{
		n = 0;
		dn = 1;
	}

	for(int i = 0; i < d-dn; ++i){
		zero_no += "0";
	}

	zero_no += RX_TO_STRING(n);

	return zero_no;
}

/*!
 * 秒数を hh:mm:ss の形式に変換
 * @param[in] sec 秒数
 * @param[in] use_msec ミリ秒まで含める(hh:mm:ss.msec)
 * @return hh:mm:ss形式の文字列
 */
inline string GenTimeString(float sec, bool use_msec = false)
{
	long value = (int)(1000*sec+0.5);	// ミリ秒

	unsigned int h = (unsigned int)(value/3600000);	// 時間
	value -= h*3600000;
	unsigned int m = (unsigned int)(value/60000);		// 分
	value -= m*60000;
	unsigned int s = (unsigned int)(value/1000);		// 秒
	value -= s*1000;
	unsigned int ms = (unsigned int)(value);			// ミリ秒

	stringstream ss;
	if(h > 0) ss << GenZeroNo(h, 2) << ":";
	ss << GenZeroNo(m, 2) << ":";
	ss << GenZeroNo(s, 2);
	if(use_msec) ss << "." << GenZeroNo(ms, 3);

	return ss.str();
}

/*!
 * 時刻を hh:mm:ss の形式に変換
 * @param[in] h,m,s 時,分,秒
 * @return hh:mm:ss形式の文字列
 */
inline string GenTimeString(int h, int m, int s)
{
	stringstream ss;
	if(h > 0) ss << GenZeroNo(h, 2) << ":";
	ss << GenZeroNo(m, 2) << ":";
	ss << GenZeroNo(s, 2);
	return ss.str();
}


//-----------------------------------------------------------------------------
// ファイル・フォルダ処理
//-----------------------------------------------------------------------------
/*!
 * ファイル，フォルダの存在確認
 * @param[in] path ファイル・フォルダパス
 */
inline int ExistFile(const string fn)
{
	FILE *fp;

	if((fp = fopen(fn.c_str(), "r")) == NULL){
		return 0;
	}

	fclose(fp);
	return 1;
}

/*!
 * フォルダ区切りの検索
 * @param[in] str ファイル・フォルダパス
 * @param[out] pos 見つかった位置
 */
inline bool FindPathBound(const string &str, string::size_type &pos)
{
	if((pos = str.find_last_of("/")) == string::npos){
		if((pos = str.find_last_of("\\")) == string::npos){
			return false;
		}
	}

	return true;
}

/*!
 * ファイル名比較関数(拡張子)
 * @param[in] fn 比較したいファイル名
 * @param[in] ext 拡張子
 * @return fnの拡張子がextと同じならtrue
 */
inline bool SearchCompExt(const string &fn, const string &ext)
{
	return (fn.find(ext, 0) != string::npos);
}


/*!
 * ファイル名生成
 * @param head : 基本ファイル名
 * @param ext  : 拡張子
 * @param n    : 連番
 * @param d    : 連番桁数
 * @return 生成したファイル名
 */
inline string CreateFileName(const string &head, const string &ext, int n, const int &d)
{
	string file_name = head;
	int dn = d-1;
	if(n > 0){
		dn = (int)(log10((float)n))+1;
	}
	else if(n == 0){
		dn = 1;
	}
	else{
		n = 0;
		dn = 1;
	}

	for(int i = 0; i < d-dn; ++i){
		file_name += "0";
	}

	file_name += RX_TO_STRING(n);
	if(!ext.empty() && ext[0] != '.') file_name += ".";
	file_name += ext;

	return file_name;
}




/*!
 * パスからファイル名のみ取り出す
 * @param[in] path パス
 * @return ファイル名
 */
inline string GetFileName(const string &path)
{
	size_t pos1;

	pos1 = path.rfind('\\');
	if(pos1 != string::npos){
		return path.substr(pos1+1, path.size()-pos1-1);
	}

	pos1 = path.rfind('/');
	if(pos1 != string::npos){
		return path.substr(pos1+1, path.size()-pos1-1);
	}

	return path;
}

/*!
 * パスから拡張子を小文字にして取り出す
 * @param[in] path ファイルパス
 * @return (小文字化した)拡張子
 */
inline string GetExtension(const string &path)
{
	string ext;
	size_t pos1 = path.rfind('.');
	if(pos1 != string::npos){
		ext = path.substr(pos1+1, path.size()-pos1);
		string::iterator itr = ext.begin();
		while(itr != ext.end()){
			*itr = tolower(*itr);
			itr++;
		}
		itr = ext.end()-1;
		while(itr != ext.begin()){	// パスの最後に\0やスペースがあったときの対策
			if(*itr == 0 || *itr == 32){
				ext.erase(itr--);
			}
			else{
				itr--;
			}
		}
	}

	return ext;
}

/*!
 * ファイルストリームを開く
 * @param[out] file ファイルストリーム
 * @param[in] path  ファイルパス
 * @param[in] rw    入出力フラグ (1:読込, 2:書込, 4:追記)
 * @return ファイルオープン成功:1, 失敗:0
 */
static inline int OpenFileStream(fstream &file, const string &path, int rw = 1)
{
	file.open(path.c_str(), (rw & 0x01 ? ios::in : 0)|(rw & 0x02 ? ios::out : 0)|(rw & 0x04 ? ios::app : 0));
	if(!file || !file.is_open() || file.bad() || file.fail()){
		return 0;
	}
	return 1;
}

/*!
 * ディレクトリ作成(多階層対応) - Windows only
 * @param[in] dir 作成ディレクトリパス
 * @return 成功で1,失敗で0 (ディレクトリがすでにある場合も1を返す)
 */
static int MkDir(string dir)
{
#ifdef WIN32
	if(_mkdir(dir.c_str()) != 0){
		char cur_dir[512];
		_getcwd(cur_dir, 512);	// カレントフォルダを確保しておく
		if(_chdir(dir.c_str()) == 0){	// chdirでフォルダ存在チェック
			cout << "MkDir : " << dir << " is already exist." << endl;
			_chdir(cur_dir);	// カレントフォルダを元に戻す
			return 1;
		}
		else{
			size_t pos = dir.find_last_of("\\/");
			if(pos != string::npos){	// 多階層の可能性有り
				int parent = MkDir(dir.substr(0, pos+1));	// 親ディレクトリを再帰的に作成
				if(parent){
					if(_mkdir(dir.c_str()) == 0){
						return 1;
					}
					else{
						return 0;
					}
				}
			}
			else{
				return 0;
			}
		}
	}
	return 1;
#else
	return 0;
#endif
}

//! ファイルパス検索用(made by 金森先生)
// [How to use]
// PathFinder p;
// p.addSearchPath("bin");
// p.addSearchPath("../bin");
// p.addSearchPath("../../bin");
// std::string filename = p.find("sample.bmp");
class PathFinder
{
public:
	void addSearchPath(const std::string& dirname) { m_SearchDirs.push_back(dirname); }

	std::string find(const std::string& filename, const std::string& separator = "/")
	{
		FILE* fp = 0;
		std::string path_name;

		for(unsigned int i = 0; i < m_SearchDirs.size(); i++)
		{
			path_name = m_SearchDirs[i] + separator + filename;
#ifdef _WIN32
			::fopen_s(&fp, path_name.c_str(), "r");
#else
			fp = fopen(path_name.c_str(), "r");
#endif

			if(fp != 0)
			{
				fclose(fp);
				return path_name;
			}
		}

		path_name = filename;
#ifdef _WIN32
		::fopen_s(&fp, path_name.c_str(), "r");
#else
		fp = fopen(path_name.c_str(), "r");
#endif

		if(fp != 0)
		{
			fclose(fp);
			return path_name;
		}

		printf("file not found: %s\n", filename.c_str());

		return "";
	}

private:
	std::vector<std::string> m_SearchDirs;
};




//-----------------------------------------------------------------------------
// グラデーション色生成
//-----------------------------------------------------------------------------
/*!
* 青->緑->赤->白と変化するグラデーション色生成
* @param[out] col 生成された色
* @param[in] x 値
* @param[in] xmin 最小値
* @param[in] xmax 最大値
*/
inline void Gradation(float col[3], float x, const float xmin = 0.0, const float xmax = 1.0)
{
	float l = xmax-xmin;
	if(fabs(l) < 1e-10) return;

	const int ncolors = 7;
	float base[ncolors][3] = { {0.0, 0.0, 0.0},
		{0.0, 0.0, 1.0},
		{0.0, 1.0, 1.0},
		{0.0, 1.0, 0.0},
		{1.0, 1.0, 0.0},
		{1.0, 0.0, 0.0},
		{1.0, 1.0, 1.0} };
	x = RX_CLAMP<float>(((x-xmin)/l), 0.0, 1.0)*(ncolors-1);
	int i = (int)x;
	float dx = x-floor(x);
	col[0] = RX_LERP(base[i][0], base[i+1][0], dx);
	col[1] = RX_LERP(base[i][1], base[i+1][1], dx);
	col[2] = RX_LERP(base[i][2], base[i+1][2], dx);
}


/*!
* グラデーション色生成
* @param[in] c1,c2 x=xmin,x=xmaxのときの色
* @param[in] x 値
* @param[in] xmin,xmax 最小値,最大値
*/
inline glm::vec3 gradation(glm::vec3 c1, glm::vec3 c2, float x, const float xmin = 0.0, const float xmax = 1.0)
{
	float l = xmax-xmin;
	if(fabs(l) < 1e-10) return c1;

	float dx = glm::clamp<float>((x-xmin)/l, 0.0f, 1.0f);
	return glm::lerp<float>(c1, c2, dx);
}



//! グラデーション色の生成
template<class T> 
inline void RX_COLOR_RAMP(T t, T *r)
{
	const int ncolors = 7;
	T c[ncolors][3] = {	{ 0.0, 0.0, 1.0 },
		{ 0.0, 0.5, 1.0 },
		{ 0.0, 1.0, 1.0 },
		{ 0.0, 1.0, 0.0 },
		{ 1.0, 1.0, 0.0 },
		{ 1.0, 0.0, 0.0 },
		{ 1.0, 0.0, 1.0 } };
	//T c[ncolors][3] = { { 1.0, 0.0, 0.0 },
	//					{ 1.0, 0.5, 0.0 },
	//					{ 1.0, 1.0, 0.0 },
	//					{ 0.0, 1.0, 0.0 },
	//					{ 0.0, 1.0, 1.0 },
	//					{ 0.0, 0.0, 1.0 },
	//					{ 1.0, 0.0, 1.0 } };
	t = t*(ncolors-1);
	int i = (int)t;
	T u = t-floor(t);
	r[0] = RX_LERP(c[i][0], c[i+1][0], u);
	r[1] = RX_LERP(c[i][1], c[i+1][1], u);
	r[2] = RX_LERP(c[i][2], c[i+1][2], u);
}




//-----------------------------------------------------------------------------
// 衝突処理/距離計算用関数
//-----------------------------------------------------------------------------
/*!
 * AABBと点の距離
 * @param[in] spos 立方体の中心を原点とした相対座標値
 * @param[in] r    半径(球の場合)
 * @param[in] sgn  立方体の内で距離が正:1,外で正:-1
 * @param[in] vMin 立方体の最小座標値(相対座標)
 * @param[in] vMax 立方体の最大座標値(相対座標)
 * @param[out] d   符号付距離値
 * @param[out] n   最近傍点の法線方向
 */
static inline bool dist_aabb_point(const glm::vec3& spos, const float& r, const int& sgn,
	const glm::vec3& vMin, const glm::vec3& vMax, float& d, glm::vec3& n)
{
	int bout = 0;
	float d0[6];
	int idx0 = -1;

	//! AABBの法線方向(内向き)
	const glm::vec3 AABB_NORMALS[6] = { glm::vec3(1.0,  0.0,  0.0),	// x-
		glm::vec3(-1.0,  0.0,  0.0),	// x+
		glm::vec3(0.0,  1.0,  0.0),	// y-
		glm::vec3(0.0, -1.0,  0.0),	// y+
		glm::vec3(0.0,  0.0,  1.0),	// z-
		glm::vec3(0.0,  0.0, -1.0) };	// z+


	// 各軸ごとに最小と最大境界外になっていないか調べる
	int c = 0;
	for(int i = 0; i < 3; ++i){
		int idx = 2*i;
		if((d0[idx] = (spos[i]-r)-vMin[i]) < 0.0){
			bout |= (1 << idx); c++;
			idx0 = idx;
		}
		idx = 2*i+1;
		if((d0[idx] = vMax[i]-(spos[i]+r)) < 0.0){
			bout |= (1 << idx); c++;
			idx0 = idx;
		}
	}

	// AABB内(全軸で境界内)
	if(bout == 0){
		float min_d = 1e10;
		int idx1 = -1;
		for(int i = 0; i < 6; ++i){
			if(d0[i] <= min_d){
				min_d = d0[i];
				idx1 = i;
			}
		}

		d = sgn*min_d;
		n = (idx1 != -1) ? float(sgn)*AABB_NORMALS[idx1] : glm::vec3(0.0);
		return true;
	}

	// AABB外
	glm::vec3 x(0.0);
	for(int i = 0; i < 3; ++i){
		if(bout & (1 << (2*i))){
			x[i] = d0[2*i];
		}
		else if(bout & (1 << (2*i+1))){
			x[i] = -d0[2*i+1];
		}
	}

	// sgn = 1:箱，-1:オブジェクト
	if(c == 1){
		// 平面近傍
		d = sgn*d0[idx0];
		n = float(sgn)*AABB_NORMALS[idx0];
	}
	else{
		// エッジ/コーナー近傍
		d = -sgn*glm::length(x);
		n = float(sgn)*(-glm::normalize(x));
	}

	return false;
}


/*!
 * 平面の陰関数値計算
 * @param[in] pos 陰関数値を計算する位置
 * @param[out] d,n 平面までの距離(陰関数値)と法線
 * @param[out] v 平面の速度(今のところ常に0)
 * @param[in] pn,pq 平面の法線と平面上の点
 */
inline bool dist_plane_point(const glm::vec3& pos, float& d, glm::vec3& n, glm::vec3& v, const glm::vec3& pn, const glm::vec3& pq)
{
	d = glm::dot(pq-pos, pn);
	n = pn;
	v = glm::vec3(0.0);
	return true;
}

/*!
* 線分(を含む直線)と点の距離
* @param[in] v0,v1 線分の両端点座標
* @param[in] p 点の座標
* @return 距離
*/
inline float dist_segment_point(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& p)
{
	glm::vec3 v = glm::normalize(v1-v0);
	glm::vec3 vp = p-v0;
	glm::vec3 vh = glm::dot(vp, v)*v;
	return glm::length(vp-vh);
}



/*!
 * 上が空いたボックス形状の陰関数値計算
 *  - 形状の中心を原点とした座標系での計算
 * @param[in] spos 立方体の中心を原点とした相対座標値
 * @param[in] r    半径(球の場合)
 * @param[in] sgn  立方体の内で距離が正:1,外で正:-1
 * @param[in] sl0  ボックスの内側サイズ(side length)
 * @param[in] sl1  ボックスの外側サイズ(side length)
 * @param[out] d   符号付距離値
 * @param[out] n   最近傍点の法線方向
 */
static inline void dist_openbox_point(const glm::vec3& spos, const float& r, const int& sgn,
	const glm::vec3& sl0, const glm::vec3& sl1, float& d, glm::vec3& n)
{
	int t = 2;
	d = RX_FEQ_INF;

	float td;
	glm::vec3 tn;

	// 底
	glm::vec3 m0, m1;
	m0 = -sl1;
	m1 = sl1;
	m1[t] = sl0[t];
	dist_aabb_point(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// 側面 -x
	m0 = glm::vec3(-sl1[0], -sl1[1], -sl0[2]);
	m1 = glm::vec3(-sl0[0], sl1[1], sl1[2]);
	dist_aabb_point(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// 側面 +x
	m0 = glm::vec3(sl0[0], -sl1[1], -sl0[2]);
	m1 = glm::vec3(sl1[0], sl1[1], sl1[2]);
	dist_aabb_point(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// 側面 -y
	m0 = glm::vec3(-sl0[0], -sl1[1], -sl0[2]);
	m1 = glm::vec3(sl0[0], -sl0[1], sl1[2]);
	dist_aabb_point(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}

	// 側面 +y
	m0 = glm::vec3(-sl0[0], sl0[1], -sl0[2]);
	m1 = glm::vec3(sl0[0], sl1[1], sl1[2]);
	dist_aabb_point(spos, r, sgn, m0, m1, td, tn);
	if(td < d){
		d = td;
		n = tn;
	}
}


/*!
 * 光線(レイ,半直線)と球の交差判定
 * @param[in] p,d レイの原点と方向
 * @param[in] c,r 球の中心と半径
 * @param[out] t1,t2 pから交点までの距離
 * @return 交点数
 */
inline int ray_sphere(const glm::vec3& p, const glm::vec3& d, const glm::vec3& sc, const float r, float& t1, float& t2)
{
	glm::vec3 q = p-sc;	// 球中心座標系での光線原点座標

	float a = glm::length(d);
	float b = 2*glm::dot(q, d);
	float c = glm::length2(q)-r*r;

	// 判別式
	float D = b*b-4*a*c;

	if(D < 0.0){ // 交差なし
		return 0;
	}
	else if(D < RX_FEQ_EPS){ // 交点数1
		t1 = -b/(2*a);
		t2 = -1;
		return 1;
	}
	else{ // 交点数2
		float sqrtD = sqrt(D);
		t1 = (-b-sqrtD)/(2*a);
		t2 = (-b+sqrtD)/(2*a);
		return 2;
	}

}

/*!
 * 三角形と球の交差判定
 * @param[in] v0,v1,v2	三角形の頂点
 * @param[in] n			三角形の法線
 * @param[in] p			最近傍点
 * @return
 */
inline bool triangle_sphere(const glm::vec3& v0, const glm::vec3& v1, const glm::vec3& v2, const glm::vec3& n,
	const glm::vec3& c, const float& r, float& dist, glm::vec3& ipoint)
{
	// ポリゴンを含む平面と球中心の距離
	float d = glm::dot(v0, n);
	float l = glm::dot(n, c)-d;

	dist = l;
	if(l > r) return false;

	// 平面との最近傍点座標
	glm::vec3 p = c-l*n;

	// 近傍点が三角形内かどうかの判定
	glm::vec3 n1 = glm::cross((v0-c), (v1-c));
	glm::vec3 n2 = glm::cross((v1-c), (v2-c));
	glm::vec3 n3 = glm::cross((v2-c), (v0-c));

	ipoint = p;
	dist = l;
	if(glm::dot(n1, n2) > 0 && glm::dot(n2, n3) > 0){		// 三角形内
		return true;
	}
	else{		// 三角形外
	 // 三角形の各エッジと球の衝突判定
		for(int e = 0; e < 3; ++e){
			glm::vec3 va0 = (e == 0 ? v0 : (e == 1 ? v1 : v2));
			glm::vec3 va1 = (e == 0 ? v1 : (e == 1 ? v2 : v0));

			float t1, t2;
			int n = ray_sphere(va0, glm::normalize(va1-va0), c, r, t1, t2);

			if(n){
				float le2 = glm::length2(va1-va0);
				if((t1 >= 0.0 && t1*t1 < le2) || (t2 >= 0.0 && t2*t2 < le2)){
					return true;
				}
			}
		}
		return false;
	}
}

/*!
 * 線分と球の交差判定
 * @param[in] s0,s1	線分の端点
 * @param[in] sc,r   球の中心座標と半径
 * @param[out] d2 線分との距離の二乗
 * @return 交差ありでtrue
 */
inline bool segment_sphere(const glm::vec3& s0, const glm::vec3& s1, const glm::vec3& sc, const float& r, float& d2)
{
	glm::vec3 v = s1-s0;
	glm::vec3 c = sc-s0;

	float vc = glm::dot(v, c);
	if(vc < 0){		// 球の中心が線分の始点s0の外にある
		d2 = glm::length2(c);
		return (d2 < r* r);	// 球中心と始点s0の距離で交差判定
	}
	else{
		float v2 = glm::length2(v);
		if(vc > v2){	// 球の中心が線分の終点s1の外にある
			d2 = glm::length2(s1-sc);
			return (d2 < r* r);	// 球中心と終点s1の距離で交差判定
		}
		else{			// 球がs0とs1の間にある
			d2 = glm::length2((vc*v)/glm::length2(v)-c);
			return (d2 < r* r);	// 直線と球中心の距離で交差判定
		}
	}

	return false;
}

/*!
* レイ/線分と三角形の交差
* @param[in] P0,P1 レイ/線分の端点orレイ上の点
* @param[in] V0,V1,V2 三角形の頂点座標
* @param[out] I 交点座標
* @retval 1 交点Iで交差 
* @retval 0 交点なし
* @retval 2 三角形の平面内
* @retval -1 三角形が"degenerate"である(面積が0，つまり，線分か点になっている)
*/
static inline int segment_triangle(glm::vec3 P0, glm::vec3 P1, glm::vec3 V0, glm::vec3 V1, glm::vec3 V2,
								   glm::vec3 &I, glm::vec3 &n, float rp)
{
	// 三角形のエッジベクトルと法線
	glm::vec3 u = V1-V0;
	glm::vec3 v = V2-V0;
	n = glm::normalize(glm::cross(u, v));
	if(glm::length(n) < 1e-6){
		return -1;	// 三角形が"degenerate"である(面積が0)
	}

	// 線分
	glm::vec3 dir = P1-P0;
	float a = glm::dot(n, P0-V0);
	float b = glm::dot(n, dir);
	if(fabs(b) < 1e-10){	// 線分と三角形平面が平行
		if(a == 0){
			return 2;	// 線分が平面上
		}
		else{
			return 0;	// 交点なし
		}
	}

	// 交点計算

	// 2端点がそれぞれ異なる面にあるかどうかを判定
	float r = -a/b;
	glm::vec3 offset = glm::vec3(0.0);
	float dn = 0;
	float sign_n = 1;
	if(a < 0){
		return 0;
	}

	if(r < 0.0){
		return 0;
	}
	else{
		if(fabs(a) > fabs(b)){
			return 0;
		}
		else{
			if(b > 0){
				return 0;
			}
		}
	}

	// 線分と平面の交点
	I = P0+r*dir;

	// 交点が三角形内にあるかどうかの判定
	float uu, uv, vv, wu, wv, D;
	uu = glm::dot(u, u);
	uv = glm::dot(u, v);
	vv = glm::dot(v, v);
	glm::vec3 w = I-V0;
	wu = glm::dot(w, u);
	wv = glm::dot(w, v);
	D = uv*uv-uu*vv;

	float s, t;
	s = (uv*wv-vv*wu)/D;
	if(s < 0.0 || s > 1.0){
		return 0;
	}

	t = (uv*wu-uu*wv)/D;
	if(t < 0.0 || (s+t) > 1.0){
		return 0;
	}

	return 1;
}

#endif // #ifndef _UTILS_H_