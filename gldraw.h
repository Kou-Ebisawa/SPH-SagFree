/*!
  @file rx_gldraw.h

  @brief OpenGLによる描画関数

  @author Makoto Fujisawa
  @date 2020-07
*/
// FILE --rx_gldraw.h--

#ifndef _RX_GLDRAW_H_
#define _RX_GLDRAW_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <cmath>

// ユーティリティ
#include "utils.h"

// 画像読み書き
#include "rx_bitmap.h"

// 文字列描画
#include "FTGL/ftgl.h"


using namespace std;




//! OpenGLのエラーチェック
// see https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/glGetError.xhtml
static inline void CheckGLError(const char *func, const char *file, int line)
{
	GLenum errCode = glGetError();

	if(errCode != GL_NO_ERROR)
	{
		std::cerr << func << ": OpenGL Error: ";

		switch(errCode)
		{
		case GL_INVALID_ENUM:
			std::cerr << "GL_INVALID_ENUM: An unacceptable value is specified for an enumerated argument";
			break;
		case GL_INVALID_VALUE:
			std::cerr << "GL_INVALID_VALUE: A numeric argument is out of range";
			break;
		case GL_INVALID_OPERATION:
			std::cerr << "GL_INVALID_OPERATION: The specified operation is not allowed in the current state";
			break;
		case GL_STACK_OVERFLOW:
			std::cerr << "GL_STACK_OVERFLOW: This command would cause a stack overflow";
			break;
		case GL_STACK_UNDERFLOW:
			std::cerr << "GL_STACK_UNDERFLOW: This command would cause a stack underflow";
			break;
		case GL_OUT_OF_MEMORY:
			std::cerr << "GL_OUT_OF_MEMORY: There is not enough memory left to execute the command";
			break;
		case GL_TABLE_TOO_LARGE:
			std::cerr << "GL_TABLE_TOO_LARGE: The specified table exceeds the implementation's maximum supported table size";
			break;
		default:
			std::cerr << "Unknown command, error code = " << std::showbase << std::hex << errCode;
			break;
		}

		std::cerr << " (file: " << file << " at line " << line << std::endl;
	}
}
#define CHECK_GL_ERROR CheckGLError(__FUNCTION__, __FILE__, __LINE__)



//-----------------------------------------------------------------------------
// OpenGLの変換関係の関数
//-----------------------------------------------------------------------------
/*!
* double版glMaterial
* @param[in] face 材質を指定する面(GL_FRONT,GL_BACK,GL_FRONT_AND_BACK)
* @param[in] pname 指定する材質の種類(GL_AMBIENT,GL_DIFFUSE,GL_SPECULAR,GL_EMISSION,GL_SHININESS,GL_AMBIENT_AND_DIFFUSE,GL_COLOR_INDEXES)
* @param[in] params 設定する材質の値
*/
inline void glMaterialdv(GLenum face, GLenum pname, const GLdouble* params)
{
	GLfloat col[4];
	col[0] = (GLfloat)params[0];
	col[1] = (GLfloat)params[1];
	col[2] = (GLfloat)params[2];
	col[3] = (GLfloat)params[3];
	glMaterialfv(face, pname, col);
}

/*!
* double版glLight
* @param[in] light 光源番号(GL_LIGHT0 ～ GL_LIGHT7)
* @param[in] pname 指定する光源パラメータの種類(GL_POSITION,GL_DIFFUSE,GL_AMBIENT,GL_SPECULAR,GL_LINEAR_ATTENUATION,GL_SPOT_DIRECTION,GL_SPOT_CUTOFF,GL_SPOT_EXPONENT)
* @param[in] params 設定する値
*/
inline void glLightdv(GLenum light, GLenum pname, const GLdouble* params)
{
	GLfloat fparams[4];
	fparams[0] = (GLfloat)params[0];
	fparams[1] = (GLfloat)params[1];
	fparams[2] = (GLfloat)params[2];
	fparams[3] = (GLfloat)params[3];
	glLightfv(light, pname, fparams);
}
inline void glLightdv3(GLenum light, GLenum pname, const GLdouble* params)
{
	GLfloat fparams[4];
	fparams[0] = (GLfloat)params[0];
	fparams[1] = (GLfloat)params[1];
	fparams[2] = (GLfloat)params[2];
	fparams[3] = 1.0f;
	glLightfv(light, pname, fparams);
}


/*!
* 3x3行列からOpenGLの変換行列を取得
* @param[in] mat 3x3行列を格納した1次元配列(mat[i*3+j]=m[i][j])
* @param[out] glmat OpenGL変換行列(4x4)を格納した1次元配列
*/
inline void GetGLMatrix(float mat[9], GLfloat glmat[16])
{
	glmat[0]  = mat[0];
	glmat[1]  = mat[3]; 
	glmat[2]  = mat[6]; 
	glmat[3]  = 0.0f; 
	glmat[4]  = mat[1]; 
	glmat[5]  = mat[4]; 
	glmat[6]  = mat[7]; 
	glmat[7]  = 0.0f; 
	glmat[8]  = mat[2]; 
	glmat[9]  = mat[5]; 
	glmat[10] = mat[8]; 
	glmat[11] = 0.0f; 
	glmat[12] = 0.0f; 
	glmat[13] = 0.0f; 
	glmat[14] = 0.0f; 
	glmat[15] = 1.0f;
}



//-----------------------------------------------------------------------------
// テクスチャ関連
//-----------------------------------------------------------------------------
/*!
 * 画像読込み -> OpenGLテクスチャ登録
 * @param[in] fn ファイル名
 * @param[inout] tex_name テクスチャ名(0なら新たに生成)
 * @param[in] mipmap ミップマップ使用フラグ
 * @param[in] compress テクスチャ圧縮使用フラグ
 */
static int loadTexture(const string& fn, GLuint& tex_name, bool mipmap, bool compress)
{
	// 画像読み込み
	int w, h, c, wstep;
	unsigned char* pimg;
	pimg = ReadBitmapFile(fn, w, h, c, wstep, true, false);
	if(!pimg){
		return 0;
	}

	GLuint iformat, format;

	// 画像フォーマット
	format = GL_RGBA;
	if(c == 1){
		format = GL_LUMINANCE;
	} else if(c == 3){
		format = GL_RGB;
	}

	// OpenGL内部の格納フォーマット
	if(compress){
		iformat = GL_COMPRESSED_RGBA_S3TC_DXT1_EXT;
		if(c == 1){
			iformat = GL_COMPRESSED_LUMINANCE_ARB;
		} else if(c == 3){
			iformat = GL_COMPRESSED_RGB_S3TC_DXT1_EXT;
		}
	} else{
		iformat = GL_RGBA;
		if(c == 1){
			iformat = GL_LUMINANCE;
		} else if(c == 3){
			iformat = GL_RGB;
		}
	}

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	// テクスチャ作成
	if(tex_name == 0){
		glGenTextures(1, &tex_name);

		glBindTexture(GL_TEXTURE_2D, tex_name);

		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (mipmap ? GL_LINEAR_MIPMAP_LINEAR : GL_LINEAR));
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

		if(mipmap){
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 6);
		}

		glTexImage2D(GL_TEXTURE_2D, 0, iformat, w, h, 0, format, GL_UNSIGNED_BYTE, pimg);

		if(mipmap){
			glGenerateMipmapEXT(GL_TEXTURE_2D);
		}
	} else{
		glBindTexture(GL_TEXTURE_2D, tex_name);
		//glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, format, GL_UNSIGNED_BYTE, pimg);
		glTexImage2D(GL_TEXTURE_2D, 0, iformat, w, h, 0, format, GL_UNSIGNED_BYTE, pimg);

		if(mipmap){
			glGenerateMipmapEXT(GL_TEXTURE_2D);
		}
	}

	glBindTexture(GL_TEXTURE_2D, 0);

	delete[] pimg;

	return 1;
}

/*!
 * 画像読込み -> OpenGLテクスチャ登録
 * @param[in] fn ファイル名
 * @param[inout] tex_name テクスチャ名(0なら新たに生成)
 * @param[in] mipmap ミップマップ使用フラグ
 * @param[in] compress テクスチャ圧縮使用フラグ
 */
static int makeCheckerBoardTexture(GLuint& tex_name, const int size = 72, bool mipmap = false, bool compress = false)
{
	// ピクセルデータ生成
	unsigned char* pimg = new unsigned char[size*size*4];
	if(!pimg) return 0;
	for(int i = 0; i < size; ++i){
		for(int j = 0; j < size; ++j){
			if((i+j%2)%2 == 0){
				// dark
				pimg[i*size*4 + j*4 + 0] = (GLubyte)195;
				pimg[i*size*4 + j*4 + 1] = (GLubyte)188;
				pimg[i*size*4 + j*4 + 2] = (GLubyte)207;
				pimg[i*size*4 + j*4 + 3] = (GLubyte)255;
			} else {
				// light
				pimg[i*size*4 + j*4 + 0] = (GLubyte)220;
				pimg[i*size*4 + j*4 + 1] = (GLubyte)213;
				pimg[i*size*4 + j*4 + 2] = (GLubyte)232;
				pimg[i*size*4 + j*4 + 3] = (GLubyte)255;
			}
		}
	}


	GLuint iformat, format;

	// 画像フォーマット
	format = GL_RGBA;

	// OpenGL内部の格納フォーマット
	if(compress){
		iformat = GL_COMPRESSED_RGBA_S3TC_DXT1_EXT;
	} else{
		iformat = GL_RGBA;
	}

	glPixelStorei(GL_UNPACK_ALIGNMENT, 4);

	// テクスチャ作成
	if(tex_name == 0){
		glGenTextures(1, &tex_name);
		glBindTexture(GL_TEXTURE_2D, tex_name);

		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, (mipmap ? GL_NEAREST_MIPMAP_NEAREST : GL_NEAREST));

		if(mipmap){
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_BASE_LEVEL, 0);
			glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 6);
		}

		glTexImage2D(GL_TEXTURE_2D, 0, iformat, size, size, 0, format, GL_UNSIGNED_BYTE, pimg);

		if(mipmap){
			glGenerateMipmapEXT(GL_TEXTURE_2D);
		}
	} else{
		glBindTexture(GL_TEXTURE_2D, tex_name);
		//glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, size, size, format, GL_UNSIGNED_BYTE, pimg);
		glTexImage2D(GL_TEXTURE_2D, 0, iformat, size, size, 0, format, GL_UNSIGNED_BYTE, pimg);

		if(mipmap){
			glGenerateMipmapEXT(GL_TEXTURE_2D);
		}
	}

	glBindTexture(GL_TEXTURE_2D, 0);

	delete[] pimg;

	return 1;
}

/*!
 * フレームバッファのRGB情報を一時的なバッファに確保
 * @retval true 保存成功
 * @retval false 保存失敗
 */
static inline bool saveFrameBuffer(string fn, int w, int h)
{
	string ext = GetExtension(fn);
	if(ext != "bmp") fn += ".bmp";
		
	int c = 3;
	int wstep = (((w+1)*c)/4)*4;
	vector<unsigned char> imm_buf(wstep*h);

	glReadPixels(0, 0, w, h, GL_RGB, GL_UNSIGNED_BYTE, &imm_buf[0]);
	WriteBitmapFile(fn, &imm_buf[0], w, h, c, RX_BMP_WINDOWS_V3, wstep, false, true);
	return true;
}


//-----------------------------------------------------------------------------
// VBO関連
//-----------------------------------------------------------------------------
struct MeshVBO
{
	GLuint vrts, vrts_attrib;
	GLuint tcds, tcds_attrib;
	GLuint nrms, nrms_attrib;
	GLuint cols, cols_attrib;
	GLuint tris;
	int nvrts, ntris;
	MeshVBO() :vrts(0),tcds(0),nrms(0),cols(0),tris(0),nvrts(0),ntris(0),vrts_attrib(0),tcds_attrib(1),nrms_attrib(2),cols_attrib(3) {}
	~MeshVBO()
	{
		glDeleteBuffers(1, &vrts);
		glDeleteBuffers(1, &tcds);
		glDeleteBuffers(1, &nrms);
		glDeleteBuffers(1, &cols);
		glDeleteBuffers(1, &tris);
	}
};

/*!
* VBOの作成
* - 最低限頂点情報があれば良い．必要のない引数は0を入れておく．
* @param[in] vrts,nvrts 頂点座標配列と頂点数
* @param[in] dim 次元(2or3)
* @param[in] tris,ntris ポリゴンを構成する頂点インデックス列を格納した配列とポリゴン数
* @param[in] nelem 1ポリゴンの頂点数(GL_TRIANGLESなら3,GL_QUADSなら4)
* @param[in] nrms,nnrms 各頂点の法線情報配列と法線数(=頂点数)
* @param[in] cols,ncols 各頂点の色情報配列と色情報数(=頂点数)
* @param[in] tcds,ntcds 各頂点のテクスチャ座標情報配列とテクスチャ座標数(=頂点数)
* @return VBO
*/
static inline int CreateVBO(MeshVBO &vbo, 
	const GLfloat* vrts, GLuint nvrts, int dim = 3,
	const int* tris = 0, GLuint ntris = 0, GLuint nelem = 3, 
	const GLfloat* nrms = 0, GLuint nnrms = 0,
	const GLfloat* cols = 0, GLuint ncols = 0,
	const GLfloat* tcds = 0, GLuint ntcds = 0)
{
	if(!nvrts) return 0;

	// VBO:頂点座標
	glGenBuffers(1, &vbo.vrts);
	glBindBuffer(GL_ARRAY_BUFFER, vbo.vrts);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * dim * nvrts, vrts, GL_STATIC_DRAW);
	vbo.nvrts = nvrts;

	// VBO:頂点テクスチャ座標
	if (tcds) {
		glGenBuffers(1, &vbo.tcds);
		glBindBuffer(GL_ARRAY_BUFFER, vbo.tcds);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 2 * ntcds, tcds, GL_STATIC_DRAW);
	}

	// VBO:頂点法線
	if (nrms) {
		glGenBuffers(1, &vbo.nrms);
		glBindBuffer(GL_ARRAY_BUFFER, vbo.nrms);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * dim * nnrms, nrms, GL_STATIC_DRAW);
	}

	// VBO:頂点描画色
	if (cols) {
		glGenBuffers(1, &vbo.cols);
		glBindBuffer(GL_ARRAY_BUFFER, vbo.cols);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * 3 * ncols, cols, GL_STATIC_DRAW);
	}

	// VBO:三角形ポリゴンインデックス
	if (tris) {
		glGenBuffers(1, &vbo.tris);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo.tris);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int) * ntris * nelem, tris, GL_STATIC_DRAW);
		vbo.ntris = ntris;
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	return 1;
}
/*!
* VBOの破棄
*/
static inline void DeleteVBO(MeshVBO vbo)
{
	if(vbo.vrts == 0) return;
	glDeleteBuffers(1, &vbo.vrts); vbo.vrts = 0;

	if(vbo.tcds != 0){ glDeleteBuffers(1, &vbo.tcds); vbo.tcds = 0; }
	if(vbo.nrms != 0){ glDeleteBuffers(1, &vbo.nrms); vbo.nrms = 0; }
	if(vbo.cols != 0){ glDeleteBuffers(1, &vbo.tcds); vbo.cols = 0; }
	if(vbo.tris != 0){ glDeleteBuffers(1, &vbo.tris); vbo.tris = 0; }

}

/*!
* VAOによるメッシュ描画
* @param[in] vao VAOオブジェクト
* @param[in] draw 描画フラグ
* @param[in] nvrts,ntris 頂点数,三角形ポリゴン数
* @param[in] col[3] 頂点,エッジ,面の色
*/
static inline void DrawMeshVBO(MeshVBO vbo, int draw, glm::vec3 col[3])
{
	if (!vbo.vrts) return;

	glBindBuffer(GL_ARRAY_BUFFER, vbo.vrts);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.nrms);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, 0);


	if (draw & 0x04) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glEnable(GL_LIGHTING);

		glColor3fv(glm::value_ptr(col[2]));
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo.tris);
		glDrawElements(GL_TRIANGLES, vbo.ntris * 3, GL_UNSIGNED_INT, 0);
	}
	if (draw & 0x02) {
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDisable(GL_LIGHTING);

		glLineWidth(1.0);
		glColor3fv(glm::value_ptr(col[1]));
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo.tris);
		glDrawElements(GL_TRIANGLES, vbo.ntris * 3, GL_UNSIGNED_INT, 0);
	}
	if (draw & 0x01) {
		glDisable(GL_LIGHTING);

		glPointSize(8.0);
		glColor3fv(glm::value_ptr(col[0]));
		glBindBuffer(GL_ARRAY_BUFFER, vbo.vrts);
		glDrawArrays(GL_POINTS, 0, vbo.nvrts);
	}
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}


/*!
* 立方体形状のポリゴンメッシュを生成(VBO用)
* - 立方体の中心は原点(0,0,0)
* - 四角形メッシュバージョン
* @param[out] nvrts,ntris 生成したメッシュの頂点数とポリゴン数を返す
* @param[in] len 辺の長さ
*/
static inline int MakeCube(int &nvrts, vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, int &ntris, vector<int> &tris,
	float len)
{
	nvrts = 0; ntris = 0;
	vrts.clear(); nrms.clear(); tris.clear();

	// 頂点座標
	float sl = 0.5f*len;	// 辺の長さの半分
	vrts.push_back(glm::vec3(-sl, -sl, -sl));
	vrts.push_back(glm::vec3( sl, -sl, -sl));
	vrts.push_back(glm::vec3( sl, -sl,  sl));
	vrts.push_back(glm::vec3(-sl, -sl,  sl));
	vrts.push_back(glm::vec3(-sl,  sl, -sl));
	vrts.push_back(glm::vec3( sl,  sl, -sl));
	vrts.push_back(glm::vec3( sl,  sl,  sl));
	vrts.push_back(glm::vec3(-sl,  sl,  sl));
	nvrts = static_cast<int>(vrts.size());

	// 頂点法線
	for(const glm::vec3 &v : vrts){
		nrms.push_back(glm::normalize(v));
	}

	// 四角形メッシュ作成
	int f[6][4] = { {0,1,2,3}, {4,7,6,5}, {0,3,7,4}, {1,5,6,2}, {0,4,5,1}, {3,2,6,7} };
	for(int j = 0; j < 6; ++j){
		for(int i = 0; i < 4; ++i){
			tris.push_back(f[j][i]);
		}
	}
	ntris = static_cast<int>(tris.size()/4);

	return 1;
}

/*!
* 立方体形状のポリゴンメッシュを生成(VBO用)
* - 立方体の中心は原点(0,0,0)
* - 面毎に頂点を分けて面法線を適用するバージョン
* @param[out] nvrts,ntris 生成したメッシュの頂点数とポリゴン数を返す
* @param[in] len 辺の長さ
*/
static inline int MakeCubeWithFaceNormal(int &nvrts, vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, int &ntris, vector<int> &tris,
	float len)
{
	nvrts = 0; ntris = 0;
	vrts.clear(); nrms.clear(); tris.clear();

	nvrts = 6*4;
	ntris = 6;
	int f[6][4] = { {0,1,2,3}, {4,7,6,5}, {0,3,7,4}, {1,5,6,2}, {0,4,5,1}, {3,2,6,7} };

	// 頂点座標
	float sl = 0.5f*len;	// 辺の長さの半分
	vector<glm::vec3> v0(8);
	v0[0] = glm::vec3(-sl, -sl, -sl);
	v0[1] = glm::vec3( sl, -sl, -sl);
	v0[2] = glm::vec3( sl, -sl,  sl);
	v0[3] = glm::vec3(-sl, -sl,  sl);
	v0[4] = glm::vec3(-sl,  sl, -sl);
	v0[5] = glm::vec3( sl,  sl, -sl);
	v0[6] = glm::vec3( sl,  sl,  sl);
	v0[7] = glm::vec3(-sl,  sl,  sl);

	// 頂点法線
	vector<glm::vec3> n0(6);
	n0[0] = glm::vec3(0, -1, 0);
	n0[1] = glm::vec3(0,  1, 0);
	n0[2] = glm::vec3(-1, 0, 0);
	n0[3] = glm::vec3( 1, 0, 0);
	n0[4] = glm::vec3(0, 0, -1);
	n0[5] = glm::vec3(0, 0,  1);

	vrts.resize(nvrts);
	nrms.resize(nvrts);
	tris.resize(ntris*4);
	for(int j = 0; j < 6; ++j){
		for(int i = 0; i < 4; ++i){
			int idx = i+4*j;
			vrts[idx] = (v0[f[j][i]]);
			nrms[idx] = n0[j];
			tris[idx] = idx;
		}
	}

	return 1;
}

/*!
* y軸上方向を法線とした平面描画のためのポリゴンメッシュを生成(VBO用)
* - 平面は原点(0,0,0)を通り，法線が(0,1,0)
* @param[out] nvrts,ntris 生成したメッシュの頂点数とポリゴン数を返す
* @param[in] len 辺の長さ
*/
static inline int MakePlaneY(int &nvrts, vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, int &ntris, vector<int> &tris,
	vector<glm::vec2> &texcoords, float len)
{
	nvrts = 0; ntris = 0;
	vrts.clear(); nrms.clear(); tris.clear();

	nvrts = 4;
	ntris = 1;

	// 頂点座標
	float sl = 0.5f*len;	// 辺の長さの半分
	vrts.resize(4);
	vrts[0] = glm::vec3(-sl, 0, -sl);
	vrts[1] = glm::vec3(-sl, 0,  sl);
	vrts[2] = glm::vec3( sl, 0,  sl);
	vrts[3] = glm::vec3( sl, 0, -sl);

	// 頂点法線
	nrms.resize(4);
	nrms[0] = glm::vec3(0, 1, 0);
	nrms[1] = glm::vec3(0, 1, 0);
	nrms[2] = glm::vec3(0, 1, 0);
	nrms[3] = glm::vec3(0, 1, 0);

	// テクスチャ座標
	texcoords.resize(4);
	texcoords[0] = glm::vec2(0, 0);
	texcoords[1] = glm::vec2(0, 1);
	texcoords[2] = glm::vec2(1, 1);
	texcoords[3] = glm::vec2(1, 0);

	// メッシュ頂点インデックス
	tris.resize(4);
	for(int i = 0; i < 4; ++i) tris[i] = i;

	return 1;
}

/*!
* 球体形状のポリゴンメッシュを生成(VBO用)
* - 球の中心は原点(0,0,0)
* @param[out] nvrts,ntris 生成したメッシュの頂点数とポリゴン数を返す
* @param[in] rad 球の半径
* @param[in] slices,stacks 緯度方向(360度)と傾度方向(180度)のポリゴン分割数
* @return 生成したVAOオブジェクト
*/
static inline int MakeSphere(int &nvrts, vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, int &ntris, vector<int> &tris,
	float rad, int slices = 16, int stacks = 8)
{
	const float pi = glm::pi<float>();
	nvrts = 0;
	ntris = 0;
	vrts.clear();
	nrms.clear();
	tris.clear();

	for(int j = 0; j <= stacks; ++j){
		float t = float(j)/float(stacks);
		float y = rad*cos(pi*t);
		float rj = rad*sin(pi*t);	// 高さyでの球の断面円半径
		for(int i = 0; i <= slices; ++i){
			float s = float(i)/float(slices);
			float x = rj*sin(2*pi*s);
			float z = rj*cos(2*pi*s);
			vrts.push_back(glm::vec3(x, y, z));
			nrms.push_back(glm::normalize(glm::vec3(x, y, z)));
		}
	}
	nvrts = static_cast<int>(vrts.size());
	// メッシュ作成
	int nx = slices+1;
	for(int j = 0; j < stacks; ++j){
		for(int i = 0; i < slices; ++i){
			tris.push_back((i)+(j)*nx);
			tris.push_back((i+1)+(j+1)*nx);
			tris.push_back((i+1)+(j)*nx);

			tris.push_back((i)+(j)*nx);
			tris.push_back((i)+(j+1)*nx);
			tris.push_back((i+1)+(j+1)*nx);
		}
	}
	ntris = static_cast<int>(tris.size()/3);

	return 1;
}

/*!
* 円筒形状のポリゴンメッシュを生成(VBO用)
* - 円筒の中心は原点(0,0,0)
* - 円筒の軸方向はz軸方向(0,0,1) - gluCylinderに合わせている
* - 法線を別にするために側面と端面の頂点を別にしている
* @param[out] nvrts,ntris 生成したメッシュの頂点数とポリゴン数を返す
* @param[in] rad,len 円筒の半径と長さ
* @param[in] slices 円筒の円に沿ったポリゴン分割数
*/
static inline int MakeCylinder(int &nvrts, vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, int &ntris, vector<int> &tris,
	float rad1, float rad2, float len, int slices = 16, bool disk = true)
{
	const float pi = glm::pi<float>();

	// 側面用頂点
	for(int i = 0; i <= slices; ++i){
		float t = float(i)/float(slices);
		float x1 = rad1*cos(2*pi*t);
		float y1 = rad1*sin(2*pi*t);
		float x2 = rad2*cos(2*pi*t);
		float y2 = rad2*sin(2*pi*t);
		vrts.push_back(glm::vec3(x1, y1, -0.5*len));
		vrts.push_back(glm::vec3(x2, y2, 0.5*len));
		glm::vec3 n = glm::vec3(x1, y1, 0.0);
		if(rad1 < 1e-6) n = glm::vec3(x2, y2, 0.0);
		if(rad1 > rad2){
			float l = glm::length(n);
			n += l*(rad1-rad2)/len*glm::vec3(0, 0, 1);
		}
		else if(rad1 < rad2){
			float l = glm::length(n);
			n += l*(rad2-rad1)/len*glm::vec3(0, 0, -1);
		}
		if(glm::length2(n) > 1e-6) n = glm::normalize(n);
		nrms.push_back(n); nrms.push_back(n);
	}
	// 端面用頂点(座標値は↑と同じだが法線が異なる)
	int voffset = vrts.size();
	for(int i = 0; i <= slices; ++i){
		vrts.push_back(vrts[2*i]);
		nrms.push_back(glm::vec3(0.0, 0.0, -1.0));
		vrts.push_back(vrts[2*i+1]);
		nrms.push_back(glm::vec3(0.0, 0.0, 1.0));
	}

	// メッシュ作成
	for(int i = 0; i < 2*slices; i += 2){
		tris.push_back(i);
		tris.push_back((i+2 >= 2*slices ? 0 : i+2));
		tris.push_back(i+1);

		tris.push_back(i+1);
		tris.push_back((i+2 >= 2*slices ? 0 : i+2));
		tris.push_back((i+2 >= 2*slices ? 1 : i+3));
	}

	// 両端面にポリゴンを貼る
	if(disk){
		// 端面中心頂点
		vrts.push_back(glm::vec3(0, 0, -0.5*len));
		nrms.push_back(glm::normalize(glm::vec3(0, 0, -1)));
		int c1 = vrts.size()-1;
		vrts.push_back(glm::vec3(0, 0,  0.5*len));
		nrms.push_back(glm::normalize(glm::vec3(0, 0,  1)));
		int c2 = vrts.size()-1;

		for(int i = 0; i < slices; ++i){
			tris.push_back(c1);
			tris.push_back(2*i+voffset);
			tris.push_back((i == slices-1 ? 0 : 2*i+2)+voffset);

			tris.push_back(c2);
			tris.push_back(2*i+1+voffset);
			tris.push_back((i == slices-1 ? 1 : 2*i+3)+voffset);
		}
	}

	nvrts = static_cast<int>(vrts.size());
	ntris = static_cast<int>(tris.size()/3);

	return 1;
}

/*!
* カプセル形状(円筒の両端が半球)のポリゴンメッシュを生成(VBO用)
* - 形状の中心は原点(0,0,0)
* - カプセル形状の軸方向はz軸方向(0,0,1) - gluCylinderに合わせている
* @param[out] nvrts,ntris 生成したメッシュの頂点数とポリゴン数を返す
* @param[in] rad,len 円筒の半径と長さ
* @param[in] slices 円筒の円に沿ったポリゴン分割数
* @return 生成したVAOオブジェクト
*/
static inline int MakeCapsule(int &nvrts, vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, int &ntris, vector<int> &tris,
	float rad, float len, int slices = 16, int stacks = 8)
{
	const float pi = glm::pi<float>();
	int offset = 0;

	// 球体の中心断面(赤道部分)に沿った頂点を2重にして，
	// その部分を円筒長さ分z方向に伸ばすことでカプセル形状を生成
	for(int j = 0; j <= stacks; ++j){
		float t = float(j)/float(stacks);
		float z = rad*cos(pi*t);
		float rj = rad*sin(pi*t);	// 高さyでの球の断面円半径
		float zlen = (j < stacks/2 ? 0.5*len : -0.5*len);	// z方向のオフセット量

		if(j == stacks/2){	// 中心断面(赤道部分)に頂点を追加
			for(int i = 0; i <= slices; ++i){
				float s = float(i)/float(slices);
				float x = rj*sin(2*pi*s);
				float y = rj*cos(2*pi*s);
				vrts.push_back(glm::vec3(x, y, z-zlen));
				nrms.push_back(glm::normalize(glm::vec3(x, y, z)));
			}
		}
		for(int i = 0; i <= slices; ++i){
			float s = float(i)/float(slices);
			float x = rj*sin(2*pi*s);
			float y = rj*cos(2*pi*s);
			vrts.push_back(glm::vec3(x, y, z+zlen));
			nrms.push_back(glm::normalize(glm::vec3(x, y, z)));
		}
	}

	// メッシュ作成
	int nx = slices+1;
	for(int j = 0; j < stacks+1; ++j){
		for(int i = 0; i < slices; ++i){
			tris.push_back((i)+(j)*nx+offset);
			tris.push_back((i+1)+(j)*nx+offset);
			tris.push_back((i+1)+(j+1)*nx+offset);

			tris.push_back((i)+(j)*nx+offset);
			tris.push_back((i+1)+(j+1)*nx+offset);
			tris.push_back((i)+(j+1)*nx+offset);
		}
	}

	nvrts = static_cast<int>(vrts.size());
	ntris = static_cast<int>(tris.size()/3);

	return 1;
}

/*!
* n×nの頂点を持つ格子状の三角形メッシュ生成(x-z平面)(VBO用)
* @param[in] c1,c2 2端点座標
* @param[out] poly ポリゴンデータ
*/
inline static int MakeGridMesh(int &nvrts, vector<glm::vec3> &vrts, vector<glm::vec3> &nrms, int &ntris, vector<int> &tris,
	glm::vec3 c1, glm::vec3 c2, int nx, int ny)
{
	// 頂点座標生成
	float dx = (c2[0]-c1[0])/static_cast<float>(nx-1);
	float dz = (c2[2]-c1[2])/static_cast<float>(ny-1);
	for(int k = 0; k < ny; ++k){
		for(int i = 0; i < nx; ++i){
			glm::vec3 pos;
			pos[0] = c1[0]+i*dx;
			pos[1] = c1[1];
			pos[2] = c1[2]+k*dz;
			vrts.push_back(pos);
			nrms.push_back(glm::vec3(0, 1, 0));
		}
	}

	// メッシュ作成
	for(int k = 0; k < ny-1; ++k){
		for(int i = 0; i < nx-1; ++i){
			tris.push_back(i+k*nx);
			tris.push_back((i+1)+(k+1)*nx);
			tris.push_back((i+1)+k*nx);

			tris.push_back(i+k*nx);
			tris.push_back(i+(k+1)*nx);
			tris.push_back((i+1)+(k+1)*nx);
		}
	}

	nvrts = static_cast<int>(vrts.size());
	ntris = static_cast<int>(tris.size()/3);

	return 1;
}


/*!
* VAOによるプリミティブ描画
* @param[in] obj 頂点数,ポリゴン数情報を含むVAO
* @param[in] draw 描画フラグ
*/
inline static void DrawPrimitive(const GLuint vao, const int nvrts, const int ntris, int draw)
{
	// エッジ描画における"stitching"をなくすためのオフセットの設定
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);

	// 図形の描画
	glDisable(GL_LIGHTING);
	glBindVertexArray(vao);
	if(draw & 0x01){
		glColor3d(1.0, 0.3, 0.3);
		glPointSize(5.0);
		glDrawArrays(GL_POINTS, 0, nvrts);
	}
	if(draw & 0x02){
		glColor3d(0.5, 0.9, 0.9);
		glLineWidth(4.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDrawElements(GL_TRIANGLES, ntris*3, GL_UNSIGNED_INT, 0);
	}
	if(draw & 0x04){
		glEnable(GL_LIGHTING);
		//glDisable(GL_CULL_FACE);
		glEnable(GL_AUTO_NORMAL);
		glEnable(GL_NORMALIZE);
		glColor3d(0.1, 0.5, 1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glDrawElements(GL_TRIANGLES, ntris*3, GL_UNSIGNED_INT, 0);
	}
	glBindVertexArray(0);
}




/*!
* 原点中心の立方体形状描画
* - 立方体の中心は原点(0,0,0)
* - 辺の長さは1で固定 (glScaleで調整する)
* - VBO使用バージョン
*/
static inline void DrawCubeVBO(void)
{
	static MeshVBO vbo;
	if(vbo.vrts == 0){
		// ポリゴンデータ作成
		int nvrts, ntris;
		vector<glm::vec3> vrts, nrms;
		vector<int> tris;
		MakeCubeWithFaceNormal(nvrts, vrts, nrms, ntris, tris, 1.0);

		// VBOの作成
		CreateVBO(vbo, (GLfloat*)&vrts[0], nvrts, 3, &tris[0], ntris, 4, (GLfloat*)&nrms[0], nvrts);
	}
	glShadeModel(GL_SMOOTH);
	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.vrts);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.nrms);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo.tris);
	glDrawElements(GL_QUADS, vbo.ntris*4, GL_UNSIGNED_INT, 0);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}
/*!
* 原点中心の球形状描画
* - 中心は原点(0,0,0)
* - 半径はrad (glScaleで調整する)
* - VBO使用バージョン
*/
static inline void DrawSphereVBO(void)
{
	static MeshVBO vbo;
	if(vbo.vrts == 0){
		// ポリゴンデータ作成
		int nvrts, ntris;
		vector<glm::vec3> vrts, nrms;
		vector<int> tris;
		MakeSphere(nvrts, vrts, nrms, ntris, tris, 0.5, 32, 16);

		// VBOの作成
		CreateVBO(vbo, (GLfloat*)&vrts[0], nvrts, 3, &tris[0], ntris, 3, (GLfloat*)&nrms[0], nvrts);
	}
	glShadeModel(GL_SMOOTH);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.vrts);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.nrms);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo.tris);
	glDrawElements(GL_TRIANGLES, vbo.ntris*3, GL_UNSIGNED_INT, 0);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

//海老沢追加--------------------------------------
static inline void DrawCollisionSphereVBO(float x,float y,float z,float rad)
{
	glPushMatrix();
	glColor3f(1.0, 1.0, 1.0);
	glTranslatef(x,y,z);
	static MeshVBO vbo;
	if (vbo.vrts == 0) {
		// ポリゴンデータ作成
		int nvrts, ntris;
		vector<glm::vec3> vrts, nrms;
		vector<int> tris;
		MakeSphere(nvrts, vrts, nrms, ntris, tris, rad, 32, 16);

		// VBOの作成
		CreateVBO(vbo, (GLfloat*)&vrts[0], nvrts, 3, &tris[0], ntris, 3, (GLfloat*)&nrms[0], nvrts);
	}
	glShadeModel(GL_SMOOTH);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.vrts);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.nrms);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo.tris);
	glDrawElements(GL_TRIANGLES, vbo.ntris * 3, GL_UNSIGNED_INT, 0);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
	glPopMatrix();
}
//--------------------------------------------------------------------

/*!
* 原点中心の円筒形状描画
* - 中心は原点(0,0,0)
* - 半径0.5(直径1)/長さ1で固定 (glScaleで調整する)
* - VBO使用バージョン
*/
static inline void DrawCylinderVBO(void)
{
	static MeshVBO vbo;
	if(vbo.vrts == 0){
		// ポリゴンデータ作成
		int nvrts, ntris;
		vector<glm::vec3> vrts, nrms;
		vector<int> tris;
		MakeCylinder(nvrts, vrts, nrms, ntris, tris, 0.5, 0.5, 1.0, 16, true);

		// VBOの作成
		CreateVBO(vbo, (GLfloat*)&vrts[0], nvrts, 3, &tris[0], ntris, 3, (GLfloat*)&nrms[0], nvrts);
	}
	glShadeModel(GL_SMOOTH);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.vrts);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.nrms);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo.tris);
	glDrawElements(GL_TRIANGLES, vbo.ntris*3, GL_UNSIGNED_INT, 0);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

/*!
* 原点中心のカプセル形状描画
* - 中心は原点(0,0,0)
* - 円筒の両端面に半球が付いた形状
* - VBO使用バージョン
*/
static inline void DrawCapsuleVBO(float rad, float len)
{
	static MeshVBO vbo_c, vbo_s;
	if(vbo_c.vrts == 0){
		// 円筒部分ポリゴンデータ/VBO作成
		int nvrts, ntris;
		vector<glm::vec3> vrts, nrms;
		vector<int> tris;
		MakeCylinder(nvrts, vrts, nrms, ntris, tris, 0.5, 0.5, 1.0, 16, 8);
		CreateVBO(vbo_c, (GLfloat*)&vrts[0], nvrts, 3, &tris[0], ntris, 3, (GLfloat*)&nrms[0], nvrts);

		// 球部分ポリゴンデータ/VBO作成
		vrts.clear(), nrms.clear(); tris.clear();
		MakeSphere(nvrts, vrts, nrms, ntris, tris, 0.5, 16, 8);
		CreateVBO(vbo_s, (GLfloat*)&vrts[0], nvrts, 3, &tris[0], ntris, 3, (GLfloat*)&nrms[0], nvrts);
	}
	glShadeModel(GL_SMOOTH);

	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);

	// 円筒部分
	glPushMatrix();
	glScalef(2*rad, 2*rad, len-2*rad);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_c.vrts);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo_c.nrms);
	glNormalPointer(GL_FLOAT, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_c.tris);
	glDrawElements(GL_TRIANGLES, vbo_c.ntris*3, GL_UNSIGNED_INT, 0);
	glPopMatrix();

	// 端面球1(z-)
	glPushMatrix();
	glTranslatef(0, 0, -0.5*len+rad);
	glScalef(2*rad, 2*rad, 2*rad);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_s.vrts);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo_s.nrms);
	glNormalPointer(GL_FLOAT, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_s.tris);
	glDrawElements(GL_TRIANGLES, vbo_s.ntris*3, GL_UNSIGNED_INT, 0);
	glPopMatrix();

	// 端面球2(z+)
	glPushMatrix();
	glTranslatef(0, 0, 0.5*len-rad);
	glScalef(2*rad, 2*rad, 2*rad);
	glBindBuffer(GL_ARRAY_BUFFER, vbo_s.vrts);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo_s.nrms);
	glNormalPointer(GL_FLOAT, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_s.tris);
	glDrawElements(GL_TRIANGLES, vbo_s.ntris*3, GL_UNSIGNED_INT, 0);
	glPopMatrix();



	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

/*!
* 原点中心の円錐形状描画
* - 中心は原点(0,0,0)
* - 半径0.5(直径1)/長さ1で固定 (glScaleで調整する)
* - 基本的にCylinderと同じ，片方の端面の半径が0の円筒と考える
* - VBO使用バージョン
*/
static inline void DrawConeVBO(void)
{
	static MeshVBO vbo;
	if(vbo.vrts == 0){
		// ポリゴンデータ作成
		int nvrts, ntris;
		vector<glm::vec3> vrts, nrms;
		vector<int> tris;
		MakeCylinder(nvrts, vrts, nrms, ntris, tris, 0.0, 0.5, 1.0, 16, true);

		// VBOの作成
		CreateVBO(vbo, (GLfloat*)&vrts[0], nvrts, 3, &tris[0], ntris, 3, (GLfloat*)&nrms[0], nvrts);
	}

	glBindBuffer(GL_ARRAY_BUFFER, vbo.vrts);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.nrms);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo.tris);
	glDrawElements(GL_TRIANGLES, vbo.ntris*3, GL_UNSIGNED_INT, 0);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

/*!
* y方向を法線とする平面ポリゴン描画
* @param[in] y,s 床の高さと水平方向の長さ
* @param[in] use_tex テクスチャON/OFF
*/
static void DrawPlaneVBO(double y = 0.0, double s = 20.0, bool use_tex = true)
{
	static MeshVBO vbo;
	if(vbo.vrts == 0){
		// ポリゴンデータ作成
		int nvrts, ntris;
		vector<glm::vec3> vrts, nrms;
		vector<glm::vec2> texcoords;
		vector<int> tris;
		MakePlaneY(nvrts, vrts, nrms, ntris, tris, texcoords, s);

		// VBOの作成
		CreateVBO(vbo, (GLfloat*)&vrts[0], nvrts, 3, &tris[0], ntris, 4, (GLfloat*)&nrms[0], nvrts, 0, 0, (GLfloat*)&texcoords[0], nvrts);
	}

	static GLuint texFloor = 0;				//!< 床のテクスチャ
	if(use_tex && texFloor == 0){
		// 床テクスチャ読み込み
		glActiveTexture(GL_TEXTURE0);
		//loadTexture("floortile.bmp", texFloor, true, false);
		makeCheckerBoardTexture(texFloor, 8, true, false);
		glBindTexture(GL_TEXTURE_2D, texFloor);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 16.0f);
	}

	if(use_tex){
		glEnable(GL_TEXTURE_2D);
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texFloor);
	}

	glPushMatrix();

	glTranslatef(0, y, 0);
	//glScalef(s, 1, s);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.vrts);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.nrms);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, 0);

	glBindBuffer(GL_ARRAY_BUFFER, vbo.tcds);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glTexCoordPointer(2, GL_FLOAT, 0, 0);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo.tris);
	glDrawElements(GL_QUADS, vbo.ntris*4, GL_UNSIGNED_INT, 0);

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

	glPopMatrix();

	if(use_tex){
		glBindTexture(GL_TEXTURE_2D, 0);
	}
}

/*!
* y方向を法線とする平面ポリゴン描画
* @param[in] y,s 床の高さと水平方向の長さ
* @param[in] use_tex テクスチャON/OFF
*/
static void DrawStaticPlane(double y = 0.0, double s = 20.0, bool use_tex = true)
{
	static GLuint texFloor = 0;				//!< 床のテクスチャ
	if(use_tex && texFloor == 0){
		// 床テクスチャ読み込み
		glActiveTexture(GL_TEXTURE0);
		//loadTexture("floortile.bmp", texFloor, true, false);
		makeCheckerBoardTexture(texFloor, 8, true, false);
		glBindTexture(GL_TEXTURE_2D, texFloor);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 16.0f);
	}

	if(use_tex){
		glEnable(GL_TEXTURE_2D);
		glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, texFloor);
	}

	// 床用ポリゴン描画
	glNormal3d(0.0, 1.0, 0.0);
	glBegin(GL_QUADS);
	{
		glTexCoord2d(0.0, 0.0);	glVertex3d(-s, y, -s);
		glTexCoord2d(0.0, 1.0);	glVertex3d(-s, y, s);
		glTexCoord2d(1.0, 1.0);	glVertex3d(s, y, s);
		glTexCoord2d(1.0, 0.0);	glVertex3d(s, y, -s);
	}
	glEnd();

	if(use_tex){
		glBindTexture(GL_TEXTURE_2D, 0);
	}
}





//-----------------------------------------------------------------------------
// VAO関連
//-----------------------------------------------------------------------------
/*!
* 頂点配列オブジェクトの作成
* - 最低限頂点情報があれば良い．必要のない引数は0を入れておく．
* @param[in] vrts,nvrts 頂点座標配列と頂点数
* @param[in] dim 次元(2or3)
* @param[in] tris,ntris 三角形ポリゴンを構成する頂点インデックス列を格納した配列と三角形ポリゴン数
* @param[in] nrms,nnrms 各頂点の法線情報配列と法線数(=頂点数)
* @param[in] cols,ncols 各頂点の色情報配列と色情報数(=頂点数)
* @param[in] tcds,ntcds 各頂点のテクスチャ座標情報配列とテクスチャ座標数(=頂点数)
* @return VAOオブジェクト
*/
static inline GLuint CreateVAO(
	const GLfloat *vrts, GLuint nvrts, int dim = 3,
	const int     *tris = 0, GLuint ntris = 0,
	const GLfloat *nrms = 0, GLuint nnrms = 0,
	const GLfloat *cols = 0, GLuint ncols = 0,
	const GLfloat *tcds = 0, GLuint ntcds = 0)
{
	// VAOの生成
	GLuint vao;
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);

	// VBO:頂点座標
	GLuint vvbo;
	glGenBuffers(1, &vvbo);
	glBindBuffer(GL_ARRAY_BUFFER, vvbo);
	glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*dim*nvrts, vrts, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, dim, GL_FLOAT, GL_FALSE, 0, 0);

	// VBO:頂点テクスチャ座標
	if(tcds){
		GLuint tvbo;
		glGenBuffers(1, &tvbo);
		glBindBuffer(GL_ARRAY_BUFFER, tvbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*2*ntcds, tcds, GL_STATIC_DRAW);
		glEnableVertexAttribArray(1);
		glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, 0);
	}

	// VBO:頂点法線
	if(nrms){
		GLuint nvbo;
		glGenBuffers(1, &nvbo);
		glBindBuffer(GL_ARRAY_BUFFER, nvbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*dim*nnrms, nrms, GL_STATIC_DRAW);
		glEnableVertexAttribArray(2);
		glVertexAttribPointer(2, dim, GL_FLOAT, GL_FALSE, 0, 0);
	}

	// VBO:頂点描画色
	if(cols){
		GLuint cvbo;
		glGenBuffers(1, &cvbo);
		glBindBuffer(GL_ARRAY_BUFFER, cvbo);
		glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat)*3*ncols, cols, GL_STATIC_DRAW);
		glEnableVertexAttribArray(3);
		glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, 0, 0);
	}

	// VBO:三角形ポリゴンインデックス
	if(tris){
		GLuint pvbo;
		glGenBuffers(1, &pvbo);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pvbo);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(int)*ntris*3, tris, GL_STATIC_DRAW);
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	return vao;
}


/*!
* VAOに関連づけられたデータの更新
* @param[in] vao VAOオブジェクト
* @param[in] index glVertexAttribPointerで関連づけられたインデックス
* @param[in] data データ配列
* @param[in] n データ数
* @return VBOオブジェクト
*/
static inline GLuint UpdateDataVAO(GLuint vao, GLuint index, const GLfloat *data, GLuint n)
{
	glBindVertexArray(vao);
	GLuint vbo;
	glGetVertexAttribIuiv(index, GL_VERTEX_ATTRIB_ARRAY_BUFFER_BINDING, &vbo);
	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(GLfloat)*n, data);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);
	return vbo;
}

/*!
* VAOの破棄
*  - リンクされているVBOなどは破棄されないようなのでそれらは個別に行う必要あり
* @param[in] vao VAOオブジェクト
*/
static inline void DeleteVAO(GLuint vao)
{
	glDeleteVertexArrays(1, &vao);
}

/*!
* VAOによるメッシュ描画
* @param[in] vao VAOオブジェクト
* @param[in] draw 描画フラグ
* @param[in] nvrts,ntris 頂点数,三角形ポリゴン数
* @param[in] col[3] 頂点,エッジ,面の色
*/
static inline void DrawMeshVAO(GLuint vao, int draw, int nvrts, int ntris, glm::vec3 col[3])
{
	if(!vao) return;

	glBindVertexArray(vao);
	if(draw & 0x04){
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		glEnable(GL_LIGHTING);

		glColor3fv(glm::value_ptr(col[2]));
		glDrawElements(GL_TRIANGLES, ntris*3, GL_UNSIGNED_INT, 0);
	}
	if(draw & 0x02){
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glDisable(GL_LIGHTING);

		glLineWidth(1.0);
		glColor3fv(glm::value_ptr(col[1]));
		glDrawElements(GL_TRIANGLES, ntris*3, GL_UNSIGNED_INT, 0);
	}
	if(draw & 0x01){
		glDisable(GL_LIGHTING);

		glPointSize(8.0);
		glColor3fv(glm::value_ptr(col[0]));
		glDrawArrays(GL_POINTS, 0, nvrts);
	}
	glBindVertexArray(0);
}





//-----------------------------------------------------------------------------
// 粒子描画関連
//-----------------------------------------------------------------------------
/*!
 * 粒子に格納されたベクトル量をGL_LINESで描画
 * @param[in] prts 粒子位置を格納した配列
 * @param[in] vels 粒子速度を格納した配列
 * @param[in] n 粒子数
 * @param[in] offset 描画する粒子の最初のインデックス
 * @param[in] c0,c1 線分の始点と終点の色
 * @param[in] len 線の長さ
 */
static inline void drawParticleVector(float* prts, float* vels, int* attr, int n, int stride, int offset = 0, float scale = 1.0f, 
									  const glm::vec3 c0 = glm::vec3(1.0), const glm::vec3 c1 = glm::vec3(0.0, 0.0, 1.0), 
									  int draw_attr = 0, double czf = 1000, double czb = -100)
{
	glBegin(GL_LINES);
	int k = offset*stride/sizeof(float);
	for(int i = 0; i < n; ++i){
		if(prts[k+2] < czf && prts[k+2] > czb && attr[i] == draw_attr){
			glColor3f(c0[0], c0[1], c0[2]);
			glVertex3f(prts[k], prts[k+1], prts[k+2]);
			glColor3f(c1[0], c1[1], c1[2]);
			glVertex3f(prts[k]+scale*vels[k], prts[k+1]+scale*vels[k+1], prts[k+2]+scale*vels[k+2]);
		}
		k += stride/sizeof(float);
	}
	glEnd();
}

/*!
* 粒子に格納されたベクトル量をGL_LINESで描画(VBO版)
* @param[in] n 粒子数
* @param[in] pbuf 粒子位置を格納したバッファ
* @param[in] vbuf 描画するベクトルを格納したバッファ
* @param[in] stride 1粒子のメモリ上のサイズ(ストライド)
* @param[in] offset 描画する粒子の最初のインデックス
* @param[in] scale ベクトル描画時のスケール
*/
static inline void drawVectors(float* pbuf, float *vbuf, int n, int stride, int offset = 0, float scale = 1.0f)
{
	if (n <= 0) return;

	// ベクトルを表す線分の端点座標の計算
	vector<glm::vec3> vs(n*2);
	for (int i = 0; i < n; ++i) {
		vs[2*i] = glm::vec3(*pbuf, *(pbuf+1), *(pbuf+2));
		vs[2*i+1] = vs[2*i]+glm::vec3(*vbuf, *(vbuf+1), *(vbuf+2))*scale;
		pbuf += stride/sizeof(float); 
		vbuf += stride/sizeof(float);
	}

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glVertexPointer(3, GL_FLOAT, sizeof(glm::vec3), &vs[0]);
	glDrawArrays(GL_LINES, offset, vs.size());
}

/*!
 * 粒子をGL_POINTSで描画
 *  - VBOがあれば用いる
 * @param[in] vbo 粒子位置を格納したVBO
 * @param[in] n 粒子数
 * @param[in] color_vbo 粒子の色を格納したVBO
 * @param[in] data 粒子座標(ホスト側メモリ，vboが0の時に使用)
 * @param[in] offset 描画する粒子の最初のインデックス
 */
static inline void drawParticlePoints(unsigned int vbo, int n, unsigned int color_vbo, float* data, int offset = 0, 
									  glm::vec3 col = glm::vec3(0.0, 0.0, 1.0), GLfloat point_size = 1.0f)
{
	if(n <= 0) return;
	glColor4d(col[0], col[1], col[2], 1.0);
	glPointSize(point_size);
	if(!vbo){
		glBegin(GL_POINTS);
		int k = offset*3;
		for(int i = 0; i < n; ++i){
			glVertex3f(data[k], data[k+1], data[k+2]);
			k += 3;
		}
		glEnd();
	} else{
		glBindBuffer(GL_ARRAY_BUFFER, vbo);
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_FLOAT, 0, 0);

		if(color_vbo){
			glBindBufferARB(GL_ARRAY_BUFFER_ARB, color_vbo);
			glColorPointer(4, GL_FLOAT, 0, 0);
			glEnableClientState(GL_COLOR_ARRAY);
		}

		glDrawArrays(GL_POINTS, offset, n);//LINESに変更で描画可能

		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
		glBindBuffer(GL_ARRAY_BUFFER, 0);
	}
}


/*!
 * 粒子をPointSpriteで描画
 *  - VBOがあれば用いる
 * @param[in] vbo 粒子位置を格納したVBO
 * @param[in] n 粒子数
 * @param[in] color_vbo 粒子の色を格納したVBO
 * @param[in] data 粒子座標(ホスト側メモリ，vboが0の時に使用)
 * @param[in] offset 描画する粒子の最初のインデックス
 * @param[in] pscale,prad PointSpriteのスケールと半径
 * @param[in] czf,czb 断面を描画する際の前方,後方クリッピング位置(z座標値)
 */
static inline void drawPointSprites(unsigned int vbo, int n, unsigned int color_vbo, float* data, int offset, 
									double pscale, double prad, double czf = 1000, double czb = -1000)
{
	static bool init = true;
	static rxGLSL glslPS;			//!< GLSLを使った描画
	if(init){
		// PointSpriteシェーダの初期化(最初の呼び出し時だけ実行)
		//glslPS = CreateGLSL(ps_vs, ps_fs, "pointsprite");
		glslPS = CreateGLSLFromFile("shaders/pointsprites.vp", "shaders/pointsprites.fp", "pointsprites");
		init = false;
	}

	// PointSpriteのための設定
	glEnable(GL_POINT_SPRITE_ARB);
	glTexEnvi(GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE);
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE_NV);
	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);

	// シェーダ変数の設定
	glUseProgram(glslPS.Prog);
	glUniform1f(glGetUniformLocation(glslPS.Prog, "pointScale"), pscale);
	glUniform1f(glGetUniformLocation(glslPS.Prog, "pointRadius"), prad);
	glUniform3f(glGetUniformLocation(glslPS.Prog, "lightDir"), 0.6, 0.6, 0.6);
	glUniform1f(glGetUniformLocation(glslPS.Prog, "zCrossFront"), czf);
	glUniform1f(glGetUniformLocation(glslPS.Prog, "zCrossBack"),  czb);

	// 点描画
	drawParticlePoints(vbo, n, color_vbo, data, offset);

	// 後処理
	glUseProgram(0);
	glDisable(GL_POINT_SPRITE_ARB);

}

/*!
 * 粒子をglutSolidSphereで描画
 * @param[in] n 粒子数
 * @param[in] offset 描画する粒子の最初のインデックス
 * @param[in] data 粒子座標(ホスト側メモリ)
 * @param[in] attr 粒子属性(ホスト側メモリ)
 * @param[in] rad  描画半径
 * @param[in] col  描画色
 * @param[in] draw_attr どの属性を描画するか
 * @param[in] czf,czb 断面を描画する際の前方,後方クリッピング位置(z座標値)
 */
static inline void drawParticleSpheres(int n, int offset, float* data, int* attr, double rad, glm::vec3 col = glm::vec3(0.4, 0.4, 1.0), 
									   int draw_attr = 0, double czf = 1000, double czb = -1000)
{
	glEnable(GL_LIGHTING);
	glColor3f(col[0], col[1], col[2]);
	int k = offset*3;
	for(int i = 0; i < n; ++i){
		if(data[k+2] < czf && data[k+2] > czb && attr[i] == draw_attr){
			glPushMatrix();
			glTranslatef(data[k], data[k+1], data[k+2]);
			glRotated(90, 1.0, 0.0, 0.0);
			glScalef(2*rad, 2*rad, 2*rad);
			DrawSphereVBO();
			glPopMatrix();
		}
		k += 3;
	}
}



/*!
* 粒子のカラーバッファ値を計算
* @param[out] cvbo カラーVBO
* @param[in] hVal ホストメモリ上の配列
* @param[in] d 配列のステップ
* @param[in] n 要素数(配列のサイズはn*d)
* @param[in] use_max 最大値で値を正規化
* @param[in] vmax 手動設定最大値
* @return 値の最大値
*/
static inline float setColorVBO(GLuint& cvbo, float* hVal, int d, int n, bool use_max, float vmax)
{
	float l = 1.0;
	float max_val = 0.0;
	for(int i = 0; i < n; ++i){
		if(hVal[d*i] > max_val) max_val = hVal[d*i];
	}

	l = max_val;
	if(!use_max) l = vmax;

	// 粒子カラーバッファ
	glBindBuffer(GL_ARRAY_BUFFER, cvbo);
	float* data = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
	//float *ptr = data;
	for(int i = 0; i < n; ++i){
		float value = hVal[d*i]/l;

		float col[3];
		Gradation(col, value, 0.0, 1.0);
		data[4*i+0] = col[0];
		data[4*i+1] = col[1];
		data[4*i+2] = col[2];
	}
	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);

	return max_val;
}


/*!
* 粒子のカラーバッファ値を設定(一定値)
* @param[out] cvbo カラーVBO
* @param[in] n 要素数
* @param[in] col 色
*/
static inline int setColorVBO(GLuint& cvbo, int n, glm::vec3 col)
{
	// カラーバッファに値を設定
	glBindBuffer(GL_ARRAY_BUFFER, cvbo);
	float* data = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
	float* ptr = data;
	for(int i = 0; i < n; ++i){
		*ptr++ = col[0];
		*ptr++ = col[1];
		*ptr++ = col[2];
	}
	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	return 0;
}

/*!
* 粒子のカラーバッファ値を設定(グラデーション)
* @param[out] cvbo カラーVBO
* @param[in] n 要素数
*/
static inline int setRampColorVBO(GLuint& cvbo, int n)
{
	// カラーバッファに値を設定
	glBindBuffer(GL_ARRAY_BUFFER, cvbo);
	float* data = (float*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY);
	float* ptr = data;
	for(int i = 0; i < n; ++i){
		float t = i/(float)n;
		RX_COLOR_RAMP<float>(t, ptr);
		ptr += 3;
		*ptr++ = 1.0f;
	}
	glUnmapBuffer(GL_ARRAY_BUFFER);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	return 0;
}





//-----------------------------------------------------------------------------
// 各種オブジェクト描画関数
//-----------------------------------------------------------------------------
/*!
 * 2点からAABBを生成して描画
 * @param[in] A,B 直方体の2辺
 */
static inline void drawAABB(glm::vec3 A, glm::vec3 B)
{
	glBegin(GL_QUADS);

	//下
	glNormal3d(0.0, -1.0, 0.0);
	glVertex3d(A[0], A[1], B[2]);
	glVertex3d(A[0], A[1], A[2]);
	glVertex3d(B[0], A[1], A[2]);
	glVertex3d(B[0], A[1], B[2]);

	//左
	glNormal3d(-1.0, 0.0, 0.0);
	glVertex3d(A[0], B[1], A[2]);
	glVertex3d(A[0], A[1], A[2]);
	glVertex3d(A[0], A[1], B[2]);
	glVertex3d(A[0], B[1], B[2]);

	//手前
	glNormal3d(0.0, 0.0, 1.0);
	glVertex3d(A[0], B[1], B[2]);
	glVertex3d(A[0], A[1], B[2]);
	glVertex3d(B[0], A[1], B[2]);
	glVertex3d(B[0], B[1], B[2]);

	//右
	glNormal3d(1.0, 0.0, 0);
	glVertex3d(B[0], B[1], B[2]);
	glVertex3d(B[0], A[1], B[2]);
	glVertex3d(B[0], A[1], A[2]);
	glVertex3d(B[0], B[1], A[2]);

	//奥
	glNormal3d(0.0, -1.0, 0);
	glVertex3d(B[0], B[1], A[2]);
	glVertex3d(B[0], A[1], A[2]);
	glVertex3d(A[0], A[1], A[2]);
	glVertex3d(A[0], B[1], A[2]);

	//上
	glNormal3d(0.0, 1.0, 0);
	glVertex3d(A[0], B[1], A[2]);
	glVertex3d(A[0], B[1], B[2]);
	glVertex3d(B[0], B[1], B[2]);
	glVertex3d(B[0], B[1], A[2]);

	glEnd();

}

/*!
 * カップ形状の直方体(一部の面が空いているAABB)を描画
 * @param[in] minp,maxp 直方体の最小・最大座標
 * @param[in] w 壁面の厚さ
 * @param[in] drawface 描画しない面を指定する(下位ビットから手前(z+),奥(z-),右(x+),左(x-),上(y+),下(y-):0x01,02,04,08,10,20)
 */
inline void drawCubeCup(glm::vec3 minp, glm::vec3 maxp, double w, int drawface = 0x3F-0x10)
{
	if(drawface & 0x01) drawAABB(glm::vec3(minp[0], minp[1], maxp[2]), glm::vec3(maxp[0], maxp[1], maxp[2]+w));//手前
	if(drawface & 0x02) drawAABB(glm::vec3(minp[0], minp[1], minp[2]-w), glm::vec3(maxp[0], maxp[1], minp[2]));//奥

	double front_offset = (drawface & 0x01) ? 0.0 : -0.8*w;
	if(drawface & 0x04) drawAABB(glm::vec3(maxp[0], minp[1], minp[2]-w), glm::vec3(maxp[0]+w, maxp[1], maxp[2]+w+front_offset));//右
	if(drawface & 0x08) drawAABB(glm::vec3(minp[0]-w, minp[1], minp[2]-w), glm::vec3(minp[0], maxp[1], maxp[2]+w+front_offset));//左

	if(drawface & 0x10) drawAABB(glm::vec3(minp[0]-w, maxp[1]-w, minp[2]-w), glm::vec3(maxp[0]+w, maxp[1], maxp[2]+w+front_offset));//上
	if(drawface & 0x20) drawAABB(glm::vec3(minp[0]-w, minp[1]-w, minp[2]-w), glm::vec3(maxp[0]+w, minp[1], maxp[2]+w+front_offset));//下
}


/*!
 * 床描画
 *  - GLSLでテクスチャ付き＆柔らかいスポットライトで照らされているような床を描画
 * @param[in] light_pos,light_color 光源位置と色
 * @param[in] y,s 床の高さと水平方向の長さ
 * @return
 */
static inline void drawFloor(glm::vec3 light_pos, glm::vec3 light_color, double y = -1.0, double s = 20.0)
{
	static bool init = true;
	static rxGLSL glslFloor;			//!< GLSLを使った描画
	static GLuint texFloor = 0;				//!< 床のテクスチャ
	static GLuint texFloorShadow = 0;		//!< 床の影用テクスチャ(まだ未実装)
	if(init){
		// シェーダの初期化(最初の呼び出し時だけ実行)
		glslFloor = CreateGLSLFromFile("shaders/floor.vp", "shaders/floor.fp", "floor");
		init = false;

		// 床テクスチャ読み込み
		glActiveTexture(GL_TEXTURE0);
		//loadTexture("floortile.bmp", texFloor, true, false);
		makeCheckerBoardTexture(texFloor, 72, true, false);
		glBindTexture(GL_TEXTURE_2D, texFloor);
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_ANISOTROPY_EXT, 16.0f);

	}

	// 視点座標系のモデルビュー行列取得と視点座標系での光源位置計算
	GLfloat eye_modelview[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, eye_modelview);
	glm::vec3 light_pos_eye = MulMat4Vec3(eye_modelview, light_pos);

	glEnable(GL_TEXTURE_2D);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

	// GLSLシェーダをセット
	glUseProgram(glslFloor.Prog);

	glUniform3f(glGetUniformLocation(glslFloor.Prog, "v3LightPosEye"), light_pos_eye[0], light_pos_eye[1], light_pos_eye[2]);
	glUniform3f(glGetUniformLocation(glslFloor.Prog, "v3LightColor"), light_color[0], light_color[1], light_color[2]);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texFloor);
	glUniform1i(glGetUniformLocation(glslFloor.Prog, "texFloor"), 0);

	glEnable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);	// 裏面をカリング

	// set shadow matrix as texture matrix
	//matrix4f shadowMatrix = renderer->getShadowMatrix();
	//glActiveTexture(GL_TEXTURE0);
	//glMatrixMode(GL_TEXTURE);
	//glLoadMatrixf((GLfloat *)shadowMatrix.get_value());

	// 床用ポリゴン描画
	glColor3d(1.0, 1.0, 1.0);
	glNormal3d(0.0, 1.0, 0.0);
	glBegin(GL_QUADS);
	{
		glTexCoord2d(0.0, 0.0);	glVertex3d(-s, y, -s);
		glTexCoord2d(0.0, 1.0);	glVertex3d(-s, y,  s);
		glTexCoord2d(1.0, 1.0);	glVertex3d( s, y,  s);
		glTexCoord2d(1.0, 0.0);	glVertex3d( s, y, -s);
	}
	glEnd();

	glUseProgram(0);
	glBindTexture(GL_TEXTURE_2D, 0);
}


/*!
 * xyz軸描画(x軸:赤,y軸:緑,z軸:青)
 * @param[in] len 軸の長さ
 */
static inline int drawAxis(double len, double line_width)
{
	glLineWidth((GLfloat)line_width);

	// x軸
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(len, 0.0, 0.0);
	glEnd();

	// y軸
	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, len, 0.0);
	glEnd();

	// z軸
	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3d(0.0, 0.0, 0.0);
	glVertex3d(0.0, 0.0, len);
	glEnd();

	return 1;
}


/*!
 * 直方体のワイヤーフレーム描画
 * @param[in] min 最小座標値
 * @param[in] max 最大座標値
 * @param[in] color 描画色
 */
static inline void drawWireCuboid(const glm::vec3& min, const glm::vec3& max, const glm::vec3& color, double line_width = 1.0)
{
	glLineWidth((GLfloat)line_width);

	glPushMatrix();
	glColor3d(color[0], color[1], color[2]);

	glBegin(GL_LINES);

	// x軸平行
	glVertex3f(min[0], min[1], min[2]);	glVertex3f(max[0], min[1], min[2]);
	glVertex3f(min[0], min[1], max[2]); glVertex3f(max[0], min[1], max[2]);
	glVertex3f(min[0], max[1], min[2]);	glVertex3f(max[0], max[1], min[2]);
	glVertex3f(min[0], max[1], max[2]);	glVertex3f(max[0], max[1], max[2]);

	// z軸平行
	glVertex3f(min[0], min[1], min[2]);	glVertex3f(min[0], min[1], max[2]);
	glVertex3f(min[0], max[1], min[2]);	glVertex3f(min[0], max[1], max[2]);
	glVertex3f(max[0], min[1], min[2]);	glVertex3f(max[0], min[1], max[2]);
	glVertex3f(max[0], max[1], min[2]);	glVertex3f(max[0], max[1], max[2]);

	// z軸平行
	glVertex3f(min[0], min[1], min[2]);	glVertex3f(min[0], max[1], min[2]);
	glVertex3f(min[0], min[1], max[2]);	glVertex3f(min[0], max[1], max[2]);
	glVertex3f(max[0], min[1], min[2]);	glVertex3f(max[0], max[1], min[2]);
	glVertex3f(max[0], min[1], max[2]);	glVertex3f(max[0], max[1], max[2]);

	glEnd();

	glPopMatrix();
}



/*!
 * 原点中心の円のワイヤーフレーム描画
 * @param rad 円の半径
 * @param n 分割数
 */
static inline void DrawWireCircle(const double& rad, const int& n)
{
	double t = 0.0;
	double dt = 2.0*RX_PI/(double)n;

	glBegin(GL_LINE_LOOP);
	do{
		glVertex3f(rad*cos(t), rad*sin(t), 0.0);
		t += dt;
	} while(t < 2.0*RX_PI);
	glEnd();
}

/*!
 * 原点中心の円のワイヤーフレーム描画(XZ平面)
 * @param rad 円の半径
 * @param n 分割数
 */
static inline void DrawWireCircleXZ(const double& rad, const int& n)
{
	double t = 0.0;
	double dt = 2.0*RX_PI/(double)n;

	glBegin(GL_LINE_LOOP);
	do{
		glVertex3f(rad*cos(t), 0.0, rad*sin(t));
		t += dt;
	} while(t < 2.0*RX_PI);
	glEnd();
}

/*!
 * 球のワイヤーフレーム描画
 * @param cen 球の中心
 * @param rad 球の半径
 * @param col 描画色
 */
static inline void DrawWireSphere(const glm::vec3& cen, const float& rad, const glm::vec3& col)
{
	glDisable(GL_LIGHTING);
	glPushMatrix();
	glTranslatef(cen[0], cen[1], cen[2]);
	glRotatef(90, 1.0, 0.0, 0.0);
	glColor3f(col[0], col[1], col[2]);

	// 緯度(x-y平面に平行)
	float z, dz;
	dz = 2.0*rad/8.0f;
	z = -(rad-dz);
	do{
		glPushMatrix();
		glTranslatef(0.0, 0.0, z);
		DrawWireCircle(sqrt(rad*rad-z*z), 32);
		glPopMatrix();
		z += dz;
	} while(z < rad);

	// 経度(z軸まわりに回転)
	float t, dt;
	t = 0.0f;
	dt = 180.0/8.0;
	do{
		glPushMatrix();
		glRotatef(t, 0.0, 0.0, 1.0);
		DrawWireCircleXZ(rad, 32);
		glPopMatrix();

		t += dt;
	} while(t < 180);

	//glutWireSphere(rad, 10, 5);
	glPopMatrix();
}

/*!
 * 球描画
 * @param cen 球の中心
 * @param rad 球の半径
 * @param col 描画色
 */
//static inline void DrawSolidSphere(const glm::vec3& cen, const float& rad, const glm::vec3& col)
//{
//	glPushMatrix();
//	glPushMatrix();
//	glTranslatef(cen[0], cen[1], cen[2]);
//	glRotatef(90, 1.0, 0.0, 0.0);
//	glColor3f(col[0], col[1], col[2]);
//	glutSolidSphere(rad, 20, 10);
//	glPopMatrix();
//}

/*!
 * 直方体の8頂点座標を中心と辺の長さの1/2(side length)から計算
 * @param[in] cn,sl 立方体の中心と辺の長さの1/2(side length)
 * @param[out] v 8頂点座標
 */
inline void SetVerticesBox(const glm::vec3& cn, const glm::vec3& sl, glm::vec3 v[8])
{
	v[0] = cn+glm::vec3(-sl[0], -sl[1], -sl[2]);
	v[1] = cn+glm::vec3(-sl[0], sl[1], -sl[2]);
	v[2] = cn+glm::vec3(-sl[0], sl[1], sl[2]);
	v[3] = cn+glm::vec3(-sl[0], -sl[1], sl[2]);

	v[4] = cn+glm::vec3(sl[0], -sl[1], -sl[2]);
	v[5] = cn+glm::vec3(sl[0], sl[1], -sl[2]);
	v[6] = cn+glm::vec3(sl[0], sl[1], sl[2]);
	v[7] = cn+glm::vec3(sl[0], -sl[1], sl[2]);
}

/*!
 * 上が空いたボックス形状(壁面厚みあり)を表すポリゴンデータの生成
 * @param[in] sl0,sl1 ボックスの内側と外側のサイズ(side length)
 * @param[in] d 空ける方向
 * @param[out] vrts,idxs ポリゴン頂点座標と接続関係
 */
inline void CreateBoxPolygon(const glm::vec3& sl0, const glm::vec3& sl1, const int& d,
							 vector<glm::vec3>& vrts, vector< vector<int> >& idxs)
{
	if(d < 0 || d > 2) return;

	double h = sl1[d]-sl0[d];
	glm::vec3 cn(0.0);

	vrts.resize(16);

	// 外側の頂点
	SetVerticesBox(cn, sl1, &vrts[0]);

	// 内側の頂点
	cn[d] += h;
	SetVerticesBox(cn, sl0, &vrts[8]);


	int idxs0[5][4] = { {0, 3, 2, 1},
						{1, 2, 6, 5},
						{5, 6, 7, 4},
						{4, 7, 3, 0},
						{0, 1, 5, 4} };

	int idxs1[4][4] = { {2, 3, 11, 10},
						{3, 7, 15, 11},
						{7, 6, 14, 15},
						{6, 2, 10, 14} };

	// 三角形作成
	idxs.resize(28);
	for(int i = 0; i < 28; ++i) idxs[i].resize(3);

	int c = 0;

	// 外側の箱
	for(int i = 0; i < 5; ++i){
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs0[i][j];
		}
		c++;
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs0[i][((j+2 > 3) ? 0 : j+2)];
		}
		c++;
	}

	// 内側の箱
	for(int i = 0; i < 5; ++i){
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs0[i][2-j]+8;
		}
		c++;
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs0[i][(((2-j)+2 > 3) ? 0 : (2-j)+2)]+8;
		}
		c++;
	}

	// 上部
	for(int i = 0; i < 4; ++i){
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs1[i][j];
		}
		c++;
		for(int j = 0; j < 3; ++j){
			idxs[c][j] = idxs1[i][((j+2 > 3) ? 0 : j+2)];
		}
		c++;
	}

}

/*!
 * 上が空いたボックス形状(壁面厚みあり)の描画
 * @param[in] sl0,sl1 内側と外側のサイズ(side length)
 */
//static inline void DrawSolidOpenBox(glm::vec3 sl0, glm::vec3 sl1)
//{
//	glPushMatrix();
//
//	glm::vec3 len0 = 2.0*sl0;
//	glm::vec3 len1 = 2.0*sl1;
//	double d = 0.5*(len1[2]-len0[2]);
//
//	glTranslatef(0.0, 0.0, 0.5*len1[2]);
//
//	vector<glm::vec3> vrts;
//	vector< vector<int> > idxs;
//	CreateBoxPolygon(sl0, sl1, 2, vrts, idxs);
//
//	// インデックスを1始まりに
//	int n = (int)idxs.size();
//	for(int i = 0; i < n; ++i){
//		for(int j = 0; j < 3; ++j){
//			idxs[i][j]++;
//		}
//	}
//
//	// ワイヤーフレーム描画
//	glDisable(GL_LIGHTING);
//	glPushMatrix();
//	glTranslatef(0.0, 0.0, d);
//	glScalef(len0[0], len0[1], len0[2]);
//	glutWireCube(1.0);
//	glPopMatrix();
//
//	glPushMatrix();
//	glScalef(len1[0], len1[1], len1[2]);
//	glutWireCube(1.0);
//	glPopMatrix();
//
//	// 面描画
//	glEnable(GL_LIGHTING);
//	rxMaterial mat;
//	mat.SetGL();
//	glColor3f(0.0, 0.0, 1.0);
//	for(int i = 0; i < (int)idxs.size(); ++i){
//		glBegin(GL_POLYGON);
//		for(int j = 0; j < (int)idxs[i].size(); ++j){
//			glVertex3dv(vrts[idxs[i][j]-1].data);
//		}
//		glEnd();
//	}
//
//	glPopMatrix();
//}

/*!
 * 上が空いたボックス形状(壁面厚みあり)の描画
 * @param[in] sl0,sl1 内側と外側のサイズ(side length)
 */
//static inline void DrawWireOpenBox(glm::vec3 sl0, glm::vec3 sl1)
//{
//	glPushMatrix();
//
//	glm::vec3 len0 = 2.0*sl0;
//	glm::vec3 len1 = 2.0*sl1;
//	double d = 0.5*(len1[2]-len0[2]);
//
//	glTranslatef(0.0, 0.0, 0.5*len1[2]);
//
//	glPushMatrix();
//	glTranslatef(0.0, 0.0, d);
//	glScalef(len0[0], len0[1], len0[2]);
//	glutWireCube(1.0);
//	glPopMatrix();
//
//	glPushMatrix();
//	glScalef(len1[0], len1[1], len1[2]);
//	glutWireCube(1.0);
//	glPopMatrix();
//
//	glPopMatrix();
//}
//


//-----------------------------------------------------------------------------
// 文字列描画関数
//-----------------------------------------------------------------------------
// FTGL
const string FONT_FILE = "Inconsolata.ttf";
static FTPixmapFont* g_font = 0;
const unsigned long g_fontsize = 14;		//!< フォントサイズ


/*!
 * 文字列描画
 * @param[in] cir_str 文字列循環バッファ
 * @param[in] static_str 静的な文字列バッファ
 * @param[in] w,h ウィンドウサイズ
 */
static inline void drawStrings(vector<string>& static_str, int w, int h, int offsetx, int offsety)
{
	if(g_font == 0){
		// フォント読み込み
		g_font = new FTPixmapFont(FONT_FILE.c_str());
		if(g_font->Error()){
			cout << "Failed to open font " << FONT_FILE << endl;
			delete g_font;
			g_font = 0;
		} else g_font->FaceSize(g_fontsize);
	}

	glDisable(GL_LIGHTING);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, w, 0, h);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	float y0 = h-offsety-(static_str.size()*(g_font->LineHeight())+5);
	glRasterPos2f(20.0f+offsetx, y0);

	if(g_font){
		// FTGLで文字列を描画
		for(int j = 0; j < (int)static_str.size(); ++j){
			glRasterPos2i(20+offsetx, y0);
			g_font->Render(static_str[j].c_str());
			y0 -= g_font->LineHeight();
		}
	} else{
		// glutBitmapStringで文字列描画
		for(int j = 0; j < (int)static_str.size(); ++j){
			glRasterPos2i(20+offsetx, y0);
			//glutBitmapString(GLUT_BITMAP_HELVETICA_12, (unsigned char*)static_str[j].c_str());
			y0 -= 20;
		}
	}

	glRasterPos2f(0, 0);
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}


/*!
 * 文字列描画
 * @param[in] cir_str 文字列循環バッファ
 * @param[in] static_str 静的な文字列バッファ
 * @param[in] w,h ウィンドウサイズ
 */
static inline void drawStringsBottom(vector<string>& static_str, int w, int h, int offsetx, int offsety)
{
	if(g_font == 0){
		// フォント読み込み
		g_font = new FTPixmapFont(FONT_FILE.c_str());
		if(g_font->Error()){
			cout << "Failed to open font " << FONT_FILE << endl;
			delete g_font;
			g_font = 0;
		} else g_font->FaceSize(g_fontsize);
	}

	glDisable(GL_LIGHTING);
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D(0, w, 0, h);
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	float x0 = 20.0f;

	if(g_font){
		float y0 = static_str.size()*(g_font->LineHeight())+5;
		glRasterPos2f(20.0f+offsetx, y0);

		// FTGLで文字列描画
		for(int j = 0; j < (int)static_str.size(); ++j){
			glRasterPos2i(20+offsetx, y0);
			g_font->Render(static_str[j].c_str());

			y0 -= g_font->LineHeight();
		}
	} else{
		float y0 = static_str.size()*20.0f;
		glRasterPos2f(20.0f+offsetx, y0);

		// glutBitmapStringで文字列描画
		for(int j = 0; j < (int)static_str.size(); ++j){
			glRasterPos2i(20+offsetx, y0);

			//glutBitmapString(GLUT_BITMAP_HELVETICA_12, (unsigned char*)static_str[j].c_str());

			y0 -= 20;
		}
	}

	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
}




#endif // #ifndef _RX_GLDRAW_H_