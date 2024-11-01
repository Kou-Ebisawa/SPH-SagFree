/*! 
 @file rx_shader.h

 @brief GLSLシェーダー
 
 @author Makoto Fujisawa
 @date 2009-11
*/
// FILE --rx_shader.h--


#ifndef _RX_SHADERS_H_
#define _RX_SHADERS_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include <cstdio>
#include <iostream>
#include <vector>
#include <string>

#include <GL/glew.h>
#include <GLFW/glfw3.h>


#define RXSTR(A) #A

using namespace std;


//-----------------------------------------------------------------------------
// HACK:GLSLシェーダ
//-----------------------------------------------------------------------------
struct rxGLSL
{
	string VertProg;
	string FragProg;
	string Name;
	GLuint Prog;
};

//-----------------------------------------------------------------------------
// MARK:GLSLコンパイル
//-----------------------------------------------------------------------------
/*!
 * GLSLプログラムのコンパイル
 * @param[in] vsource vertexシェーダプログラム冁E��
 * @param[in] fsource pixel(fragment)シェーダプログラム冁E��
 * @return GLSLプログラム番号
 */
static GLuint CompileProgram(const char *vsource, const char *fsource)
{
	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

	glShaderSource(vertexShader, 1, &vsource, 0);
	glShaderSource(fragmentShader, 1, &fsource, 0);
	
	glCompileShader(vertexShader);
	glCompileShader(fragmentShader);

	GLuint program = glCreateProgram();

	glAttachShader(program, vertexShader);
	glAttachShader(program, fragmentShader);

	glLinkProgram(program);

	// check if program linked
	GLint success = 0;
	glGetProgramiv(program, GL_LINK_STATUS, &success);

	if (!success) {
		char temp[256];
		glGetProgramInfoLog(program, 256, 0, temp);
		printf("Failed to link program:\n%s\n", temp);
		glDeleteProgram(program);
		program = 0;
	}

	return program;
}

/*!
 * GLSLシェーダコンパイル
 * @param[in] target ターゲチE��(GL_VERTEX_SHADER,GL_FRAGMENT_SHADER)
 * @param[in] shader シェーダコーチE
 * @return GLSLオブジェクチE
 */
inline GLuint CompileGLSLShader(GLenum target, const char* shader)
{
	// GLSLオブジェクト作�E
	GLuint object = glCreateShader(target);

	if(!object) return 0;

	glShaderSource(object, 1, &shader, NULL);
	glCompileShader(object);

	// コンパイル状態�E確誁E
	GLint compiled = 0;
	glGetShaderiv(object, GL_COMPILE_STATUS, &compiled);

	if(!compiled){
		char temp[256] = "";
		glGetShaderInfoLog( object, 256, NULL, temp);
		fprintf(stderr, " Compile failed:\n%s\n", temp);

		glDeleteShader(object);
		return 0;
	}

	return object;
}

/*!
 * GLSLシェーダコンパイル
 * @param[in] target ターゲチE��(GL_VERTEX_SHADER,GL_FRAGMENT_SHADER)
 * @param[in] fn シェーダファイルパス
 * @return GLSLオブジェクチE
 */
inline GLuint CompileGLSLShaderFromFile(GLenum target, const char* fn)
{
	FILE *fp;

	// バイナリとしてファイル読み込み
	fp = fopen(fn, "rb");
	if(fp == NULL) return 0;

	// ファイルの末尾に移動し現在位置(ファイルサイズ)を取征E
	fseek(fp, 0, SEEK_END);
	long size = ftell(fp);

	fseek(fp, 0, SEEK_SET);

	// シェーダの冁E��格紁E
	char *text = new char[size+1];
	fread(text, size, 1, fp);
	text[size] = '\0';

	//printf("%s\n", text);


	fclose(fp);

	// シェーダコンパイル
	printf("Compile %s\n", fn);
	GLuint object = CompileGLSLShader(target, text);

	delete [] text;

	return object;
}

/*!
 * バ�EチE��クスとフラグメントシェーダで構�EされるGLSLプログラム作�E
 * @param[in] vs バ�EチE��クスシェーダオブジェクチE
 * @param[in] fs フラグメントシェーダオブジェクチE
 * @return GLSLプログラムオブジェクチE
 */
inline GLuint LinkGLSLProgram(GLuint vs, GLuint fs)
{
	// プログラムオブジェクト作�E
	GLuint program = glCreateProgram();

	// シェーダオブジェクトを登録
	glAttachShader(program, vs);
	glAttachShader(program, fs);

	// プログラムのリンク
	glLinkProgram(program);

	// エラー出劁E
	GLint charsWritten, infoLogLength;
	glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLogLength);

	char * infoLog = new char[infoLogLength];
	glGetProgramInfoLog(program, infoLogLength, &charsWritten, infoLog);
	printf("%s\n", infoLog);
	delete [] infoLog;

	// リンカチE��チE
	GLint linkSucceed = GL_FALSE;
	glGetProgramiv(program, GL_LINK_STATUS, &linkSucceed);
	if(linkSucceed == GL_FALSE){
		glDeleteProgram(program);
		return 0;
	}

	return program;
}


/*!
 * バ�EチE��クス/ジオメトリ/フラグメントシェーダで構�EされるGLSLプログラム作�E
 * @param[in] vs バ�EチE��クスシェーダオブジェクチE
 * @param[in] gs ジオメトリシェーダオブジェクチE
 * @param[in] inputType ジオメトリシェーダへの入力タイチE
 * @param[in] vertexOut バ�EチE��クスの出劁E
 * @param[in] outputType ジオメトリシェーダからの出力タイチE
 * @param[in] fs フラグメントシェーダオブジェクチE
 * @return GLSLプログラムオブジェクチE
 */
inline GLuint LinkGLSLProgram(GLuint vs, GLuint gs, GLint inputType, GLint vertexOut, GLint outputType, GLuint fs)
{
	// プログラムオブジェクト作�E
	GLuint program = glCreateProgram();

	// シェーダオブジェクトを登録
	glAttachShader(program, vs);
	glAttachShader(program, gs);

	glProgramParameteriEXT(program, GL_GEOMETRY_INPUT_TYPE_EXT, inputType);
	glProgramParameteriEXT(program, GL_GEOMETRY_VERTICES_OUT_EXT, vertexOut);
	glProgramParameteriEXT(program, GL_GEOMETRY_OUTPUT_TYPE_EXT, outputType);
	glAttachShader(program, fs);

	// プログラムのリンク
	glLinkProgram(program);

	// エラー出劁E
	GLint charsWritten, infoLogLength;
	glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLogLength);

	char * infoLog = new char[infoLogLength];
	glGetProgramInfoLog(program, infoLogLength, &charsWritten, infoLog);
	printf("%s\n", infoLog);
	delete [] infoLog;

	// リンカチE��チE
	GLint linkSucceed = GL_FALSE;
	glGetProgramiv(program, GL_LINK_STATUS, &linkSucceed);
	if(linkSucceed == GL_FALSE){
		glDeleteProgram(program);
		return 0;
	}

	return program;
}


/*!
 * GLSLのコンパイル・リンク(ファイルより)
 * @param[in] vs 頂点シェーダファイルパス
 * @param[in] fs フラグメントシェーダファイルパス
 * @param[in] name プログラム吁E
 * @return GLSLオブジェクチE
 */
inline rxGLSL CreateGLSLFromFile(const string &vs, const string &fs, const string &name)
{
	rxGLSL gs;
	gs.VertProg = vs;
	gs.FragProg = fs;
	gs.Name = name;
	GLuint v, f;
	printf("compile the vertex shader : %s\n", name.c_str());
	if(!(v = CompileGLSLShaderFromFile(GL_VERTEX_SHADER, vs.c_str()))){
		// skip the first three chars to deal with path differences
		v = CompileGLSLShaderFromFile(GL_VERTEX_SHADER, &vs.c_str()[3]);
	}

	printf("compile the fragment shader : %s\n", name.c_str());
	if(!(f = CompileGLSLShaderFromFile(GL_FRAGMENT_SHADER, fs.c_str()))){
		f = CompileGLSLShaderFromFile(GL_FRAGMENT_SHADER, &fs.c_str()[3]);
	}

	gs.Prog = LinkGLSLProgram(v, f);
	//gs.Prog = GLSL_CreateShaders(gs.VertProg.c_str(), gs.FragProg.c_str());

	return gs;
}

/*!
 * #versionなどのプリプロセチE��を文字�Eとして書かれたシェーダ中に含む場合，改行がぁE��くいかなぁE�Eで�E�E
 *  #version 110 @ のように最後に@を付け�E�改行に変換する
 * @param[in] s  シェーダ斁E���E
 * @param[in] vs 変換後�Eシェーダ斁E���E
 */
inline void CreateGLSLShaderString(const char* s, vector<char> &vs)
{
	int idx = 0;
	char c = s[0];
	while(c != '\0'){
		if(c == '@') c = '\n'; // #versionなどを可能にするために@を改行に変換

		vs.push_back(c);
		idx++;
		c = s[idx];
	}
	vs.push_back('\0');
}

/*!
 * GLSLのコンパイル・リンク(斁E���Eより)
 * @param[in] vs 頂点シェーダ冁E��
 * @param[in] fs フラグメントシェーダ冁E��
 * @param[in] name プログラム吁E
 * @return GLSLオブジェクチE
 */
inline rxGLSL CreateGLSL(const char* vs, const char* fs, const string &name)
{
	rxGLSL gs;
	gs.VertProg = "from char";
	gs.FragProg = "from char";
	gs.Name = name;

	vector<char> vs1, fs1;
	CreateGLSLShaderString(vs, vs1);
	CreateGLSLShaderString(fs, fs1);
	
	//printf("vertex shader : %d\n%s\n", vs1.size(), &vs1[0]);
	//printf("pixel shader  : %d\n%s\n", fs1.size(), &fs1[0]);

	GLuint v, f;
	printf("compile the vertex shader : %s\n", name.c_str());
	v = CompileGLSLShader(GL_VERTEX_SHADER, &vs1[0]);
	printf("compile the fragment shader : %s\n", name.c_str());
	f = CompileGLSLShader(GL_FRAGMENT_SHADER, &fs1[0]);
	gs.Prog = LinkGLSLProgram(v, f);

	return gs;
}



#endif // #ifndef _RX_SHADERS_H_