/*!
  @file scene.cpp

  @brief GLFWによるOpenGL描画

  @author Makoto Fujisawa
  @date   2021-11
*/


#pragma warning (disable: 4996)
#pragma warning (disable: 4819)


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "scene.h"

// テクスチャ・画像
#include "rx_bitmap.h"
#include "rx_texture.h"

// OpenGL描画関連
#include "gldraw.h"
#include "rx_trackball.h"	// 視点変更用トラックボールクラス
#include "rx_shaders.h" 

// SPHによる流体シミュレーション
#include "sph.h"

// MC法による陰関数場からの三角形メッシュ生成
#include "mc.h"

#include <random>
#include <time.h>
#include <cmath>

//髪型アセットの作成
#include "rx_obj.h"

#include <ctime>
#include <cstdlib>

#include <windows.h>


//-----------------------------------------------------------------------------
// 定数・変数
//-----------------------------------------------------------------------------
const glm::vec4 RX_LIGHT0_POS(2.0f, 4.0f, 2.0f, 0.0f); //SPHの場合(0.0f,0.0f,1.0f,1.0f)
const glm::vec4 RX_LIGHT1_POS(-1.0f, -10.0f, -1.0f, 0.0f);
const glm::vec4 RX_LIGHT_AMBI(0.3f, 0.3f, 0.3f, 1.0f);
const glm::vec4 RX_LIGHT_DIFF(1.0f, 1.0f, 1.0f, 1.0f);
const glm::vec4 RX_LIGHT_SPEC(1.0f, 1.0f, 1.0f, 1.0f);

//海老沢追加
const GLfloat BrownAmbi[] = { 0.00f, 0.00f, 0.00f, 0.00f };
const glm::vec4 BrownDiff(0.25f, 0.15f, 0.12f, 1.00f);
const glm::vec4 BrownSpec(0.72f, 0.52f, 0.34f, 1.00f);
const float BrownShine(0.00f);

#define PI_F 3.14159265359f

//PBDのために配列に変換
const GLfloat RX_FOV = 45.0f;

// 流体シミュ関連
SPH *g_sim = 0;//そのまま利用
float g_dt = 0.001;		//!< 時間ステップ幅
int g_colortype = SPH::C_CONSTANT;
float g_bnd_scale = 1.0f;

rxPolygons g_poly;

//デバック用の表示する毛髪をランダムに選択
int g_display_hair = -1;

//サイズ感を表すパラメータ(12/20追加)
float g_simulation_size_scene = 2.5f;

//サイズ感を変えたことによる重力の変更
glm::vec3 g_user_gravity = glm::vec3(0.f, -9.81f, 0.f) / 10.f;

//画像を毎ステップ保存するかどうかのフラグ
bool g_spacing_flag = false;

//風を適用するかどうかのフラグ
bool g_Apply_flag = false;

//連番OBJを保存するかどうかのフラグ
bool g_obj_Animation_flag = false;

//実験時の風の強さの設定
glm::vec3 g_experiment_power = glm::vec3(0.f);

//-----------------------------------------------------------------------------
// ScenePBDクラスのstatic変数の定義と初期化
//-----------------------------------------------------------------------------
int ScenePBD::m_winw = 1000;					//!< 描画ウィンドウの幅
int ScenePBD::m_winh = 1000;					//!< 描画ウィンドウの高さ
rxTrackball ScenePBD::m_view;					//!< トラックボール
float ScenePBD::m_bgcolor[3] = { 1, 1, 1 };	//!< 背景色
bool ScenePBD::m_animation_on = false;			//!< アニメーションON/OFF

int ScenePBD::m_draw = 0;						//!< 描画フラグ
int ScenePBD::m_currentstep = 0;				//!< 現在のステップ数
int ScenePBD::m_simg_spacing = -1;				//!< 画像保存間隔(=-1なら保存しない)

ElasticPBD* ScenePBD::m_elasticbody = 0;
//弾性体の配列を設定
vector<ElasticPBD*> ScenePBD::m_elasticbodies;

//CPUからGPUにデータを送るための配列
vector<glm::vec3> ScenePBD::m_NewPos;
vector<glm::vec3> ScenePBD::m_NewVelocity;

int ScenePBD::m_num_elasticbodies = 0;//弾性体の数
int ScenePBD::m_allnum = 0;//全ての粒子数の合計

int ScenePBD::m_picked = -1;
float ScenePBD::m_pickdist = 1.0;	//!< ピックされた点までの距離

//海老沢追加
SceneConfig ScenePBD::sph_scene;
int ScenePBD::sph_scene_idx;

//海老沢追加
void ScenePBD::drawHairObject(unsigned int vbo, int n, unsigned int color_vbo,unsigned int tang_vbo, float* data, int offset, double pscale, double prad, double czf, double czb) {
	static bool init = true;
	static rxGLSL glslPS;			//!< GLSLを使った描画
	if (init) {
		// PointSpriteシェーダの初期化(最初の呼び出し時だけ実行)
		//glslPS = CreateGLSL(ps_vs, ps_fs, "pointsprite");
		glslPS = CreateGLSLFromFile("shaders/shading.vp", "shaders/shading.fp", "hair_shading");
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
	/*glUniform1f(glGetUniformLocation(glslPS.Prog, "pointScale"), pscale);
	glUniform1f(glGetUniformLocation(glslPS.Prog, "pointRadius"), prad);
	glUniform3f(glGetUniformLocation(glslPS.Prog, "lightDir"), 0.6, 0.6, 0.6);
	glUniform1f(glGetUniformLocation(glslPS.Prog, "zCrossFront"), czf);
	glUniform1f(glGetUniformLocation(glslPS.Prog, "zCrossBack"), czb);*/

	//draw particle point-------------------------------
	//glm::vec3 col(1.0, 0.0, 0.0);
	GLfloat line_width = 3.0;
	//glColor4d(col[0], col[1], col[2], 1.0);
	glLineWidth(line_width);

	glBindBuffer(GL_ARRAY_BUFFER, vbo);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, 0, 0);

	//法線を接線として保存-----------------------------
	glBindBuffer(GL_ARRAY_BUFFER, tang_vbo);
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, 0, 0);
	//-------------------------------------------------
	if (color_vbo) {
		glBindBufferARB(GL_ARRAY_BUFFER_ARB, color_vbo);
		glColorPointer(4, GL_FLOAT, 0, 0);
		glEnableClientState(GL_COLOR_ARRAY);
	}

	//glBindVertexArray(g_sim->m_vao_pos);
	int start = 0;
	for (int i = 0; i < m_num_elasticbodies; i++) {
		int count = m_elasticbodies[i]->GetNumParticles();
		if (g_display_hair<0||i==g_display_hair) {
			//初期状態では全ての毛髪を表示，Nキーを一度押すと，ランダムに一本の毛髪のみを表示
			glDrawArrays(GL_LINE_STRIP, start, count);
		}
		start += count;
	}

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	//--------------------------------------------------

	glUseProgram(0);
	glDisable(GL_POINT_SPRITE_ARB);
}

void ScenePBD::Init(int argc, char* argv[])
{
	// GLEWの初期化
	GLenum err = glewInit();
	if (err != GLEW_OK) cout << "GLEW Error : " << glewGetErrorString(err) << endl;

	// マルチサンプル設定の確認
	//GLint buf, sbuf;
	//glGetIntegerv(GL_SAMPLE_BUFFERS, &buf);
	//cout << "number of sample buffers is " << buf << endl;
	//glGetIntegerv(GL_SAMPLES, &sbuf);
	//cout << "number of samples is " << sbuf << endl;
	glEnable(GL_MULTISAMPLE);

	glEnable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);

	glEnable(GL_POINT_SMOOTH);

	// 光源&材質描画設定
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHTING);
	glLightfv(GL_LIGHT0, GL_POSITION, glm::value_ptr(RX_LIGHT0_POS));
	glLightfv(GL_LIGHT0, GL_DIFFUSE, glm::value_ptr(RX_LIGHT_DIFF));
	glLightfv(GL_LIGHT0, GL_SPECULAR, glm::value_ptr(RX_LIGHT_SPEC));
	glLightfv(GL_LIGHT0, GL_AMBIENT, glm::value_ptr(RX_LIGHT_AMBI));

	glShadeModel(GL_SMOOTH);

	//海老沢追加---------------------------------------------------------
	//glMaterialfv(GL_FRONT, GL_AMBIENT, BrownAmbi);
	//glMaterialfv(GL_FRONT, GL_DIFFUSE, glm::value_ptr(BrownDiff));
	//glMaterialfv(GL_FRONT, GL_SPECULAR, glm::value_ptr(BrownSpec));
	//glMaterialf(GL_FRONT, GL_SHININESS, BrownShine);
	//--------------------------------------------------------------------

	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);

	// 視点初期化
	resetview();

	// 描画フラグ初期化
	m_draw = 0;
	m_draw |= RXD_BBOX;
	m_draw |= RXD_FLOOR;

	//利用するファイルの指定
	//char* filename = "Assets/1024-32/Curly(22).obj";
	//テスト用
	//char* filename = "AssetsNotUsed/Test_Case_hair_1024.obj";

	//Bob
	char* bob_filename = "Assets/LargeSize(1024-32)/Bob(24).obj";
	int bob_type = BOB_STYLE;
	//bob_filename = "AssetsNotUsed/Test_Case_hair.obj";//Kajiya-Kayモデルのスクショに利用
	//Curl
	char* curl_filename = "Assets/LargeSize(1024-32)/Curl.obj";
	int curl_type = CURL_STYLE;
	//Wavy
	char* wavy_filename = "Assets/LargeSize(1024-32)/Windy.obj";
	int wavy_type = WAVY_STYLE;
	//実験用
	char* experiment_filename = "AssetsTimeExperiment/Curl1024.obj";


	// PBD初期設定
	//initStraightRod();
	//initCenterSpiralRod();
	//initNaturalSpiralRod();
	//initExampleRod();
	initMoreRod(curl_filename, true, curl_type);//SagFreeがされない
	//initMoreRod(wavy_filename, true, wavy_type);//SagFreeがされる
}


/*!
 * 再描画イベント処理関数
 */
void ScenePBD::Draw(void)
{
	// ビューポート,透視変換行列,モデルビュー変換行列の設定
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glm::mat4 mp = glm::perspective(RX_FOV, (float)m_winw / m_winh, 0.2f, 1000.0f);
	glMultMatrixf(glm::value_ptr(mp));
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// 描画バッファのクリア
	glClearColor((GLfloat)m_bgcolor[0], (GLfloat)m_bgcolor[1], (GLfloat)m_bgcolor[2], 1.0f);
	glClearDepth(1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glPushMatrix();

	// マウスによる回転・平行移動の適用
	m_view.Apply();

	// シミュレーション空間情報
	glm::vec3 cen(0.0), dim(2.0);

	// 床描画
	glEnable(GL_LIGHTING);
	glm::vec3 lightpos(RX_LIGHT0_POS[0] * g_simulation_size_scene, RX_LIGHT0_POS[1] * g_simulation_size_scene, RX_LIGHT0_POS[2] * g_simulation_size_scene);//大きいサイズ感
	//glm::vec3 lightpos(RX_LIGHT0_POS[0],RX_LIGHT0_POS[1],RX_LIGHT0_POS[2]);//小さいサイズ感
	glm::vec3 lightcol(RX_LIGHT_DIFF[0], RX_LIGHT_DIFF[1], RX_LIGHT_DIFF[2]);
	if (m_draw & RXD_FLOOR) drawFloor(lightpos, lightcol, -1.5 * g_simulation_size_scene * dim[1], 30.0 * g_simulation_size_scene);//大きいサイズ感
	//if (m_draw & RXD_FLOOR) drawFloor(lightpos, lightcol, -1.5 * dim[1], 30.0);//小さいサイズ感

	// 境界,軸描画
	glDisable(GL_LIGHTING);
	if (m_draw & RXD_BBOX) drawWireCuboid(cen - 0.5f * g_simulation_size_scene * dim, cen + 0.5f * g_simulation_size_scene * dim, glm::vec3(0.5, 1.0, 0.5), 2.0);//大きいサイズ感
	//if (m_draw & RXD_BBOX) drawWireCuboid(cen - 0.5f * dim, cen + 0.5f * dim, glm::vec3(0.5, 1.0, 0.5), 2.0);//小さいサイズ感
	if (m_draw & RXD_AXIS) drawAxis(0.6 * dim[0], 3.0);

	// オブジェクト描画
	glEnable(GL_LIGHTING);
	glPushMatrix();

	//複数の弾性体
	/*for (int i = 0; i < m_elasticbodies.size(); i++) {
		if (m_elasticbodies[i]) m_elasticbodies[i]->Draw(m_draw);
	}*/
	DrawCollisionSphereVBO(g_sim->m_center.x, g_sim->m_center.y, g_sim->m_center.z, g_sim->m_rad);

	//SPH法側の描画処理-------------------------------------------------------------------------------------------
	if (g_sim) {
		int pnum = g_sim->GetNumParticles();	// 全粒子数
		unsigned int pvbo = g_sim->GetCurrentReadBuffer();			// 粒子座標VBO
		unsigned int cvbo = g_sim->GetColorBuffer();					// 粒子描画色VBO
		unsigned int tang_vbo = g_sim->GetTangBuffer(); //接線vbo
		float prad = (float)g_sim->m_params.particle_radius;				// 粒子半径
		float pscale = 1.2 * m_winh / tanf(RX_FOV * 0.5f * (float)RX_PI / 180.0f); // 描画スケール(PointSprite用)
		drawHairObject(pvbo, pnum, cvbo, tang_vbo, 0, 0, pscale, prad);
		//drawPointSprites(pvbo, pnum, cvbo, 0, 0, pscale, prad);
		//drawParticlePoints(pvbo, pnum, cvbo, 0, 0, glm::vec3(0, 0, 1), 4.0f);
	}
	//--------------------------------------------------------------------------------------------------------------------------------

	glPopMatrix();

	glPopMatrix();
}

/*!
 * タイマーコールバック関数
 */
void ScenePBD::Timer(void)
{
	if (m_animation_on) {
		// 描画を画像ファイルとして保存
		if (m_simg_spacing > 0 && m_currentstep % m_simg_spacing == 0) savedisplay(-1);

		//連番Objの保存
		if (g_obj_Animation_flag) {
			vector<glm::vec3> vrts(g_sim->GetNumParticles());
			vector<int> fix_array(g_sim->GetNumParticles());
			g_sim->GetArrayFromDevice(g_sim->P_POSITION, &vrts[0], g_sim->GetNumParticles());
			g_sim->GetArrayFromDevice(g_sim->P_FIX, &fix_array[0], g_sim->GetNumParticles());

			static int nsave = 1;
			string filename = CreateFileName("../AnimationObj/hair_", ".obj", nsave++, 5);
			SaveAnimationObj(filename, vrts, fix_array);
		}


		//時間計測---------------------------------------------------------
		//LARGE_INTEGER t1, t2;
		//LARGE_INTEGER f;
		//QueryPerformanceFrequency((LARGE_INTEGER*)&f);	// 高分解能パフォーマンスカウンタの周波数を取得

		//cudaDeviceSynchronize();
		//QueryPerformanceCounter((LARGE_INTEGER*)&t1);

		g_sim->Update(DT, m_currentstep);

		//cudaDeviceSynchronize();
		//QueryPerformanceCounter((LARGE_INTEGER*)&t2);

		//if (m_currentstep < 100)printf("%f \n", (double)(t2.QuadPart - t1.QuadPart)*1000 / (double)(f.QuadPart));//msecに合わせる
		//-----------------------------------------------------------------

		//ステップ数の出力(ちゃんと動いているかの確認
		cout << m_currentstep << " steps" << endl;

		//動画用のあるタイムステップで風を与える処理
		if(m_currentstep==399){
			g_sim->m_wind_flag = true;
			g_sim->m_wind_power = make_float3(g_experiment_power.x, g_experiment_power.y, g_experiment_power.z);
		}
		else if (m_currentstep == 799) {
			g_sim->m_wind_flag = false;
			g_sim->m_wind_power = make_float3(0.f,0.f,0.f);
		}

		//実験用の終了
		if (m_currentstep == 1000) {
			switchanimation(-1);
		}

		if (m_currentstep > MAX_STEPS) m_currentstep = 0;
		m_currentstep++;
	}
}

/*!
 * マウスイベント処理関数
 * @param[in] button マウスボタン(GLFW_MOUSE_BUTTON_LEFT,GLFW_MOUSE_BUTTON_MIDDLE,GLFW_MOUSE_BUTTON_RIGHT)
 * @param[in] action マウスボタンの状態(GLFW_PRESS, GLFW_RELEASE)
 * @param[in] mods 修飾キー(CTRL,SHIFT,ALT) -> https://www.glfw.org/docs/latest/group__mods.html
 */
void ScenePBD::Mouse(GLFWwindow* window, int button, int action, int mods)
{
	double x, y;
	glfwGetCursorPos(window, &x, &y);
	glm::vec2 mpos(x / (float)m_winw, (m_winh - y - 1.0) / (float)m_winh);

	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		if (action == GLFW_PRESS) {
			// マウスピック
			if (m_elasticbody && !(mods & 0x02)) {
				// マウスクリックした方向のレイ
				glm::vec3 ray_from, ray_to;
				glm::vec3 init_pos = glm::vec3(0.0);
				m_view.CalLocalPos(ray_from, init_pos);
				m_view.GetRayTo(x, y, RX_FOV, ray_to);
				glm::vec3 dir = glm::normalize(ray_to - ray_from);	// 視点からマウス位置へのベクトル

				// レイと各頂点(球)の交差判定
				float d;
				int v = m_elasticbody->IntersectRay(ray_from, dir, d);
				if (v >= 0) {
					//cout << "vertex " << v << " is selected - " << d << endl;
					m_picked = v;
					m_pickdist = d;
				}
			}

			if (m_picked == -1) {
				// マウスドラッグによる視点移動
				m_view.Start(x, y, mods);
			}

		}
		else if (action == GLFW_RELEASE) {
			m_view.Stop(x, y);
			if (m_picked != -1) {
				if (mods & 0x01) {	// SHIFTを押しながらボタンを離すことでその位置に固定
					m_elasticbody->FixVertex(m_picked);
				}
				else {
					m_elasticbody->UnFixVertex(m_picked);
					m_picked = -1;
				}
			}
		}
	}
}
/*!
 * モーションイベント処理関数(マウスボタンを押したままドラッグ)
 * @param[in] x,y マウス座標(スクリーン座標系)
 */
void ScenePBD::Cursor(GLFWwindow* window, double x, double y)
{
	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE &&
		glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_RELEASE &&
		glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE) {
		return;
	}

	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
		if (m_picked == -1) {
			m_view.Motion(x, y);
		}
		else {
			glm::vec3 ray_from, ray_to;
			glm::vec3 init_pos = glm::vec3(0.0);
			m_view.CalLocalPos(ray_from, init_pos);
			m_view.GetRayTo(x, y, RX_FOV, ray_to);
			glm::vec3 dir = glm::normalize(ray_to - ray_from);	// 視点からマウス位置へのベクトル

			glm::vec3 cur_pos = m_elasticbody->GetVertexPos(m_picked);
			glm::vec3 new_pos = ray_from + dir * static_cast<float>(m_pickdist);
			m_elasticbody->FixVertex(m_picked, new_pos);
		}
	}
}

/*!
 * キーボードイベント処理関数
 * @param[in] key キーの種類 -> https://www.glfw.org/docs/latest/group__keys.html
 * @param[in] mods 修飾キー(CTRL,SHIFT,ALT) -> https://www.glfw.org/docs/latest/group__mods.html
 */
void ScenePBD::Keyboard(GLFWwindow* window, int key, int mods)
{
	switch (key) {
	case GLFW_KEY_S: // シミュレーションスタート/ストップ
		switchanimation(-1);
		break;

	case GLFW_KEY_C:	// 固定点の全解除
		m_elasticbody->UnFixAllVertex();
		break;

	case GLFW_KEY_N:
		//ランダムに一本の毛髪を表示
		srand((unsigned)time(NULL));
		g_display_hair = rand() % m_num_elasticbodies;
		break;

	default:
		break;
	}
}


/*!
 * リサイズイベント処理関数
 * @param[in] w キャンバス幅(ピクセル数)
 * @param[in] h キャンバス高さ(ピクセル数)
 */
void ScenePBD::Resize(GLFWwindow* window, int w, int h)
{
	m_winw = w; m_winh = h;
	m_view.SetRegion(w, h);
	glViewport(0, 0, m_winw, m_winh);
}

/*!
 * ImGUIのウィジット配置
 */
void ScenePBD::ImGui(GLFWwindow* window)
{
	ImGui::Text("menu:");

	if (ImGui::Button("start/stop")) { switchanimation(-1); }
	if (ImGui::Button("run a step")) { m_animation_on = true; Timer(); m_animation_on = false; }
	if (ImGui::Button("reset viewpos")) { resetview(); }
	if (ImGui::Button("save screenshot")) { savedisplay(-1); }
	ImGui::Separator();
	ImGui::CheckboxFlags("vertex", &m_draw, RXD_VERTEX);
	ImGui::CheckboxFlags("edge", &m_draw, RXD_EDGE);
	ImGui::CheckboxFlags("face", &m_draw, RXD_FACE);
	ImGui::CheckboxFlags("tetrahedral(wire)", &m_draw, RXD_TETRA);
	ImGui::CheckboxFlags("tetrahedral(cs)", &m_draw, RXD_TETCS);
	ImGui::CheckboxFlags("aabb", &m_draw, RXD_BBOX);
	ImGui::CheckboxFlags("axis", &m_draw, RXD_AXIS);
	ImGui::CheckboxFlags("floor", &m_draw, RXD_FLOOR);
	ImGui::Separator();
	//実験的に毛髪を設定する際に利用
	//if (ImGui::Button("HairInteraction")) { initExampleRod(); }
	//ImGui::Separator();
	//各髪型のファイル定義
	char* BobHairFile = "Assets/LargeSize(1024-32)/Bob(24).obj";
	char* CurlHairFile = "Assets/LargeSize(1024-32)/Curl.obj";
	char* WavyHairFile = "Assets/LargeSize(1024-32)/Windy.obj";
	//各髪型のタイプ定義
	int bob_type = BOB_STYLE;
	int curl_type = CURL_STYLE;
	int wavy_type = WAVY_STYLE;
	//髪型を読み込みSagFree処理をする
	if (ImGui::Button("BobHair(SagFree)")) { initMoreRod(BobHairFile, true, bob_type); }
	//SagFree処理をしない
	if (ImGui::Button("BobHair")) { initMoreRod(BobHairFile, false, bob_type); }
	ImGui::Separator();
	if (ImGui::Button("CurlHair(SagFree)")) { initMoreRod(CurlHairFile, true, curl_type); }
	if (ImGui::Button("CurlHair")) { initMoreRod(CurlHairFile, false, curl_type); }
	ImGui::Separator();
	if (ImGui::Button("WavyHair(SagFree)")) { initMoreRod(WavyHairFile, true, wavy_type); }
	if (ImGui::Button("WavyHair")) { initMoreRod(WavyHairFile, false, wavy_type); }
	ImGui::Separator();
	//風の設定
	//刻み幅
	float step = 0.125;
	string step_str = "SetWindPower (step:" + to_string(step) + ")";
	ImGui::Text(step_str.c_str());
	if (ImGui::SliderFloat("", &g_sim->m_wind_power.x, 0.0f, 5.0f)) {
		//値を刻み幅で丸める
		g_sim->m_wind_power.x = roundf(g_sim->m_wind_power.x / step) * step;
	}
	if (ImGui::Checkbox("ApplyWind", &g_Apply_flag)) {
		if (g_Apply_flag) {
			g_sim->m_wind_flag = true;
		}
		else {
			g_sim->m_wind_flag = false;
		}
	}
	//画像を毎ステップ保存
	if (ImGui::Checkbox("SaveAllStep", &g_spacing_flag)) {
		if (g_spacing_flag) {
			m_simg_spacing = 1;
		}
		else {
			m_simg_spacing = -1;
		}
	}
	//Objを毎ステップ保存
	if (ImGui::Checkbox("SaveAnimationObj", &g_obj_Animation_flag)) {

	}

	ImGui::Separator();
	if (ImGui::Button("quit")) { glfwSetWindowShouldClose(window, GLFW_FALSE); }

}

void ScenePBD::Destroy()
{
	clearArray();

	if (m_elasticbody) delete m_elasticbody;
	//弾性体の配列を削除
	for (int i = 0; i < m_elasticbodies.size(); i++) {
		if (m_elasticbodies[i]) delete m_elasticbodies[i];
	}
}


//海老沢追加 
//複数弾性体の配列を確保
void ScenePBD::initArray(int num_particles) {
	//SPHとXPBDとのデータのやり取りを管理する配列の作成(弾性体の追加)
	for (int i = 0; i < num_particles; i++) {
		m_NewPos.push_back(glm::vec3(0.0));
		m_NewVelocity.push_back(glm::vec3(0.0));
	}
}

//海老沢追加
//複数弾性体の配列をクリア
void ScenePBD::clearArray(void) {
	m_NewPos.clear();
	m_NewVelocity.clear();
}



/*
 * アニメーションN/OFF
 * @param[in] on trueでON, falseでOFF
 */
bool ScenePBD::switchanimation(int on)
{
	m_animation_on = (on == -1) ? !m_animation_on : (on ? true : false);
	return m_animation_on;
}

/*!
 * 現在の画面描画を画像ファイルとして保存(連番)
 * @param[in] stp 現在のステップ数(ファイル名として使用)
 */
void ScenePBD::savedisplay(const int& stp)
{
	static int nsave = 1;
	string fn = CreateFileName("img_", ".bmp", (stp == -1 ? nsave++ : stp), 5);
	saveFrameBuffer(fn, m_winw, m_winh);
	std::cout << "saved the screen image to " << fn << std::endl;
}

/*!
 * 視点の初期化
 */
void ScenePBD::resetview(void)//視点、カメラの初期位置変更中
{
	//通常の設定
	double q[4] = { 1, 0, 0, 0 };
	//横からの視点に回転させる場合
	//double q[4] = { 1/sqrt(2), 0, 1 / sqrt(2), 0};

	m_view.SetQuaternion(q);
	m_view.SetScaling(-3.0*g_simulation_size_scene);//大きいサイズ感
	//m_view.SetScaling(-3.0);//小さいサイズ感
	m_view.SetTranslation(0.0, 0.0);//元(0.0,0.0)
}

/*!
 * 衝突処理関数
 * @param[in] p 現在の座標
 * @param[out] np 衝突後の座標
 * @param[in] v 速度
 * @param[in] obj オブジェクト番号
 */
void Collision(glm::vec3& p, glm::vec3& np, glm::vec3& v, int obj)
{
	// 床以外の物体との衝突処理をしたい場合にここに書く

	//	// 球体との衝突判定
	//	glm::vec3 rpos = p-g_v3Cen;
	//	float d = norm(rpos)-g_fRad;
	//	if(d < 0.0){
	//		glm::vec3 n = Unit(rpos);
	//		np = g_v3Cen+n*g_fRad;
	//	}
}

//objファイルの読み込み
bool ScenePBD::readObjFile(const char* filename,vector<glm::vec3> &PosArray,vector<glm::ivec2> &IndexArray,vector<int> &FixArray) {
	//file open
	ifstream f_in(filename);
	if (!f_in) {
		cout << "Error! cannot open file \"" << filename << "\"" << endl;
		f_in.close();
		return false;
	}

	int id_e = 0;
	glm::ivec2 prev_index(0);
	
	while (!f_in.eof()) {
		char buf[64];
		f_in.getline(buf, 64);

		//vertex
		if (buf[0] == 'v') {
			glm::vec3 pos;
			sscanf_s(buf, "v %f %f %f", &pos.x, &pos.y, &pos.z);
			//大きいサイズ感
			//pos *= 10.0f;
			PosArray.push_back(pos);
		}
		if (buf[0] == 'l') {
			glm::ivec2 index;
			sscanf_s(buf, "l %d %d", &index.x, &index.y);
			IndexArray.push_back(index - glm::ivec2(1));

			if (prev_index.y != index.x) {//新しい毛髪
				FixArray.push_back(id_e);
			}
			prev_index = index;
			id_e++;
		}
	}
}


//髪型用にObjファイルを編集
bool ScenePBD::MakeHairObjFile(const char* In_filename,const char* Out_filename) {
	//file open
	ifstream f_in(In_filename);
	ofstream f_out(Out_filename);
	if (!f_in) {
		cout << "Error! cannot open file \"" << In_filename << "\"" << endl;
		f_in.close();
		return false;
	}

	if (!f_out) {
		cout << "Error! cannot open file \"" << Out_filename << "\"" << endl;
		f_out.close();
		return false;
	}

	int id_e = 0;
	glm::ivec2 prev_index(0);

	while (!f_in.eof()) {
		char buf[64];
		f_in.getline(buf, 64);

		if (buf[0] == 'l') {
			glm::ivec2 index;
			sscanf_s(buf, "l %d %d", &index.x, &index.y);

			//連続したインデックスでないなら毛髪として表現不能
			if (index.y-index.x!=1) {
				continue;
			}
		}

		f_out << buf << endl;
	}
}

bool ScenePBD::SaveAnimationObj(string filename, const vector<glm::vec3>& vrts, const vector<int>& fix) {
	ofstream file;

	file.open(filename.c_str());
	if (!file || !file.is_open() || file.bad() || file.fail()) {
		cout << "rxOBJ::Save : Invalid file specified" << endl;
		return false;
	}

	//頂点数
	int nv = (int)vrts.size();
	for (int i = 0; i < nv; i++) {
		file << "v " << Vec3ToString(vrts[i]) << endl;
	}

	//インデックスの出力(頂点数と固定点の識別のサイズは同じことを前提)
	int idx = 1;
	for (int i = 0; i < nv-1; i++) {//最後の粒子での出力は不要
		if (fix[i + 1] == 0) {
			file << "l " << idx << " " << idx + 1 << endl;
		}
		else if(fix[i+1]==1){
			//特に処理はない
		}
		else {
			std::cerr << "fix data is incorrect!!" << endl;
		}
		idx++;
	}

	file.close();
	return true;
}

//SPHの初期化
void ScenePBD::initSPH(int max_particles,int num_elastic,float mass,float3 sphere_center,float sphere_rad, bool use_bp)
{
	// 描画スケール(PointSprite用)
	float view_scale = 1.2 * m_winh / tanf(RX_FOV * 0.5f * (float)RX_PI / 180.0f);

	// シミュレーションクラスのオブジェクト生成
	if (g_sim) delete g_sim;
	g_sim = new SPH();
	sph_scene.Clear();
	sph_scene.Set(g_sim);
	sph_scene.ReadSceneFiles();
	sph_scene_idx = sph_scene.GetCurrentSceneIdx();
	sph_scene.SetCurrentScene(sph_scene_idx);
	// シーン全体情報の読み込み
	if (sph_scene.LoadSpaceFromFile()) {
		SceneParameter env = sph_scene.GetEnv();

		//海老沢追加---------------------------------------------------------------
		env.max_particles = max_particles;//全体の粒子数の変更
		env.gravity.y = g_user_gravity.y;//海老沢追加 重力加速度を無理やり変更
		env.mass = mass;
		//cout << "mass " << mass << endl;
		env.dt = DT;
		cout << "DT " << DT << endl;
		//-------------------------------------------------------------------------

		g_sim->Initialize(env, num_elastic, sphere_center, sphere_rad);
		//g_dt = DT;
	}
	m_currentstep = 0;

	// シーンの個別設定(粒子，固体)の読み込み
	g_sim->Reset();
	cout << "load a scene from file." << endl;
	sph_scene.LoadSceneFromFile();
	if (use_bp) {
		g_sim->m_use_bparticles = true;
		g_sim->SetBoundaryParticles(true);
	}
	else {
		g_sim->m_use_bparticles = false;
	}
}

//海老沢追加
//配列で弾性体を扱う
//num_elasticbodies: 弾性体の数
//num_particles: 一つの弾性体当たりの粒子の数(全ての弾性体で粒子数は同じ)
//kss: 伸び剛性
//kbt: 曲げ剛性
void ScenePBD::initElasticbodies(int num_elasticbodies, int num_particles,float kss,float kbt,int type) {
	//速度と位置の配列を初期化
	clearArray();
	//ScenePBD側の変数に格納
	m_num_elasticbodies = num_elasticbodies;
	m_allnum = m_num_elasticbodies * num_particles;

	//既に値が格納されている場合には、これを初期化
	for (int i = 0; i < m_elasticbodies.size(); i++) {
		if (m_elasticbodies[i]) delete m_elasticbodies[i];
	}
	m_elasticbodies.clear();

	glm::vec3 env_min(-3.0, GOUND_HEIGHT - 10.0, -3.0);//海老沢変更中
	glm::vec3 env_max(3.0, GOUND_HEIGHT + 3.0, 3.0);

	srand((unsigned int)time(NULL));//球上に配置する際のためのシード値の設定
	//衝突させる球の設定(この球の上半分が毛髪の根元とする)
	glm::vec3 center(0.0, 0.0, 0.0);
	float rad = 0.35 * g_simulation_size_scene;
	//新しい弾性体を設定
	for (int i = 0; i < num_elasticbodies; i++) {
		m_elasticbodies.push_back(0);
		m_elasticbodies[i] = new ElasticPBD(0);

		m_elasticbodies[i]->SetCoefficient(kss, kbt);
		m_elasticbodies[i]->SetSimulationSpace(env_min, env_max);
		m_elasticbodies[i]->SetCollisionFunc(Collision);
		m_elasticbodies[i]->SetCollisionSphere(center, rad);

		//とりあえず、弾性体の配置は横一列に並べるとする。
		switch (type) {
		default:
		case STRAIGHT:
			//y軸
			//m_elasticbodies[i]->GenerateStrand(glm::vec3(0.015 * i, 0.50, 0.f), glm::vec3(0.015 * i, -0.50, 0.f), num_particles - 1);
			//x軸逆向き
			m_elasticbodies[i]->GenerateStrand(glm::vec3(0.0, 0.50, i * 0.015), glm::vec3(std::pow(-1.0,i), 0.50, i * 0.015), num_particles - 1);//変更中
			//z軸
			//m_elasticbodies[i]->GenerateStrand(glm::vec3(0.015 * i, 0.50, 0.0), glm::vec3(0.015 * i, 0.50, 1.0), num_particles - 1);
			break;
		case CENTER_SPIRAL:
			m_elasticbodies[i]->GenerateCenterSpiral(glm::vec3(i * 0.1, 0.50, 0.0), glm::vec3(i * 0.1, -0.50, 0.0), num_particles - 1);
			break;
		case NATURAL_SPIRAL:
			m_elasticbodies[i]->GenerateNaturalSpiral(glm::vec3(i * 0.1, 0.50, 0.0), glm::vec3(i * 0.1, -0.50, 0.0), num_particles - 1);
			break;
		case EXAMPLE:
			m_elasticbodies[i]->GenerateExampleRod(glm::vec3(i * 0.1, 0.50, 0.0), glm::vec3(i * 0.1, -0.50, 0.0), num_particles - 1);
			break;
		case ON_SPHERE:
			const float PI = 3.14159265359;
			//極座標からランダムな位置に配置
			float theta = rand() * (PI + 1.0) / (1.0 + RAND_MAX);
			float phi = rand() * (2.0 * PI + 1.0) / (1.0 + RAND_MAX);
			
			//ベクトル表現可能な形に
			float x = rad * std::cos(phi) * std::cos(theta);
			float y = rad * std::sin(phi) * std::cos(theta);
			if (y < 0)y *= -1;
			//if (y < rad - 0.2&&y*y<rad*rad-x*x) y += 0.1;
			//float z = rad * std::sin(theta);
			float z = rad * rad - x * x - y * y;

			if (x < 0)x *= -1;
			//if (z < 0)z *= -1;
			//弾性体の生成(長さは1,0固定としている)
			//cout << "Init Pos" << glm::to_string(glm::vec3(x, y, z)) << endl;
			//cout << "Length" << glm::length(glm::vec3(x, y, z));
			glm::vec3 start_pos(x, y, z);
			glm::vec3 d = start_pos - center;
			d = glm::normalize(d)*g_simulation_size_scene;
			m_elasticbodies[i]->GenerateStrand(glm::vec3(x, y, z), start_pos + d, num_particles - 1);
		}
		m_elasticbodies[i]->FixVertex(0);

		//SagFree処理の実行(CPU)
		//m_elasticbodies[i]->SagFree(kss, kbt);
	}
}

//海老沢追加
//Objファイルから弾性体の生成(それぞれの配列はobjからとってきたもの)
//PosArray:粒子位置の配列
//IndexArray:エッジを構成する2粒子のインデックスの配列
//FixArray:固定点に連結するエッジのインデックスの配列(このエッジから新しい弾性体が始まる)
//ks:伸び剛性
//kbt:曲げ剛性
//def_mass:質量
//num_elastic:弾性体の数を返す
//all_particle:合計の粒子数を返す
void ScenePBD::initElasticFromObj(vector<glm::vec3> PosArray, vector<glm::ivec2> IndexArray, vector<int> FixArray,float ks,float kbt,float def_mass,int& num_elastic,int& all_particle,int type) {
	clearArray();
	int num_elasticbodies = FixArray.size();
	num_elastic = num_elasticbodies;
	all_particle = PosArray.size();
	m_num_elasticbodies = num_elasticbodies;
	m_allnum = all_particle;
	//既に値が格納されている場合には、これを初期化
	for (int i = 0; i < m_elasticbodies.size(); i++) {
		if (m_elasticbodies[i]) delete m_elasticbodies[i];
	}
	m_elasticbodies.clear();

	glm::vec3 env_min(-3.0, GOUND_HEIGHT - 10.0, -3.0);//海老沢変更中
	glm::vec3 env_max(3.0, GOUND_HEIGHT + 3.0, 3.0);

	glm::vec3 center(0.0, 0.0, 0.0);
	float rad = 0.25;
	float mass = def_mass;
	
	int fix_index = 0;
	for (int i = 0; i < num_elasticbodies; i++) {
		m_elasticbodies.push_back(0);
		m_elasticbodies[i] = new ElasticPBD(0);

		m_elasticbodies[i]->SetCoefficient(ks, kbt);
		m_elasticbodies[i]->SetSimulationSpace(env_min, env_max);
		m_elasticbodies[i]->SetCollisionFunc(Collision);
		m_elasticbodies[i]->SetCollisionSphere(center, rad);

		//m_elasticbodies[i]->Clear();//頂点の削除
		int iter = 0;//エッジの追加の回数(下の反復の回数)
		int last;
		if (fix_index == num_elasticbodies - 1) last = IndexArray.size();
		else last = FixArray[fix_index + 1];

		for (int j = FixArray[fix_index]; j < last; j++) {
			int index_pos0 = IndexArray[j].x;
			int index_pos1 = IndexArray[j].y;

			glm::vec3 pos0 = PosArray[index_pos0];
			glm::vec3 pos1 = PosArray[index_pos1];

			if (iter == 0)m_elasticbodies[i]->AddVertex(pos0, mass);

			m_elasticbodies[i]->AddVertex(pos1, mass);
			m_elasticbodies[i]->AddEdge(iter, iter + 1);
			if (iter > 0)m_elasticbodies[i]->AddDarbouxVector(iter);
			iter++;
		}
		//両端の点を比較し，y座標が大きい方を固定点とする．(基本は[0]が固定点とする)
		if (m_elasticbodies[i]->GetVertexPos(0).y < m_elasticbodies[i]->GetVertexPos(iter).y) {//最後の頂点をiterでとれているかは少し不安
			m_elasticbodies[i]->ArrayReverse();
		}

		fix_index++;
		m_elasticbodies[i]->FixVertex(0);
		
		
		//SagFree処理(CPU)
		//m_elasticbodies[i]->SagFree(ks, kbt);
	}

	cout << "num_particles " << m_allnum << endl;
	cout << "num_edges " << IndexArray.size() << endl;
	cout << "num_elasticbodies" << num_elastic << endl;
}

//海老沢追加
//XPBD(CPU)からSPH(GPU)にデータを転送するための配列にデータを格納
//num_elasticbodies: 弾性体の数
//num_particles: 一つの弾性体当たりの粒子の数(全ての弾性体で粒子数は同じ)
void ScenePBD::XPBDtoSPH(void) {
	vector<glm::vec3> NewPos, NewVel,NewTang;
	for (int i = 0; i < m_num_elasticbodies; i++) {
		int particles = m_elasticbodies[i]->GetNumParticles();
		//cout << "Num Particle" << particles << endl;
		for(int j=0;j<particles;j++){
			NewVel.push_back(m_elasticbodies[i]->GetVertexVelPointer()[j]);//m_num_particlesからparticlesに変更
			NewPos.push_back(m_elasticbodies[i]->GetVertexPosPointer()[j]);
			//接線を保存------------------------------------------------------------------------
			if (j == 0) NewTang.push_back(glm::normalize(m_elasticbodies[i]->GetTangPointer()[j]));//正規化をして代入
			else NewTang.push_back(glm::normalize(m_elasticbodies[i]->GetTangPointer()[j - 1]));
		}
	}

	//データの転送
	if (g_sim) g_sim->Set(NewPos, NewVel,NewTang, m_allnum);
}

//海老沢追加
//SPH(GPU)からXPBD(CPU)にデータを転送する
//num_elasticbodies: 弾性体の数
//num_particles: 一つの弾性体当たりの粒子の数(全ての弾性体で粒子数は同じ)
void ScenePBD::SPHtoXPBD(void) {
	vector<float> acc(m_allnum * DIM), pos(m_allnum * DIM);
	if (g_sim) {
		//加速度を渡す
		g_sim->GetArrayFromDevice(g_sim->P_ACC, &acc[0], m_allnum);
		//位置を渡す
		g_sim->GetArrayFromDevice(g_sim->P_POSITION, &pos[0], m_allnum);
	}
	int start = 0;
	for (int i = 0; i < m_num_elasticbodies; i++) {
		//加速度設定
		m_elasticbodies[i]->SetVertexAcc(acc, start, m_elasticbodies[i]->GetNumParticles());
		//位置設定(CPUで時間積分をする場合)
		//m_elasticbodies[i]->SetVertexPos(pos, start, m_elasticbodies[i]->GetNumParticles());
		//GPUで時間積分をする場合
		m_elasticbodies[i]->SetVertexCurPos(pos, start, m_elasticbodies[i]->GetNumParticles());
		start += m_elasticbodies[i]->GetNumParticles();
	}
}

//海老沢追加
//デバック用に弾性体の配列を削除
void ScenePBD::DestroyElasticArray() {
	for (int i = 0; i < m_elasticbodies.size(); i++) {
		if (m_elasticbodies[i]) delete m_elasticbodies[i];
	}
	m_elasticbodies.clear();
}

//海老沢追加
//デバイスにXPBDのパラメータを渡す処理をまとめる
void ScenePBD::XPBDParamsToDevice(int num_elasticbodies) {
	vector<float> Mass_array, Length_array, Kss_array, Kbt_array;
	vector<glm::quat> Quat_array, Darboux_array;
	vector<int> Fix_array,Last_index;
	for (int i = 0; i < num_elasticbodies; i++) {
		int particles = m_elasticbodies[i]->GetNumParticles();
		for (int j = 0; j < particles - 1; j++) {//num_particlesからparticlesに変更
			Mass_array.push_back(m_elasticbodies[i]->GetVertexMass()[j]);
			Length_array.push_back(m_elasticbodies[i]->GetInitial_Length()[j]);
			Kss_array.push_back(m_elasticbodies[i]->GetInitial_Kss()[j]);
			//cout << "Kss" << m_elasticbodies[i]->GetInitial_Kss()[j] << endl;
			Kbt_array.push_back(m_elasticbodies[i]->GetInitial_Kbt()[j]);
			Quat_array.push_back(m_elasticbodies[i]->GetQuat()[j]);
			if (j != particles - 2) Darboux_array.push_back(m_elasticbodies[i]->GetRestDarboux()[j]);//num_particlesからparticlesに変更
			else Darboux_array.push_back(glm::quat(0.f, glm::vec3(0.0)));
			if (m_elasticbodies[i]->GetVertexFix()[j]) Fix_array.push_back(1);
			else Fix_array.push_back(0);

			/*if (j == 0) {
				cout << "quat[0]:" << m_elasticbodies[i]->GetQuat()[j][0] << " quat[1]:" << m_elasticbodies[i]->GetQuat()[j][1] << " quat[2]:" << m_elasticbodies[i]->GetQuat()[j][2] << " quat.w" << m_elasticbodies[i]->GetQuat()[j][3] << endl;
				cout << "quat.x" << m_elasticbodies[i]->GetQuat()[j].x << " quat.y:" << m_elasticbodies[i]->GetQuat()[j].y << " quat.z:" << m_elasticbodies[i]->GetQuat()[j].z << " quat.w" << m_elasticbodies[i]->GetQuat()[j].w << endl;
			}*/
		}
		Mass_array.push_back(m_elasticbodies[i]->GetVertexMass()[particles-1]);
		Length_array.push_back(0.f);
		Kss_array.push_back(0.f);
		Kbt_array.push_back(0.f);
		Quat_array.push_back(glm::quat(0.f, glm::vec3(0.0)));
		Darboux_array.push_back(glm::quat(0.f, glm::vec3(0.0)));
		if (m_elasticbodies[i]->GetVertexFix()[particles - 1]) Fix_array.push_back(1);//num_particlesからparticleに変更
		else Fix_array.push_back(0);

		Last_index.push_back(Fix_array.size() - 1);
	}

	g_sim->SetXPBD_Params(m_allnum, &Mass_array[0], &Length_array[0], &Kss_array[0], &Kbt_array[0], &Quat_array[0], &Darboux_array[0], &Fix_array[0],&Last_index[0]);
	//SPHの質量設定
	//g_sim->m_params.mass = Mass_array[0];
}

void ScenePBD::initMoreRod(char* filename,bool sag_free_flag ,int type) {
	switchanimation(false);

	m_draw |= RXD_VERTEX;
	m_draw |= RXD_EDGE;

	m_currentstep = 0;

	//後からobjファイル内で定義するため，この初期値に意味はない
	int num_elasticbodies = 10;//弾性体の数
	int num_particles = 11;//弾性体一つ当たりの粒子数
	int all_particles = 0;

	//弾性体の初期化
	//initElasticbodies(num_elasticbodies, num_particles, ks, kbt, ON_SPHERE);//ON_SPHEREで球上に配置

	//髪型用にObjファイルの作成
	//髪型アセットの作成----------------------------------------------------------------------------------------------------------------------------
	//Objファイルの読み込み
	// 
	// 
	// 
	// objファイルを作成する場合



	//rxOBJ obj;
	//vector <glm::vec3> vrts,vnms;//頂点情報,法線情報
	//vector <rxFace> plys;//面情報
	//rxMTL mats;//材質情報
	//obj.Read("AssetsNotUsed/Curl.obj", vrts, vnms, plys, mats);

	////頂点列をクランプする
	//FitVertices(vrts);

	////ポリゴンへの設定
	//rxPolygonsE poly;
	//poly.vertices = vrts;
	//poly.faces = plys;

	////髪型アセットの作成
	//int num = 32;//目標は32
	////閾値T
	////小さいサイズ感
	////float T = 0.25;
	////大きいサイズ感
	//float T = 0.25 * g_simulation_size_scene;
	//MakeHairObj("Test_Case", poly, T , num);



	//エッジの探索
	/*int edge_count=SearchEdge(poly);
	cout << "edge_count" << edge_count << endl;
	
	SearchBoundaryEdge(poly);*/


	//-----------------------------------------------------------------------------------------------------------------------------------------------



	//objファイルから弾性体の生成
	vector<glm::vec3> PosArray;
	vector<glm::ivec2> IndexArray;
	vector<int> FixArray;

	//Objファイルを読む
	readObjFile(filename, PosArray, IndexArray, FixArray);
	//test.objの場合、回転をさせる
	if (strcmp(filename, "Assets/3000.obj")) {
		for (int i = 0; i < PosArray.size(); i++) {
			float tmp = PosArray[i].x;
			PosArray[i].x = -PosArray[i].z;
			PosArray[i].z = -tmp;
		}
	}
	//パラメータ設定
	float mass, ks, kbt;
	//質量
	mass = 1.0e-1;//5.0e-3
	//大きいサイズ感により，変更
	ks = 20.f * g_simulation_size_scene;//伸び剛性
	kbt = 5.f * g_simulation_size_scene;//曲げ剛性
	//球のサイズ
	float3 sphere_center = make_float3(-5.0 * g_simulation_size_scene, 0.0, 0.0);
	float sphere_rad = 0.35 * g_simulation_size_scene;
	//LocalTorqueStepの調整パラメータ
	float K_min = 0.005f;
	//ボブヘア
	if (type == BOB_STYLE) {
		//球の位置やサイズ，質量，剛性以外にも変更すべきパラメータがあった場合には，typeを受けるg_sim->hair_typeなどを作って渡し，sph.cppのUpdate関数で参照できるように
		//球の中心と半径
		sphere_center = make_float3(0.0, 0.175 * g_simulation_size_scene, 0.0);
		sphere_rad = 0.4 * g_simulation_size_scene;
		//剛性
		ks = 20.f;//20.f
		kbt = 75.f;//7.5f 一時的に変更
		//質量
		mass = 1.0e-1;
		//LocalTorqueStepの調整パラメータ
		//卒論用
		//K_min = 0.01f;//0.0075f
		//中間発表用の動画
		K_min = 0.0075f;

		//風の強さ 0.75(等しい実験)
		g_experiment_power = glm::vec3(1.0f, 0.f, 0.f);//0.75(卒論) 最終発表(1.0f)
		//SagFree下で無理やり通常時と同じくらい動かしたい場合(卒論)
		//g_experiment_power = glm::vec3(3.25f, 0.f, 0.f);//2.75f(卒論) 最終発表(3.25f)
	}
	//カールヘア
	else if (type == CURL_STYLE) {
		//球の中心と半径
		sphere_center = make_float3(0.0, 0.55 * g_simulation_size_scene, -0.2 * g_simulation_size_scene);
		sphere_rad = 0.33 * g_simulation_size_scene;
		//剛性
		ks = 120.f;
		kbt = 30.f;
		//質量
		mass = 1.0e-1;
		//LocalTorqueStepの調整パラメータ(卒論)
		//K_min = 0.005f;
		//中間発表用の動画
		K_min = 0.0035f;

		//風の強さ 1.25
		//g_experiment_power = glm::vec3(1.85f, 0.f, 0.f);//卒論-最終発表(1.85f)
	}
	//ウェーブヘア
	else if (type == WAVY_STYLE) {
		//球の中心と半径
		sphere_center = make_float3(-0.05 * g_simulation_size_scene, 0.55 * g_simulation_size_scene, 0.35 * g_simulation_size_scene);
		sphere_rad = 0.325 * g_simulation_size_scene;
		//剛性
		ks = 150.f;
		kbt = 50.f;
		//質量
		mass = 1.0e-1;
		//LocalTorqueStepの調整パラメータ(卒論)
		K_min = 0.00175f;//0.0025
		//中間発表用の動画

		//風の強さ 1.25
		//g_experiment_power = glm::vec3(1.75f, 0.f, 0.f);
	}
	else {
		std::cerr << "No HairStyle defined!!" << endl;
	}

	initElasticFromObj(PosArray, IndexArray, FixArray, ks, kbt, mass, num_elasticbodies, all_particles, type);
	
	m_picked = -1;
	//SPH法の初期化
	initSPH(all_particles, num_elasticbodies, mass, sphere_center, sphere_rad);//元はm_allnum(initElasticbodiesで設定)

	//フラグをオンにしてみる
	//g_sim->m_example_flag = true;

	//データ転送に用いる配列の初期化
	//initArray(all_particles);//元はm_allnum(initElasticbodiesで設定)

	//XPBDからSPH(CPU->GPU)に位置、速度のデータを移す
	XPBDtoSPH();

	//XPBDのパラメータ
	XPBDParamsToDevice(num_elasticbodies);

	//時間計測---------------------------------------------------------
	LARGE_INTEGER t1, t2;
	LARGE_INTEGER f;
	QueryPerformanceFrequency((LARGE_INTEGER*)&f);	// 高分解能パフォーマンスカウンタの周波数を取得

	cudaDeviceSynchronize();
	QueryPerformanceCounter((LARGE_INTEGER*)&t1);

	//SagFree処理(GPU)
	if (sag_free_flag) g_sim->SagFree(g_user_gravity, K_min);

	cudaDeviceSynchronize();
	QueryPerformanceCounter((LARGE_INTEGER*)&t2);

	printf("\n %f [msec] \n \n", (double)(t2.QuadPart - t1.QuadPart) * 1000 / (double)(f.QuadPart));//msecに合わせる
	printf("QueryPerformanceFrequency %f \n", (double)(f.QuadPart));
	//-----------------------------------------------------------------

	//ファイルにGPUに移したパラメータ出力
	//g_sim->OutputParticles("debug.txt");
}


/*!
 * ストランドで初期化
 */
void ScenePBD::initStraightRod(void)
{
	//シミュレーションをストップ
	switchanimation(false);

	m_draw |= RXD_VERTEX;
	m_draw |= RXD_EDGE;

	m_currentstep = 0;
	float ks = 10.f;
	float kbt = 2.0f;

	int num_elasticbodies = 2;
	int num_particles = 11;
	//int all_particles = num_elasticbodies * num_particles;

	initElasticbodies(num_elasticbodies, num_particles, ks, kbt, STRAIGHT);

	m_picked = -1;
	float mass = 1.0e-1;

	//SPH法の設定
	//initSPH(m_allnum,num_elasticbodies,mass);
	initArray(m_allnum);
	
	XPBDtoSPH();
	XPBDParamsToDevice(num_elasticbodies);

	//デバイスメモリのパラメータを出力
	//g_sim->OutputParticles("debug.txt");
}

void ScenePBD::initCenterSpiralRod(void)
{
	//シミュレーションをストップ
	switchanimation(false);

	m_draw |= RXD_VERTEX;
	m_draw |= RXD_EDGE;

	m_currentstep = 0;
	
	//剛性の設定
	//螺旋型のシーン用
	float ks = 120.f;
	float kbt = 120.f;//元100000.f

	int num_elasticbodies = 2;
	int num_particles = 18;

	initElasticbodies(num_elasticbodies, num_particles, ks, kbt, CENTER_SPIRAL);

	m_picked = -1;

	float mass = 1.0e-1;
	//initSPH(m_allnum,num_elasticbodies,mass);
	initArray(m_allnum);

	XPBDtoSPH();
	XPBDParamsToDevice(num_elasticbodies);
}

/*!
 * 四面体メッシュで構成されたボール形状で初期化
 */
void ScenePBD::initNaturalSpiralRod(void)
{
	//シミュレーションをストップ
	switchanimation(false);

	m_draw |= RXD_VERTEX;
	m_draw |= RXD_EDGE;

	m_currentstep = 0;
	glm::vec3 env_min(-3.0, GOUND_HEIGHT - 10.0, -3.0);//海老沢変更中
	glm::vec3 env_max(3.0, GOUND_HEIGHT + 3.0, 3.0);

	//剛性の設定
	//螺旋型のシーン用
	float ks = 120.f;
	float kbt = 200.f;//元100000.f
	
	int num_elasticbodies = 2;
	int num_particles = 22;

	initElasticbodies(num_elasticbodies, num_particles, ks, kbt, NATURAL_SPIRAL);

	m_picked = -1;
	float mass = 1.0e-1;

	//initSPH(m_allnum,num_elasticbodies,mass);
	initArray(m_allnum);

	XPBDtoSPH();
	XPBDParamsToDevice(num_elasticbodies);
}

void ScenePBD::initExampleRod(void) {
	//シミュレーションをストップ
	switchanimation(false);

	m_draw |= RXD_VERTEX;
	m_draw |= RXD_EDGE;

	m_currentstep = 0;
	
	//剛性の設定
	//実験用のシーン用
	float ks = 20.f;
	float kbt = 20.f;

	//球のサイズ
	float3 sphere_center = make_float3(0.0, 0.0, 0.0);
	float sphere_rad = 0.35 * g_simulation_size_scene;

	//パラメータ設定
	//Shear&Stretching Constraint  w = 1.0 / (length); wq = 1.0 / (1.0e-5) * length;
	//球 center=(0.0,0.0,0.0) rad=0.35
	//SPH rest_dens=5.0e2
	//g=-9.81

	int num_elasticbodies = 400;
	int num_particles = 18;

	initElasticbodies(num_elasticbodies, num_particles, ks, kbt, ON_SPHERE);

	m_picked = -1;

	float mass = 1.0e-1;
	initSPH(m_allnum, num_elasticbodies, mass, sphere_center, sphere_rad);
	initArray(m_allnum);

	//SPHの結果を表示する際には，他と少し変える----
	//g_sim->m_center = make_float3(0.0, 0.0, 0.0);
	//g_sim->m_rad = 0.35;
	g_sim->m_example_flag = true;
	//---------------------------------------------

	XPBDtoSPH();
	XPBDParamsToDevice(num_elasticbodies);
}

//SPHの方に風を与える
void ScenePBD::ChangeWindPower(float3 wind) {
	g_sim->m_wind_power = wind;
	return;
}

/*!
 * マウスピックの解除
 */
void ScenePBD::clearPick(void)
{
	if (m_elasticbody && m_picked != -1) m_elasticbody->UnFixVertex(m_picked);
	m_picked = -1;
}

//座標群を(0,1)にフィットさせる
void ScenePBD::FitVertices(vector<glm::vec3> &vertices) {
	int size = vertices.size();
	float x_num, y_num, z_num;
	x_num = y_num = z_num = 0;
	float x_max, x_min, y_max, y_min, z_max, z_min;
	//平均値のための計算
	for (int i = 0; i < size; i++) {
		glm::vec3 v = vertices[i];
		x_num += v.x;
		y_num += v.y;
		z_num += v.z;
	}

	float x_ave = x_num / size;
	float y_ave = y_num / size;
	float z_ave = z_num / size;

	//平均値を原点に持ってくる
	for (int i = 0; i < size; i++) {
		//全ての点がほぼ同じ位置にある場合などは想定していない
		glm::vec3 v = vertices[i];
		v -= glm::vec3(x_ave, y_ave, z_ave);

		if (i == 0) {
			x_max = x_min = v.x;
			y_max = y_min = v.y;
			z_max = z_min = v.z;
		}
		else {
			//最大値の更新
			if (x_max < v.x) x_max = v.x;
			if (y_max < v.y) y_max = v.y;
			if (z_max < v.z) z_max = v.z;

			//最小値の更新
			if (v.x < x_min) x_min = v.x;
			if (v.y < y_min) y_min = v.y;
			if (v.z < z_min) z_min = v.z;
		}
		vertices[i] = v;
	}

	float range_x = abs(x_max - x_min);//念のため絶対値に
	float range_y = abs(y_max - y_min);
	float range_z = abs(z_max - z_min);

	float max_range = max(max(range_x, range_y), range_z);

	for (int i = 0; i < size; i++) {
		glm::vec3 v = vertices[i];
		float x_new = (v.x - x_min) / max_range;
		float y_new = (v.y - y_min) / max_range;
		float z_new = (v.z - z_min) / max_range;

		//形状によって少し変える必要場ある．
		//アセット用
		glm::vec3 scale(1.35, 1.f, 1.35);
		//CurlHair
		vertices[i] = glm::vec3((x_new - 0.1) * scale[0], (y_new - 0.5) * scale[1], (z_new - 0.25) * scale[2]) * 2.0f;//x:-0.1 z:-0.25
		//WavyHair
		//vertices[i] = glm::vec3((x_new - 0.3) * scale[0], (y_new - 0.5) * scale[1], (z_new - 0.3) * scale[2]) * 2.0f;//x:-0.1 z:-0.25
		//頂いたObjファイル用
		//vertices[i] = glm::vec3(x_new - 0.5, y_new - 0.5, z_new - 0.4) * 1.5f;

		//大きいサイズ感
		vertices[i] *= g_simulation_size_scene;
	}

}

//計算点の数を制限する



