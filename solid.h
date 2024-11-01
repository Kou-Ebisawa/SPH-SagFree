/*!
  @file rx_sph_solid.h
	
  @brief SPH用固体定義
 
  @author Makoto Fujisawa
  @date 2008-12
*/

#ifndef _SPH_SOLID_H_
#define _SPH_SOLID_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "utils.h"
#include "gldraw.h"


//-----------------------------------------------------------------------------
// 定義・定数
//-----------------------------------------------------------------------------

const int RX_DISPMAP_N = 128;


//-----------------------------------------------------------------------------
// 衝突点情報記述用クラス
//-----------------------------------------------------------------------------
class rxCollisionInfo
{
protected:
	glm::vec3 m_contact_point;	//!< 衝突点
	glm::vec3 m_normal;			//!< (衝突点での)法線
	float m_penetration;		//!< めり込み量

	glm::vec3 m_velocity;		//!< 衝突点の速度

public:
	//! デフォルトコンストラクタ
	rxCollisionInfo()
	  : m_contact_point(glm::vec3(0.0)), 
		m_normal(glm::vec3(0.0)), 
		m_penetration(0.0), 
		m_velocity(glm::vec3(0.0))
	{
	}

	//! コンストラクタ
	rxCollisionInfo(const glm::vec3 &contact_point, 
					const glm::vec3 &normal = glm::vec3(0.0), 
					const float &penetration = 0.0, 
					const glm::vec3 &veloc = glm::vec3(0.0))
	  : m_contact_point(contact_point), 
		m_normal(normal), 
		m_penetration(penetration), 
		m_velocity(veloc)
	{
	}

	//! デストラクタ
	~rxCollisionInfo(){}

public:
	const glm::vec3& Contact() const { return m_contact_point; }
	glm::vec3& Contact(){ return m_contact_point; }

	const glm::vec3& Normal() const { return m_normal; }
	glm::vec3& Normal(){ return m_normal; }

	const float& Penetration() const { return m_penetration; }
	float& Penetration(){ return m_penetration; }

	const glm::vec3& Velocity() const { return m_velocity; }
	glm::vec3& Velocity(){ return m_velocity; }
};




//-----------------------------------------------------------------------------
// rxSolid : 固体オブジェクト基底クラス
//-----------------------------------------------------------------------------
class rxSolid
{
public:
	// 固体種別識別用
	enum{
		RXS_SPHERE, 
		RXS_AABB, 
		RXS_OPEN_BOX, 
		RXS_POLYGON, 
		RXS_IMPLICIT, 
		RXS_OTHER = -1, 
	};

protected:
	glm::vec3 m_mass_center;			//!< 重心座標
	glm::vec3 m_velocity;				//!< 固体速度
	glm::quat m_quat;					//!< 姿勢を表す四元数

	int m_name;							//!< 固体名(RXS_SPHERE,RXS_AABBなど)
	float m_offset;						//!< 固体表面陰関数のオフセット(境界粒子配置時に設定)
	
	bool m_fix;							//!< 固定フラグ
	bool m_gen_bparticles;				//!< 境界粒子生成フラグ

	int m_sgn;							//!< 箱:-1, オブジェクト:1

public:
	rxSolid() : m_name(RXS_OTHER)
	{
		m_fix = true;
		m_mass_center = glm::vec3(0.0);
		m_velocity = glm::vec3(0.0);
		m_sgn = 1;
		m_offset = (float)(0.0);
		m_gen_bparticles = true;

		m_quat = glm::quat(1.0f, 0.0f, 0.0f, 0.0f);
	}

	//
	// 仮想関数
	//
	/*!
	 * 距離値/曲率値計算
	 * @param[in] pos,pos0,pos1 グローバル座標での位置(pos0,pos1は前ステップと現ステップでの位置)
	 * @param[in] r 球の半径
	 * @param[out] col 距離などの情報(衝突情報)
	 */
	virtual bool GetDistance(const glm::vec3 &pos, rxCollisionInfo &col) = 0;	//!< 距離関数計算
	virtual bool GetDistanceR(const glm::vec3 &pos, const float &r, rxCollisionInfo &col) = 0;
	virtual bool GetDistance(const glm::vec3 &pos0, const glm::vec3 &pos1, rxCollisionInfo &col) = 0;	//!< 距離関数計算
	virtual bool GetDistanceR(const glm::vec3 &pos0, const glm::vec3 &pos1, const float &r, rxCollisionInfo &col) = 0;
	virtual bool GetCurvature(const glm::vec3 &pos, float &k) = 0;	//!< 距離関数の曲率計算

	/*!
	 * OpenGLによる描画
	 * @param[in] drw 描画フラグ(drw&2 == 1でワイヤフレーム描画)
	 */
	virtual void Draw(int drw) = 0;				//!< OpenGLでの描画

	/*!
	 * 固体の大きさ(領域)取得/設定関数
	 */
	virtual glm::vec3 GetMin(void) = 0;
	virtual glm::vec3 GetMax(void) = 0;
	virtual void SetMin(const glm::vec3& p){}
	virtual void SetMax(const glm::vec3& p){}

	//!< OpenGL変換行列の適用
	virtual void SetGLMatrix(void)
	{
		glTranslatef(m_mass_center[0], m_mass_center[1], m_mass_center[2]); 
		glMultMatrixf(glm::value_ptr(glm::mat4_cast(m_quat)));
	}

	//
	// 陰関数値(衝突判定,境界粒子生成に使う)
	//
	/*!
	 * 陰関数値とその勾配の計算
	 * @param[in] x,y,z 計算位置
	 * @return 勾配と陰関数値を格納した4次元ベクトル(0〜2:勾配,3:値)
	 */
	static glm::vec4 GetImplicitG_s(float x, float y, float z, void* ptr){ return ((rxSolid*)ptr)->GetImplicitG(glm::vec3(x, y, z)); }
	inline glm::vec4 GetImplicitG(glm::vec3 pos)
	{ 
		rxCollisionInfo col; 
		GetDistance(pos, col); 
		return glm::vec4(col.Penetration()+m_offset, col.Normal()[0], col.Normal()[1], col.Normal()[2]);
	}

	/*!
	 * 陰関数値の計算
	 * @param[in] x,y,z 計算位置
	 * @return 陰関数値
	 */
	static float GetImplicit_s(void* ptr, float x, float y, float z){ return ((rxSolid*)ptr)->GetImplicit(glm::vec3(x, y, z)); }
	inline float GetImplicit(glm::vec3 pos){ rxCollisionInfo col; GetDistance(pos, col); return col.Penetration()+m_offset; }
	static bool GetDistance_s(glm::vec3 pos, rxCollisionInfo &col, void* x){ return ((rxSolid*)x)->GetDistance(pos, col); }

	inline void SetOffset(float offset){ m_offset = offset; }

	//
	// 取得・設定関数
	//
	//! グローバルから固体ローカルへの座標変換
	inline glm::vec3 CalLocalCoord(const glm::vec3& pos){ return glm::conjugate(m_quat)*(pos-m_mass_center); }
	//! 固体ローカルからグローバルへの座標変換
	inline glm::vec3 CalGlobalCoord(const glm::vec3 &pos){ return m_quat*pos+m_mass_center; }
	//! 固体の回転を任意のベクトルに適用
	inline glm::vec3 ApplyRot(const glm::vec3& v){ return m_quat*v; }

	inline glm::vec3 GetPosition(void){	return m_mass_center; }				//!< 固体重心位置の取得
	inline void SetPosition(const glm::vec3& pos){ m_mass_center = pos; }	//!< 固体重心位置の設定

	inline glm::quat GetRotation(void){ return m_quat; }					//!< 回転を表す四元数の取得
	inline void SetRotation(const glm::quat& q){ m_quat = q; }				//!< 回転を表す四元数の設定

	inline glm::vec3 GetVelocityAtGrobal(const glm::vec3& pos){ return m_velocity; }	//!< 指定した座標値の固体速度の取得
	inline void SetVelocity(const glm::vec3& vec){ m_velocity = vec; }		//!< 固体中心速度の設定

	inline bool GetFix(void) const { return m_fix; }			//!< 固定フラグの取得
	inline void SetFix(bool fix){ m_fix = fix; }				//!< 固定フラグの設定

	inline bool& IsSolidParticles(void){ return m_gen_bparticles; }	//!< 固体粒子フラグの取得

	inline bool RigidSimulation(const float &dt)				//!< 剛体シミュレーション
	{
		m_mass_center += dt*m_velocity;
		return true;
	}

	inline const int& Name() const { return m_name; }
	inline int& Name(){ return m_name; }

protected:
	inline bool calCurvature(const glm::vec3& pos, float& k, bool (dstfunc)(glm::vec3, rxCollisionInfo&, void*), void* fp = 0);
};



//-----------------------------------------------------------------------------
// MARK:rxSolidBox : 直方体
//-----------------------------------------------------------------------------
class rxSolidBox : public rxSolid
{
protected:
	glm::vec3 m_max, m_min;	//!< 最大座標，最小座標(中心からの相対値)

public:
	// コンストラクタ
	rxSolidBox(glm::vec3 minp, glm::vec3 maxp, int sgn)
	{
		glm::vec3 sl  = 0.5f*(maxp-minp);
		glm::vec3 ctr = 0.5f*(maxp+minp);

		m_min = -sl;
		m_max =  sl;
		m_mass_center = ctr;

		m_sgn = sgn;

		m_name = RXS_AABB;
	}

	virtual glm::vec3 GetMin(void)
	{
		glm::vec3 minp, maxp;
		getRotMinMax(minp, maxp);
		return minp;
		//return m_mass_center+m_min; 
	}

	virtual glm::vec3 GetMax(void)
	{
		glm::vec3 minp, maxp;
		getRotMinMax(minp, maxp);
		return maxp;
		//return m_mass_center+m_max; 
	}

	virtual bool GetDistance(const glm::vec3 &pos, rxCollisionInfo &col){ return GetDistanceR(pos, 0.0, col); }
	virtual bool GetDistanceR(const glm::vec3 &pos, const float &r, rxCollisionInfo &col)
	{
		float d;
		glm::vec3 n, spos = CalLocalCoord(pos);
		dist_aabb_point(spos, r, -m_sgn, m_min, m_max, d, n);

		col.Velocity() = GetVelocityAtGrobal(pos);
		col.Penetration() = d; col.Normal() = ApplyRot(n); col.Contact() = pos+n*fabs(d);
		return (col.Penetration() <= 0.0);
	}
	virtual bool GetDistance(const glm::vec3 &pos0, const glm::vec3 &pos1, rxCollisionInfo &col){ return GetDistance(pos1, col); }
	virtual bool GetDistanceR(const glm::vec3 &pos0, const glm::vec3 &pos1, const float &r, rxCollisionInfo &col){ return GetDistanceR(pos1, r, col); }

	virtual bool GetCurvature(const glm::vec3 &pos, float &k){ return calCurvature(pos, k, &rxSolidBox::GetDistance_s, this); }
	virtual void Draw(int drw = 4)
	{
		glPushMatrix();
		SetGLMatrix();
		glm::vec3 len = glm::abs(m_max-m_min);
		glScalef(len[0], len[1], len[2]);
		if(drw & 4){ glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); DrawCubeVBO(); }
		else if(drw & 2){ glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); glLineWidth(3.0); DrawCubeVBO(); }
		glPopMatrix();
	}

	virtual void SetMin(const glm::vec3& p){ m_min = p; }
	virtual void SetMax(const glm::vec3& p){ m_max = p; }


protected:
	//! 回転後の最小最大座標取得
	void getRotMinMax(glm::vec3 &minp, glm::vec3 &maxp)
	{
		// 回転後の8頂点座標
		glm::vec3 cpos[8] = {
			glm::vec3(m_min[0], m_min[1], m_min[2]),
			glm::vec3(m_min[0], m_min[1], m_max[2]),
			glm::vec3(m_min[0], m_max[1], m_min[2]),
			glm::vec3(m_min[0], m_max[1], m_max[2]),
			glm::vec3(m_max[0], m_min[1], m_min[2]),
			glm::vec3(m_max[0], m_min[1], m_max[2]),
			glm::vec3(m_max[0], m_max[1], m_min[2]),
			glm::vec3(m_max[0], m_max[1], m_max[2]) };
		for(int i = 0; i < 8; ++i){
			cpos[i] = CalGlobalCoord(cpos[i]);
		}

		// 回転後の最大最小座標
		minp = maxp = m_mass_center;
		for(int i = 0; i < 8; ++i){
			for(int d = 0; d < 3; ++d){
				if(cpos[i][d] < minp[d]) minp[d] = cpos[i][d];
				if(cpos[i][d] > maxp[d]) maxp[d] = cpos[i][d];
			}
		}

	}

};

//-----------------------------------------------------------------------------
// rxSolidOpenBox : 直方体(開)
//-----------------------------------------------------------------------------
class rxSolidOpenBox : public rxSolid
{
protected:
	glm::vec3 m_slen_inside, m_slen_outside;	// 箱の内部/外部の辺の長さの1/2

public:
	// コンストラクタ
	rxSolidOpenBox(glm::vec3 ctr, glm::vec3 sl_in, glm::vec3 sl_out, int sgn)
	{
		m_slen_inside  = sl_in;
		m_slen_outside = sl_out;
		m_mass_center = ctr;

		//cout << "SLenIn  " << glm::to_string(m_slen_inside) << endl;
		//cout << "SLenOut " << glm::to_string(m_slen_outside) << endl;

		m_sgn = sgn;

		m_name = RXS_OPEN_BOX;
	}

	glm::vec3 GetInMin(void) const { return -m_slen_inside; }
	glm::vec3 GetInMax(void) const { return  m_slen_inside; }
	glm::vec3 GetOutMin(void) const { return -m_slen_outside; }
	glm::vec3 GetOutMax(void) const { return  m_slen_outside; }
	glm::vec3 GetInSideLength(void) const { return m_slen_inside; }
	glm::vec3 GetOutSideLength(void) const { return m_slen_outside; }

	virtual glm::vec3 GetMin(void){ return -m_slen_outside; }
	virtual glm::vec3 GetMax(void){ return  m_slen_outside; }

	virtual bool GetDistance(const glm::vec3 &pos, rxCollisionInfo &col){ return GetDistanceR(pos, 0.0, col); }
	virtual bool GetDistanceR(const glm::vec3 &pos, const float &r, rxCollisionInfo &col)
	{
		glm::vec3 n, spos = CalLocalCoord(pos);
		float d = RX_FEQ_INF;
		dist_openbox_point(spos, r, -m_sgn, m_slen_inside, m_slen_outside, d, n);

		col.Velocity() = GetVelocityAtGrobal(pos);//glm::vec3(0.0);
		col.Penetration() = d; col.Normal() = n; col.Contact() = pos+n*fabs(d);
		return (col.Penetration() <= 0.0);
	}

	virtual bool GetDistance(const glm::vec3 &pos0, const glm::vec3 &pos1, rxCollisionInfo &col){ return GetDistance(pos1, col); }
	virtual bool GetDistanceR(const glm::vec3 &pos0, const glm::vec3 &pos1, const float &r, rxCollisionInfo &col){ return GetDistanceR(pos1, r, col); }

	virtual bool GetCurvature(const glm::vec3 &pos, float &k){ return calCurvature(pos, k, &rxSolidOpenBox::GetDistance_s, this); }
	virtual void Draw(int drw = 4)
	{
		//glPushMatrix();
		//SetGLMatrix();
		//if(drw & 2) DrawWireOpenBox(m_slen_inside, m_slen_outside);
		//else DrawSolidOpenBox(m_slen_inside, m_slen_outside);
		//glPopMatrix();
	}
};



//-----------------------------------------------------------------------------
// rxSolidSphere : 球
//-----------------------------------------------------------------------------
class rxSolidSphere : public rxSolid
{
protected:
	float m_radius;		//!< 半径
	float m_radius_sqr;	//!< 半径の自乗

public:
	// コンストラクタ
	rxSolidSphere(glm::vec3 ctr, float rad, int sgn)
		: m_radius(rad)
	{
		m_sgn = sgn;
		m_mass_center = ctr;
		m_radius_sqr = rad*rad;

		m_name = RXS_SPHERE;
	}

	glm::vec3 GetCenter(void) const { return m_mass_center; }
	float GetRadius(void) const { return m_radius; }

	virtual glm::vec3 GetMin(void){ return m_mass_center-glm::vec3(m_radius); }
	virtual glm::vec3 GetMax(void){ return m_mass_center+glm::vec3(m_radius); }

	virtual bool GetDistance(const glm::vec3 &pos, rxCollisionInfo &col){ return GetDistanceR(pos, 0.0, col); }
	virtual bool GetDistanceR(const glm::vec3 &pos, const float &r, rxCollisionInfo &col)
	{
		glm::vec3 rpos = pos-m_mass_center;
		float d = m_sgn*(glm::length(rpos)-m_radius);
		if(d < r){
			glm::vec3 n = glm::normalize(rpos);

			col.Penetration() = d-r;
			col.Contact() = m_mass_center+n*(m_radius+m_sgn*r);
			col.Normal() = (float)(m_sgn)*n;

			col.Velocity() = GetVelocityAtGrobal(pos);
		} else{
			return false;
		}

		return (col.Penetration() <= 0.0);
	}


	virtual bool GetDistance(const glm::vec3 &pos0, const glm::vec3 &pos1, rxCollisionInfo &col){ return GetDistance(pos1, col); }
	virtual bool GetDistanceR(const glm::vec3 &pos0, const glm::vec3 &pos1, const float &r, rxCollisionInfo &col){ return GetDistanceR(pos1, r, col); }

	virtual bool GetCurvature(const glm::vec3 &pos, float &k){ return false; }
	virtual void Draw(int drw = 4)
	{
		glPushMatrix();
		SetGLMatrix();	// 中心位置へ座標系を移動
		glScalef(2*m_radius, 2*m_radius, 2*m_radius);
		if(drw & 4){ glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); DrawSphereVBO(); }
		else if(drw & 2){ glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); glLineWidth(3.0); DrawSphereVBO(); }
		glPopMatrix();
	}

};






//-----------------------------------------------------------------------------
// rxSolidクラスの実装
//-----------------------------------------------------------------------------
/*!
 * 距離関数から曲率を計算
 * @param[in] pos 計算点
 * @param[out] k 曲率
 * @param[in] dstfunc 距離関数
 */
inline bool rxSolid::calCurvature(const glm::vec3& pos, float& k, bool (dstfunc)(glm::vec3, rxCollisionInfo&, void*), void* fp)
{
	k = 0.0;

	float h = 0.005;
	float x0, y0, z0;
	float p[3][3][3];
	rxCollisionInfo col;

	x0 = pos[0]-0.5*h;
	y0 = pos[1]-0.5*h;
	z0 = pos[2]-0.5*h;

	//	dstfunc(glm::vec3(x0-h, y0-h, z0-h), col, fp); p[0][0][0] = col.Penetration();
	dstfunc(glm::vec3(x0-h, y0-h, z0), col, fp); p[0][0][1] = col.Penetration();
	//	dstfunc(glm::vec3(x0-h, y0-h, z0+h), col, fp); p[0][0][2] = col.Penetration();
	dstfunc(glm::vec3(x0-h, y0, z0-h), col, fp); p[0][1][0] = col.Penetration();
	dstfunc(glm::vec3(x0-h, y0, z0), col, fp); p[0][1][1] = col.Penetration();
	dstfunc(glm::vec3(x0-h, y0, z0+h), col, fp); p[0][1][2] = col.Penetration();
	//	dstfunc(glm::vec3(x0-h, y0+h, z0-h), col, fp); p[0][2][0] = col.Penetration();
	dstfunc(glm::vec3(x0-h, y0+h, z0), col, fp); p[0][2][1] = col.Penetration();
	//	dstfunc(glm::vec3(x0-h, y0+h, z0+h), col, fp); p[0][2][2] = col.Penetration();

	dstfunc(glm::vec3(x0, y0-h, z0-h), col, fp); p[1][0][0] = col.Penetration();
	dstfunc(glm::vec3(x0, y0-h, z0), col, fp); p[1][0][1] = col.Penetration();
	dstfunc(glm::vec3(x0, y0-h, z0+h), col, fp); p[1][0][2] = col.Penetration();
	dstfunc(glm::vec3(x0, y0, z0-h), col, fp); p[1][1][0] = col.Penetration();
	dstfunc(glm::vec3(x0, y0, z0), col, fp); p[1][1][1] = col.Penetration();
	dstfunc(glm::vec3(x0, y0, z0+h), col, fp); p[1][1][2] = col.Penetration();
	dstfunc(glm::vec3(x0, y0+h, z0-h), col, fp); p[1][2][0] = col.Penetration();
	dstfunc(glm::vec3(x0, y0+h, z0), col, fp); p[1][2][1] = col.Penetration();
	dstfunc(glm::vec3(x0, y0+h, z0+h), col, fp); p[1][2][2] = col.Penetration();

	//	dstfunc(glm::vec3(x0+h, y0-h, z0-h), col, fp); p[2][0][0] = col.Penetration();
	dstfunc(glm::vec3(x0+h, y0-h, z0), col, fp); p[2][0][1] = col.Penetration();
	//	dstfunc(glm::vec3(x0+h, y0-h, z0+h), col, fp); p[2][0][2] = col.Penetration();
	dstfunc(glm::vec3(x0+h, y0, z0-h), col, fp); p[2][1][0] = col.Penetration();
	dstfunc(glm::vec3(x0+h, y0, z0), col, fp); p[2][1][1] = col.Penetration();
	dstfunc(glm::vec3(x0+h, y0, z0+h), col, fp); p[2][1][2] = col.Penetration();
	//	dstfunc(glm::vec3(x0+h, y0+h, z0-h), col, fp); p[2][2][0] = col.Penetration();
	dstfunc(glm::vec3(x0+h, y0+h, z0), col, fp); p[2][2][1] = col.Penetration();
	//	dstfunc(glm::vec3(x0+h, y0+h, z0+h), col, fp); p[2][2][2] = col.Penetration();

	float px, py, pz, pxx, pyy, pzz, pxy, pyz, pxz, np;
	px = (p[2][1][1]-p[0][1][1])/(2.0*h);
	py = (p[1][2][1]-p[1][0][1])/(2.0*h);
	pz = (p[1][1][2]-p[1][1][0])/(2.0*h);

	pxx = (p[2][1][1]-2.0*p[1][1][1]+p[0][1][1])/(h*h);
	pyy = (p[1][2][1]-2.0*p[1][1][1]+p[1][0][1])/(h*h);
	pzz = (p[1][1][2]-2.0*p[1][1][1]+p[1][1][0])/(h*h);

	pxy = (p[0][0][1]+p[2][2][1]-p[0][2][1]-p[2][0][1])/(4.0*h*h);
	pxz = (p[0][1][0]+p[2][1][2]-p[0][1][2]-p[2][1][0])/(4.0*h*h);
	pyz = (p[1][0][0]+p[1][2][2]-p[1][0][2]-p[1][2][0])/(4.0*h*h);

	np = px*px+py*py+pz*pz;
	if(np > RX_FEQ_EPS){
		np = sqrt(np);

		// 曲率の計算
		k = (px*px*pyy-2.0*px*py*pxy+py*py*pxx+px*px*pzz-2.0*px*pz*pxz+pz*pz*pxx+py*py*pzz-2.0*py*pz*pyz+pz*pz*pyy)/(np*np*np);
	}

	k = -k;

	return true;
}

	




#endif	// _RX_SPH_SOLID_H_
