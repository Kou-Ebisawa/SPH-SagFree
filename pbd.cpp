/*! 
  @file pbd.cpp
	
  @brief PBDによる弾性体シミュレーション
 
  @author Makoto Fujisawa
  @date 2021
*/

//cloth,volumeに関するConstraintとInitCloth,InitBallそれぞれコメントアウト

//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "pbd.h"
#include "msh.h"
#include "sph.h"

#include <iostream>
#include <memory>
#include <ostream>

//海老沢追加-------------------------------------
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <Eigen/SparseLU>
#include <Eigen/SVD>

Eigen::Vector3d Vec_glmToEigen(const glm::vec3& glmVec) {
	return Eigen::Vector3d(glmVec.x, glmVec.y, glmVec.z);
}

glm::vec3 Vec_eigenToGlm(const Eigen::Vector3d& eigenVec) {
	return glm::vec3(eigenVec.x(), eigenVec.y(), eigenVec.z());
}

Eigen::Quaterniond Quat_glmToEigen(const glm::quat glmQuat) {
	return Eigen::Quaterniond(glmQuat.w, glmQuat.x, glmQuat.y, glmQuat.z);
}

glm::quat Quat_eigenToGlm(const Eigen::Quaterniond eigenQuat) {
	return glm::quat(eigenQuat.w(), eigenQuat.x(), eigenQuat.y(), eigenQuat.z());
}

float g_wq = 1.0e-6;

//-----------------------------------------------------------------------------
// 課題用関数
//-----------------------------------------------------------------------------
/*!
 * PBDのstretching constraintの適用
 *  - エッジ間の長さを保つようにm_vNewPosを更新
 * @param[in] ks 合成パラメータ(スライド中のk,現在の反復回数で変化するので引数として与える)
 */

//CPUでSagFreeの処理をする場合
//行列を用いて個々の毛髪ごとに線形システムを解くこととする
void ElasticPBD::GlobalForceStep(void) {

	using ScalarType = double;
	using IndexType = int64_t;
	using Triplet = Eigen::Triplet<ScalarType, IndexType>;

	using ScalarType = double;  // 以前に指定していたら書く必要はない
	using SparseMatrix = Eigen::SparseMatrix<ScalarType>;

	std::vector<Triplet> triplets;

	const rxEdge& e = m_poly.edges[m_iNumEdge-1];
	int v1 = e.v[0];	// エッジ
	int v2 = e.v[1];
	glm::vec3 p1 = m_poly.vertices[e.v[0]];
	glm::vec3 p2 = m_poly.vertices[e.v[1]];
	float mass = m_vMass[1];//(一律であること前提
	mass = 5.0e-3;

	//固定点を除く頂点で未定乗数を求める
	SparseMatrix A(2 * (m_iNumEdge) * 3, 2 * (m_iNumEdge) * 3);//Num*3でfssの微分とNum*3で未定乗数で微分
	Eigen::VectorXd x(2 * (m_iNumEdge) * 3);
	Eigen::VectorXd b(2 * m_iNumEdge * 3);

	//bの値を初期化
	b = Eigen::VectorXd::Zero(2 * (m_iNumEdge) * 3);

	//行列の値を設定
	for (int i = 0; i < 2*(m_iNumEdge) * 3; i++) {
		//fssへの要素についての偏微分
		if (i < (m_iNumEdge) * 3) {//fssの各要素の値に2を代入
			triplets.emplace_back(i, i, 2);
			triplets.emplace_back(i, (m_iNumEdge) * 3 + i, -1);
			if (i > 2)triplets.emplace_back(i, (m_iNumEdge) * 3 + i - 3, 1);
		}
		else{//未定乗数での微分に1を代入
			if (i < (2 * m_iNumEdge - 1) * 3) triplets.emplace_back(i, i + 3 - 3 * (m_iNumEdge), 1);
			triplets.emplace_back(i, i - 3 * (m_iNumEdge), -1);
		}
	}

	//Aに値を設定
	A.setFromTriplets(triplets.begin(), triplets.end());

	//右辺値を設定
	for (int i = 0; i < (m_iNumEdge); i++) {
		b[3 * (m_iNumEdge) + 3 * i] = mass*m_v3Gravity.x;
		b[3 * (m_iNumEdge) + 3 * i + 1] = mass * m_v3Gravity.y;
		b[3 * (m_iNumEdge) + 3 * i + 2] = mass * m_v3Gravity.z;
		if (i == m_iNumEdge - 1) {
			b[3 * (m_iNumEdge) + 3 * i] = mass * m_v3Gravity.x;
			b[3 * (m_iNumEdge) + 3 * i + 1] = mass * m_v3Gravity.y;
			b[3 * (m_iNumEdge) + 3 * i + 2] = mass * m_v3Gravity.z;
		}
	}

	Eigen::SparseLU<SparseMatrix>solver;
	solver.compute(A);
	x = solver.solve(b);

	for (int i = 0; i < (m_iNumEdge) * 3; i++) {
		m_Fss[i / 3][i % 3] = x[i];
	}
}

void ElasticPBD::LocalForceStep(float ks) {
	for (int i = 0; i < m_iNumEdge; i++) {
		if (i == 0) continue;//最初の四元数を固定した場合
		float l1, l2;
		const rxEdge& e = m_poly.edges[i];
		glm::vec3 p1 = m_poly.vertices[e.v[0]];//xp
		glm::vec3 p2 = m_poly.vertices[e.v[1]];//xp+
		float l0 = m_vLengths[i];

		glm::vec3 ds = glm::normalize(m_vQuat[i] * glm::vec3(0, 0, 1));
		ds = glm::normalize(ds);

		glm::vec3 Fss_i = m_Fss[i];//m_Fss[i]を置き換え
		float fs_Len = glm::length(Fss_i);

		//ksが判別式を満たすかどうかを判断
		float tmp_ks = ks;
		//判別式を満たさない場合にksに追加する量
		float delta = 10.f;
		
		//判別式を満たすまでのループ
		while (1) {
			//float B = glm::dot((p1 - p2), m_Fss[i]) / (tmp_ks * tmp_ks) + 1;//論文にある通り
			float B = glm::dot((p1 - p2), Fss_i) / (tmp_ks)+1;//自分の感じた方
			float AC = glm::length2(p1 - p2) * glm::length2(Fss_i) / (tmp_ks * tmp_ks);
			//判別式
			float discrim = B * B - 4 * AC;
			if (discrim >= 0) break;
			tmp_ks += delta;
		}

		//更新後のkssを配列に保存
		m_NewKss[i] = tmp_ks;

		//||ds,3||=1(Eq.14)から求める------------------------------------
		float a = fs_Len * fs_Len;
		float b = -tmp_ks * (2.f * glm::dot(Fss_i, p1 - p2) + tmp_ks);
		float c = tmp_ks * tmp_ks * glm::length2(p1 - p2);

		l1 = sqrt((-b + sqrt(abs(pow(b, 2.0) - 4.f * a * c))) / (2.f * a));
		l2 = sqrt((-b - sqrt(abs(pow(b, 2.0) - 4.f * a * c))) / (2.f * a));
		//---------------------------------------------------------------

		//2次解のうち、適切な長さを選択する------------------------------------------------------------
		float Cur_l = m_vLengths[i];//現在の長さ
		
		//元の長さの出力と更新後の長さの候補
		//if (i == 0) cout << "rest_length " << Cur_l << "-----------------------------------------------------------------------------------" << endl;
		//最も長さが近いものを選んで基準長の更新
		if (abs(l1 - Cur_l) > abs(l2 - Cur_l) && abs(l2) > 1.0e-10) m_vLengths[i] = l2;
		else if (abs(l1) > 1.0e-10)m_vLengths[i] = l1;
		else {
			cout << "Not Updated In LocalForceStep!!" << endl;
		}
		//----------------------------------------------------------------------------------------------

		//姿勢の更新--------------------------------------------------------------------------------------------------------
		glm::vec3 new_ds = (m_vLengths[i] * Fss_i) / tmp_ks - (p1 - p2) / m_vLengths[i];//直感と+-逆だが、これで直線はよし(Eq.14上の式のds3の符号はおそらく+)
		Eigen::Quaterniond new_qs= Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(), Vec_glmToEigen(new_ds));
		//------------------------------------------------------------------------------------------------------------------

		//結果の確認-----------------------------------------------------------
		//glm::vec3 new_Fss = (tmp_ks / m_vLengths[i]) * ((p1 - p2) / m_vLengths[i] + new_ds);//最後の姿勢を+に変更
		//cout << "NewFss" << i << ": " << glm::to_string(new_Fss) << endl;
		//---------------------------------------------------------------------

		m_vQuat[i] = Quat_eigenToGlm(new_qs);
	}
}

//GlobalTorqueStepを線形システムで解く
void ElasticPBD::VideoGlobalTorqueStep(float kbt) {
	using ScalarType = double;
	using IndexType = int64_t;
	using Triplet = Eigen::Triplet<ScalarType, IndexType>;

	using ScalarType = double;  // 以前に指定していたら書く必要はない
	using SparseMatrix = Eigen::SparseMatrix<ScalarType>;

	std::vector<Triplet> triplets;

	//全てのエッジで四元数を解く
	SparseMatrix A(m_iNumEdge * 4, 4*m_iNumEdge);//Num*3でfssの微分とNum*3で未定乗数で微分
	Eigen::VectorXd x(m_iNumEdge * 4);
	Eigen::VectorXd b(m_iNumEdge * 4);

	//bの値を初期化
	b = Eigen::VectorXd::Zero(m_iNumEdge * 4);

	for (int i = 0; i < m_iNumEdge; i++) {
		const rxEdge& e = m_poly.edges[i];
		glm::vec3 p1 = m_poly.vertices[e.v[0]];//xp
		glm::vec3 p2 = m_poly.vertices[e.v[1]];//xp+
		float l0 = m_vLengths[i];//フォースステップでの更新後の長さを利用
		glm::quat qs = m_vQuat[i];
		float ks = m_NewKss[i];

		glm::vec4 torqueSS;
		glm::vec3 V = (p1 - p2) / l0;
		torqueSS[0] = 4 * (qs.x * qs.x * qs.x + qs.x * qs.y * qs.y + qs.x * qs.z * qs.z + qs.x * qs.w * qs.w + V.z * qs.x - V.x * qs.z + V.y * qs.w);//x
		torqueSS[1] = 4 * (qs.y * qs.y * qs.y + qs.y * qs.x * qs.x + qs.y * qs.z * qs.z + qs.y * qs.w * qs.w + V.z * qs.y - V.x * qs.w - V.y * qs.z);//y
		torqueSS[2] = 4 * (qs.z * qs.z * qs.z + qs.z * qs.x * qs.x + qs.z * qs.y * qs.y + qs.z * qs.w * qs.w - V.z * qs.z - V.x * qs.x - V.y * qs.y);//z
		torqueSS[3] = 4 * (qs.w * qs.w * qs.w + qs.w * qs.x * qs.x + qs.w * qs.y * qs.y + qs.w * qs.z * qs.z - V.z * qs.w - V.x * qs.y + V.y * qs.x);//w
		torqueSS = 1.f / 2.f * ks * torqueSS;

		float K = 2 * kbt / (l0);//kに当たる部分
		float s0, s1;//曲げ剛性の方向についてとりあえずpositiveだけを考える(本来はエッジごとに指定されるため、Omega_{i,0}に掛けられるsは常に同じ
		s0 = s1 = 1;

		//Aの設定だが、トルクに逆元をかける場合、ここは方向の設定のみ(方向はpositiveに限定)
		if (i == 0) {
			for (int j = 0; j < 4; j++) {
				triplets.emplace_back(4 * i + j, 4 * i + j, s0);//四元数の要素の数だけs0を追加
			}
		}
		else {
			for (int j = 0; j < 4; j++) {
				triplets.emplace_back(4 * i + j, 4 * (i - 1) + j, s0);//四元数の要素の数だけs0を追加
				triplets.emplace_back(4 * i + j, 4 * i + j, s1);
			}
		}

		//右辺値を設定
		if (i == 0) {
			glm::quat kq = K * m_vQuat[i];//定数と四元数の積
			glm::quat kq_inv = glm::inverse(kq);//\tauに欠ける逆元
			glm::quat quat_torqueSS(torqueSS[3], torqueSS[0], torqueSS[1], torqueSS[2]);

			//直接利用する変数
			glm::quat Omega_Next = glm::conjugate(glm::conjugate(m_vQuat[i]) * m_vQuat[i + 1]);//\bar{\bar{q_i}q_{i+1}}に対応
			glm::quat New_tau = kq_inv * quat_torqueSS;//\tau_iをkqの逆数で割ったもの

			b[4 * i] = Omega_Next.x - New_tau.x;
			b[4 * i + 1] = Omega_Next.y - New_tau.y;
			b[4 * i + 2] = Omega_Next.z - New_tau.z;
			b[4 * i + 3] = Omega_Next.w - New_tau.w;
		}
		else {
			glm::quat kq = K * m_vQuat[i];//定数と四元数の積
			glm::quat kq_inv = glm::inverse(kq);//\tauに欠ける逆元
			glm::quat quat_torqueSS(torqueSS[3], torqueSS[0], torqueSS[1], torqueSS[2]);

			//直接利用する変数
			glm::quat Omega_Prev = glm::conjugate(m_vQuat[i - 1]) * m_vQuat[i];//\bar{q_{i-1}}q_iに対応
			glm::quat Omega_Next = glm::conjugate(glm::conjugate(m_vQuat[i]) * m_vQuat[i + 1]);//\bar{\bar{q_i}q_{i+1}}に対応
			glm::quat New_tau = kq_inv * quat_torqueSS;//\tau_iをkqの逆数で割ったもの

			b[4 * i] = Omega_Prev.x + Omega_Next.x - New_tau.x;
			b[4 * i + 1] = Omega_Prev.y + Omega_Next.y - New_tau.y;
			b[4 * i + 2] = Omega_Prev.z + Omega_Next.z - New_tau.z;
			b[4 * i + 3] = Omega_Prev.w + Omega_Next.w - New_tau.w;
		}
	}

	//Aに値を設定
	A.setFromTriplets(triplets.begin(), triplets.end());
	Eigen::MatrixXd A_dense = Eigen::MatrixXd(A);
	Eigen::MatrixXd A_inv = A_dense.inverse();

	Eigen::SparseLU<SparseMatrix>solver;
	solver.compute(A);
	x = solver.solve(b);

	glm::quat New_Quat;
	for (int i = 0; i < m_iNumEdge; i++) {
		New_Quat.x = x[i * 4];
		New_Quat.y = x[i * 4 + 1];
		New_Quat.z = x[i * 4 + 2];
		New_Quat.w = x[i * 4 + 3];

		m_eOrgOmega[i] = New_Quat;
	}
}

glm::vec4 FourDtorqueSolver(glm::quat q1, glm::quat q2, glm::quat Darboux) {//1セグメント前の曲げ・ねじれのトルクを計算する(2つ目以降のτの右辺値)
	glm::vec4 torque;
	float sum = q2.x * q2.x + q2.y * q2.y + q2.z * q2.z + q2.w * q2.w;
	torque[0] = sum * q1.x - q2.x * Darboux.w + q2.y * Darboux.z - q2.z * Darboux.y + q2.w * Darboux.x;
	torque[1] = sum * q1.y - q2.x * Darboux.z - q2.y * Darboux.w + q2.z * Darboux.x + q2.w * Darboux.y;
	torque[2] = sum * q1.z + q2.x * Darboux.y - q2.y * Darboux.x - q2.z * Darboux.w + q2.w * Darboux.z;
	torque[3] = sum * q1.w - q2.x * Darboux.x - q2.y * Darboux.y - q2.z * Darboux.z - q2.w * Darboux.w;

	return torque;
}

glm::quat SecondFourDtorqueSolver(glm::quat q1, glm::quat q2, glm::vec4 torqueSS) {
	//手計算から求める
	//bに右辺値を設定(torqueSSに他の右辺値も入れておく)
	Eigen::MatrixXd A(4, 4);//(4,3)
	Eigen::VectorXd x(4);
	Eigen::VectorXd b(4);

	float sum = q2.x * q2.x + q2.y * q2.y + q2.z * q2.z + q2.w * q2.w;
	//新たな基準ダルボーベクトルにかかわる係数を行列とする------
	A << q2.w, -q2.z, q2.y, -q2.x,
		q2.z, q2.w, -q2.x, -q2.y,
		-q2.y, q2.x, q2.w, -q2.z,
		-q2.x, -q2.y, -q2.z, -q2.w;
	//-----------------------------------------------------------

	//右辺値を設定-----------------------------------------------
	b = Eigen::VectorXd::Zero(4);
	for (int i = 0; i < 4; i++) {
		b[i] = torqueSS[i];
	}

	b[0] -= sum * q1.x;
	b[1] -= sum * q1.y;
	b[2] -= sum * q1.z;
	b[3] -= sum * q1.w;
	//------------------------------------------------------------

	x = A.partialPivLu().solve(b);//Eigenソルバにより、行列を解く


	glm::quat quat;
	quat = glm::quat(x[3], x[0], x[1], x[2]);

	return quat;
}

void ElasticPBD::GlobalTorqueStep(float kbt) {
	for (int i = 0; i < m_iNumEdge-1; i++) {
		const rxEdge& e = m_poly.edges[i];
		glm::vec3 p1 = m_poly.vertices[e.v[0]];//xp
		glm::vec3 p2 = m_poly.vertices[e.v[1]];//xp+
		float l0 = m_vLengths[i];//フォースステップでの更新後の長さを利用
		glm::quat qs = m_vQuat[i];//更新後の四元数
		float ks = m_NewKss[i];//更新後の剛性

		glm::vec4 torqueSS;
		glm::vec3 V = (p1 - p2)/l0;
		torqueSS[0] = 4 * (qs.x * qs.x * qs.x + qs.x * qs.y * qs.y + qs.x * qs.z * qs.z + qs.x * qs.w * qs.w + V.z * qs.x - V.x * qs.z + V.y * qs.w);//x
		torqueSS[1] = 4 * (qs.y * qs.y * qs.y + qs.y * qs.x * qs.x + qs.y * qs.z * qs.z + qs.y * qs.w * qs.w + V.z * qs.y - V.x * qs.w - V.y * qs.z);//y
		torqueSS[2] = 4 * (qs.z * qs.z * qs.z + qs.z * qs.x * qs.x + qs.z * qs.y * qs.y + qs.z * qs.w * qs.w - V.z * qs.z - V.x * qs.x - V.y * qs.y);//z
		torqueSS[3] = 4 * (qs.w * qs.w * qs.w + qs.w * qs.x * qs.x + qs.w * qs.y * qs.y + qs.w * qs.z * qs.z - V.z * qs.w - V.x * qs.y + V.y * qs.x);//w
		torqueSS = 1.f / 2.f * ks * torqueSS;

		//BTのエネルギー,トルクを求める
		if (i == 0) {
			glm::quat c_Omega = glm::conjugate(m_vQuat[i]) * m_vQuat[i + 1];//現在のダルボーベクトル
			glm::quat c_Darboux = m_eOrgOmega[i];//基準ダルボーベクトル、φを求めないと制約が出せないので、ひとまずこれを利用
			//Eigenメソッドを使えるように
			Eigen::Quaterniond Current_eigen = Quat_glmToEigen(c_Omega);
			Eigen::Quaterniond Omega_eigen = Quat_glmToEigen(c_Darboux);
			//制約に使うダルボーベクトルの方向を求める---------------------
			Eigen::Quaterniond W1, W2;
			W1.coeffs() = Current_eigen.coeffs() - Omega_eigen.coeffs();
			W2.coeffs() = Current_eigen.coeffs() + Omega_eigen.coeffs();

			float s;
			if (W1.squaredNorm() < W2.squaredNorm()) {
				s = 1;
			}
			else {
				s = -1;
			}
			//s = 1;//sを除く
			//------------------------------------------------------------
			glm::quat New_Darboux=SecondFourDtorqueSolver(m_vQuat[i], m_vQuat[i+1], -torqueSS / kbt);

			//更新した値を計算に用いる場合
			m_eOrgOmega[i] = (s * New_Darboux);
		}
		else {
			//更新前のダルボーベクトルも含めてデータを持ってくる---------------------------------------
			glm::quat p_Omega = glm::conjugate(m_vQuat[i - 1]) * m_vQuat[i];//既知
			glm::quat p_Darboux = m_eOrgOmega[i - 1];//既知

			glm::quat c_Omega = glm::conjugate(m_vQuat[i]) * m_vQuat[i + 1];//既知
			glm::quat c_Darboux = m_eOrgOmega[i];//φを求めないと制約が出せないので、更新前の値を利用
			//-----------------------------------------------------------------------------------------

			//ひとつ前のダルボーベクトルを求め、姿勢の方向を導出---------------------------------------
			Eigen::Quaterniond p_Omega_eigen = Quat_glmToEigen(p_Omega);
			Eigen::Quaterniond p_Darboux_eigen = Quat_glmToEigen(p_Darboux);

			Eigen::Quaterniond W1, W2;
			W1.coeffs() = p_Omega_eigen.coeffs() - p_Darboux_eigen.coeffs();
			W2.coeffs() = p_Omega_eigen.coeffs() + p_Darboux_eigen.coeffs();

			float s1;
			if (W1.squaredNorm() < W2.squaredNorm()) {
				s1 = 1;
			}
			else {
				s1 = -1;
			}
			//s1 = 1;
			//-----------------------------------------------------------------------------------------

			//ひとつ前のダルボーベクトルに関するEbtを求める--------------------------------------------
			glm::vec4 torqueBT_prev;//1/2は掛けない
			glm::vec4 tmp;//Ωc-φΩ0を代入
			glm::quat qs = m_vQuat[i-1];

			glm::vec4 prev_torqueBT = FourDtorqueSolver(m_vQuat[i - 1], m_vQuat[i], s1*p_Darboux);
			prev_torqueBT *= kbt;
			//--------------------------------------------------------------------------------------------
			
			//現在のダルボーベクトルを求める------------------------------------
			Eigen::Quaterniond c_Omega_eigen = Quat_glmToEigen(c_Omega);
			Eigen::Quaterniond c_Darboux_eigen = Quat_glmToEigen(c_Darboux);

			Eigen::Quaterniond W3, W4;
			W3.coeffs() = c_Omega_eigen.coeffs() - c_Darboux_eigen.coeffs();
			W4.coeffs() = c_Omega_eigen.coeffs() + c_Darboux_eigen.coeffs();

			float s2;
			if (W3.squaredNorm() < W4.squaredNorm()) {
				s2 = 1;
			}
			else {
				s2 = -1;
			}
			//s2 = 1;
			//-------------------------------------------------------------------

			//右辺値を設定-------------------------------------------------------
			glm::vec4 rightForm;
			rightForm[0] = -torqueSS[0] - prev_torqueBT[0];
			rightForm[1] = -torqueSS[1] - prev_torqueBT[1];
			rightForm[2] = -torqueSS[2] - prev_torqueBT[2];
			rightForm[3] = -torqueSS[3] - prev_torqueBT[3];
			//-------------------------------------------------------------------
			glm::quat New_Darboux = SecondFourDtorqueSolver(m_vQuat[i], m_vQuat[i + 1], rightForm/kbt);

			//更新後の値を計算に利用
			m_eOrgOmega[i] = s2*New_Darboux;
		}
	}
}

//LocalTorqueStepにおける基準ダルボーベクトルを求め、曲げ剛性を上げる処理
glm::quat solveInverseRot(glm::quat quat, glm::quat tar, float& bendK) {
	const float SAFETY_FACTOR = glm::min(quat.w, 0.002f);//0.00002f,0.002f
	glm::vec4 quat_vec = glm::vec4(quat.x, quat.y, quat.z, quat.w);
	glm::vec4 tar_vec = glm::vec4(tar.x, tar.y, tar.z, tar.w);

	//cout << "reminder Length" << glm::length(glm::vec3(quat.x, quat.y, quat.z)) << endl;

	tar_vec -= quat_vec * glm::dot(tar_vec, quat_vec);
	glm::vec4 Omega_vec = -tar_vec / bendK;

	if (glm::length2(Omega_vec) > SAFETY_FACTOR * SAFETY_FACTOR) {
		//新しい曲げ剛性を求める
		bendK = glm::length(tar_vec) / SAFETY_FACTOR;

		//Ωの再計算
		Omega_vec = -tar_vec / bendK;
		//cout << "Omega_updated!" << endl;
	}

	glm::vec4 d = sqrt(1 - glm::length2(Omega_vec)) * quat_vec;
	Omega_vec += d;
	return glm::quat(Omega_vec[3], Omega_vec[0], Omega_vec[1], Omega_vec[2]);
}

void ElasticPBD::LocalTorqueStep(float kbt) {
	//cout << "LastDarboux" << endl;
	for (int i = 0; i < m_iNumEdge-1; i++) {
		glm::quat q1 = m_vQuat[i];
		glm::quat q2 = m_vQuat[i + 1];
		glm::quat qs = glm::conjugate(q1) * q2;//現在のダルボーベクトル(Omegaに当たる)
		//cout << "qs w: " << qs.w << " x:" << qs.x << " y:" << qs.y << " z:" << qs.z << endl;
		glm::quat Omega = m_eOrgOmega[i];//グローバルステップ更新後の基準ダルボーベクトル
	
		float tmp_kbt = kbt;
		glm::quat Last_Omega = solveInverseRot(qs, Omega, tmp_kbt);
		m_eOrgOmega[i] = Last_Omega;
		m_NewKbt[i] = tmp_kbt;
	}
}

void ElasticPBD::SetCoefficient(float ks, float kbt) {
	PBD_ks = ks;
	PBD_kbt = kbt;

	return;
}

void ElasticPBD::projectStretchingConstraint(float dt)
{
	if(m_iNumEdge <= 1) return;

	for (int i = 0; i < m_iNumEdge; ++i) {
		// 四面体を使うときの内部エッジかどうかの判定＆内部エッジを使うかどうかのフラグチェック
		if (m_vInEdge[i] && !m_bUseInEdge) continue;

		// エッジ情報の取得とエッジ両端の頂点番号および質量の取得(固定点の質量は大きくする)
		const rxEdge& e = m_poly.edges[i];
		int v1 = e.v[0];	// エッジ
		int v2 = e.v[1];

		float m1 = m_vMass[v1];
		float m2 = m_vMass[v2];

		//-------------------------------------
		if (m1 < glm::epsilon<float>() || m2 < glm::epsilon<float>()) continue;

		// 2頂点の位置ベクトル
		glm::vec3 p1 = m_vNewPos[v1];
		glm::vec3 p2 = m_vNewPos[v2];

		// 計算点間の元の長さ(制約条件)
		float d = m_vLengths[i];

		glm::vec3 dp1, dp2;
		glm::vec3 Xdp1, Xdp2;
		dp1, dp2 = glm::vec3(0);

		float w1 = 1 / m1;
		float w2 = 1 / m2;
		float wq = 1.0e-6;//q(姿勢)における仮定の重さ

		//(strech&shear_Constraint)---------------------------------------------------------
		//XPBD---------------------------------------------------------------------------------------
		float DT = dt;//時間ステップ
		float A;//コンプライアンス(大きくなるほど柔らかくなる)
		float tmp_ks = m_NewKss[i];//実際に用いる値
		A = 1 / (tmp_ks * DT * DT);//α(チルダ)ss

		//e3の回転
		glm::vec3 ds3 = m_vQuat[i] * glm::vec3(0, 0, 1);

		float weight = w1 + w2 + d * d * (4 * wq + A);//分母
		glm::vec3 mole = d * (p1 - p2 + d * ds3 - d * A * m_LamdaSS[i]);//lamda(分子)
		glm::vec3 lamda = mole / weight;//Δλ

		m_LamdaSS[i] += lamda;//ラグランジュ乗数の更新

		//変形の方向(論文中(Eq.3の少し下)の方法とは符号を反転させている。
		Xdp1 = -lamda / d;
		Xdp2 = lamda / d;

		glm::quat quat2 = -2.f * glm::quat(0.f, lamda) * m_vQuat[i] * glm::conjugate(glm::quat(0, 0, 0, 1));
		//XPBDによる四元数の更新
		m_vQuat[i] = glm::normalize(m_vQuat[i] + quat2);

		//XPBDの位置の更新
		if (!m_vFix[v1]) m_vNewPos[v1] += Xdp1;
		if (!m_vFix[v2]) m_vNewPos[v2] += Xdp2;
		//---------------------------------------------------------------------------------------------
	}
}

void ElasticPBD::projectBendingConstraint(float dt) {
	if (m_iNumEdge <= 1) return;

	//bend&twist_Constraint-------------------------------------------------------------------------
	for (int i = 0; i < m_iNumEdge - 1; ++i) {

		//XPBDの曲げ・ねじれ制約
		if (m_eOrgOmega.size() > 0) {
			// エッジ情報の取得とエッジ両端の頂点番号および質量の取得(固定点の質量は大きくする)
			const rxEdge& e = m_poly.edges[i];
			int v1 = e.v[0];	// エッジ
			int v2 = e.v[1];

			float m1 = m_vMass[v1];
			float m2 = m_vMass[v2];

			float wq = 2 / (m1 + m2);//q(姿勢)における仮定の重さ
			wq = 1.0e-6;

			glm::quat q1 = m_vQuat[i];
			glm::quat q2 = m_vQuat[i + 1];

			//-------------------------------------
			if (m1 < glm::epsilon<float>() || m2 < glm::epsilon<float>()) continue;

			glm::quat Current_glm = glm::conjugate(q1) * q2;
			glm::quat Omega_glm = m_eOrgOmega[i];

			Eigen::Quaterniond Current_eigen = Quat_glmToEigen(Current_glm);
			Eigen::Quaterniond Omega_eigen = Quat_glmToEigen(Omega_glm);

			Eigen::Quaterniond W1, W2;
			W1.coeffs() = Current_eigen.coeffs() - Omega_eigen.coeffs();
			W2.coeffs() = Current_eigen.coeffs() + Omega_eigen.coeffs();

			glm::quat s;
			if (W1.squaredNorm() < W2.squaredNorm()) {
				s = Quat_eigenToGlm(W1);
			}
			else {
				s = Quat_eigenToGlm(W2);
			}
			

			float DT = dt;//時間ステップ
			float A;//コンプライアンス(大きくなるほど柔らかくなる)
			float tmp_kbt = m_NewKbt[i];//更新された曲げ剛性を利用
			
			A = 1 / (tmp_kbt * DT * DT);//α(チルダ)bt
			float weight = 2 * wq + A;//エッジの仮定の重みが一律であると仮定

			//lamdaをquatで指定
			glm::quat lamda;
			lamda = (-s + (-A * m_LamdaBT[i])) / weight;

			m_LamdaBT[i] = m_LamdaBT[i] + lamda;

			//C^btによる四元数の更新-----------------------------------------------------------
			m_vQuat[i] = glm::normalize(m_vQuat[i] + q2 * glm::conjugate(lamda));
			m_vQuat[i + 1] = glm::normalize(m_vQuat[i + 1] + q1 * lamda);

			//実部を0にする場合
			glm::quat delta_q1 = q2 * glm::conjugate(lamda);
			glm::quat delta_q2 = q1 * lamda;

			delta_q1.w = delta_q2.w = 0;

			m_vQuat[i] = glm::normalize(m_vQuat[i] + delta_q1);
			m_vQuat[i + 1] = glm::normalize(m_vQuat[i + 1] + delta_q2);
		}
	}
}

//海老沢追加
//球との衝突を扱う
//center:球の中心
//rad:球の半径
void ElasticPBD::projectCollisionConstraint(float dt) {
	if (m_iNumVertices < 1) return;

	for (int i = 0; i < m_iNumVertices; i++) {
		if (m_vFix[i]) continue;

		glm::vec3 pos = m_vNewPos[i];
		glm::vec3 d = pos - m_center;

		float l = m_rad - glm::length(d);
		if (l <= 0) continue;

		glm::vec3 n = glm::normalize(d);
		m_vVel[i] -= (glm::dot(n, m_vVel[i])) * n;
		m_vVel[i].y = 0;//球との衝突によるy座標方向への変形は抑えてしまう
		m_vNewPos[i] += l * n;
	}
}

//海老沢追加
//衝突する球の設定
//center:球の中心
//rad:球の半径
void ElasticPBD::SetCollisionSphere(glm::vec3 c, float r) {
	m_center = c;
	m_rad = r;
}

void ElasticPBD::SagFree(float ks, float kbt) {
	GlobalForceStep();
	LocalForceStep(ks);
	//GlobalTorqueStep(kbt);
	VideoGlobalTorqueStep(kbt);
	LocalTorqueStep(kbt);
}

//普通の位置ベース法
void ElasticPBD::PBDConstraint(float ks) {
	for (int i = 0; i < m_iNumEdge; ++i) {
		// 四面体を使うときの内部エッジかどうかの判定＆内部エッジを使うかどうかのフラグチェック
		if (m_vInEdge[i] && !m_bUseInEdge) continue;

		// エッジ情報の取得とエッジ両端の頂点番号および質量の取得(固定点の質量は大きくする)
		const rxEdge& e = m_poly.edges[i];
		int v1 = e.v[0];	// エッジ
		int v2 = e.v[1];

		float m1 = m_vMass[v1];
		float m2 = m_vMass[v2];

		//-------------------------------------
		if (m1 < glm::epsilon<float>() || m2 < glm::epsilon<float>()) continue;

		// 2頂点の位置ベクトル
		glm::vec3 p1 = m_vNewPos[v1];
		glm::vec3 p2 = m_vNewPos[v2];

		// 計算点間の元の長さ(制約条件)
		float d = m_vLengths[i];
		float w1 = 1 / m1;
		float w2 = 1 / m2;
		float wq = 1.0e-6;
		
		glm::vec3 e3(0, 0, 1);
		glm::vec3 d3 = e3 * m_vQuat[i];

		//calc delta pos
		glm::vec3 gamma = (p2 - p1) / d - d3;
		gamma /= (w2 + w1) / d + 4.0f * wq * d;
		glm::vec3 delta_p1 = w1 * gamma;
		glm::vec3 delta_p2 = -w2 * gamma;

		//calc delta q
		glm::quat e3q = glm::quat(0.0, e3);
		glm::quat q_e3_bar = m_vQuat[i] * glm::conjugate(e3q);
		glm::quat gammaq = glm::quat(0.0, gamma);
		glm::quat delta_q = wq * d * gammaq * q_e3_bar;

		if (!m_vFix[v1]) m_vNewPos[v1] += ks*delta_p1;
		if (!m_vFix[v2]) m_vNewPos[v2] += ks*delta_p2;

		m_vQuat[i] = glm::normalize(m_vQuat[i] + delta_q);
	}
}

void ElasticPBD::onlyStretch(float ks) {
	for (int i = 0; i < m_iNumEdge; ++i) {
		// 四面体を使うときの内部エッジかどうかの判定＆内部エッジを使うかどうかのフラグチェック
		if (m_vInEdge[i] && !m_bUseInEdge) continue;

		// エッジ情報の取得とエッジ両端の頂点番号および質量の取得(固定点の質量は大きくする)
		const rxEdge& e = m_poly.edges[i];
		int v1 = e.v[0];	// エッジ
		int v2 = e.v[1];

		float m1 = m_vMass[v1];
		float m2 = m_vMass[v2];

		//-------------------------------------
		if (m1 < glm::epsilon<float>() || m2 < glm::epsilon<float>()) continue;

		// 2頂点の位置ベクトル
		glm::vec3 p1 = m_vNewPos[v1];
		glm::vec3 p2 = m_vNewPos[v2];

		// 計算点間の元の長さ(制約条件)
		float d = m_vLengths[i];
		float w1 = 1 / m1;
		float w2 = 1 / m2;
		float wq = 1.0e-6;

		glm::vec3 n = p1 - p2;
		float l = glm::length(n);

		if (l < 1.0e-6) continue;
		n /= l;

		// 頂点の移動量
		glm::vec3 dp1 = -w1 / (w1 + w2) * (l - d) * n;
		glm::vec3 dp2 = w2 / (w1 + w2) * (l - d) * n;

		// ----課題ここまで----

		// 頂点位置を修正
		if (!m_vFix[v1]) m_vNewPos[v1] += dp1;
		if (!m_vFix[v2]) m_vNewPos[v2] += dp2;
	}
}


/*!
 * シミュレーションステップを進める
 */
void ElasticPBD::Update(float dt)
{
	// 外力で速度を更新＆予測位置p'の計算
	calExternalForces(dt);
	// 衝突処理(p'を更新)
	//genCollConstraints(dt);

	if (PBDstep == 0) {
		//SagFree(PBD_ks, PBD_kbt);
	}

	for (int j = 0; j < m_iNumEdge; ++j) {//反復前にラグランジュ乗数を初期化
		m_LamdaSS[j] = glm::vec3(0.0);
		m_LamdaBT[j] = glm::quat(0.f,glm::vec3(0.f));
	}

	// 制約条件を満たすための反復処理
	for(int i = 0; i < m_iNmax; ++i){//m_iNmax
		//ともに配列に入れられたグローバル変数を利用
		float ks = 1.0f - pow(1.0f - m_fK, 1.0f / static_cast<float>((i + 1)));
		projectStretchingConstraint(dt);
		projectBendingConstraint(dt);
		//projectCollisionConstraint(dt);

	}

	// 速度と位置の更新
	integrate(dt);
	PBDstep++;
	//ステップの表示
	//cout << PBDstep << "step ---------------------------------------------------------" << endl;
}


//-----------------------------------------------------------------------------
// PBDクラスの実装
//-----------------------------------------------------------------------------
/*!
 * デフォルトコンストラクタ
 */
ElasticPBD::ElasticPBD(int obj)
{
	m_iObjectNo = obj;
	m_bUseInEdge = true;

	m_v3Gravity = glm::vec3(0.0, -9.81, 0.0);
	//m_v3Gravity = glm::vec3(0.0);//重力なし
	//m_v3Gravity = 20.f*glm::vec3(0.0, -0.5, 0.0);//元の重力設定
	m_fWind = 0.0f;

	//m_v3Min = glm::vec3(-1.0,-10.0,-10.0);//海老沢変更 当たり判定を下に広げる
	m_v3Min = glm::vec3(-100.0);//衝突判定を無視(大きいサイズ感であると衝突の可能性あり)
	m_v3Max = glm::vec3(100.0);//

	m_fK = 0.9;
	m_iNmax = 64;//最大反復回数(元は3)
	m_fpCollision = 0;

	PBDstep = 0;//海老沢追加、ステップ数0に設定

	Clear();
}

/*!
 * デストラクタ
 */
ElasticPBD::~ElasticPBD()
{
}

/*!
 * 全頂点を消去
 */
void ElasticPBD::Clear()
{
	m_iNumVertices = 0;
	m_iNumEdge = 0;
	m_iNumTris = 0;
	m_vCurPos.clear();
	m_vNewPos.clear();

	m_vQuat.clear(); //海老沢追加
	m_eOrgOmega.clear();//海老沢追加
	m_vTang.clear();//海老沢追加

	m_vMass.clear();
	m_vVel.clear();
	m_vFix.clear();
	m_poly.Clear();
	m_vLengths.clear();
	m_vInEdge.clear();

	//剛性の配列の削除
	m_NewKss.clear();
	m_NewKbt.clear();
}

/*!
 * 描画関数
 * @param[in] drw 描画フラグ
 */
void ElasticPBD::Draw(int drw)
{
	if(drw & 0x02){
		// エッジ描画における"stitching"をなくすためのオフセットの設定
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
	}

	// 頂点描画
	if(drw & 0x0001){
		glDisable(GL_LIGHTING);
		glPointSize(5.0);
		glColor3d(1.0, 0.3, 0.3);
		for(int i = 0; i < m_iNumVertices; ++i){
			glm::vec3 v = m_vCurPos[i];
			glPushName(i);
			glBegin(GL_POINTS);
			glVertex3d(v[0], v[1], v[2]);
			glEnd();
			glPopName();
		}
	}

	// エッジ描画
	if(drw & 0x0002){
		glEnable(GL_LIGHTING);
		glColor3d(0.5, 0.9, 0.9);
		glLineWidth(4.0);//エッジの太さを変更(元1.0) 海老沢変更
		glBegin(GL_LINES);
		for(int i = 0; i < m_iNumEdge; ++i){
			const rxEdge &e = m_poly.edges[i];
			glm::vec3 v0 = m_vCurPos[e.v[0]];
			glm::vec3 v1 = m_vCurPos[e.v[1]];
			glVertex3d(v0[0], v0[1], v0[2]);
			glVertex3d(v1[0], v1[1], v1[2]);
		}
		glEnd();
	}
}


/*!
 * 頂点を追加
 * @param[in] pos 頂点位置
 * @param[out] mass 頂点質量
 */
void ElasticPBD::AddVertex(const glm::vec3 &pos, float mass)
{
	m_poly.vertices.push_back(pos);
	m_vCurPos.push_back(pos);
	m_vNewPos.push_back(pos);
	m_vMass.push_back(mass);
	m_vVel.push_back(glm::vec3(0.0));
	m_vAcc.push_back(glm::vec3(0.0));//海老沢 加速度を追加
	m_vFix.push_back(false);
	m_iNumVertices++;

	//海老沢追加
	m_LamdaSS.push_back(glm::vec3(0.0));//ラグランジュ乗数の初期値を0にする
	//cout << glm::to_string(pos) << endl;
}

/*!
 * エッジを追加
 *  - ポリゴンとの関係などを計算する場合はMakeEdge関数を使うこと
 *  - ポリゴンがない1次元弾性用
 * @param[in] v0,v1 エッジを構成する頂点番号
 */
void ElasticPBD::AddEdge(int v0, int v1)
{
	rxEdge e;
	e.v[0] = v0; e.v[1] = v1;
	m_poly.edges.push_back(e);

	//初期状態
	m_vLengths.push_back(glm::length(m_poly.vertices[e.v[0]] - m_poly.vertices[e.v[1]]));

	m_vInEdge.push_back(0);
	m_iNumEdge++;

	//海老沢追加 四元数の追加-----------------------------------------------------
	glm::vec3 p1 = m_vNewPos[v0];
	glm::vec3 p2 = m_vNewPos[v1];

	Eigen::Quaterniond qinit = Eigen::Quaterniond::FromTwoVectors(Eigen::Vector3d::UnitZ(), Vec_glmToEigen(p2 - p1));
	qinit.normalize();
	m_vQuat.push_back(glm::quat(qinit.w(), qinit.x(), qinit.y(), qinit.z()));
	//cout << "x:" << m_vQuat[0][0] << " y:" << m_vQuat[0][1] << " z:" << m_vQuat[0][2] <<" w:"<<m_vQuat[0][3] << endl;
	//----------------------------------------------------------------------------

	//接線の保存------------------------------------------------------------------
	//m_vTang.push_back(glm::vec3(0.0,1.0,0.0));
	m_vTang.push_back(p2-p1);
	//----------------------------------------------------------------------------

	m_LamdaBT.push_back(glm::quat(0.f, glm::vec3(0.f)));
	m_Fss.push_back(glm::vec3(0.0));

	m_NewKss.push_back(PBD_ks);
	m_NewKbt.push_back(PBD_kbt);
}

//ダルボーベクトルを追加(ダルボーベクトル配列[i]はセグメントiとi+1
void ElasticPBD::AddDarbouxVector(int i) {
	glm::quat Darboux = glm::conjugate(m_vQuat[i-1]) * m_vQuat[i];
	Darboux = glm::normalize(Darboux);
	m_eOrgOmega.push_back(Darboux);//i=1から始まり、i-1番目に代入
	//cout << "m_eOrgOmega" << i << " " << "w:" << Darboux.w << " x:" << Darboux.x << " y:" << Darboux.y << " z:" << Darboux.z << endl;
}

//! 三角体の面積計算(このバージョンでは使わないが面積を保存するような拘束条件もあるので念のため残しておく)
static inline float calArea(const glm::vec3 &p0, const glm::vec3 &p1, const glm::vec3 &p2)
{
	return (1.0f/2.0f)*glm::length(glm::cross(p1-p0, p2-p0));
}

static inline long calKey(vector<int> v, int n = 10000)
{
	if(v[0] > v[1]) RX_SWAP(v[0], v[1]);
	return v[0]+v[1]*n;
}

/*!
 * nの頂点を持つ1次元弾性体生成
 * @param[in] c1,c2 2端点座標
 * @param[in] n エッジ分割数(頂点数ではないので注意.頂点数はn+1になる)
 */
void ElasticPBD::GenerateStrand(glm::vec3 c1, glm::vec3 c2, int n)
{
	Clear();

	glm::vec3 dir = c2-c1;

	//cout << endl << "Init Straight Rod----------------------------" << endl;

	// 頂点座標生成
	double mass = 0.01;// 各頂点の質量
	for(int i = 0; i < n+1; ++i){
		//直線のシーン
		AddVertex(c1+dir*(static_cast<float>(i)/static_cast<float>(n)), mass);//海老沢変更
	}
	for(int i = 0; i < n; ++i){
		AddEdge(i, i+1);//AddEdge内で四元数を追加
		if (i > 0)AddDarbouxVector(i);
	}
}

//重心が中心に来るような徐々に大きくなる螺旋型
void ElasticPBD::GenerateCenterSpiral(glm::vec3 c1, glm::vec3 c2, int n)
{
	Clear();

	glm::vec3 dir = c2 - c1;

	cout << endl << "Init Center Spiral----------------------------" << endl;

	// 頂点座標生成
	double mass = 1.f;// 各頂点の質量
	for (int i = 0; i < n + 1; ++i) {
		//螺旋形の実装
		float x, y, z, rad, pi;
		pi = 3.141592653589;
		y = c1.y + dir.y * (static_cast<float>(i) / static_cast<float>(n));
		//rad = i;//大きいサイズ感
		rad = 0.01 * i; //小さいサイズ感
		x = c1.x + rad * cos(i * pi / 5);
		z = c1.z + rad * sin(i * pi / 5);
		AddVertex(glm::vec3(x, y, z), mass);
	}
	for (int i = 0; i < n; ++i) {
		AddEdge(i, i + 1);//AddEdge内で四元数を追加
		if (i > 0)AddDarbouxVector(i);
	}
}

//一般的な螺旋型
void ElasticPBD::GenerateNaturalSpiral(glm::vec3 c1, glm::vec3 c2, int n)
{
	Clear();

	glm::vec3 dir = c2 - c1;

	cout << endl << "Init Natural Spiral----------------------------" << endl;

	// 頂点座標生成
	double mass = 1.f;// 各頂点の質量
	for (int i = 0; i < n + 1; ++i) {
		//螺旋形の実装
		float x, y, z, rad, pi;
		pi = 3.141592653589;
		y = c1.y + dir.y * (static_cast<float>(i) / static_cast<float>(n));
		//rad = 15;//大きいサイズ感
		rad = 0.15;//小さいサイズ感
		x = c1.x + rad * cos(i * pi / 6);
		z = c1.z + rad * sin(i * pi / 6);
		AddVertex(glm::vec3(x, y, z), mass);
	}
	for (int i = 0; i < n; ++i) {
		AddEdge(i, i + 1);//AddEdge内で四元数を追加
		if (i > 0)AddDarbouxVector(i);
	}
}

//実験用の関数で任意の形を設定(ここでは、ジグザグの形を設定
void ElasticPBD::GenerateExampleRod(glm::vec3 c1, glm::vec3 c2, int n)
{
	Clear();

	glm::vec3 dir = c2 - c1;

	cout << endl << "Init Example Rod----------------------------" << endl;

	// 頂点座標生成
	double mass = 1.f;// 各頂点の質量
	for (int i = 0; i < n + 1; ++i) {
		float x, y, z;
		glm::vec3 edge;
		//直線のシーン
		if (i % 2 == 1) {
			//x = z = 7.5;//大きいサイズ感
			x = c1.x + 0.075;//小さいサイズ感
			z = c1.z + 0.075;
		}
		else {
			x = c1.x;
			z = c1.z;
		}
		//edge = glm::vec3(x, -7.5 * i, z);//大きいサイズ感
		edge=glm::vec3(x,-0.075*i,z);//小さいサイズ感
		AddVertex(c1 + edge, mass);//海老沢変更
	}
	for (int i = 0; i < n; ++i) {
		AddEdge(i, i + 1);//AddEdge内で四元数を追加
		if (i > 0)AddDarbouxVector(i);
	}
}


/*!
 * 固定頂点を設定(位置変更も含む)
 * @param[in] i 頂点インデックス
 * @param[in] pos 固定位置
 */
void ElasticPBD::FixVertex(int i, const glm::vec3 &pos)
{
	m_vCurPos[i] = pos;
	m_vNewPos[i] = pos;
	m_vFix[i] = true;
}

/*!
 * 固定頂点を設定
 * @param[in] i 頂点インデックス
 */
void ElasticPBD::FixVertex(int i)
{
	m_vFix[i] = true;
}


/*!
 * 頂点の固定を解除
 * @param[in] i 頂点インデックス
 */
void ElasticPBD::UnFixVertex(int i)
{
	m_vFix[i] = false;
}

/*!
 * 全頂点の固定を解除
 */
void ElasticPBD::UnFixAllVertex(void)
{
	for(int i = 0; i < m_iNumVertices; ++i)	m_vFix[i] = false;
}

/*!
 * 頂点選択(レイと頂点(球)の交差判定)
 * @param[in]  ray_origin,ray_dir レイ(光線)の原点と方向ベクトル
 * @param[out] t 交差があったときの原点から交差点までの距離(媒介変数の値)
 * @param[in]  rad 球の半径(これを大きくすると頂点からマウスクリック位置が多少離れていても選択されるようになる)
 * @return 交差していればその頂点番号，交差していなければ-1を返す
 */
int ElasticPBD::IntersectRay(const glm::vec3 &ray_origin, const glm::vec3 &ray_dir, float &t, float rad)
{
	int v = -1;
	float min_t = 1.0e6;
	float rad2 = rad*rad;
	float a = glm::length2(ray_dir);
	if(a < 1.0e-6) return -1;

	for(int i = 0; i < m_iNumVertices; ++i){
		glm::vec3 cen = m_vCurPos[i];
		glm::vec3 s = ray_origin-cen;
		float b = 2.0f*glm::dot(s, ray_dir);
		float c = glm::length2(s)-rad2;

		float D = b*b-4.0f*a*c;
		if(D < 0.0f) continue;	// 交差なし

		float t0 = (-b-sqrt(D))/(2.0*a);
		float t1 = (-b+sqrt(D))/(2.0*a);
		if(t0 > 0.0 && t1 > 0.0 && t0 < min_t){	// 2交点がある場合
			v = i; min_t = t0;
		}
		else if(t0 < 0.0 && t1 > 0.0 && t1 < min_t){	// 1交点のみの場合(光線の始点が球内部にある)
			v = i; min_t = t1;
			v = i; min_t = t1;
		}
	}
	if(v != -1) t = glm::length(ray_origin-m_vCurPos[v]);
	return v;
}



/*!
 * 外力
 *  - 重力と境界壁からの力の影響
 */
void ElasticPBD::calExternalForces(float dt)
{

	for (int i = 0; i < m_iNumVertices; i++) {
		if (m_vFix[i]) continue;
		float mass = m_vMass[i];

		//位置ベース法に基づき、速度の更新(既に位置の差をタイムステップで割ったものが、m_vVelに代入図済み)
		m_vVel[i] += m_v3Gravity * dt + glm::vec3(m_fWind, 0.0f, 0.0f);

		//SPH法から加速度を設定
		//m_vVel[i] += m_vAcc[i] * dt + glm::vec3(m_fWind, 0.0f, 0.0f);
		//m_fWind = 1.5f*0.001;
		//重力を質点にかけない
		//m_vVel[i] += glm::vec3(m_fWind, 0.0f, 0.f);

		m_vNewPos[i] = m_vCurPos[i]+m_vVel[i]*dt;
	}
}

/*!
* 衝突判定
* @param[in] dt タイムステップ幅
*/
void ElasticPBD::genCollConstraints(float dt)
{
	// 境界壁の影響
	float res = 0.9;	// 反発係数
	for(int i = 0; i < m_iNumVertices; ++i){
		if(m_vFix[i]) continue;
		glm::vec3 &p = m_vCurPos[i];
		glm::vec3 &np = m_vNewPos[i];
		glm::vec3 &v = m_vVel[i];
		if(np[0] < m_v3Min[0] || np[0] > m_v3Max[0]){
			np[0] = p[0]-v[0]*dt*res;
			np[1] = p[1];
			np[2] = p[2];
		}
		if(np[1] < m_v3Min[1] || np[1] > m_v3Max[1]){
			np[1] = p[1]-v[1]*dt*res;
			np[0] = p[0];
			np[2] = p[2];
		}
		if(np[2] < m_v3Min[2] || np[2] > m_v3Max[2]){
			np[2] = p[2]-v[2]*dt*res;
			np[0] = p[0];
			np[1] = p[1];
		}
		clamp(m_vNewPos[i]);
	}

	// 他のオブジェクトとの衝突
	if(m_fpCollision != 0){
		for(int i = 0; i < m_iNumVertices; ++i){
			if(m_vFix[i]) continue;
			glm::vec3 &p = m_vCurPos[i];
			glm::vec3 &np = m_vNewPos[i];
			glm::vec3 &v = m_vVel[i];
			m_fpCollision(p, np, v, m_iObjectNo);
		}
	}
}


/*!
 * 速度と位置の更新
 *  - 新しい位置と現在の位置座標から速度を算出
 * @param[in] dt タイムステップ幅
 */
void ElasticPBD::integrate(float dt)
{
	float dt1 = 1.0f/dt;
	for(int i = 0; i < m_iNumVertices; ++i){
		m_vVel[i] = (m_vNewPos[i]-m_vCurPos[i])*dt1;
		m_vCurPos[i] = m_vNewPos[i];
	}
}

//海老沢追加
//SPH法から加速度を設定
//start:
void ElasticPBD::SetVertexAcc(vector<float> data,int start, int num) {
	if (num > m_iNumVertices) return;

	for (int i = 0; i < num; ++i) {
		m_vAcc[i][0] = data[DIM * (start + i)];
		m_vAcc[i][1] = data[DIM * (start + i) + 1];
		m_vAcc[i][2] = data[DIM * (start + i) + 2];
	}
	return;
}

//海老沢追加
//SPH法から位置を設定
void ElasticPBD::SetVertexPos(vector<float> data, int start, int num) {
	if (num > m_iNumVertices) return;

	for (int i = 0; i < num; i++) {
		m_vNewPos[i].x = data[DIM * (start + i)];
		m_vNewPos[i].y = data[DIM * (start + i) + 1];
		m_vNewPos[i].z = data[DIM * (start + i) + 2];
	}
}

void ElasticPBD::SetVertexCurPos(vector<float> data, int start, int num) {
	if (num > m_iNumVertices) return;

	for (int i = 0; i < num; i++) {
		m_vCurPos[i].x = data[DIM * (start + i)];
		m_vCurPos[i].y = data[DIM * (start + i) + 1];
		m_vCurPos[i].z = data[DIM * (start + i) + 2];
	}
}

//全ての頂点を出力
void ElasticPBD::PrintAllVertices(void) {
	for (int i = 0; i < m_iNumVertices; i++) {
		cout << "Number " << i << " pos" << glm::to_string(m_vCurPos[i]) << endl;
	}
}

bool ElasticPBD::judgePos(void) {
	float fixPos = m_vCurPos[0].y;
	float lastPos = m_vCurPos[m_iNumVertices - 1].y;
	return fixPos > lastPos;
}

//毛髪を反転させる
void ElasticPBD::ArrayReverse(void) {
	std::reverse(m_vMass.begin(), m_vMass.end());//質量
	std::reverse(m_vLengths.begin(), m_vLengths.end());//基準長
	std::reverse(m_NewKss.begin(), m_NewKss.end());//伸び剛性
	std::reverse(m_NewKbt.begin(), m_NewKbt.end());//曲げ剛性
	std::reverse(m_vQuat.begin(), m_vQuat.end());//姿勢
	std::reverse(m_eOrgOmega.begin(), m_eOrgOmega.end());//基準ダルボーベクトル
}
