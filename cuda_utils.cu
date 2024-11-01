/*! 
@file cuda_utils.cu

@brief CUDA���ʃf�o�C�X�֐�
	- CUDA�f�o�C�X(__device__)�֐����L�q����t�@�C��
	- cu�t�@�C������C���N���[�h����
	- �r���h����͏��O����

@author Makoto Fujisawa
@date 2023-02
*/

#ifndef _CUDA_UTILS_CU_
#define _CUDA_UTILS_CU_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include <stdio.h>
#include <math.h>

#include "helper_math.h"
#include <math_constants.h>

#include "cuda_utils.h"

//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include <cstdio>
#include <GL/glew.h>
#if __APPLE__
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <GL/gl.h>
	#include <GL/glu.h>
#endif

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>


//-----------------------------------------------------------------------------
// device�֐� - �f�o�C�X(GPU)�Ŏ��s�E�f�o�C�X�֐�����̂݌Ăяo����
//-----------------------------------------------------------------------------
/*!
* [a,b]�ɃN�����v
* @param[in] x �N�����v���������l
* @param[in] a,b �N�����v���E
* @return �N�����v���ꂽ���l
*/
__device__ 
inline float CxClamp(float x, float a, float b)
{
    return max(a, min(b, x));
}
__device__ 
inline int CxClamp(int x, int a, int b)
{
    return max(a, min(b, x));
}

/*!
* a/b�̌v�Z���ʂ�؂�グ
* @param[in] a,b a/b
* @return �؂�グ�����Z����
*/
__device__
inline uint DivCeil(uint a, uint b)
{
    return (a % b != 0) ? (a / b + 1) : (a / b);
}

/*!
* �[������ for float3
* @param[in] v �l
*/
__device__
inline int CuIsZero(float3 v)
{
    if(fabsf(v.x) < 1.0e-10 && fabsf(v.y) < 1.0e-10 && fabsf(v.z) < 1.0e-10){
        return 1;
    }
    else{
        return 0;
    }
}

/*!
* �s��ƃx�N�g���̐�
* @param[in] m 3x3�s��(float3�̑傫��3�̔z��)
* @param[in] v 3D�x�N�g��
* @return �ς̌���
*/
__device__
inline float3 CuMulMV(float3 *m, float3 v)
{
    return make_float3(dot(m[0], v), dot(m[1], v), dot(m[2], v));
}



//-----------------------------------------------------------------------------
// �O���b�h
//-----------------------------------------------------------------------------
/*!
* 1D�C���f�b�N�X����3D�C���f�b�N�X�ւ̕ϊ�(�O���b�h���͔C��)
* @param[in] i 1D�C���f�b�N�X
* @param[in] gridSize �O���b�h��
* @return 3D�C���f�b�N�X
*/
__device__
inline int3 calcGridPos(int i, int3 ngrid)
{
	int3 gridPos;
	int w = i%(ngrid.x*ngrid.y);
	gridPos.x = w%ngrid.x;
	gridPos.y = w/ngrid.x;
	gridPos.z = i/(ngrid.x*ngrid.y);
	return gridPos;
}
/*!
* 3D�C���f�b�N�X����1D�C���f�b�N�X�ւ̕ϊ�(�O���b�h���͔C��)
* @param[in] p 3D�C���f�b�N�X
* @param[in] gridSize �O���b�h��
* @return 1D�C���f�b�N�X
*/
__device__
inline uint calcGridIndex(int3 p, int3 ngrid)
{
	p.x = min(p.x, ngrid.x-1);
	p.y = min(p.y, ngrid.y-1);
	p.z = min(p.z, ngrid.z-1);
	return (p.z*ngrid.x*ngrid.y)+(p.y*ngrid.x)+p.x;
}

//----------------------
// �Փ˔���
//----------------------
/*!
* AABB�Ɠ_�̋���
* @param[in] p �_���W
* @param[in] box_cen AABB�̒��S
* @param[in] box_ext AABB�̊e�ӂ̒�����1/2
* @param[out] cp AABB�\�ʂ̍ŋߖT�_
* @param[out] d ���s��AABB�̋���
* @param[out] n ��_�ɂ�����P�ʖ@���x�N�g��
*/
__device__
inline int distPointAABB(float3 p, float3 box_cen, float3 box_ext, float3& cp, float& d, float3& n)
{
    cp = p-box_cen;

    float3 tmp = fabs(cp)-box_ext;
    float res = ((tmp.x > tmp.y && tmp.x > tmp.z) ? tmp.x : (tmp.y > tmp.z ? tmp.y : tmp.z));

    float sgn = (res > 0.0) ? -1.0 : 1.0;

    int coli = 0;
    n = make_float3(0.0f);

    if(cp.x > box_ext.x){
        cp.x = box_ext.x;
        n.x -= 1.0;
        coli++;
    } else if(cp.x < -box_ext.x){
        cp.x = -box_ext.x;
        n.x += 1.0;
        coli++;
    }

    if(cp.y > box_ext.y){
        cp.y = box_ext.y;
        n.y -= 1.0;
        coli++;
    } else if(cp.y < -box_ext.y){
        cp.y = -box_ext.y;
        n.y += 1.0;
        coli++;
    }

    if(cp.z > box_ext.z){
        cp.z = box_ext.z;
        n.z -= 1.0;
        coli++;
    } else if(cp.z < -box_ext.z){
        cp.z = -box_ext.z;
        n.z += 1.0;
        coli++;
    }

    n = normalize(n);

    cp += box_cen;
    d = sgn*length(cp-p);

    return coli;
}

/*!
* �����Ɖ~�̌�������(2D, A��)
* @param[in] A,B �����̗��[�_���W
* @param[in] C �~�̒��S
* @param[in] r �~�̔��a
* @param[out] P ��_���W
* @return ��_��
*/
__device__ 
static int CuLineCircleIntersection(float2 A, float2 B, float2 C, float r, float2 P[2], float t[2])
{
	float rr = r*r;
	float2 AC = C-A;
	float2 BC = C-B;

	float2 v = B-A;
	float l = length(v);
	v /= l;

	float td = dot(v, AC);
	float2 D = A+td*v;
	float dd = dot(D-C, D-C);

	if(dd < rr){
		float dt = sqrtf(rr-dd);

		float da = rr-dot(AC, AC);
		float db = rr-dot(BC, BC);

		int inter = 0;
		float t1 = td-dt;
		float t2 = td+dt;
		if(t1 >= 0 && t1 <= l){
			P[inter] = A+t1*v;
			t[inter] = t1;
			inter++;
		}
		if(t2 >= 0 && t2 <= l){
			P[inter] = A+t2*v;
			t[inter] = t2;
			inter++;
		}

		return inter;
	}
	else{
		return 0;
	}
}


/*!
* AABB�Ƌ��̋���
* @param[in] spos �����S
* @param[in] r �����a
* @param[in] sgn
* @param[in] box_min,box_max AABB�ŏ��C�ő���W�l
* @param[out] cp AABB�\�ʂ̍ŋߖT�_
* @param[out] d ���s��AABB�̋���
* @param[out] n ��_�ɂ�����P�ʖ@���x�N�g��
*/
__device__
inline int collisionSphereAABB(float3 spos, float r, int sgn, float3 box_min, float3 box_max, float3 &cp, float &d, float3 &n)
{
	float3 dist_min;	// box_min�Ƃ̋���
	float3 dist_max;	// box_max�Ƃ̋���
	float d0 = 0.0f;
	float3 n0 = make_float3(0.0f, 0.0f, 0.0f);
	int bout = 0;
	int count = 0;

	// �e�����Ƃɍŏ��ƍő勫�E�O�ɂȂ��Ă��Ȃ������ׂ�
	if((dist_min.x = (spos.x-r)-box_min.x) < 0.0){ bout |= 0x0001; count++; d0 = dist_min.x; n0 = make_float3( 1.0,  0.0,  0.0);}
	if((dist_min.y = (spos.y-r)-box_min.y) < 0.0){ bout |= 0x0002; count++; d0 = dist_min.y; n0 = make_float3( 0.0,  1.0,  0.0);}
	if((dist_min.z = (spos.z-r)-box_min.z) < 0.0){ bout |= 0x0004; count++; d0 = dist_min.z; n0 = make_float3( 0.0,  0.0,  1.0);}
	if((dist_max.x = box_max.x-(spos.x+r)) < 0.0){ bout |= 0x0008; count++; d0 = dist_max.x; n0 = make_float3(-1.0,  0.0,  0.0);}
	if((dist_max.y = box_max.y-(spos.y+r)) < 0.0){ bout |= 0x0010; count++; d0 = dist_max.y; n0 = make_float3( 0.0, -1.0,  0.0);}
	if((dist_max.z = box_max.z-(spos.z+r)) < 0.0){ bout |= 0x0020; count++; d0 = dist_max.z; n0 = make_float3( 0.0,  0.0, -1.0);}

	// �����̓�(�S���ŋ��E��)
	if(bout == 0){
		float min_d = 1e10;
		if(dist_min.x < min_d){ min_d = dist_min.x; n = make_float3( 1.0,  0.0,  0.0); }
		if(dist_min.y < min_d){ min_d = dist_min.y; n = make_float3( 0.0,  1.0,  0.0); }
		if(dist_min.z < min_d){ min_d = dist_min.z; n = make_float3( 0.0,  0.0,  1.0); }

		if(dist_max.x < min_d){ min_d = dist_max.x; n = make_float3(-1.0,  0.0,  0.0); }
		if(dist_max.y < min_d){ min_d = dist_max.y; n = make_float3( 0.0, -1.0,  0.0); }
		if(dist_max.z < min_d){ min_d = dist_max.z; n = make_float3( 0.0,  0.0, -1.0); }

		d = (float)sgn*min_d;
		n *= (float)sgn;
		cp = spos+n*fabs(d);
		return 1;
	}

	// �����̊O
	// sgn = 1:���C-1:�I�u�W�F�N�g
	if(count == 1){
		// ���ʋߖT
		d = (float)sgn*d0;
		n = (float)sgn*n0;
		cp = spos+n*fabs(d);
	}
	else{
		// �G�b�W/�R�[�i�[�ߖT
		float3 x = make_float3(0.0f, 0.0f, 0.0f);
		if(bout & 0x0001) x.x =  dist_min.x;
		if(bout & 0x0002) x.y =  dist_min.y;
		if(bout & 0x0004) x.z =  dist_min.z;
		if(bout & 0x0008) x.x = -dist_max.x;
		if(bout & 0x0010) x.y = -dist_max.y;
		if(bout & 0x0020) x.z = -dist_max.z;

		d = length(x);
		n = normalize(x);

		d *= -(float)sgn;
		n *= -(float)sgn;

		cp = spos+n*fabs(d);

		float3 disp = make_float3(0.00001);
		//Random(disp, 0, 0.00001);
		disp = disp*n;
		cp += disp;
	}

	return 0;
}


/*!
* AABB�Ɠ_�̋���
* @param[in] p �_���W
* @param[in] box_cen AABB�̒��S
* @param[in] box_ext AABB�̊e�ӂ̒�����1/2
* @param[out] cp AABB�\�ʂ̍ŋߖT�_
* @param[out] d ���s��AABB�̋���
* @param[out] n ��_�ɂ�����P�ʖ@���x�N�g��
*/
__device__
inline int collisionPointAABB(float3 p, float3 box_cen, float3 box_ext, float3 &cp, float &d, float3 &n)
{
	cp = p-box_cen;

	float3 tmp = fabs(cp)-box_ext;
	float res = ((tmp.x > tmp.y && tmp.x > tmp.z) ? tmp.x : (tmp.y > tmp.z ? tmp.y : tmp.z));

	float sgn = (res > 0.0) ? -1.0 : 1.0;

	int coli = 0;
	n = make_float3(0.0f);

	if(cp.x > box_ext.x){
		cp.x = box_ext.x;
		n.x -= 1.0;
		coli++;
	}
	else if(cp.x < -box_ext.x){
		cp.x = -box_ext.x;
		n.x += 1.0;
		coli++;
	}

	if(cp.y > box_ext.y){
		cp.y = box_ext.y;
		n.y -= 1.0;
		coli++;
	}
	else if(cp.y < -box_ext.y){
		cp.y = -box_ext.y;
		n.y += 1.0;
		coli++;
	}

	if(cp.z > box_ext.z){
		cp.z = box_ext.z;
		n.z -= 1.0;
		coli++;
	}
	else if(cp.z < -box_ext.z){
		cp.z = -box_ext.z;
		n.z += 1.0;
		coli++;
	}

	n = normalize(n);

	//if(coli > 1){
	//	float3 disp;
	//	Random(disp, 0, 0.00001);
	//	disp = disp*n;
	//	cp += disp;
	//}

	cp += box_cen;
	d = sgn*length(cp-p);

	return 0;
}


/*!
* �_��BOX�̋���
* @param[in] p �_���W
* @param[in] box_cen BOX�̒��S
* @param[in] box_ext BOX�̊e�ӂ̒�����1/2
* @param[in] box_rot BOX�̕����s��(3x3��]�s��)
* @param[in] box_inv_rot BOX�̕����s��̋t�s��(3x3)
* @param[out] cp BOX�\�ʂ̍ŋߖT�_
* @param[out] d �_��BOX�̋���
* @param[out] n ��_�ɂ�����P�ʖ@���x�N�g��
*/
__device__
inline int collisionPointBox(float3 p, float3 box_cen, float3 box_ext, float3 box_rot[3], float3 box_inv_rot[3], float3 &cp, float &d, float3 &n)
{
	cp = p-box_cen;
	cp = CuMulMV(box_rot, cp);

	float3 tmp = fabs(cp)-box_ext;

	int coli = 0;
	n = make_float3(0.0f);

	if(tmp.x < 0.0 && tmp.y < 0.0 && tmp.z < 0.0){
		tmp = fabs(tmp);

		if(tmp.x <= tmp.y && tmp.x <= tmp.z){	// x���ʂɋ߂�
			if(cp.x > 0){
				cp.x = box_ext.x;
				n.x += 1.0;
			}
			else{
				cp.x = -box_ext.x;
				n.x -= 1.0;
			}
		}
		else if(tmp.y <= tmp.x && tmp.y <= tmp.z){ // y���ʂɋ߂�
			if(cp.y > 0){
				cp.y = box_ext.y;
				n.y += 1.0;
			}
			else{
				cp.y = -box_ext.y;
				n.y -= 1.0;
			}
		}
		else{ // z���ʂɋ߂�
			if(cp.z > 0){
				cp.z = box_ext.z;
				n.z += 1.0;
			}
			else{
				cp.z = -box_ext.z;
				n.z -= 1.0;
			}
		}

		coli++;
	}

	cp = CuMulMV(box_inv_rot, cp);
	n  = CuMulMV(box_inv_rot, n);

	n = normalize(n);
	cp += box_cen;

	float sgn = (coli) ? -1.0 : 1.0;
	d = sgn*(length(cp-p));

	return 0;
}

/*!
* �_�Ƌ��̋���
* @param[in] p �_���W
* @param[in] sphere_cen ���̒��S
* @param[in] sphere_rad ���̔��a
* @param[out] cp �_�Ƌ����S�����Ԑ����Ƌ��̌�_
* @param[out] d �_�Ƌ��\�ʂ̋���
* @param[out] n �����S����_�ւ̒P�ʃx�N�g��
*/
__device__
inline int collisionPointSphere(float3 p, float3 sphere_cen, float sphere_rad, float3 &cp, float &d, float3 &n)
{
	n = make_float3(0.0f);

	float3 l = p-sphere_cen;
	float ll = length(l);

	d = ll-sphere_rad;
	if(d < 0.0){
		n = normalize(p-sphere_cen);
		cp = sphere_cen+n*sphere_rad;
	}

	return 0;
}

/*!
* �_�ƕ��ʂ̋���
* @param[in] v  �_�̍��W
* @param[in] px ���ʏ�̓_
* @param[in] pn ���ʂ̖@��
* @return ����
*/
__device__ 
inline float distPointPlane(float3 v, float3 px, float3 pn)
{
	return dot((v-px), pn)/length(pn);
}

/*!
* �O�p�`�Ɠ_�̋����ƍŋߖT�_
* @param[in] v0,v1,v2	�O�p�`�̒��_
* @param[in] n			�O�p�`�̖@��
* @param[in] p			�_
* @return 
*/
__device__ 
inline int distPointTriangle(float3 v0, float3 v1, float3 v2, float3 n, float3 p, float &dist, float3 &p0)
{
	// �|���S�����܂ޕ��ʂƓ_�̋���
	float l = distPointPlane(p, v0, n);

	// ���ʂƂ̍ŋߖT�_���W
	float3 np = p-l*n;

	// �ߖT�_���O�p�`�����ǂ����̔���
	float3 n1 = cross((v0-p), (v1-p));
	float3 n2 = cross((v1-p), (v2-p));
	float3 n3 = cross((v2-p), (v0-p));

	if(dot(n1, n2) > 0 && dot(n2, n3) > 0){
		// �O�p�`��
		dist = l;
		p0 = np;
		return 1;
	}
	else{
		// �O�p�`�O
		return 0;
	}
}


/*!
* ���C/�����ƎO�p�`�̌���
* @param[in] P0,P1 ���C/�����̒[�_or���C��̓_
* @param[in] V0,V1,V2 �O�p�`�̒��_���W
* @param[out] I ��_���W
* @retval 1 ��_I�Ō��� 
* @retval 0 ��_�Ȃ�
* @retval 2 �O�p�`�̕��ʓ�
* @retval -1 �O�p�`��"degenerate"�ł���(�ʐς�0�C�܂�C�������_�ɂȂ��Ă���)
*/
inline __device__ 
int intersectSegmentTriangle(float3 P0, float3 P1, 
							 float3 V0, float3 V1, float3 V2, 
							 float3 &I, float3 &n, float rp = 0.01)
{
	// �O�p�`�̃G�b�W�x�N�g���Ɩ@��
	float3 u = V1-V0;		
	float3 v = V2-V0;			
	n = normalize(cross(u, v));
	if(CuIsZero(n)){
		return -1;	// �O�p�`��"degenerate"�ł���(�ʐς�0)
	}

	// ����
	float3 dir = P1-P0;
	float a = dot(n, P0-V0);
	float b = dot(n, dir);
	if(fabs(b) < 1e-10){	// �����ƎO�p�`���ʂ����s
		if(a == 0){
			return 2;	// ���������ʏ�
		}
		else{
			return 0;	// ��_�Ȃ�
		}
	}


	// ��_�v�Z

	// 2�[�_�����ꂼ��قȂ�ʂɂ��邩�ǂ����𔻒�
	float r = -a/b;
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

	// �����ƕ��ʂ̌�_
	I = P0+r*dir;

	// ��_���O�p�`���ɂ��邩�ǂ����̔���
	float uu, uv, vv, wu, wv, D;
	uu = dot(u, u);
	uv = dot(u, v);
	vv = dot(v, v);
	float3 w = I-V0;
	wu = dot(w, u);
	wv = dot(w, v);
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




/*!
* �����Ƌ��̌�������
* @param[in] s0,s1	�����̒[�_
* @param[in] sc,r   ���̒��S���W�Ɣ��a
* @param[out] d2 �����Ƃ̋����̓��
* @return ���������true
*/
__device__ 
inline bool segment_sphere(const float3 &s0, const float3 &s1, const float3 &sc, const float &r, float &d2)
{
	float3 v = s1-s0;
	float3 c = sc-s0;

	float vc = dot(v, c);
	if(vc < 0){		// ���̒��S�������̎n�_s0�̊O�ɂ���
		d2 = dot(c, c);
		return (d2 < r*r);	// �����S�Ǝn�_s0�̋����Ō�������
	}
	else{
		float v2 = dot(v, v);
		if(vc > v2){	// ���̒��S�������̏I�_s1�̊O�ɂ���
			d2 = dot(s1-sc, s1-sc);
			return (d2 < r*r);	// �����S�ƏI�_s1�̋����Ō�������
		}
		else{			// ����s0��s1�̊Ԃɂ���
			float3 a = (vc*v)/dot(v, v)-c;
			d2 = dot(a, a);
			return (d2 < r*r);	// �����Ƌ����S�̋����Ō�������
		}
	}
}

/*!
* ����(���܂ޒ���)�Ɠ_�̋���
* @param[in] v0,v1 �����̗��[�_���W
* @param[in] p �_�̍��W
* @return ����
*/
__device__ 
inline double segment_point_dist(const float3 &v0, const float3 &v1, const float3 &p)
{
	float3 v = normalize(v1-v0);
	float3 vp = p-v0;
	float3 vh = dot(vp, v)*v;
	return length(vp-vh);
}


/*!
* ����(���C,������)�Ƌ��̌�������
* @param[in] p,d ���C�̌��_�ƕ���
* @param[in] c,r ���̒��S�Ɣ��a
* @param[out] t1,t2 p�����_�܂ł̋���
* @return ��_��
*/
__device__
inline int ray_sphere(const float3 &p, const float3 &d, const float3 &sc, const float r, float &t1, float &t2)
{
	float3 q = p-sc;	// �����S���W�n�ł̌������_���W

	float a = dot(d, d);
	float b = 2*dot(q, d);
	float c = dot(q, q)-r*r;

	// ���ʎ�
	float D = b*b-4*a*c;

	if(D < 0.0){ // �����Ȃ�
		return 0;
	}
	else if(D < 1e-8){ // ��_��1
		t1 = -b/(2*a);
		t2 = -1;
		return 1;
	}
	else{ // ��_��2
		float sqrtD = sqrt(D);
		t1 = (-b-sqrtD)/(2*a);
		t2 = (-b+sqrtD)/(2*a);
		return 2;
	}

}
/*!
* �O�p�`�Ƌ��̌�������
* @param[in] v0,v1,v2	�O�p�`�̒��_
* @param[in] n			�O�p�`�̖@��
* @param[in] p			�ŋߖT�_
* @return 
*/
__device__
inline bool triangle_sphere(const float3 &v0, const float3 &v1, const float3 &v2, const float3 &n, 
							const float3 &c, const float &r, float &dist, float3 &ipoint)
{
	// �|���S�����܂ޕ��ʂƋ����S�̋���
	float d = dot(v0, n);
	float l = dot(n, c)-d;

	dist = l;
	if(l > r) return false;

	// ���ʂƂ̍ŋߖT�_���W
	float3 p = c-l*n;

	// �ߖT�_���O�p�`�����ǂ����̔���
	float3 n1 = cross((v0-c), (v1-c));
	float3 n2 = cross((v1-c), (v2-c));
	float3 n3 = cross((v2-c), (v0-c));

	ipoint = p;
	dist = l;
	if(dot(n1, n2) > 0 && dot(n2, n3) > 0){		// �O�p�`��
		return true;
	}
	else{		// �O�p�`�O
				// �O�p�`�̊e�G�b�W�Ƌ��̏Փ˔���
		for(int e = 0; e < 3; ++e){
			float3 va0 = (e == 0 ? v0 : (e == 1 ? v1 : v2));
			float3 va1 = (e == 0 ? v1 : (e == 1 ? v2 : v0));

			float t1, t2;
			int n = ray_sphere(va0, normalize(va1-va0), c, r, t1, t2);

			if(n){
				float le2 = dot(va1-va0, va1-va0);
				if((t1 >= 0.0 && t1*t1 < le2) || (t2 >= 0.0 && t2*t2 < le2)){
					return true;
				}
			}
		}
		return false;
	}
}


#endif // #ifndef _CUDA_UTILS_CU_