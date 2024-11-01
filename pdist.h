/*!
  @file pdist.h

  @brief ���q�z�u�����p�̊֐��Q

  @author Makoto Fujisawa
  @date 2020-06,2022-07
*/

#ifndef _PARTICLE_DISTRIBUTION_H_
#define _PARTICLE_DISTRIBUTION_H_

#include <vector>
#include <deque>

#include "utils.h"

#include "mc.h"
#include "particle_on_surf.h"


//! ���q�z�u�p�\����
struct rxPdist
{
	std::vector<glm::vec3> pdist;	//!< �ǉ����q�z�u
	std::deque<int> steps;	//!< �ǉ��X�e�b�v(��������)
	glm::vec3 vel;			//!< �ǉ����q�̏������x
	int attr;			//!< �ǉ����q�̑���
};

/*!
 * �{�b�N�X�`����ɕ��ׂ�ꂽ���q�z�u�̍쐬
 * @param[in] cen �{�b�N�X�̒��S���W
 * @param[in] ext �{�b�N�X�̑傫��(�e�ӂ̒�����1/2)
 * @param[in] vel �������x
 * @param[in] spacing ���q�z�u�Ԋu
 * @param[in] attr ���q����(�����l0)
 * @param[in] step �ǉ��X�e�b�v��(�����l0)
 */
inline static rxPdist MakeBox(glm::vec3 cen, glm::vec3 ext, glm::vec3 vel, float spacing, int attr = 0, int step = 0)
{
	rxPdist pd;

	// �z�u�����ʒu(���_)�ƃ����_���ȗh�炬�ʂ̐ݒ�
	int sx = (int)(ext[0]/spacing)-1;
	int sy = (int)(ext[1]/spacing)-1;
	int sz = (int)(ext[2]/spacing)-1;
	srand((unsigned)time(NULL));
	float jitter = spacing*0.001f;

	// 3�����i�q��ɗ��q��z�u
	int count = 0;
	for(int z = -sz; z <= sz; ++z){
		for(int y = -sy; y <= sy; ++y){
			for(int x = -sx; x <= sx; ++x){
				glm::vec3 p = cen+glm::vec3(x, y, z)*spacing;
				for(int l = 0; l < 3; ++l) p[l] += (RX_FRAND()*2.0f-1.0f)*jitter;
				pd.pdist.push_back(p);
				count++;
			}
		}
	}

	// ���q���x�E�����̐ݒ�
	pd.vel = vel;
	pd.attr = attr;

	// �ǉ��^�C�~���O(�X�e�b�v��)
	pd.steps.push_back(step);

	return pd;
}


/*!
 * �{�b�N�X�`����ɕ��ׂ�ꂽ���q�z�u�̍쐬
 * @param[in] cen �{�b�N�X�̒��S���W
 * @param[in] ext �{�b�N�X�̑傫��(�e�ӂ̒�����1/2)
 * @param[in] vel �������x
 * @param[in] spacing ���q�z�u�Ԋu
 * @param[in] attr ���q����(�����l0)
 * @param[in] step �ǉ��X�e�b�v��(�����l0)
 */
inline static rxPdist MakeOne(glm::vec3 cen, glm::vec3 vel, int attr = 0, int step = 0)
{
	rxPdist pd;
	pd.pdist.push_back(cen);

	// ���q���x�E�����̐ݒ�
	pd.vel = vel;
	pd.attr = attr;

	// �ǉ��^�C�~���O(�X�e�b�v��)
	pd.steps.push_back(step);

	return pd;
}


/*!
 * ����̃{�b�N�X�`����ɕ��ׂ�ꂽ���q�z�u�̍쐬
 * @param[in] cen �{�b�N�X�̒��S���W
 * @param[in] ext0 �{�b�N�X�̓����̑傫��(�e�ӂ̒�����1/2)
 * @param[in] ext1 �{�b�N�X�̊O���̑傫��(�e�ӂ̒�����1/2)
 * @param[in] vel �������x
 * @param[in] spacing ���q�z�u�Ԋu
 * @param[in] attr ���q����(�����l0)
 * @param[in] step �ǉ��X�e�b�v��(�����l0)
 */
inline static rxPdist MakeHollowBox(glm::vec3 cen, glm::vec3 ext0, glm::vec3 ext1, glm::vec3 vel, float spacing, int attr = 0, int step = 0)
{
	rxPdist pd;

	// �z�u�����ʒu(���_)�ƃ����_���ȗh�炬�ʂ̐ݒ�
	int sx = (int)(ext1[0]/spacing)-1;
	int sy = (int)(ext1[1]/spacing)-1;
	int sz = (int)(ext1[2]/spacing)-1;
	srand((unsigned)time(NULL));
	float jitter = spacing*0.001f;

	// 3�����i�q��ɗ��q��z�u
	int count = 0;
	for(int z = -sz; z <= sz; ++z){
		for(int y = -sy; y <= sy; ++y){
			for(int x = -sx; x <= sx; ++x){
				glm::vec3 p = glm::vec3(x, y, z)*spacing;
				if(fabs(p[0]) < ext0[0] && fabs(p[1]) < ext0[1] && fabs(p[2]) < ext0[2]) continue;

				p += cen;
				for(int l = 0; l < 3; ++l) p[l] += (RX_FRAND()*2.0f-1.0f)*jitter;
				pd.pdist.push_back(p);
				count++;
			}
		}
	}

	// ���q���x�E�����̐ݒ�
	pd.vel = vel;
	pd.attr = attr;

	// �ǉ��^�C�~���O(�X�e�b�v��)
	pd.steps.push_back(step);

	return pd;
}


/*!
 * ���`����ɕ��ׂ�ꂽ���q�z�u�̍쐬
 * @param[in] cen ���̒��S���W
 * @param[in] ext ���̔��a
 * @param[in] vel �������x
 * @param[in] spacing ���q�z�u�Ԋu
 * @param[in] attr ���q����(�����l0)
 * @param[in] step �ǉ��X�e�b�v��(�����l0)
 */
inline static rxPdist MakeSphere(glm::vec3 cen, float r, glm::vec3 vel, float spacing, int attr, int step)
{
	rxPdist pd;

	// �����_���ȗh�炬�ʂ̐ݒ�
	srand((unsigned)time(NULL));
	float jitter = spacing*0.01f;

	// �O���b�h����̔��a
	int sr = (int)(r/spacing)+1;

	// 3�����i�q��(����)�ɗ��q��z�u
	int count = 0;
	for(int z = -sr; z <= sr; ++z){
		for(int y = -sr; y <= sr; ++y){
			for(int x = -sr; x <= sr; ++x){
				float dx[3];
				dx[0] = x*spacing;
				dx[1] = y*spacing;
				dx[2] = z*spacing;
				float l = sqrtf(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
				if((l <= r)) {
					glm::vec3 p(0.0);
					for(int j = 0; j < 3; ++j){
						p[j] = cen[j]+dx[j]+(RX_FRAND()*2.0f-1.0f)*jitter;
					}
					pd.pdist.push_back(p);
					count++;
				}
			}
		}
	}


	// ���q���x�E�����̐ݒ�
	pd.vel = vel;
	pd.attr = attr;

	// �ǉ��^�C�~���O(�X�e�b�v��)
	pd.steps.push_back(step);

	return pd;
}

/*!
 * �t�̗��q�̗������C���̐ݒ�
 * @param[in] cen ���̒��S���W
 * @param[in] ext ���̔��a
 * @param[in] vel �������x
 * @param[in] spacing ���q�z�u�Ԋu
 * @param[in] attr ���q����(�����l0)
 * @param[in] step �ǉ��X�e�b�v��(�����l0)
 */
inline static rxPdist MakeInletLine(glm::vec3 pos1, glm::vec3 pos2, glm::vec3 vel, glm::vec3 up, int accum, int span, float spacing, int attr, int steps)
{
	rxPdist pd;

	glm::vec3 rel;	// �[�_�Ԃ̑��Έʒu�x�N�g��
	float l = 0.0;		// �[�_�Ԃ̋���
	for(int i = 0; i < 3; ++i){
		rel[i] = pos2[i]-pos1[i];
		l += rel[i]*rel[i];
	}
	l = sqrt(l);

	int n = l/(spacing);	// ���ׂ闱�q�̐�
	if(!n) return pd;

	// �����_���ȗh�炬�ʂ̐ݒ�
	srand((unsigned)time(NULL));
	float jitter = spacing*0.01f;

	// 3�����i�q��(����)�ɗ��q��z�u
	int count = 0;
	for(int j = 0; j < accum; ++j){
		for(int i = 0; i < n; ++i){
			glm::vec3 p;
			for(int k = 0; k < 3; ++k){
				p[k] = pos1[k]+rel[k]*(i+0.5)/n+(RX_FRAND()*2.0f-1.0f)*jitter;
			}
			pd.pdist.push_back(p);
			count++;
		}
		pos1 += up*spacing;
	}


	// ���q���x�E�����̐ݒ�
	pd.vel = vel;
	pd.attr = attr;

	// �ǉ��^�C�~���O(�X�e�b�v��)
	for(int i = 0; i < steps; i += span){
		pd.steps.push_back(i);
	}

	return pd;
}



/*!
 * �\�ʗ��q����
 *  - �A�֐������͂Ƃ��āCMC�@�ŏ������q�𐶐��C���q�ԂɃo�l��z�u���Ĉʒu�𒲐�����
 *  - -�Q�l: A. Witkin and P. Heckbert, "Using particles to sample and control implicit surfaces", SIGGRAPH1994.
 * @param[in] func �A�֐��l�擾�֐��|�C���^(��1�����Ƀ��[�U�[�ϐ�,��2-4������3�������W��n���C�@���ƉA�֐��l(Vec4)��Ԃ��֐�)
 * @param[in] func_ptr func�̑�1�����ɓn�����[�U�ϐ�
 * @param[in] minp,maxp �v�Z�̈�
 * @param[in] rad ���q���a
 * @param[out] ppos �����������q���i�[����z��(���g�͂��ׂč폜����ĐV������񂪊i�[�����̂Œ���)
 * @param[in] dim ppos�ֈʒu���i�[����Ƃ��̎�����
 *
 */
static int GenerateParticlesOnSurf(glm::vec4 (*func)(float, float, float, void*), void* func_ptr, glm::vec3 minp, glm::vec3 maxp,
								   float rad, vector<float> &ppos, int dim = 3)
{
	int num_particles = 0;

	vector<glm::vec3> vrts, nrms;
	vector<int> tris;

	minp -= glm::vec3(6.0*rad);
	maxp += glm::vec3(6.0*rad);

	// �����z�u���q�v�Z:���b�V�������i�q�̌���
	float h = 2.0f*rad;
	float thr = 0.0f;
	int n[3];
	for(int i = 0; i < 3; ++i){
		n[i] = (int)((maxp[i]-minp[i])/h)+1;
	}

	// �����z�u���q�v�Z:Marching Cubes�œ��l�ʃ��b�V������
	//  -> �K�v�Ȃ̂͒��_���W�����Ȃ̂Ń��b�V���C���f�b�N�X�쐬���Ȃ��Ώ����͑����Ȃ肻��
	MCMesh mc;
	mc.CreateMesh(func, func_ptr, minp, h, n, thr, vrts, nrms, tris);

	// ���_��
	num_particles = (int)vrts.size();

	if(num_particles){
		cout << num_particles << " particles are generated for the boundary." << endl;

		ParticleOnSurf sp;

		// ���q�����z�u
		sp.Initialize(vrts, rad, minp, maxp, func, func_ptr);

		// ���q�ʒu�C�����̍ő唽���v�Z�񐔂Ƌ��e�덷
		int iter = 50;
		float eps = 1.0e-4;

		// ���q�ʒu�C��
		sp.Update(0.01, iter, eps);

		cout << iter << " iterations and the error value is " << eps << endl;

		// ���q���̎擾(�ʒu�C�����ɕK�v�ɉ����ė��q���폜����邱�Ƃ�����̂ŕK���Ď擾���Ă�������)
		num_particles = sp.GetNumParticles();

		// �ő̗��q�����i�[���郁�����̈�̊m��(nmax��-1�̏ꍇ�͐V�K�m�ہC�����łȂ���Ίm�ۍςݔz��ɒǉ�)
		if(!ppos.empty()) ppos.clear();
		ppos.resize(dim*num_particles, 0.0f);

		// �ő̗��q���̎��o��
		float* spos = sp.GetPositionArray();
		for(int i = 0; i < num_particles; ++i){
			ppos[dim*(i)+0] = spos[dim*i+0];
			ppos[dim*(i)+1] = spos[dim*i+1];
			ppos[dim*(i)+2] = spos[dim*i+2];
		}
	}

	return num_particles;
}



#endif // #ifdef _RX_PARTICLE_DISTRIBUTION_H_
