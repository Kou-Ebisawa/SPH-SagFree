/*! 
  @file sph.cu
	
  @brief CUDA : SPH�@

  @author Makoto Fujisawa
  @date 2023-02
*/


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

#include "sph_kernel.cu"

#include <thrust/device_vector.h>
#include <thrust/scan.h>
#include <thrust/sort.h>

//-----------------------------------------------------------------------------
// CUDA�֐�
//-----------------------------------------------------------------------------
extern "C"
{
/*!
 * �p�����[�^��GPU�֓]��
 * @param[in] hparams �z�X�g(CPU)�������Ɋi�[���ꂽ�p�����[�^
 */
void CuSetParameters(const SceneParameter* hparams)
{
	CUCHECK(cudaMemcpyToSymbol(params, hparams, sizeof(SceneParameter), 0, cudaMemcpyHostToDevice));
}

/*!
 * �X���b�h������u���b�N/�O���b�h���̌v�Z(�K�v�X���b�h��n�ȏ�ɂȂ�悤�ɐݒ�)
 * @param[in] n �K�v�X���b�h��
 * @param[out] block,grid �u���b�N��/�O���b�h��
 */
void CuCalGridN(int n, dim3& block, dim3& grid)
{
	// �X���b�h���̐ݒ�(n�ȏ�ɂȂ�悤�ɐݒ�)
	block = dim3(THREAD_NUM, 1, 1); // 1�u���b�N����X���b�h��
	grid = dim3((n+block.x-1)/block.x, 1, 1); // 1�O���b�h������̃u���b�N��
}

/*!
 * CUDA�J�[�l�����g���ĕ���v�Z:���q���x�̌v�Z
 *  - ���q�ʒu��cell.dSortedPos����擾����̂ň����Ƃ��ēn���K�v�Ȃ�
 * @param[out] ddens ���q���x(�f�o�C�X������)
 * @param[in] dvol ���q�̐�(�f�o�C�X������)
 * @param[in] n ���q��
 */
void CuSphDensity(float* drestdens,float* ddens, float* dvol,float* dmas, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxSphDensity<<<grid, block>>>(drestdens,ddens, dvol,dmas, n);	// �J�[�l�����s
	cudaThreadSynchronize();
}

/*!
 * CUDA�J�[�l�����g���ĕ���v�Z:���q���͒l�𖧓x����v�Z
 * @param[out] dpres ���q����(�f�o�C�X������)
 * @param[in] ddens ���q���x(�f�o�C�X������)
 * @param[in] n ���q��
 */
void CuSphPressure(float* drestdens,float* dpres, float* ddens, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxSphPressure<<<grid, block>>>(drestdens,dpres, ddens, n);	// �J�[�l�����s
	cudaThreadSynchronize();
}

/*!
* CUDA�J�[�l�����g���ĕ���v�Z:�e���q�̉Q�x�v�Z[Macklin2013]
*  - M. Macklin & M. M\"{u}ller, "Position Based Fluids", ACM ToG, 32(4), pp.104:1-104:12, 2013.
*  - �O���b�h�@�����̌��̎�@�� 
*    R. Fedkiw; J. Stam & H. Jensen, "Visual simulation of smoke", Proc. SIGGRAPH 2001, pp.15-22, 2001.
* @param[out] dvort �e���q�̉Q�x�x�N�g��(�f�o�C�X������)
* @param[in] dvel ���q���x�z��(�f�o�C�X������)
* @param[in] ddens ���q���x(�f�o�C�X������)
* @param[in] dvol ���q�̐�(�f�o�C�X������)
* @param[in] datt ���q����(0�ŗ���,1�ŋ��E)(�f�o�C�X������)
* @param[in] n ���q��
*/
void CuSphVorticity(float* dvort, float* dvel, float* ddens, float* dvol, int* datt, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxSphVorticity<<<grid, block>>>(dvort,dvel, ddens, dvol, datt, n);	// �J�[�l�����s
	cudaThreadSynchronize();
}

/*!
 * CUDA�J�[�l�����g���ĕ���v�Z:���q�ɓ�����(���͍�&�O�͍�)�̌v�Z
 * @param[out] dacc �e���q�ɓ�����(�����xdv/dt)(�f�o�C�X������)
 * @param[in] dvel ���q���x�z��(�f�o�C�X������)
 * @param[in] ddens ���q���x(�f�o�C�X������)
 * @param[in] dpres ���q����(�f�o�C�X������)
 * @param[in] dvort �e���q�̉Q�x�x�N�g��(�f�o�C�X������)
 * @param[in] dvol  ���q�̐�(�f�o�C�X������)
 * @param[in] datt ���q����(0�ŗ���,1�ŋ��E)(�f�o�C�X������)
 * @param[in] n ���q��
 */
void CuSphForces(float* drestdens,float* dacc, float* dvel, float* ddens, float* dpres, float* dvort, float* dvol,float* dmas, int* datt,float3 power,float* dfss, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxSphForces<<<grid, block>>>(drestdens,dacc, dvel, ddens, dpres, dvort, dvol,dmas, datt,power,dfss, n);	// �J�[�l�����s
	cudaThreadSynchronize();
}
/*!
* CUDA�J�[�l�����g���ĕ���v�Z:���q�ɓ�����(�S����)�̌v�Z[Becker2007]
*  - M. Becker & M. Teschner, "Weakly Compressible SPH for Free Surface Flows", Proc. SCA2007, pp.209-217, 2007.
* @param[inout] dacc �e���q�ɓ�����(�����xdv/dt)(�f�o�C�X������)
* @param[in] dvel ���q���x�z��(�f�o�C�X������)
* @param[in] ddens ���q���x(�f�o�C�X������)
* @param[in] dvol ���q�̐�(�f�o�C�X������)
* @param[in] datt ���q����(0�ŗ��̗��q�C����ȊO�ŋ��E���q)(�f�o�C�X������)
* @param[in] n ���q��
*/
void CuSphViscosityForces(float* drestdens,float* dacc, float* dvel, float* ddens, float* dvol,float* dmas, int* datt, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxSphViscosity<<<grid, block>>>(drestdens,dacc, dvel, ddens, dvol,dmas, datt, n);	// �J�[�l�����s
	cudaThreadSynchronize();
}

/*!
* CUDA�J�[�l�����g���ĕ���v�Z:XSPH�l�H�S���̌v�Z[Schechter2012]
*  - H. Schechter & R. Bridson, "Ghost SPH for animating water", ACM ToG, 31(4), pp.61:1-61:8, 2012.
*  - ��(�����x)�Ƃ��ĔS������������̂ł͂Ȃ��C���x�𒼐ڍX�V����
*  - ���̂̕����I�����Ƃ��Ă̔S���ł͂Ȃ��v�Z���萫�̂��߂̌v�Z�Ƃ����������悳����
* @param[inout] dvel ���q���x�z��(�f�o�C�X������)
* @param[in] ddens ���q���x(�f�o�C�X������)
* @param[in] dvol ���q�̐�(�f�o�C�X������)
* @param[in] datt ���q����(0�ŗ��̗��q�C����ȊO�ŋ��E���q)(�f�o�C�X������)
* @param[in] n ���q��
*/
void CuSphXSPHViscosity(float* dvel, float* ddens, float* dvol,float* dmas, int* datt, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxSphXSPHViscosity<<<grid, block>>>(dvel, ddens, dvol,dmas, datt, n);	// �J�[�l�����s
	cudaThreadSynchronize();
}


/*!
 * CUDA�J�[�l�����g���ĕ���v�Z:�����x�ɏ]���Ĉʒu�Ƒ��x���X�V
 * @param[inout] dpos ���q�ʒu�z��(�f�o�C�X������)
 * @param[inout] dvel ���q���x�z��(�f�o�C�X������)
 * @param[in] dacc �e���q�ɓ�����(�����x)���i�[�����z��(�f�o�C�X������)
 * @param[in] dvol ���q�̐�(�f�o�C�X������)
 * @param[in] datt ���q����(0�ŗ��̗��q�C����ȊO�ŋ��E���q)(�f�o�C�X������)
 * fix:�C�V��ǉ� �Œ�_��\��
 * @param[in] n ���q��
 */
void CuSphIntegrate(float* dpos, float* dvel, float* dacc, int* datt, int* dfix,int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxSphIntegrate<<<grid, block>>>(dpos, dvel, dacc, datt,dfix, n);	// �J�[�l�����s
	cudaThreadSynchronize();
}


/*!
* CUDA�J�[�l�����g���ĕ���v�Z:�����x�ɏ]���đ��x�݂̂��X�V
*  - XSPH�p
* @param[inout] dvel ���q���x�z��(�f�o�C�X������)
* @param[in] dacc �e���q�ɓ�����(�����x)���i�[�����z��(�f�o�C�X������)
* @param[in] dvol ���q�̐�(�f�o�C�X������)
* @param[in] datt ���q����(0�ŗ��̗��q�C����ȊO�ŋ��E���q)(�f�o�C�X������)
* @param[in] n ���q��
*/
void CuSphIntegrateV(float* dvel, float* dacc, int* datt, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxSphIntegrateVelocity<<<grid, block>>>(dvel, dacc, datt, n);	// �J�[�l�����s
	cudaThreadSynchronize();
}
/*!
* CUDA�J�[�l�����g���ĕ���v�Z:���x���]���Ĉʒu���X�V
*  - XSPH�p
* @param[inout] dpos ���q�ʒu�z��(�f�o�C�X������)
* @param[in] dvel ���q���x�z��(�f�o�C�X������)
* @param[in] dvol ���q�̐�(�f�o�C�X������)
* @param[in] datt ���q����(0�ŗ��̗��q�C����ȊO�ŋ��E���q)(�f�o�C�X������)
* @param[in] n ���q��
*/
void CuSphIntegrateP(float* dpos, float* dvel, int* datt,int*dfix, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxSphIntegratePosition<<<grid, block>>>(dpos, dvel, datt,dfix, n);	// �J�[�l�����s
	cudaThreadSynchronize();
}



/*!
* CUDA�J�[�l�����g���ĕ���v�Z:���q�̐ς̌v�Z
*  - ���q�ʒu��cell.dSortedPos����擾����̂ň����Ƃ��ēn���K�v�Ȃ�
* @param[out] dvol ���q�̐�(�f�o�C�X������)
* @param[in] datt ���q����(0�ŗ��̗��q�C����ȊO�ŋ��E���q)(�f�o�C�X������)
* @param[in] n �������闱�q��(offset����̑��ΓI�Ȉʒu)
* @param[in] v ���̗��q�̏ꍇ�̗��q�̐ϒl
*/
void CuSphCalVolume(float* dvol, int *datt, int n, float v)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxSphCalVolume<<<grid, block>>>(dvol, datt, n, v);	// �J�[�l�����s
	cudaThreadSynchronize();
}

/*!
* CUDA�J�[�l�����g���ĕ���v�Z:���q���x�̌v�Z
*  - ���q�ʒu��cell.dSortedPos����擾����̂ň����Ƃ��ēn���K�v�Ȃ�
* @param[out] dF ���x�l���i�[����O���b�h�Z���z��(�f�o�C�X������)
* @param[in] dvol ���q�̐�(�f�o�C�X������)
* @param[in] datt ���q����(0�ŗ���,1�ŋ��E)(�f�o�C�X������)
* @param[in] n ���q��
* @param[in] gnum �O���b�h��
* @param[in] gmin �O���b�h�ŏ����W
* @param[in] glen �O���b�h��
*/
void CuSphDensityInGrid(float* dF, float* dvol, int* datt, int n, int3 gnum, float3 gmin, float3 glen)
{
	// ���O���b�h�Z����=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	int numcell = gnum.x*gnum.y*gnum.z;
	dim3 block, grid;
	CuCalGridN(numcell, block, grid);	
	CxSphDensityAtCell<<<grid, block>>>(dF, dvol, datt, n, gnum, gmin, glen);	// �J�[�l�����s
	cudaThreadSynchronize();
}


/*!
* CUDA�J�[�l�����g���ĕ���v�Z:���q�̕`��F�𖧓x����v�Z
* @param[out] dcol ���q�F�z��(�f�o�C�X������)
* @param[in] dval  ���q�����ʔz��(�f�o�C�X������)
* @param[in] n ���q��
* @param[in] c1,c2 �����ʂ��ŏ�,�ő�̂Ƃ��̐F(�Ԃ̐F�͐��`��Ԃŋ��߂���)
* @param[in] range x�v�f�ɍŏ��l�Cy�v�f�ɍő�l
*/
void CuColorScalar(float* dcol, int* datt, float* dval, int n, float3 c1, float3 c2, float2 range)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxColorScalar<<<grid, block>>>(dcol, datt, dval, n, c1, c2, range);	// �J�[�l�����s
	cudaThreadSynchronize();
}
void CuColorVector(float* dcol, int* datt, float* dval, int n, float3 c1, float3 c2, float2 range)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxColorVector<<<grid, block>>>(dcol, datt, dval, n, c1, c2, range);	// �J�[�l�����s
	cudaThreadSynchronize();
}
/*!
* CUDA�J�[�l�����g���ĕ���v�Z:���q�̕`��F�ݒ� - ���̐F
* @param[out] dcol ���q�F�z��(�f�o�C�X������)
* @param[in] col �`��F
* @param[in] n ���q��
*/
void CuColorConstant(float* dcol, int* datt, float3 col, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxColorConstant<<<grid, block>>>(dcol, datt, col, n);	// �J�[�l�����s
	cudaThreadSynchronize();
}


/*!
 * �e���q�̃O���b�h�n�b�V���l�v�Z(�ߖT�T���p)
 * @param[out] dhash �e���q�̃O���b�h�n�b�V���l���i�[�����z��
 * @param[out] dsortedidx �e���q�̃C���f�b�N�X���i�[�����z��(�ォ��n�b�V���l�Ń\�[�g����� -> �����_�ł͂܂��\�[�h�ς݂ł͂Ȃ�)
 * @param[in] dpos ���q�ʒu���i�[�����z��
 * @param[in] n ���q��
 */
void CuCalcHash(uint* dhash, uint* dindex, float* dpos, int n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxCalcHash<<<grid, block>>>(dhash, dindex, dpos, n);	// �J�[�l�����s
	cudaThreadSynchronize();
}

/*!
 * thrust::sort_by_key�ɂ��n�b�V���l�Ɋ�Â��\�[�g
 * @param[in] dhash �n�b�V���l
 * @param[in] dindex �C���f�b�N�X(�p�[�e�B�N���C�|���S���Ȃ�)
 * @param[in] n �f�[�^��
 */
void CuSort(unsigned int* dhash, uint* dindex, uint n)
{
	thrust::sort_by_key(thrust::device_ptr<unsigned int>(dhash),
					    thrust::device_ptr<unsigned int>(dhash+n),
					    thrust::device_ptr<unsigned int>(dindex));
	cudaThreadSynchronize();
}

/*!
 * �p�[�e�B�N���z����\�[�g���ꂽ���Ԃɕ��ёւ��C�e�Z���̎n�܂�ƏI���̃C���f�b�N�X������
 * @param[in] cell �ߖT�T���p�O���b�h�f�[�^
 * @param[in] dpos ���q�ʒu
 * @param[in] dvel ���q���x
 */
void CuReorderDataAndFindCellStart(Cell cell, float* dpos, float* dvel, uint n)
{
	dim3 block, grid;
	CuCalGridN(n, block, grid);	// ���q��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z

	// �Z���X�^�[�g�ʒu�z��̏�����
	CUCHECK(cudaMemset(cell.dCellStart, 0xffffffff, cell.uNumCells*sizeof(uint)));

	// �V�F�A�[�h�������̃T�C�Y
	uint smemSize = sizeof(uint)*(THREAD_NUM+1);

	// �J�[�l�����s
	CxReorderDataAndFindCellStartD<<<grid, block, smemSize>>>(cell, dpos, dvel, n);
	cudaThreadSynchronize();
}

//�C�V��ǉ�-----------------------------------------------------------------------------------------------------
//XPBD�̐L�сE����f�C�Ȃ��E�˂��ꐧ��̏���
//dpos:�ʒu
//dmas:����
//dlen:���
//dkss:�L�э���
//dkbt:�Ȃ�����
//dquat:�p��(�l����)
//domega:��_���{�[�x�N�g��
//dlamb_ss:XPBD�̐L�сE����f����ɗp�����
//dlamb_bt:XPBD�̋Ȃ��E�˂��ꐧ��ɗp�����
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//dt:�^�C���X�e�b�v
//n:���q��
//iter:������
//example_flag:�`��ɂ���āC�������ꕔ�ς���
void CuXPBDConstraint(float* dpos,float* dcurpos,float* dmas, float* dlen, float* dkss,float* dkbt, float* dquat,float* dcurquat, float* domega, float* dlamb_ss,float* dlamb_bt,int* dfix, float dt,int n,int iter,bool example_flag) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	//XPBD�̏����̂��߂ɁC�ɂ�0�ɂ���
	CxSetLambdaZero << <grid, block >> > (dlamb_ss,dlamb_bt, n);
	cudaThreadSynchronize();
	//�L�ѐ���̔���
	for (int i = 0; i < iter; i++) {
		//�S�Ă̐���𓯎��Ɏ��s����ƁC�Փ˂��������邽�߁C��Ƌ����ɕ����Ď��s����
		//�����Ԗڂ�id�����s
		CxStretchingShearConstraint << <grid, block >> > (dpos, dcurpos, dmas, dlen, dkss, dquat, dcurquat, dlamb_ss, dfix, dt, n, 0, i, example_flag);
		//��Ԗڂ�id�����s
		CxStretchingShearConstraint << <grid, block >> > (dpos, dcurpos, dmas, dlen, dkss, dquat, dcurquat, dlamb_ss, dfix, dt, n, 1, i, example_flag);
		
		CxBendTwistConstraint << <grid, block >> > (dmas, dquat, dcurquat, domega, dkbt, dlamb_bt, dlen, dfix, dt, n, 0, i, example_flag);
		CxBendTwistConstraint << <grid, block >> > (dmas, dquat, dcurquat, domega, dkbt, dlamb_bt, dlen, dfix, dt, n, 1, i, example_flag);
	}
}

//�Փː���
//�����ł́C��Ԏ������e�Ղȋ��Ƃ̏Փ˂���������
//dpos:�ʒu
//dvel:���x
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//center:�є��Ƃ̏Փ˂������������̒��S
//rad:�є��Ƃ̏Փ˂������������̔��a
//dt:�^�C���X�e�b�v
//n:���q��
void CuCollisionConstraint(float* dpos, float* dvel, int* dfix, float3 center, float rad, float dt, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxCollisionConstraint << <grid, block >> > (dpos, dvel, dfix, center, rad, dt, n);
	cudaThreadSynchronize();
}

//�C�V��ǉ�
//���Ԑϕ�
//�ʒu�x�[�X�@�ɏ]���C���݂̈ʒu�ƈʒu�C����̈ʒu���瑬�x�����߁C�ʒu���X�V
//dpos:�ʒu(�ʒu�C����)
//dcurpos:�O�X�e�b�v�̈ʒu(�ʒu�C���O)
//dvel:���x
//dt:�^�C���X�e�b�v
//n:���q��
//vel_control:���ȉ��̑��x�̏ꍇ�ɐ؂�̂Ă��s�����ǂ������w��
void CuIntegrate(float* dpos,float* dcurpos,float* dvel,float dt,int n,bool vel_control) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxIntegrate << <grid, block >> > (dpos, dcurpos, dvel, dt, n, vel_control);
	cudaThreadSynchronize();
}

//�C�V��ǉ�
//����d�͂��C���[�W�����O�͌v�Z
//�f�o�b�N�p
void CuCalExternalForces(float* dpos,float*dvel,float* dmass,int* dfix,float3 gravity, float3 wind, float dt, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxCalExternalForces << <grid, block >> > (dpos, dvel, dmass, dfix, gravity, wind, dt, n);
	cudaThreadSynchronize();
}

//�C�V��ǉ�
//�ʒu�x�[�X�@(�g���ʒu�x�[�X�@�łȂ��C�L�сE����f����̂�)
//�f�o�b�N�p
void CuPBDStretchingConstraint(float* dpos, float* dmas, float* dlen, float* dkss, float* dquat, int* dfix, int n, int iter) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	for (int i = 0; i < iter; i++) {
		CxStretchingConstraint<<<grid,block>>>(dpos, dmas, dlen, dkss, dquat, dfix, n, 0);
		cudaThreadSynchronize();
		CxStretchingConstraint<<<grid,block>>>(dpos, dmas, dlen, dkss, dquat, dfix, n, 1);
	}
}

//�f�o�b�N�p�̔z�����ʒu�̏o��
void CuPrint3Dfloat(float* dpos,float* dvel,float* dacc,int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxPrint3Dfloat << <grid, block >> > (dpos, dvel, dacc, n);
	cudaThreadSynchronize();
}

//�ڐ��̍X�V
//kajiya-kay���f���ł̃����_�����O�ɗ��p
//dpos:�ʒu
//dtang:�G�b�W���Ƃ̐ڐ�
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//n:���q��
void CuTangUpdate(float* dpos, float* dtang, int* dfix, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxTangUpdate << <grid, block >> > (dpos, dtang, dfix, n);
	cudaThreadSynchronize();
}

//�p���x�ȂǏ�������f�o�C�X�������ɐݒ�����Ă�����̂̏����l��0�ɐݒ�
//dangvel:�p���x
//dfss:�G�b�W���Ƃɂ������(GlobalForceStep�ŋ��߂�)
//dpbf_lambda:���x����̌v�Z�ߒ��ɕK�v�ȃɂ��������m��
//n:���q��
void CuSetParametersZero(float* dangvel, float* dfss, float* dpbf_lambda, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxSetParametersZero << <grid, block >> > (dangvel,dfss,dpbf_lambda, n);
	cudaThreadSynchronize();
}

//�p�����x�̍X�V
//dangvel:�p���x
//dquat:�p��(�l����)
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//dt:�^�C���X�e�b�v
//n:���q��
void CuAngVelUpdate(float* dangvel, float* dquat,int* dfix,float dt, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxAngVelUpdate << <grid, block >> > (dangvel, dquat,dfix, dt, n);
	cudaThreadSynchronize();
}

//�e�����x�̎��Ԑϕ�
//dangvel:�p���x
//dcurquat:�O�X�e�b�v�̎p��(�ʒu�C���O)
//dquat:���݂̎p��(�ʒu�C����)
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//dt:�^�C���X�e�b�v
//n:���q��
//vel_control:�p���x�����ȉ��Ȃ�؂�̂Ă��s�����ǂ����𔻒�
void CuAngVelIntegrate(float* dangvel,float* dcurquat, float* dquat,int* dfix,float dt, int n,bool vel_control) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxAngVelIntegrate << <grid, block >> > (dangvel, dcurquat, dquat, dfix, dt, n, vel_control);
	cudaThreadSynchronize();
}

//��ƂȂ閧�x�̐ݒ�
//dpos:�ʒu
//dRestDens:���q���Ƃɐݒ肷�����x
//dvol:�̐�
//n:���q��
void CuRestDensSet(float* dpos,float* dRestDens, float* dvol,float* dmas, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxRestDensSet << <grid, block >> > (dpos, dRestDens, dvol, dmas, n);
	cudaThreadSynchronize();
}

//�ꗥ�̊�ƂȂ閧�x�̐ݒ�
//�f�o�b�N�p
void CuRestTotalDens(float* drestdens,float dens, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxRestTotalDens << <grid, block >> > (drestdens, dens, n);
	cudaThreadSynchronize();
}

//�O���[�o���t�H�[�X�X�e�b�v
//�d�͂Ȃǂɂ��G�b�W�ɂ�����͂����߂�
//dfss:�G�b�W���Ƃɂ������
//dmass:����
//last_index:�є����Ƃ̍Ō�̗��q�̃C���f�b�N�X���i�[
//gravity:�d��
//num_elastic:�����ł́C�є����Ƃɕ���v�Z���邽�߁C�є��̐���n��
void CuGlobalForceStep(float* dpos,float* dfss,float* dmass, int* last_index, float3 gravity,float* ddens,float* drestdens,float* dvol, int num_elastic) {
	dim3 block, grid;
	CuCalGridN(num_elastic, block, grid);
	CxGlobalForceStep << <grid, block >> > (dpos,dfss, dmass, last_index, gravity, ddens, drestdens, dvol, num_elastic);
	cudaThreadSynchronize();
}

//���[�J���t�H�[�X�X�e�b�v
//�O���[�o���t�H�[�X�X�e�b�v�ŋ��߂��G�b�W���Ƃ̗͂���C�ό`��h�����߂̊����p�������߂�
//dpos:�ʒu
//dlen:���
//dquat:�p��
//dcurquat:�O�X�e�b�v�̎p��(�V�~�����[�V�����J�n�O��curquat��quat�ƈ�v���邽�߁C�X�V��̒l����)
//dkss:�L�э���
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//n:���q��(�G�b�W���Ƃɕ���v�Z)
void CuLocalForceStep(float* dpos, float* dlen, float* dquat,float* dcurquat, float* dkss, float* dfss, int* dfix, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxLocalForceStep << <grid, block >> > (dpos, dlen, dquat, dcurquat, dkss, dfss, dfix, n);
	cudaThreadSynchronize();
}

//�O���[�o���g���N�X�e�b�v
//�є����ƂɃt�H�[�X�X�e�b�v�Ő������g���N��ł�������_���{�[�x�N�g�������߂�
//dpos:�ʒu
//dquat:�p��
//domega:��_���{�[�x�N�g��
//dlen:���
//dkss:�L�э���
//dkbt:�Ȃ�����
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//last_index:�є����Ƃ̍Ō�̗��q�̃C���f�b�N�X���i�[
//num_elastic:�����ł́C�є����Ƃɕ���v�Z���邽�߁C�є��̐���n��
void CuGlobalTorqueStep(float* dpos, float* dquat, float* domega, float* dlen, float* dkss, float* dkbt, int* dfix, int* dlast_index, int num_elastic) {
	dim3 block, grid;
	CuCalGridN(num_elastic, block, grid);
	//�����狁�߂�
	CxGlobalTorqueStep << <grid, block >> > (dpos, dquat, domega, dlen, dkss, dkbt, dfix, dlast_index, num_elastic);
	//�ォ�狁�߂�
	//CxGlobalTorqueStep_Upstair << <grid, block >> > (dpos, dquat, domega, dlen, dkss, dkbt, dfix, dlast_index, num_elastic);
	cudaThreadSynchronize();
}

//���[�J���g���N�X�e�b�v
//��_���{�[�x�N�g����K�؂Ȍ`�Ő��K������
//dquat:�p��
//domega:��_���{�[�x�N�g��
//deln:���
//dkbt:�Ȃ�����
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//n:���q��(��_���{�[�x�N�g�����Ƃɕ���v�Z)
void CuLocalTorqueStep(float* dquat,float* domega, float* dlen, float* dkbt, int* dfix, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxLocalTorqueStep << <grid, block >> > (dquat, domega, dlen, dkbt, dfix, n);
	cudaThreadSynchronize();
}

//���x����̌v�Z
//pos:�ʒu
//ddens:���݂̖��x
//drestdens:����x
//dpbf_lambda:����ɗp�����
//dvol:�̐�
//n:���q��
void CuPbfConstraint(float* dpos,float* ddens,float* drestdens,float*dpbf_lambda,float*dvol,float* dmas,int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxSphDensity << <grid, block >> > (drestdens, ddens, dvol,dmas, n);//���x�v�Z
	CxPbfLambda << <grid, block >> > (ddens,drestdens, dpbf_lambda, dvol,dmas, n);//����ɗp����ɂ����߂�
	cudaThreadSynchronize();
	CxPbfConstraint << <grid, block >> > (dpos, drestdens, dpbf_lambda, dvol,dmas, n);//���񏈗�
	cudaThreadSynchronize();
}

//PBF�ŉ����ꍇ�̊O�͍��̌v�Z
//dacc:�����x
//datt:���q����(0�ŗ���,1�ŋ��E)
//power:���Ȃǂ̗�
//n:���q��
void CuPbfExternalForces(float* dacc, int* datt, float3 power,bool wind_flag, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxPbfExternalForces << <grid, block >> > (dacc, datt, power, wind_flag, n);
	cudaThreadSynchronize();
}

//���C����
void CuFrictionConstraint(float* dpos, float* dcurpos, float* drestdens, float* dvol, float* ddens, int* dfix, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	//�e���q�Ƃ̖��C�͂����v������C�Î~���C���ǂ����𔻒�
	//CxFrictionConstraint << <grid, block >> > (dpos, dcurpos, drestdens, dvol, ddens, dfix, n);
	//�e���q�ƐÎ~���C���𔻒肵����C���v
	CxFrictionAllParticlesConstraint << <grid, block >> > (dpos, dcurpos, drestdens, dvol, ddens, dfix, n);
	cudaThreadSynchronize();
}

//���C����̌�C�p�����C������
void CuFrictionConstraint_withQuat(float* dpos, float* dcurpos, float* drestdens, float* dvol, float* ddens, float* dquat, float* dlen, int* dfix, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxFrictionConstraint_withQuat << <grid, block >> > (dpos, dcurpos, drestdens, dvol, ddens, dquat, dlen, dfix, n);
	cudaThreadSynchronize();
}

//2���_����p����ݒ�
void CuQuatSet(float* dpos, float* dquat, int* dfix, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxQuatSet << <grid, block >> > (dpos, dquat, dfix, n);
	cudaThreadSynchronize();
}

//�g���N���v�Z������
void CuCalcTorque(float* dpos,float* dmas, float* dquat, float* dfss, float* dlength,float* dkss, int* dfix, float3 gravity, int n) {
	dim3 block, grid;
	CuCalGridN(n, block, grid);
	CxCalcTorque << <grid, block >> > (dpos, dmas, dquat, dfss, dlength, dkss, dfix, gravity, n);
	cudaThreadSynchronize();
}

//--------------------------------------------------------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// GPU�����⏕�֐�
//-----------------------------------------------------------------------------
/*!
 * CUDA�f�o�C�X�̐ݒ� - id�𒼐ڎw��
 * @param[in] id �f�o�C�XID
 */
void CuSetDevice(int id)
{
	int device_count = 0;
	cudaGetDeviceCount(&device_count);
	if(id < 0 || id >= device_count){
		id = 0;
	}
	cudaSetDevice(id);

	cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, id);

	std::cout << " ---- GPU Info ----" << std::endl;
	std::cout << " Device        : " << prop.name << std::endl;
	std::cout << " Global mem    : " << prop.totalGlobalMem << " Byte" << std::endl;
	std::cout << " Constant mem  : " << prop.totalConstMem  << " Byte" << std::endl;
	std::cout << " Thresds/Block : " << prop.maxThreadsPerBlock << std::endl;
	std::cout << std::endl;
	THREAD_NUM = prop.maxThreadsPerBlock; // 1�u���b�N������̃X���b�h���ő�l
}

/*!
 * CUDA�f�o�C�X�̐ݒ�
 *  - �R�}���h���C�������Ɋ�Â�CUDA�f�o�C�X��ݒ�((��)-device 0)
 * @param[in] argc �R�}���h���C�������̐�
 * @param[in] argv �R�}���h���C���������X�g(argv[0]�͎��s�t�@�C����)
 */
void CuInit()
{
	CuSetDevice(0);
}


/*!
 * �f�o�C�X�������̊m��
 * @param[out] dPtr �f�o�C�X�������ւ̃|�C���^
 * @param[in] size �m�ۃT�C�Y(��������̃T�C�Y)
 */
void CuAllocateArray(void** dPtr, size_t size)
{
	CUCHECK(cudaMalloc(dPtr, size));
}

/*!
 * �f�o�C�X�������̉��
 * @param[in] devPtr �f�o�C�X�������ւ̃|�C���^
 */
void CuFreeArray(void* dPtr)
{
	CUCHECK(cudaFree(dPtr));
}

/*!
 * �f�o�C�X�������̈�̏�����
 * @param[in] dPtr �f�o�C�X�������ւ̃|�C���^
 * @param[in] val �����l
 * @param[in] size ����������̈�̃T�C�Y(��������̃T�C�Y)
 */
void CuSetArrayValue(void* dPtr, int val, size_t size)
{
	CUCHECK(cudaMemset(dPtr, val, size));
}

/*!
 * �f�o�C�X�������ԃR�s�[
 * @param[in] dDst �R�s�[��
 * @param[in] dSrc �R�s�[��
 * @param[in] size �R�s�[�T�C�Y(��������̃T�C�Y)
 */
void CuCopyArrayD2D(void* dDst, void* dSrc, int size)
{
	CUCHECK(cudaMemcpy(dDst, dSrc, size, cudaMemcpyDeviceToDevice));
}


/*!
 * VBO���}�b�s���O
 * @param[in] vbo VBO,PBO��
 */
void* CuMapGLBufferObject(cudaGraphicsResource** resource)
{
	void* ptr;
	CUCHECK(cudaGraphicsMapResources(1, resource, 0));
	size_t num_bytes;
	CUCHECK(cudaGraphicsResourceGetMappedPointer((void**)&ptr, &num_bytes, *resource));
	return ptr;
}

/*!
 * VBO���A���}�b�v
 * @param[in] vbo VBO,PBO��
 */
void CuUnmapGLBufferObject(cudaGraphicsResource* resource)
{
	CUCHECK(cudaGraphicsUnmapResources(1, &resource, 0));
}

/*!
 * PBO,VBO�o�b�t�@��CUDA�ɓo�^
 * @param[in] vbo VBO,PBO��
 */
void CuRegisterGLBufferObject(unsigned int vbo, cudaGraphicsResource** resource)
{
	CUCHECK(cudaGraphicsGLRegisterBuffer(resource, vbo, cudaGraphicsMapFlagsNone));
}

/*!
 * PBO,VBO�o�b�t�@��CUDA����폜
 * @param[in] vbo VBO,PBO��
 */
void CuUnregisterGLBufferObject(cudaGraphicsResource* resource)
{
	CUCHECK(cudaGraphicsUnregisterResource(resource));
}

/*!
 * �f�o�C�X����z�X�g�������ւ̃R�s�[
 * @param[in] hDst �R�s�[��z�X�g������(�Œ�size���m�ۂ���Ă��邱��)
 * @param[in] dSrc �R�s�[���f�o�C�X������
 * @param[in] vbo dSrc��VBO�̏ꍇ�CVBO��ID�D�����łȂ��ꍇ��0���w��
 * @param[in] size �R�s�[�T�C�Y(��������̃T�C�Y)
 */
void CuCopyArrayFromDevice(void* hDst, void* dSrc, cudaGraphicsResource** resource, int offset, int size)
{
	if(resource) dSrc = CuMapGLBufferObject(resource);

	CUCHECK(cudaMemcpy(hDst, (char*)dSrc+offset, size, cudaMemcpyDeviceToHost));

	if(resource) CuUnmapGLBufferObject(*resource);
}

/*!
 * �z�X�g����f�o�C�X�������ւ̃R�s�[
 * @param[in] dDst �R�s�[��f�o�C�X������(�Œ�size���m�ۂ���Ă��邱��)
 * @param[in] hSrc �R�s�[���z�X�g������
 * @param[in] offset �R�s�[��I�t�Z�b�g
 * @param[in] size �R�s�[�T�C�Y(��������̃T�C�Y)
 */
void CuCopyArrayToDevice(void* dDst, const void* hSrc, int offset, int size)
{
	CUCHECK(cudaMemcpy((char*)dDst+offset, hSrc, size, cudaMemcpyHostToDevice));
}

/*!
 * �X���b�h����
 */
void CuThreadSync(void)
{
	CUCHECK(cudaThreadSynchronize());
}

/*!
 * �f�o�C�X�v���p�e�B�̕\��
 */
void CuDeviceProp(void)
{
	int n;	//�f�o�C�X��
	CUCHECK(cudaGetDeviceCount(&n));

	for(int i = 0; i < n; ++i){
		cudaDeviceProp dev;

		// �f�o�C�X�v���p�e�B�擾
		CUCHECK(cudaGetDeviceProperties(&dev, i));

		printf("device %d\n", i);
		printf(" device name : %s\n", dev.name);
		printf(" total global memory : %d (MB)\n", (int)dev.totalGlobalMem/1024/1024);
		printf(" shared memory / block : %d (KB)\n", (int)dev.sharedMemPerBlock/1024);
		printf(" register / block : %d\n", dev.regsPerBlock);
		printf(" warp size : %d\n", dev.warpSize);
		printf(" max pitch : %d (B)\n", (int)dev.memPitch);
		printf(" max threads / block : %d\n", dev.maxThreadsPerBlock);
		printf(" max size of each dim. of block : (%d, %d, %d)\n", dev.maxThreadsDim[0], dev.maxThreadsDim[1], dev.maxThreadsDim[2]);
		printf(" max size of each dim. of grid  : (%d, %d, %d)\n", dev.maxGridSize[0], dev.maxGridSize[1], dev.maxGridSize[2]);
		printf(" clock rate : %d (MHz)\n", dev.clockRate/1000);
		printf(" total constant memory : %d (KB)\n", (int)dev.totalConstMem/1024);
		printf(" compute capability : %d.%d\n", dev.major, dev.minor);
		printf(" alignment requirement for texture : %d\n", (int)dev.textureAlignment);
		printf(" device overlap : %s\n", (dev.deviceOverlap ? "ok" : "not"));
		printf(" num. of multiprocessors : %d\n", dev.multiProcessorCount);
		printf(" kernel execution timeout : %s\n", (dev.kernelExecTimeoutEnabled ? "on" : "off"));
		printf(" integrated : %s\n", (dev.integrated ? "on" : "off"));
		printf(" host memory mapping : %s\n", (dev.canMapHostMemory ? "on" : "off"));

		printf(" compute mode : ");
		if(dev.computeMode == cudaComputeModeDefault) printf("default mode (multiple threads can use) \n");
		else if(dev.computeMode == cudaComputeModeExclusive) printf("exclusive mode (only one thread will be able to use)\n");
		else if(dev.computeMode == cudaComputeModeProhibited) printf("prohibited mode (no threads can use)\n");

	}

	//printf("Device with Maximum GFLOPS : %d\n", gpuGetMaxGflopsDeviceId());
}

/*!
 * thrust::exclusive_scan�̌Ăяo��
 * @param[out] dScanData scan��̃f�[�^
 * @param[in] dData ���f�[�^
 * @param[in] num �f�[�^��
 */
void CuScan(float* dScanData, float* dData, int num)
{
	thrust::exclusive_scan(thrust::device_ptr<float>(dData),
		thrust::device_ptr<float>(dData + num),
		thrust::device_ptr<float>(dScanData));
}

/*!
 * �f�o�b�O�p : �X�J���[�l���������z��̒l�̕��ϒl�����߂ĕԂ�
 * @param[out] dcol ���q�F�z��(�f�o�C�X������)
 * @param[in] ddens ���q���x�z��(�f�o�C�X������)
 * @param[in] n ���q��
 * @param[in] userparam ���[�U�[�p�����[�^(�C��)
 */
float CuCalAverage(float* data, int n)
{
	if (n == 0) return 0;
	float avg = 0.0f;

	float* data_scan = 0;
	CuAllocateArray((void**)&data_scan, n * sizeof(float));

	// ���v�l�����߂邽�߂�scan(prefix sum)���v�Z
	CuScan(data_scan, data, n);

	// Exclusive scan (�Ō�̗v�f��0�Ԗڂ���n-2�Ԗڂ܂ł̍��v�ɂȂ��Ă���)�Ȃ̂ŁC
	// Scan�O�z��̍Ō�(n-1�Ԗ�)�̗v�f�ƍ��v���邱�ƂŃ|���S�������v�Z
	float lval, lsval;
	CUCHECK(cudaMemcpy((void*)&lval, (void*)(data + n - 1), sizeof(float), cudaMemcpyDeviceToHost));
	CUCHECK(cudaMemcpy((void*)&lsval, (void*)(data_scan + n - 1), sizeof(float), cudaMemcpyDeviceToHost));
	float total = lval + lsval;
	avg = total / n;

	if (data_scan != 0) CuFreeArray(data_scan);

	return avg;
}
float CuCalAverageV(float* vdata, int n)
{
	if (n == 0) return 0;
	float avg = 0.0f;

	float* data = 0;
	float* data_scan = 0;
	CuAllocateArray((void**)&data, n * sizeof(float));
	CuAllocateArray((void**)&data_scan, n * sizeof(float));

	dim3 block, grid;
	CuCalGridN(n, block, grid);	// �f�[�^��=�X���b�h���Ƃ��ău���b�N/�O���b�h�T�C�Y���v�Z
	CxVectorToScalar<<<grid, block>>>(vdata, data, n);	// �J�[�l�����s
	cudaThreadSynchronize();

	// ���v�l�����߂邽�߂�scan(prefix sum)���v�Z
	CuScan(data_scan, data, n);

	// Exclusive scan (�Ō�̗v�f��0�Ԗڂ���n-2�Ԗڂ܂ł̍��v�ɂȂ��Ă���)�Ȃ̂ŁC
	// Scan�O�z��̍Ō�(n-1�Ԗ�)�̗v�f�ƍ��v���邱�ƂŃ|���S�������v�Z
	float lval, lsval;
	CUCHECK(cudaMemcpy((void*)&lval, (void*)(data + n - 1), sizeof(float), cudaMemcpyDeviceToHost));
	CUCHECK(cudaMemcpy((void*)&lsval, (void*)(data_scan + n - 1), sizeof(float), cudaMemcpyDeviceToHost));
	float total = lval + lsval;
	avg = total / n;

	if (data_scan != 0) CuFreeArray(data_scan);

	return avg;
}

}   // extern "C"
