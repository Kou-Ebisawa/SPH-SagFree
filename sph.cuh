/*! 
  @file sph.cuh
	
  @brief CUDA�֐��̐錾
		 - CUDA���Ăяo��C++�̃R�[�h�͊�{�I�ɂ��̃t�@�C�����C���N���[�h����
 
  @author Makoto Fujisawa
  @date 2023-02
*/

#ifndef _SPH_CUH_
#define _SPH_CUH_


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "cuda_utils.h"

//-----------------------------------------------------------------------------
// CUDA�֐�
//-----------------------------------------------------------------------------
extern "C"
{

//-----------------------------------------------------------------------------
// ���q����
void CuSphDensity(float* drestdens,float* ddens, float* dvol, int n);
void CuSphPressure(float* drestdens,float* dpres, float* ddens, int n);

void CuSphVorticity(float* dvort, float* dvel, float* ddens, float* dvol, int* datt, int n);
void CuSphForces(float* drestdens,float* dacc, float* dvel, float* ddens, float* dpres, float* dvort, float* dvol, int* datt,float3 power, int n);
void CuSphViscosityForces(float* drestdens,float* dacc, float* dvel, float* ddens, float* dvol, int* datt, int n);
void CuSphXSPHViscosity(float* dvel, float* ddens, float* dvol, int* datt, int n);
void CuSphIntegrate(float* dpos, float* dvel, float* dacc, int* datt, int* fix,int n);//�C�V��ύX fix��ǉ�
void CuSphIntegrateV(float* dvel, float* dacc, int* datt, int n);
void CuSphIntegrateP(float* dpos, float* dvel, int* datt, int* dfix, int n);

//�C�V��ǉ�--------------------------------------------------------------------------------
//XPBD�̐���
void CuXPBDConstraint(float* dpos, float* dmas, float* dlen, float* dkss, float* dkbt,float* dquat,float* domega, float* dlamb_ss,float* dlamb_bt,int* dfix, float dt, int n,int iter);
//�Փː���
void CuCollisionConstraint(float* dpos, float* dvel, int* dfix, float3 center, float rad, float dt, int n);
//���Ԑϕ�
void CuIntegrate(float* dpos, float* dcurpos, float* dvel, float dt, int n);
//�O�͌v�Z
void CuCalExternalForces(float* dpos, float* dmass, float* dvel, int* dfix,float3 gravity, float3 wind, float dt, int n);
//PBD�̈ʒu�x�[�X�@
void CuPBDStretchingConstraint(float* dpos, float* dmas, float* dlen, float* dkss, float* dquat, int* dfix, int n, int iter);
//���݂̂���z��̏o��(�f�o�b�N�p)
void CuPrint3Dfloat(float* dpos,float* dvel,float* dacc,int n);
//�ڐ��̍X�V
void CuTangUpdate(float* dpos, float* dtang, int* dfix, int n);
//�d�݂͂̂̌v�Z
void CuOnlyGravity(float* dvel, float* dmass, int* dfix,float dt, int n);
//�������̃p�����[�^��0����
void CuSetParametersZero(float* dangvel,float* dfss,float*dpbf_lambda, int n);
//�p�����x�̍X�V
void CuAngVelUpdate(float* dangvel, float* dquat,int* dfix,float dt, int n);
//�e�����x�̎��Ԑϕ�
void CuAngVelIntegrate(float* dangvel,float* dcurquat, float* dquat,int* dfix,float dt, int n);
//�l�����̐ݒ�
void CuQuatSet(float* dcurquat, float* dquat, int* dfix, int n);
//��ƂȂ�p���̕ύX
void CuRestDensSet(float* dpos, float* dRestDens, float* dvol, int n);
//�ꗥ�̊���x�̐ݒ�
void CuRestTotalDens(float* drestdens,float dens, int n);
//SagFree����-------
//�O���[�o���t�H�[�X�X�e�b�v
void CuGlobalForceStep(float* dfss,float* dmass,int* dlast_index,float3 gravity,int num_elastic);
//���[�J���t�H�[�X�X�e�b�v
void CuLocalForceStep(float* dpos, float* dlen, float* dquat,float* dcurquat, float* dkss, float* dfss, int* dfix, int n);
//�O���[�o���g���N�X�e�b�v
void CuGlobalTorqueStep(float* dpos, float* dquat, float* domega, float* dlen, float* dkss, float* dkbt, int* dfix, int* dlast_index, int num_elastic);
//�O���[�o���g���N�X�e�b�v(Video���Q�l�ɂ�����)
void CuVideoGlobalTorqueStep(float* dpos, float* dquat, float* domega, float* dlen, float* dkss, float* dkbt, int* dfix, int* dlast_index, int num_elastic);
//���[�J���g���N�X�e�b�v
void CuLocalTorqueStep(float* dquat, float* domega, float* dlen, float* dkbt, int* dfix, int n);
//-----------------
//���x����
void CuPbfConstraint(float* dpos, float* ddens, float* drestdens, float* dpbf_lambda, float* dvol, int n);
//PBF�ŉ����ꍇ�̊O�͌v�Z
void CuPbfExternalForces(float* dacc, int* datt, float3 power, int n);
//------------------------------------------------------------------------------------------

// ���q�̐όv�Z
void CuSphCalVolume(float* dvol, int *datt, int n, float v);

// ���q���x�ꐶ��(�\�ʃ��b�V�������p)
void CuSphDensityInGrid(float* dF, float* dvol, int* datt, int n, int3 gnum, float3 gmin, float3 glen);

// ���q�`��F�v�Z
void CuColorScalar(float* dcol, int* datt, float* dval, int n, float3 c1, float3 c2, float2 range);
void CuColorVector(float* dcol, int* datt, float* dval, int n, float3 c1, float3 c2, float2 range);
void CuColorConstant(float* dcol, int* datt, float3 col, int n);

// �ߖT���q�T���p
void CuCalcHash(uint* dhash, uint* dindex, float* dpos, int n);
void CuSort(unsigned int* dhash, uint* dindex, uint n);
void CuReorderDataAndFindCellStart(Cell cell, float* dpos, float* dvel, uint n);

// �p�����[�^��GPU�������ɓ]��
void CuSetParameters(const SceneParameter* hparams);

// �f�o�b�O�p
float CuCalAverage(float* data, int n);		// �X�J���[�l���������z��̒l�̕��ϒl�����߂ĕԂ�
float CuCalAverageV(float* data, int n);	// �x�N�g���l���������z��̃x�N�g���̒����̕��ϒl�����߂ĕԂ�
void CuScan(float* dScanData, float* dData, int num);

//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// GPU�����⏕
void CuInit();
void CuSetDevice(int id);
void CuDeviceProp(void);

void CuAllocateArray(void** devPtr, size_t size);
void CuSetArrayValue(void* devPtr, int val, size_t size);
void CuFreeArray(void* devPtr);

void CuCopyArrayD2D(void* dDst, void* dSrc, int size);
void CuCopyArrayFromDevice(void* host, void* device, cudaGraphicsResource** resource, int offset, int size);
void CuCopyArrayToDevice(void* device, const void* host, int offset, int size);

void CuThreadSync(void);

void CuRegisterGLBufferObject(unsigned int vbo, cudaGraphicsResource** resource);
void CuUnregisterGLBufferObject(cudaGraphicsResource* resource);
void* CuMapGLBufferObject(cudaGraphicsResource** resource);
void CuUnmapGLBufferObject(cudaGraphicsResource* resource);
//-----------------------------------------------------------------------------

} // extern "C"


#endif // #ifdef _SPH_CUH_