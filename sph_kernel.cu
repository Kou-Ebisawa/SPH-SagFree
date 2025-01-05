/*!
  @file sph_kernel.cu

  @brief CUDA : SPH�@
         - CUDA�J�[�l���y�уf�o�C�X�֐����L�q
         - �z�X�g�֐��������ꂽ*.cu�t�@�C������̂݃C���N���[�h(������cpp����C���N���[�h���Ȃ��悤��)

  @author Makoto Fujisawa
  @date 2023-02
 */


//-----------------------------------------------------------------------------
// �C���N���[�h�t�@�C��
//-----------------------------------------------------------------------------
#include "helper_math.h"
#include <math_constants.h>

#include "cuda_utils.h"
#include "cuda_utils.cu"


//-----------------------------------------------------------------------------
// �萔(�f�o�C�X������)
//-----------------------------------------------------------------------------
__device__ __constant__ SceneParameter params;	// �e��p�����[�^


//-----------------------------------------------------------------------------
// device�֐� - �f�o�C�X(GPU)�Ŏ��s�E�f�o�C�X�֐�����̂݌Ăяo����
//-----------------------------------------------------------------------------
/*!
* �O���b�h�ʒu�v�Z
* @param[in] p ���W
* @return �O���b�h���W
*/
__device__
inline int3 calcGridPos(float3 p)
{
    int3 grid;
    grid.x = floor((p.x-params.cell.WorldOrigin.x)/params.cell.CellWidth.x);
    grid.y = floor((p.y-params.cell.WorldOrigin.y)/params.cell.CellWidth.y);
    grid.z = floor((p.z-params.cell.WorldOrigin.z)/params.cell.CellWidth.z);

    grid.x = min(max(grid.x, 0), params.cell.GridSize.x-1);
    grid.y = min(max(grid.y, 0), params.cell.GridSize.y-1);
    grid.z = min(max(grid.z, 0), params.cell.GridSize.z-1);

    return grid;
}

/*!
* �O���b�h���W����1�����z�񒆂ł̈ʒu���v�Z
* @param[in] gridPos �O���b�h���W
* @return �A�h���X
*/
__device__
inline uint calcGridHash(int3 gridPos)
{
    return __umul24(__umul24(gridPos.z, params.cell.GridSize.y), params.cell.GridSize.x)+__umul24(gridPos.y, params.cell.GridSize.x)+gridPos.x;
}

/*!
* ���q�̏Փˏ����m�ِ�
* @param[inout] p,v ���q�ʒu,���x
* @param[in] dt �^�C���X�e�b�v��
*/
__device__
void collision(float3 &p, float3 &v, float dt)
{
    float d;
    float3 nrm, cp;
    float res = params.res;

    // �{�b�N�X�`��̃I�u�W�F�N�g�Ƃ̏Փ�
#if MAX_BOX_NUM
    for(int i = 0; i < params.num_box; ++i){
        if(params.box[i].flg == 0) continue;
        collisionPointBox(p, params.box[i].cen, params.box[i].ext+make_float3(params.particle_radius), params.box[i].rot, params.box[i].inv_rot, cp, d, nrm);
        if(d < 0.0){
            res = (res > 0) ? (res*fabs(d)/(dt*length(v))) : 0.0f;
            v -= (1+res)*nrm*dot(nrm, v);
            p = cp;
        }
    }
#endif

    // ���`��̃I�u�W�F�N�g�Ƃ̏Փ�
#if MAX_SPHERE_NUM
    for(int i = 0; i < params.num_sphere; ++i){
        if(params.sphere[i].flg == 0) continue;
        collisionPointSphere(p, params.sphere[i].cen, params.sphere[i].rad+params.particle_radius, cp, d, nrm);
        if(d < 0.0){
            res = (res > 0) ? (res*fabs(d)/(dt*length(v))) : 0.0f;
            v -= (1+res)*nrm*dot(nrm, v);
            p = cp;
        }
    }
#endif

    // �V�~�����[�V������Ԃ̋��E(AABB)�Ƃ̏Փ�
    float3 l0 = params.boundary_min;
    float3 l1 = params.boundary_max;
    if(distPointAABB(p, 0.5*(l1+l0), 0.5*(l1-l0), cp, d, nrm)){
        res = (res > 0) ? (res*fabs(d)/(dt*length(v))) : 0.0f;
        v -= (1+res)*nrm*dot(nrm, v);
        p = cp;
    }
}


//-----------------------------------------------------------------------------
// global�֐� - �f�o�C�X(GPU)�Ŏ��s�E�z�X�g�֐�����̂݌Ăяo����
//-----------------------------------------------------------------------------
/*!
 * SPH�@�ɂ�闱�q���x�̌v�Z(Poly6�J�[�l��)
 * @param[out] ddens ���q���x
 * @param[in] dvol ���q�̐�
 * @param[in] n ���q��
 */
__global__ 
void CxSphDensity(float*drestdens,float* ddens, float* dvol, float* dmas,int n)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    float3 pos0 = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    float m = params.mass;
    float a = params.aw;
    float rest_dens = params.rest_dens;
    //�C���f�b�N�X�̌v�Z
    uint sid = params.cell.dSortedIndex[id];
    //�C�V��SPH�ǉ�
    rest_dens = drestdens[sid];

    // ���q�𒆐S�Ƃ��Ĕ��ah���Ɋ܂܂��O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0-make_float3(h));
    grid_pos1 = calcGridPos(pos0+make_float3(h));

    // ���͂̃O���b�h�Z�����܂߂ċߖT�T�����Ė��x�v�Z
    float dens = 0.0f;
    for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
        for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
            for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // �O���b�h�n�b�V���l

                // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                uint startIndex = params.cell.dCellStart[ghash];
                if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
                    // �Z�����̃p�[�e�B�N���Ŕ���
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for(uint j = startIndex; j < endIndex; ++j){
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 rij = pos0-pos1;
                        float r = length(rij);
                        if(r <= h){
                            // Poly6�J�[�l���Ŗ��x���v�Z (rho = �� m Wij)
                            float q = h*h-r*r;

                            float m = params.mass;

                            dens += m*a*q*q*q;
                        }
                    }
                }
            }
        }
    }
    ddens[sid] = dens;
}

/*!
 * ���q�ɓ����͂̌v�Z
 *  - �d�́C���͂Ȃ�
 * @param[out] dpres ���q����
 * @param[in] ddens ���q���x
 * @param[in] n ���q��
 */
__global__ 
void CxSphPressure(float* drestdens,float* dpres, float* ddens, int n)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    float p = 0.0f;

    // ���z�C�̂̏�ԕ������Ɋ�Â��v�Z[Muller2003]
    //p = params.gas_k * (ddens[id] - params.rest_dens);

    // Tait�������Ɋ�Â��v�Z:WCSPH[Becker2007]
    float rdens = ddens[id]/params.rest_dens;
    //�C�V��SPH�ǉ�
    rdens = ddens[id] / drestdens[id];

    //p = params.B * (powf(rdens, params.gamma)-1.0f);
    // gamma=7�Œ�̏ꍇ
    p = params.B * (rdens*rdens*rdens*rdens*rdens*rdens*rdens-1.0f);

    // �����̏ꍇ��0�ɂ���
    p = clamp(p, 0.0, 1.0e6);

    // ���q���͂̍X�V(�O���[�o��������������������)
    dpres[id] = p;
}



/*!
* ���q�Q�x�̌v�Z
*  - vorticity confinement�ɂ�闐���\���̂��߂̉Q�x�v�Z
* @param[out] dvort �e���q�̉Q�x�x�N�g��(�f�o�C�X������)
* @param[in] dvel ���q���x
* @param[in] ddens ���q���x
* @param[in] dvol ���q�̐�
* @param[in] datt ���q����(0�ŗ���,1�ŋ��E)
* @param[in] n ���q��
*/
__global__ 
void CxSphVorticity(float* dvort, float* dvel, float* ddens, float* dvol, int* datt, int n)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    uint sid = params.cell.dSortedIndex[id];
    if(datt[sid] != 0){  // ���E���q�̏ꍇ�͉Q�x=0
        float3 v0 = make_float3(0.0f);
        dvort[DIM*sid+0] = v0.x;  dvort[DIM*sid+1] = v0.y; dvort[DIM*sid+2] = v0.z;
        return;
    }

    // ���qi�̕ϐ��l
    float3 pos0 = params.cell.dSortedPos[id];
    float3 vel0 = make_float3(dvel[DIM*sid], dvel[DIM*sid+1], dvel[DIM*sid+2]);

    float h = params.h;
    float m = params.mass;
    float a = params.ag;

    // �p�[�e�B�N�����͂̃O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0-make_float3(h));
    grid_pos1 = calcGridPos(pos0+make_float3(h));

    // ���͂̃O���b�h�Z�����܂߂ċߖT�T�����Ĉ��͂ɂ���
    float3 f = make_float3(0.0f);
    for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
        for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
            for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);

                // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                uint startIndex = params.cell.dCellStart[ghash];
                if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
                    // �Z�����̃p�[�e�B�N���Ŕ���
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for(uint j = startIndex; j < endIndex; ++j){
                        uint sj = params.cell.dSortedIndex[j];
                        if(sj == sid) continue;

                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 vel1 = make_float3(dvel[DIM*sj], dvel[DIM*sj+1], dvel[DIM*sj+2]);

                        float3 rij = pos0-pos1;
                        float3 vji = vel1-vel0;
                        float dens1 = ddens[sj];

                        float r = length(rij);
                        if(r <= h && r > 0.0001){
                            float q = h-r;
                            float3 gw = a*q*q*rij/r;
                            f += m/dens1*cross(vji, gw);
                        }
                    }
                }
            }
        }
    }

    // ���q�̉Q�x�x�N�g���̍X�V(�O���[�o��������������������)
    dvort[DIM*sid+0] = f.x;
    dvort[DIM*sid+1] = f.y;
    dvort[DIM*sid+2] = f.z;
}


/*!
 * ���q�ɓ����͂̌v�Z
 *  - ���x�����ɂ���悤�Ȉ���
 *  - �d��
 *  - vorticity confinement
 * @param[out] dacc ���q�ɓ�����(�����x)
 * @param[in] dvel ���q���x�z��
 * @param[in] ddens ���q���x
 * @param[in] dpres ���q����
 * @param[in] dvort �e���q�̉Q�x�x�N�g��
 * @param[in] dvol  ���q�̐�
 * @param[in] datt ���q����(0�ŗ���,1�ŋ��E)
 * @param[in] n ���q��
 */
//�C�V��power�ǉ�
//�C�V��dfss�ǉ�
__global__ 
void CxSphForces(float* drestdens,float* dacc, float* dvel, float* ddens, float* dpres, float* dvort, float* dvol,float* dmas, int* datt,float3 power,float* dfss, int n)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    uint sid = params.cell.dSortedIndex[id];
    if(datt[sid] != 0){  // ���E���q�̏ꍇ�͗��q�ɂ������=0
        float3 v0 = make_float3(0.0f);
        dacc[DIM*sid+0] = v0.x;  dacc[DIM*sid+1] = v0.y; dacc[DIM*sid+2] = v0.z;
        return;
    }

    // ���qi�̕ϐ��l
    float3 pos0 = params.cell.dSortedPos[id];
    float3 omega0 = make_float3(dvort[DIM*sid], dvort[DIM*sid+1], dvort[DIM*sid+2]);
    float dens0 = ddens[sid];
    float pres0 = dpres[sid];
    //int3 grid = calcGridPos(pos0);
    float prsi = pres0/(dens0*dens0);

    float h = params.h;
    //float m = params.mass;
    float a = params.ag;

    // �p�[�e�B�N�����͂̃O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0-make_float3(h));
    grid_pos1 = calcGridPos(pos0+make_float3(h));

    // ���͂̃O���b�h�Z�����܂߂ċߖT�T�����Ĉ��͂ɂ���
    float3 f = make_float3(0.0f);
    float3 eta = make_float3(0.0f);
    for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
        for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
            for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);

                // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                uint startIndex = params.cell.dCellStart[ghash];
                if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
                    // �Z�����̃p�[�e�B�N���Ŕ���
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for(uint j = startIndex; j < endIndex; ++j){
                        uint sj = params.cell.dSortedIndex[j];
                        if(sj == sid) continue;

                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 omega1 = make_float3(dvort[DIM*sj], dvort[DIM*sj+1], dvort[DIM*sj+2]);

                        float3 rij = pos0-pos1;

                        float dens1 = ddens[sj];
                        float pres1 = dpres[sj];
                        float prsj = pres1/(dens1*dens1);

                        float r = length(rij);
                        if(r <= h && r > 0.0001){
                            float q = h-r;
                            float3 gw = a*q*q*rij/r;

                            float m = params.mass;

                            f += -m*(prsi+prsj)*gw; // ���͍��̌v�Z
                            eta += (m/dens1)*length(omega1)*gw;
                        }
                    }
                }
            }
        }
    }
    
    // �d��
    f += params.gravity+power;
    //f += power;

    // Vorticity Confinement
    if(length(eta) > 1e-3){
        f += params.vorticity*dens0*cross(normalize(eta), omega0);
    }

    // ���q�ɂ�����O��(�����x)�̍X�V(�O���[�o��������������������)
    dacc[DIM*sid+0] = f.x; 
    dacc[DIM*sid+1] = f.y;
    dacc[DIM*sid+2] = f.z;
}

/*!
* ���q�ɓ����S���͂̌v�Z
*  - XSPH�łȂ��S�������(�����x)�Ƃ��Čv�Z������@[Becker2007]
* @param[inout] dacc ���q�ɓ�����(�����x)
* @param[in] dvel ���q���x
* @param[in] ddens ���q���x
* @param[in] dvol ���q�̐�
* @param[in] datt ���q����(0�ŗ��̗��q�C����ȊO�ŋ��E���q)
* @param[in] n ���q��
*/
__global__ 
void CxSphViscosity(float* drestdens,float* dacc, float* dvel, float* ddens, float* dvol,float* dmas, int* datt, int n)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    uint sid = params.cell.dSortedIndex[id];
    if(datt[sid] != 0) return;

    // ���qi�̕ϐ��l
    float3 pos0 = params.cell.dSortedPos[id];
    float3 vel0 = make_float3(dvel[DIM*sid], dvel[DIM*sid+1], dvel[DIM*sid+2]);
    float dens0 = ddens[sid];

    float h = params.h;
    float a = params.ag;
    float alpha = params.viscosity;  // �S���萔�D�_������[0.08, 0.5]�ƍs���Ă��邪���ꂾ�Ƒ傫������...
    float cs = 88.5;
    float eps = 0.001*h*h;
    float rest_dens = params.rest_dens;
    //�C�V��SPH�ǉ�
    rest_dens = drestdens[sid];

    // �p�[�e�B�N�����͂̃O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0-make_float3(h));
    grid_pos1 = calcGridPos(pos0+make_float3(h));

    // ���͂̃O���b�h�Z�����܂߂ċߖT�T�����Ĉ��͂ɂ���
    float3 f = make_float3(0.0f);
    for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
        for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
            for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);

                // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                uint startIndex = params.cell.dCellStart[ghash];
                if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
                    // �Z�����̃p�[�e�B�N���Ŕ���
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for(uint j = startIndex; j < endIndex; ++j){
                        uint sj = params.cell.dSortedIndex[j];
                        if(sj == sid) continue;

                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 vel1 = make_float3(dvel[DIM*sj], dvel[DIM*sj+1], dvel[DIM*sj+2]);
                        float3 rij = pos0-pos1;

                        float dens1 = ddens[sj];
                        float m = dmas[j];
                        float vx = dot(vel0-vel1, rij);

                        float r = length(rij);
                        if(r <= h && r > 0.0001 && vx < 0.0f){
                            float nu = (2.0f*alpha*h*cs)/(dens0+dens1);
                            float visc = -nu*(vx/(r*r+eps));
                            float q = h-r;
                            f += -m*visc*a*q*q*rij/r;
                        }
                    }
                }
            }
        }
    }

    // ���q�ɂ�����O��(�����x)�̍X�V(�O���[�o��������������������)
    dacc[DIM*sid+0] += f.x; 
    dacc[DIM*sid+1] += f.y;
    dacc[DIM*sid+2] += f.z;
}

/*!
* ���q�ɓ����͂̌v�Z
*  - XSPH Artificial Viscosity�ɂ�鑬�x�X�V[Schechter2012]
* @param[inout] dvel ���q���x
* @param[in] ddens ���q���x
* @param[in] dvol ���q�̐�
* @param[in] datt ���q����(0�ŗ��̗��q�C����ȊO�ŋ��E���q)
* @param[in] n ���q��
*/
__global__ 
void CxSphXSPHViscosity(float* dvel, float* ddens, float* dvol,float* dmas, int* datt, int n)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    uint sid = params.cell.dSortedIndex[id];
    if(datt[sid] != 0) return;

    // ���qi�̕ϐ��l
    float3 pos0 = params.cell.dSortedPos[id];
    float3 vel0 = make_float3(dvel[DIM*sid], dvel[DIM*sid+1], dvel[DIM*sid+2]);

    float h = params.h;
    float a = params.aw;
    float eps = params.viscosity;
    float rest_dens = params.rest_dens;

    // �p�[�e�B�N�����͂̃O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0-make_float3(h));
    grid_pos1 = calcGridPos(pos0+make_float3(h));

    // ���͂̃O���b�h�Z�����܂߂ċߖT�T�����Ĉ��͂ɂ���
    float3 dv = make_float3(0.0f);
    for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
        for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
            for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);

                // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                uint startIndex = params.cell.dCellStart[ghash];
                if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
                    // �Z�����̃p�[�e�B�N���Ŕ���
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for(uint j = startIndex; j < endIndex; ++j){
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 vel1 = make_float3(dvel[DIM*sj], dvel[DIM*sj+1], dvel[DIM*sj+2]);
                        float3 rij = pos0-pos1;

                        float dens1 = ddens[sj];
                        float m = dmas[sj];

                        float r = length(rij);
                        if(r <= h && r > 0.0001){
                            float q = h*h-r*r;
                            //float m = rest_dens*dvol[sj];
                            dv += (m/dens1)*(vel1-vel0)*a*q*q*q;
                        }
                    }
                }
            }
        }
    }

    dv *= eps;

    // ���q���x�̍X�V(�O���[�o��������������������)
    dvel[DIM*sid+0] += dv.x; 
    dvel[DIM*sid+1] += dv.y;
    dvel[DIM*sid+2] += dv.z;
}


/*!
 * ���q��O�i�I�C���[�@�ňړ�������
 *  - �ʒu�̑��x�ɂ��ϕ�
 *  - ���E�������܂�
 * @param[inout] dpos ���q�ʒu
 * @param[inout] dvel ���q���x
 * @param[in] dacc ���q�ɓ�����(�����x)
 * @param[in] datt ���q����(0�ŗ��̗��q�C����ȊO�ŋ��E���q)
 * fix:�C�V��ǉ� 1�Ȃ�ΌŒ�_
* @param[in] n ���q��
 */
__global__ 
void CxSphIntegrate(float* dpos, float* dvel, float* dacc, int* datt,int* dfix, int n)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)
    if(datt[id] != 0){  // ���E���q�̏ꍇ�͑��x��0�ɂ��Ĉʒu�͕ς��Ȃ�
        float3 v0 = make_float3(0.0f);
        dvel[DIM*id+0] = v0.x;  dvel[DIM*id+1] = v0.y; dvel[DIM*id+2] = v0.z;
        return;
    }
    //�C�V��ǉ�
    //�Œ�_�Ȃ�΁A���x��0�ɂ��ď������X�L�b�v
    if (dfix[id] == 1||id==0 || dfix[id - 1] == 1) {// 
        float3 v0 = make_float3(0.0f);
        dvel[DIM * id + 0] = v0.x;  dvel[DIM * id + 1] = v0.y; dvel[DIM * id + 2] = v0.z;
        return;
    }
    // ���q�ʒu,���x,��
    float3 p = make_float3(dpos[DIM*id+0], dpos[DIM*id+1], dpos[DIM*id+2]);
    float3 v = make_float3(dvel[DIM*id+0], dvel[DIM*id+1], dvel[DIM*id+2]);
    float3 a = make_float3(dacc[DIM*id+0], dacc[DIM*id+1], dacc[DIM*id+2]);
    float dt = params.dt;

    // �X�V�ʒu�C���x�̍X�V
    v += dt*a;
    p += dt*v;

    // ���͋��E�Ƃ̏Փˏ���
    collision(p, v, dt);

    // ���q�ʒu�E���x�̍X�V(�O���[�o��������������������)
    dpos[DIM*id+0] = p.x;  dpos[DIM*id+1] = p.y; dpos[DIM*id+2] = p.z;
    dvel[DIM*id+0] = v.x;  dvel[DIM*id+1] = v.y; dvel[DIM*id+2] = v.z;
}

/*!
* ���q��O�i�I�C���[�@�ňړ�������
*  - ���x�̎��Ԑϕ��̂�, XSPH�̂݁C���E�����Ȃ�
* @param[inout] dvel ���q���x
* @param[in] dacc ���q�ɓ�����(�����x)
* @param[in] datt ���q����(0�ŗ��̗��q�C����ȊO�ŋ��E���q)
* @param[in] n ���q��
*/
__global__ 
void CxSphIntegrateVelocity(float* dvel, float* dacc, int* datt, int n)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)
    if(datt[id] != 0){  // ���E���q�̏ꍇ�͑��x��0�ɂ��Ĉʒu�͕ς��Ȃ�
        float3 v0 = make_float3(0.0f);
        dvel[DIM*id+0] = v0.x;  dvel[DIM*id+1] = v0.y; dvel[DIM*id+2] = v0.z;
        return;
    }

    // ���q�ʒu,���x,��
    float3 v = make_float3(dvel[DIM*id+0], dvel[DIM*id+1], dvel[DIM*id+2]);
    float3 a = make_float3(dacc[DIM*id+0], dacc[DIM*id+1], dacc[DIM*id+2]);
    float dt = params.dt;

    // �X�V�ʒu�C���x�̍X�V
    v += dt*a;

    // ���q�ʒu�E���x�̍X�V(�O���[�o��������������������)
    dvel[DIM*id+0] = v.x;  dvel[DIM*id+1] = v.y; dvel[DIM*id+2] = v.z;
}
/*!
* ���q��O�i�I�C���[�@�ňړ�������
*  - �ʒu�̑��x�ɂ��ϕ��̂݁CXSPH�p�C���E�������܂�
* @param[inout] dpos ���q�ʒu
* @param[inout] dvel ���q���x
* @param[in] datt ���q����(0�ŗ��̗��q�C����ȊO�ŋ��E���q)
* @param[in] n ���q��
*/
//�C�V�� dfix�ǉ�
__global__ 
void CxSphIntegratePosition(float* dpos, float* dvel, int* datt,int* dfix, int n)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n || datt[id] != 0) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    //�C�V��ǉ�
    //�Œ�_�Ȃ�΁A���x��0�ɂ��ď������X�L�b�v
    if (dfix[id] == 1 || id == 0 || dfix[id - 1] == 1) { 
        float3 v0 = make_float3(0.0f);
        dvel[DIM * id + 0] = v0.x;  dvel[DIM * id + 1] = v0.y; dvel[DIM * id + 2] = v0.z;
        return;
    }

    // ���q�ʒu,���x,��
    float3 p = make_float3(dpos[DIM*id+0], dpos[DIM*id+1], dpos[DIM*id+2]);
    float3 v = make_float3(dvel[DIM*id+0], dvel[DIM*id+1], dvel[DIM*id+2]);
    float dt = params.dt;

    // �X�V�ʒu�C���x�̍X�V
    p += dt*v;

    // ���͋��E�Ƃ̏Փˏ���
    collision(p, v, dt);

    // ���q�ʒu�E���x�̍X�V(�O���[�o��������������������)
    dpos[DIM*id+0] = p.x;  dpos[DIM*id+1] = p.y; dpos[DIM*id+2] = p.z;
    dvel[DIM*id+0] = v.x;  dvel[DIM*id+1] = v.y; dvel[DIM*id+2] = v.z;
}


/*!
* ���E���q�����̂��߂̗��q�̐όv�Z
*  - ���E���q�� "Versatile Rigid-Fluid Coupling for Incompressible SPH", 2.2 ��(3)�̏��V_bi �Ōv�Z
*  - ���̗��q�� V=mass/rest_dens
* @param[out] dvol ���q�̐�
* @param[in] datt ���q����(0�ŗ��̗��q�C����ȊO�ŋ��E���q)
* @param[in] n �������闱�q��(offset����̑��ΓI�Ȉʒu)
* @param[in] v ���̗��q�̏ꍇ�̗��q�̐ϒl
*/
__global__ 
void CxSphCalVolume(float* dvol, int *datt, int n, float v)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    uint sid = params.cell.dSortedIndex[id];
    int att = datt[sid];
    if(att == 0){   // ���̗��q�̏ꍇ
        dvol[sid] = v;
        return;
    }

    float3 pos0 = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    float m = params.mass;
    float a = params.aw;

    // ���q�𒆐S�Ƃ��Ĕ��ah���Ɋ܂܂��O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0-make_float3(h));
    grid_pos1 = calcGridPos(pos0+make_float3(h));

    // ���͂̃O���b�h�Z�����܂߂ċߖT�T�����đ̐όv�Z
    float mw = 0.0f;    // ��mW
    for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
        for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
            for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // �O���b�h�n�b�V���l
                uint startIndex = params.cell.dCellStart[ghash];                // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
                    // �Z�����̃p�[�e�B�N���Ŕ���
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for(uint j = startIndex; j < endIndex; ++j){
                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 rij = pos0-pos1;
                        float r = length(rij);
                        if(r <= h){
                            float q = h*h-r*r;
                            mw += m*a*q*q*q;
                        }
                    }
                }
            }
        }
    }

    // �v�Z�����̐ς��O���[�o���������ɏ�������
    dvol[sid] = m/mw;
}


/*!
* �O���b�h��ł̖��x���v�Z(�\�ʃ��b�V�������p)
* @param[out] dF ���x�l���i�[����O���b�h�Z���z��
* @param[in] dvol ���q�̐�
* @param[in] datt ���q����(0�ŗ���,1�ŋ��E)
* @param[in] n ���q��
* @param[in] gnum �O���b�h��
* @param[in] gmin �O���b�h�ŏ����W
* @param[in] glen �O���b�h��
*/
__global__
void CxSphDensityAtCell(float* dF, float* dvol, int* datt, int n, 
                        int3 gnum, float3 gmin, float3 glen)
{
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    int3 gridPos = calcGridPos(id, gnum);

    if(gridPos.x < gnum.x && gridPos.y < gnum.y && gridPos.z < gnum.z){
        float3 pos0;    // �O���b�h�Z�����S���W
        pos0.x = gmin.x+(gridPos.x)*glen.x;
        pos0.y = gmin.y+(gridPos.y)*glen.y;
        pos0.z = gmin.z+(gridPos.z)*glen.z;

        float h = params.effective_radius;
        float m = params.mass;
        float a = params.aw;

        int3 grid_pos0, grid_pos1;
        grid_pos0 = calcGridPos(pos0-make_float3(h));
        grid_pos1 = calcGridPos(pos0+make_float3(h));

        float dens = 0.0f;
        for(int z = grid_pos0.z; z <= grid_pos1.z; ++z){
            for(int y = grid_pos0.y; y <= grid_pos1.y; ++y){
                for(int x = grid_pos0.x; x <= grid_pos1.x; ++x){
                    int3 ngrid = make_int3(x, y, z);
                    uint ghash = calcGridHash(ngrid);   // �O���b�h�n�b�V���l

                    // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                    uint startIndex = params.cell.dCellStart[ghash];
                    if(startIndex != 0xffffffff){	// �Z������łȂ����̃`�F�b�N
                        // �Z�����̃p�[�e�B�N���Ŕ���
                        uint endIndex = params.cell.dCellEnd[ghash];
                        for(uint j = startIndex; j < endIndex; ++j){
                            uint sj = params.cell.dSortedIndex[j];
                            if(datt[sj] != 0) continue; // ���E���q�͕\�ʃ��b�V�������ɂ͎g��Ȃ�
                            float3 pos1 = params.cell.dSortedPos[j];
                            float3 rij = pos0-pos1;
                            float r = length(rij);
                            if(r <= h){
                                // Poly6�J�[�l���Ŗ��x���v�Z (rho = �� m Wij)
                                float q = h*h-r*r;
                                dens += m*a*q*q*q;
                            }
                        }
                    }
                }
            }
        }

        dF[gridPos.x+gridPos.y*gnum.x+gridPos.z*gnum.x*gnum.y] = dens;
    }

}


//-----------------------------------------------------------------------------
// ���q�̕`��F�̌v�Z
//-----------------------------------------------------------------------------
/*!
* ���q�̕`��F�𗱎q�̎������ʂ���v�Z
* @param[out] dcol ���q�F�z��(�f�o�C�X������)
* @param[in] dval ���q�����ʔz��(�f�o�C�X������)
* @param[in] n ���q��
*/
__global__ 
void CxColorScalar(float* dcol, int* datt, float* dval, int n, float3 c1, float3 c2, float2 range)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    // �����ʂƑ���
    float val = dval[id];
    int att = datt[id];

    // ���q�F�̌v�Z
    float l = range.y-range.x;
    float t = clamp((val-range.x)/l, 0.0f, 1.0f);
    float3 col = lerp(c1, c2, t);

    // ���q�F�̃O���[�o���������ւ̊i�[
    dcol[4*id+0] = col.x;  dcol[4*id+1] = col.y; dcol[4*id+2] = col.z;
    dcol[4*id+3] = (1.0f-(float)(att));
}
__global__ 
void CxColorVector(float* dcol, int* datt, float* dval, int n, float3 c1, float3 c2, float2 range)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    // �����ʂƑ���
    float val = length(make_float3(dval[DIM*id], dval[DIM*id+1], dval[DIM*id+2]));
    int att = datt[id];

    // ���q�F�̌v�Z
    float l = range.y-range.x;
    float t = clamp((val-range.x)/l, 0.0f, 1.0f);
    float3 col = lerp(c1, c2, t);

    // ���q�F�̃O���[�o���������ւ̊i�[
    dcol[4*id+0] = col.x;  dcol[4*id+1] = col.y; dcol[4*id+2] = col.z;
    dcol[4*id+3] = (1.0f-(float)(att));
}

/*!
* ���q�̕`��F�ݒ� : ���ׂē����F
* @param[out] dcol ���q�F�z��(�f�o�C�X������)
* @param[in] ddens ���q���x�z��(�f�o�C�X������)
* @param[in] n ���q��
*/
__global__ 
void CxColorConstant(float* dcol, int* datt, float3 col, int n)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)
    
    // ���q����
    int att = datt[id];

    // ���q�F�̃O���[�o���������ւ̊i�[
    dcol[4*id+0] = col.x;  dcol[4*id+1] = col.y; dcol[4*id+2] = col.z;
    dcol[4*id+3] = (1.0f-(float)(att));
}

/*!
* �f�o�b�O�p : �x�N�g���z�񂩂�x�N�g���̑傫��(�X�J���[�l)�̔z����v�Z
* @param[in] dvdata �x�N�g���l�z��(�f�o�C�X������)
* @param[out] dsdata �X�J���[�l�z��(�f�o�C�X������)
* @param[in] n ���q��
*/
__global__ 
void CxVectorToScalar(float* dvdata, float* dsdata, int n)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    // �x�N�g���̑傫�������߂�
    float val = length(make_float3(dvdata[DIM*id], dvdata[DIM*id+1], dvdata[DIM*id+2]));

    // ���ʂ̃O���[�o���������ւ̊i�[
    dsdata[id] = val;
}


//-----------------------------------------------------------------------------
// �n�b�V��
//-----------------------------------------------------------------------------
/*!
 * �e���q�̃O���b�h�n�b�V���l�v�Z
 * @param[out] dhash �e���q�̃O���b�h�n�b�V���l���i�[�����z��
 * @param[out] dsortedidx �e���q�̃C���f�b�N�X���i�[�����z��(�ォ��n�b�V���l�Ń\�[�g����� -> �����_�ł͂܂��\�[�h�ς݂ł͂Ȃ�)
 * @param[in] dpos ���q�ʒu���i�[�����z��
 * @param[in] n ���q��
 */
__global__
void CxCalcHash(uint* dhash, uint* dsortedidx, float* dpos, uint n)
{
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    // ���q�ʒu
    float3 p = make_float3(dpos[DIM*id+0], dpos[DIM*id+1], dpos[DIM*id+2]);

    // ���q�ʒu����܂܂��O���b�h�Z���̃n�b�V���l���v�Z
    int3 grid = calcGridPos(p);
    uint hash = calcGridHash(grid);

    dhash[id] = hash;
    dsortedidx[id] = id;
}

/*!
 * �p�[�e�B�N���f�[�^���\�[�g���āC�n�b�V�����̊e�Z���̍ŏ��̃A�h���X������
 * @param[in] cell �p�[�e�B�N���O���b�h�f�[�^
 * @param[in] dpos ���q�ʒu�z��
 * @param[in] dvel ���q���x�z��
 */
__global__
void CxReorderDataAndFindCellStartD(Cell cell, float* dpos, float* dvel, uint n)
{
    // �V�F�A�[�h������
    extern __shared__ uint sharedHash[];	// �T�C�Y : blockSize+1

    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    uint hash;
    if(id < n){ // �V�F�A�[�h�������g�p�̂��߂ɃX���b�h�������s���̂�(id >= n)�̎���return�͂��Ȃ�
        hash = cell.dGridParticleHash[id];	// �n�b�V���l
        sharedHash[threadIdx.x+1] = hash;	// �n�b�V���l���V�F�A�[�h�������Ɋi�[

        if(id > 0 && threadIdx.x == 0){
            // �e�V�F�A�[�h�������̍ŏ��ׂ͗̃O���b�h�̃p�[�e�B�N���̃n�b�V���l���i�[
            sharedHash[0] = cell.dGridParticleHash[id-1];
        }
    }

    __syncthreads();    // �X���b�h����(���̃X���b�h���V�F�A�[�h�������i�[�������I����܂ő҂�)

    if(id < n){
        // �C���f�b�N�X0�ł���C�������́C��O�̃p�[�e�B�N���̃O���b�h�n�b�V���l���قȂ�ꍇ�C
        // �p�[�e�B�N���͕����̈�̍ŏ�
        if(id == 0 || hash != sharedHash[threadIdx.x]){
            cell.dCellStart[hash] = id;
            if(id > 0){
                // ��O�̃p�[�e�B�N���́C��O�̕����̈�̍Ō�
                cell.dCellEnd[sharedHash[threadIdx.x]] = id;
            }
        }

        // �C���f�b�N�X���Ō�Ȃ�΁C�����̈�̍Ō�
        if(id == cell.uNumParticles-1){
            cell.dCellEnd[hash] = id+1;
        }

        // �ʒu�Ƒ��x�̃f�[�^����ёւ�
        // �\�[�g�����C���f�b�N�X�ŎQ�Ƃ��\�����T�����̃O���[�o���������A�N�Z�X���ɗ͗}���邽�߂Ƀf�[�^���̂��̂���ёւ���
        uint sid = cell.dSortedIndex[id];
        float3 pos = make_float3(dpos[DIM*sid+0], dpos[DIM*sid+1], dpos[DIM*sid+2]);
        float3 vel = make_float3(dvel[DIM*sid+0], dvel[DIM*sid+1], dvel[DIM*sid+2]);

        cell.dSortedPos[id] = pos;
        cell.dSortedVel[id] = vel;
    }
}

//---------------------------------------------
//�ȉ��A�C�V��ǉ� �l������float4(x,y,z,w)�Ƃ��Ĉ���(�o������Q�l)

//�o������̃R�[�h����
//��̎l�����̐ς����
__host__ __device__ float4 quatProduct(float4 a, float4 b) {
    return make_float4(
        a.x * b.w + a.w * b.x - a.z * b.y + a.y * b.z,
        a.y * b.w + a.z * b.x + a.w * b.y - a.x * b.z,
        a.z * b.w - a.y * b.x + a.x * b.y + a.w * b.z,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z
    );
}

//�o������̃R�[�h����
//�l�����̋��������
__host__ __device__ float4 quatConjugate(float4 quat) {
    return make_float4(-quat.x, -quat.y, -quat.z, quat.w);
}

//�o������̃R�[�h����
//3�����x�N�g���Ǝl�����̐ς����
__host__ __device__ float3 rotVecByQuat(float3 vec, float4 quat) {
    float4 vecq = make_float4(vec, 0.0f);
    float4 vecq_dash = quatProduct(quatProduct(quat, vecq), quatConjugate(quat));
    return make_float3(vecq_dash);
};

//3�����x�N�g���̒��������߂�
__host__ __device__ float Length(float3 vec) {
    return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

//4�����x�N�g���̒��������߂�
__host__ __device__ float Length(float4 vec) {
    return sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z + vec.w * vec.w);
}

//XPBD�ɗp����ɂ�S��0�ɂ���
//dlamb_ss:�L�ѐ���̃�
//dlamb_bt:�Ȃ�����̃�
__global__
void CxSetLambdaZero(float* dlamb_ss,float* dlamb_bt,int n) {
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    dlamb_ss[DIM * id] = dlamb_ss[DIM * id + 1] = dlamb_ss[DIM * id + 2] = 0.f;
    dlamb_bt[QUAT * id] = dlamb_bt[QUAT * id + 1] = dlamb_bt[QUAT * id + 2] = dlamb_bt[QUAT * id + 3] = 0.f;
}

//XPBD�̐L�сE����f����
//dpos:�ʒu
//dcurpos:�ʒu�X�V�O�̈ʒu(�����p�ɑ��x���l����ꍇ)
//dmas:����
//dlen:���
//dkss:�Ȃ�����
//dquat:�p��
//dcurquat:�X�V�O�̎p��(�����p�Ɋp���x���l����ꍇ)
//dlamb_ss:XPBD�̃�
//fix:�Œ�_���ǂ���(1�Ȃ�ΌŒ�_)
//dt:�^�C���X�e�b�v
//n:���q��
//odd_even:�����̃X���b�hID����̃X���b�hID���𔻕�
//example_flag:�`��ɂ���āC�������ꕔ�ύX����
__global__
void CxStretchingShearConstraint(float* dpos,float* dcurpos, float* dmas, float* dlen, float* dkss, float* dquat,float* dcurquat, float* dlamb_ss, int* dfix, float dt, int n,int odd_even,int iter,bool example_flag) {
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    id = id * 2 + odd_even;
    
    if (id >= n-1) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�) �Ō�̗��q�̓X�L�b�v
    if (dfix[id] == 1 || dfix[id + 1] == 1) return;//�e���̖��ɕӂ̐������s�� 

    float3 pos0 = make_float3(dpos[DIM * id], dpos[DIM * id + 1], dpos[DIM * id + 2]);
    float3 pos1 = make_float3(dpos[DIM * (id + 1)], dpos[DIM * (id + 1) + 1], dpos[DIM * (id + 1) + 2]);

    float mass = dmas[id];//���݂͗��q�̏d���͑S�ē������Ƃ��Ă���
    float length = dlen[id];
    float kss = dkss[id];
    float4 quat = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);//�l������(x,y,z,w)�ɂȂ�悤�Ɏ󂯎��
    float3 lambda_ss = make_float3(dlamb_ss[DIM * id], dlamb_ss[DIM * id + 1], dlamb_ss[DIM * id + 2]);

    float w1 = 1 / mass;//�d��
    float w2 = 1 / mass;//�d��
    float wq = 1.0;//�ӂ̏d��
    //�d�݂̓K�؂Ȑݒ�ɂ��Ă͂܂���܂��Ă��Ȃ�
    //�d�݂�2�̒��_�̎��ʂ̘a�Ƃ��Đݒ�
    wq = (w1 + w2) ;
    //�����𗘗p���Ă݂�
    //wq = (w1 + w2) / length / length * 4;
    //wq = w1 * w2 / (w1 + w2);
    //wq = 7.5e3;
    //Example Rod
    if (example_flag) {
        w1 = w2 = 1/length;
        wq = 1.0e5*length;
    }

    //�Œ�_�̏ꍇ�C�d�݂�0�ɂ���
    if (dfix[id-1] == 1) {
        w1 = 0;
        wq = 0;
    }

    float alpha = 1 / (kss * dt * dt);
    float3 ds3 = rotVecByQuat(make_float3(0.f, 0.f, 1.f), quat);//�p���̊�ƂȂ�x�N�g����z������O��Ƃ��Ă���

    //�������l����ꍇ-----------------------------
    //float damping = 0.1;//�����W��
    //float beta = damping * dt * dt;
    //float gamma = alpha * beta / dt;
    //float weight_with_damping = (1 + gamma) * ((1 / (length * length) * (w1 + w2)) + 4 * wq) + alpha;//(1+gamma)*(\nablaC^2 w1+\nablaC^2 w2+\nablaC^2 wq)+alpha(�ꕔ�ȗ������ăR�����g)

    //float3 curpos0 = make_float3(dcurpos[DIM * id], dcurpos[DIM * id + 1], dcurpos[DIM * id + 2]);
    //float3 curpos1 = make_float3(dcurpos[DIM * (id + 1)], dcurpos[DIM * (id + 1) + 1], dcurpos[DIM * (id + 1) + 2]);
    //float4 curquat = make_float4(dcurquat[QUAT * id], dcurquat[QUAT * id + 1], dcurquat[QUAT * id + 2], dcurquat[QUAT * id + 3]);

    //float3 tmp_v0 = pos0 - curpos0;
    //float3 tmp_v1 = pos1 - curpos1;

    //float3 tmp_angvel = 2.f * make_float3(quatProduct(quatConjugate(curquat), quat));
    //float4 nablaC_q = -2.f * quatProduct(quat, make_float4(0.f, 0.f, -1.f, 0.f));//q_s*\bar{e_3}

    //float3 sum = tmp_v0 * (-1 / length) + tmp_v1 * (1 / length) + make_float3(quatProduct(nablaC_q,make_float4(tmp_angvel,0.f)));//\nablaC*v���v�Z�C�ŏI�e��\nablaC*(x_i-x^n)���l�����ɒu��������

    //float3 mole_with_damping = (pos1 - pos0) / length - ds3 + alpha * lambda_ss + gamma * sum;//�}�C�i�X�������Ă��Ȃ�
    //---------------------------------------------

    float weight = w1 + w2 + length * length * (4 * wq + alpha);//����
    float3 mole = length * (pos0 - pos1 + length * ds3 - alpha * length * lambda_ss);//���q
    float3 lambda = mole / weight;//����

    //�������l����ꍇ------------------------------------------
    //lambda = -mole_with_damping / weight_with_damping;
    //----------------------------------------------------------

    //�d�ݒǉ�
    float3 delta_pos0 = -w1*lambda / length;//��x0(�_���ƕ����t)
    float3 delta_pos1 = w2*lambda / length;//��x1(�_���ƕ����t
    float4 q_e3_bar = quatProduct(quat, make_float4(0.f, 0.f, -1.f, 0.f));
    
    float4 inter_quat = quatProduct(make_float4(lambda, 0.f), q_e3_bar);
    float4 delta_quat = -wq * 2.f * inter_quat;
    float4 new_quat = quat + delta_quat;//�X�V���qs
    new_quat = normalize(new_quat);

    //�ʒu�X�V
    if (dfix[id-1] == 0) {
        dpos[id * DIM] += delta_pos0.x;
        dpos[id * DIM + 1] += delta_pos0.y;
        dpos[id * DIM + 2] += delta_pos0.z;
    }

    dpos[(id + 1) * DIM] += delta_pos1.x;
    dpos[(id + 1) * DIM + 1] += delta_pos1.y;
    dpos[(id + 1) * DIM + 2] += delta_pos1.z;
    
    //�p���̍X�V
    if (dfix[id] == 0) {
        dquat[id * QUAT] = new_quat.x;
        dquat[id * QUAT + 1] = new_quat.y;
        dquat[id * QUAT + 2] = new_quat.z;
        dquat[id * QUAT + 3] = new_quat.w;
    }

    //XPBD�ł̃ɂ̐ݒ�
    dlamb_ss[id * DIM] += lambda.x;
    dlamb_ss[id * DIM + 1] += lambda.y;
    dlamb_ss[id * DIM + 2] += lambda.z;
}

//�Ȃ��˂��ꐧ��̕���������
//delta_omega=cur_omega-rest_omega
//delta_omega_plus=cur_omega+rest_omega
__host__ __device__ 
int deltaOmegaSign(float4 delta_omega, float4 delta_omega_plus) {
    if (dot(delta_omega, delta_omega) > dot(delta_omega_plus, delta_omega_plus)) return -1;

    return 1;
}

//�Ȃ��˂��ꐧ��̒ǉ�
//dmas:����
//dquat:�ӂ̎p��(�l����)
//dcurquat:�X�V�O�̎p��(�����p�ɑ��x���l����)
//domega:��_���{�[�x�N�g��
//dkbt:�Ȃ�����
//dlamb_bt:�Ȃ�����ɗp�����
//dfix:�Œ�_���ǂ�����\��(1�Ȃ�ΌŒ�_)
//dt:�^�C���X�e�b�v
//n:���q��
//odd_even:�����̃X���b�hID����̃X���b�hID���𔻕�
//example_flag:�`��ɂ���āC�������ꕔ�ύX����
__global__ 
void CxBendTwistConstraint(float* dmas,float* dquat,float* dcurquat, float* domega, float* dkbt, float* dlamb_bt, float* dlength, int* dfix, float dt, int n, int odd_even, int iter,bool example_flag) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    id = id * 2 + odd_even;
    if (id >= n - 2) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�) �Ō��2���q�̓X�L�b�v
    if (dfix[id + 1] == 1 || dfix[id + 2] == 1) return;//���[�̗��q�Ƃ��̂ЂƂO�̗��q�ł́A�G�b�W������Ȃ�

    float kbt = dkbt[id];
    float4 quat1 = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);
    float4 quat2 = make_float4(dquat[QUAT * (id + 1)], dquat[QUAT * (id + 1) + 1], dquat[QUAT * (id + 1) + 2], dquat[QUAT * (id + 1) + 3]);
    float4 rest_omega = make_float4(domega[QUAT * id], domega[QUAT * id + 1], domega[QUAT * id + 2], domega[QUAT * id + 3]);
    float4 lambda_bt = make_float4(dlamb_bt[QUAT * id], dlamb_bt[QUAT * id + 1], dlamb_bt[QUAT * id + 2], dlamb_bt[QUAT * id + 3]);

    float dlen1 = dlength[id];
    float dlen2 = dlength[id + 1];

    float wq1 = 1.0f;//1/(1.0e-3)*dlen1
    float wq2 = 1.0f;
    float alpha = 1 / (kbt * dt * dt);
    //�d�݂̐ݒ�ɗ��p(�S�Ă̎��ʂ͓������Ɖ���)
    float mass = dmas[id];
    //�d�݂̓K�؂Ȑݒ�͂܂���܂��Ă��Ȃ�
    //wq1 = 1.0f / dlen1;
    //wq2 = 1.0f / dlen2;
    wq1 = wq2 = 2 / mass * 10.f;//*10.f
    //wq1 = 2 / mass / dlen1 / dlen1 * 4;
    //wq2 = 2 / mass / dlen2 / dlen2 * 4;
    if (example_flag) {
        wq1 = wq2 = 1.0e5 * dlen1;
    }

    //�Œ肷��G�b�W�Ȃ�d�݂�0�ɂ���
    if (dfix[id] == 1) {
        wq1 = 0;
    }

    float weight = wq1 + wq2 + alpha;//����
    float4 cur_omega = quatProduct(quatConjugate(quat1), quat2);
    int s = deltaOmegaSign(cur_omega - rest_omega, cur_omega + rest_omega);//���������߂�

    //�������l����ꍇ-------------------------------------------------
    //float damping = 0.05;//�����W��
    //float beta = damping * dt * dt;
    //float gamma = alpha * beta / dt;
    //float weight_with_damping = (1 + gamma) * (wq1 + wq2) + alpha;

    //float4 curquat1 = make_float4(dcurquat[QUAT * id], dcurquat[QUAT * id + 1], dcurquat[QUAT * id + 2], dcurquat[QUAT * id + 3]);
    //float4 curquat2 = make_float4(dcurquat[QUAT * (id + 1)], dcurquat[QUAT * (id + 1) + 1], dcurquat[QUAT * (id + 1) + 2], dcurquat[QUAT * (id + 1) + 3]);

    //float4 tmp_angvel1 = 2.f * quatProduct(quatConjugate(curquat1), quat1);
    //float4 tmp_angvel2 = 2.f * quatProduct(quatConjugate(curquat2), quat2);

    //float4 sum = tmp_angvel1 * (-quat2) + tmp_angvel2 * quat1;

    //float4 mole_with_damping = -(cur_omega - s * rest_omega) - alpha * lambda_bt - gamma * sum;//�}�C�i�X�����łɂ����Ă���
    //-----------------------------------------------------------------

    float4 delta_omega = cur_omega - (s * rest_omega);
    //delta_omega.w = 0.f;//omega.w=0�Ƃ��Ă���
    float4 lambda = (-delta_omega - alpha * lambda_bt) / weight;
    //lambda.w = 0.f;

    //�������l����ꍇ------------------------------------------------
    //lambda = mole_with_damping / weight_with_damping;
    //----------------------------------------------------------------

    //�d�ݒǉ�
    float4 delta_quat1 = wq1 * quatProduct(quat2, quatConjugate(lambda));
    float4 delta_quat2 = wq2 * quatProduct(quat1, lambda);
    float4 new_quat1 = normalize(quat1 + delta_quat1);
    float4 new_quat2 = normalize(quat2 + delta_quat2);

    //quat1�̍X�V
    if (dfix[id] == 0) {
        dquat[QUAT * id] = new_quat1.x;
        dquat[QUAT * id + 1] = new_quat1.y;
        dquat[QUAT * id + 2] = new_quat1.z;
        dquat[QUAT * id + 3] = new_quat1.w;
    }
    //quat2�̍X�V
    dquat[QUAT * (id + 1)] = new_quat2.x;
    dquat[QUAT * (id + 1) + 1] = new_quat2.y;
    dquat[QUAT * (id + 1) + 2] = new_quat2.z;
    dquat[QUAT * (id + 1) + 3] = new_quat2.w;

    //lambda�̍X�V
    dlamb_bt[QUAT * id] += lambda.x;
    dlamb_bt[QUAT * id + 1] += lambda.y;
    dlamb_bt[QUAT * id + 2] += lambda.z;
    dlamb_bt[QUAT * id + 3] += lambda.w;
}

//�Փː���̎���
//dpos:�ʒu
//dvel:���x
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//center:�є��Ƃ̏Փ˂������������̒��S
//rad:�є��Ƃ̏Փ˂������������̔��a
//dt:�^�C���X�e�b�v
//n:���q��
__global__
void CxCollisionConstraint(float* dpos, float* dvel, int* dfix, float3 center, float rad, float dt, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n)return;
    if (dfix[id] == 1) return;

    float3 pos = make_float3(dpos[DIM * id], dpos[DIM * id + 1], dpos[DIM * id + 2]);
    float3 vel = make_float3(dvel[DIM * id], dvel[DIM * id + 1], dvel[DIM * id + 2]);

    float3 d = pos - center;
    float l = rad - length(d);
    if (l <= 0) return;

    float3 norm = normalize(d);

    dpos[DIM * id] += l * norm.x;
    dpos[DIM * id + 1] += l * norm.y;
    dpos[DIM * id + 2] += l * norm.z;

    float3 addVel = -dot(norm, vel) * norm;

    dvel[DIM * id] += addVel.x;
    dvel[DIM * id + 1] += addVel.y;
    dvel[DIM * id + 2] += addVel.z;
}

//���Ԑϕ�
//dpos:�X�V��̈ʒu
//dcurpos:�X�V�O�̈ʒu
//dvel:���x(�X�V�������x����)
//dt:�^�C���X�e�b�v
//n:���q��
__global__
void CxIntegrate(float* dpos, float* dcurpos,float* dvel,float dt,int n,bool vel_control) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    float3 pos = make_float3(dpos[DIM * id], dpos[DIM * id + 1], dpos[DIM * id + 2]);
    float3 cur_pos = make_float3(dcurpos[DIM * id], dcurpos[DIM * id + 1], dcurpos[DIM * id + 2]);
    float3 vel = (pos-cur_pos)/dt;
    
    //���Ԑϕ��̍ۂɁC��葬�x�ȉ��̗��q�ɂ��Ă͕ω��Ȃ��Ƃ��āC�Œ肵�Ă��܂�
    if (vel_control&&length(vel) < VEL_EPSILON) {
        //���x��0�ɌŒ�
        vel = make_float3(0.f);
        //�ʒu���X�V�O�̒l�ɖ߂�
        pos = cur_pos;
        dpos[DIM * id] = pos.x;
        dpos[DIM * id + 1] = pos.y;
        dpos[DIM * id + 2] = pos.z;
    }

    //if (length(vel) > 1.0e-2) printf("id %d vel x:%f,y:%f,z:%f\n", id, vel.x, vel.y, vel.z);

    dvel[DIM * id] = vel.x;
    dvel[DIM * id + 1] = vel.y;
    dvel[DIM * id + 2] = vel.z;

    dcurpos[DIM * id] = pos.x;
    dcurpos[DIM * id + 1] = pos.y;
    dcurpos[DIM * id + 2] = pos.z;
}

//�O�͌v�Z
//����d�͂��C���[�W���������x��^����
//�S�Ă̗��q�ɓ��������x��^����
//�f�o�b�N�p
__global__
void CxCalExternalForces(float* dpos,float* dvel,float* dmass,int* dfix,float3 gravity, float3 wind,float dt, int n){
int id = blockDim.x * blockIdx.x + threadIdx.x;
if (id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�){
if (dfix[id] == 1) return;
if (dfix[id - 1] == 1)return;

float mass = dmass[id];

dvel[DIM * id] += (gravity.x+wind.x) * dt;
dvel[DIM * id + 1] += (gravity.y+wind.y) * dt;
dvel[DIM * id + 2] += (gravity.z+wind.z) * dt;

dpos[DIM * id] += dvel[DIM * id] * dt;
dpos[DIM * id + 1] += dvel[DIM * id + 1] * dt;
dpos[DIM * id + 2] += dvel[DIM * id + 2] * dt;
}

//�ʒu�x�[�X�@�̐L�ѐ���
//�f�o�b�N�p
__global__
void CxStretchingConstraint(float* dpos, float* dmas, float* dlen, float* dkss, float* dquat, int* dfix, int n,int odd_even) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    id = id * 2 + odd_even;

    if (id >= n - 1) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�) �Ō�̗��q�̓X�L�b�v
    if (dfix[id + 1] == 1) return;//�e���̖��ɕӂ̐������s��

    float3 pos0 = make_float3(dpos[DIM * id], dpos[DIM * id + 1], dpos[DIM * id + 2]);
    float3 pos1 = make_float3(dpos[DIM * (id + 1)], dpos[DIM * (id + 1) + 1], dpos[DIM * (id + 1) + 2]);

    float mass = dmas[id];//���݂͗��q�̏d���͑S�ē������Ƃ��Ă���
    float length = dlen[id];
    float kss = dkss[id];
    float4 quat = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);//�l������(x,y,z,w)�ɂȂ�悤�Ɏ󂯎��

    float3 e3 = make_float3(0, 0, 1);
    float3 d3 = rotVecByQuat(e3, quat);

    float w = 1.0 / mass;//�����p
    float wq = 1.0e5;//�����p
    wq = 1 / length;

    //��pos���v�Z
    float3 gamma = (pos1 - pos0) / length - d3;
    gamma /= (2*w) / length + 4.0f * wq * length;
    float ks = 1.0;
    gamma *= ks;

    float3 delta_pos0 = w * gamma;
    float3 delta_pos1 = -w * gamma;

    // calc delta_q
    float4 e3q = make_float4(e3, 0.0f);
    float4 q_e3_bar = quatProduct(quat, quatConjugate(e3q)); // calc q*e3_bar
    float4 gammaq = make_float4(gamma, 0.0f);
    float4 inter_quat = quatProduct(gammaq, q_e3_bar);
    float4 delta_quat = wq * length * quatProduct(gammaq, q_e3_bar);//2.0�ǉ�

    float4 new_quat = quat + delta_quat;//�X�V���qs
    new_quat = normalize(new_quat);

    if (dfix[id] == 0) {
        //printf("dpos %f", Length(delta_pos0));
        dpos[id * DIM] += delta_pos0.x;
        dpos[id * DIM + 1] += delta_pos0.y;
        dpos[id * DIM + 2] += delta_pos0.z;
    }

    dpos[(id + 1) * DIM] += delta_pos1.x;
    dpos[(id + 1) * DIM + 1] += delta_pos1.y;
    dpos[(id + 1) * DIM + 2] += delta_pos1.z;

    dquat[id * QUAT] = new_quat.x;
    dquat[id * QUAT + 1] = new_quat.y;
    dquat[id * QUAT + 2] = new_quat.z;
    dquat[id * QUAT + 3] = new_quat.w;
}

//�ʒu�A���x�A�����x���o��
//�f�o�b�N�p
__global__
void CxPrint3Dfloat(float* dpos, float* dvel, float* dacc, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return;

    printf("id %d\n", id);
    printf("pos x:%f,y:%f,z:%f\n", dpos[id * DIM], dpos[id * DIM + 1], dpos[id * DIM + 2]);
    printf("vel x:%f,y:%f,z:%f\n", dvel[id * DIM], dvel[id * DIM + 1], dvel[id * DIM + 2]);
    printf("acc x:%f,y:%f,z:%f\n", dacc[id * DIM], dacc[id * DIM + 1], dacc[id * DIM + 2]);
}

//�ڐ��̍X�V
//dpos:�ʒu
//dtang:�G�b�W���Ƃ̐ڐ�
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//n:���q��
__global__
void CxTangUpdate(float* dpos, float* dtang, int* dfix, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return;
    if (dfix[id] == 1)return;

    dtang[DIM * id] = dpos[DIM * id] - dpos[DIM * (id - 1)];
    dtang[DIM * id + 1] = dpos[DIM * id + 1] - dpos[DIM * (id - 1) + 1];
    dtang[DIM * id + 2] = dpos[DIM * id + 2] - dpos[DIM * (id - 1) + 2];
}

//�p���x�ȂǏ�������f�o�C�X�������ɐݒ�����Ă�����̂̏����l��0�ɐݒ�
//dangvel:�p���x
//dfss:�G�b�W���Ƃɂ������(GlobalForceStep�ŋ��߂�)
//dpbf_lambda:���x����̌v�Z�ߒ��ɕK�v�ȃɂ��������m��
//n:���q��
__global__
void CxSetParametersZero(float* dangvel,float* dfss,float* dpbf_lambda, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return;
    dangvel[DIM * id] = dangvel[DIM * id + 1] = dangvel[DIM * id + 2] = 0;
    dfss[DIM * id] = dfss[DIM * id + 1] = dfss[DIM * id + 2] = 0;
    dpbf_lambda[id] = 0.f;
}

//�p�����x�̍X�V
//dangvel:�p���x
//dquat:�p��(�l����)
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//dt:�^�C���X�e�b�v
//n:���q��
__global__
void CxAngVelUpdate(float* dangvel, float* dquat,int* dfix,float dt, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n-1) return;
    if (dfix[id + 1] == 1) return;

    float4 quat = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);
    float3 avel = make_float3(dangvel[DIM * id], dangvel[DIM * id + 1], dangvel[DIM * id + 2]);

    float4 avelq = make_float4(avel, 0.0f);
    quat += 0.5f * quatProduct(quat, avelq) * dt;
    quat = normalize(quat);

    dquat[QUAT * id] = quat.x;
    dquat[QUAT * id + 1] = quat.y;
    dquat[QUAT * id + 2] = quat.z;
    dquat[QUAT * id + 3] = quat.w;
}

//�e�����x�̎��Ԑϕ�
//dangvel:�p���x
//dcurquat:�O�X�e�b�v�̎p��(�ʒu�C���O)
//dquat:���݂̎p��(�ʒu�C����)
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//dt:�^�C���X�e�b�v
//n:���q��
//vel_control:�p���x�����ȉ��Ȃ�؂�̂Ă��s�����ǂ����𔻒�
__global__
void CxAngVelIntegrate(float* dangvel,float* dcurquat, float* dquat,int* dfix,float dt, int n,bool vel_control) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n-1) return;
    if (dfix[id + 1] == 1) return;

    float4 quat = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);
    float4 cur_quat = make_float4(dcurquat[QUAT * id], dcurquat[QUAT * id + 1], dcurquat[QUAT * id + 2], dcurquat[QUAT * id + 3]);

    float4 delta_rot = quatProduct(quatConjugate(cur_quat), quat);
    float3 new_AngVel = 2.0f * make_float3(delta_rot.x, delta_rot.y, delta_rot.z) / dt;

    //�p���x�����ȉ��Ȃ�C�����Ă��Ȃ��Ƃ��ČŒ肷��
    if (vel_control&&length(new_AngVel) < ANGVEL_EPSILON) {
        new_AngVel = make_float3(0.f);
        quat = cur_quat;
        dquat[QUAT * id] = quat.x;
        dquat[QUAT * id + 1] = quat.y;
        dquat[QUAT * id + 2] = quat.z;
        dquat[QUAT * id + 3] = quat.w;
    }

    //if (length(new_AngVel) > 1.0e-2)printf("id %d angvel x:%f,y:%f,z:%f\n", id, new_AngVel.x, new_AngVel.y, new_AngVel.z);

    dangvel[DIM * id] = new_AngVel.x;
    dangvel[DIM * id + 1] = new_AngVel.y;
    dangvel[DIM * id + 2] = new_AngVel.z;

    dcurquat[QUAT * id] = quat.x;
    dcurquat[QUAT * id + 1] = quat.y;
    dcurquat[QUAT * id + 2] = quat.z;
    dcurquat[QUAT * id + 3] = quat.w;
}

//����x���ŏ��Ɍv�Z�������x�ɐݒ�
//dpos:�ʒu
//dRestDens:���q���Ƃɐݒ肷�����x
//dvol:�̐�
//dmas:����
//n:���q��
__global__
void CxRestDensSet(float* dpos,float* dRestDens,float* dvol, float*dmas,int n) {
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    float3 pos0 = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    //float m = params.mass;
    float a = params.aw;
    float rest_dens = params.rest_dens;

    // ���q�𒆐S�Ƃ��Ĕ��ah���Ɋ܂܂��O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos0 - make_float3(h));
    grid_pos1 = calcGridPos(pos0 + make_float3(h));

    // ���͂̃O���b�h�Z�����܂߂ċߖT�T�����Ė��x�v�Z
    float dens = 0.0f;
    for (int z = grid_pos0.z; z <= grid_pos1.z; ++z) {
        for (int y = grid_pos0.y; y <= grid_pos1.y; ++y) {
            for (int x = grid_pos0.x; x <= grid_pos1.x; ++x) {
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // �O���b�h�n�b�V���l

                // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                uint startIndex = params.cell.dCellStart[ghash];
                if (startIndex != 0xffffffff) {	// �Z������łȂ����̃`�F�b�N
                    // �Z�����̃p�[�e�B�N���Ŕ���
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for (uint j = startIndex; j < endIndex; ++j) {
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos1 = params.cell.dSortedPos[j];
                        float3 rij = pos0 - pos1;
                        float r = length(rij);
                        if (r <= h) {
                            // Poly6�J�[�l���Ŗ��x���v�Z (rho = �� m Wij)
                            float q = h * h - r * r;
                            // ���̗��q��rest_dens*dvol[sj] = m�ƂȂ�悤�ɐݒ肵�Ă���
                            // ���E���q�͑z�̐ςƏ������x���畡���w���E���q���������ꍇ�̉��z���ʃ�=��0*Vb�����߂Ďg�� 
                            //float m = rest_dens * dvol[sj];
                            //printf("mass %f\n", m);
                            float m = params.mass;
                            dens += m * a * q * q * q;
                        }
                    }
                }
            }
        }
    }


    // �v�Z�������x���f�o�C�X�������ɏ�������
    uint sid = params.cell.dSortedIndex[id];
    //�Œ���̖��x��ݒ�
    dRestDens[sid] = fmaxf(dens,rest_dens);//�o������̐ݒ���@
}

//�ꗥ�̊�ƂȂ閧�x�̐ݒ�
//�f�o�b�N�p
__global__
void CxRestTotalDens(float* drestdens, float dens, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    drestdens[id] = dens;
}

//�ȉ�SagFree�̏������ڐA--------------------------------------------------------------------------------------------------------------------------------
//LU������p���āA�A��1��������������
//1�Ȃ琬��,0�Ȃ玸�s
//���݁C�g�p���Ă��Ȃ�
__device__ __host__
int LUDecomp(float A[][4], int n) {
    if (n <= 0) return 0;

    for (int i = 0; i < n; i++) {
        //l_ij�̌v�Z(i>=j)
        for (int j = 0; j <= i; ++j) {
            float lu = A[i][j];
            for (int k = 0; k < j; k++) {
                lu -= A[i][k] * A[k][j];//l_ik*u_kj
            }
            A[i][j] = lu;
        }

        //u_ij�̌v�Z(i<j)
        for (int j = i + 1; j < n; ++j) {
            float lu = A[i][j];
            for (int k = 0; k < i; ++k) {
                lu -= A[i][k] * A[k][j];
            }
            A[i][j] = lu / A[i][i];
        }
    }

    return 1;
}

//A:LU�������ꂽ�s��
//b:�E�Ӄx�N�g��
//x:���ʃx�N�g��
//n:�s��̑傫��
//���݁C�g�p���Ă��Ȃ�
__device__ __host__
int LUSolver(const float A[][4], const float b[], float x[], int n) {
    if (n <= 0) return 0;

    //�O�i���
    //LY=b����Y���v�Z
    for (int i = 0; i < n; ++i) {
        float bly = b[i];
        for (int j = 0; j < i; ++j) {
            bly -= A[i][j] * x[j];
        }
        //if (A[i][i] < 1.0e-6) printf("A trace error!\n");
        x[i] = bly / A[i][i];
    }

    //��ޑ��
    //UX=Y����X���v�Z
    for (int i = n - 1; i >= 0; --i) {
        float yux = x[i];
        for (int j = i + 1; j < n; ++j) {
            yux -= A[i][j] * x[j];
        }
        x[i] = yux;
    }

    return 1;
}

//3�����x�N�g������l���������Ƃ߂�
//z��������ƂȂ���x�N�g���Ƃ���
__device__ __host__
float4 quatFromDirector(float3 d3) {
    d3 = normalize(d3);
    float3 e3 = make_float3(0.f, 0.f, 1.0f);//z������O��
    float3 w = cross(e3, d3);
    float4 q = make_float4(w, dot(e3, d3));
    q.w += length(q);
    return normalize(q);
}

//��̃x�N�g������l���������߂�
//�Е��͊�ƂȂ���
__device__ __host__
float4 quatFromTwoVectors(float3 a, float3 b) {
    float3 v0 = normalize(a);
    float3 v1 = normalize(b);
    float c = dot(v0, v1);

    float3 axis = cross(v0, v1);
    float s = sqrt((1 + c) * 2);
    float3 vec = axis / s;
    float w = s * 0.5;

    float4 quat = make_float4(vec.x,vec.y,vec.z,w);
    return normalize(quat);
}

//�Î~���C�ƂȂ�͂����߂�
//�O���[�o���X�e�b�v�ŗp����\��
//���ꂼ��̗��q������������������Ȃ��Ɩ��C������������������܂�Ȃ����߁C�K���ȕ����������Ƃ���
//���݂́C�ЂƂ܂��C�S�Ă̗��q�����ΓI�ɉ����������ɓ������̂Ƃ���(����Ȃ��Ƃ͂��蓾�Ȃ���)
__device__
float3 calcFrictionForce(float3 dir_i, float3 dir_j, float* ddens, float* drestdens, float* dvol, int id) {
    float3 pos_i = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    //�C���f�b�N�X�̌v�Z
    uint sid = params.cell.dSortedIndex[id];
    //�Î~���C�W��(�����C�W���͐Î~���C�W����0.1�{�Ƃ���)
    float mu = MU;

    // ���q�𒆐S�Ƃ��Ĕ��ah���Ɋ܂܂��O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos_i - make_float3(h));
    grid_pos1 = calcGridPos(pos_i + make_float3(h));

    //�ŏI�I�Ȗ��C��
    float3 x_fric = make_float3(0.f);

    for (int z = grid_pos0.z; z <= grid_pos1.z; ++z) {
        for (int y = grid_pos0.y; y <= grid_pos1.y; ++y) {
            for (int x = grid_pos0.x; x <= grid_pos1.x; ++x) {
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // �O���b�h�n�b�V���l
                // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                uint startIndex = params.cell.dCellStart[ghash];
                if (startIndex != 0xffffffff) {	// �Z������łȂ����̃`�F�b�N
                    // �Z�����̃p�[�e�B�N���Ŕ���
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for (uint j = startIndex; j < endIndex; ++j) {
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos_j = params.cell.dSortedPos[j];

                        //j�̗��q�̏������x
                        float restdens_j = drestdens[sj];
                        //j�̑̐�
                        float vol_j = dvol[sj];
                        //����
                        float m = restdens_j * vol_j;

                        float3 r_ij = pos_i - pos_j;
                        float r = length(r_ij);
                        if (r <= 1.0e-3) continue;
                        if (r < h) {
                            float3 v_ij = dir_i - dir_j;
                            //v_ij = make_float3(0, 1, 0);//����������Ɉ��������Ă���Ƃ���
                            v_ij = dir_i;

                            r_ij = normalize(r_ij);
                            //�Փ˖@���ɑ΂��Đ����Ȑ��������߂�
                            float3 delxn = v_ij - r_ij * dot(v_ij, r_ij);//delta x_��=v_ij-x_||

                            float q = h * h - r * r;//(h^2-||rij||^2)
                            x_fric += m / ddens[sj] * MU * delxn * params.aw * q * q * q;//aw*q^3
                        }
                    }
                }
            }
        }
    }
    return x_fric;
}

__device__ __host__
float3 CalcNormalTorque(float3 pos0, float3 pos1, float4 quat, float3 fss, float len, float mass0, float mass1, float3 gravity) {
    float3 mid = (pos0 + pos1) / 2.f;
    //�����x�N�g�������߂�
    float3 dir = rotVecByQuat(make_float3(0, 0, 1), quat);
    //0,1�̎��_�ɑ΂���d�S(�G�b�W����)����̃x�N�g��
    float3 r0 = normalize(-dir + mid) * len / 2;
    float3 r1 = normalize(dir + mid) * len / 2;

    //���͂ɂ��G�b�W�ɂ�����͂𗼕��̎��_�ƊO�ς����
    float3 tau_int0 = cross(r0, fss);
    float3 tau_int1 = cross(r1, fss);

    //���͂ɂ��g���N���o��
    float3 tau_internal = tau_int0 + tau_int1;
    //printf("id %d tau_internal x:%f,y:%f,z:%f\n",id,tau_internal.x,tau_internal.y,tau_internal.z);

    //�O��(�d��)�ɂ��g���N�����߂�
    float3 tau_ext0 = cross(r0, mass0 * gravity);
    float3 tau_ext1 = cross(r1, mass1 * gravity);

    //�O�͂ɂ��g���N���o��
    float3 tau_external = tau_ext0 + tau_ext1;
    //printf("id %d tau_external x:%f,y:%f,z:%f\n", id, tau_external.x, tau_external.y, tau_external.z);

    float3 total_torque = tau_internal + tau_external;
    //if (length(total_torque) > 1.0e-3) printf("id %d total_torque x:%f,y:%f,z:%f\n", id, total_torque.x, total_torque.y, total_torque.z);

    return total_torque;
}

//�O���[�o���t�H�[�X�X�e�b�v
//�ЂƂ܂����C�����l�����ɁA�P���ɏd�͂��狁�߂邱�ƂƂ���
//�G�b�W����ɖ�������
//dfss:�G�b�W���Ƃɂ������
//dmass:����
//last_index:�є����Ƃ̍Ō�̗��q�̃C���f�b�N�X���i�[
//gravity:�d��
//num_elastic:�����ł́C�є����Ƃɕ���v�Z���邽�߁C�є��̐���n��
__global__
void CxGlobalForceStep(float* dpos,float* dfss,float* dmass,int* dlast_ind,float3 gravity,float* ddens,float* drestdens,float* dvol,int num_elastic) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    if (id >= num_elastic) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    int min;
    if (id == 0) min = 1;
    else min = dlast_ind[id - 1] + 2;//�Œ�_�ɗאڂ���G�b�W�͌v�Z�Ɋ܂܂Ȃ�
    int max = dlast_ind[id] - 1;
    for (int i = max; i > min - 1; i--) {//i=0�̎��͗�0
        float mass = dmass[i+1];

        float3 pos0 = make_float3(dpos[i * DIM], dpos[i * DIM + 1], dpos[i * DIM + 2]);
        float3 pos1 = make_float3(dpos[(i + 1) * DIM], dpos[(i + 1) * DIM + 1], dpos[(i + 1) * DIM + 2]);

        float3 dir = normalize(pos1 - pos0);

        float3 prev_fss;
        if (i == max) prev_fss = make_float3(0.f);
        else prev_fss = make_float3(dfss[(i + 1) * DIM], dfss[(i + 1) * DIM + 1], dfss[(i + 1) * DIM + 2]);

        dfss[i * DIM] = -(mass * gravity.x - prev_fss.x);
        dfss[i * DIM + 1] = -(mass * gravity.y - prev_fss.y);
        dfss[i * DIM + 2] = -(mass * gravity.z - prev_fss.z);
        
        //11/16�ǉ�
        //�G�b�W�̕����ɂ݈̂��������Ă݂�--------------------
       /* float a = -(mass * gravity.y - prev_fss.y) / dir.y;

        dfss[i * DIM] = a * dir.x;
        dfss[i * DIM + 1] = a * dir.y;
        dfss[i * DIM + 2] = a * dir.z;*/
        //----------------------------------------------------

        //float3 dir_i = make_float3(0, -1, 0);
        //float3 dir_j = make_float3(0, -1, 0);
        //����,dir_i��dir_j�ɈӖ��͂Ȃ�
        //float3 friction_force = 0.5*calcFrictionForce(dir_i, dir_j, ddens, drestdens, dvol,id);
        
        //���C���܂߂čl����
        /*dfss[i * DIM] = -(mass * gravity.x - prev_fss.x - friction_force.x);
        dfss[i * DIM + 1] = -(mass * gravity.y - prev_fss.y - friction_force.y);
        dfss[i * DIM + 2] = -(mass * gravity.z - prev_fss.z - friction_force.z);*/

        //�l��ς���ꍇ
        /*if (i == max) {
            dfss[i * DIM] *= 20;
            dfss[i * DIM + 1] *= 20;
            dfss[i * DIM + 1] *= 20;
        }*/
    }
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
__global__
void CxLocalForceStep(float* dpos, float* dlen, float* dquat,float* dcurquat, float* dkss, float* dfss, int* dfix, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    if (id >= n - 1) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�) �Ō�̗��q�̓X�L�b�v
    if (dfix[id] == 1 || dfix[id + 1] == 1) return;//�e���̖��ɕӂ̐������s���A�ŏ��̃G�b�W�͌Œ肵�Ĉ���

    float3 pos1 = make_float3(dpos[DIM * id], dpos[DIM * id + 1], dpos[DIM * id + 2]);
    float3 pos2 = make_float3(dpos[DIM * (id + 1)], dpos[DIM * (id + 1) + 1], dpos[DIM * (id + 1) + 2]);

    float3 fss = make_float3(dfss[DIM * id], dfss[DIM * id + 1], dfss[DIM * id + 2]);

    float fs_Len = Length(fss);

    float l0 = dlen[id];
    float tmp_ks = dkss[id];
    float4 quat = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);//�l������(x,y,z,w)�ɂȂ�悤�Ɏ󂯎��

    float3 e3 = make_float3(0, 0, 1);
    float3 d3 = normalize(rotVecByQuat(e3, quat));

    //���ʎ��𖞂������ǂ����ɂ���āC�L�э���������������----------------
    float delta = 10.f;
    float B = dot((pos1 - pos2), fss) / (tmp_ks)+1;
    float AC = dot(pos1 - pos2, pos1 - pos2) * fs_Len * fs_Len / (tmp_ks * tmp_ks);//Length2��dot�ŕ\��
    float discrim = B * B - 4 * AC;
    //�����ŋ��߂�ꍇ
    /*while (1) {
        if (discrim > 1.0e-2) break;
        tmp_ks += delta;
    }*/

    //���ʎ��𖞂����Ȃ��ꍇ�A�莮�I��kss�����߂�
    if (discrim < 0) {
        //2������苁�߂Ă���A�{��+-sqrt(...)�ł��邪�Akss>0���A+�݂̂��l����
        tmp_ks = - abs(dot(fss, pos1 - pos2)) + 2 * fs_Len * Length(pos1 - pos2);
        printf("id %d discrim<0!!\n", id);
    }
    
    //�L�э������X�V
    dkss[id] = tmp_ks;

    //��̌��ƂȂ钷�������߂�(float�^�̂܂܂��Ɛ��x�ɉe�������邽�߁Adouble�^�ɕύX)
    double a = fs_Len * fs_Len;
    double b = -tmp_ks * (2.f * dot(fss, pos1 - pos2) + tmp_ks);
    double c = tmp_ks * tmp_ks * dot(pos1 - pos2, pos1 - pos2);//Length2��dot�ŕ\��

    //�����2�̌��
    double l1 = sqrt((-b + sqrt(abs(b * b - 4.f * a * c))) / (2.f * a));
    double l2 = sqrt((-b - sqrt(abs(b * b - 4.f * a * c))) / (2.f * a));

    //2�̌��̂����A��茻�݂̒����ɋ߂����̂�I��
    float rest_length;
    if (abs(l1 - l0) > abs(l2 - l0) && abs(l2) > 1.0e-10) rest_length = l2;
    else if (abs(l1) > 1.0e-10) rest_length = l1;
    else {
        printf("Error Occured! in LocalForceStep!");
    }

    //�p���̍X�V
    float3 new_ds = (rest_length * fss) / tmp_ks - (pos1 - pos2) / rest_length;

    //�l�������x�N�g�����狁�߂�(��̎�@�����邪�A�قƂ�Ǖς��Ȃ��Ɛ���)
    float4 new_qs = quatFromDirector(new_ds);
    //float4 new_qs = quatFromTwoVectors(e3, new_ds);

    //�p���x�N�g���̊m�F
    float3 d3_from_qs = rotVecByQuat(e3, new_qs);

    //�͂̊m�F�p(Eq.14�̏�ł́Cds3��-�Œ�`����Ă��邪�C�����+���ƍl������D)
    float3 new_Fss = (tmp_ks / rest_length) * ((pos1 - pos2) / rest_length + d3_from_qs);

    //printf("id %d old_quat x:%f,y:%f,z:%f,w:%f new_quat x:%f,y:%f,z:%f,w:%f\n", id, quat.x, quat.y, quat.z, quat.w, new_qs.x, new_qs.y, new_qs.z, new_qs.w);

    //����̍X�V
    dlen[id] = rest_length;
    //�p���̍X�V
    dquat[id * QUAT] = dcurquat[id * QUAT] = new_qs.x;
    dquat[id * QUAT + 1] = dcurquat[id * QUAT + 1] = new_qs.y;
    dquat[id * QUAT + 2] = dcurquat[id * QUAT + 2] = new_qs.z;
    dquat[id * QUAT + 3] = dcurquat[id * QUAT + 3] = new_qs.w;
}

//�L�сE����f����̃g���N�����߂�
//qs:���݁C���ڂ��Ă���G�b�W�̎p��
//pos1,pos2:�G�b�W�̗��[�̗��q�̈ʒu
//len:���
//kss:�L�э���
__device__ __host__
float4 StretchingShearTorque(float4 qs, float3 pos1, float3 pos2, float len, float kss) {
    float3 V = (pos1 - pos2) / len;

    float4 torqueSS;
    torqueSS.x = 4 * (qs.x * qs.x * qs.x + qs.x * qs.y * qs.y + qs.x * qs.z * qs.z + qs.x * qs.w * qs.w + V.z * qs.x - V.x * qs.z + V.y * qs.w);//x
    torqueSS.y = 4 * (qs.y * qs.y * qs.y + qs.y * qs.x * qs.x + qs.y * qs.z * qs.z + qs.y * qs.w * qs.w + V.z * qs.y - V.x * qs.w - V.y * qs.z);//y
    torqueSS.z = 4 * (qs.z * qs.z * qs.z + qs.z * qs.x * qs.x + qs.z * qs.y * qs.y + qs.z * qs.w * qs.w - V.z * qs.z - V.x * qs.x - V.y * qs.y);//z
    torqueSS.w = 4 * (qs.w * qs.w * qs.w + qs.w * qs.x * qs.x + qs.w * qs.y * qs.y + qs.w * qs.z * qs.z - V.z * qs.w - V.x * qs.y + V.y * qs.x);//w
    torqueSS = 1.f / 2.f * kss * torqueSS;

   /* float a = qs.x;
    float b = qs.y;
    float c = qs.z;
    float d = qs.w;

    torqueSS.x = 2 * (V.x - 2 * a * c - 2 * b * d) * (-2 * c) + 2 * (V.y + 2 * a * d - 2 * b * c) * (2 * d) + 2 * (V.z - d * d + a * a + b * b + c * c) * 2 * a;
    torqueSS.y = 2 * (V.x - 2 * a * c - 2 * b * d) * (-2 * d) + 2 * (V.y + 2 * a * d - 2 * b * c) * (-2 * c) + 2 * (V.z - d * d + a * a + b * b + c * c) * 2 * b;
    torqueSS.z = 2 * (V.x - 2 * a * c - 2 * b * d) * (-2 * a) + 2 * (V.y + 2 * a * d - 2 * b * c) * (-2 * b) + 2 * (V.z - d * d + a * a + b * b + c * c) * 2 * c;
    torqueSS.w = 2 * (V.x - 2 * a * c - 2 * b * d) * (-2 * b) + 2 * (V.y + 2 * a * d - 2 * b * c) * (2 * a) + 2 * (V.z - d * d + a * a + b * b + c * c) * (-2 * d);

    torqueSS *= 1.f / 2.f * kss;*/

    return torqueSS;
}

//�Ȃ��E�˂��ꐧ��̃g���N�����߂�(���ݗ��p����)
__device__ __host__
float4 BendTwistTorque(float4 q1, float4 q2, float4 Darboux, float kbt) {
    float4 torqueBT;
    float sum = q2.x * q2.x + q2.y * q2.y + q2.z * q2.z + q2.w * q2.w;
    torqueBT.x = sum * q1.x - q2.x * Darboux.w + q2.y * Darboux.z - q2.z * Darboux.y + q2.w * Darboux.x;
    torqueBT.y = sum * q1.y - q2.x * Darboux.z - q2.y * Darboux.w + q2.z * Darboux.x + q2.w * Darboux.y;
    torqueBT.z = sum * q1.z + q2.x * Darboux.y - q2.y * Darboux.x - q2.z * Darboux.w + q2.w * Darboux.z;
    torqueBT.w = sum * q1.w - q2.x * Darboux.x - q2.y * Darboux.y - q2.z * Darboux.z - q2.w * Darboux.w;

    return kbt * torqueBT;
}

//�g���N���狁�߂�ꍇ�ɗ��p(���ݗ��p����)
__device__ __host__
float4 SolveTorqueSolver(float4 q1, float4 q2, float4 torqueSS, float4 torqueBT, float kbt) {
    float4 rightForm = -(torqueSS + torqueBT);
    rightForm = rightForm / kbt;

    float4 New_Darboux;

    float sum = q2.x * q2.x + q2.y * q2.y + q2.z * q2.z + q2.w * q2.w;
    float b[4];
    b[0] = rightForm.x - sum * q1.x;
    b[1] = rightForm.y - sum * q1.y;
    b[2] = rightForm.z - sum * q1.z;
    b[3] = rightForm.w - sum * q1.w;

    float A[4][4];
    A[0][0] = q2.w;
    A[0][1] = -q2.z;
    A[0][2] = q2.y;
    A[0][3] = -q2.x;

    A[1][0] = q2.z;
    A[1][1] = q2.w;
    A[1][2] = -q2.x;
    A[1][3] = -q2.y;

    A[2][0] = -q2.y;
    A[2][1] = q2.x;
    A[2][2] = q2.w;
    A[2][3] = -q2.z;

    A[3][0] = -q2.x;
    A[3][1] = -q2.y;
    A[3][2] = -q2.z;
    A[3][3] = -q2.w;

    float x[4];
    int n1=LUDecomp(A, 4);
    //if (n1 == 0) printf("LUDecomp failure!!\n");
    int n2=LUSolver(A, b, x, 4);
    //if (n2 == 0)printf("LUSolver failure!!\n");

    New_Darboux = make_float4(x[0], x[1], x[2], x[3]);
    //printf("%d:New_Darboux x:%f,y:%f,z:%f,w:%\nf", New_Darboux.x, New_Darboux.y, New_Darboux.z, New_Darboux.w);
    return New_Darboux;
}

//�l�����̋t����Ԃ�
__device__ __host__
float4 QuatInverse(float4 quat) {
    return quatConjugate(quat) / dot(quat, quat);
}

//����̕����Q�l�ɐ��`�V�X�e���I�ɉ������C���ۂɂ́A��ԉ��̃G�b�W���珇�Ԃɉ����Ă���
//��ڂ̃G�b�W�͊��S�ɌŒ肷��̂ŁA��ڂ̃G�b�W�܂ł��l����
//�אڃG�b�W�̍X�V��̒l���g���K�v������̂ŁA�є��P�ʂł̕��񉻂ɂȂ�
//dpos:�ʒu
//dquat:�p��
//domega:��_���{�[�x�N�g��
//dlen:���
//dkss:�L�э���
//dkbt:�Ȃ�����
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//last_index:�є����Ƃ̍Ō�̗��q�̃C���f�b�N�X���i�[
//num_elastic:�����ł́C�є����Ƃɕ���v�Z���邽�߁C�є��̐���n��
__global__
void CxGlobalTorqueStep(float* dpos, float* dquat, float* domega, float* dlen, float* dkss, float* dkbt, int* dfix, int* dlast_index, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    if (id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�) �Ō�̗��q�̓X�L�b�v
    int last_edge = dlast_index[id] - 1;//�Ō�̗��q-1���Ō�̃G�b�W
    float min;
    if (id == 0) min = 0;
    else min = dlast_index[id - 1]+1;//�ŏ��̗��q

    //��ԉ��̃G�b�W(i==last_edge)�̏����͑��ƈႢ�A��݂̂̊�_���{�[�x�N�g�����ւ��B��x����Ȃ̂ŁA���[�v���g�킸�ɏ������C�����͈ȍ~�̃��[�v�Ƃقړ���
    float3 init_pos1 = make_float3(dpos[last_edge * DIM], dpos[last_edge * DIM + 1], dpos[last_edge * DIM + 2]);
    float3 init_pos2 = make_float3(dpos[(last_edge + 1) * DIM], dpos[(last_edge + 1) * DIM + 1], dpos[(last_edge + 1) * DIM + 2]);

    float4 init_quat0= make_float4(dquat[(last_edge - 1) * QUAT], dquat[(last_edge - 1) * QUAT + 1], dquat[(last_edge - 1) * QUAT + 2], dquat[(last_edge - 1) * QUAT + 3]);
    float4 init_quat1= make_float4(dquat[last_edge * QUAT], dquat[last_edge * QUAT + 1], dquat[last_edge * QUAT + 2], dquat[last_edge * QUAT + 3]);

    float init_l0 = dlen[last_edge];
    float init_kss = dkss[last_edge];
    float init_kbt = dkbt[last_edge];

    float init_K = init_kbt;
    float4 init_kq_inv = QuatInverse(init_K * init_quat1);

    float4 init_torqueSS = StretchingShearTorque(init_quat1, init_pos1, init_pos2, init_l0, init_kss);
    init_torqueSS = quatProduct(init_torqueSS, init_kq_inv);//�K���ɒ萔�������Ă݂�

    float4 init_Cur_Omega_Prev = quatProduct(quatConjugate(init_quat0), init_quat1);
    float4 init_Rest_Omega_Prev = init_Cur_Omega_Prev - init_torqueSS;
    
    //printf("id %d init_torqueSS x:%f,y:%f,z:%f,w:%f\n", id, init_torqueSS.x, init_torqueSS.y, init_torqueSS.z, init_torqueSS.w);
    //printf("id %d init_Cur_Omega_Prev x:%f,y:%f,z:%f,w:%f\n", id, init_Cur_Omega_Prev.x, init_Cur_Omega_Prev.y, init_Cur_Omega_Prev.z, init_Cur_Omega_Prev.w);
    //printf("id %d init_Rest_Omega_Prev x:%f,y:%f,z:%f,w:%f\n", id, init_Rest_Omega_Prev.x, init_Rest_Omega_Prev.y, init_Rest_Omega_Prev.z, init_Rest_Omega_Prev.w);

    //init_Rest_Omega_Prev = quatConjugate(init_Rest_Omega_Prev);
    //�_���{�[�x�N�g���̕������Ȃ��˂��ꐧ��̂悤�ɋ��߂�K�v�����邪�A�ق�1�ł���Ɛ����ł��邽�߁A1�ň���(s*Omega�����̂܂܃_���{�[�x�N�g���̔z��ɓ����)
    domega[(last_edge - 1) * QUAT] = init_Rest_Omega_Prev.x;
    domega[(last_edge - 1) * QUAT + 1] = init_Rest_Omega_Prev.y;
    domega[(last_edge - 1) * QUAT + 2] = init_Rest_Omega_Prev.z;
    domega[(last_edge - 1) * QUAT + 3] = init_Rest_Omega_Prev.w;
    //-----------------------------------------------------------------------------------------------------------------------------

    //���̃G�b�W�̏���
    for (int i = last_edge-1; i > min; i--) {//i�̂ЂƂO�̊�_���{�[�x�N�g�������߂�
        float3 pos1 = make_float3(dpos[i * DIM], dpos[i * DIM + 1], dpos[i * DIM + 2]);
        float3 pos2 = make_float3(dpos[(i + 1) * DIM], dpos[(i + 1) * DIM + 1], dpos[(i + 1) * DIM + 2]);

        float4 quat0 = make_float4(dquat[(i - 1) * QUAT], dquat[(i - 1) * QUAT + 1], dquat[(i - 1) * QUAT + 2], dquat[(i - 1) * QUAT + 3]);
        float4 quat1 = make_float4(dquat[i * QUAT], dquat[i * QUAT + 1], dquat[i * QUAT + 2], dquat[i * QUAT + 3]);
        float4 quat2 = make_float4(dquat[(i + 1) * QUAT], dquat[(i + 1) * QUAT + 1], dquat[(i + 1) * QUAT + 2], dquat[(i + 1) * QUAT + 3]);

        float l0 = dlen[i];
        float kss = dkss[i];
        float kbt = dkbt[i];

        //�L�сE����f����̃g���N������Ŋ���
        float K = kbt;
        float4 kq_inv = QuatInverse(K * quat1);

        //�L�сE����f����̃g���N�����߂�
        float4 torqueSS = StretchingShearTorque(quat1, pos1, pos2, l0, kss);
        torqueSS = quatProduct(torqueSS, kq_inv);//100�{

        float4 Cur_Omega_Prev = quatProduct(quatConjugate(quat0), quat1);//���݂̃G�b�W�ƂЂƂO�̃G�b�W�̃_���{�[�x�N�g��
        float4 Cur_Omega_Next = quatConjugate(quatProduct(quatConjugate(quat1), quat2));//���݂̃G�b�W�ƈ�ۂ̃G�b�W�̃_���{�[�x�N�g��
        //�����@���ŋ��߂�(�v�Z���ʂ͕ς��Ȃ���)
        //Cur_Omega_Next = quatProduct(quatConjugate(quat2), quat1);
        
        float4 Rest_Omega_Next = make_float4(domega[i * QUAT], domega[i * QUAT + 1], domega[i * QUAT + 2], domega[i * QUAT + 3]);
        //�ŏI�I�ɋ��߂��_���{�[�x�N�g��
        float4 Rest_Omega_Prev = Cur_Omega_Next + Cur_Omega_Prev - Rest_Omega_Next - torqueSS;//�ЂƂO�̃G�b�W�Ƃ̊Ԃ̊�_���{�[�x�N�g��(Appendix������torqueSS��+�ɕύX)
        Rest_Omega_Prev = quatConjugate(Rest_Omega_Prev);
        //���ʂ��o��
        /*if (i == min + 1) {
            printf("id %d torqueSS x:%f,y:%f,z:%f,w:%f\n", id, torqueSS.x, torqueSS.y, torqueSS.z, torqueSS.w);
            printf("id %d Cur_Omega_Prev x:%f,y:%f,z:%f,w:%f\n", id, Cur_Omega_Prev.x, Cur_Omega_Prev.y, Cur_Omega_Prev.z, Cur_Omega_Prev.w);
            printf("id %d Rest_Omega_Prev x:%f,y:%f,z:%f,w:%f\n", id, Rest_Omega_Prev.x, Rest_Omega_Prev.y, Rest_Omega_Prev.z, Rest_Omega_Prev.w);
        }*/

        //�_���{�[�x�N�g���̕������Ȃ��˂��ꐧ��̂悤�ɋ��߂�K�v�����邪�A�ق�1�ł���Ɛ����ł��邽�߁A1�ň���(s*Omega�����̂܂܃_���{�[�x�N�g���̔z��ɓ����)
        domega[(i - 1) * QUAT] = Rest_Omega_Prev.x;
        domega[(i - 1) * QUAT + 1] = Rest_Omega_Prev.y;
        domega[(i - 1) * QUAT + 2] = Rest_Omega_Prev.z;
        domega[(i - 1) * QUAT + 3] = Rest_Omega_Prev.w;
    }
}

//���܂łƂ͋t�ɏォ������悤�ɂ���
//
__global__
void CxGlobalTorqueStep_Upstair(float* dpos, float* dquat, float* domega, float* dlen, float* dkss, float* dkbt, int* dfix, int* dlast_index, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    if (id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�) �Ō�̗��q�̓X�L�b�v
    int last_edge = dlast_index[id] - 1;//�Ō�̗��q-1���Ō�̃G�b�W
    int min;
    if (id == 0) min = 1;
    else min = dlast_index[id - 1] + 2;//�ŏ��̗��q

    //��ԉ��̃G�b�W(i==last_edge)�̏����͑��ƈႢ�A��݂̂̊�_���{�[�x�N�g�����ւ��B��x����Ȃ̂ŁA���[�v���g�킸�ɏ������C�����͈ȍ~�̃��[�v�Ƃقړ���
    float3 init_pos1 = make_float3(dpos[min * DIM], dpos[min * DIM + 1], dpos[min * DIM + 2]);
    float3 init_pos2 = make_float3(dpos[(min + 1) * DIM], dpos[(min + 1) * DIM + 1], dpos[(min + 1) * DIM + 2]);

    float4 init_quat0 = make_float4(dquat[(min - 1) * QUAT], dquat[(min - 1) * QUAT + 1], dquat[(min - 1) * QUAT + 2], dquat[(min - 1) * QUAT + 3]);
    float4 init_quat1 = make_float4(dquat[min * QUAT], dquat[min * QUAT + 1], dquat[min * QUAT + 2], dquat[min * QUAT + 3]);

    float init_l0 = dlen[min];
    float init_kss = dkss[min];
    float init_kbt = dkbt[min];

    float init_K = init_kbt;
    float4 init_kq_inv = QuatInverse(init_K * init_quat1);

    float4 init_torqueSS = StretchingShearTorque(init_quat1, init_pos1, init_pos2, init_l0, init_kss);
    init_torqueSS = quatProduct(init_torqueSS, init_kq_inv);//�K���ɒ萔�������Ă݂�

    float4 init_Cur_Omega_Prev = quatProduct(quatConjugate(init_quat0), init_quat1);
    float4 init_Rest_Omega_Prev = init_Cur_Omega_Prev - init_torqueSS;

    //printf("id %d init_torqueSS x:%f,y:%f,z:%f,w:%f\n", id, init_torqueSS.x, init_torqueSS.y, init_torqueSS.z, init_torqueSS.w);
    //printf("id %d init_Cur_Omega_Prev x:%f,y:%f,z:%f,w:%f\n", id, init_Cur_Omega_Prev.x, init_Cur_Omega_Prev.y, init_Cur_Omega_Prev.z, init_Cur_Omega_Prev.w);
    //printf("id %d init_Rest_Omega_Prev x:%f,y:%f,z:%f,w:%f\n", id, init_Rest_Omega_Prev.x, init_Rest_Omega_Prev.y, init_Rest_Omega_Prev.z, init_Rest_Omega_Prev.w);

    //init_Rest_Omega_Prev = quatConjugate(init_Rest_Omega_Prev);
    //�_���{�[�x�N�g���̕������Ȃ��˂��ꐧ��̂悤�ɋ��߂�K�v�����邪�A�ق�1�ł���Ɛ����ł��邽�߁A1�ň���(s*Omega�����̂܂܃_���{�[�x�N�g���̔z��ɓ����)
    domega[min * QUAT] = init_Rest_Omega_Prev.x;
    domega[min * QUAT + 1] = init_Rest_Omega_Prev.y;
    domega[min * QUAT + 2] = init_Rest_Omega_Prev.z;
    domega[min * QUAT + 3] = init_Rest_Omega_Prev.w;
    //-----------------------------------------------------------------------------------------------------------------------------

    //���̃G�b�W�̏���
    for (int i = min+1; i < last_edge-1; i++) {//i�̂ЂƂO�̊�_���{�[�x�N�g�������߂�
        float3 pos1 = make_float3(dpos[i * DIM], dpos[i * DIM + 1], dpos[i * DIM + 2]);
        float3 pos2 = make_float3(dpos[(i + 1) * DIM], dpos[(i + 1) * DIM + 1], dpos[(i + 1) * DIM + 2]);

        float4 quat0 = make_float4(dquat[(i - 1) * QUAT], dquat[(i - 1) * QUAT + 1], dquat[(i - 1) * QUAT + 2], dquat[(i - 1) * QUAT + 3]);
        float4 quat1 = make_float4(dquat[i * QUAT], dquat[i * QUAT + 1], dquat[i * QUAT + 2], dquat[i * QUAT + 3]);
        float4 quat2 = make_float4(dquat[(i + 1) * QUAT], dquat[(i + 1) * QUAT + 1], dquat[(i + 1) * QUAT + 2], dquat[(i + 1) * QUAT + 3]);

        float l0 = dlen[i];
        float kss = dkss[i];
        float kbt = dkbt[i];

        //�L�сE����f����̃g���N������Ŋ���
        float K = kbt;
        float4 kq_inv = QuatInverse(K * quat1);

        //�L�сE����f����̃g���N�����߂�
        float4 torqueSS = StretchingShearTorque(quat1, pos1, pos2, l0, kss);
        torqueSS = quatProduct(torqueSS, kq_inv);

        float4 Cur_Omega_Prev = quatProduct(quatConjugate(quat0), quat1);//���݂̃G�b�W�ƂЂƂO�̃G�b�W�̃_���{�[�x�N�g��
        float4 Cur_Omega_Next = quatConjugate(quatProduct(quatConjugate(quat1), quat2));//���݂̃G�b�W�ƈ�ۂ̃G�b�W�̃_���{�[�x�N�g��
        //�����@���ŋ��߂�(�v�Z���ʂ͕ς��Ȃ���)
        //Cur_Omega_Next = quatProduct(quatConjugate(quat2), quat1);

        float4 Rest_Omega_Next = make_float4(domega[(i - 1) * QUAT], domega[(i - 1) * QUAT + 1], domega[(i - 1) * QUAT + 2], domega[(i - 1) * QUAT + 3]);
        //�ŏI�I�ɋ��߂��_���{�[�x�N�g��
        float4 Rest_Omega_Prev = Cur_Omega_Next + Cur_Omega_Prev - Rest_Omega_Next - torqueSS;//�ЂƂO�̃G�b�W�Ƃ̊Ԃ̊�_���{�[�x�N�g��(Appendix������torqueSS��+�ɕύX)
        Rest_Omega_Prev = quatConjugate(Rest_Omega_Prev);
        //���ʂ��o��
        /*if (i == min + 1) {
            printf("id %d torqueSS x:%f,y:%f,z:%f,w:%f\n", id, torqueSS.x, torqueSS.y, torqueSS.z, torqueSS.w);
            printf("id %d Cur_Omega_Prev x:%f,y:%f,z:%f,w:%f\n", id, Cur_Omega_Prev.x, Cur_Omega_Prev.y, Cur_Omega_Prev.z, Cur_Omega_Prev.w);
            printf("id %d Rest_Omega_Prev x:%f,y:%f,z:%f,w:%f\n", id, Rest_Omega_Prev.x, Rest_Omega_Prev.y, Rest_Omega_Prev.z, Rest_Omega_Prev.w);
        }*/

        //�_���{�[�x�N�g���̕������Ȃ��˂��ꐧ��̂悤�ɋ��߂�K�v�����邪�A�ق�1�ł���Ɛ����ł��邽�߁A1�ň���(s*Omega�����̂܂܃_���{�[�x�N�g���̔z��ɓ����)
        domega[i * QUAT] = Rest_Omega_Prev.x;
        domega[i * QUAT + 1] = Rest_Omega_Prev.y;
        domega[i * QUAT + 2] = Rest_Omega_Prev.z;
        domega[i * QUAT + 3] = Rest_Omega_Prev.w;
    }
}


//LocalTorqueStep�������̂Ɏg��
//cur_omega:���݂̓�̎p�����狁�߂�_���{�[�x�N�g��
//rest_omega:GlobalTorqueStep�ŋ��߂����K���O�̊�_���{�[�x�N�g��
//bendK:�Ȃ�����
//K_min:LocalTorqueStep�ł̒����p�����[�^
__device__ __host__
float4 solveInverseRot(float4 cur_omega, float4 rest_omega, float& bendK,float K_min) {
    const float SAFETY_FACTOR = min(abs(cur_omega.w), K_min);//0.00002f,0.002f,�ŏI�I�ɂ�0.005f
    //const float SAFETY_FACTOR = min(length(make_float3(cur_omega.x, cur_omega.y, cur_omega.z)), 0.2f);
    
    rest_omega -= dot(rest_omega, cur_omega) * cur_omega;

    //printf("omega_orth x:%f,y:%f,z:%f,w:%f\n", rest_omega.x, rest_omega.y, rest_omega.z, rest_omega.w);
    //printf("omega_orth length %f\n", length(rest_omega));

    float4 Omega = -rest_omega / bendK;

    //if (SAFETY_FACTOR > 0.2 + 1.0e-3 || SAFETY_FACTOR < 0.2 - 1.0e-3)printf("SAFETY_FACTOR %f\n", SAFETY_FACTOR);

    if (dot(Omega,Omega) > SAFETY_FACTOR * SAFETY_FACTOR) {//Length2��dot�ɒu������
        bendK = Length(rest_omega) / SAFETY_FACTOR;
        Omega = -rest_omega / bendK;
    }

    float4 d = sqrt(1 - dot(Omega, Omega)) * cur_omega;//Length2��dot�ɒu������
    Omega += d;

    return Omega;
}

//���[�J���g���N�X�e�b�v
//��_���{�[�x�N�g����K�؂Ȍ`�Ő��K������
//dquat:�p��
//domega:��_���{�[�x�N�g��
//deln:���
//dkbt:�Ȃ�����
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//K_min:LocalTorqueStep�ł̒����p�����[�^
//n:���q��(��_���{�[�x�N�g�����Ƃɕ���v�Z)
__global__
void CxLocalTorqueStep(float* dquat, float* domega, float* dlen,float* dkbt, int* dfix,float K_min,int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    if (id >= n - 2) return; // ���q���𒴂���X���b�hID�̃`�F�b�N
    if (dfix[id + 1] == 1 || dfix[id + 2] == 1) return;//��_���{�[�x�N�g�����Ȃ������̓X�L�b�v

    float4 quat1 = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);
    float4 quat2 = make_float4(dquat[QUAT * (id + 1)], dquat[QUAT * (id + 1) + 1], dquat[QUAT * (id + 1) + 2], dquat[QUAT * (id + 1) + 3]);
    float4 cur_omega = quatProduct(quatConjugate(quat1), quat2);//���݂̃_���{�[�x�N�g��

    float4 rest_omega = make_float4(domega[QUAT * id], domega[QUAT * id + 1], domega[QUAT * id + 2], domega[QUAT * id + 3]);

    float length = dlen[id];
    float tmp_kbt = dkbt[id];

    float4 last_omega = solveInverseRot(cur_omega, rest_omega, tmp_kbt, K_min);

    //�Ȃ��˂��ꐧ��̍����̍X�V
    dkbt[id] = tmp_kbt;
    //��_���{�[�x�N�g���̍X�V
    domega[QUAT * id] = last_omega.x;
    domega[QUAT * id + 1] = last_omega.y;
    domega[QUAT * id + 2] = last_omega.z;
    domega[QUAT * id + 3] = last_omega.w;
}

//���x����
//���x����ɗp����lambda�̌v�Z
//ddens:���݂̖��x
//drestdens:����x
//dpbf_lambda:���x����ɗp�����
//dvol:�̐�
//n:���q��
__global__
void CxPbfLambda(float* ddens,float* drestdens,float* dpbf_lambda,float* dvol,float* dmas,int n) {
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    float3 pos_i = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    //float m = params.mass;
    float a = params.aw;
    float rest_dens = params.rest_dens;
    //�C���f�b�N�X�̌v�Z
    uint sid = params.cell.dSortedIndex[id];
    //�C�V��SPH�ǉ�
    rest_dens = drestdens[sid];
    float dens = ddens[sid];

    // ���q�𒆐S�Ƃ��Ĕ��ah���Ɋ܂܂��O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos_i - make_float3(h));
    grid_pos1 = calcGridPos(pos_i + make_float3(h));
    //---------------------------------------------

    //���x���琧��̌v�Z�ɗp����ɂ����߂�
    float C = dens / rest_dens - 1.0f;
    if (C > 0.f) {
        float lambda_denom_i = 0.f;
        float3 grad_i_C = make_float3(0.f);
        for (int z = grid_pos0.z; z <= grid_pos1.z; ++z) {
            for (int y = grid_pos0.y; y <= grid_pos1.y; ++y) {
                for (int x = grid_pos0.x; x <= grid_pos1.x; ++x) {
                    int3 ngrid = make_int3(x, y, z);
                    uint ghash = calcGridHash(ngrid);   // �O���b�h�n�b�V���l

                    // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                    uint startIndex = params.cell.dCellStart[ghash];
                    if (startIndex != 0xffffffff) {	// �Z������łȂ����̃`�F�b�N
                        // �Z�����̃p�[�e�B�N���Ŕ���
                        uint endIndex = params.cell.dCellEnd[ghash];
                        for (uint j = startIndex; j < endIndex; ++j) {
                            uint sj = params.cell.dSortedIndex[j];
                            float3 pos_j = params.cell.dSortedPos[j];
                            float rest_dens_j = drestdens[sj];

                            //�L�����a���ɑ��݂��邩���m�F
                            float3 r_ij = pos_i - pos_j;
                            float r = length(r_ij);

                            if (r <= 1.0e-3f) continue;
                            if (r < h) {
                                float q = h - r;
                                //float m = rest_dens * dvol[sj];//�������x�Ƒ̐ς���v�Z
                                float m = dmas[sj];
                                //Grad_J_C�����߂�
                                float3 grad_j_C = -m / rest_dens_j * (-params.ag * q * q * r_ij / r);//params.ag�ɃJ�[�l���̌��z�萔(spiky�J�[�l��)���i�[����Ă���.m/\rho_rest*W(�J�[�l��)
                                lambda_denom_i += dot(grad_j_C, grad_j_C);
                                grad_i_C += grad_j_C;
                            }
                        }
                    }
                }
            }
        }
        lambda_denom_i += dot(grad_i_C, grad_i_C);
        //�ɂɒl���i�[(float�^)
        dpbf_lambda[sid] = -C / (lambda_denom_i + 1.0e-6f);//�Â�1.0e-6f�Œ�`
    }

    else {
        dpbf_lambda[sid] = 0.f;
    }
}

//���x����ɂ��ʒu�C��
//dpos:�ʒu
//drestdens:����x
//dpbf_lambda:���x����ɗp�����
//dvol:�̐�
//n:���q��
__global__
void CxPbfConstraint(float*dpos,float* drestdens,float* dpbf_lambda,float* dvol,float* dmas,int n) {
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    float3 pos_i = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    //�C���f�b�N�X�̌v�Z
    uint sid = params.cell.dSortedIndex[id];
    //pbf�̌v�Z�ɕK�v�ȃɂ������Ă���
    float pbf_lambda_i = dpbf_lambda[sid];
    //�������x
    float restdens_i = drestdens[id];

    // ���q�𒆐S�Ƃ��Ĕ��ah���Ɋ܂܂��O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos_i - make_float3(h));
    grid_pos1 = calcGridPos(pos_i + make_float3(h));

    // ���͂̃O���b�h�Z�����܂߂ċߖT�T�����Ė��x�v�Z
    float3 delta_pos_i = make_float3(0.f);
    for (int z = grid_pos0.z; z <= grid_pos1.z; ++z) {
        for (int y = grid_pos0.y; y <= grid_pos1.y; ++y) {
            for (int x = grid_pos0.x; x <= grid_pos1.x; ++x) {
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // �O���b�h�n�b�V���l
                // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                uint startIndex = params.cell.dCellStart[ghash];
                if (startIndex != 0xffffffff) {	// �Z������łȂ����̃`�F�b�N
                    // �Z�����̃p�[�e�B�N���Ŕ���
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for (uint j = startIndex; j < endIndex; ++j) {
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos_j = params.cell.dSortedPos[j];
                        float restdens_j = drestdens[sj];
                        float pbf_lambda_j = dpbf_lambda[sj];

                        float3 r_ij = pos_i - pos_j;
                        float r = length(r_ij);
                        if (r <= 1.0e-3) continue;
                        if (r < h) {
                            float q;
                            q = h - r;
                            float m = restdens_j * dvol[sj];//�������x�Ƒ̐ς��玿�ʂ��v�Z
                            //printf("pbf mass %f\n", m);
                            m = dmas[sid];
                            delta_pos_i += m / restdens_i * (pbf_lambda_i + pbf_lambda_j) * (params.ag * q * q * r_ij / r);
                        }
                    }
                }
            }
        }
    }

    //���ʂ̊i�[
    dpos[DIM * sid] += delta_pos_i.x;
    dpos[DIM * sid + 1] += delta_pos_i.y;
    dpos[DIM * sid + 2] += delta_pos_i.z;
}

//SPH��PBF�ŉ����ꍇ�̈��͌v�Z���������ꍇ
//dacc:SPH�ł̈ʒu���X�V����ۂɗp��������x
//datt:���q����(0�ŗ���,1�ŋ��E)
//power:���Ȃǂ̗�
//n:���q��
__global__
void CxPbfExternalForces(float* dacc, int* datt, float3 power, bool m_wind_flag, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return;

    uint sid = params.cell.dSortedIndex[id];
    if (datt[sid] != 0) {  // ���E���q�̏ꍇ�͗��q�ɂ������=0
        float3 v0 = make_float3(0.0f);
        dacc[DIM * sid + 0] = v0.x;  dacc[DIM * sid + 1] = v0.y; dacc[DIM * sid + 2] = v0.z;
        return;
    }

    float3 acc = make_float3(dacc[DIM * id], dacc[DIM * id + 1], dacc[DIM * id + 2]);
    acc = params.gravity;

    if (m_wind_flag) {
        acc += power;
    }

    dacc[DIM * id] = acc.x;
    dacc[DIM * id + 1] = acc.y;
    dacc[DIM * id + 2] = acc.z;
}

//���C����̎���
//���C����͂Ƃ肠�����C���x�X�V�̍ۂɁC��x�̂ݍs���悤�ɐݒ肷��D
//��������ꍇ�ɂ�XPBD�ɂ��Ȃ��ƁC�����Ɉˑ����ďC���ʂ��ω�����ƍl�����邪�C���C����͐������C�����킯�ł͂Ȃ��D
//���ӗ��q���厖�Ȃ̂ŁCSort���ꂽ�ʒu�𗘗p����
//dpos:�ʒu
//dcurpos:�ʒu�C���O�̈ʒu
//drestdens:����x
//dvol:�̐�(���ʂ����z���ʂ����_0*V�Œ�`)
//ddens:���݂̖��x
//dfix:dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//n:���q��
__global__
void CxFrictionConstraint(float* dpos, float* dcurpos,float* drestdens,float* dvol,float*ddens, int* dfix, int n) {
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    float3 pos_i = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    //�C���f�b�N�X�̌v�Z
    uint sid = params.cell.dSortedIndex[id];

    if (dfix[sid] == 1) return;//�Œ�_�ł���΁C�X�L�b�v

    //�O�X�e�b�v�̈ʒu
    float3 cur_pos_i = make_float3(dcurpos[sid * DIM], dcurpos[sid * DIM + 1], dcurpos[sid * DIM + 2]);
    //�ʒu�C���Ȃǂɂ��ړ���
    float3 v_i = pos_i - cur_pos_i;
    //�Î~���C�W��
    float mu = MU;

    // ���q�𒆐S�Ƃ��Ĕ��ah���Ɋ܂܂��O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos_i - make_float3(h));
    grid_pos1 = calcGridPos(pos_i + make_float3(h));

    //�ŏI�I�Ȗ��C��
    float3 x_fric = make_float3(0.f);

    for (int z = grid_pos0.z; z <= grid_pos1.z; ++z) {
        for (int y = grid_pos0.y; y <= grid_pos1.y; ++y) {
            for (int x = grid_pos0.x; x <= grid_pos1.x; ++x) {
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // �O���b�h�n�b�V���l
                // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                uint startIndex = params.cell.dCellStart[ghash];
                if (startIndex != 0xffffffff) {	// �Z������łȂ����̃`�F�b�N
                    // �Z�����̃p�[�e�B�N���Ŕ���
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for (uint j = startIndex; j < endIndex; ++j) {
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos_j = params.cell.dSortedPos[j];

                        //�O�X�e�b�v�̈ʒu
                        float3 cur_pos_j = make_float3(dcurpos[sj * DIM], dcurpos[sj * DIM + 1], dcurpos[sj * DIM + 2]);
                        //j�̗��q�̏������x
                        float restdens_j = drestdens[sj];
                        //j�̑̐�
                        float vol_j = dvol[sj];
                        //����
                        float m = restdens_j * vol_j;

                        float3 r_ij = pos_i - pos_j;
                        float r = length(r_ij);
                        if (r <= 1.0e-3) continue;
                        if (r < h) {
                            //�����ɖ��C���������
                            float3 v_j = pos_j - cur_pos_j;

                            float3 v_ij = v_i - v_j;

                            r_ij = normalize(r_ij);
                            //�Փ˖@���ɑ΂��Đ����Ȑ��������߂�
                            float3 delxn = v_ij - r_ij * dot(v_ij, r_ij);//delta x_��=v_ij-x_||

                            float q = h * h - r * r;//(h^2-||rij||^2)
                            x_fric += m / ddens[sj] * MU * delxn * params.aw * q * q * q;//aw*q^3
                        }
                    }
                }
            }
        }
    }

    //i�̈ړ��ʂɂ�����x_friction�����̐��������߂�
    float3 norm_x_fric = normalize(x_fric);
    float3 dir_i_fric = norm_x_fric * dot(v_i, norm_x_fric);

    //printf("id %d friction delta x:%f,y:%f,z:%f\n",id, x_fric.x, x_fric.y, x_fric.z);

    float3 delta_x;
    if (length(dir_i_fric) <= length(x_fric)) {//�Î~���C�͂Ƃ��Ĉ����p�^�[��
        delta_x = -dir_i_fric;
    }
    else {//�����C�͂Ƃ��Ĉ����p�^�[��
        //�����炭��肠��
        delta_x = -x_fric * min(MU / length(x_fric), 1.0f);
        //delta_x = -dir_i_fric * 0.1;
        //delta_x = x_fric * 0.1;
        //delta_x = make_float3(0.f);
    }

    //�ʒu�C���ɂ��X�V
    dpos[DIM * sid] += delta_x.x;
    dpos[DIM * sid + 1] += delta_x.y;
    dpos[DIM * sid + 2] += delta_x.z;
}

//���C����̎���
//�p����L�сE����f����̈ʒu�C���ʂɊ�Â��āC��x�ɂ���ĕω������Ă݂�
__global__
void CxFrictionConstraint_withQuat(float* dpos, float* dcurpos, float* drestdens, float* dvol, float* ddens,float* dquat,float* dlen, int* dfix, int n) {
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    if (id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)
    //if (dfix[id] == 1) return;//�Œ�_�ł���΁C�X�L�b�v

    float3 pos_i = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    //�C���f�b�N�X�̌v�Z
    uint sid = params.cell.dSortedIndex[id];
    //�O�X�e�b�v�̈ʒu
    float3 cur_pos_i = make_float3(dcurpos[sid * DIM], dcurpos[sid * DIM + 1], dcurpos[sid * DIM + 2]);

    if (dfix[sid] == 1) return;
    //�ЂƂO��quat----------------------
    float4 quat1 = make_float4(dquat[QUAT * (sid - 1)], dquat[QUAT * (sid - 1) + 1], dquat[QUAT * (sid - 1) + 2], dquat[QUAT * (sid - 1) + 3]);
    //����quat
    float4 quat2 = make_float4(dquat[QUAT * sid], dquat[QUAT * sid + 1], dquat[QUAT * sid + 2], dquat[QUAT * sid + 3]);
    float len = dlen[sid];
    //-----------------------------------

    //�ʒu�C���Ȃǂɂ��ړ���
    float3 v_i = pos_i - cur_pos_i;
    //�Î~���C�W��(�����C�W���͐Î~���C�W����0.1�{�Ƃ���)
    float mu = MU;

    // ���q�𒆐S�Ƃ��Ĕ��ah���Ɋ܂܂��O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos_i - make_float3(h));
    grid_pos1 = calcGridPos(pos_i + make_float3(h));

    //�ŏI�I�Ȗ��C��
    float3 x_fric = make_float3(0.f);

    for (int z = grid_pos0.z; z <= grid_pos1.z; ++z) {
        for (int y = grid_pos0.y; y <= grid_pos1.y; ++y) {
            for (int x = grid_pos0.x; x <= grid_pos1.x; ++x) {
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // �O���b�h�n�b�V���l
                // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                uint startIndex = params.cell.dCellStart[ghash];
                if (startIndex != 0xffffffff) {	// �Z������łȂ����̃`�F�b�N
                    // �Z�����̃p�[�e�B�N���Ŕ���
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for (uint j = startIndex; j < endIndex; ++j) {
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos_j = params.cell.dSortedPos[j];

                        //�O�X�e�b�v�̈ʒu
                        float3 cur_pos_j = make_float3(dcurpos[sj * DIM], dcurpos[sj * DIM + 1], dcurpos[sj * DIM + 2]);
                        //j�̗��q�̏������x
                        float restdens_j = drestdens[sj];
                        //j�̑̐�
                        float vol_j = dvol[sj];
                        //����
                        float m = restdens_j * vol_j;

                        float3 r_ij = pos_i - pos_j;
                        float r = length(r_ij);
                        if (r <= 1.0e-3) continue;
                        if (r < h) {
                            //�����ɖ��C���������
                            float3 v_j = pos_j - cur_pos_j;

                            float3 v_ij = v_i - v_j;

                            r_ij = normalize(r_ij);
                            //�Փ˖@���ɑ΂��Đ����Ȑ��������߂�
                            float3 delxn = v_ij - r_ij * dot(v_ij, r_ij);//delta x_��=v_ij-x_||

                            float q = h * h - r * r;//(h^2-||rij||^2)
                            x_fric += m / ddens[sj] * MU * delxn * params.aw * q * q * q;//aw*q^3
                        }
                    }
                }
            }
        }
    }

    //i�̈ړ��ʂɂ�����x_friction�����̐��������߂�
    float3 norm_x_fric = normalize(x_fric);
    float3 dir_i_fric = norm_x_fric * dot(v_i, norm_x_fric);

    //printf("id %d friction delta x:%f,y:%f,z:%f\n",id, x_fric.x, x_fric.y, x_fric.z);

    float3 delta_x;
    if (length(dir_i_fric) <= length(x_fric)) {//�Î~���C�͂Ƃ��Ĉ����p�^�[��
        delta_x = -dir_i_fric;
    }
    else {//�����C�͂Ƃ��Ĉ����p�^�[��
        //delta_x = -dir_i_fric * min(length(x_fric) / length(dir_i_fric), 1.f);//[Macklin 2014]���Q�l�ɓK���ɋ��߂�
        delta_x = -x_fric * min(MU / length(x_fric), 1.0f);
        //delta_x = make_float3(0.f);
    }

    //�ʒu�C���ɂ��X�V
    dpos[DIM * sid] += delta_x.x;
    dpos[DIM * sid + 1] += delta_x.y;
    dpos[DIM * sid + 2] += delta_x.z;

    //�L�сE����f����̃�x=����/l0���C�������p���ɊҌ�����D��̃G�b�W�Ɖ��̃G�b�W�̗����ɊҌ�
    //���̂��߂ɁC�{��sid����Ƌ����ɕ�����K�v���邪�C�ǂݏo�����珑�����݂܂łɏ������������߁C���񉻂��Ă����v�ł��낤�Ɛ����D

    //�L�сE����f����ł̃���
    float3 lambda = delta_x * len;

    float4 delta_quat1 = 2.f * quatProduct(make_float4(lambda, 0.f), quatProduct(quat1, make_float4(0.f, 0.f, -1.f, 0.f)));
    float4 delta_quat2 = -2.f * quatProduct(make_float4(lambda, 0.f), quatProduct(quat2, make_float4(0.f, 0.f, -1.f, 0.f)));

    float4 new_quat1 = normalize(quat1 + delta_quat1);
    float4 new_quat2 = normalize(quat2 + delta_quat2);

    if (dfix[sid - 1] == 1) {
        dquat[QUAT * (sid - 1)] = new_quat1.x;
        dquat[QUAT * (sid - 1) + 1] = new_quat1.y;
        dquat[QUAT * (sid - 1) + 2] = new_quat1.z;
        dquat[QUAT * (sid - 1) + 3] = new_quat1.w;
    }
    if (dfix[sid + 1] == 1) {
        dquat[QUAT * sid] = new_quat2.x;
        dquat[QUAT * sid + 1] = new_quat2.y;
        dquat[QUAT * sid + 2] = new_quat2.z;
        dquat[QUAT * sid + 3] = new_quat2.w;
    }
}

//�X�̗��q�Ƃ̖��C���l���C�X�ɐÎ~���C���l����
//���C����͂Ƃ肠�����C���x�X�V�̍ۂɁC��x�̂ݍs���悤�ɐݒ肷��D
//��������ꍇ�ɂ�XPBD�ɂ��Ȃ��ƁC�����Ɉˑ����ďC���ʂ��ω�����ƍl�����邪�C���C����͐������C�����킯�ł͂Ȃ��D
//���ӗ��q���厖�Ȃ̂ŁCSort���ꂽ�ʒu�𗘗p����
//dpos:�ʒu
//dcurpos:�ʒu�C���O�̈ʒu
//drestdens:����x
//dvol:�̐�(���ʂ����z���ʂ����_0*V�Œ�`)
//ddens:���݂̖��x
//dfix:dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//n:���q��
__global__
void CxFrictionAllParticlesConstraint(float* dpos, float* dcurpos, float* drestdens, float* dvol, float* ddens, int* dfix, int n) {
    // �O���b�h,�u���b�N���̃X���b�h�ʒu�𗱎q�C���f�b�N�X�Ƃ���
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if (id >= n) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�)

    float3 pos_i = params.cell.dSortedPos[id];
    float h = params.effective_radius;
    //�C���f�b�N�X�̌v�Z
    uint sid = params.cell.dSortedIndex[id];

    if (dfix[sid] == 1) return;//�Œ�_�ł���΁C�X�L�b�v

    //�O�X�e�b�v�̈ʒu
    float3 cur_pos_i = make_float3(dcurpos[sid * DIM], dcurpos[sid * DIM + 1], dcurpos[sid * DIM + 2]);
    //�ʒu�C���Ȃǂɂ��ړ���
    float3 v_i = pos_i - cur_pos_i;
    //�Î~���C�W��
    float mu = MU;

    // ���q�𒆐S�Ƃ��Ĕ��ah���Ɋ܂܂��O���b�h(caclGridPos���ŋ��E��������)
    int3 grid_pos0, grid_pos1;
    grid_pos0 = calcGridPos(pos_i - make_float3(h));
    grid_pos1 = calcGridPos(pos_i + make_float3(h));

    //�ŏI�I�Ȗ��C��
    float3 x_fric = make_float3(0.f);

    for (int z = grid_pos0.z; z <= grid_pos1.z; ++z) {
        for (int y = grid_pos0.y; y <= grid_pos1.y; ++y) {
            for (int x = grid_pos0.x; x <= grid_pos1.x; ++x) {
                int3 ngrid = make_int3(x, y, z);
                uint ghash = calcGridHash(ngrid);   // �O���b�h�n�b�V���l
                // �Z�����̃p�[�e�B�N���̃X�^�[�g�C���f�b�N�X
                uint startIndex = params.cell.dCellStart[ghash];
                if (startIndex != 0xffffffff) {	// �Z������łȂ����̃`�F�b�N
                    // �Z�����̃p�[�e�B�N���Ŕ���
                    uint endIndex = params.cell.dCellEnd[ghash];
                    for (uint j = startIndex; j < endIndex; ++j) {
                        uint sj = params.cell.dSortedIndex[j];
                        float3 pos_j = params.cell.dSortedPos[j];

                        //�O�X�e�b�v�̈ʒu
                        float3 cur_pos_j = make_float3(dcurpos[sj * DIM], dcurpos[sj * DIM + 1], dcurpos[sj * DIM + 2]);
                        //j�̗��q�̏������x
                        float restdens_j = drestdens[sj];
                        //j�̑̐�
                        float vol_j = dvol[sj];
                        //����
                        float m = restdens_j * vol_j;

                        float3 r_ij = pos_i - pos_j;
                        float r = length(r_ij);
                        if (r <= 1.0e-3) continue;
                        if (r < h) {
                            //�����ɖ��C���������
                            float3 v_j = pos_j - cur_pos_j;

                            float3 v_ij = v_i - v_j;

                            r_ij = normalize(r_ij);

                            //�Փ˖@���ɑ΂��Đ����Ȑ��������߂�
                            float3 delxn = v_ij - r_ij * dot(v_ij, r_ij);//delta x_��=v_ij-x_||

                            float q = h * h - r * r;//(h^2-||rij||^2)
                            float3 tmp_x_fric = m / ddens[sj] * MU * delxn * params.aw * q * q * q;//aw*q^3

                            x_fric -= tmp_x_fric;

                            //�Î~���C�̎c�蕨
                            ////���C�͂𐳋K�����ĕ����x�N�g����
                            //float3 norm_tmp_x_fric = normalize(tmp_x_fric);
                            ////v_i�̂����C���C�͂̕����̐��������o��
                            //float3 dir_i_fric = norm_tmp_x_fric * dot(v_i, norm_tmp_x_fric);

                            ////�Î~���C�͂Ȃ炻����̕����̐�����ł�����
                            //if (length(dir_i_fric) <= length(tmp_x_fric)) {
                            //    x_fric -= dir_i_fric;
                            //}
                            ////�����C�Ȃ�C���̂܂ܓK�p���邱�ƂƂ���
                            //else {
                            //    x_fric -= tmp_x_fric;
                            //}
                        }
                    }
                }
            }
        }
    }

    //�ʒu�C���ɂ��X�V
    dpos[DIM * sid] += x_fric.x;
    dpos[DIM * sid + 1] += x_fric.y;
    dpos[DIM * sid + 2] += x_fric.z;
}

//�V����2���_����Ԃ̃G�b�W�̎p�������߂�
//dpos:�ʒu
//dquat:�p��
//dfix:�Œ�_(�є��̊J�n�_)�������z��(1�Ȃ�Œ�_,0�Ȃ炻��ȊO)
//n:���q��
__global__
void CxQuatSet(float* dpos, float* dquat, int* dfix, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    if (id >= n - 1) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�) �Ō�̗��q�̓X�L�b�v
    if (dfix[id] == 1 || dfix[id + 1] == 1) return;//�e���̖��ɕӂ̐������s��

    float3 pos0 = make_float3(dpos[DIM * id], dpos[DIM * id + 1], dpos[DIM * id + 2]);
    float3 pos1 = make_float3(dpos[DIM * (id + 1)], dpos[DIM * (id + 1) + 1], dpos[DIM * (id + 1) + 2]);
    //2���_�̊Ԃ̕����x�N�g��
    float3 dir = pos1 - pos0;
    //�G�b�W�̎p�������߂�
    float4 quat = quatFromDirector(dir);

    dquat[QUAT * id] = quat.x;
    dquat[QUAT * id + 1] = quat.y;
    dquat[QUAT * id + 2] = quat.z;
    dquat[QUAT * id + 3] = quat.w;
}

//���͂ɂ��g���N�̌v�Z������
__global__
void CxCalcTorque(float* dpos,float* dmas,float* dquat, float* dfss,float* dlength,float* dkss, int* dfix,float3 gravity, int n) {
    int id = blockDim.x * blockIdx.x + threadIdx.x;

    if (id >= n - 1) return; // ���q���𒴂���X���b�hID�̃`�F�b�N(�]�肪�o�Ȃ��悤�Ƀu���b�N���Ȃǂ��ݒ�ł���Ȃ�K�v�Ȃ�) �Ō�̗��q�̓X�L�b�v
    if (dfix[id] == 1 || dfix[id + 1] == 1) return;//�ŏ��̓�̗��q�͌Œ�_�Ƃ��Ĉ���

    float3 pos0 = make_float3(dpos[DIM * (id - 1)], dpos[DIM * (id - 1) + 1], dpos[DIM * (id - 1) + 2]);
    float3 pos1 = make_float3(dpos[DIM * id], dpos[DIM * id + 1], dpos[DIM * id + 2]);
    float3 mid = (pos0 + pos1) / 2.f;

    float4 quat = make_float4(dquat[QUAT * id], dquat[QUAT * id + 1], dquat[QUAT * id + 2], dquat[QUAT * id + 3]);
    float3 fss = make_float3(dfss[DIM * id], dfss[DIM * id + 1], dfss[DIM * id + 2]);
    //�p����e3����0->1�̕����x�N�g�������߁C����͒P�ʃx�N�g���Ȃ̂ŁC���������߂�K�v������D
    float len = dlength[id];

    float mass0 = dmas[id-1];
    float mass1 = dmas[id];

    float3 torque_calc = CalcNormalTorque(pos0, pos1, quat, fss, len, mass0, mass1, make_float3(0.f, -9.81, 0.f));
    printf("id %d calcTorque x:%f,y:%f,z:%f\n", id, torque_calc.x, torque_calc.y, torque_calc.z);

    float4 torque = StretchingShearTorque(quat, pos0, pos1, len, dkss[id]);
    //printf("id %d StretchingShearTorque x:%f,y:%f,z:%f,w:%f\n",id, torque.x, torque.y, torque.z, torque.w);
}
