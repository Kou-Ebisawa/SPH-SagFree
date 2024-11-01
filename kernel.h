/*!
  @file kernel.h
	
  @brief SPH法のカーネル計算
 
  @author Makoto Fujisawa
  @date 2012-12
*/

#ifndef _SPH_KERNEL_H_
#define _SPH_KERNEL_H_


//-----------------------------------------------------------------------------
// インクルードファイル
//-----------------------------------------------------------------------------
#include "utils.h"

//-----------------------------------------------------------------------------
// Poly6カーネル
//-----------------------------------------------------------------------------
/*!
 * カーネル係数
 * @param[in] h 有効半径
 * @param[in] d 次元(1,2,3)
 * @param[in] type ノーマル:1，勾配:2，ラプラシアン:3
 * @return カーネル係数値
 */
inline static float KernelCoefPoly6(float h, int d, int type)
{
	float a = 1.0;
	if(d < 1) d = 1;
	if(d > 3) d = 3;
	switch(type){
	case 1:	// ノーマル
		switch(d){
		case 2: a = 4.0/(RX_PI*pow((float)h, (float)8.0));			break;
		case 3:	a = 315.0/(64.0*RX_PI*pow((float)h, (float)9.0));		break;
		}
		break;

	case 2:	// 勾配
		switch(d){
		case 2:	a = -24.0/(RX_PI*pow((float)h, (float)8.0));			break;
		case 3:	a = -945.0/(32.0*RX_PI*pow((float)h, (float)9.0));	break;
		}
		break;

	case 3:	// ラプラシアン
		switch(d){
		case 2: a = -24.0/(RX_PI*pow((float)h, (float)8.0));			break;
		case 3: a = -945.0/(32.0*RX_PI*pow((float)h, (float)9.0));	break;
		}
		break;

	default:
		break;
	}

	return a;
}


/*!
 * Poly6カーネル関数値の計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @param[in] a カーネル係数
 * @return 関数値
 */
inline float KernelPoly6(float r, float h, float a)
{
	if(r >= 0.0 && r <= h){
		float q = h*h-r*r;
		return a*q*q*q;
	}
	else{
		return 0.0;
	}
}

/*!
 * Poly6カーネル関数勾配値の計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @param[in] rij 相対位置ベクトル
 * @param[in] a カーネル係数
 * @return 勾配値
 */
template<class T> 
inline T KernelPoly6G(float r, float h, float a, T rij)
{
	if(r >= 0.0 && r <= h){
		float q = h*h-r*r;
		return  a*q*q*rij;
	}
	else{
		return T(0.0);
	}
}
 
/*!
 * Splineカーネル関数ラプラシアンの計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @param[in] a カーネル係数
 * @param[in] d 次元(1,2,3)
 * @return ラプラシアンの値
 */
inline float KernelPoly6L(float r, float h, float a, float d)
{
	if(r >= 0.0 && r <= h){
		float q = h*h-r*r;
		return a*(3.0*q*q-4.0*r*r*q);
	}
	else{
		return 0.0;
	}
}


/*!
 * Poly6カーネルの値を計算
 * @param[in] r カーネル中心までの距離
 * @param[in] h 有効半径
 * @param[in] ptr 関数呼び出しポインタ
 * @return カーネル関数値
 */
static inline float CalKernelPoly6(float r, float h, void* ptr = 0)
{
	static float a = 0.0;
	if(a == 0.0) a = KernelCoefPoly6(h, 3, 1);
	return KernelPoly6(r, h, a);
}

/*!
 * Poly6カーネルの値を計算(最大値が1になるように正規化)
 * @param[in] r カーネル中心までの距離
 * @param[in] h 有効半径
 * @param[in] ptr 関数呼び出しポインタ
 * @return カーネル関数値
 */
static inline float CalKernelPoly6r(float r, float h, void* ptr = 0)
{
	static float a = 0.0;
	static float b = 0.0;
	if(a == 0.0) a = KernelCoefPoly6(h, 3, 1);
	if(b == 0.0) b = KernelPoly6(0.0, h, a);
	return KernelPoly6(r, h, a)/b;
}




//-----------------------------------------------------------------------------
// Spikyカーネル
//-----------------------------------------------------------------------------
/*!
 * カーネル係数
 * @param[in] h 有効半径
 * @param[in] d 次元(1,2,3)
 * @param[in] type ノーマル:1，勾配:2，ラプラシアン:3
 * @return カーネル係数値
 */
inline static float KernelCoefSpiky(float h, int d, int type)
{
	float a = 1.0;
	if(d < 1) d = 1;
	if(d > 3) d = 3;
	switch(type){
	case 1:	// ノーマル
		switch(d){
		case 2: a = 10.0/(RX_PI*pow((float)h, (float)5.0));	break;
		case 3:	a = 15.0/(RX_PI*pow((float)h, (float)6.0));	break;
		}
		break;

	case 2:	// 勾配
		switch(d){
		case 2:	a = -30.0/(RX_PI*pow((float)h, (float)5.0));	break;
		case 3:	a = -45.0/(RX_PI*pow((float)h, (float)6.0));	break;
		}
		break;

	case 3:	// ラプラシアン
		switch(d){
		case 2: a = -60.0/(RX_PI*pow((float)h, (float)5.0));	break;
		case 3: a = -90.0/(RX_PI*pow((float)h, (float)6.0));	break;
		}
		break;

	default:
		break;
	}

	return a;
}

/*!
 * Spikyカーネル関数値の計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @param[in] a カーネル係数
 * @return 関数値
 */
inline float KernelSpiky(float r, float h, float a)
{
	if(r >= 0.0 && r <= h){
		float q = h-r;
		return a*q*q*q;
	}
	else{
		return 0.0;
	}
}

/*!
 * Spikyカーネル関数勾配値の計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @param[in] a カーネル係数
 * @param[in] rij 相対位置ベクトル
 * @return 勾配値
 */
template<class T> 
inline T KernelSpikyG(float r, float h, float a, T rij)
{
	if(r > 0.0 && r <= h){
		float q = h-r;
		return  a*q*q*rij/r;
	}
	else{
		return T(0.0);
	}
}
 
/*!
 * Splineカーネル関数ラプラシアンの計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @param[in] a カーネル係数
 * @param[in] d 次元(1,2,3)
 * @return ラプラシアンの値
 */
inline float KernelSpikyL(float r, float h, float a, float d)
{
	if(r > 0.0 && r <= h){
		float q = h-r;
		return a*(q*q/r-q);
	}
	else{
		return 0.0;
	}
}



//-----------------------------------------------------------------------------
// Viscカーネル
//-----------------------------------------------------------------------------
/*!
 * カーネル係数
 * @param[in] h 有効半径
 * @param[in] d 次元(1,2,3)
 * @param[in] type ノーマル:1，勾配:2，ラプラシアン:3
 * @return カーネル係数値
 */
inline static float KernelCoefVisc(float h, int d, int type)
{
	float a = 1.0;
	if(d < 1) d = 1;
	if(d > 3) d = 3;
	switch(type){
	case 1:	// ノーマル
		switch(d){
		case 2: a = 10.0/(3.0*RX_PI*pow((float)h, (float)2.0));	break;
		case 3:	a = 15.0/(2.0*RX_PI*pow((float)h, (float)3.0));	break;
		}
		break;

	case 2:	// 勾配
		switch(d){
		case 2:	a = 10.0/(3.0*RX_PI*pow((float)h, (float)4.0));	break;
		case 3:	a = 15.0/(2.0*RX_PI*pow((float)h, (float)5.0));	break;
		}
		break;

	case 3:	// ラプラシアン
		switch(d){
		case 2: a = 20.0/(3.0*RX_PI*pow((float)h, (float)5.0));	break;
		case 3: a = 45.0/(RX_PI*pow((float)h, (float)6.0));		break;
		}
		break;

	default:
		break;
	}

	return a;
}

/*!
 * Viscカーネル関数値の計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @param[in] a カーネル係数
 * @return 関数値
 */
inline float KernelVisc(float r, float h, float a)
{
	if(r > 0.0 && r <= h){
		float q = r/h;
		return a*(-q*q*q/2.0+q*q+2.0/q-1.0);
	}
	else{
		return 0.0;
	}
}

/*!
 * Viscカーネル関数勾配値の計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @param[in] a カーネル係数
 * @param[in] rij 相対位置ベクトル
 * @return 勾配値
 */
template<class T> 
inline T KernelViscG(float r, float h, float a, T rij)
{
	if(r > 0.0 && r <= h){
		float q = r/h;
		return  a*(-1.5/q+2.0-q*q*q/2.0)*rij;
	}
	else{
		return T(0.0);
	}
}
 
/*!
 * Splineカーネル関数ラプラシアンの計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @param[in] a カーネル係数
 * @param[in] d 次元(1,2,3)
 * @return ラプラシアンの値
 */
inline float KernelViscL(float r, float h, float a, float d)
{
	if(r > 0.0 && r <= h){
		return a*(h-r);
	}
	else{
		return 0.0;
	}
}


//-----------------------------------------------------------------------------
// Splineカーネル
//-----------------------------------------------------------------------------
/*!
 * カーネル係数
 * @param[in] h 有効半径
 * @param[in] d 次元(1,2,3)
 * @param[in] type ノーマル:1，勾配:2，ラプラシアン:3
 * @return カーネル係数値
 */
inline static float KernelCoefSpline(float h, int d, int type)
{
	float a = 1.0;
	if(d < 1) d = 1;
	if(d > 3) d = 3;
	switch(type){
	case 1:	// ノーマル
		switch(d){
		case 1: a = 2.0/(3.0*h);				break;
		case 2: a = 10.0/(7.0*RX_PI*h*h);		break;
		case 3:	a = 1.0/(RX_PI*h*h*h);			break;
		}
		break;

	case 2:	// 勾配
		switch(d){
		case 1:	a = 3.0/(2.0*h*h*h);			break;
		case 2:	a = 45.0/(14.0*RX_PI*h*h*h*h);	break;
		case 3:	a = 9.0/(4.0*RX_PI*h*h*h*h*h);	break;
		}
		break;

	case 3:	// ラプラシアン
		switch(d){
		case 1: a = 1.0/(2.0*RX_PI*h*h*h);		break;
		case 2: a = 45.0/(42.0*RX_PI*h*h*h*h);	break;
		case 3: a = 3.0/(4.0*RX_PI*h*h*h*h*h);	break;
		}
		break;

	default:
		break;
	}

	return a;
}


/*!
 * Splineカーネル関数値の計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @param[in] a カーネル係数
 * @return 関数値
 */
inline float KernelSpline(float r, float h, float a)
{
	float q = r/h;
	if(q >= 0.0 && q <= 1.0){
		return a*(1.0-1.5*q*q+0.75*q*q*q);
	}
	else if(q > 1.0 && q <= 2.0){
		return a*0.25*(2.0-q)*(2.0-q)*(2.0-q);
	}
	else{
		return 0.0;
	}
}

/*!
 * Splineカーネル関数勾配値の計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @param[in] rij 相対位置ベクトル
 * @param[in] a カーネル係数
 * @return 勾配値
 */
template<class T> 
inline T KernelSplineG(float r, float h, float a, T rij)
{
	float q = r/h;
	if(q >= 0.0 && q < 1.0){
		return  a*(q-4.0/3.0)*rij;
	}
	else if(q >= 1.0 && q < 2.0){
		return -a*(2.0-q)*(2.0-q)*rij/q/3.0;
	}
	else{
		return T(0.0);
	}
}
 
/*!
 * Splineカーネル関数ラプラシアンの計算
 * @param[in] r 距離
 * @param[in] h 有効半径
 * @param[in] a カーネル係数
 * @param[in] d 次元(1,2,3)
 * @return ラプラシアンの値
 */
inline float KernelSplineL(float r, float h, float a, float d)
{
	float q = r/h;
	if(q >= 0.0 && q < 1.0){
		return a*(3.0*(d+1.0)*q-4.0*d);
	}
	else if(q >= 1.0 && q < 2.0){
		return a*((1.0-d)*(2.0-q)*(2.0-q)/q+2.0*(2.0-q));
	}
	else{
		return 0.0;
	}
}



#endif	// _RX_KERNEL_H_