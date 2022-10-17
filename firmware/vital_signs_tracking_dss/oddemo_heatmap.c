/**
 * oddemo_heatmap.c
 *
 * Functions relating to the Occupancy Detection Heatmap
 *
 * Copyright (C) 2019 Texas Instruments Incorporated - http://www.ti.com/ 
 * 
 * 
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions 
 *  are met:
 *
 *    Redistributions of source code must retain the above copyright 
 *    notice, this list of conditions and the following disclaimer.
 *
 *    Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the 
 *    documentation and/or other materials provided with the   
 *    distribution.
 *
 *    Neither the name of Texas Instruments Incorporated nor the names of
 *    its contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
 *  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
 *  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 *  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 *  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
 *  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <math.h>

/* C674x mathlib */
#include <ti/mathlib/mathlib.h>
#include <ti/alg/mmwavelib/mmwavelib.h>

#include "radar_c674x.h"

#include "../common/oddemo_common.h"
#include <float.h>
#include <oddemo_heatmap.h>

//global variables
uint16_t oddemo_row_init = 0;
float oddemo_heatmap_max;
float oddemo_row_noise[ODDEMO_MAX_RANGE];
ODDemo_Heatmap_Stats oddemo_stats[ODDEMO_MAX_RANGE];
uint8_t peak_positions[9];

extern ODDEMO_Parms oddemo_parms;

/**
 *  \fn     void MATRIX_4x4_BWInversionfp(
 *                            IN      cplxf_t  * RESTRICT Input,
 *                            OUT     cplxf_t  * RESTRICT output);
 *
 *  \brief   4x4 matrix inversion using block wise method.
 *
 *  \param[in]    Input
 *              Input 4x4 matrix that needs to be inversed, stored sequentially
 *              in format (1,1), (1,2).... (1,4), (2,1), (2,2)...(4,4).
 *
 *  \param[out]   output
 *              Output 4x4 matrix. Stored sequentially as the input.
 *
 *  \pre
 *
 *  \post
 *
 *  \sa
 *
 */
void MATRIX_4x4_BWInversionfp(cplxf_t *Input, cplxf_t *output)
{

    float A0, A3, D0, D3, F0, F3;
    float detA, oneOverdetA;
    float ftemp1;
    __float2_t A1, D1, F1, C0, C1, C2, C3;
    __float2_t B0, B1, B2, B3;
    __float2_t dtemp, dtemp1, dtemp2;

    /* load A */
    dtemp = _amem8_f2(&Input[0]);
    A0 = _hif2(dtemp);
    dtemp = _amem8_f2(&Input[5]);
    A3 = _hif2(dtemp);
    A1 = _amem8_f2(&Input[1]);

    /* Calculate D = inv(D) */
    dtemp = _amem8_f2(&Input[10]);
    D0 = _hif2(dtemp);
    dtemp = _amem8_f2(&Input[15]);
    D3 = _hif2(dtemp);
    D1 = _amem8_f2(&Input[11]);
    dtemp = _dmpysp(D1, D1);
    detA = D0 * D3 - _hif2(dtemp) - _lof2(dtemp);
    oneOverdetA = _rcpsp(detA);
    oneOverdetA = oneOverdetA * (2.f - detA * oneOverdetA);
    oneOverdetA = oneOverdetA * (2.f - detA * oneOverdetA);

    ftemp1 = D0 * oneOverdetA;
    D0 = D3 * oneOverdetA;
    D3 = ftemp1;
    D1 = _dmpysp(D1, _ftof2(-oneOverdetA, -oneOverdetA));

    /* load B */
    B0 = _amem8_f2(&Input[2]);
    B1 = _amem8_f2(&Input[3]);
    B2 = _amem8_f2(&Input[6]);
    B3 = _amem8_f2(&Input[7]);

    /* calculate C = B*inv(D) */
    //results           =   _cmpysp(D1, B1);
    //C0                =   _dsubsp(_hif2_128(results), _lof2_128(results));
    C0 = _complex_conjugate_mpysp(D1, B1);
    C0 = _daddsp(C0, _dmpysp(B0, _ftof2(D0, D0)));
    //results           =   _cmpysp(D1, B0);
    //C1                =   _daddsp(_hif2_128(results), _lof2_128(results));
    C1 = _complex_mpysp(D1, B0);
    dtemp1 = C1;
    C1 = _daddsp(C1, _dmpysp(B1, _ftof2(D3, D3)));
    //results           =   _cmpysp(D1, B3);
    //C2                =   _dsubsp(_hif2_128(results), _lof2_128(results));
    C2 = _complex_conjugate_mpysp(D1, B3);
    C2 = _daddsp(C2, _dmpysp(B2, _ftof2(D0, D0)));
    //results           =   _cmpysp(D1, B2);
    //C3                =   _daddsp(_hif2_128(results), _lof2_128(results));
    C3 = _complex_mpysp(D1, B2);
    dtemp2 = C3;
    C3 = _daddsp(C3, _dmpysp(B3, _ftof2(D3, D3)));

    /* calculate F = A - B *inv(D) * conj(B) -- Hermitian */
    dtemp = _dmpysp(B0, B0);
    F0 = A0 - D0 * (_hif2(dtemp) + _lof2(dtemp));
    dtemp = _dmpysp(B1, B1);
    F0 -= D3 * (_hif2(dtemp) + _lof2(dtemp));
    //results           =   _cmpysp(B1, dtemp1);
    //dtemp         =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp = _complex_conjugate_mpysp(B1, dtemp1);
    F0 -= 2.f * _hif2(dtemp);

    dtemp = _dmpysp(B2, B2);
    F3 = A3 - D0 * (_hif2(dtemp) + _lof2(dtemp));
    dtemp = _dmpysp(B3, B3);
    F3 -= D3 * (_hif2(dtemp) + _lof2(dtemp));
    //results           =   _cmpysp(B3, dtemp2);
    //dtemp         =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp = _complex_conjugate_mpysp(B3, dtemp2);
    F3 -= 2.f * _hif2(dtemp);

    //results           =   _cmpysp(B2, C0);
    //F1                =   _dsubsp(A1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    F1 = _dsubsp(A1, _complex_conjugate_mpysp(B2, C0));
    //results           =   _cmpysp(B3, C1);
    //F1                =   _dsubsp(F1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    F1 = _dsubsp(F1, _complex_conjugate_mpysp(B3, C1));

    /* Calculate F = inv(F) */
    dtemp = _dmpysp(F1, F1);
    detA = F0 * F3 - _hif2(dtemp) - _lof2(dtemp);
    oneOverdetA = _rcpsp(detA);
    oneOverdetA = oneOverdetA * (2.f - detA * oneOverdetA);
    oneOverdetA = oneOverdetA * (2.f - detA * oneOverdetA);

    ftemp1 = F0 * oneOverdetA;
    F0 = F3 * oneOverdetA;
    F3 = ftemp1;
    F1 = _dmpysp(F1, _ftof2(-oneOverdetA, -oneOverdetA));

    /* NW output */
    _amem8_f2(&output[0]) = _ftof2(F0, 0.f);
    _amem8_f2(&output[1]) = F1;
    _amem8_f2(&output[5]) = _ftof2(F3, 0.f);
    //_amem8_f2(&output[4]) =   _lltof2(_f2toll(F1) ^ 0x0000000080000000);
    _amem8_f2(&output[4]) = _ftof2(_hif2(F1), -_lof2(F1));

    /* NE output = - F * C, SW = conj(NW)*/
    //results           =   _cmpysp(F1, C2);
    //dtemp         =   _daddsp(_hif2_128(results), _lof2_128(results));
    dtemp = _complex_mpysp(F1, C2);
    dtemp1 = dtemp;
    dtemp = _daddsp(dtemp, _dmpysp(C0, _ftof2(F0, F0)));
    _amem8_f2(&output[2]) = _ftof2(-_hif2(dtemp), -_lof2(dtemp));
    _amem8_f2(&output[8]) = _ftof2(-_hif2(dtemp), _lof2(dtemp));
    //_amem8_f2(&output[2])     =   _lltof2(_f2toll(dtemp) ^ 0x8000000080000000);
    //_amem8_f2(&output[8])     =   _lltof2(_f2toll(dtemp) ^ 0x8000000000000000);

    //results           =   _cmpysp(F1, C3);
    //dtemp         =   _daddsp(_hif2_128(results), _lof2_128(results));
    dtemp = _complex_mpysp(F1, C3);
    dtemp2 = dtemp;
    dtemp = _daddsp(dtemp, _dmpysp(C1, _ftof2(F0, F0)));
    B1 = dtemp;
    _amem8_f2(&output[3]) = _ftof2(-_hif2(dtemp), -_lof2(dtemp));
    _amem8_f2(&output[12]) = _ftof2(-_hif2(dtemp), _lof2(dtemp));
    //_amem8_f2(&output[3])     =   _lltof2(_f2toll(dtemp) ^ 0x8000000080000000);
    //_amem8_f2(&output[12])    =   _lltof2(_f2toll(dtemp) ^ 0x8000000000000000);

    //results           =   _cmpysp(F1, C0);
    //dtemp         =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp = _complex_conjugate_mpysp(F1, C0);
    dtemp = _daddsp(dtemp, _dmpysp(C2, _ftof2(F3, F3)));
    _amem8_f2(&output[6]) = _ftof2(-_hif2(dtemp), -_lof2(dtemp));
    _amem8_f2(&output[9]) = _ftof2(-_hif2(dtemp), _lof2(dtemp));
    //_amem8_f2(&output[6])     =   _lltof2(_f2toll(dtemp) ^ 0x8000000080000000);
    //_amem8_f2(&output[9])     =   _lltof2(_f2toll(dtemp) ^ 0x8000000000000000);

    //results           =   _cmpysp(F1, C1);
    //dtemp         =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp = _complex_conjugate_mpysp(F1, C1);
    dtemp = _daddsp(dtemp, _dmpysp(C3, _ftof2(F3, F3)));
    B3 = dtemp;
    _amem8_f2(&output[7]) = _ftof2(-_hif2(dtemp), -_lof2(dtemp));
    _amem8_f2(&output[13]) = _ftof2(-_hif2(dtemp), _lof2(dtemp));
    //_amem8_f2(&output[7])     =   _lltof2(_f2toll(dtemp) ^ 0x8000000080000000);
    //_amem8_f2(&output[13])    =   _lltof2(_f2toll(dtemp) ^ 0x8000000000000000);

    /* SE output */
    /* inv(D) - conj(C) * inv(F) * C, whrer C = B * inv(D) */
    dtemp = _dmpysp(C0, C0);
    A0 = D0 + F0 * (_hif2(dtemp) + _lof2(dtemp));
    dtemp = _dmpysp(C2, C2);
    A0 += F3 * (_hif2(dtemp) + _lof2(dtemp));
    //results           =   _cmpysp(C0, dtemp1);
    //dtemp         =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp = _complex_conjugate_mpysp(C0, dtemp1);
    A0 += 2.f * _hif2(dtemp);

    dtemp = _dmpysp(C1, C1);
    A3 = D3 + F0 * (_hif2(dtemp) + _lof2(dtemp));
    dtemp = _dmpysp(C3, C3);
    A3 += F3 * (_hif2(dtemp) + _lof2(dtemp));
    //results           =   _cmpysp(C1, dtemp2);
    //dtemp         =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp = _complex_conjugate_mpysp(C1, dtemp2);
    A3 += 2.f * _hif2(dtemp);
    _amem8_f2(&output[10]) = _ftof2(A0, 0.f);
    _amem8_f2(&output[15]) = _ftof2(A3, 0.f);

    //results           =   _cmpysp(C0, B1);
    //dtemp         =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp = _complex_conjugate_mpysp(C0, B1);
    A1 = _daddsp(D1, dtemp);
    //results           =   _cmpysp(C2, B3);
    //A1                =   _daddsp(A1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    A1 = _daddsp(A1, _complex_conjugate_mpysp(C2, B3));

    _amem8_f2(&output[11]) = A1;
    _amem8_f2(&output[14]) = _ftof2(_hif2(A1), -_lof2(A1));
    //_amem8_f2(&output[14])=   _lltof2(_f2toll(A1) ^ 0x0000000080000000);
}

void MATRIX_Mult4x4fp(cplxf_t *A, cplxf_t *B, cplxf_t *C)

{

#if GENERICCODE
    int32_t jj, kk, mm;
    float   tempre, tempim;

    for ( jj = 0; jj < 4; jj++ )
    {
        for ( kk = 0; kk < 4; kk++ )
        {
            tempre      =   0.f;
            tempim      =   0.f;
            for (mm = 0; mm < 4; mm++)
            {
                tempre      +=  A[4 * jj + mm].real * B[4 * mm + kk].real - A[4 * jj + mm].imag * B[4 * mm + kk].imag;
                tempim      +=  A[4 * jj + mm].imag * B[4 * mm + kk].real + A[4 * jj + mm].real * B[4 * mm + kk].imag;
            }
            C[4 * jj + kk].real     =   tempre;
            C[4 * jj + kk].imag     =   tempim;
        }
    }
#else
    int32_t kk, mm;
    __float2_t *A1, *A2, *A3, *A4;
    __float2_t *C1, *C2, *C3, *C4;
    __float2_t dtemp, dtemp1, dtemp2, dtemp3;

    A1 = (__float2_t*) &A[0];
    A2 = (__float2_t*) &A[4];
    A3 = (__float2_t*) &A[8];
    A4 = (__float2_t*) &A[12];
    C1 = (__float2_t*) &C[0];
    C2 = (__float2_t*) &C[4];
    C3 = (__float2_t*) &C[8];
    C4 = (__float2_t*) &C[12];
    for (kk = 0; kk < 4; kk++)
    {
        dtemp = _ftof2(0.f, 0.f);
        dtemp1 = dtemp;
        dtemp2 = dtemp;
        dtemp3 = dtemp;
#ifdef _TMS320C6x
        #pragma UNROLL(4);
        #endif
        for (mm = 0; mm < 4; mm++)
        {
            dtemp = _daddsp(dtemp, _complex_mpysp(_amem8_f2(&A1[mm]),
            _amem8_f2(&B[4 * mm + kk])));
            dtemp1 = _daddsp(dtemp1, _complex_mpysp(_amem8_f2(&A2[mm]),
            _amem8_f2(&B[4 * mm + kk])));
            dtemp2 = _daddsp(dtemp2, _complex_mpysp(_amem8_f2(&A3[mm]),
            _amem8_f2(&B[4 * mm + kk])));
            dtemp3 = _daddsp(dtemp3, _complex_mpysp(_amem8_f2(&A4[mm]),
            _amem8_f2(&B[4 * mm + kk])));
        }
        _amem8_f2(&C1[kk]) = dtemp;
        _amem8_f2(&C2[kk]) = dtemp1;
        _amem8_f2(&C3[kk]) = dtemp2;
        _amem8_f2(&C4[kk]) = dtemp3;
    }
#endif
}

/**
 *  \fn     void MATRIX_single8x8MatInv(
 *                            IN      Cplx32  * restrict Input,
 *                            IN      int32_t * RESTRICT scratch,
 *                            OUT     cplx32_t  * Output);
 *
 *  \brief   Single 8x8 matrix (complex positive semi-definitive Hermitian) inversion using block wise method.
 *
 *  \param[in]    Input
 *              Input 8x8 matrix that needs to be inversed, stored sequentially
 *              in format (1,1), (1,2).... (1,8), (2,1), (2,2)...(8,8).
 *              Must be aligned to 8-byte boundary.
 *
 *  \param[out]    scratch
 *              Input pointer to the scratch memory. Must be of size 7 * 2 * 16 = 288 32-bit words.
 *              Must be aligned to 8-byte boundary.
 *
 *  \param[out]   Output
 *              Output 8x8 matrix. Stored sequentially as the input.
 *              Must be aligned to 8-byte boundary.
 *
 *  \pre     Input matrix must be complex positive semi-definitive Hermitian matrix
 *
 *  \post
 *
 *  \sa
 *
 */
void MATRIX_single8x8MatInv(cplxf_t *Input, int32_t *scratch, cplxf_t *output)
{
    cplxf_t *A;
    cplxf_t *B;
    cplxf_t *C;
    cplxf_t *D;
    cplxf_t *T;
    cplxf_t *invT;
    cplxf_t *invA;
    int32_t jj, offset;
    int32_t scratchIndx;
    __float2_t dtemp1;

    scratchIndx = 0;
    A = (cplxf_t*) &scratch[scratchIndx];
    scratchIndx += 16 * 2;
    B = (cplxf_t*) &scratch[scratchIndx];
    scratchIndx += 16 * 2;
    C = (cplxf_t*) &scratch[scratchIndx];
    scratchIndx += 16 * 2;
    D = (cplxf_t*) &scratch[scratchIndx];
    scratchIndx += 16 * 2;
    T = (cplxf_t*) &scratch[scratchIndx];
    scratchIndx += 16 * 2;
    invA = (cplxf_t*) &scratch[scratchIndx];
    scratchIndx += 16 * 2;
    invT = (cplxf_t*) &scratch[scratchIndx];
    scratchIndx += 16 * 2;

    /*Copy A */
    _amem8_f2(&A[0]) = _amem8_f2(&Input[0]);
    dtemp1 = _amem8_f2(&Input[1]);
    _amem8_f2(&A[1]) = dtemp1;
    _amem8_f2(&A[4]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    dtemp1 = _amem8_f2(&Input[2]);
    _amem8_f2(&A[2]) = dtemp1;
    _amem8_f2(&A[8]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    dtemp1 = _amem8_f2(&Input[3]);
    _amem8_f2(&A[3]) = dtemp1;
    _amem8_f2(&A[12]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    _amem8_f2(&A[5]) = _amem8_f2(&Input[8]);
    dtemp1 = _amem8_f2(&Input[9]);
    _amem8_f2(&A[6]) = dtemp1;
    _amem8_f2(&A[9]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    dtemp1 = _amem8_f2(&Input[10]);
    _amem8_f2(&A[7]) = dtemp1;
    _amem8_f2(&A[13]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    _amem8_f2(&A[10]) = _amem8_f2(&Input[15]);
    dtemp1 = _amem8_f2(&Input[16]);
    _amem8_f2(&A[11]) = dtemp1;
    _amem8_f2(&A[14]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    _amem8_f2(&A[15]) = _amem8_f2(&Input[21]);

    /*copy C */
    offset = 4;
    for (jj = 0; jj < 4; jj++)
    {
        dtemp1 = _amem8_f2(&Input[jj + offset]);
        _amem8_f2(&C[4 * jj]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    }
    offset = 11;
    for (jj = 0; jj < 4; jj++)
    {
        dtemp1 = _amem8_f2(&Input[jj + offset]);
        _amem8_f2(&C[4 * jj + 1]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    }
    offset = 17;
    for (jj = 0; jj < 4; jj++)
    {
        dtemp1 = _amem8_f2(&Input[jj + offset]);
        _amem8_f2(&C[4 * jj + 2]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    }
    offset = 22;
    for (jj = 0; jj < 4; jj++)
    {
        dtemp1 = _amem8_f2(&Input[jj + offset]);
        _amem8_f2(&C[4 * jj + 3]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    }

    /*Copy D */
    _amem8_f2(&D[0]) = _amem8_f2(&Input[26]);
    dtemp1 = _amem8_f2(&Input[27]);
    _amem8_f2(&D[1]) = dtemp1;
    _amem8_f2(&D[4]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    dtemp1 = _amem8_f2(&Input[28]);
    _amem8_f2(&D[2]) = dtemp1;
    _amem8_f2(&D[8]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    dtemp1 = _amem8_f2(&Input[29]);
    _amem8_f2(&D[3]) = dtemp1;
    _amem8_f2(&D[12]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    _amem8_f2(&D[5]) = _amem8_f2(&Input[30]);
    dtemp1 = _amem8_f2(&Input[31]);
    _amem8_f2(&D[6]) = dtemp1;
    _amem8_f2(&D[9]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    dtemp1 = _amem8_f2(&Input[32]);
    _amem8_f2(&D[7]) = dtemp1;
    _amem8_f2(&D[13]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    _amem8_f2(&D[10]) = _amem8_f2(&Input[33]);
    dtemp1 = _amem8_f2(&Input[34]);
    _amem8_f2(&D[11]) = dtemp1;
    _amem8_f2(&D[14]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    _amem8_f2(&D[15]) = _amem8_f2(&Input[35]);

    /* calculate inv(A)*/
    MATRIX_4x4_BWInversionfp(A, invA);

    /* calculate tempC = C*inv(A) */
    //tsc_in            =   clock();
    MATRIX_Mult4x4fp(C, invA, B);
    //tsc_out           =   clock();
    //CycleCounts       =   tsc_out - tsc_in;

    /* calculate T = D - C*inv(A)*B */
    //results           =   _cmpysp(_amem8_f2(&C[0]), _amem8_f2(&B[0]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&C[0]), _amem8_f2(&B[0]));

    //results           =   _cmpysp(_amem8_f2(&C[1]), _amem8_f2(&B[1]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[1]), _amem8_f2(&B[1])));

    //results           =   _cmpysp(_amem8_f2(&C[2]), _amem8_f2(&B[2]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[2]), _amem8_f2(&B[2])));

    //results           =   _cmpysp(_amem8_f2(&C[3]), _amem8_f2(&B[3]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[3]), _amem8_f2(&B[3])));

    dtemp1 = _dsubsp(_amem8_f2(&D[0]), dtemp1);
    _amem8_f2(&T[0]) = dtemp1;

    //results           =   _cmpysp(_amem8_f2(&C[4]), _amem8_f2(&B[0]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&C[4]), _amem8_f2(&B[0]));

    //results           =   _cmpysp(_amem8_f2(&C[5]), _amem8_f2(&B[1]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[5]), _amem8_f2(&B[1])));

    //results           =   _cmpysp(_amem8_f2(&C[6]), _amem8_f2(&B[2]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[6]), _amem8_f2(&B[2])));

    //results           =   _cmpysp(_amem8_f2(&C[7]), _amem8_f2(&B[3]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[7]), _amem8_f2(&B[3])));

    dtemp1 = _dsubsp(_amem8_f2(&D[1]), dtemp1);
    _amem8_f2(&T[1]) = dtemp1;
    _amem8_f2(&T[4]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    //_amem8_f2(&T[4])  =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

    //results           =   _cmpysp(_amem8_f2(&C[4]), _amem8_f2(&B[4]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&C[4]), _amem8_f2(&B[4]));

    //results           =   _cmpysp(_amem8_f2(&C[5]), _amem8_f2(&B[5]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[5]), _amem8_f2(&B[5])));

    //results           =   _cmpysp(_amem8_f2(&C[6]), _amem8_f2(&B[6]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[6]), _amem8_f2(&B[6])));

    //results           =   _cmpysp(_amem8_f2(&C[7]), _amem8_f2(&B[7]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[7]), _amem8_f2(&B[7])));

    dtemp1 = _dsubsp(_amem8_f2(&D[5]), dtemp1);
    _amem8_f2(&T[5]) = dtemp1;

    //results           =   _cmpysp(_amem8_f2(&C[8]), _amem8_f2(&B[0]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&C[8]), _amem8_f2(&B[0]));

    //results           =   _cmpysp(_amem8_f2(&C[9]), _amem8_f2(&B[1]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[9]), _amem8_f2(&B[1])));

    //results           =   _cmpysp(_amem8_f2(&C[10]), _amem8_f2(&B[2]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[10]), _amem8_f2(&B[2])));

    //results           =   _cmpysp(_amem8_f2(&C[11]), _amem8_f2(&B[3]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[11]), _amem8_f2(&B[3])));

    dtemp1 = _dsubsp(_amem8_f2(&D[2]), dtemp1);
    _amem8_f2(&T[2]) = dtemp1;
    _amem8_f2(&T[8]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    //_amem8_f2(&T[8])  =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

    //results           =   _cmpysp(_amem8_f2(&C[8]), _amem8_f2(&B[4]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&C[8]), _amem8_f2(&B[4]));

    //results           =   _cmpysp(_amem8_f2(&C[9]), _amem8_f2(&B[5]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[9]), _amem8_f2(&B[5])));

    //results           =   _cmpysp(_amem8_f2(&C[10]), _amem8_f2(&B[6]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[10]), _amem8_f2(&B[6])));

    //results           =   _cmpysp(_amem8_f2(&C[11]), _amem8_f2(&B[7]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[11]), _amem8_f2(&B[7])));

    dtemp1 = _dsubsp(_amem8_f2(&D[6]), dtemp1);
    _amem8_f2(&T[6]) = dtemp1;
    _amem8_f2(&T[9]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    //_amem8_f2(&T[9])  =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

    //results           =   _cmpysp(_amem8_f2(&C[8]), _amem8_f2(&B[8]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&C[8]), _amem8_f2(&B[8]));

    //results           =   _cmpysp(_amem8_f2(&C[9]), _amem8_f2(&B[9]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[9]), _amem8_f2(&B[9])));

    //results           =   _cmpysp(_amem8_f2(&C[10]), _amem8_f2(&B[10]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[10]), _amem8_f2(&B[10])));

    //results           =   _cmpysp(_amem8_f2(&C[11]), _amem8_f2(&B[11]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[11]), _amem8_f2(&B[11])));

    dtemp1 = _dsubsp(_amem8_f2(&D[10]), dtemp1);
    _amem8_f2(&T[10]) = dtemp1;

    //results           =   _cmpysp(_amem8_f2(&C[12]), _amem8_f2(&B[0]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&C[12]), _amem8_f2(&B[0]));

    //results           =   _cmpysp(_amem8_f2(&C[13]), _amem8_f2(&B[1]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[13]), _amem8_f2(&B[1])));

    //results           =   _cmpysp(_amem8_f2(&C[14]), _amem8_f2(&B[2]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[14]), _amem8_f2(&B[2])));

    //results           =   _cmpysp(_amem8_f2(&C[15]), _amem8_f2(&B[3]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[15]), _amem8_f2(&B[3])));

    dtemp1 = _dsubsp(_amem8_f2(&D[3]), dtemp1);
    _amem8_f2(&T[3]) = dtemp1;
    _amem8_f2(&T[12]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    //_amem8_f2(&T[12]) =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

    //results           =   _cmpysp(_amem8_f2(&C[12]), _amem8_f2(&B[4]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&C[12]), _amem8_f2(&B[4]));

    //results           =   _cmpysp(_amem8_f2(&C[13]), _amem8_f2(&B[5]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[13]), _amem8_f2(&B[5])));

    //results           =   _cmpysp(_amem8_f2(&C[14]), _amem8_f2(&B[6]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[14]), _amem8_f2(&B[6])));

    //results           =   _cmpysp(_amem8_f2(&C[15]), _amem8_f2(&B[7]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[15]), _amem8_f2(&B[7])));

    dtemp1 = _dsubsp(_amem8_f2(&D[7]), dtemp1);
    _amem8_f2(&T[7]) = dtemp1;
    _amem8_f2(&T[13]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    //_amem8_f2(&T[13]) =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

    //results           =   _cmpysp(_amem8_f2(&C[12]), _amem8_f2(&B[8]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&C[12]), _amem8_f2(&B[8]));

    //results           =   _cmpysp(_amem8_f2(&C[13]), _amem8_f2(&B[9]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[13]), _amem8_f2(&B[9])));

    //results           =   _cmpysp(_amem8_f2(&C[14]), _amem8_f2(&B[10]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[14]), _amem8_f2(&B[10])));

    //results           =   _cmpysp(_amem8_f2(&C[15]), _amem8_f2(&B[11]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[15]), _amem8_f2(&B[11])));

    dtemp1 = _dsubsp(_amem8_f2(&D[11]), dtemp1);
    _amem8_f2(&T[11]) = dtemp1;
    _amem8_f2(&T[14]) = _ftof2(_hif2(dtemp1), -_lof2(dtemp1));
    //_amem8_f2(&T[14]) =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

    //results           =   _cmpysp(_amem8_f2(&C[12]), _amem8_f2(&B[12]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&C[12]), _amem8_f2(&B[12]));

    //results           =   _cmpysp(_amem8_f2(&C[13]), _amem8_f2(&B[13]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[13]), _amem8_f2(&B[13])));

    //results           =   _cmpysp(_amem8_f2(&C[14]), _amem8_f2(&B[14]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[14]), _amem8_f2(&B[14])));

    //results           =   _cmpysp(_amem8_f2(&C[15]), _amem8_f2(&B[15]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&C[15]), _amem8_f2(&B[15])));

    dtemp1 = _dsubsp(_amem8_f2(&D[15]), dtemp1);
    _amem8_f2(&T[15]) = dtemp1;

    /* Calculate T = inv(T) */
    MATRIX_4x4_BWInversionfp(T, invT);

    /* output SE 4x4 matrix */
    _amem8_f2(&output[26]) = _amem8_f2(&invT[0]);
    _amem8_f2(&output[27]) = _amem8_f2(&invT[1]);
    _amem8_f2(&output[28]) = _amem8_f2(&invT[2]);
    _amem8_f2(&output[29]) = _amem8_f2(&invT[3]);
    _amem8_f2(&output[30]) = _amem8_f2(&invT[5]);
    _amem8_f2(&output[31]) = _amem8_f2(&invT[6]);
    _amem8_f2(&output[32]) = _amem8_f2(&invT[7]);
    _amem8_f2(&output[33]) = _amem8_f2(&invT[10]);
    _amem8_f2(&output[34]) = _amem8_f2(&invT[11]);
    _amem8_f2(&output[35]) = _amem8_f2(&invT[15]);

    /* output SW and NE 4x4 matrix = - invA*B*invT */
    MATRIX_Mult4x4fp(invT, B, C);
    jj = 0;
    dtemp1 = _amem8_f2(&C[jj * 4 + 0]);
    _amem8_f2(&output[4]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    dtemp1 = _amem8_f2(&C[jj * 4 + 1]);
    _amem8_f2(&output[11]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    dtemp1 = _amem8_f2(&C[jj * 4 + 2]);
    _amem8_f2(&output[17]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    dtemp1 = _amem8_f2(&C[jj * 4 + 3]);
    _amem8_f2(&output[22]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    jj = 1;
    dtemp1 = _amem8_f2(&C[jj * 4 + 0]);
    _amem8_f2(&output[5]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    dtemp1 = _amem8_f2(&C[jj * 4 + 1]);
    _amem8_f2(&output[12]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    dtemp1 = _amem8_f2(&C[jj * 4 + 2]);
    _amem8_f2(&output[18]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    dtemp1 = _amem8_f2(&C[jj * 4 + 3]);
    _amem8_f2(&output[23]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    jj = 2;
    dtemp1 = _amem8_f2(&C[jj * 4 + 0]);
    _amem8_f2(&output[6]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    dtemp1 = _amem8_f2(&C[jj * 4 + 1]);
    _amem8_f2(&output[13]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    dtemp1 = _amem8_f2(&C[jj * 4 + 2]);
    _amem8_f2(&output[19]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    dtemp1 = _amem8_f2(&C[jj * 4 + 3]);
    _amem8_f2(&output[24]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    jj = 3;
    dtemp1 = _amem8_f2(&C[jj * 4 + 0]);
    _amem8_f2(&output[7]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    dtemp1 = _amem8_f2(&C[jj * 4 + 1]);
    _amem8_f2(&output[14]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    dtemp1 = _amem8_f2(&C[jj * 4 + 2]);
    _amem8_f2(&output[20]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));
    dtemp1 = _amem8_f2(&C[jj * 4 + 3]);
    _amem8_f2(&output[25]) = _ftof2(-_hif2(dtemp1), _lof2(dtemp1));

#if 0
    zeros       =   _ftof2(0.f, 0.f);
    for ( jj = 0; jj < 4; jj++ )
    {
        dtemp1                              =   _dsubsp(zeros, _amem8_f2(&C[jj * 4 + 0]));
        _amem8_f2(&output[(jj + 4) * 8 + 0])    =   dtemp1;
        _amem8_f2(&output[jj + 4])      =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

        dtemp1                              =   _dsubsp(zeros, _amem8_f2(&C[jj * 4 + 1]));
        _amem8_f2(&output[(jj + 4) * 8 + 1])    =   dtemp1;
        _amem8_f2(&output[jj + 12])     =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

        dtemp1                              =   _dsubsp(zeros, _amem8_f2(&C[jj * 4 + 2]));
        _amem8_f2(&output[(jj + 4) * 8 + 2])    =   dtemp1;
        _amem8_f2(&output[jj + 20])     =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

        dtemp1                              =   _dsubsp(zeros, _amem8_f2(&C[jj * 4 + 3]));
        _amem8_f2(&output[(jj + 4) * 8 + 3])    =   dtemp1;
        _amem8_f2(&output[jj + 28])     =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);
    }
#endif

    /* output NW 4x4 matrix = invA + invA*B*invT*C*invA */
    /* output NW 4x4 matrix = A + conj(B)*C */
    //results           =   _cmpysp(_amem8_f2(&B[0]), _amem8_f2(&C[0]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&B[0]), _amem8_f2(&C[0]));

    //results           =   _cmpysp(_amem8_f2(&B[4]), _amem8_f2(&C[4]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[4]), _amem8_f2(&C[4])));

    //results           =   _cmpysp(_amem8_f2(&B[8]), _amem8_f2(&C[8]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[8]), _amem8_f2(&C[8])));

    //results           =   _cmpysp(_amem8_f2(&B[12]), _amem8_f2(&C[12]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[12]), _amem8_f2(&C[12])));

    dtemp1 = _daddsp(_amem8_f2(&invA[0]), dtemp1);
    _amem8_f2(&output[0]) = dtemp1;

    //results           =   _cmpysp(_amem8_f2(&B[0]), _amem8_f2(&C[1]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&B[0]), _amem8_f2(&C[1]));

    //results           =   _cmpysp(_amem8_f2(&B[4]), _amem8_f2(&C[5]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[4]), _amem8_f2(&C[5])));

    //results           =   _cmpysp(_amem8_f2(&B[8]), _amem8_f2(&C[9]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[8]), _amem8_f2(&C[9])));

    //results           =   _cmpysp(_amem8_f2(&B[12]), _amem8_f2(&C[13]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[12]), _amem8_f2(&C[13])));

    dtemp1 = _daddsp(_amem8_f2(&invA[1]), dtemp1);
    _amem8_f2(&output[1]) = dtemp1;
    //_amem8_f2(&output[8]) =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

    //results           =   _cmpysp(_amem8_f2(&B[1]), _amem8_f2(&C[1]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&B[1]), _amem8_f2(&C[1]));

    //results           =   _cmpysp(_amem8_f2(&B[5]), _amem8_f2(&C[5]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[5]), _amem8_f2(&C[5])));

    //results           =   _cmpysp(_amem8_f2(&B[9]), _amem8_f2(&C[9]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[9]), _amem8_f2(&C[9])));

    //results           =   _cmpysp(_amem8_f2(&B[13]), _amem8_f2(&C[13]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[13]), _amem8_f2(&C[13])));

    dtemp1 = _daddsp(_amem8_f2(&invA[5]), dtemp1);
    //_amem8_f2(&output[9]) =   dtemp1;
    _amem8_f2(&output[8]) = dtemp1;

    //results           =   _cmpysp(_amem8_f2(&B[0]), _amem8_f2(&C[2]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&B[0]), _amem8_f2(&C[2]));

    //results           =   _cmpysp(_amem8_f2(&B[4]), _amem8_f2(&C[6]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[4]), _amem8_f2(&C[6])));

    //results           =   _cmpysp(_amem8_f2(&B[8]), _amem8_f2(&C[10]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[8]), _amem8_f2(&C[10])));

    //results           =   _cmpysp(_amem8_f2(&B[12]), _amem8_f2(&C[14]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[12]), _amem8_f2(&C[14])));

    dtemp1 = _daddsp(_amem8_f2(&invA[2]), dtemp1);
    _amem8_f2(&output[2]) = dtemp1;
    //_amem8_f2(&output[16])    =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

    //results           =   _cmpysp(_amem8_f2(&B[1]), _amem8_f2(&C[2]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&B[1]), _amem8_f2(&C[2]));

    //results           =   _cmpysp(_amem8_f2(&B[5]), _amem8_f2(&C[6]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[5]), _amem8_f2(&C[6])));

    //results           =   _cmpysp(_amem8_f2(&B[9]), _amem8_f2(&C[10]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[9]), _amem8_f2(&C[10])));

    //results           =   _cmpysp(_amem8_f2(&B[13]), _amem8_f2(&C[14]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[13]), _amem8_f2(&C[14])));

    dtemp1 = _daddsp(_amem8_f2(&invA[6]), dtemp1);
    _amem8_f2(&output[9]) = dtemp1;
    //_amem8_f2(&output[10])    =   dtemp1;
    //_amem8_f2(&output[17])    =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

    //results           =   _cmpysp(_amem8_f2(&B[2]), _amem8_f2(&C[2]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&B[2]), _amem8_f2(&C[2]));

    //results           =   _cmpysp(_amem8_f2(&B[6]), _amem8_f2(&C[6]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[6]), _amem8_f2(&C[6])));

    //results           =   _cmpysp(_amem8_f2(&B[10]), _amem8_f2(&C[10]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[10]), _amem8_f2(&C[10])));

    //results           =   _cmpysp(_amem8_f2(&B[14]), _amem8_f2(&C[14]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[14]), _amem8_f2(&C[14])));

    dtemp1 = _daddsp(_amem8_f2(&invA[10]), dtemp1);
    _amem8_f2(&output[15]) = dtemp1;
    //_amem8_f2(&output[18])    =   dtemp1;

    //results           =   _cmpysp(_amem8_f2(&B[0]), _amem8_f2(&C[3]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&B[0]), _amem8_f2(&C[3]));

    //results           =   _cmpysp(_amem8_f2(&B[4]), _amem8_f2(&C[7]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[4]), _amem8_f2(&C[7])));

    //results           =   _cmpysp(_amem8_f2(&B[8]), _amem8_f2(&C[11]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[8]), _amem8_f2(&C[11])));

    //results           =   _cmpysp(_amem8_f2(&B[12]), _amem8_f2(&C[15]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[12]), _amem8_f2(&C[15])));

    dtemp1 = _daddsp(_amem8_f2(&invA[3]), dtemp1);
    _amem8_f2(&output[3]) = dtemp1;
    //_amem8_f2(&output[24])    =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

    //results           =   _cmpysp(_amem8_f2(&B[1]), _amem8_f2(&C[3]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&B[1]), _amem8_f2(&C[3]));

    //results           =   _cmpysp(_amem8_f2(&B[5]), _amem8_f2(&C[7]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[5]), _amem8_f2(&C[7])));

    //results           =   _cmpysp(_amem8_f2(&B[9]), _amem8_f2(&C[11]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[9]), _amem8_f2(&C[11])));

    //results           =   _cmpysp(_amem8_f2(&B[13]), _amem8_f2(&C[15]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[13]), _amem8_f2(&C[15])));

    dtemp1 = _daddsp(_amem8_f2(&invA[7]), dtemp1);
    _amem8_f2(&output[10]) = dtemp1;
    //_amem8_f2(&output[11])    =   dtemp1;
    //_amem8_f2(&output[25])    =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

    //results           =   _cmpysp(_amem8_f2(&B[2]), _amem8_f2(&C[3]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&B[2]), _amem8_f2(&C[3]));

    //results           =   _cmpysp(_amem8_f2(&B[6]), _amem8_f2(&C[7]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[6]), _amem8_f2(&C[7])));

//  results         =   _cmpysp(_amem8_f2(&B[10]), _amem8_f2(&C[11]));
//  dtemp1          =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[10]), _amem8_f2(&C[11])));

//  results         =   _cmpysp(_amem8_f2(&B[14]), _amem8_f2(&C[15]));
//  dtemp1          =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[14]), _amem8_f2(&C[15])));

    dtemp1 = _daddsp(_amem8_f2(&invA[11]), dtemp1);
    _amem8_f2(&output[16]) = dtemp1;
    //_amem8_f2(&output[19])    =   dtemp1;
    //_amem8_f2(&output[26])    =   _lltof2(_f2toll(dtemp1) ^ 0x0000000080000000);

    //results           =   _cmpysp(_amem8_f2(&B[3]), _amem8_f2(&C[3]));
    //dtemp1            =   _dsubsp(_hif2_128(results), _lof2_128(results));
    dtemp1 = _complex_conjugate_mpysp(_amem8_f2(&B[3]), _amem8_f2(&C[3]));

    //results           =   _cmpysp(_amem8_f2(&B[7]), _amem8_f2(&C[7]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[7]), _amem8_f2(&C[7])));

    //results           =   _cmpysp(_amem8_f2(&B[11]), _amem8_f2(&C[11]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[11]), _amem8_f2(&C[11])));

    //results           =   _cmpysp(_amem8_f2(&B[15]), _amem8_f2(&C[15]));
    //dtemp1            =   _daddsp(dtemp1, _dsubsp(_hif2_128(results), _lof2_128(results)));
    dtemp1 = _daddsp(
            dtemp1,
            _complex_conjugate_mpysp(_amem8_f2(&B[15]), _amem8_f2(&C[15])));

    dtemp1 = _daddsp(_amem8_f2(&invA[15]), dtemp1);
    _amem8_f2(&output[21]) = dtemp1;
    //_amem8_f2(&output[27])    =   dtemp1;
}

void ODDemo_Heatmap_steeringVecGen(ODDemo_DataPathObj *obj)
{

    uint32_t i, j;
    double ftemp1, freal1, fimag1, frealJ, fimagJ;

    // Ant0's steeringVec is 1 for all angle possiblities, so we don't save them
    for (i = 0; i < obj->steeringVecSize; i++)
    {
        ftemp1 = (double) sin(
                (-obj->estAngleRange + (double) i * obj->estAngleResolution)
                        * (double) RADARDEMO_AOAESTBF_PIOVER180);
        freal1 = (double) cos(-RADARDEMO_AOAESTBF_PI * ftemp1);
        fimag1 = (double) sin(-RADARDEMO_AOAESTBF_PI * ftemp1);
        frealJ = freal1;
        fimagJ = fimag1;
        obj->steeringVec[(obj->nRxAnt - 1) * i + 0].real = (float) frealJ;
        obj->steeringVec[(obj->nRxAnt - 1) * i + 0].imag = (float) fimagJ;
        for (j = 2; j < obj->nRxAnt; j++)
        {
            ftemp1 = frealJ;
            frealJ = frealJ * freal1 - fimagJ * fimag1;
            fimagJ = ftemp1 * fimag1 + fimagJ * freal1;
            obj->steeringVec[(obj->nRxAnt - 1) * i + j - 1].real =
                    (float) frealJ;
            obj->steeringVec[(obj->nRxAnt - 1) * i + j - 1].imag =
                    (float) fimagJ;
        }
    }
}

void ODDemo_Heatmap_aoaEstCaponBF_range(ODDemo_DataPathObj *obj)
{
    /*Calculate covariance matrix and invert */
    ODDemo_Heatmap_aoaEstCaponBF_covInv(
            (uint8_t) (obj->fallBackToConvBFFlag ^ 1),
            (uint8_t) obj->clutterRemovalFlag, oddemo_parms.gamma,
            (int32_t) obj->nRxAnt, (int32_t) obj->nChirps,
            (int32_t) obj->log2Chirps, (int32_t*) &obj->scratchPad[0],
            (cplx16_t*) obj->inputAntSamples, (cplxf_t*) &obj->invRnMatrix[0]);

    /* Capon beamforming */
    ODDemo_Heatmap_aoaEstCaponBF_heatmap(
            (uint8_t) (obj->fallBackToConvBFFlag ^ 1), (int32_t) obj->nRxAnt,
            (int32_t) obj->steeringVecSize, (cplxf_t*) obj->steeringVec,
            (cplxf_t*) &obj->invRnMatrix[0], (float*) obj->rangeAzimuthHeatMap);
}

#if 0
void ODDemo_Heatmap_aoaEstCaponBF_doppler(ODDemo_DataPathObj *obj, uint8_t rangeIndx, uint8_t azimuthIndx)
{
    ODDEMO_Heatmap_aoaEstCaponBF_dopplerEstInput(
            (uint8_t)           1,
            (int32_t)           obj->nRxAnt,
            (int32_t)           obj->nChirps,
            //(cplx16_t *)        obj->inputAntSamples,
            (cplx16_t *)        obj->swapBuf,
            (cplxf_t *)         &obj->steeringVec[azimuthIndx * (obj->nRxAnt - 1)],
            (cplxf_t  *)        &obj->invRnMatrix[rangeIndx * ODDEMO_INVRNMATRIX_BUF_SIZE],
            (int32_t *)         &obj->scratchPad[2* obj->nChirps],
            (float)             1.0,
            (float *)           obj->bfOutput);
}
#endif

/*!
 *   \fn     RADARDEMO_aoaEstCaponBF_covInv
 *
 *   \brief   Per range bin, estimate the covariance matrices from input 1D FFT results, and calculate the inverse of these matrices.
 *
 *   \param[in]    invFlag
 *               Flag to indicate matrix inversion will be performed.
 *               If set to 1, output invRnMatrices will contain inversion of covariance matrices.
 *               If set to 0, output invRnMatrices will contain covariance matrices without inversion.
 *
 *   \param[in]    clutterRmFlag
 *               Flag to indicate clutter removal will be performed if set to 1. Otherwise, disabled.
 *
 *   \param[in]    gamma
 *               Scaling factor for diagnal loading.
 *
 *   \param[in]    nRxAnt
 *               number of antenna
 *
 *   \param[in]    nChirps
 *               number of input chirps
 *
 *   \param[in]    scratch
 *               scratch memory, must be of size TBD!!!
 *               Must be aligned to 8-byte boundary.
 *
 *   \param[in]    inputAntSamples
 *               input samples from radar cube (1D FFT output) for the current (one) range bin to be processed
 *               Must be aligned to 8-byte boundary.
 *
 *   \param[out]    invRnMatrices
 *               Output inverse of covariance matrices for the current range bin, in order of upper triangle of nRxAnt x nRxAnt Hermitian matrix.
 *               Must be aligned to 8-byte boundary.
 *
 *   \ret       none
 *
 *   \pre       none
 *
 *   \post      none
 *
 *   IMPORTANT NOTE: Previous versions of this function expected input as 16-bit
 *   complex with imag in the low 16-bits. This version expects real in the low
 *   16-bits as it now does real/imag swap inline.
 */

//made these global for debug inspection:
int32_t sumVal[2];
int32_t meanVal;

void ODDemo_Heatmap_aoaEstCaponBF_covInv(uint8_t invFlag, uint8_t clutterRmFlag,
                                         float gamma, int32_t nRxAnt,
                                         int32_t nChirps, int32_t log2Chirps,
                                         int32_t *scratch,
                                         cplx16_t *inputAntSamples,
                                         cplxf_t *invRnMatrices)
{
    int32_t antIdx, chirpIdx, i, j, scratchOffset, rnIdx;
    cplx16_t *input1;
    cplx16_t *input2;
    __float2_t *Rn;
    int64_t lltemp, llinput1, llinput2;
    __float2_t acc, acc1, acc2, acc3, scale2;
    int32_t itemp1;
    int32_t itemp2;
    int32_t roundup;
    cplxf_t *invRn;
    float ftemp;

    scratchOffset = 0;
    Rn = (__float2_t*) &scratch[scratchOffset];
    scratchOffset = scratchOffset + 2 * nRxAnt * (1 + (nRxAnt >> 1)); /*72 32-bit word for 8 antennas*/

    ftemp = _rcpsp((float) nChirps);
    scale2 = _ftof2(ftemp, ftemp);

    /* Do clutter removal, summing to 32-bits. Real/imag order does not matter here
     * as real and imag are averaged independently. */
    if (clutterRmFlag)
    {
        roundup = 1 << (log2Chirps - 1);

        for (antIdx = 0; antIdx < nRxAnt; antIdx++)
        {
            cmplx32ImRe_t *pSumVal = (cmplx32ImRe_t*) sumVal;
            cmplx16ImRe_t *pMeanVal = (cmplx16ImRe_t*) &meanVal;
            int16_t *input = (int16_t*) &inputAntSamples[antIdx * nChirps];

            mmwavelib_vecsum((int16_t*) input, (int32_t*) sumVal,
                             (int32_t) nChirps);

            pMeanVal->real = (pSumVal->real + roundup) >> log2Chirps;
            pMeanVal->imag = (pSumVal->imag + roundup) >> log2Chirps;

            mmwavelib_vecsubc((int16_t*) input, (int16_t*) input,
                              (uint32_t) meanVal, (int32_t) nChirps);
        }
    }

    /* Rn estimation */
    rnIdx = 0;
    for (antIdx = 0; antIdx < nRxAnt; antIdx++)
    {
        input1 = (cplx16_t*) &inputAntSamples[antIdx * nChirps];

        //i = antIdx case
        acc = _ftof2(0.f, 0.f);
        acc1 = _ftof2(0.f, 0.f);
        acc2 = _ftof2(0.f, 0.f);
        acc3 = _ftof2(0.f, 0.f);
        for (chirpIdx = 0; chirpIdx < nChirps; chirpIdx += 8)
        {
            llinput1 = _amem8(&input1[chirpIdx]);       //load 2 complex symbols
            itemp2 = _hill(llinput1);                      //Grab the 1st symbol
            itemp2 = _packlh2(itemp2, itemp2);             //Swap real/imag
            itemp1 = _packhl2(itemp2, _ssub2(0, itemp2));  //Negate imag
            lltemp = _cmpy(itemp2, itemp1);                //Complex multiply

            itemp2 = _loll(llinput1);                      //Grab the 2nd symbol
            itemp2 = _packlh2(itemp2, itemp2);             //Swap real/imag
            itemp1 = _packhl2(itemp2, _ssub2(0, itemp2));  //Negate imag
            lltemp = _dadd(lltemp, _cmpy(itemp2, itemp1)); //Complex mult and add 1st symbol
            acc = _daddsp(acc, _dintsp(lltemp));   //add both symbols to the sum

            llinput1 = _amem8(&input1[chirpIdx + 2]); //load next 2 complex symbols
            itemp2 = _hill(llinput1);
            itemp2 = _packlh2(itemp2, itemp2);
            itemp1 = _packhl2(itemp2, _ssub2(0, itemp2));
            lltemp = _cmpy(itemp2, itemp1);

            itemp2 = _loll(llinput1);
            itemp2 = _packlh2(itemp2, itemp2);
            itemp1 = _packhl2(itemp2, _ssub2(0, itemp2));
            lltemp = _dadd(lltemp, _cmpy(itemp2, itemp1));
            acc2 = _daddsp(acc2, _dintsp(lltemp));

            llinput1 = _amem8(&input1[chirpIdx + 4]); //load next 2 complex symbols
            itemp2 = _hill(llinput1);
            itemp2 = _packlh2(itemp2, itemp2);
            itemp1 = _packhl2(itemp2, _ssub2(0, itemp2));
            lltemp = _cmpy(itemp2, itemp1);

            itemp2 = _loll(llinput1);
            itemp2 = _packlh2(itemp2, itemp2);
            itemp1 = _packhl2(itemp2, _ssub2(0, itemp2));
            lltemp = _dadd(lltemp, _cmpy(itemp2, itemp1));
            acc1 = _daddsp(acc1, _dintsp(lltemp));

            llinput1 = _amem8(&input1[chirpIdx + 6]); //load next 2 complex symbols
            itemp2 = _hill(llinput1);
            itemp2 = _packlh2(itemp2, itemp2);
            itemp1 = _packhl2(itemp2, _ssub2(0, itemp2));
            lltemp = _cmpy(itemp2, itemp1);

            itemp2 = _loll(llinput1);
            itemp2 = _packlh2(itemp2, itemp2);
            itemp1 = _packhl2(itemp2, _ssub2(0, itemp2));
            lltemp = _dadd(lltemp, _cmpy(itemp2, itemp1));
            acc3 = _daddsp(acc3, _dintsp(lltemp));

        }
        acc = _daddsp(acc, acc1);
        acc = _daddsp(acc, acc2);
        acc = _daddsp(acc, acc3);
        acc = _dmpysp(acc, scale2);
        _amem8_f2(&Rn[rnIdx++]) = _ftof2(_hif2(acc), 0.f);

        for (i = antIdx + 1; i < nRxAnt; i++)
        {
            input2 = (cplx16_t*) &inputAntSamples[i * nChirps];

            acc = _ftof2(0.f, 0.f);
            acc1 = _ftof2(0.f, 0.f);
            for (chirpIdx = 0; chirpIdx < nChirps; chirpIdx += 4)
            {
                llinput1 = _amem8(&input1[chirpIdx]);
                llinput2 = _amem8(&input2[chirpIdx]);
                itemp2 = _hill(llinput2);
                itemp2 = _packlh2(itemp2, itemp2); //swap
                itemp1 = _packhl2(itemp2, _ssub2(0, itemp2));
                itemp2 = _hill(llinput1);
                itemp2 = _packlh2(itemp2, itemp2); //swap
                lltemp = _cmpy(itemp2, itemp1);

                itemp2 = _loll(llinput2);
                itemp2 = _packlh2(itemp2, itemp2); //swap
                itemp1 = _packhl2(itemp2, _ssub2(0, itemp2));
                itemp2 = _loll(llinput1);
                itemp2 = _packlh2(itemp2, itemp2); //swap
                lltemp = _dadd(lltemp, _cmpy(itemp2, itemp1));
                acc = _daddsp(acc, _dintsp(lltemp));

                llinput1 = _amem8(&input1[chirpIdx + 2]);
                llinput2 = _amem8(&input2[chirpIdx + 2]);
                itemp2 = _hill(llinput2);
                itemp2 = _packlh2(itemp2, itemp2); //swap
                itemp1 = _packhl2(itemp2, _ssub2(0, itemp2));
                itemp2 = _hill(llinput1);
                itemp2 = _packlh2(itemp2, itemp2); //swap
                lltemp = _cmpy(itemp2, itemp1);

                itemp2 = _loll(llinput2);
                itemp2 = _packlh2(itemp2, itemp2); //swap
                itemp1 = _packhl2(itemp2, _ssub2(0, itemp2));
                itemp2 = _loll(llinput1);
                itemp2 = _packlh2(itemp2, itemp2); //swap
                lltemp = _dadd(lltemp, _cmpy(itemp2, itemp1));
                acc1 = _daddsp(acc1, _dintsp(lltemp));
            }
            acc = _daddsp(acc, acc1);
            acc = _dmpysp(acc, scale2);
            _amem8_f2(&Rn[rnIdx++]) = acc;
        }
    }

    if (invFlag)
    {
        /*add diagnal loading */
        j = 0;
        ftemp = 0;
        itemp1 = nRxAnt;
        for (i = 0; i < nRxAnt; i++)
        {
            ftemp += _hif2(_amem8_f2(&Rn[j]));
            j += itemp1;
            itemp1--;
        }
        if (nRxAnt == 8)
            ftemp *= 0.125f;
        else if (nRxAnt == 4)
            ftemp *= 0.25f;
        j = 0;
        ftemp *= gamma;
        acc = _ftof2(ftemp, 0.f);
        itemp1 = nRxAnt;
        for (i = 0; i < nRxAnt; i++)
        {
            _amem8_f2(&Rn[j]) = _daddsp(_amem8_f2(&Rn[j]), acc);
            j += itemp1;
            itemp1--;
        }

        /* matrix inversion */
        if (nRxAnt == 8)
        {
            MATRIX_single8x8MatInv((cplxf_t*) Rn, &scratch[scratchOffset],
                                   invRnMatrices);
        }
        else if (nRxAnt == 4)
        {
            cplxf_t *inputRn;
            __float2_t f2temp1;

            inputRn = (cplxf_t*) &scratch[scratchOffset]; // size 2 * 16
            invRn = (cplxf_t*) &scratch[scratchOffset + 2 * 16]; // size 2 * 16

            rnIdx = 0;
            for (i = 0; i < nRxAnt; i++)
            {
                for (j = i; j < nRxAnt; j++)
                {
                    _amem8_f2(&inputRn[i * nRxAnt + j]) = _amem8_f2(
                            &Rn[rnIdx++]);
                }
            }
            for (i = 1; i < nRxAnt; i++)
            {
                for (j = 0; j < i; j++)
                {
                    f2temp1 = _amem8_f2(&inputRn[j * nRxAnt + i]);

                    _amem8_f2(&inputRn[i * nRxAnt + j]) = _ftof2(
                            _hif2(f2temp1), -_lof2(f2temp1));
                }
            }

            MATRIX_4x4_BWInversionfp(inputRn, invRn);
            rnIdx = 0;
            for (i = 0; i < nRxAnt; i++)
            {
                for (j = i; j < nRxAnt; j++)
                {
                    _amem8_f2(&invRnMatrices[rnIdx++]) = _amem8_f2(
                            &invRn[i * nRxAnt + j]);
                }
            }
        }
    }
    else
    {

        if (nRxAnt == 8)
            rnIdx = 36;
        else if (nRxAnt == 4)
            rnIdx = 10;

        for (i = 0; i < rnIdx; i++)
        {
            _amem8_f2(&invRnMatrices[i]) = _amem8_f2(&Rn[i]);
        }
    }
}

/*!
 *   \fn     RADARDEMO_aoaEstCaponBF_heatmap
 *
 *   \brief   Use Capon beamforming to generate range azimuth heatmap per range bin.
 *
 *   \param[in]    bfFlag
 *               Flag to indicate which covariance matrix based beamforming will be performed.
 *               If set to 1, Capon BF.
 *               If set to 0, conventional BF.
 *
 *   \param[in]    nRxAnt
 *               number of antenna
 *
 *   \param[in]    numAzimuthBins
 *               number of Azimuth bins
 *
 *   \param[in]    steeringVec
 *              steering vector for beamforming.
 *
 *   \param[in]    scratch
 *               scratch memory, must be of size of 360 32-bit words for 8 antennas, and 82 32-bit words for 4 antennas.
 *               Must be aligned to 8-byte boundary.
 *
 *   \param[in]    invRnMatrices
 *               Output inverse of covariance matrices or the current range bin, in order of upper triangle of nRxAnt x nRxAnt Hermitian matrix.
 *               Must be aligned to 8-byte boundary.
 *
 *   \param[out]    rangeAzimuthHeatMap
 *               Output range azimuth heatmap, in the format of numInputRangeBins by numAzimuthBins
 *               Must be aligned to 8-byte boundary.
 *
 *   \ret       none
 *
 *   \pre       none
 *
 *   \post      none
 *
 *
 */
void ODDemo_Heatmap_aoaEstCaponBF_heatmap(uint8_t bfFlag, int32_t nRxAnt,
                                          int32_t numAzimuthBins,
                                          cplxf_t *steeringVec,
                                          cplxf_t *invRnMatrices,
                                          float *rangeAzimuthHeatMap)
{
    int32_t azimIdx;
    __float2_t *steeringVecPtr;
    __float2_t *invRnMatPtr;
    float *heatMapPtr;
    __float2_t f2temp, steerVecIn;
    float output, result;

    if (nRxAnt == 8)
    {
        steeringVecPtr = (__float2_t*) &steeringVec[0];
        invRnMatPtr = (__float2_t*) &invRnMatrices[0];

        //this may need to change depending on which dimension to search first.
        heatMapPtr = (float*) &rangeAzimuthHeatMap[0];
        for (azimIdx = 0; azimIdx < numAzimuthBins; azimIdx++)
        {
            output = _hif2(_amem8_f2(&invRnMatPtr[0]));
            output += _hif2(_amem8_f2(&invRnMatPtr[8]));
            output += _hif2(_amem8_f2(&invRnMatPtr[15]));
            output += _hif2(_amem8_f2(&invRnMatPtr[21]));
            output += _hif2(_amem8_f2(&invRnMatPtr[26]));
            output += _hif2(_amem8_f2(&invRnMatPtr[30]));
            output += _hif2(_amem8_f2(&invRnMatPtr[33]));
            output += _hif2(_amem8_f2(&invRnMatPtr[35]));

            f2temp = _amem8_f2(&invRnMatPtr[1]);
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[9]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[16]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[22]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[27]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[31]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[34]));
            steerVecIn = _amem8_f2(steeringVecPtr++);
            f2temp = _dmpysp(f2temp, steerVecIn);
            output += 2.f * (_hif2(f2temp) + _lof2(f2temp));

            f2temp = _amem8_f2(&invRnMatPtr[2]);
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[10]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[17]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[23]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[28]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[32]));
            steerVecIn = _amem8_f2(steeringVecPtr++);
            f2temp = _dmpysp(f2temp, steerVecIn);
            output += 2.f * (_hif2(f2temp) + _lof2(f2temp));

            f2temp = _amem8_f2(&invRnMatPtr[3]);
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[11]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[18]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[24]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[29]));
            steerVecIn = _amem8_f2(steeringVecPtr++);
            f2temp = _dmpysp(f2temp, steerVecIn);
            output += 2.f * (_hif2(f2temp) + _lof2(f2temp));

            f2temp = _amem8_f2(&invRnMatPtr[4]);
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[12]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[19]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[25]));
            steerVecIn = _amem8_f2(steeringVecPtr++);
            f2temp = _dmpysp(f2temp, steerVecIn);
            output += 2.f * (_hif2(f2temp) + _lof2(f2temp));

            f2temp = _amem8_f2(&invRnMatPtr[5]);
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[13]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[20]));
            steerVecIn = _amem8_f2(steeringVecPtr++);
            f2temp = _dmpysp(f2temp, steerVecIn);
            output += 2.f * (_hif2(f2temp) + _lof2(f2temp));

            f2temp = _amem8_f2(&invRnMatPtr[6]);
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[14]));
            steerVecIn = _amem8_f2(steeringVecPtr++);
            f2temp = _dmpysp(f2temp, steerVecIn);
            output += 2.f * (_hif2(f2temp) + _lof2(f2temp));

            f2temp = _amem8_f2(&invRnMatPtr[7]);
            steerVecIn = _amem8_f2(steeringVecPtr++);
            f2temp = _dmpysp(f2temp, steerVecIn);
            output += 2.f * (_hif2(f2temp) + _lof2(f2temp));

            result = _rcpsp(output);
            result = result * (2.f - output * result);
            result = result * (2.f - output * result);
            //result        =   result * result;

            if (!bfFlag)
                result = output;
            heatMapPtr[azimIdx] = result;
        }
    }
    else if (nRxAnt == 4)
    {
        steeringVecPtr = (__float2_t*) &steeringVec[0];
        invRnMatPtr = (__float2_t*) &invRnMatrices[0];

        //this may need to change depending on which dimension to search first.
        heatMapPtr = (float*) &rangeAzimuthHeatMap[0];
        for (azimIdx = 0; azimIdx < numAzimuthBins; azimIdx++)
        {
            output = _hif2(_amem8_f2(&invRnMatPtr[0]));
            output += _hif2(_amem8_f2(&invRnMatPtr[4]));
            output += _hif2(_amem8_f2(&invRnMatPtr[7]));
            output += _hif2(_amem8_f2(&invRnMatPtr[9]));

            f2temp = _amem8_f2(&invRnMatPtr[1]);
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[5]));
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[8]));
            steerVecIn = _amem8_f2(steeringVecPtr++);
            f2temp = _dmpysp(f2temp, steerVecIn);
            output += 2.f * (_hif2(f2temp) + _lof2(f2temp));

            f2temp = _amem8_f2(&invRnMatPtr[2]);
            f2temp = _daddsp(f2temp, _amem8_f2(&invRnMatPtr[6]));
            steerVecIn = _amem8_f2(steeringVecPtr++);
            f2temp = _dmpysp(f2temp, steerVecIn);
            output += 2.f * (_hif2(f2temp) + _lof2(f2temp));

            f2temp = _amem8_f2(&invRnMatPtr[3]);
            steerVecIn = _amem8_f2(steeringVecPtr++);
            f2temp = _dmpysp(f2temp, steerVecIn);
            output += 2.f * (_hif2(f2temp) + _lof2(f2temp));

            result = _rcpsp(output);
            result = result * (2.f - output * result);
            result = result * (2.f - output * result);
            //result        =   result * result;

            if (!bfFlag)
                result = output;
            heatMapPtr[azimIdx] = result;
        }
    }
}

/*!
 *   \fn     ODDEMO_Heatmap_aoaEstCaponBF_dopplerEstInput
 *
 *   \brief   Calculate the beam forming output for doppler estimation, for the current range bin and azimuth bin.
 *
 *   \param[in]    bfFlag
 *               Flag to indicate which covariance matrix based beamforming will be performed.
 *               If set to 1, Capon BF.
 *               If set to 0, conventional BF.
 *
 *   \param[in]    nRxAnt
 *               number of antenna
 *
 *   \param[in]    nChirps
 *               number of chirps
 *
 *   \param[in]    inputAntSamples
 *              Input 1D FFT results for the current range bin.
 *              Must be aligned to 8-byte boundary.
 *
 *   \param[in]    steeringVec
 *              steering vector for beamforming for the current azimuth bin.
 *
 *   \param[in]    invRnMatrices
 *               Output inverse of covariance matrices for the current range bin, in order of numInputRangeBins by upper triangle of nRxAnt x nRxAnt Hermitian matrix.
 *               Must be aligned to 8-byte boundary.
 *
 *   \param[in]    scratch
 *               Output inverse of covariance matrices for the current range bin, in order of numInputRangeBins by upper triangle of nRxAnt x nRxAnt Hermitian matrix.
 *               Must be aligned to 8-byte boundary.
 *
 *   \param[in]    rangeAzimuthHeatMap
 *               input range azimuth heatmap for current range bin and azimuth bin
 *
 *   \param[out]    bfOutput
 *               Beamforming output for the current range bin and azimuth bin.
 *               Must be in the order of real0, imag0, real1, imag1... as required by DSP LIB single precision floating-point FFT.
 *               Must be aligned to 8-byte boundary.
 *
 *   \ret       none
 *
 *   \pre       none
 *
 *   \post      none
 *
 *
 */
void ODDEMO_Heatmap_aoaEstCaponBF_dopplerEstInput(uint8_t bfFlag,
                                                  int32_t nRxAnt,
                                                  int32_t nChirps,
                                                  cplx16_t *inputAntSamples,
                                                  cplxf_t *steeringVec,
                                                  cplxf_t *invRnMatrices,
                                                  int32_t *scratch,
                                                  float rangeAzimuthHeatMap,
                                                  float *bfOutput)
{
    int32_t chirpIdx, i;
    __float2_t *restrict bweights;
    __float2_t *restrict rnPtr;
    __float2_t *restrict steerVecPtr;
    __float2_t f2temp1, acc2, scale2;

    /* beamforming weights = A*Rn, beamforming output = input(1x8)*A*Rn */
    scale2 = _ftof2(rangeAzimuthHeatMap, rangeAzimuthHeatMap);
    bweights = (__float2_t*) &scratch[0];
    rnPtr = (__float2_t*) &invRnMatrices[0];
    steerVecPtr = (__float2_t*) &steeringVec[0];

    if (nRxAnt == 8)
    {
        int32_t *restrict inPtr1;
        int32_t *restrict inPtr2;
        int32_t *restrict inPtr3;
        int32_t *restrict inPtr4;
        int32_t *restrict inPtr5;
        int32_t *restrict inPtr6;
        int32_t *restrict inPtr7;
        int32_t *restrict inPtr8;

        inPtr1 = (int32_t*) &inputAntSamples[0 * nChirps];
        inPtr2 = (int32_t*) &inputAntSamples[1 * nChirps];
        inPtr3 = (int32_t*) &inputAntSamples[2 * nChirps];
        inPtr4 = (int32_t*) &inputAntSamples[3 * nChirps];
        inPtr5 = (int32_t*) &inputAntSamples[4 * nChirps];
        inPtr6 = (int32_t*) &inputAntSamples[5 * nChirps];
        inPtr7 = (int32_t*) &inputAntSamples[6 * nChirps];
        inPtr8 = (int32_t*) &inputAntSamples[7 * nChirps];

        if (bfFlag)
        {
            //bweights[0]
            acc2 = _amem8_f2(&rnPtr[0]);
            for (i = 1; i < nRxAnt; i++)
            {
                f2temp1 = _complex_conjugate_mpysp(
                        _amem8_f2(&rnPtr[i]), _amem8_f2(&steerVecPtr[i - 1]));
                acc2 = _daddsp(acc2, f2temp1);
            }
            _amem8_f2(&bweights[0]) = _dmpysp(acc2, scale2);

            //bweights[1]
            acc2 = _amem8_f2(&rnPtr[1]);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[8]),
            _amem8_f2(&steerVecPtr[0]));
            acc2 = _daddsp(acc2, f2temp1);
            for (i = 2; i < nRxAnt; i++)
            {
                f2temp1 = _complex_conjugate_mpysp(
                        _amem8_f2(&rnPtr[i + 7]),
                        _amem8_f2(&steerVecPtr[i - 1]));
                acc2 = _daddsp(acc2, f2temp1);
            }
            _amem8_f2(&bweights[1]) = _dmpysp(acc2, scale2);

            //bweights[2]
            acc2 = _amem8_f2(&rnPtr[2]);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[9]),
            _amem8_f2(&steerVecPtr[0]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[15]),
            _amem8_f2(&steerVecPtr[1]));
            acc2 = _daddsp(acc2, f2temp1);
            for (i = 3; i < nRxAnt; i++)
            {
                f2temp1 = _complex_conjugate_mpysp(
                        _amem8_f2(&rnPtr[i + 13]),
                        _amem8_f2(&steerVecPtr[i - 1]));
                acc2 = _daddsp(acc2, f2temp1);
            }
            _amem8_f2(&bweights[2]) = _dmpysp(acc2, scale2);

            //bweights[3]
            acc2 = _amem8_f2(&rnPtr[3]);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[10]),
            _amem8_f2(&steerVecPtr[0]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[16]),
            _amem8_f2(&steerVecPtr[1]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[21]),
            _amem8_f2(&steerVecPtr[2]));
            acc2 = _daddsp(acc2, f2temp1);
            for (i = 4; i < nRxAnt; i++)
            {
                f2temp1 = _complex_conjugate_mpysp(
                        _amem8_f2(&rnPtr[i + 18]),
                        _amem8_f2(&steerVecPtr[i - 1]));
                acc2 = _daddsp(acc2, f2temp1);
            }
            _amem8_f2(&bweights[3]) = _dmpysp(acc2, scale2);

            //bweights[4]
            acc2 = _amem8_f2(&rnPtr[4]);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[11]),
            _amem8_f2(&steerVecPtr[0]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[17]),
            _amem8_f2(&steerVecPtr[1]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[22]),
            _amem8_f2(&steerVecPtr[2]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[26]),
            _amem8_f2(&steerVecPtr[3]));
            acc2 = _daddsp(acc2, f2temp1);
            for (i = 5; i < nRxAnt; i++)
            {
                f2temp1 = _complex_conjugate_mpysp(
                        _amem8_f2(&rnPtr[i + 22]),
                        _amem8_f2(&steerVecPtr[i - 1]));
                acc2 = _daddsp(acc2, f2temp1);
            }
            _amem8_f2(&bweights[4]) = _dmpysp(acc2, scale2);

            //bweights[5]
            acc2 = _amem8_f2(&rnPtr[5]);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[12]),
            _amem8_f2(&steerVecPtr[0]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[18]),
            _amem8_f2(&steerVecPtr[1]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[23]),
            _amem8_f2(&steerVecPtr[2]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[27]),
            _amem8_f2(&steerVecPtr[3]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[30]),
            _amem8_f2(&steerVecPtr[4]));
            acc2 = _daddsp(acc2, f2temp1);
            for (i = 6; i < nRxAnt; i++)
            {
                f2temp1 = _complex_conjugate_mpysp(
                        _amem8_f2(&rnPtr[i + 25]),
                        _amem8_f2(&steerVecPtr[i - 1]));
                acc2 = _daddsp(acc2, f2temp1);
            }
            _amem8_f2(&bweights[5]) = _dmpysp(acc2, scale2);

            //bweights[6]
            acc2 = _amem8_f2(&rnPtr[6]);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[13]),
            _amem8_f2(&steerVecPtr[0]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[19]),
            _amem8_f2(&steerVecPtr[1]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[24]),
            _amem8_f2(&steerVecPtr[2]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[28]),
            _amem8_f2(&steerVecPtr[3]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[31]),
            _amem8_f2(&steerVecPtr[4]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[33]),
            _amem8_f2(&steerVecPtr[5]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_conjugate_mpysp(_amem8_f2(&rnPtr[34]),
            _amem8_f2(&steerVecPtr[6]));
            acc2 = _daddsp(acc2, f2temp1);
            _amem8_f2(&bweights[6]) = _dmpysp(acc2, scale2);

            //bweights[7]
            acc2 = _amem8_f2(&rnPtr[7]);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[14]),
            _amem8_f2(&steerVecPtr[0]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[20]),
            _amem8_f2(&steerVecPtr[1]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[25]),
            _amem8_f2(&steerVecPtr[2]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[29]),
            _amem8_f2(&steerVecPtr[3]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[32]),
            _amem8_f2(&steerVecPtr[4]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[34]),
            _amem8_f2(&steerVecPtr[5]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[35]),
            _amem8_f2(&steerVecPtr[6]));
            acc2 = _daddsp(acc2, f2temp1);
            _amem8_f2(&bweights[7]) = _dmpysp(acc2, scale2);
        }
        else
        {
            _amem8_f2(&bweights[0]) = _ftof2(1.f, 0.f);
            for (i = 0; i < nRxAnt - 1; i++)
            {
                _amem8_f2(&bweights[i + 1]) = _amem8_f2(&steerVecPtr[i]);
            }
        }

        for (chirpIdx = 0; chirpIdx < nChirps; chirpIdx++)
        {
            acc2 = _complex_mpysp(_dinthsp(_amem4(&inPtr1[chirpIdx])),
            _amem8_f2(&bweights[0]));
            acc2 = _daddsp(acc2,
                           _complex_mpysp(_dinthsp(_amem4(&inPtr2[chirpIdx])),
                           _amem8_f2(&bweights[1])));
            acc2 = _daddsp(acc2,
                           _complex_mpysp(_dinthsp(_amem4(&inPtr3[chirpIdx])),
                           _amem8_f2(&bweights[2])));
            acc2 = _daddsp(acc2,
                           _complex_mpysp(_dinthsp(_amem4(&inPtr4[chirpIdx])),
                           _amem8_f2(&bweights[3])));
            acc2 = _daddsp(acc2,
                           _complex_mpysp(_dinthsp(_amem4(&inPtr5[chirpIdx])),
                           _amem8_f2(&bweights[4])));
            acc2 = _daddsp(acc2,
                           _complex_mpysp(_dinthsp(_amem4(&inPtr6[chirpIdx])),
                           _amem8_f2(&bweights[5])));
            acc2 = _daddsp(acc2,
                           _complex_mpysp(_dinthsp(_amem4(&inPtr7[chirpIdx])),
                           _amem8_f2(&bweights[6])));
            acc2 = _daddsp(acc2,
                           _complex_mpysp(_dinthsp(_amem4(&inPtr8[chirpIdx])),
                           _amem8_f2(&bweights[7])));
            _amem8_f2(&bfOutput[2 * chirpIdx]) = _ftof2(_lof2(acc2),
            _hif2(acc2));
        }
    }
    else if (nRxAnt == 4)
    {
        int32_t *restrict inPtr1;
        int32_t *restrict inPtr2;
        int32_t *restrict inPtr3;
        int32_t *restrict inPtr4;

        inPtr1 = (int32_t*) &inputAntSamples[0 * nChirps];
        inPtr2 = (int32_t*) &inputAntSamples[1 * nChirps];
        inPtr3 = (int32_t*) &inputAntSamples[2 * nChirps];
        inPtr4 = (int32_t*) &inputAntSamples[3 * nChirps];
        if (bfFlag)
        {
            //bweights[0]
            acc2 = _amem8_f2(&rnPtr[0]);
            for (i = 1; i < nRxAnt; i++)
            {
                f2temp1 = _complex_conjugate_mpysp(
                        _amem8_f2(&rnPtr[i]), _amem8_f2(&steerVecPtr[i - 1]));
                acc2 = _daddsp(acc2, f2temp1);
            }
            _amem8_f2(&bweights[0]) = _dmpysp(acc2, scale2);

            //bweights[1]
            acc2 = _amem8_f2(&rnPtr[1]);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[4]),
            _amem8_f2(&steerVecPtr[0]));
            acc2 = _daddsp(acc2, f2temp1);
            for (i = 2; i < nRxAnt; i++)
            {
                f2temp1 = _complex_conjugate_mpysp(
                        _amem8_f2(&rnPtr[i + 3]),
                        _amem8_f2(&steerVecPtr[i - 1]));
                acc2 = _daddsp(acc2, f2temp1);
            }
            _amem8_f2(&bweights[1]) = _dmpysp(acc2, scale2);

            //bweights[2]
            acc2 = _amem8_f2(&rnPtr[2]);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[5]),
            _amem8_f2(&steerVecPtr[0]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[7]),
            _amem8_f2(&steerVecPtr[1]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_conjugate_mpysp(_amem8_f2(&rnPtr[8]),
            _amem8_f2(&steerVecPtr[2]));
            acc2 = _daddsp(acc2, f2temp1);
            _amem8_f2(&bweights[2]) = _dmpysp(acc2, scale2);

            //bweights[3]
            acc2 = _amem8_f2(&rnPtr[3]);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[6]),
            _amem8_f2(&steerVecPtr[0]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[8]),
            _amem8_f2(&steerVecPtr[1]));
            acc2 = _daddsp(acc2, f2temp1);
            f2temp1 = _complex_mpysp(_amem8_f2(&rnPtr[9]),
            _amem8_f2(&steerVecPtr[2]));
            acc2 = _daddsp(acc2, f2temp1);
            _amem8_f2(&bweights[3]) = _dmpysp(acc2, scale2);
        }
        else
        {
            _amem8_f2(&bweights[0]) = _ftof2(1.f, 0.f);
            for (i = 0; i < nRxAnt - 1; i++)
            {
                _amem8_f2(&bweights[i + 1]) = _amem8_f2(&steerVecPtr[i]);
            }
        }

        for (chirpIdx = 0; chirpIdx < nChirps; chirpIdx++)
        {
            acc2 = _complex_mpysp(_dinthsp(_amem4(&inPtr1[chirpIdx])),
            _amem8_f2(&bweights[0]));
            acc2 = _daddsp(acc2,
                           _complex_mpysp(_dinthsp(_amem4(&inPtr2[chirpIdx])),
                           _amem8_f2(&bweights[1])));
            acc2 = _daddsp(acc2,
                           _complex_mpysp(_dinthsp(_amem4(&inPtr3[chirpIdx])),
                           _amem8_f2(&bweights[2])));
            acc2 = _daddsp(acc2,
                           _complex_mpysp(_dinthsp(_amem4(&inPtr4[chirpIdx])),
                           _amem8_f2(&bweights[3])));
            _amem8_f2(&bfOutput[chirpIdx]) = _ftof2(_lof2(acc2), _hif2(acc2));
        }

    }
}

//This function captures min, max and average for each heatmap row.
//It currently also creates average empty FOV power during the first several
//seconds only if "rowNoise" CLI messages have not been received.
void ODDemo_Heatmap_get_stats(float *heatmap)
{
    uint16_t rng_idx, az_idx;
    float min, max, sum;
    float *heat = heatmap;

    oddemo_heatmap_max = 0.0;

    for (rng_idx = 0; rng_idx < ODDEMO_MAX_RANGE; rng_idx++)
    {
        min = 999999.9;
        max = 0.0;
        sum = 0.0;

        //Find the min and max of all cells on this row
        for (az_idx = 0; az_idx < ODDEMO_MAX_AZIMUTH; az_idx++)
        {
            if (heat[az_idx] > max)
                max = heat[az_idx];

            if (heat[az_idx] < min)
                min = heat[az_idx];

            sum += heat[az_idx];
        }

        if (max > oddemo_heatmap_max) //keep the overall heatmap maximum
            oddemo_heatmap_max = max;

        oddemo_stats[rng_idx].min = min;
        oddemo_stats[rng_idx].max = max;
        oddemo_stats[rng_idx].avg = sum / ODDEMO_MAX_AZIMUTH;

        heat += ODDEMO_MAX_AZIMUTH; //point to the next row
    }

    //Do startup row calibration. If oddemo_row_init is greater than this value
    //then we received the values from CLI and do not need to do it again.
    if (oddemo_row_init < ODDEMO_ROW_NOISE_FRAMES)
    {
        if (oddemo_row_init == 0)
            memset(oddemo_row_noise, 0, sizeof(oddemo_row_noise));

        if (oddemo_row_init >= 4) //ignore the first 4 frames
        {
            for (rng_idx = 0; rng_idx < ODDEMO_MAX_RANGE; rng_idx++)
                oddemo_row_noise[rng_idx] += oddemo_stats[rng_idx].avg;
        }

        oddemo_row_init++;

        if (oddemo_row_init == ODDEMO_ROW_NOISE_FRAMES) //done collecting data
        {
            for (rng_idx = 0; rng_idx < ODDEMO_MAX_RANGE; rng_idx++)
                oddemo_row_noise[rng_idx] /= (ODDEMO_ROW_NOISE_FRAMES - 4.0);
        }
    }
}

//This function scans the heatmap row by row looking for cases where the minimum
//value is greater than the saved noise floor for the row.  If true, the minimum
//is subtracted from the row and then rescaled to restore the peak maximum to its
//original value.
void ODDemo_Heatmap_arc_removal(float *heatmap)
{
    uint16_t rng_idx, az_idx;
    float min;
    float row_range, factor;
    float temp, noise_level;
    float *heat = heatmap;

    ODDemo_Heatmap_get_stats(heatmap);

    for (rng_idx = 0; rng_idx < ODDEMO_MAX_RANGE; rng_idx++)
    {
        noise_level = oddemo_row_noise[rng_idx];
        min = oddemo_stats[rng_idx].min;

        if (min > noise_level)
        {
            row_range = oddemo_stats[rng_idx].max - min;
            factor = 1.0 + ((min - noise_level) / row_range);

            //This loop will expand the range of values on this row so that the min value is
            //moved to the noise floor, and the max of any peak is kept at the same value.
            for (az_idx = 0; az_idx < ODDEMO_MAX_AZIMUTH; az_idx++)
            {
                temp = (heat[az_idx] - min) * factor;
                heat[az_idx] = temp + noise_level;
            }
        }

        heat += ODDEMO_MAX_AZIMUTH; //point to the next row
    }
}

//This function scales a heatmap and converts it to 16-bit integer data.
//It uses the maximum cell value found during arc removal.
void ODDemo_Heatmap_scale_heatmap16(float *heatin, uint16_t *heatout)
{
    uint16_t idx;
    float scale;
    uint16_t *out = heatout;
    float *ptr = heatin;

    ptr = heatin;

    if (oddemo_heatmap_max <= (float) ODDEMO_HEATMAP_MAX16) //no scaling needed
    {
        for (idx = 0; idx < ODDEMO_ANGLEHEATMAP_BUF_SIZE; idx++)
            *out++ = (uint16_t) *ptr++;
    }
    else //scale the heatmap to the max value
    {
        scale = (float) ODDEMO_HEATMAP_MAX16 / oddemo_heatmap_max;

        for (idx = 0; idx < ODDEMO_ANGLEHEATMAP_BUF_SIZE; idx++)
            *out++ = (uint16_t) (*ptr++ * scale);
    }
}

//This function scales a heatmap and converts it to 8-bit integer data.
//It uses the maximum cell value found during arc removal.
void ODDemo_Heatmap_scale_heatmap8(float *heatin, uint8_t *heatout)
{
    uint16_t idx;
    float scale;
    uint8_t *out = heatout;
    float *ptr = heatin;

    ptr = heatin;

    if (oddemo_heatmap_max <= (float) ODDEMO_HEATMAP_MAX8) //no scaling needed
    {
        for (idx = 0; idx < ODDEMO_ANGLEHEATMAP_BUF_SIZE; idx++)
            *out++ = (uint8_t) *ptr++;
    }
    else //scale the heatmap to the max value
    {
        scale = (float) ODDEMO_HEATMAP_MAX8 / oddemo_heatmap_max;

        for (idx = 0; idx < ODDEMO_ANGLEHEATMAP_BUF_SIZE; idx++)
            *out++ = (uint8_t) (*ptr++ * scale);
    }
}

/***
 * This function gets the peaks from the heatmap from a certain threshold.
 * Two peaks close to each other will be grouped. Max of 4 biggest peaks
 * will be returned.
 */
void heatmapGetPeaks(float *heatmap)
{
    uint16_t rng_idx, az_idx;
    float *heat = &heatmap[8 * ODDEMO_MAX_AZIMUTH];
    float maxPeaks[4];
    uint8_t numPeaks = 0;

    for (rng_idx = 8; rng_idx < ODDEMO_MAX_RANGE; rng_idx++)
    {
        for (az_idx = 0; az_idx < ODDEMO_MAX_AZIMUTH; az_idx++)
        {
//            float surroundingHeat = 0;
            if (heat[az_idx] < 2000.0f)
            {
                continue;
            }

            if (rng_idx > 0)
            {
                if (heat[az_idx - ODDEMO_MAX_AZIMUTH] > heat[az_idx])
                {
                    continue;
                }
//                surroundingHeat += heat[az_idx - ODDEMO_MAX_AZIMUTH];
                if (az_idx > 0)
                {
                    if (heat[az_idx - ODDEMO_MAX_AZIMUTH - 1] > heat[az_idx])
                    {
                        continue;
                    }
//                    surroundingHeat += heat[az_idx - ODDEMO_MAX_AZIMUTH - 1];
                }
                if (az_idx < ODDEMO_MAX_AZIMUTH - 1)
                {
                    if (heat[az_idx - ODDEMO_MAX_AZIMUTH + 1] > heat[az_idx])
                    {
                        continue;
                    }
//                    surroundingHeat += heat[az_idx - ODDEMO_MAX_AZIMUTH + 1];
                }
            }
            if (rng_idx < ODDEMO_MAX_RANGE - 1)
            {
                if (heat[az_idx + ODDEMO_MAX_AZIMUTH] > heat[az_idx])
                {
                    continue;
                }
//                surroundingHeat += heat[az_idx + ODDEMO_MAX_AZIMUTH];
                if (az_idx > 0)
                {
                    if (heat[az_idx + ODDEMO_MAX_AZIMUTH - 1] > heat[az_idx])
                    {
                        continue;
                    }
//                    surroundingHeat += heat[az_idx + ODDEMO_MAX_AZIMUTH - 1];
                }
                if (az_idx < ODDEMO_MAX_AZIMUTH - 1)
                {
                    if (heat[az_idx + ODDEMO_MAX_AZIMUTH + 1] > heat[az_idx])
                    {
                        continue;
                    }
//                    surroundingHeat += heat[az_idx + ODDEMO_MAX_AZIMUTH + 1];
                }
            }
            if (az_idx > 0)
            {
                if (heat[az_idx - 1] > heat[az_idx])
                {
                    continue;
                }
//                surroundingHeat += heat[az_idx - 1];
            }
            if (az_idx < ODDEMO_MAX_AZIMUTH - 1)
            {
                if (heat[az_idx + 1] > heat[az_idx])
                {
                    continue;
                }
//                surroundingHeat += heat[az_idx + 1];
            }
//            if ((surroundingHeat / 8) < heat[az_idx] / 4) {
//                continue;
//            }
            // It is a peak
            addPeakToArray(rng_idx, az_idx, heat[az_idx], maxPeaks, &numPeaks);
        }

        heat += ODDEMO_MAX_AZIMUTH;
    }
    if (numPeaks > 1)
    {
        uint16_t i;
        for (i = 0; i < numPeaks - 1; i++)
        {
            uint16_t j = 0;
            for (j = i + 1; j < numPeaks; j++)
            {
                if (peak_positions[i * 2 + 1] < peak_positions[j * 2 + 1])
                {
                    uint16_t temp = peak_positions[i * 2];
                    peak_positions[i * 2] = peak_positions[j * 2];
                    peak_positions[j * 2] = temp;
                    temp = peak_positions[i * 2 + 1];
                    peak_positions[i * 2 + 1] = peak_positions[j * 2 + 1];
                    peak_positions[j * 2 + 1] = temp;
                }
            }
        }
    }
    if (numPeaks > 0) {
        uint8_t i;
        for (i = numPeaks * 2; i < 8; i++)
        {
            peak_positions[i] = 0;
        }

        peak_positions[8] = numPeaks;
    }
}

void addPeakToArray(uint8_t rng_idx, uint8_t az_idx, float power,
                    float *maxPeaks, uint8_t *numPeaks)
{
    uint8_t i;
    uint8_t new_locations[8];
    uint8_t new_loc_idx = 0;
    float newMaxPeaks[4];
    uint8_t notBiggest = 0;
    for (i = 0; i < *numPeaks; i++)
    {
        if (notBiggest == 1)
        {
            new_locations[new_loc_idx * 2] = peak_positions[i * 2];
            new_locations[new_loc_idx * 2 + 1] = peak_positions[i * 2 + 1];
            newMaxPeaks[new_loc_idx] = maxPeaks[i];
            new_loc_idx++;
            continue;
        }
        float dis = sqrt(
                pow((int16_t) peak_positions[i * 2 + 1] - (int16_t) az_idx, 2)
                        + pow((int16_t) peak_positions[i * 2]
                                      - (int16_t) rng_idx,
                              2) * 1.0);
        if (dis < 10)
        {
            // One peak must go
            if (power < maxPeaks[i])
            {
                // Newly added peak is too close and smaller than existing peak in list
                new_locations[new_loc_idx * 2] = peak_positions[i * 2];
                new_locations[new_loc_idx * 2 + 1] = peak_positions[i * 2 + 1];
                newMaxPeaks[new_loc_idx] = maxPeaks[i];
                new_loc_idx++;
                notBiggest = 1;
                continue;
            }
        }
        else
        {
            // Peaks are added separately
            new_locations[new_loc_idx * 2] = peak_positions[i * 2];
            new_locations[new_loc_idx * 2 + 1] = peak_positions[i * 2 + 1];
            newMaxPeaks[new_loc_idx] = maxPeaks[i];
            new_loc_idx++;
        }
    }
    if (notBiggest == 0)
    {
        // Just add the peak to the list
        if (new_loc_idx < 4)
        {
            newMaxPeaks[new_loc_idx] = power;
            new_locations[new_loc_idx * 2] = rng_idx;
            new_locations[(new_loc_idx * 2) + 1] = az_idx;
            new_loc_idx++;
        }
        else
        {
            uint16_t j;
            for (j = 0; j < new_loc_idx; j++)
            {
                if (newMaxPeaks[j] >= power)
                {
                    continue;
                }
                newMaxPeaks[j] = power;
                new_locations[j * 2] = rng_idx;
                new_locations[(j * 2) + 1] = az_idx;
                break;
            }
        }
    }
    memcpy(&peak_positions[0], &new_locations[0], sizeof(uint8_t) * 8);
    memcpy(&maxPeaks[0], &newMaxPeaks[0], sizeof(float) * 4);
    *numPeaks = new_loc_idx;
}
