/**
 * oddemo_decision.c
 *
 * Functions to enable the decision making for Occupancy Detection
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

#include "../common/oddemo_common.h"
#include <math.h>
#include <float.h>
#include <oddemo_decision.h>

uint16_t oddemo_decision[ODDEMO_ZONE_PAIR];


void DSPF_sp_mat_vec_mul(float *x1, const int r1, const int c1,
                         float *x2, float *y)
{
    int i, j;
    float sum;

    // Multiply each row in x1 by vector x2.
    // The product of row m in x1 and vector x2 is placed
    // in position (m) in the result vector.
    for (i = 0; i < r1; i++)
    {
        sum = 0;

        for (j = 0; j < c1; j++)
            sum += x1[j + i * c1] * x2[j];

        y[i] = sum;
    }
}


/**
 *  This routine finds out the index of maximum number in the input array. This
 *  function returns the index of the greatest value.
 *
 *         @param x   Pointer to input array.
 *         @param nx  Number of inputs in the input array.
 *
 * @par Algorithm:
 * DSPF_sp_maxidx_cn.c is the natural C equivalent of the optimized
 * linear assembly code without restrictions. Note that the linear
 * assembly code is optimized and restrictions may apply.
 *
 * @par Assumptions:
 *   nx is a multiple of 4 and >= 4. <BR>
 *
 * @par Implementation Notes:
 * @b Interruptibility: The code is interruptible. <BR>
 * @b Endian Support: The code supports both big and little endian modes. <BR>
 *
 */
int DSPF_sp_maxidx(const float* x, const int nx)
{
    int   i, idx1, idx2, idx3, idx4;
    float max1, max2, max3, max4;
    __float2_t x_1_2, x_3_4;

    max1 = -FLT_MAX;
    max2 = -FLT_MAX;
    max3 = -FLT_MAX;
    max4 = -FLT_MAX;

//    _nassert((int)x % 8 == 0);
//    _nassert(nx % 4 == 0);
//    _nassert(nx > 0);

    #pragma MUST_ITERATE(1,,1)
    for (i = 0; i < nx; i+=4)
    {
        x_1_2 = _amem8_f2((void *)&x[i]);
        x_3_4 = _amem8_f2((void *)&x[i+2]);

        if (_lof2(x_1_2) > max1) {
            max1 = _lof2(x_1_2);
            idx1 = i;
        }

        if (_hif2(x_1_2) > max2) {
            max2 = _hif2(x_1_2);
            idx2 = i + 1;
        }

        if (_lof2(x_3_4) > max3) {
            max3 = _lof2(x_3_4);
            idx3 = i + 2;
        }

        if (_hif2(x_3_4) > max4) {
            max4 = _hif2(x_3_4);
            idx4 = i + 3;
        }
    }

    if (max2 > max1)
    {
        max1 = max2;
        idx1 = idx2;
    }
    if (max4 > max3)
    {
        max3 = max4;
        idx3 = idx4;
    }
    if (max3 > max1)
        return idx3;
    else
        return idx1;
}


void ODDemo_Decision_process(float          *coeffMatrix,
                             ODDEMO_Feature *featureVect,
                             int8_t         *decision)
{
    uint32_t idx;
    uint32_t result;
    float vec[ODDEMO_MATRIX_ROW_SIZE];
    float prob[ODDEMO_MATRIX_NUM_ROWS];

    //create a column vector
    vec[0] = 1.0;
    vec[1] = featureVect->powerMA[0];
    vec[2] = featureVect->powerMA[1];
    vec[3] = featureVect->powRatio[0];
    vec[4] = featureVect->powRatio[1];
    vec[5] = featureVect->crossCorr;

    DSPF_sp_mat_vec_mul(coeffMatrix,
                        ODDEMO_MATRIX_NUM_ROWS,
                        ODDEMO_MATRIX_ROW_SIZE,
                        vec,
                        prob); //4 element vector

    //Perform a sigmoid function on the result vector
    for (idx = 0; idx < ODDEMO_MATRIX_NUM_ROWS; idx++)
        prob[idx] = 1.0 / (1.0 + exp(-prob[idx]));

    //The maximum entry in the resulting vector indicates the decision index.
    //This needs to be converted into a bit-mapped value, 1 bit per zone.
    result = DSPF_sp_maxidx(prob,
                            ODDEMO_MATRIX_NUM_ROWS);

    for (idx = 0; idx < 2; idx++)
        decision[idx] = (result & (1 << idx)) >> idx;
}
