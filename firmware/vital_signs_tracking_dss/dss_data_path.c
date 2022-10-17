/**
 * dss_data_path.c
 *
 * Implements Data path processing functionality
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

/**************************************************************************
 *************************** Include Files ********************************
 **************************************************************************/

/* Standard Include Files. */
#include <stdint.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

/* BIOS/XDC Include Files. */
#include <xdc/std.h>
#include <xdc/cfg/global.h>
#include <xdc/runtime/IHeap.h>
#include <xdc/runtime/System.h>
#include <xdc/runtime/Error.h>
#include <xdc/runtime/Memory.h>
#include <ti/sysbios/BIOS.h>
#include <ti/sysbios/knl/Task.h>
#include <ti/sysbios/knl/Event.h>
#include <ti/sysbios/knl/Semaphore.h>
#include <ti/sysbios/knl/Clock.h>
#include <ti/sysbios/heaps/HeapBuf.h>
#include <ti/sysbios/heaps/HeapMem.h>
#include <ti/sysbios/knl/Event.h>
#if defined (SUBSYS_DSS)
#include <ti/sysbios/family/c64p/Hwi.h>
#include <ti/sysbios/family/c64p/EventCombiner.h>
#endif
#define DebugP_ASSERT_ENABLED 1
#include <ti/drivers/osal/DebugP.h>
#include <assert.h>
#include <ti/common/sys_common.h>
#include <ti/drivers/osal/SemaphoreP.h>
#include <ti/drivers/edma/edma.h>
#include <ti/drivers/esm/esm.h>
#include <ti/drivers/soc/soc.h>
#include <ti/utils/cycleprofiler/cycle_profiler.h>

#include <ti/alg/mmwavelib/mmwavelib.h>

/* C64P dsplib (fixed point part for C674X) */
#include "DSP_fft32x32.h"
#include "DSP_fft16x16.h"

/* C674x mathlib */
/* Suppress the mathlib.h warnings
 *  #48-D: incompatible redefinition of macro "TRUE"
 *  #48-D: incompatible redefinition of macro "FALSE"
 */
#pragma diag_push
#pragma diag_suppress 48
#include <ti/mathlib/mathlib.h>
#pragma diag_pop

// Occupancy_Detection
#include <./common/mmw_config.h>
#include "dss_mmw.h"
#include "dss_data_path.h"
#include "./common/oddemo_common.h"
#include "oddemo_heatmap.h"
#include "oddemo_feature.h"
#include "oddemo_decision.h"
#include "dss_config_edma_util.h"
#include "dss_resources.h"
//#include <ti/demo/utils/rx_ch_bias_measure.h>

// Vital Signs
#include "dss_vitalSignsDemo_utilsFunc.h"

extern ODDemo_DataPathObj oddemo_dataPathObj;
//extern VitalSigns_DataPathObj vitalsigns_dataPathObj;
extern VitalSigns_DataPathObj vitalsigns_dataPathObj_zone_1,
        vitalsigns_dataPathObj_zone_2, vitalsigns_dataPathObj_zone_3, vitalsigns_dataPathObj_zone_4;
extern MmwDemo_DSS_MCB gMmwDssMCB;

extern uint16_t oddemo_num_zones;
extern uint16_t oddemo_zone_pairs;
extern ODDEMO_Feature oddemo_feature[ODDEMO_ZONE_PAIR];
extern int8_t oddemo_decision[ODDEMO_MAX_ZONES];
extern ODDEMO_Zone oddemo_zone[ODDEMO_MAX_ZONES];
extern VS_Feature vital_signs;
extern ODDEMO_Parms oddemo_parms;
extern uint8_t peak_positions[9];

/** @brief Lookup table for Twiddle table generation and DFT single bin DFT calculation.
 * It contains 256 complex exponentials e(k) = cos(2*pi*k/1024)+j*sin(2*pi*k/1024), k=0,...,255,
 * Imaginary values are in even, and real values are in odd locations. Values are in Q31 format,
 * saturated to +/- 2147483647
 *
 * Following Matlab code was used to generate this table:
 *
 * %Generates lookup table for fast twiddle table generation for DSP lib FFT
 * %functions 16x16 and 32x32 for FFT sizes up to 1024. It saturates the
 * %amplitude to +/- 2147483647. The values are saved to file as imaginary
 * %(sine) to even, and real (cosine) to odd index positions.
 * M = 2147483647.5;
 * nMax = 1024;
 * table=round(M*exp(1j*2*pi/1024*[0:1*nMax/4-1]'));
 * table=max(min(real(table),2147483647),-2147483647) + 1j* max(min(imag(table),2147483647),-2147483647);
 * fid = fopen('twiddle_table_all.dat','w');
 *
 * tableRI = [imag(table) real(table)]';
 * tableRI = tableRI(:);
 * tableRI(tableRI<0) = tableRI(tableRI<0) + 4294967296;
 *
 * fprintf(fid, [repmat(' 0x%08x,', 1, 8) '\n'], tableRI);
 * fclose(fid);
 *
 **/
#pragma DATA_ALIGN(twiddleTableCommon, 8)
const int32_t twiddleTableCommon[2 * 256] = { 0x00000000, 0x7fffffff,
                                              0x00c90f88, 0x7fff6216,
                                              0x01921d20, 0x7ffd885a,
                                              0x025b26d7, 0x7ffa72d1,
                                              0x03242abf, 0x7ff62182,
                                              0x03ed26e6, 0x7ff09477,
                                              0x04b6195d, 0x7fe9cbbf,
                                              0x057f0035, 0x7fe1c76b,
                                              0x0647d97c, 0x7fd8878d,
                                              0x0710a345, 0x7fce0c3e,
                                              0x07d95b9e, 0x7fc25596,
                                              0x08a2009a, 0x7fb563b2,
                                              0x096a9049, 0x7fa736b4,
                                              0x0a3308bd, 0x7f97cebc,
                                              0x0afb6805, 0x7f872bf2,
                                              0x0bc3ac35, 0x7f754e7f,
                                              0x0c8bd35e, 0x7f62368f,
                                              0x0d53db92, 0x7f4de450,
                                              0x0e1bc2e4, 0x7f3857f5,
                                              0x0ee38766, 0x7f2191b3,
                                              0x0fab272b, 0x7f0991c3,
                                              0x1072a048, 0x7ef0585f,
                                              0x1139f0cf, 0x7ed5e5c6,
                                              0x120116d5, 0x7eba3a39,
                                              0x12c8106e, 0x7e9d55fc,
                                              0x138edbb1, 0x7e7f3956,
                                              0x145576b1, 0x7e5fe493,
                                              0x151bdf86, 0x7e3f57fe,
                                              0x15e21444, 0x7e1d93e9,
                                              0x16a81305, 0x7dfa98a7,
                                              0x176dd9de, 0x7dd6668e,
                                              0x183366e9, 0x7db0fdf7,
                                              0x18f8b83c, 0x7d8a5f3f,
                                              0x19bdcbf3, 0x7d628ac5,
                                              0x1a82a026, 0x7d3980ec,
                                              0x1b4732ef, 0x7d0f4218,
                                              0x1c0b826a, 0x7ce3ceb1,
                                              0x1ccf8cb3, 0x7cb72724,
                                              0x1d934fe5, 0x7c894bdd,
                                              0x1e56ca1e, 0x7c5a3d4f,
                                              0x1f19f97b, 0x7c29fbee,
                                              0x1fdcdc1b, 0x7bf88830,
                                              0x209f701c, 0x7bc5e28f,
                                              0x2161b3a0, 0x7b920b89,
                                              0x2223a4c5, 0x7b5d039d,
                                              0x22e541af, 0x7b26cb4f,
                                              0x23a6887e, 0x7aef6323,
                                              0x24677757, 0x7ab6cba3,
                                              0x25280c5e, 0x7a7d055b,
                                              0x25e845b6, 0x7a4210d8,
                                              0x26a82186, 0x7a05eead,
                                              0x27679df4, 0x79c89f6d,
                                              0x2826b928, 0x798a23b1,
                                              0x28e5714b, 0x794a7c11,
                                              0x29a3c485, 0x7909a92c,
                                              0x2a61b101, 0x78c7aba1,
                                              0x2b1f34eb, 0x78848413,
                                              0x2bdc4e6f, 0x78403328,
                                              0x2c98fbba, 0x77fab988,
                                              0x2d553afb, 0x77b417df,
                                              0x2e110a62, 0x776c4edb,
                                              0x2ecc681e, 0x77235f2d,
                                              0x2f875262, 0x76d94988,
                                              0x3041c760, 0x768e0ea5,
                                              0x30fbc54d, 0x7641af3c,
                                              0x31b54a5d, 0x75f42c0a,
                                              0x326e54c7, 0x75a585cf,
                                              0x3326e2c2, 0x7555bd4b,
                                              0x33def287, 0x7504d345,
                                              0x34968250, 0x74b2c883,
                                              0x354d9057, 0x745f9dd1,
                                              0x36041ad9, 0x740b53fa,
                                              0x36ba2014, 0x73b5ebd0,
                                              0x376f9e46, 0x735f6626,
                                              0x382493b0, 0x7307c3d0,
                                              0x38d8fe93, 0x72af05a6,
                                              0x398cdd32, 0x72552c84,
                                              0x3a402dd2, 0x71fa3948,
                                              0x3af2eeb7, 0x719e2cd2,
                                              0x3ba51e29, 0x71410804,
                                              0x3c56ba70, 0x70e2cbc6,
                                              0x3d07c1d6, 0x708378fe,
                                              0x3db832a6, 0x70231099,
                                              0x3e680b2c, 0x6fc19385,
                                              0x3f1749b8, 0x6f5f02b1,
                                              0x3fc5ec98, 0x6efb5f12,
                                              0x4073f21d, 0x6e96a99c,
                                              0x4121589a, 0x6e30e349,
                                              0x41ce1e64, 0x6dca0d14,
                                              0x427a41d0, 0x6d6227fa,
                                              0x4325c135, 0x6cf934fb,
                                              0x43d09aec, 0x6c8f351c,
                                              0x447acd50, 0x6c242960,
                                              0x452456bd, 0x6bb812d1,
                                              0x45cd358f, 0x6b4af278,
                                              0x46756828, 0x6adcc964,
                                              0x471cece6, 0x6a6d98a4,
                                              0x47c3c22f, 0x69fd614a,
                                              0x4869e665, 0x698c246c,
                                              0x490f57ee, 0x6919e320,
                                              0x49b41533, 0x68a69e81,
                                              0x4a581c9d, 0x683257ab,
                                              0x4afb6c98, 0x67bd0fbc,
                                              0x4b9e038f, 0x6746c7d7,
                                              0x4c3fdff3, 0x66cf811f,
                                              0x4ce10034, 0x66573cbb,
                                              0x4d8162c4, 0x65ddfbd3,
                                              0x4e210617, 0x6563bf92,
                                              0x4ebfe8a4, 0x64e88926,
                                              0x4f5e08e3, 0x646c59bf,
                                              0x4ffb654d, 0x63ef328f,
                                              0x5097fc5e, 0x637114cc,
                                              0x5133cc94, 0x62f201ac,
                                              0x51ced46e, 0x6271fa69,
                                              0x5269126e, 0x61f1003f,
                                              0x53028517, 0x616f146b,
                                              0x539b2aef, 0x60ec3830,
                                              0x5433027d, 0x60686cce,
                                              0x54ca0a4a, 0x5fe3b38d,
                                              0x556040e2, 0x5f5e0db3,
                                              0x55f5a4d2, 0x5ed77c89,
                                              0x568a34a9, 0x5e50015d,
                                              0x571deef9, 0x5dc79d7c,
                                              0x57b0d256, 0x5d3e5236,
                                              0x5842dd54, 0x5cb420df,
                                              0x58d40e8c, 0x5c290acc,
                                              0x59646497, 0x5b9d1153,
                                              0x59f3de12, 0x5b1035cf,
                                              0x5a82799a, 0x5a82799a,
                                              0x5b1035cf, 0x59f3de12,
                                              0x5b9d1153, 0x59646497,
                                              0x5c290acc, 0x58d40e8c,
                                              0x5cb420df, 0x5842dd54,
                                              0x5d3e5236, 0x57b0d256,
                                              0x5dc79d7c, 0x571deef9,
                                              0x5e50015d, 0x568a34a9,
                                              0x5ed77c89, 0x55f5a4d2,
                                              0x5f5e0db3, 0x556040e2,
                                              0x5fe3b38d, 0x54ca0a4a,
                                              0x60686cce, 0x5433027d,
                                              0x60ec3830, 0x539b2aef,
                                              0x616f146b, 0x53028517,
                                              0x61f1003f, 0x5269126e,
                                              0x6271fa69, 0x51ced46e,
                                              0x62f201ac, 0x5133cc94,
                                              0x637114cc, 0x5097fc5e,
                                              0x63ef328f, 0x4ffb654d,
                                              0x646c59bf, 0x4f5e08e3,
                                              0x64e88926, 0x4ebfe8a4,
                                              0x6563bf92, 0x4e210617,
                                              0x65ddfbd3, 0x4d8162c4,
                                              0x66573cbb, 0x4ce10034,
                                              0x66cf811f, 0x4c3fdff3,
                                              0x6746c7d7, 0x4b9e038f,
                                              0x67bd0fbc, 0x4afb6c98,
                                              0x683257ab, 0x4a581c9d,
                                              0x68a69e81, 0x49b41533,
                                              0x6919e320, 0x490f57ee,
                                              0x698c246c, 0x4869e665,
                                              0x69fd614a, 0x47c3c22f,
                                              0x6a6d98a4, 0x471cece6,
                                              0x6adcc964, 0x46756828,
                                              0x6b4af278, 0x45cd358f,
                                              0x6bb812d1, 0x452456bd,
                                              0x6c242960, 0x447acd50,
                                              0x6c8f351c, 0x43d09aec,
                                              0x6cf934fb, 0x4325c135,
                                              0x6d6227fa, 0x427a41d0,
                                              0x6dca0d14, 0x41ce1e64,
                                              0x6e30e349, 0x4121589a,
                                              0x6e96a99c, 0x4073f21d,
                                              0x6efb5f12, 0x3fc5ec98,
                                              0x6f5f02b1, 0x3f1749b8,
                                              0x6fc19385, 0x3e680b2c,
                                              0x70231099, 0x3db832a6,
                                              0x708378fe, 0x3d07c1d6,
                                              0x70e2cbc6, 0x3c56ba70,
                                              0x71410804, 0x3ba51e29,
                                              0x719e2cd2, 0x3af2eeb7,
                                              0x71fa3948, 0x3a402dd2,
                                              0x72552c84, 0x398cdd32,
                                              0x72af05a6, 0x38d8fe93,
                                              0x7307c3d0, 0x382493b0,
                                              0x735f6626, 0x376f9e46,
                                              0x73b5ebd0, 0x36ba2014,
                                              0x740b53fa, 0x36041ad9,
                                              0x745f9dd1, 0x354d9057,
                                              0x74b2c883, 0x34968250,
                                              0x7504d345, 0x33def287,
                                              0x7555bd4b, 0x3326e2c2,
                                              0x75a585cf, 0x326e54c7,
                                              0x75f42c0a, 0x31b54a5d,
                                              0x7641af3c, 0x30fbc54d,
                                              0x768e0ea5, 0x3041c760,
                                              0x76d94988, 0x2f875262,
                                              0x77235f2d, 0x2ecc681e,
                                              0x776c4edb, 0x2e110a62,
                                              0x77b417df, 0x2d553afb,
                                              0x77fab988, 0x2c98fbba,
                                              0x78403328, 0x2bdc4e6f,
                                              0x78848413, 0x2b1f34eb,
                                              0x78c7aba1, 0x2a61b101,
                                              0x7909a92c, 0x29a3c485,
                                              0x794a7c11, 0x28e5714b,
                                              0x798a23b1, 0x2826b928,
                                              0x79c89f6d, 0x27679df4,
                                              0x7a05eead, 0x26a82186,
                                              0x7a4210d8, 0x25e845b6,
                                              0x7a7d055b, 0x25280c5e,
                                              0x7ab6cba3, 0x24677757,
                                              0x7aef6323, 0x23a6887e,
                                              0x7b26cb4f, 0x22e541af,
                                              0x7b5d039d, 0x2223a4c5,
                                              0x7b920b89, 0x2161b3a0,
                                              0x7bc5e28f, 0x209f701c,
                                              0x7bf88830, 0x1fdcdc1b,
                                              0x7c29fbee, 0x1f19f97b,
                                              0x7c5a3d4f, 0x1e56ca1e,
                                              0x7c894bdd, 0x1d934fe5,
                                              0x7cb72724, 0x1ccf8cb3,
                                              0x7ce3ceb1, 0x1c0b826a,
                                              0x7d0f4218, 0x1b4732ef,
                                              0x7d3980ec, 0x1a82a026,
                                              0x7d628ac5, 0x19bdcbf3,
                                              0x7d8a5f3f, 0x18f8b83c,
                                              0x7db0fdf7, 0x183366e9,
                                              0x7dd6668e, 0x176dd9de,
                                              0x7dfa98a7, 0x16a81305,
                                              0x7e1d93e9, 0x15e21444,
                                              0x7e3f57fe, 0x151bdf86,
                                              0x7e5fe493, 0x145576b1,
                                              0x7e7f3956, 0x138edbb1,
                                              0x7e9d55fc, 0x12c8106e,
                                              0x7eba3a39, 0x120116d5,
                                              0x7ed5e5c6, 0x1139f0cf,
                                              0x7ef0585f, 0x1072a048,
                                              0x7f0991c3, 0x0fab272b,
                                              0x7f2191b3, 0x0ee38766,
                                              0x7f3857f5, 0x0e1bc2e4,
                                              0x7f4de450, 0x0d53db92,
                                              0x7f62368f, 0x0c8bd35e,
                                              0x7f754e7f, 0x0bc3ac35,
                                              0x7f872bf2, 0x0afb6805,
                                              0x7f97cebc, 0x0a3308bd,
                                              0x7fa736b4, 0x096a9049,
                                              0x7fb563b2, 0x08a2009a,
                                              0x7fc25596, 0x07d95b9e,
                                              0x7fce0c3e, 0x0710a345,
                                              0x7fd8878d, 0x0647d97c,
                                              0x7fe1c76b, 0x057f0035,
                                              0x7fe9cbbf, 0x04b6195d,
                                              0x7ff09477, 0x03ed26e6,
                                              0x7ff62182, 0x03242abf,
                                              0x7ffa72d1, 0x025b26d7,
                                              0x7ffd885a, 0x01921d20,
                                              0x7fff6216, 0x00c90f88 };

/**
 *  @b Description
 *  @n
 *      Following function is equivalent of the dsplib's gen_twiddle_fft32x32() function
 *      optimized for speed to allow quick reconfiguration when switching sub-frames
 *      in advanced frame mode. The alternative is to store tables for each sub-frame
 *      which is very high in memory consumption. The maximum error with respect to the
 *      dsplib function is in the LSB (+/- 1).
 */
int32_t MmwDemo_gen_twiddle_fft32x32_fast(int32_t *w, int32_t n, double scale)
{
    int32_t i, j, k;
    int32_t log2n = 30 - _norm(n); //Note n is always power of 2
    int32_t step = 1024 >> log2n;
    int32_t step6 = 3 * step;
    int32_t step4 = 2 * step;
    int32_t step2 = 1 * step;
    int32_t ind, indLsb, indMsb;
    int64_t *restrict table = (int64_t*) twiddleTableCommon;
    int64_t *restrict wd = (int64_t*) w;
    int32_t xRe;
    int32_t xIm;

    for (j = 1, k = 0; j < n >> 2; j = j << 2)
    {
        for (i = 0; i < n >> 2; i += j)
        {
            ind = step2 * (i);
            indLsb = ind & 0xFF;
            indMsb = (ind >> 8) & 0x3;
            xRe = _hill(table[indLsb]);
            xIm = _loll(table[indLsb]);
            if (indMsb == 0)
            {
                wd[k + 0] = _itoll(xRe, xIm);
            }
            if (indMsb == 1)
            {
                wd[k + 0] = _itoll(-xIm, xRe);
            }
            if (indMsb == 2)
            {
                wd[k + 0] = _itoll(-xRe, -xIm);
            }

            ind = step4 * (i);
            indLsb = ind & 0xFF;
            indMsb = (ind >> 8) & 0x3;
            xRe = _hill(table[indLsb]);
            xIm = _loll(table[indLsb]);
            if (indMsb == 0)
            {
                wd[k + 1] = _itoll(xRe, xIm);
            }
            if (indMsb == 1)
            {
                wd[k + 1] = _itoll(-xIm, xRe);
            }
            if (indMsb == 2)
            {
                wd[k + 1] = _itoll(-xRe, -xIm);
            }

            ind = step6 * (i);
            indLsb = ind & 0xFF;
            indMsb = (ind >> 8) & 0x3;
            xRe = _hill(table[indLsb]);
            xIm = _loll(table[indLsb]);
            if (indMsb == 0)
            {
                wd[k + 2] = _itoll(xRe, xIm);
            }
            if (indMsb == 1)
            {
                wd[k + 2] = _itoll(-xIm, xRe);
            }
            if (indMsb == 2)
            {
                wd[k + 2] = _itoll(-xRe, -xIm);
            }

            k += 3;
        }
    }
    return 2 * k;
}

/**
 *  @b Description
 *  @n
 *      Following function is equivalent of the dsplib's gen_twiddle_fft16x16() function
 *      optimized for speed to allow quick reconfiguration when switching sub-frames
 *      in advanced frame mode. The alternative is to store tables for each sub-frame
 *      which is very high in memory consumption. The maximum error with respect to the
 *     dsplib function is in the LSB (+/- 1).
 */
int32_t MmwDemo_gen_twiddle_fft16x16_fast(short *w, int32_t n)
{
    int32_t i, j, k;
    int32_t log2n = 30 - _norm(n); //Note n is always power of 2
    int32_t step = 1024 >> log2n;
    int32_t step6 = 3 * step;
    int32_t step4 = 2 * step;
    int32_t step2 = 1 * step;
    int32_t ind, indLsb, indMsb;
    int64_t *restrict table = (int64_t*) twiddleTableCommon;
    uint32_t *restrict wd = (uint32_t*) w;
    int32_t xRe;
    int32_t xIm;

    for (j = 1, k = 0; j < n >> 2; j = j << 2)
    {
        for (i = 0; i < n >> 2; i += j << 1)
        {
            ind = step2 * (i + j);
            indLsb = ind & 0xFF;
            indMsb = (ind >> 8) & 0x3;
            xRe = ((int32_t) _sadd(_hill(table[indLsb]), 0x00008000)) >> 16;
            xIm = ((int32_t) _sadd(_loll(table[indLsb]), 0x00008000)) >> 16;
            if (indMsb == 0)
            {
                wd[k + 1] = _pack2(xRe, xIm);
            }
            if (indMsb == 1)
            {
                wd[k + 1] = _pack2(-xIm, xRe);
            }
            if (indMsb == 2)
            {
                wd[k + 1] = _pack2(-xRe, -xIm);
            }

            ind = step2 * (i);
            indLsb = ind & 0xFF;
            indMsb = (ind >> 8) & 0x3;
            xRe = ((int32_t) _sadd(_hill(table[indLsb]), 0x00008000)) >> 16;
            xIm = ((int32_t) _sadd(_loll(table[indLsb]), 0x00008000)) >> 16;
            if (indMsb == 0)
            {
                wd[k + 0] = _pack2(xRe, xIm);
            }
            if (indMsb == 1)
            {
                wd[k + 0] = _pack2(-xIm, xRe);
            }
            if (indMsb == 2)
            {
                wd[k + 0] = _pack2(-xRe, -xIm);
            }

            ind = step4 * (i + j);
            indLsb = ind & 0xFF;
            indMsb = (ind >> 8) & 0x3;
            xRe = ((int32_t) _sadd(_hill(table[indLsb]), 0x00008000)) >> 16;
            xIm = ((int32_t) _sadd(_loll(table[indLsb]), 0x00008000)) >> 16;
            if (indMsb == 0)
            {
                wd[k + 3] = _pack2(-xRe, -xIm);
            }
            if (indMsb == 1)
            {
                wd[k + 3] = _pack2(xIm, -xRe);
            }
            if (indMsb == 2)
            {
                wd[k + 3] = _pack2(xRe, xIm);
            }

            ind = step4 * (i);
            indLsb = ind & 0xFF;
            indMsb = (ind >> 8) & 0x3;
            xRe = ((int32_t) _sadd(_hill(table[indLsb]), 0x00008000)) >> 16;
            xIm = ((int32_t) _sadd(_loll(table[indLsb]), 0x00008000)) >> 16;
            if (indMsb == 0)
            {
                wd[k + 2] = _pack2(-xRe, -xIm);
            }
            if (indMsb == 1)
            {
                wd[k + 2] = _pack2(xIm, -xRe);
            }
            if (indMsb == 2)
            {
                wd[k + 2] = _pack2(xRe, xIm);
            }

            ind = step6 * (i + j);
            indLsb = ind & 0xFF;
            indMsb = (ind >> 8) & 0x3;
            xRe = ((int32_t) _sadd(_hill(table[indLsb]), 0x00008000)) >> 16;
            xIm = ((int32_t) _sadd(_loll(table[indLsb]), 0x00008000)) >> 16;
            if (indMsb == 0)
            {
                wd[k + 5] = _pack2(xRe, xIm);
            }
            if (indMsb == 1)
            {
                wd[k + 5] = _pack2(-xIm, xRe);
            }
            if (indMsb == 2)
            {
                wd[k + 5] = _pack2(-xRe, -xIm);
            }

            ind = step6 * (i);
            indLsb = ind & 0xFF;
            indMsb = (ind >> 8) & 0x3;
            xRe = ((int32_t) _sadd(_hill(table[indLsb]), 0x00008000)) >> 16;
            xIm = ((int32_t) _sadd(_loll(table[indLsb]), 0x00008000)) >> 16;
            if (indMsb == 0)
            {
                wd[k + 4] = _pack2(xRe, xIm);
            }
            if (indMsb == 1)
            {
                wd[k + 4] = _pack2(-xIm, xRe);
            }
            if (indMsb == 2)
            {
                wd[k + 4] = _pack2(-xRe, -xIm);
            }

            k += 6;
        }
    }
    return 2 * k;
}

/*! @brief L2 heap used for allocating buffers in L2 SRAM,
 mostly scratch buffers */
#define MMW_L2_HEAP_SIZE    0xC000U
//#define MMW_L2_HEAP_SIZE    0xA000U

/*! @brief L1 heap used for allocating buffers in L1D SRAM,
 mostly scratch buffers */
//#define MMW_L1_HEAP_SIZE    0x4000U
#define MMW_L1_HEAP_SIZE    0x3800U

/*! L3 RAM buffer */
#pragma DATA_SECTION(gMmwL3, ".l3data");
#pragma DATA_ALIGN(gMmwL3, 8);
uint8_t gMmwL3[ODDEMO_L3_SIZE];

/*! L2 Heap */
#pragma DATA_SECTION(gMmwL2, ".l2data");
#pragma DATA_ALIGN(gMmwL2, 8);
uint8_t gMmwL2[MMW_L2_HEAP_SIZE];

/*! L1 Heap */
#pragma DATA_SECTION(gMmwL1, ".l1data");
#pragma DATA_ALIGN(gMmwL1, 8);
uint8_t gMmwL1[MMW_L1_HEAP_SIZE];

// Occupancy_Detection Demo
#pragma DATA_SECTION(oddemo_steeringVec, ".l3data");
#pragma DATA_ALIGN(oddemo_steeringVec, 8);
cplxf_t oddemo_steeringVec[ODDEMO_STEERINGVEC_L1_BUF_SIZE];

#pragma DATA_SECTION(oddemo_scratchPad, ".l1data");
#pragma DATA_ALIGN(oddemo_scratchPad, 8);
uint32_t oddemo_scratchPad[ODDEMO_SCRATCH_L1_BUF_SIZE];

#pragma DATA_SECTION(oddemo_invRnMatrix, ".l1data");
#pragma DATA_ALIGN(oddemo_invRnMatrix, 8);
//cplxf_t oddemo_invRnMatrix[ODDEMO_INVRNMATRIX_BUFFER];
cplxf_t oddemo_invRnMatrix[ODDEMO_INVRNMATRIX_BUF_SIZE];

#pragma DATA_SECTION(oddemo_invRnMatrix_VS, ".l1data");
#pragma DATA_ALIGN(oddemo_invRnMatrix, 8);
cplxf_t oddemo_invRnMatrix_VS[ODDEMO_INVRNMATRIX_BUF_SIZE];

//oddemo_rangeAzimuthHeatMap is allocated
//in L3 to avoid copying to HSM Buffer by MSS when transmitted over UART.

#pragma DATA_SECTION(oddemo_rangeAzimuthHeatMap, ".l3data");
#pragma DATA_ALIGN(oddemo_rangeAzimuthHeatMap, 8);
float oddemo_rangeAzimuthHeatMap[ODDEMO_ANGLEHEATMAP_BUF_SIZE];

#pragma DATA_SECTION(oddemo_shortHeatMap, ".l3data");
#pragma DATA_ALIGN(oddemo_shortHeatMap, 8);
uint16_t oddemo_shortHeatMap[ODDEMO_ANGLEHEATMAP_BUF_SIZE];

#pragma DATA_SECTION(oddemo_scratch_output, ".l3data");
#pragma DATA_ALIGN(oddemo_scratch_output, 8);
float oddemo_scratch_output[72];

float oddemo_coeffMatrix[ODDEMO_ZONE_PAIR][ODDEMO_MATRIX_SIZE];

/*! Types of FFT window */
/*! FFT window 16 - samples format is int16_t */
#define FFT_WINDOW_INT16 0
/*! FFT window 32 - samples format is int32_t */
#define FFT_WINDOW_INT32 1

/* FFT Window */
/*! Hanning window */
#define MMW_WIN_HANNING  0
/*! Blackman window */
#define MMW_WIN_BLACKMAN 1
/*! Rectangular window */
#define MMW_WIN_RECT     2

/*! If MMW_USE_SINGLE_POINT_DFT is defined azimuth calculation uses single
 * point DFT, otherwise FFT function from DSP lib*/
#define MMW_USE_SINGLE_POINT_DFT

/* Vital Signs Params */

static uint32_t gFrameCount;

float tempFloat;

void MmwDemo_genWindow(void *win, uint32_t windowDatumType, uint32_t winLen,
                       uint32_t winGenLen, int32_t oneQformat, uint32_t winType);

#define MMW_EDMA_TRIGGER_ENABLE  1
#define MMW_EDMA_TRIGGER_DISABLE 0

extern volatile cycleLog_t gCycleLog;

/**
 *  @b Description
 *  @n
 *      Resets the Doppler line bit mask. Doppler line bit mask indicates Doppler
 *      lines (bins) on wich the CFAR in Doppler direction detected objects.
 *      After the CFAR in Doppler direction is completed for all range bins, the
 *      CFAR in range direction is performed on indicated Doppler lines.
 *      The array dopplerLineMask is uint32_t array. The LSB bit of the first word
 *      corresponds to Doppler line (bin) zero.
 *
 */
void MmwDemo_resetDopplerLines(MmwDemo_1D_DopplerLines_t *ths)
{
    memset((void*) ths->dopplerLineMask, 0,
           ths->dopplerLineMaskLen * sizeof(uint32_t));
    ths->currentIndex = 0;
}

/**
 *  @b Description
 *  @n
 *      Sets the bit in the Doppler line bit mask dopplerLineMask corresponding to Doppler
 *      line on which CFAR in Doppler direction detected object. Indicating the Doppler
 *      line being active in observed frame. @sa MmwDemo_resetDopplerLines
 */
void MmwDemo_setDopplerLine(MmwDemo_1D_DopplerLines_t *ths,
                            uint16_t dopplerIndex)
{
    uint32_t word = dopplerIndex >> 5;
    uint32_t bit = dopplerIndex & 31;

    ths->dopplerLineMask[word] |= (0x1 << bit);
}

/**
 *  @b Description
 *  @n
 *      Checks whetehr Doppler line is active in the observed frame. It checks whether the bit
 *      is set in the Doppler line bit mask corresponding to Doppler
 *      line on which CFAR in Doppler direction detected object.
 *      @sa MmwDemo_resetDopplerLines
 */
uint32_t MmwDemo_isSetDopplerLine(MmwDemo_1D_DopplerLines_t *ths,
                                  uint16_t index)
{
    uint32_t dopplerLineStat;
    uint32_t word = index >> 5;
    uint32_t bit = index & 31;

    if (ths->dopplerLineMask[word] & (0x1 << bit))
    {
        dopplerLineStat = 1;
    }
    else
    {
        dopplerLineStat = 0;
    }
    return dopplerLineStat;
}

/**
 *  @b Description
 *  @n
 *      Gets the Doppler index from the Doppler line bit mask, starting from the
 *      smallest active Doppler lin (bin). Subsequent calls return the next
 *      active Doppler line. @sa MmwDemo_resetDopplerLines
 *
 */
int32_t MmwDemo_getDopplerLine(MmwDemo_1D_DopplerLines_t *ths)
{
    uint32_t index = ths->currentIndex;
    uint32_t word = index >> 5;
    uint32_t bit = index & 31;

    while (((ths->dopplerLineMask[word] >> bit) & 0x1) == 0)
    {
        index++;
        bit++;
        if (bit == 32)
        {
            word++;
            bit = 0;
            if (word >= ths->dopplerLineMaskLen)
            {
                MmwDemo_dssAssert(0);
            }
        }
    }
    ths->currentIndex = index + 1;
    return index;
}

/**
 *  @b Description
 *  @n
 *      Power of 2 round up function.
 */
uint32_t MmwDemo_pow2roundup(uint32_t x)
{
    uint32_t result = 1;
    while (x > result)
    {
        result <<= 1;
    }
    return result;
}

/**
 *  @b Description
 *  @n
 *      Waits for 1D FFT data to be transferrd to input buffer.
 *      This is a blocking function.
 *
 *  @param[in] obj  Pointer to data path object
 *  @param[in] pingPongId ping-pong id (ping is 0 and pong is 1)
 *
 *  @retval
 *      NONE
 */
void MmwDemo_dataPathWait1DInputData(MmwDemo_DSS_DataPathObj *obj,
                                     uint32_t pingPongId)
{
    MmwDemo_DSS_dataPathContext_t *context = obj->context;

#ifdef EDMA_1D_INPUT_BLOCKING
    Bool       status;

    status = Semaphore_pend(context->EDMA_1D_InputDone_semHandle[pingPongId], BIOS_WAIT_FOREVER);
    if (status != TRUE)
    {
        System_printf("Error: Semaphore_pend returned %d\n",status);
    }
#else
    /* wait until transfer done */
    volatile bool isTransferDone;
    uint8_t chId;
    if (pingPongId == 0)
    {
        chId = MMW_EDMA_CH_1D_IN_PING;
    }
    else
    {
        chId = MMW_EDMA_CH_1D_IN_PONG;
    }
    do
    {
        if (EDMA_isTransferComplete(
                context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE], chId,
                (bool*) &isTransferDone) != EDMA_NO_ERROR)
        {
            MmwDemo_dssAssert(0);
        }
    }
    while (isTransferDone == false);
#endif
}

/**
 *  @b Description
 *  @n
 *      Waits for 1D FFT data to be transferred to output buffer.
 *      This is a blocking function.
 *
 *  @param[in] obj  Pointer to data path object
 *  @param[in] pingPongId ping-pong id (ping is 0 and pong is 1)
 *
 *  @retval
 *      NONE
 */
void MmwDemo_dataPathWait1DOutputData(MmwDemo_DSS_DataPathObj *obj,
                                      uint32_t pingPongId)
{
    MmwDemo_DSS_dataPathContext_t *context = obj->context;

#ifdef EDMA_1D_OUTPUT_BLOCKING
    Bool       status;

    status = Semaphore_pend(context->EDMA_1D_OutputDone_semHandle[pingPongId], BIOS_WAIT_FOREVER);
    if (status != TRUE)
    {
        System_printf("Error: Semaphore_pend returned %d\n",status);
    }
#else
    /* wait until transfer done */
    volatile bool isTransferDone;
    uint8_t chId;
    if (pingPongId == 0)
    {
        chId = MMW_EDMA_CH_1D_OUT_PING;
    }
    else
    {
        chId = MMW_EDMA_CH_1D_OUT_PONG;
    }
    do
    {
        if (EDMA_isTransferComplete(
                context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE], chId,
                (bool*) &isTransferDone) != EDMA_NO_ERROR)
        {
            MmwDemo_dssAssert(0);
        }
    }
    while (isTransferDone == false);
#endif
}

void MmwDemo_EDMA_CQTransferCompletionCallbackFxn(
        uintptr_t arg, uint8_t transferCompletionCode)
{
    MmwDemo_DSS_DataPathObj *obj = (MmwDemo_DSS_DataPathObj*) arg;
    uint8_t chirpIdx;

    switch (transferCompletionCode)
    {
    /* Validate numSlices from data buffer for every chirp.
     Note1: If the configured number of slices covers more time period than the sampling time,
     the actual reported number of slices maybe smaller.

     Note 2: If there are extensive CQ data processing,
     then the processing code should be moved to a task
     */

    case MMW_EDMA_CH_SIGIMG_MON:
        obj->datapathCQ.sigImgEdmaCnt++;
        for (chirpIdx = 0; chirpIdx < obj->numChirpsPerChirpEvent; chirpIdx++)
        {
            uint8_t *sigImgData = obj->datapathCQ.sigImgData;

            if (sigImgData[chirpIdx * obj->datapathCQ.sigImgMonDataSizePerChirp]
                    > obj->datapathCQ.sigImgMonCfg->numSlices)
            {
                obj->datapathCQ.sigImgErrCnt++;
                MmwDemo_dssAssert(0);
            }
        }
        break;

    case MMW_EDMA_CH_RX_SATURATION_MON:
        obj->datapathCQ.rxSatEdmaCnt++;
        for (chirpIdx = 0; chirpIdx < obj->numChirpsPerChirpEvent; chirpIdx++)
        {
            uint8_t *rxSatData = obj->datapathCQ.rxSatData;
            if (rxSatData[chirpIdx * obj->datapathCQ.satMonDataSizePerChirp]
                    > obj->datapathCQ.rxSatMonCfg->numSlices)
            {
                obj->datapathCQ.rxSatErrCnt++;
                MmwDemo_dssAssert(0);
            }
        }
        break;

    default:
        MmwDemo_dssAssert(0)
        ;
        break;
    }
}

/**
 *  @b Description
 *  @n
 *      Configures all EDMA channels and param sets used in data path processing
 *  @param[in] obj  Pointer to data path object
 *
 *  @retval
 *      -1 if error, 0 for no error
 */
int32_t MmwDemo_dataPathConfigCQEdma(MmwDemo_DSS_DataPathObj *obj)
{
    uint32_t eventQueue;
    int32_t retVal = 0;
    MmwDemo_DSS_dataPathContext_t *context = obj->context;

    /*****************************************************
     * EDMA configuration for getting CQ data from CQ buffer
     * to Datapath CQ storage
     *****************************************************/
    eventQueue = 0U;
    if (obj->datapathCQ.anaMonCfg->sigImgMonEn)
    {
        retVal = EDMAutil_configType1(
                context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
                (uint8_t*) (SOC_translateAddress(
                        (uint32_t) obj->datapathCQ.sigImgMonAddr,
                        SOC_TranslateAddr_Dir_TO_EDMA, NULL)),
                (uint8_t*) (SOC_translateAddress(
                        (uint32_t) obj->datapathCQ.sigImgData,
                        SOC_TranslateAddr_Dir_TO_EDMA, NULL)),
                MMW_EDMA_CH_SIGIMG_MON,
                false,
                MMW_EDMA_CH_SIGIMG_MON, obj->datapathCQ.sigImgMonTotalSize, 1,
                0, 0, eventQueue, MmwDemo_EDMA_CQTransferCompletionCallbackFxn,
                (uintptr_t) obj);
        if (retVal < 0)
        {
            return -1;
        }
    }

    if (obj->datapathCQ.anaMonCfg->rxSatMonEn)
    {
        retVal = EDMAutil_configType1(
                context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
                (uint8_t*) (SOC_translateAddress(
                        (uint32_t) obj->datapathCQ.satMonAddr,
                        SOC_TranslateAddr_Dir_TO_EDMA, NULL)),
                (uint8_t*) (SOC_translateAddress(
                        (uint32_t) obj->datapathCQ.rxSatData,
                        SOC_TranslateAddr_Dir_TO_EDMA, NULL)),
                MMW_EDMA_CH_RX_SATURATION_MON,
                false,
                MMW_EDMA_CH_RX_SATURATION_MON, obj->datapathCQ.satMonTotalSize,
                1, 0, 0, eventQueue,
                MmwDemo_EDMA_CQTransferCompletionCallbackFxn, (uintptr_t) obj);
        if (retVal < 0)
        {
            return -1;
        }
    }
    return 0;
}

/**
 *  @b Description
 *  @n
 *      Function to configure CQ.
 *  @param[in] ptrDataPathObj Pointer to data path object.
 *
 *  @retval
 *      0 if no error, else error (there will be system prints for these).
 */
void MmwDemo_dssDataPathStartCQEdma(MmwDemo_DSS_DataPathObj *ptrDataPathObj)
{
    MmwDemo_DSS_dataPathContext_t *context = ptrDataPathObj->context;

    /* Manually start EDMA transfer */
    if (ptrDataPathObj->datapathCQ.anaMonCfg->sigImgMonEn)
    {
        /* Kick off DMA to fetch data from CQ1 buffer */
        EDMA_startDmaTransfer(context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
        MMW_EDMA_CH_SIGIMG_MON);
    }

    if (ptrDataPathObj->datapathCQ.anaMonCfg->rxSatMonEn)
    {
        /* Kick off DMA to fetch data from CQ2 buffer */
        EDMA_startDmaTransfer(context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
        MMW_EDMA_CH_RX_SATURATION_MON);
    }
}

/**
 *  @b Description
 *  @n
 *      Configures all EDMA channels and param sets used in data path processing
 *  @param[in] obj  Pointer to data path object
 *
 *  @retval
 *      -1 if error, 0 for no error
 */
int32_t MmwDemo_dataPathConfigEdma(MmwDemo_DSS_DataPathObj *obj)
{
    uint32_t eventQueue;
    uint16_t shadowParam = EDMA_NUM_DMA_CHANNELS;
    int32_t retVal = 0;
    MmwDemo_DSS_dataPathContext_t *context = obj->context;
    //MmwDemo_AnaMonitorCfg*      ptrAnaMonitorCfg;

    /*****************************************************
     * EDMA configuration for getting ADC data from ADC buffer
     * to L2 (prior to 1D FFT)
     * For ADC Buffer to L2 use EDMA-A TPTC =1
     *****************************************************/
    eventQueue = 0U;
    /* Ping - copies chirp samples from even antenna numbers (e.g. RxAnt0 and RxAnt2) */
    retVal = EDMAutil_configType1(
            context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
            (uint8_t*) (&obj->ADCdataBuf[0]),
            (uint8_t*) (SOC_translateAddress((uint32_t) &obj->adcDataIn[0],
                                             SOC_TranslateAddr_Dir_TO_EDMA,
                                             NULL)),
            MMW_EDMA_CH_1D_IN_PING,
            false,
            shadowParam++,
            obj->numAdcSamples * BYTES_PER_SAMP_1D,
            MAX(obj->numRxAntennas / 2, 1) * obj->numChirpsPerChirpEvent,
            (obj->numAdcSamples * BYTES_PER_SAMP_1D * 2)
                    * obj->numChirpsPerChirpEvent,
            0, eventQueue,
#ifdef EDMA_1D_INPUT_BLOCKING
            MmwDemo_EDMA_transferCompletionCallbackFxn,
#else
            NULL,
#endif
            (uintptr_t) obj);
    if (retVal < 0)
    {
        return -1;
    }

    /* Pong - copies chirp samples from odd antenna numbers (e.g. RxAnt1 and RxAnt3) */
    retVal = EDMAutil_configType1(
            context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
            (uint8_t*) (&obj->ADCdataBuf[obj->numAdcSamples
                    * obj->numChirpsPerChirpEvent]),
            (uint8_t*) (SOC_translateAddress(
                    (uint32_t) (&obj->adcDataIn[obj->numRangeBins]),
                    SOC_TranslateAddr_Dir_TO_EDMA, NULL)),
            MMW_EDMA_CH_1D_IN_PONG,
            false,
            shadowParam++,
            obj->numAdcSamples * BYTES_PER_SAMP_1D,
            MAX(obj->numRxAntennas / 2, 1) * obj->numChirpsPerChirpEvent,
            (obj->numAdcSamples * BYTES_PER_SAMP_1D * 2)
                    * obj->numChirpsPerChirpEvent,
            0, eventQueue,
#ifdef EDMA_1D_INPUT_BLOCKING
            MmwDemo_EDMA_transferCompletionCallbackFxn,
#else
            NULL,
#endif
            (uintptr_t) obj);
    if (retVal < 0)
    {
        return -1;
    }

    //OD DEMO: Replaced 1.2.0.5 EDMAs with 1.0.0.5 version to keep heatmap format
    eventQueue = 1U;
    /*****************************************************
     * EDMA configuration for storing 1d fft output in transposed manner to L3.
     * It copies all Rx antennas of the chirp per trigger event.
     *****************************************************/
    /* Ping - Copies from ping FFT output (even chirp indices)  to L3 */
    retVal = EDMAutil_configType2a(
            context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
            (uint8_t*) (SOC_translateAddress((uint32_t) (&obj->fftOut1D[0]),
                                             SOC_TranslateAddr_Dir_TO_EDMA,
                                             NULL)),
            (uint8_t*) (&obj->radarCube[0]),
            MMW_EDMA_CH_1D_OUT_PING,
            false, shadowParam++,
            BYTES_PER_SAMP_1D,
            obj->numRangeBins, obj->numTxAntennas, obj->numRxAntennas,
            obj->numDopplerBins, eventQueue,
            NULL,
            (uintptr_t) obj);
    if (retVal < 0)
    {
        return -1;
    }

    /* Ping - Copies from pong FFT output (odd chirp indices)  to L3 */
    retVal = EDMAutil_configType2a(
            context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
            (uint8_t*) (SOC_translateAddress(
                    (uint32_t) (&obj->fftOut1D[obj->numRxAntennas
                            * obj->numRangeBins]),
                    SOC_TranslateAddr_Dir_TO_EDMA, NULL)),
            (uint8_t*) (&obj->radarCube[0]),
            MMW_EDMA_CH_1D_OUT_PONG,
            false, shadowParam++,
            BYTES_PER_SAMP_1D,
            obj->numRangeBins, obj->numTxAntennas, obj->numRxAntennas,
            obj->numDopplerBins, eventQueue,
            NULL,
            (uintptr_t) obj);
    if (retVal < 0)
    {
        return -1;
    }

    return (0);
}

/**
 *  @b Description
 *  @n
 *    Outputs magnitude squared float array of input complex32 array
 *
 *  @retval
 *      Not Applicable.
 */
void MmwDemo_magnitudeSquared(cmplx32ReIm_t *restrict inpBuff,
                              float *restrict magSqrdBuff, uint32_t numSamples)
{
    uint32_t i;
    for (i = 0; i < numSamples; i++)
    {
        magSqrdBuff[i] = (float) inpBuff[i].real * (float) inpBuff[i].real
                + (float) inpBuff[i].imag * (float) inpBuff[i].imag;
    }
}

#define pingPongId(x) ((x) & 0x1U)
#define isPong(x) (pingPongId(x))

/**
 *  @b Description
 *  @n
 *    Compensation of DC range antenna signature
 *
 *
 *  @retval
 *      Not Applicable.
 */
void MmwDemo_dcRangeSignatureCompensation(MmwDemo_DSS_DataPathObj *obj,
                                          uint8_t chirpPingPongId)
{
    uint32_t rxAntIdx, binIdx;
    uint32_t ind;
    int32_t chirpPingPongOffs;
    int32_t chirpPingPongSize;
    MmwDemo_CalibDcRangeSigCfg *calibDcCfg = &obj->cliCfg->calibDcRangeSigCfg;

    chirpPingPongSize = obj->numRxAntennas
            * (calibDcCfg->positiveBinIdx - calibDcCfg->negativeBinIdx + 1);
    if (obj->dcRangeSigCalibCntr == 0)
    {
        memset(obj->dcRangeSigMean, 0,
               obj->numTxAntennas * chirpPingPongSize * sizeof(cmplx32ImRe_t));
    }

    chirpPingPongOffs = chirpPingPongId * chirpPingPongSize;

    /* Calibration */
    if (obj->dcRangeSigCalibCntr
            < (calibDcCfg->numAvgChirps * obj->numTxAntennas))
    {
        /* Accumulate */
        ind = 0;
        for (rxAntIdx = 0; rxAntIdx < obj->numRxAntennas; rxAntIdx++)
        {
            uint32_t chirpInOffs = chirpPingPongId
                    * (obj->numRxAntennas * obj->numRangeBins)
                    + (obj->numRangeBins * rxAntIdx);
            int64_t *meanPtr =
                    (int64_t*) &obj->dcRangeSigMean[chirpPingPongOffs];
            uint32_t *fftPtr = (uint32_t*) &obj->fftOut1D[chirpInOffs];
            int64_t meanBin;
            uint32_t fftBin;
            int32_t Re, Im;
            for (binIdx = 0; binIdx <= calibDcCfg->positiveBinIdx; binIdx++)
            {
                meanBin = _amem8(&meanPtr[ind]);
                fftBin = _amem4(&fftPtr[binIdx]);
                Im = _loll(meanBin) + _ext(fftBin, 0, 16);
                Re = _hill(meanBin) + _ext(fftBin, 16, 16);
                _amem8(&meanPtr[ind]) = _itoll(Re, Im);
                ind++;
            }

            chirpInOffs = chirpPingPongId
                    * (obj->numRxAntennas * obj->numRangeBins)
                    + (obj->numRangeBins * rxAntIdx) + obj->numRangeBins
                    + calibDcCfg->negativeBinIdx;
            fftPtr = (uint32_t*) &obj->fftOut1D[chirpInOffs];
            for (binIdx = 0; binIdx < -calibDcCfg->negativeBinIdx; binIdx++)
            {
                meanBin = _amem8(&meanPtr[ind]);
                fftBin = _amem4(&fftPtr[binIdx]);
                Im = _loll(meanBin) + _ext(fftBin, 0, 16);
                Re = _hill(meanBin) + _ext(fftBin, 16, 16);
                _amem8(&meanPtr[ind]) = _itoll(Re, Im);
                ind++;
            }
        }
        obj->dcRangeSigCalibCntr++;

        if (obj->dcRangeSigCalibCntr
                == (calibDcCfg->numAvgChirps * obj->numTxAntennas))
        {
            /* Divide */
            int64_t *meanPtr = (int64_t*) obj->dcRangeSigMean;
            int32_t Re, Im;
            int64_t meanBin;
            int32_t divShift = obj->log2NumAvgChirps;
            for (ind = 0; ind < (obj->numTxAntennas * chirpPingPongSize); ind++)
            {
                meanBin = _amem8(&meanPtr[ind]);
                Im = _sshvr(_loll(meanBin), divShift);
                Re = _sshvr(_hill(meanBin), divShift);
                _amem8(&meanPtr[ind]) = _itoll(Re, Im);
            }
        }
    }
    else
    {
        /* fftOut1D -= dcRangeSigMean */
        ind = 0;
        for (rxAntIdx = 0; rxAntIdx < obj->numRxAntennas; rxAntIdx++)
        {
            uint32_t chirpInOffs = chirpPingPongId
                    * (obj->numRxAntennas * obj->numRangeBins)
                    + (obj->numRangeBins * rxAntIdx);
            int64_t *meanPtr =
                    (int64_t*) &obj->dcRangeSigMean[chirpPingPongOffs];
            uint32_t *fftPtr = (uint32_t*) &obj->fftOut1D[chirpInOffs];
            int64_t meanBin;
            uint32_t fftBin;
            int32_t Re, Im;
            for (binIdx = 0; binIdx < calibDcCfg->positiveBinIdx; binIdx++)
            {
                meanBin = _amem8(&meanPtr[ind]);
                fftBin = _amem4(&fftPtr[binIdx]);
                Im = _ext(fftBin, 0, 16) - _loll(meanBin);
                Re = _ext(fftBin, 16, 16) - _hill(meanBin);
                _amem4(&fftPtr[binIdx]) = _pack2(Im, Re);
                ind++;
            }

            chirpInOffs = chirpPingPongId
                    * (obj->numRxAntennas * obj->numRangeBins)
                    + (obj->numRangeBins * rxAntIdx) + obj->numRangeBins
                    + calibDcCfg->negativeBinIdx;
            fftPtr = (uint32_t*) &obj->fftOut1D[chirpInOffs];
            for (binIdx = 0; binIdx < -calibDcCfg->negativeBinIdx; binIdx++)
            {
                meanBin = _amem8(&meanPtr[ind]);
                fftBin = _amem4(&fftPtr[binIdx]);
                Im = _ext(fftBin, 0, 16) - _loll(meanBin);
                Re = _ext(fftBin, 16, 16) - _hill(meanBin);
                _amem4(&fftPtr[binIdx]) = _pack2(Im, Re);
                //_amem4(&fftPtr[binIdx]) = _packh2(_sshvl(Im,16) , _sshvl(Re, 16));
                ind++;
            }
        }
    }
}

void VS_Processing_Chain(VitalSigns_DataPathObj *obj_VS, uint8_t zone_number)
{

    float doppler_input_real, doppler_input_imag;
    uint8_t zoneIdx = zone_number - 1;

    //Extract the first [real,imaginary] pair from bfOutput as input for the Vital Signs Chain
    //Potentially average several pairs from bfOutput to calculate the phase value
    doppler_input_real = *oddemo_dataPathObj.bfOutput;
    oddemo_dataPathObj.bfOutput++;

    doppler_input_imag = *oddemo_dataPathObj.bfOutput;
    oddemo_dataPathObj.bfOutput--;

    //Beginning of vital signs processing chain
    /*Vital Signs Processing Chain*/
//        gFrameCount++;                         // Increment the Global Frame Count
//        static uint16_t frameCountLocal;       // Local circular count
    uint16_t loopIndexBuffer;              // Index that loops over the buffers

    /* Obtain GUI-related Flags sent through the CLI */
//        VitalSignsDemo_GuiMonSel *pGuiMonSel;
//        pGuiMonSel = (VitalSignsDemo_GuiMonSel *) &(gMmwDssMCB.cliCfg[0].vitalSigns_GuiMonSel);
    /* Variables for Impulse Noise Removal */
    static float dataCurr[] = { 0, 0, 0, 0 };
    static float dataPrev2[] = { 0, 0, 0, 0 };
    static float dataPrev1[] = { 0, 0, 0, 0 };

    /* Variables for Clutter Removal */
//        float maxValClutter;
//        uint16_t guiFlag_ClutterRemoval  = pGuiMonSel->guiFlag_ClutterRemoval;
//        uint16_t rangeBinMaxClutter;
    /* Variables for Phase Unwrapping */
    static float phasePrevFrame[] = { 0, 0, 0, 0 }; // Phase value of Previous frame (For phase unwrapping)
    static float diffPhaseCorrectionCum[] = { 0, 0, 0, 0 }; // Phase correction cumulative (For phase unwrapping)
    static float phaseUsedComputationPrev[] = { 0, 0, 0, 0 }; // Phase values used for the Previous frame
    float phaseUsedComputation;    // Unwrapped Phase value used for computation

    /* Variables for Detecting Motion Corrupted Segments */
    uint16_t guiFlag_MotionDetection =
            gMmwDssMCB.cliCfg[0].motionDetectionParamsCfg.enabled; // GUI Flag. Set from the configuration File
    uint16_t guiFlag_GainControl =
            gMmwDssMCB.cliCfg[0].motionDetectionParamsCfg.gainControl; // GUI Flag. Set from the configuration File
    uint16_t indexMotionDetection;      // Temporary Index
    float sumEnergy; // Energy in the data-segment checked for presence of large-scale motion

    /* Vital Signs Waveform */
    float outputFilterBreathOut;      // Breathing waveform after the IIR-Filter
    float outputFilterHeartOut;         // Cardiac waveform after the IIR-Filter

    /* Variables for FFT-based Spectral Estimation */
    uint16_t pPeakSortOutIndex[MAX_NUM_PEAKS_SPECTRUM]; // Sorted Peaks in the Spectrum
    uint16_t numPeaks_BreathSpectrum; // Number of Peaks in the Breathing Spectrum
    uint16_t numPeaks_heartSpectrum;  // Number of Peaks in the Cardiac Spectrum
    uint16_t maxIndexHeartBeatSpect, maxIndexBreathSpect; // Indices corresponding to the max peak in the Breathing and Cardiac Spectrum
    uint16_t maxIndexHeartBeatSpect_4Hz; // Indices corresponding to the max peak from [1.6 - 4.0] Hz
    float breathingRateEst_FFT, heartRateEst_FFT; // Vital Signs Estimate based on the FFT
    float heartRateEst_FFT_4Hz;         // Vital Signs Estimate based on the FFT

    /* Confidence Metric associated with the estimates */
    float confidenceMetricBreath[MAX_NUM_PEAKS_SPECTRUM]; // Confidence Metric associated with each Breathing Spectrum Peak
    float confidenceMetricHeart[MAX_NUM_PEAKS_SPECTRUM]; // Confidence Metric associated with each Cardiac Spectrum Peak
    float confidenceMetricBreathOut, confidenceMetricHeartOut; // Confidence Metric associated with the estimates
    float confidenceMetricHeartOut_4Hz; // Confidence Metric for the 1st Heart beat Harmonic
    float confidenceMetricHeartOut_xCorr; // Confidence Metric for the Autocorrelation
    float confidenceMetricBreathOut_xCorr; // Confidence Metric for the Autocorrelation based breathing-rate estimate
    float peakValueBreathSpect;

//this is too much for the stack:
    /* Variables for peak-counting */
    uint16_t pPeakLocsHeart[MAX_PEAKS_ALLOWED_WFM]; // Peak locations (indices) of the Cardiac Waveform
    uint16_t pPeakLocsBreath[MAX_PEAKS_ALLOWED_WFM]; // Peak locations (indices) of the Breathing Waveform
    uint16_t pPeakLocsValid[MAX_PEAKS_ALLOWED_WFM]; // Peak locations after only retaining the valid peaks
    uint16_t numPeaksBreath, numPeaksHeart; // Number of peaks in the time-domain filtered waveform
    float breathingRateEst_peakCount, heartRateEst_peakCount; // Vital Signs Estimate based on peak-Interval
    float heartRateEst_peakCount_filtered; // Heart-rate peak-interval based estimate after filtering

    /* Exponential smoothing filter */
    static float breathWfmOutUpdated[] = { 0, 0, 0, 0 };
    static float heartWfmOutUpdated[] = { 0, 0, 0, 0 }; // Updated values after exponential smoothing
    float breathWfmOutPrev, heartWfmOutPrev; // Exponential smoothing values at time instance (t-1)
    float sumEnergyBreathWfm, sumEnergyHeartWfm; // These values are used to make a decision if the energy in the waveform is sufficient for it to be classfied as a valid waveform

    /* Variables for Auto-Correlation */
    float heartRateEst_xCorr; // Heart-rate estimate from the Autocorrelation Method
    float breathRateEst_xCorr; // Breathing-rate estimate from the Autocorrelation Method

    /* For FIR Filtering */
    static float pDataIn[][FIR_FILTER_SIZE] = { { 0 }, { 0 }, { 0 }, { 0 } };

    /* Variables for Extracting the Range-FFT output  */
//        uint16_t rangeBinIndex;             // Range-bin Index
    float rangeBinPhase;       // Phase of the Range-bin selected for processing
    static uint16_t rangeBinIndexPhase[4]; // Index of the Range Bin for which the phase is computed
//        uint16_t rangeBinMax;               // Index of the Strongest Range-Bin

    uint16_t indexTemp, indexNumPeaks;  // Temporary Indices
//        int16_t temp_real, temp_imag;       // Temporary variables storing the Real and Imaginary part of a range-bin in the Range-FFT Output
//        float absVal;                       // Absolute value based on the Real and Imaginary values
    float maxVal; // Maximum Value of the Range-Bin in the current range-profile

    maxVal = 0;
    rangeBinPhase = atan2(doppler_input_imag, doppler_input_real);

    obj_VS->unwrapPhasePeak = unwrap(rangeBinPhase, phasePrevFrame[zoneIdx],
                                     &diffPhaseCorrectionCum[zoneIdx]);
    phasePrevFrame[zoneIdx] = rangeBinPhase;

#if TEST_TONE  // Generates a Test-Tone

        float testTone_ampBreath_mm  = 1.0;
        float testTone_freqBreath_Hz = 0.45;
        float testTone_ampHeart_mm   = 0.01;
        float testTone_freqHeart_Hz  = 1.8;

        obj_VS->unwrapPhasePeak = (4 * PI_/WAVELENGTH_MM) *(
                               (testTone_ampBreath_mm * sin(2*PI_*testTone_freqBreath_Hz* gFrameCount/obj_VS->samplingFreq_Hz))+
                               (testTone_ampHeart_mm  * sin(2*PI_*testTone_freqHeart_Hz*gFrameCount/obj_VS->samplingFreq_Hz))
                               );
         //                     +0.5*(testTone_ampBreath_mm *sin(2*PI_*(2*testTone_freqBreath_Hz)* gFrameCount/obj_VS->samplingFreq_Hz)); // Harmonic signal
    #endif

    // Computes the phase differences between successive phase samples
    if (FLAG_COMPUTE_PHASE_DIFFERENCE)
    {
        phaseUsedComputation = obj_VS->unwrapPhasePeak
                - phaseUsedComputationPrev[zoneIdx];
        phaseUsedComputationPrev[zoneIdx] = obj_VS->unwrapPhasePeak;
    }
    else
    {
        phaseUsedComputation = obj_VS->unwrapPhasePeak;
    }

    // Removes impulse like noise from the waveforms
    if (FLAG_REMOVE_IMPULSE_NOISE)
    {
        dataPrev2[zoneIdx] = dataPrev1[zoneIdx];
        dataPrev1[zoneIdx] = dataCurr[zoneIdx];
        dataCurr[zoneIdx] = phaseUsedComputation;
        phaseUsedComputation = filter_RemoveImpulseNoise(
                dataPrev2[zoneIdx], dataPrev1[zoneIdx], dataCurr[zoneIdx],
                obj_VS->noiseImpulse_Thresh);
    }

    // IIR Filtering
    outputFilterBreathOut = filter_IIR_BiquadCascade(
            phaseUsedComputation, obj_VS->pFilterCoefsBreath,
            obj_VS->pScaleValsBreath, obj_VS->pDelayBreath,
            IIR_FILTER_BREATH_NUM_STAGES);
    outputFilterHeartOut = filter_IIR_BiquadCascade(
            phaseUsedComputation, obj_VS->pFilterCoefsHeart_4Hz,
            obj_VS->pScaleValsHeart_4Hz, obj_VS->pDelayHeart,
            IIR_FILTER_HEART_NUM_STAGES);

    // Copies the "Breathing Waveform" in a circular Buffer
    for (loopIndexBuffer = 1;
            loopIndexBuffer < obj_VS->circularBufferSizeBreath;
            loopIndexBuffer++)
    {
        obj_VS->pVitalSigns_Breath_CircularBuffer[loopIndexBuffer - 1] =
                obj_VS->pVitalSigns_Breath_CircularBuffer[loopIndexBuffer];
    }
    obj_VS->pVitalSigns_Breath_CircularBuffer[obj_VS->circularBufferSizeBreath
            - 1] = outputFilterBreathOut;

    // Detection of Motion corrupted Segments
    if (guiFlag_MotionDetection == 1)
    {
        // Update the Motion Removal Circular Buffer
        for (loopIndexBuffer = 1;
                loopIndexBuffer < obj_VS->motionDetection_BlockSize;
                loopIndexBuffer++)
        {
            obj_VS->pMotionCircularBuffer[loopIndexBuffer - 1] =
                    obj_VS->pMotionCircularBuffer[loopIndexBuffer];
        }
        obj_VS->pMotionCircularBuffer[obj_VS->motionDetection_BlockSize - 1] =
                outputFilterHeartOut;
        indexMotionDetection = gFrameCount % obj_VS->motionDetection_BlockSize;

        // Only perform these steps for every obj_VS->motionDetection_BlockSize sample
        if (indexMotionDetection == 0)
        {
            // Check if the current segment is "Noisy"
            sumEnergy = 0;
            for (loopIndexBuffer = 0;
                    loopIndexBuffer < obj_VS->motionDetection_BlockSize;
                    loopIndexBuffer++)
            {
                sumEnergy += (obj_VS->pMotionCircularBuffer[loopIndexBuffer]
                        * obj_VS->pMotionCircularBuffer[loopIndexBuffer]);
            }

            if (sumEnergy > obj_VS->motionDetection_Thresh)
            {
                obj_VS->motionDetected = 1; // Temporary variable to send to the GUI
            }
            else
            {
                obj_VS->motionDetected = 0; // Temporary variable to send to the GUI
            }

            // If NO motion detected in the current segment
            if (obj_VS->motionDetected == 0)
            {
                uint16_t tempEndIndex;
                //  Shift the current contents of the circular Buffer
                for (loopIndexBuffer = obj_VS->motionDetection_BlockSize;
                        loopIndexBuffer < obj_VS->circularBufferSizeHeart;
                        loopIndexBuffer++)
                {
                    obj_VS->pVitalSigns_Heart_CircularBuffer[loopIndexBuffer
                            - obj_VS->motionDetection_BlockSize] =
                            obj_VS->pVitalSigns_Heart_CircularBuffer[loopIndexBuffer];
                }
                // Copy the current data segment to the end of the Circular Buffer
                for (loopIndexBuffer = 0;
                        loopIndexBuffer < obj_VS->motionDetection_BlockSize;
                        loopIndexBuffer++)
                {
                    tempEndIndex = obj_VS->circularBufferSizeHeart
                            - obj_VS->motionDetection_BlockSize;
                    obj_VS->pVitalSigns_Heart_CircularBuffer[tempEndIndex
                            + loopIndexBuffer] =
                            obj_VS->pMotionCircularBuffer[loopIndexBuffer];
                }
            }
        }
    }
    // If Motion DETECTED then don't UPDATE or SHIFT the values in the buffer
    else // Regular processing
    {
        // Copies the "Cardiac Waveform" in a circular Buffer
        for (loopIndexBuffer = 1;
                loopIndexBuffer < obj_VS->circularBufferSizeHeart;
                loopIndexBuffer++)
        {
            obj_VS->pVitalSigns_Heart_CircularBuffer[loopIndexBuffer - 1] =
                    obj_VS->pVitalSigns_Heart_CircularBuffer[loopIndexBuffer];
        }
        obj_VS->pVitalSigns_Heart_CircularBuffer[obj_VS->circularBufferSizeHeart
                - 1] = outputFilterHeartOut;
    }

    /* Spectral Estimation based on the Inter-Peaks Distance */
    numPeaksHeart = find_Peaks(obj_VS->pVitalSigns_Heart_CircularBuffer,
                               float_type, pPeakLocsHeart, obj_VS->pPeakValues,
                               0, obj_VS->circularBufferSizeHeart - 1);
    if (numPeaksHeart != 0)
    {
        numPeaksHeart = filterPeaksWfm(pPeakLocsHeart, pPeakLocsValid,
                                       numPeaksHeart,
                                       obj_VS->peakDistanceHeart_Min,
                                       obj_VS->peakDistanceHeart_Max);
    }
    heartRateEst_peakCount = CONVERT_HZ_BPM
            * ((numPeaksHeart * obj_VS->samplingFreq_Hz)
                    / obj_VS->circularBufferSizeHeart);

    for (loopIndexBuffer = 1; loopIndexBuffer < FIR_FILTER_SIZE;
            loopIndexBuffer++)
    {
        pDataIn[zoneIdx][loopIndexBuffer - 1] =
                pDataIn[zoneIdx][loopIndexBuffer];
    }
    pDataIn[zoneIdx][FIR_FILTER_SIZE - 1] = heartRateEst_peakCount;
    heartRateEst_peakCount_filtered = filter_FIR(pDataIn[zoneIdx],
                                                 obj_VS->pFilterCoefs,
                                                 FIR_FILTER_SIZE);

    numPeaksBreath = find_Peaks(obj_VS->pVitalSigns_Breath_CircularBuffer,
                                int32_type, pPeakLocsBreath,
                                obj_VS->pPeakValues, 0,
                                obj_VS->circularBufferSizeBreath - 1);
    if (numPeaksBreath != 0)
    {
        numPeaksBreath = filterPeaksWfm(pPeakLocsBreath, pPeakLocsValid,
                                        numPeaksBreath,
                                        obj_VS->peakDistanceBreath_Min,
                                        obj_VS->peakDistanceBreath_Max);
    }

    breathingRateEst_peakCount = CONVERT_HZ_BPM
            * ((numPeaksBreath * obj_VS->samplingFreq_Hz)
                    / obj_VS->circularBufferSizeHeart);
    heartRateEst_peakCount = CONVERT_HZ_BPM
            * ((numPeaksHeart * obj_VS->samplingFreq_Hz)
                    / obj_VS->circularBufferSizeBreath);

    // Input to the FFT needs to be complex
    memset((uint8_t*) obj_VS->pVitalSignsBuffer_Cplx, 0,
           obj_VS->breathingWfm_Spectrum_FftSize * sizeof(cmplx32ReIm_t));
    for (loopIndexBuffer = 0;
            loopIndexBuffer < obj_VS->circularBufferSizeBreath;
            loopIndexBuffer++)
    {
        tempFloat = obj_VS->scale_breathingWfm
                * obj_VS->pVitalSigns_Breath_CircularBuffer[loopIndexBuffer];
        //tempFloat = BREATH_WAVEFORM_SCALAR*obj_VS->pVitalSigns_Breath_CircularBuffer[loopIndexBuffer];
        obj_VS->pVitalSignsBuffer_Cplx[loopIndexBuffer].real =
                (int32_t) tempFloat;
        obj_VS->pVitalSignsBuffer_Cplx[loopIndexBuffer].imag = 0;
    }

    // Input is overwritten by the DSP_fft32x32 function
    DSP_fft32x32((int32_t*) obj_VS->pVitalSignsSpectrumTwiddle32x32,
                 obj_VS->breathingWfm_Spectrum_FftSize,
                 (int32_t*) obj_VS->pVitalSignsBuffer_Cplx,
                 (int32_t*) obj_VS->pVitalSigns_SpectrumCplx);

    MmwDemo_magnitudeSquared(obj_VS->pVitalSigns_SpectrumCplx,
                             obj_VS->pVitalSigns_Breath_AbsSpectrum,
                             obj_VS->breathingWfm_Spectrum_FftSize);

    memset((uint8_t*) obj_VS->pVitalSignsBuffer_Cplx, 0,
           obj_VS->heartWfm_Spectrum_FftSize * sizeof(cmplx32ReIm_t));

    // Pre-Processing Steps for the Cardiac Waveform
    // Perform Automatic Gain Control if enabled from the GUI
    if (guiFlag_GainControl == 1)
    {
        computeAGC(obj_VS->pVitalSigns_Heart_CircularBuffer,
                   obj_VS->circularBufferSizeHeart,
                   obj_VS->motionDetection_BlockSize,
                   obj_VS->motionDetection_Thresh);
    }

    if (guiFlag_MotionDetection == 1)
    {
        outputFilterHeartOut =
                obj_VS->pMotionCircularBuffer[obj_VS->motionDetection_BlockSize
                        - 1];
    }
    else
    {
        outputFilterHeartOut =
                obj_VS->pVitalSigns_Heart_CircularBuffer[obj_VS->circularBufferSizeHeart
                        - 1];
    }

    // Perform Autocorrelation on the Waveform
    if (PERFORM_XCORR)
    {
        float temp;
        // Perform Autocorrelation on the Cardiac-Waveform
        uint16_t xCorr_numPeaks;
        uint16_t maxIndex_lag;
        computeAutoCorrelation(obj_VS->pVitalSigns_Heart_CircularBuffer,
                               obj_VS->circularBufferSizeHeart, obj_VS->pXcorr,
                               obj_VS->xCorr_minLag, obj_VS->xCorr_maxLag);
        xCorr_numPeaks = find_Peaks(obj_VS->pXcorr, float_type,
                                    obj_VS->pPeakIndex, obj_VS->pPeakValues,
                                    obj_VS->xCorr_minLag, obj_VS->xCorr_maxLag);
        maxIndex_lag = computeMaxIndex((float*) obj_VS->pXcorr,
                                       obj_VS->xCorr_minLag,
                                       obj_VS->xCorr_maxLag);
        temp = (float) (1.0) / (maxIndex_lag / obj_VS->samplingFreq_Hz);
        heartRateEst_xCorr = (float) CONVERT_HZ_BPM * temp;

        if (xCorr_numPeaks == 0)
        {
            confidenceMetricHeartOut_xCorr = 0;
        }
        else
        {
            confidenceMetricHeartOut_xCorr = obj_VS->pXcorr[maxIndex_lag];
        }

        // Auto-correlation on the Breathing Waveform
        computeAutoCorrelation(obj_VS->pVitalSigns_Breath_CircularBuffer,
                               obj_VS->circularBufferSizeBreath, obj_VS->pXcorr,
                               obj_VS->xCorr_Breath_minLag,
                               obj_VS->xCorr_Breath_maxLag);
        xCorr_numPeaks = find_Peaks(obj_VS->pXcorr, float_type,
                                    obj_VS->pPeakIndex, obj_VS->pPeakValues,
                                    obj_VS->xCorr_Breath_minLag,
                                    obj_VS->xCorr_Breath_maxLag);
        maxIndex_lag = computeMaxIndex((float*) obj_VS->pXcorr,
                                       obj_VS->xCorr_Breath_minLag,
                                       obj_VS->xCorr_Breath_maxLag);
        temp = (float) (1.0) / (maxIndex_lag / obj_VS->samplingFreq_Hz);
        breathRateEst_xCorr = (float) CONVERT_HZ_BPM * temp;

        if (xCorr_numPeaks == 0)
        {
            confidenceMetricBreathOut_xCorr = 0;
        }
        else
        {
            confidenceMetricBreathOut_xCorr = obj_VS->pXcorr[maxIndex_lag];
        }
    }

    // Apply Window on the Cardiac Waveform prior to FFT-based spectral estimation
    // and copies the Pre-processed data to pCircularBufferHeart
    if (FLAG_APPLY_WINDOW)
    {
        uint16_t index_win;
        uint16_t index_WinEnd;

        //float tempFloat;
        index_WinEnd = obj_VS->circularBufferSizeHeart - 1;
        for (index_win = 0; index_win < DOPPLER_WINDOW_SIZE; index_win++)
        {
            tempFloat = obj_VS->pDopplerWindow[index_win];
            obj_VS->pVitalSignsBuffer_Cplx[index_win].real =
                    (int32_t) obj_VS->scale_heartWfm * tempFloat
                            * obj_VS->pVitalSigns_Heart_CircularBuffer[index_win];
            obj_VS->pVitalSignsBuffer_Cplx[index_WinEnd].real =
                    (int32_t) obj_VS->scale_heartWfm * tempFloat
                            * obj_VS->pVitalSigns_Heart_CircularBuffer[index_WinEnd];
            obj_VS->pVitalSignsBuffer_Cplx[index_win].imag = 0;
            obj_VS->pVitalSignsBuffer_Cplx[index_WinEnd].imag = 0;
            index_WinEnd--;
        }
        for (loopIndexBuffer = DOPPLER_WINDOW_SIZE;
                loopIndexBuffer
                        < obj_VS->circularBufferSizeHeart - DOPPLER_WINDOW_SIZE;
                loopIndexBuffer++)
        {
            tempFloat = obj_VS->scale_heartWfm
                    * obj_VS->pVitalSigns_Heart_CircularBuffer[loopIndexBuffer];
            //tempFloat = HEART_WAVEFORM_SCALAR*obj_VS->pVitalSigns_Heart_CircularBuffer[loopIndexBuffer];
            obj_VS->pVitalSignsBuffer_Cplx[loopIndexBuffer].real =
                    (int32_t) tempFloat;
            obj_VS->pVitalSignsBuffer_Cplx[loopIndexBuffer].imag = 0;
        }
    }
    else
    {
        for (loopIndexBuffer = 0;
                loopIndexBuffer < obj_VS->circularBufferSizeHeart;
                loopIndexBuffer++)
        {
            //obj_VS->pVitalSignsBuffer_Cplx[loopIndexBuffer].real = (int32_t) (obj_VS->scale_heartWfm)*obj_VS->pVitalSigns_Heart_CircularBuffer[loopIndexBuffer];
            tempFloat = obj_VS->scale_heartWfm
                    * obj_VS->pVitalSigns_Heart_CircularBuffer[loopIndexBuffer];
            //tempFloat = HEART_WAVEFORM_SCALAR*obj_VS->pVitalSigns_Heart_CircularBuffer[loopIndexBuffer];
            obj_VS->pVitalSignsBuffer_Cplx[loopIndexBuffer].real =
                    (int32_t) tempFloat;
            obj_VS->pVitalSignsBuffer_Cplx[loopIndexBuffer].imag = 0;
        }
    }

    // FFT of the Cardiac Waveform
    DSP_fft32x32((int32_t*) obj_VS->pVitalSignsSpectrumTwiddle32x32,
                 obj_VS->heartWfm_Spectrum_FftSize,
                 (int32_t*) obj_VS->pVitalSignsBuffer_Cplx,
                 (int32_t*) obj_VS->pVitalSigns_SpectrumCplx);

    MmwDemo_magnitudeSquared(obj_VS->pVitalSigns_SpectrumCplx,
                             obj_VS->pVitalSigns_Heart_AbsSpectrum,
                             obj_VS->heartWfm_Spectrum_FftSize);

    // Pick the Peaks in the Breathing Spectrum
    numPeaks_BreathSpectrum = find_Peaks(obj_VS->pVitalSigns_Breath_AbsSpectrum,
                                         float_type, obj_VS->pPeakIndex,
                                         obj_VS->pPeakValues,
                                         obj_VS->breath_startFreq_Index,
                                         obj_VS->breath_endFreq_Index);
    indexNumPeaks =
            (numPeaks_BreathSpectrum < MAX_NUM_PEAKS_SPECTRUM) ?
                    numPeaks_BreathSpectrum : MAX_NUM_PEAKS_SPECTRUM;

    if (indexNumPeaks != 0)
    {
        heapsort_index(obj_VS->pPeakValues, numPeaks_BreathSpectrum,
                       obj_VS->pPeakIndexSorted);
        for (indexTemp = 0; indexTemp < indexNumPeaks; indexTemp++)
        {
            pPeakSortOutIndex[indexTemp] =
                    obj_VS->pPeakIndex[obj_VS->pPeakIndexSorted[numPeaks_BreathSpectrum
                            - indexTemp - 1]];
            confidenceMetricBreath[indexTemp] = computeConfidenceMetric(
                    obj_VS->pVitalSigns_Breath_AbsSpectrum,
                    obj_VS->confMetric_spectrumBreath_IndexStart,
                    obj_VS->confMetric_spectrumBreath_IndexEnd,
                    pPeakSortOutIndex[indexTemp],
                    obj_VS->confMetric_numIndexAroundPeak_breath);
        }
        maxIndexBreathSpect = pPeakSortOutIndex[0]; // The maximum peak
        confidenceMetricBreathOut = confidenceMetricBreath[0];
    }
    else
    {
        maxIndexBreathSpect = computeMaxIndex(
                (float*) obj_VS->pVitalSigns_Breath_AbsSpectrum,
                obj_VS->breath_startFreq_Index, obj_VS->breath_endFreq_Index);
        confidenceMetricBreathOut = computeConfidenceMetric(
                obj_VS->pVitalSigns_Breath_AbsSpectrum, 0,
                PHASE_FFT_SIZE / 4,
                maxIndexBreathSpect,
                obj_VS->confMetric_numIndexAroundPeak_breath);
    }
    peakValueBreathSpect =
            obj_VS->pVitalSigns_Breath_AbsSpectrum[maxIndexBreathSpect]
                    / (10 * obj_VS->scale_breathingWfm);

    // Pick the Peaks in the Heart Spectrum [1.6 - 4.0 Hz]
    numPeaks_heartSpectrum = find_Peaks(obj_VS->pVitalSigns_Heart_AbsSpectrum,
                                        uint32_type, obj_VS->pPeakIndex,
                                        obj_VS->pPeakValues,
                                        obj_VS->heart_startFreq_Index_1p6Hz,
                                        obj_VS->heart_endFreq_Index_4Hz);
    indexNumPeaks =
            (numPeaks_heartSpectrum < MAX_NUM_PEAKS_SPECTRUM) ?
                    numPeaks_heartSpectrum : MAX_NUM_PEAKS_SPECTRUM;
    if (indexNumPeaks != 0)
    {
        heapsort_index(obj_VS->pPeakValues, numPeaks_heartSpectrum,
                       obj_VS->pPeakIndexSorted);
        for (indexTemp = 0; indexTemp < indexNumPeaks; indexTemp++)
        {
            pPeakSortOutIndex[indexTemp] =
                    obj_VS->pPeakIndex[obj_VS->pPeakIndexSorted[numPeaks_heartSpectrum
                            - indexTemp - 1]];
        }
        maxIndexHeartBeatSpect_4Hz = pPeakSortOutIndex[0]; // The maximum peak
        confidenceMetricHeartOut_4Hz = computeConfidenceMetric(
                obj_VS->pVitalSigns_Heart_AbsSpectrum,
                obj_VS->confMetric_spectrumHeart_IndexStart_1p6Hz,
                obj_VS->confMetric_spectrumHeart_IndexStart_4Hz,
                maxIndexHeartBeatSpect_4Hz,
                obj_VS->confMetric_numIndexAroundPeak_heart);
    }
    else
    {
        maxIndexHeartBeatSpect_4Hz = computeMaxIndex(
                (float*) obj_VS->pVitalSigns_Heart_AbsSpectrum,
                obj_VS->heart_startFreq_Index_1p6Hz,
                obj_VS->heart_endFreq_Index_4Hz);
        confidenceMetricHeartOut_4Hz = computeConfidenceMetric(
                obj_VS->pVitalSigns_Heart_AbsSpectrum, 0,
                PHASE_FFT_SIZE / 4,
                maxIndexHeartBeatSpect_4Hz,
                obj_VS->confMetric_numIndexAroundPeak_heart);
    }
    heartRateEst_FFT_4Hz = (float) CONVERT_HZ_BPM * maxIndexHeartBeatSpect_4Hz
            * (obj_VS->freqIncrement_Hz);

    // If a peak is within [1.6 2.0] Hz then check if a harmonic is present is the cardiac spectrum region [0.8 - 2.0] Hz
    if (heartRateEst_FFT_4Hz < MAX_HEART_RATE_BPM)
    {
        for (indexTemp = 1; indexTemp < numPeaks_heartSpectrum; indexTemp++)

            if (abs(heartRateEst_FFT_4Hz
                    - CONVERT_HZ_BPM * (obj_VS->freqIncrement_Hz)
                            * pPeakSortOutIndex[indexTemp]) < HEART_HAMRONIC_THRESH_BPM)
            {
                heartRateEst_FFT_4Hz = CONVERT_HZ_BPM
                        * (obj_VS->freqIncrement_Hz)
                        * pPeakSortOutIndex[indexTemp];
                break;
            }
    }

    // Pick the Peaks in the Cardiac Spectrum
    numPeaks_heartSpectrum = find_Peaks(obj_VS->pVitalSigns_Heart_AbsSpectrum,
                                        float_type, obj_VS->pPeakIndex,
                                        obj_VS->pPeakValues,
                                        obj_VS->heart_startFreq_Index,
                                        obj_VS->heart_endFreq_Index);
    indexNumPeaks =
            (numPeaks_heartSpectrum < MAX_NUM_PEAKS_SPECTRUM) ?
                    numPeaks_heartSpectrum : MAX_NUM_PEAKS_SPECTRUM;

    if (indexNumPeaks != 0)
    {
        heapsort_index(obj_VS->pPeakValues, numPeaks_heartSpectrum,
                       obj_VS->pPeakIndexSorted);
        for (indexTemp = 0; indexTemp < indexNumPeaks; indexTemp++)
        {
            pPeakSortOutIndex[indexTemp] =
                    obj_VS->pPeakIndex[obj_VS->pPeakIndexSorted[numPeaks_heartSpectrum
                            - indexTemp - 1]];
            confidenceMetricHeart[indexTemp] = computeConfidenceMetric(
                    obj_VS->pVitalSigns_Heart_AbsSpectrum,
                    obj_VS->confMetric_spectrumHeart_IndexStart,
                    obj_VS->confMetric_spectrumHeart_IndexEnd,
                    pPeakSortOutIndex[indexTemp],
                    obj_VS->confMetric_numIndexAroundPeak_heart);
        }
        maxIndexHeartBeatSpect = pPeakSortOutIndex[0]; // The maximum peak
        confidenceMetricHeartOut = confidenceMetricHeart[0];
    }
    else
    {
        maxIndexHeartBeatSpect = computeMaxIndex(
                (float*) obj_VS->pVitalSigns_Heart_AbsSpectrum,
                obj_VS->heart_startFreq_Index, obj_VS->heart_endFreq_Index);
        confidenceMetricHeartOut = computeConfidenceMetric(
                obj_VS->pVitalSigns_Heart_AbsSpectrum, 0,
                PHASE_FFT_SIZE / 4,
                maxIndexHeartBeatSpect,
                obj_VS->confMetric_numIndexAroundPeak_heart);
    }

    //Remove the First Breathing Harmonic (if present in the cardiac Spectrum)
    if (FLAG_HARMONIC_CANCELLATION)
    {
        float diffIndex = abs(
                maxIndexHeartBeatSpect
                        - BREATHING_HARMONIC_NUM * maxIndexBreathSpect);
        if (diffIndex * (obj_VS->freqIncrement_Hz) * CONVERT_HZ_BPM
                < BREATHING_HAMRONIC_THRESH_BPM) // Only cancel the 2nd Breathing Harmonic
        {
            maxIndexHeartBeatSpect = pPeakSortOutIndex[1]; // Pick the 2nd Largest peak in the cardiac-spectrum
            confidenceMetricHeartOut = confidenceMetricHeart[1];
        }
    }

    heartRateEst_FFT = (float) CONVERT_HZ_BPM * maxIndexHeartBeatSpect
            * (obj_VS->freqIncrement_Hz);
    breathingRateEst_FFT = (float) CONVERT_HZ_BPM * maxIndexBreathSpect
            * (obj_VS->freqIncrement_Hz);

#ifdef HARMONICS_ENERGY
    float heartRateEst_HarmonicEnergy;
    float breathRateEst_HarmonicEnergy;
    uint16_t maxIndexHeartBeatSpect_temp;
    uint16_t maxIndexBreathSpect_temp;

    memset((void*) obj_VS->pDataOutTemp, 0,
           obj_VS->heartWfm_Spectrum_FftSize * sizeof(float));
    computeEnergyHarmonics(obj_VS->pVitalSigns_Heart_AbsSpectrum,
                           obj_VS->pDataOutTemp,
                           obj_VS->confMetric_spectrumHeart_IndexStart,
                           obj_VS->confMetric_spectrumHeart_IndexEnd,
                           obj_VS->confMetric_numIndexAroundPeak_heart);
    maxIndexHeartBeatSpect_temp = computeMaxIndex(
            obj_VS->pDataOutTemp, obj_VS->confMetric_spectrumHeart_IndexStart,
            obj_VS->confMetric_spectrumHeart_IndexEnd);
    heartRateEst_HarmonicEnergy = (float) CONVERT_HZ_BPM
            * maxIndexHeartBeatSpect_temp * (obj_VS->freqIncrement_Hz);

    memset((void*) obj_VS->pDataOutTemp, 0,
           obj_VS->heartWfm_Spectrum_FftSize * sizeof(float));
    computeEnergyHarmonics(obj_VS->pVitalSigns_Breath_AbsSpectrum,
                           obj_VS->pDataOutTemp,
                           obj_VS->confMetric_spectrumBreath_IndexStart,
                           obj_VS->confMetric_spectrumBreath_IndexEnd,
                           obj_VS->confMetric_numIndexAroundPeak_breath);
    maxIndexBreathSpect_temp = computeMaxIndex(
            obj_VS->pDataOutTemp, obj_VS->confMetric_spectrumBreath_IndexStart,
            obj_VS->confMetric_spectrumBreath_IndexEnd);
    breathRateEst_HarmonicEnergy = (float) CONVERT_HZ_BPM
            * maxIndexBreathSpect_temp * (obj_VS->freqIncrement_Hz);

#endif

    //  Median Value for Heart Rate and Breathing Rate based on 'MEDIAN_WINDOW_LENGTH' previous estimates
    if (FLAG_MEDIAN_FILTER)
    {
        for (loopIndexBuffer = 1; loopIndexBuffer < MEDIAN_WINDOW_LENGTH;
                loopIndexBuffer++)
        {
            obj_VS->pBufferHeartRate[loopIndexBuffer - 1] =
                    obj_VS->pBufferHeartRate[loopIndexBuffer];
            obj_VS->pBufferBreathingRate[loopIndexBuffer - 1] =
                    obj_VS->pBufferBreathingRate[loopIndexBuffer];
            obj_VS->pBufferHeartRate_4Hz[loopIndexBuffer - 1] =
                    obj_VS->pBufferHeartRate_4Hz[loopIndexBuffer];
        }
        obj_VS->pBufferHeartRate[MEDIAN_WINDOW_LENGTH - 1] = heartRateEst_FFT;
        obj_VS->pBufferBreathingRate[MEDIAN_WINDOW_LENGTH - 1] =
                breathingRateEst_FFT;
        obj_VS->pBufferHeartRate_4Hz[MEDIAN_WINDOW_LENGTH - 1] =
                heartRateEst_FFT_4Hz;

        heapsort_index(obj_VS->pBufferHeartRate, MEDIAN_WINDOW_LENGTH,
                       obj_VS->pPeakSortTempIndex);
        heartRateEst_FFT =
                obj_VS->pBufferHeartRate[obj_VS->pPeakSortTempIndex[MEDIAN_WINDOW_LENGTH
                        / 2 + 1]];

        heapsort_index(obj_VS->pBufferHeartRate_4Hz, MEDIAN_WINDOW_LENGTH,
                       obj_VS->pPeakSortTempIndex);
        heartRateEst_FFT_4Hz =
                obj_VS->pBufferHeartRate_4Hz[obj_VS->pPeakSortTempIndex[MEDIAN_WINDOW_LENGTH
                        / 2 + 1]];

        heapsort_index(obj_VS->pBufferBreathingRate, MEDIAN_WINDOW_LENGTH,
                       obj_VS->pPeakSortTempIndex);
        breathingRateEst_FFT =
                obj_VS->pBufferBreathingRate[obj_VS->pPeakSortTempIndex[MEDIAN_WINDOW_LENGTH
                        / 2 + 1]];
    }

    // Exponential Smoothing
    breathWfmOutPrev = breathWfmOutUpdated[zoneIdx];
    breathWfmOutUpdated[zoneIdx] = (obj_VS->alpha_breathing)
            * (outputFilterBreathOut * outputFilterBreathOut)
            + (1 - (obj_VS->alpha_breathing)) * breathWfmOutPrev; // Exponential Smoothing
    sumEnergyBreathWfm = breathWfmOutUpdated[zoneIdx] * 10000;

    heartWfmOutPrev = heartWfmOutUpdated[zoneIdx];
    heartWfmOutUpdated[zoneIdx] = (obj_VS->alpha_heart)
            * (outputFilterHeartOut * outputFilterHeartOut)
            + (1 - (obj_VS->alpha_heart)) * heartWfmOutPrev; // Exponential Smoothing
    sumEnergyHeartWfm = heartWfmOutUpdated[zoneIdx] * 10000;

    // Output Values
    obj_VS->VitalSigns_Output.unwrapPhasePeak_mm = obj_VS->unwrapPhasePeak;
    obj_VS->VitalSigns_Output.outputFilterBreathOut = outputFilterBreathOut;
    obj_VS->VitalSigns_Output.outputFilterHeartOut = outputFilterHeartOut;
    obj_VS->VitalSigns_Output.rangeBinIndexPhase = rangeBinIndexPhase[zoneIdx]; //frameCountLocal;
    obj_VS->VitalSigns_Output.maxVal = maxVal;
    obj_VS->VitalSigns_Output.sumEnergyHeartWfm = sumEnergyHeartWfm;
    obj_VS->VitalSigns_Output.sumEnergyBreathWfm = sumEnergyBreathWfm
            * peakValueBreathSpect;

    obj_VS->VitalSigns_Output.confidenceMetricBreathOut =
            confidenceMetricBreathOut;
    obj_VS->VitalSigns_Output.confidenceMetricHeartOut =
            confidenceMetricHeartOut; // Confidence Metric associated with the estimates
    obj_VS->VitalSigns_Output.confidenceMetricHeartOut_4Hz =
            confidenceMetricHeartOut_4Hz;
    obj_VS->VitalSigns_Output.confidenceMetricHeartOut_xCorr =
            confidenceMetricHeartOut_xCorr;

    obj_VS->VitalSigns_Output.breathingRateEst_FFT = breathingRateEst_FFT;
    obj_VS->VitalSigns_Output.breathingRateEst_peakCount =
            breathingRateEst_peakCount;
    obj_VS->VitalSigns_Output.heartRateEst_peakCount_filtered =
            heartRateEst_peakCount_filtered;
    obj_VS->VitalSigns_Output.heartRateEst_xCorr = heartRateEst_xCorr;
    obj_VS->VitalSigns_Output.heartRateEst_FFT_4Hz = heartRateEst_FFT_4Hz;
    obj_VS->VitalSigns_Output.heartRateEst_FFT = heartRateEst_FFT;

    obj_VS->VitalSigns_Output.processingCyclesOut =
            gMmwDssMCB.dataPathObj[0].timingInfo.interFrameProcCycles
                    / DSP_CLOCK_MHZ;
    obj_VS->VitalSigns_Output.rangeBinStartIndex = obj_VS->rangeBinStartIndex;
    obj_VS->VitalSigns_Output.rangeBinEndIndex = obj_VS->rangeBinEndIndex;
    obj_VS->VitalSigns_Output.motionDetectedFlag = obj_VS->motionDetected;
    obj_VS->VitalSigns_Output.breathingRateEst_xCorr = breathRateEst_xCorr; //breathRateEst_HarmonicEnergy;
    obj_VS->VitalSigns_Output.confidenceMetricBreathOut_xCorr =
            confidenceMetricBreathOut_xCorr;     //breathRateEst_HarmonicEnergy;
    obj_VS->VitalSigns_Output.breathingRateEst_harmonicEnergy =
            breathRateEst_HarmonicEnergy;         //heartRateEst_HarmonicEnergy;
    obj_VS->VitalSigns_Output.heartRateEst_harmonicEnergy =
            heartRateEst_HarmonicEnergy;

    //End of processing chain

    //Copy results into vital_signs structure based on zone_number
    vital_signs.unwrapped_waveform[zone_number - 1] = obj_VS->unwrapPhasePeak;
    vital_signs.heart_waveform[zone_number - 1] = outputFilterHeartOut;
    vital_signs.breathing_waveform[zone_number - 1] = outputFilterBreathOut;
    vital_signs.heart_rate[zone_number - 1] = heartRateEst_FFT;
    vital_signs.breathing_rate[zone_number - 1] = breathingRateEst_FFT;
}

/**
 *  @b Description
 *  @n
 *    Interchirp processing. It is executed per chirp event, after ADC
 *    buffer is filled with chirp samples.
 *
 *  @retval
 *      Not Applicable.
 */
void MmwDemo_interChirpProcessing(MmwDemo_DSS_DataPathObj *obj,
                                  uint8_t chirpPingPongId)
{
    uint32_t antIndx, waitingTime;
    volatile uint32_t startTime;
    volatile uint32_t startTime1;
    MmwDemo_DSS_dataPathContext_t *context = obj->context;

    waitingTime = 0;
    startTime = Cycleprofiler_getTimeStamp();

    /* Kick off DMA to fetch data from ADC buffer for first channel */
    EDMA_startDmaTransfer(context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
    MMW_EDMA_CH_1D_IN_PING);

    /* 1d fft for first antenna, followed by kicking off the DMA of fft output */
    for (antIndx = 0; antIndx < obj->numRxAntennas; antIndx++)
    {
        /* kick off DMA to fetch data for next antenna */
        if (antIndx < (obj->numRxAntennas - 1))
        {
            if (isPong(antIndx))
            {
                EDMA_startDmaTransfer(
                        context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
                        MMW_EDMA_CH_1D_IN_PING);
            }
            else
            {
                EDMA_startDmaTransfer(
                        context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
                        MMW_EDMA_CH_1D_IN_PONG);
            }
        }

        /* verify if DMA has completed for current antenna */
        startTime1 = Cycleprofiler_getTimeStamp();
        MmwDemo_dataPathWait1DInputData(obj, pingPongId(antIndx));
        waitingTime += Cycleprofiler_getTimeStamp() - startTime1;

        mmwavelib_windowing16x16(
                (int16_t*) &obj->adcDataIn[pingPongId(antIndx)
                        * obj->numRangeBins],
                (int16_t*) obj->window1D, obj->numAdcSamples);
        memset((void*) &obj->adcDataIn[pingPongId(antIndx) * obj->numRangeBins
                       + obj->numAdcSamples],
               0,
               (obj->numRangeBins - obj->numAdcSamples)
                       * sizeof(cmplx16ReIm_t));

        DSP_fft16x16(
                (int16_t*) obj->twiddle16x16_1D,
                obj->numRangeBins,
                (int16_t*) &obj->adcDataIn[pingPongId(antIndx)
                        * obj->numRangeBins],
                (int16_t*) &obj->fftOut1D[chirpPingPongId
                        * (obj->numRxAntennas * obj->numRangeBins)
                        + (obj->numRangeBins * antIndx)]);

    }

    if (obj->cliCfg->calibDcRangeSigCfg.enabled)
    {
        MmwDemo_dcRangeSignatureCompensation(obj, chirpPingPongId);
    }

    gCycleLog.interChirpProcessingTime += Cycleprofiler_getTimeStamp()
            - startTime - waitingTime;
    gCycleLog.interChirpWaitTime += waitingTime;
}

void Re16bitIm16bit_swap(cmplx16ReIm_t *input, cmplx16ImRe_t *output, int n)
{

    int i;

    for (i = 0; i < n; i++)
    {
        output[i].real = input[i].real;
        output[i].imag = input[i].imag;
    }
}

/**
 *  @b Description
 *  @n
 *    Interframe processing. It is called from MmwDemo_dssDataPathProcessEvents
 *    after all chirps of the frame have been received and 1D FFT processing on them
 *    has been completed.
 *
 *  @retval
 *      Not Applicable.
 */
cmplx16ImRe_t oddemo_swapBuf[8 * 256]; //max virtual ant * max chirp loops

#ifdef ODDEMO_PROFILE
int t_start, t_stop, t_overhead, t_opt;
#endif

void MmwDemo_interFrameProcessing(MmwDemo_DSS_DataPathObj *obj)
{
    uint32_t rangeIdx;
    uint32_t chirpStep;
//    uint16_t pair;
    uint8_t rangeIndx, azimuthIndx;
//    float doppler_input_real, doppler_input_imag;
    uint8_t vital_signs_zone_number;
//    MmwDemo_GuiMonSel *guiMonSel = &gMmwDssMCB.cliCfg[0].guiMonSel;

    VitalSigns_DataPathObj *obj_VS_zone_1;
    VitalSigns_DataPathObj *obj_VS_zone_2;
    VitalSigns_DataPathObj *obj_VS_zone_3;
    VitalSigns_DataPathObj *obj_VS_zone_4;
    obj_VS_zone_1 = &vitalsigns_dataPathObj_zone_1;
    obj_VS_zone_2 = &vitalsigns_dataPathObj_zone_2;
    obj_VS_zone_3 = &vitalsigns_dataPathObj_zone_3;
    obj_VS_zone_4 = &vitalsigns_dataPathObj_zone_4;

    gFrameCount++;                         // Increment the Global Frame Count
    static uint16_t frameCountLocal;       // Local circular count

    frameCountLocal = gFrameCount % RESET_LOCAL_COUNT_VAL; // Checks for new positions each x frames

    // Beginning of Occupancy_Detection Processing Chain

#ifdef ODDEMO_PROFILE
    TSCL    = 0;
    t_start = _itoll(TSCH, TSCL);
    t_stop  = _itoll(TSCH, TSCL);

    t_overhead = t_stop - t_start;
    t_start = _itoll(TSCH, TSCL);
    #endif

    oddemo_dataPathObj.log2Chirps = obj->log2NumDopplerBins;
    chirpStep = oddemo_dataPathObj.nRxAnt * oddemo_dataPathObj.nChirps;

    cmplx16ReIm_t *DVS_InputSamples;

    for (rangeIdx = 0; rangeIdx < ODDEMO_MAX_RANGE; rangeIdx++)
    {
        ODDemo_Heatmap_aoaEstCaponBF_range(&oddemo_dataPathObj);

        oddemo_dataPathObj.inputAntSamples += chirpStep;
        oddemo_dataPathObj.rangeAzimuthHeatMap += ODDEMO_MAX_AZIMUTH;
    }

    oddemo_dataPathObj.inputAntSamples = (cmplx16ReIm_t*) &gMmwL3[0];
    oddemo_dataPathObj.rangeAzimuthHeatMap = &oddemo_rangeAzimuthHeatMap[0];

    ODDemo_Heatmap_arc_removal(oddemo_rangeAzimuthHeatMap);

//    if (guiMonSel->decision == 1)
//    {
//        for (pair = 0; pair < oddemo_zone_pairs; pair++)
//        {
//            ODDemo_Feature_extract(pair, oddemo_rangeAzimuthHeatMap,
//                                   &oddemo_feature[pair]);
//
//            if (guiMonSel->decision == 1)
//            {
//                int8_t *dptr = &oddemo_decision[pair * 2];
//
//                ODDemo_Decision_process(&oddemo_coeffMatrix[pair][0],
//                                        &oddemo_feature[pair], dptr);
//            }
//        }
//    }

    // End of Occupancy_Detection Processing Chain

    // Beginning of dynamic zone selection
    if (frameCountLocal == 1)
    {
        heatmapGetPeaks(oddemo_rangeAzimuthHeatMap);
    }
    // End of dynamic zone selection

    // Beginning of dynamic multiple vital sign estimation
    uint8_t personsDetected = peak_positions[8];
    float *dopplerFFTInput;
    uint32_t scratchOffset = 0;
    if (personsDetected > 0)
    {
        //Beginning of Vital Signs Processing Chain - Zone 1
        vital_signs_zone_number = 1;
        //Extract range and index bins based on zone definitions
        rangeIndx = peak_positions[0];
        azimuthIndx = peak_positions[1];

        //Index to the correct starting point for the Input Samples
        DVS_InputSamples = (cmplx16ReIm_t*) &gMmwL3[0];
        DVS_InputSamples += (chirpStep * rangeIndx);
        //    Re16bitIm16bit_swap(DVS_InputSamples, oddemo_swapBuf, chirpStep);

        //    /*Calculate covariance matrix and invert */
        ODDemo_Heatmap_aoaEstCaponBF_covInv(
                (uint8_t) (1),                       //Conventional Capon BF = 1
                (uint8_t) (0),                             //Clutter Removal = 0
                oddemo_parms.gamma,
                (int32_t) oddemo_dataPathObj.nRxAnt,
                (int32_t) (1),                            //Number of Chirps = 1
                (int32_t) (1),                //not used without clutter removal
                (int32_t*) &oddemo_dataPathObj.scratchPad[0],
                (cplx16_t*) DVS_InputSamples,
                (cplxf_t*) &oddemo_dataPathObj.invRnMatrix_VS[0]);

        //Capon BF Doppler Estimation
        //Used as input into Vital Signs Processing Chain
        scratchOffset = 0;

        dopplerFFTInput =
                (float*) &oddemo_dataPathObj.scratchPad[scratchOffset];
        oddemo_dataPathObj.bfOutput = dopplerFFTInput;

        Re16bitIm16bit_swap(DVS_InputSamples, oddemo_swapBuf, chirpStep);

        ODDEMO_Heatmap_aoaEstCaponBF_dopplerEstInput(
                (uint8_t) 1,
                (int32_t) oddemo_dataPathObj.nRxAnt,
                (int32_t) 1,
                (cplx16_t*) oddemo_swapBuf,
                (cplxf_t*) &oddemo_dataPathObj.steeringVec[azimuthIndx
                        * (oddemo_dataPathObj.nRxAnt - 1)],
                (cplxf_t*) &oddemo_dataPathObj.invRnMatrix_VS[0],
                (int32_t*) &oddemo_dataPathObj.scratchPad[2 * 1], (float) 1.0,
                (float*) oddemo_dataPathObj.bfOutput);

        VS_Processing_Chain(obj_VS_zone_1, vital_signs_zone_number);

        //End of Vital Signs Processing Chain - Zone 1
    }

    if (personsDetected > 1)
    {
        //Beginning of Vital Signs Processing Chain - Zone 2
        vital_signs_zone_number = 2;
        //Extract range and index bins based on zone definitions
        rangeIndx = peak_positions[2];
        azimuthIndx = peak_positions[3];

        //Index to the correct starting point for the Input Samples
        DVS_InputSamples = (cmplx16ReIm_t*) &gMmwL3[0];
        DVS_InputSamples += (chirpStep * rangeIndx);
        //    Re16bitIm16bit_swap(DVS_InputSamples, oddemo_swapBuf, chirpStep);

        //    /*Calculate covariance matrix and invert */
        ODDemo_Heatmap_aoaEstCaponBF_covInv(
                (uint8_t) (1),                       //Conventional Capon BF = 1
                (uint8_t) (0),                             //Clutter Removal = 0
                oddemo_parms.gamma,
                (int32_t) oddemo_dataPathObj.nRxAnt,
                (int32_t) (1),                            //Number of Chirps = 1
                (int32_t) (1),                //not used without clutter removal
                (int32_t*) &oddemo_dataPathObj.scratchPad[0],
                (cplx16_t*) DVS_InputSamples,
                (cplxf_t*) &oddemo_dataPathObj.invRnMatrix_VS[0]);

        //Capon BF Doppler Estimation
        //Used as input into Vital Signs Processing Chain
        scratchOffset = 0;

        dopplerFFTInput =
                (float*) &oddemo_dataPathObj.scratchPad[scratchOffset];
        oddemo_dataPathObj.bfOutput = dopplerFFTInput;

        Re16bitIm16bit_swap(DVS_InputSamples, oddemo_swapBuf, chirpStep);

        ODDEMO_Heatmap_aoaEstCaponBF_dopplerEstInput(
                (uint8_t) 1,
                (int32_t) oddemo_dataPathObj.nRxAnt,
                (int32_t) 1,
                (cplx16_t*) oddemo_swapBuf,
                (cplxf_t*) &oddemo_dataPathObj.steeringVec[azimuthIndx
                        * (oddemo_dataPathObj.nRxAnt - 1)],
                (cplxf_t*) &oddemo_dataPathObj.invRnMatrix_VS[0],
                (int32_t*) &oddemo_dataPathObj.scratchPad[2 * 1], (float) 1.0,
                (float*) oddemo_dataPathObj.bfOutput);

        VS_Processing_Chain(obj_VS_zone_2, vital_signs_zone_number);

        //End of Vital Signs Processing Chain - Zone 2
    }

    if (personsDetected > 2)
    {
        //Beginning of Vital Signs Processing Chain - Zone 3
        vital_signs_zone_number = 3;
        //Extract range and index bins based on zone definitions
        rangeIndx = peak_positions[4];
        azimuthIndx = peak_positions[5];

        //Index to the correct starting point for the Input Samples
        DVS_InputSamples = (cmplx16ReIm_t*) &gMmwL3[0];
        DVS_InputSamples += (chirpStep * rangeIndx);
        //    Re16bitIm16bit_swap(DVS_InputSamples, oddemo_swapBuf, chirpStep);

        //    /*Calculate covariance matrix and invert */
        ODDemo_Heatmap_aoaEstCaponBF_covInv(
                (uint8_t) (1),                       //Conventional Capon BF = 1
                (uint8_t) (0),                             //Clutter Removal = 0
                oddemo_parms.gamma,
                (int32_t) oddemo_dataPathObj.nRxAnt,
                (int32_t) (1),                            //Number of Chirps = 1
                (int32_t) (1),                //not used without clutter removal
                (int32_t*) &oddemo_dataPathObj.scratchPad[0],
                (cplx16_t*) DVS_InputSamples,
                (cplxf_t*) &oddemo_dataPathObj.invRnMatrix_VS[0]);

        //Capon BF Doppler Estimation
        //Used as input into Vital Signs Processing Chain
        scratchOffset = 0;

        dopplerFFTInput =
                (float*) &oddemo_dataPathObj.scratchPad[scratchOffset];
        oddemo_dataPathObj.bfOutput = dopplerFFTInput;

        Re16bitIm16bit_swap(DVS_InputSamples, oddemo_swapBuf, chirpStep);

        ODDEMO_Heatmap_aoaEstCaponBF_dopplerEstInput(
                (uint8_t) 1,
                (int32_t) oddemo_dataPathObj.nRxAnt,
                (int32_t) 1,
                (cplx16_t*) oddemo_swapBuf,
                (cplxf_t*) &oddemo_dataPathObj.steeringVec[azimuthIndx
                        * (oddemo_dataPathObj.nRxAnt - 1)],
                (cplxf_t*) &oddemo_dataPathObj.invRnMatrix_VS[0],
                (int32_t*) &oddemo_dataPathObj.scratchPad[2 * 1], (float) 1.0,
                (float*) oddemo_dataPathObj.bfOutput);

        VS_Processing_Chain(obj_VS_zone_3, vital_signs_zone_number);

        //End of Vital Signs Processing Chain - Zone 3
    }

    if (personsDetected > 3)
    {
        //Beginning of Vital Signs Processing Chain - Zone 4
        vital_signs_zone_number = 4;
        //Extract range and index bins based on zone definitions
        rangeIndx = peak_positions[6];
        azimuthIndx = peak_positions[7];

        //Index to the correct starting point for the Input Samples
        DVS_InputSamples = (cmplx16ReIm_t*) &gMmwL3[0];
        DVS_InputSamples += (chirpStep * rangeIndx);
        //    Re16bitIm16bit_swap(DVS_InputSamples, oddemo_swapBuf, chirpStep);

        //    /*Calculate covariance matrix and invert */
        ODDemo_Heatmap_aoaEstCaponBF_covInv(
                (uint8_t) (1),                       //Conventional Capon BF = 1
                (uint8_t) (0),                             //Clutter Removal = 0
                oddemo_parms.gamma,
                (int32_t) oddemo_dataPathObj.nRxAnt,
                (int32_t) (1),                            //Number of Chirps = 1
                (int32_t) (1),                //not used without clutter removal
                (int32_t*) &oddemo_dataPathObj.scratchPad[0],
                (cplx16_t*) DVS_InputSamples,
                (cplxf_t*) &oddemo_dataPathObj.invRnMatrix_VS[0]);

        //Capon BF Doppler Estimation
        //Used as input into Vital Signs Processing Chain
        scratchOffset = 0;

        dopplerFFTInput =
                (float*) &oddemo_dataPathObj.scratchPad[scratchOffset];
        oddemo_dataPathObj.bfOutput = dopplerFFTInput;

        Re16bitIm16bit_swap(DVS_InputSamples, oddemo_swapBuf, chirpStep);

        ODDEMO_Heatmap_aoaEstCaponBF_dopplerEstInput(
                (uint8_t) 1,
                (int32_t) oddemo_dataPathObj.nRxAnt,
                (int32_t) 1,
                (cplx16_t*) oddemo_swapBuf,
                (cplxf_t*) &oddemo_dataPathObj.steeringVec[azimuthIndx
                        * (oddemo_dataPathObj.nRxAnt - 1)],
                (cplxf_t*) &oddemo_dataPathObj.invRnMatrix_VS[0],
                (int32_t*) &oddemo_dataPathObj.scratchPad[2 * 1], (float) 1.0,
                (float*) oddemo_dataPathObj.bfOutput);

        VS_Processing_Chain(obj_VS_zone_4, vital_signs_zone_number);

        //End of Vital Signs Processing Chain - Zone 4
    }
    // End of dynamic multiple vital sign estimation

    oddemo_scratch_output[50] = vital_signs.unwrapped_waveform[0];
    oddemo_scratch_output[51] = vital_signs.heart_waveform[0];
    oddemo_scratch_output[52] = vital_signs.breathing_waveform[0];
    oddemo_scratch_output[53] = vital_signs.heart_rate[0];
    oddemo_scratch_output[54] = vital_signs.breathing_rate[0];

    oddemo_scratch_output[55] = vital_signs.unwrapped_waveform[1];
    oddemo_scratch_output[56] = vital_signs.heart_waveform[1];
    oddemo_scratch_output[57] = vital_signs.breathing_waveform[1];
    oddemo_scratch_output[58] = vital_signs.heart_rate[1];
    oddemo_scratch_output[59] = vital_signs.breathing_rate[1];

    oddemo_scratch_output[60] = vital_signs.unwrapped_waveform[2];
    oddemo_scratch_output[61] = vital_signs.heart_waveform[2];
    oddemo_scratch_output[62] = vital_signs.breathing_waveform[2];
    oddemo_scratch_output[63] = vital_signs.heart_rate[2];
    oddemo_scratch_output[64] = vital_signs.breathing_rate[2];

    oddemo_scratch_output[65] = vital_signs.unwrapped_waveform[3];
    oddemo_scratch_output[66] = vital_signs.heart_waveform[3];
    oddemo_scratch_output[67] = vital_signs.breathing_waveform[3];
    oddemo_scratch_output[68] = vital_signs.heart_rate[3];
    oddemo_scratch_output[69] = vital_signs.breathing_rate[3];

#ifdef ODDEMO_PROFILE
    t_stop = _itoll(TSCH, TSCL);
    t_opt  = (t_stop - t_start) - t_overhead;
    #endif
}

/**
 *  @b Description
 *  @n
 *    Chirp processing. It is called from MmwDemo_dssDataPathProcessEvents. It
 *    is executed per chirp
 *
 *  @retval
 *      Not Applicable.
 */
void MmwDemo_processChirp(MmwDemo_DSS_DataPathObj *obj,
                          uint16_t chirpIndxInMultiChirp)
{
    volatile uint32_t startTime;
    MmwDemo_DSS_dataPathContext_t *context = obj->context;

    startTime = Cycleprofiler_getTimeStamp();
    if (obj->chirpCount > 1) //verify if ping(or pong) buffer is free for odd(or even) chirps
    {
        MmwDemo_dataPathWait1DOutputData(obj, pingPongId(obj->chirpCount));
    }
    gCycleLog.interChirpWaitTime += Cycleprofiler_getTimeStamp() - startTime;

    MmwDemo_interChirpProcessing(obj, pingPongId(obj->chirpCount));

    /*Modify destination address in Param set and DMA for sending 1DFFT output (for all antennas) to L3  */
    if (isPong(obj->chirpCount))
    {
        EDMAutil_triggerType3(
                context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
                (uint8_t*) NULL,
                (uint8_t*) (&obj->radarCube[(obj->numDopplerBins
                        * obj->numRxAntennas * (obj->numTxAntennas - 1))
                        + obj->dopplerBinCount]),
                (uint8_t) MMW_EDMA_CH_1D_OUT_PONG,
                (uint8_t) MMW_EDMA_TRIGGER_ENABLE);
    }
    else
    {
        EDMAutil_triggerType3(
                context->edmaHandle[MMW_DATA_PATH_EDMA_INSTANCE],
                (uint8_t*) NULL,
                (uint8_t*) (&obj->radarCube[obj->dopplerBinCount]),
                (uint8_t) MMW_EDMA_CH_1D_OUT_PING,
                (uint8_t) MMW_EDMA_TRIGGER_ENABLE);
    }

    obj->chirpCount++;
    obj->txAntennaCount++;
    if (obj->txAntennaCount == obj->numTxAntennas)
    {
        obj->txAntennaCount = 0;
        obj->dopplerBinCount++;
        if (obj->dopplerBinCount == obj->numDopplerBins)
        {
            obj->dopplerBinCount = 0;
            obj->chirpCount = 0;
        }
    }
}

/**
 *  @b Description
 *  @n
 *  Wait for transfer of data corresponding to the last 2 chirps (ping/pong)
 *  to the radarCube matrix before starting interframe processing.
 *  @retval
 *      Not Applicable.
 */
void MmwDemo_waitEndOfChirps(MmwDemo_DSS_DataPathObj *obj)
{
    volatile uint32_t startTime;

    startTime = Cycleprofiler_getTimeStamp();
    /* Wait for transfer of data corresponding to last 2 chirps (ping/pong) */
    MmwDemo_dataPathWait1DOutputData(obj, 0);
    MmwDemo_dataPathWait1DOutputData(obj, 1);

    gCycleLog.interChirpWaitTime += Cycleprofiler_getTimeStamp() - startTime;
}

void MmwDemo_edmaErrorCallbackFxn(EDMA_Handle handle,
                                  EDMA_errorInfo_t *errorInfo)
{
    MmwDemo_dssAssert(0);
}

void MmwDemo_edmaTransferControllerErrorCallbackFxn(
        EDMA_Handle handle, EDMA_transferControllerErrorInfo_t *errorInfo)
{
    MmwDemo_dssAssert(0);
}

void MmwDemo_dataPathObjInit(MmwDemo_DSS_DataPathObj *obj,
                             MmwDemo_DSS_dataPathContext_t *context,
                             MmwDemo_CliCfg_t *cliCfg,
                             MmwDemo_CliCommonCfg_t *cliCommonCfg,
                             MmwDemo_Cfg *cfg)
{
    memset((void*) obj, 0, sizeof(MmwDemo_DSS_DataPathObj));
    obj->context = context;
    obj->cliCfg = cliCfg;
    obj->cliCommonCfg = cliCommonCfg;
    obj->cfg = cfg;
}

void MmwDemo_dataPathInit1Dstate(MmwDemo_DSS_DataPathObj *obj)
{
    obj->chirpCount = 0;
    obj->dopplerBinCount = 0;
    obj->txAntennaCount = 0;

    /* reset profiling logs before start of frame */
    memset((void*) &gCycleLog, 0, sizeof(cycleLog_t));
}

void MmwDemo_dataPathDeleteSemaphore(MmwDemo_DSS_dataPathContext_t *context)
{
#ifdef EDMA_1D_INPUT_BLOCKING
    Semaphore_delete(&context->EDMA_1D_InputDone_semHandle[0]);
    Semaphore_delete(&context->EDMA_1D_InputDone_semHandle[1]);
#endif
#ifdef EDMA_1D_OUTPUT_BLOCKING
    Semaphore_delete(&context->EDMA_1D_OutputDone_semHandle[0]);
    Semaphore_delete(&context->EDMA_1D_OutputDone_semHandle[1]);
#endif
}

int32_t MmwDemo_dataPathInitEdma(MmwDemo_DSS_dataPathContext_t *context)
{
    Semaphore_Params semParams;
    uint8_t numInstances;
    int32_t errorCode;
    EDMA_Handle handle;
    EDMA_errorConfig_t errorConfig;
    uint32_t instanceId;
    EDMA_instanceInfo_t instanceInfo;

    Semaphore_Params_init(&semParams);
    semParams.mode = Semaphore_Mode_BINARY;
#ifdef EDMA_1D_INPUT_BLOCKING
    context->EDMA_1D_InputDone_semHandle[0] = Semaphore_create(0, &semParams, NULL);
    context->EDMA_1D_InputDone_semHandle[1] = Semaphore_create(0, &semParams, NULL);
#endif
#ifdef EDMA_1D_OUTPUT_BLOCKING
    context->EDMA_1D_OutputDone_semHandle[0] = Semaphore_create(0, &semParams, NULL);
    context->EDMA_1D_OutputDone_semHandle[1] = Semaphore_create(0, &semParams, NULL);
#endif

    numInstances = EDMA_getNumInstances();

    /* Initialize the edma instance to be tested */
    for (instanceId = 0; instanceId < numInstances; instanceId++)
    {
        EDMA_init(instanceId);

        handle = EDMA_open(instanceId, &errorCode, &instanceInfo);
        if (handle == NULL)
        {
            System_printf(
                    "Error: Unable to open the edma Instance, erorCode = %d\n",
                    errorCode);
            return -1;
        }
        context->edmaHandle[instanceId] = handle;

        errorConfig.isConfigAllEventQueues = true;
        errorConfig.isConfigAllTransferControllers = true;
        errorConfig.isEventQueueThresholdingEnabled = true;
        errorConfig.eventQueueThreshold = EDMA_EVENT_QUEUE_THRESHOLD_MAX;
        errorConfig.isEnableAllTransferControllerErrors = true;
        errorConfig.callbackFxn = MmwDemo_edmaErrorCallbackFxn;
        errorConfig.transferControllerCallbackFxn =
                MmwDemo_edmaTransferControllerErrorCallbackFxn;
        if ((errorCode = EDMA_configErrorMonitoring(handle, &errorConfig))
                != EDMA_NO_ERROR)
        {
            System_printf(
                    "Debug: EDMA_configErrorMonitoring() failed with errorCode = %d\n",
                    errorCode);
            return -1;
        }
    }
    return 0;
}

void MmwDemo_printHeapStats(char *name, uint32_t heapUsed, uint32_t heapSize)
{
    System_printf("Heap %s : size %d (0x%x), free %d (0x%x)\n", name, heapSize,
                  heapSize, heapSize - heapUsed, heapSize - heapUsed);
}

#define SOC_MAX_NUM_RX_ANTENNAS SYS_COMMON_NUM_RX_CHANNEL
#define SOC_MAX_NUM_TX_ANTENNAS SYS_COMMON_NUM_TX_ANTENNAS

void MmwDemo_dataPathComputeDerivedConfig(MmwDemo_DSS_DataPathObj *obj)
{
    obj->log2NumDopplerBins = MmwDemo_floorLog2(obj->numDopplerBins);

//    /* check for numDopplerBins to be exact power of 2 */
//    if ((1U << obj->log2NumDopplerBins) != obj->numDopplerBins)
//    {
//        System_printf("Number of doppler bins must be a power of 2\n");
//        MmwDemo_dssAssert(0);
//    }
}

void MmwDemo_dataPathConfigBuffers(MmwDemo_DSS_DataPathObj *obj,
                                   uint32_t adcBufAddress,
                                   VitalSigns_DataPathObj *obj_VS_1,
                                   VitalSigns_DataPathObj *obj_VS_2,
                                   VitalSigns_DataPathObj *obj_VS_3,
                                   VitalSigns_DataPathObj *obj_VS_4)
{
    /* below defines for debugging purposes, do not remove as overlays can be hard to debug */
//#define NO_L1_ALLOC /* don't allocate from L1D, use L2 instead */
//#define NO_OVERLAY  /* do not overlay */
#define ALIGN(x,a)  (((x)+((a)-1))&~((a)-1))

#ifdef NO_OVERLAY
#define MMW_ALLOC_BUF(name, nameType, startAddr, alignment, size) \
        obj->name = (nameType *) ALIGN(prev_end, alignment); \
        prev_end = (uint32_t)obj->name + (size) * sizeof(nameType);
#else
#define MMW_ALLOC_BUF(name, nameType, startAddr, alignment, size) \
        obj->name = (nameType *) ALIGN(startAddr, alignment); \
        uint32_t name##_end = (uint32_t)obj->name + (size) * sizeof(nameType);

#define MMW_ALLOC_BUF_NO_END(name, nameType, startAddr, alignment) \
        obj->name = (nameType *) ALIGN(startAddr, alignment);

#define MMW_ALLOC_BUF_VS_1(name, nameType, startAddr, alignment, size) \
        obj_VS_1->name = (nameType *) ALIGN(startAddr, alignment); \
        uint32_t name##_vs1_end = (uint32_t)obj_VS_1->name + (size) * sizeof(nameType);

#define MMW_ALLOC_BUF_VS_2(name, nameType, startAddr, alignment, size) \
        obj_VS_2->name = (nameType *) ALIGN(startAddr, alignment); \
        uint32_t name##_vs2_end = (uint32_t)obj_VS_2->name + (size) * sizeof(nameType);

#define MMW_ALLOC_BUF_VS_3(name, nameType, startAddr, alignment, size) \
        obj_VS_3->name = (nameType *) ALIGN(startAddr, alignment); \
        uint32_t name##_vs3_end = (uint32_t)obj_VS_3->name + (size) * sizeof(nameType);

#define MMW_ALLOC_BUF_VS_4(name, nameType, startAddr, alignment, size) \
        obj_VS_4->name = (nameType *) ALIGN(startAddr, alignment); \
        uint32_t name##_vs4_end = (uint32_t)obj_VS_4->name + (size) * sizeof(nameType);

#define MMW_ALLOC_VS4_NO_END(name, nameType, startAddr, alignment) \
        obj_VS_4->name = (nameType *) ALIGN(startAddr, alignment);
#endif

    uint32_t heapUsed;
    uint32_t heapL1start = (uint32_t) &gMmwL1[0];
    uint32_t heapL2start = (uint32_t) &gMmwL2[0];
    uint32_t heapL3start = (uint32_t) &gMmwL3[0];

    /* L3 is overlaid with one-time only accessed code. Although heap is not
     required to be initialized to 0, it may help during debugging when viewing memory
     in CCS */
    memset((void*) heapL3start, 0, sizeof(gMmwL3));

    /* L1 allocation

     Buffers are overlayed in the following order. Notation "|" indicates parallel
     and "+" means cascade

     { 1D
     (adcDataIn)
     } |
     { 2D
     (dstPingPong +  fftOut2D) +
     (windowingBuf2D | log2Abs) + sumAbs
     } |
     { CFAR
     detObj2DRaw
     } |
     { 3D
     (
     #ifdef MMW_USE_SINGLE_POINT_DFT
     azimuthIn (must be at least beyond dstPingPong) + azimuthOut + azimuthMagSqr)
     #else
     azimuthIn (must be at least beyond windowingBuf2D) + azimuthOut + azimuthMagSqr)
     #endif
     }
     */

    heapUsed = 2 * obj->numRangeBins * sizeof(cmplx16ReIm_t);

    MMW_ALLOC_BUF_NO_END(adcDataIn, cmplx16ReIm_t, heapL1start,
                         MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN);

    memset((void*) obj->adcDataIn, 0, heapUsed);

    MmwDemo_printHeapStats("L1", heapUsed, MMW_L1_HEAP_SIZE);

    /* L2 allocation
     {
     { 1D
     (fftOut1D)
     } |
     { 2D + 3D
     (cfarDetObjIndexBuf + detDopplerLines.dopplerLineMask) + sumAbsRange
     }
     } +
     {
     twiddle16x16_1D +
     window1D +
     twiddle32x32_2D +
     window2D +
     detObj2D +
     detObj2dAzimIdx +
     azimuthTwiddle32x32 +
     azimuthModCoefs
     }
     */

    MMW_ALLOC_BUF(fftOut1D, cmplx16ReIm_t, heapL2start,
                  MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                  2 * obj->numRxAntennas * obj->numRangeBins);

    MMW_ALLOC_BUF(twiddle16x16_1D, cmplx16ReIm_t, fftOut1D_end,
                  MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN, obj->numRangeBins);

    MMW_ALLOC_BUF(window1D, int16_t, twiddle16x16_1D_end,
                  MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                  obj->numAdcSamples / 2);

    MMW_ALLOC_BUF(twiddle32x32_2D, cmplx32ReIm_t, window1D_end,
                  MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN, obj->numDopplerBins);

    MMW_ALLOC_BUF_NO_END(window2D, int32_t, twiddle32x32_2D_end,
                         MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN);

    heapUsed = twiddle32x32_2D_end - heapL2start;

    MmwDemo_printHeapStats("L2", heapUsed, MMW_L2_HEAP_SIZE);

    /* L3 allocation:
     ADCdataBuf (for unit test) +
     radarCube +
     azimuthStaticHeatMap +
     detMatrix
     */

    MMW_ALLOC_BUF(ADCdataBuf, cmplx16ReIm_t, heapL3start,
                  MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                  obj->numRangeBins * obj->numRxAntennas * obj->numTxAntennas);

    if (adcBufAddress != NULL)
    {
        obj->ADCdataBuf = (cmplx16ReIm_t*) adcBufAddress;
        ADCdataBuf_end = heapL3start;
    }

    MMW_ALLOC_BUF(
            radarCube,
            cmplx16ReIm_t,
            ADCdataBuf_end,
            MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
            obj->numRangeBins * obj->numDopplerBins * obj->numRxAntennas
                    * obj->numTxAntennas);

    MMW_ALLOC_BUF_VS_1(pVitalSigns_Breath_CircularBuffer, float, radarCube_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_1->circularBufferSizeBreath);

    MMW_ALLOC_BUF_VS_1(pVitalSigns_Heart_CircularBuffer, float,
                       pVitalSigns_Breath_CircularBuffer_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_1->circularBufferSizeHeart);

    MMW_ALLOC_BUF_VS_1(pMotionCircularBuffer, float,
                       pVitalSigns_Heart_CircularBuffer_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_1->motionDetection_BlockSize);

    MMW_ALLOC_BUF_VS_1(pDopplerWindow, float, pMotionCircularBuffer_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       DOPPLER_WINDOW_SIZE);

    MMW_ALLOC_BUF_VS_1(pVitalSigns_SpectrumCplx, cmplx32ReIm_t,
                       pDopplerWindow_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_1->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_1(pVitalSignsBuffer_Cplx, cmplx32ReIm_t,
                       pVitalSigns_SpectrumCplx_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_1->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_1(pVitalSignsSpectrumTwiddle32x32, cmplx32ReIm_t,
                       pVitalSignsBuffer_Cplx_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_1->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_1(pVitalSigns_Breath_AbsSpectrum, float,
                       pVitalSignsSpectrumTwiddle32x32_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_1->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_1(pVitalSigns_Heart_AbsSpectrum, float,
                       pVitalSigns_Breath_AbsSpectrum_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_1->heartWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_1(
            pFilterCoefsBreath, float, pVitalSigns_Heart_AbsSpectrum_vs1_end,
            MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
            IIR_FILTER_BREATH_NUM_STAGES * IIR_FILTER_COEFS_SECOND_ORDER);

    MMW_ALLOC_BUF_VS_1(pScaleValsBreath, float, pFilterCoefsBreath_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       IIR_FILTER_BREATH_NUM_STAGES + 1);

    MMW_ALLOC_BUF_VS_1(
            pFilterCoefsHeart_4Hz, float, pScaleValsBreath_vs1_end,
            MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
            IIR_FILTER_HEART_NUM_STAGES * IIR_FILTER_COEFS_SECOND_ORDER);

    MMW_ALLOC_BUF_VS_1(pScaleValsHeart_4Hz, float,
                       pFilterCoefsHeart_4Hz_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       IIR_FILTER_HEART_NUM_STAGES + 1);

    MMW_ALLOC_BUF_VS_1(pXcorr, float, pScaleValsHeart_4Hz_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN, XCORR_NUM_LAGS);

    MMW_ALLOC_BUF_VS_1(pTempReal_Prev, float, pXcorr_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj->numRangeBins);

    MMW_ALLOC_BUF_VS_1(pTempImag_Prev, float, pTempReal_Prev_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj->numRangeBins);

    MMW_ALLOC_BUF_VS_1(pRangeProfileClutterRemoved, float,
                       pTempImag_Prev_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj->numRangeBins);

    MMW_ALLOC_BUF_VS_1(pDataOutTemp, float, pRangeProfileClutterRemoved_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_1->heartWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_2(pVitalSigns_Breath_CircularBuffer, float,
                       pDataOutTemp_vs1_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_2->circularBufferSizeBreath);

    MMW_ALLOC_BUF_VS_2(pVitalSigns_Heart_CircularBuffer, float,
                       pVitalSigns_Breath_CircularBuffer_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_2->circularBufferSizeHeart);

    MMW_ALLOC_BUF_VS_2(pMotionCircularBuffer, float,
                       pVitalSigns_Heart_CircularBuffer_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_2->motionDetection_BlockSize);

    MMW_ALLOC_BUF_VS_2(pDopplerWindow, float, pMotionCircularBuffer_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       DOPPLER_WINDOW_SIZE);

    MMW_ALLOC_BUF_VS_2(pVitalSigns_SpectrumCplx, cmplx32ReIm_t,
                       pDopplerWindow_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_2->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_2(pVitalSignsBuffer_Cplx, cmplx32ReIm_t,
                       pVitalSigns_SpectrumCplx_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_2->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_2(pVitalSignsSpectrumTwiddle32x32, cmplx32ReIm_t,
                       pVitalSignsBuffer_Cplx_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_2->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_2(pVitalSigns_Breath_AbsSpectrum, float,
                       pVitalSignsSpectrumTwiddle32x32_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_2->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_2(pVitalSigns_Heart_AbsSpectrum, float,
                       pVitalSigns_Breath_AbsSpectrum_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_2->heartWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_2(
            pFilterCoefsBreath, float, pVitalSigns_Heart_AbsSpectrum_vs2_end,
            MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
            IIR_FILTER_BREATH_NUM_STAGES * IIR_FILTER_COEFS_SECOND_ORDER);

    MMW_ALLOC_BUF_VS_2(pScaleValsBreath, float, pFilterCoefsBreath_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       IIR_FILTER_BREATH_NUM_STAGES + 1);

    MMW_ALLOC_BUF_VS_2(
            pFilterCoefsHeart_4Hz, float, pScaleValsBreath_vs2_end,
            MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
            IIR_FILTER_HEART_NUM_STAGES * IIR_FILTER_COEFS_SECOND_ORDER);

    MMW_ALLOC_BUF_VS_2(pScaleValsHeart_4Hz, float,
                       pFilterCoefsHeart_4Hz_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       IIR_FILTER_HEART_NUM_STAGES + 1);

    MMW_ALLOC_BUF_VS_2(pXcorr, float, pScaleValsHeart_4Hz_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN, XCORR_NUM_LAGS);

    MMW_ALLOC_BUF_VS_2(pTempReal_Prev, float, pXcorr_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj->numRangeBins);

    MMW_ALLOC_BUF_VS_2(pTempImag_Prev, float, pTempReal_Prev_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj->numRangeBins);

    MMW_ALLOC_BUF_VS_2(pRangeProfileClutterRemoved, float,
                       pTempImag_Prev_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj->numRangeBins);

    MMW_ALLOC_BUF_VS_2(pDataOutTemp, float, pRangeProfileClutterRemoved_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_2->heartWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_3(pVitalSigns_Breath_CircularBuffer, float,
                       pDataOutTemp_vs2_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_3->circularBufferSizeBreath);

    MMW_ALLOC_BUF_VS_3(pVitalSigns_Heart_CircularBuffer, float,
                       pVitalSigns_Breath_CircularBuffer_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_3->circularBufferSizeHeart);

    MMW_ALLOC_BUF_VS_3(pMotionCircularBuffer, float,
                       pVitalSigns_Heart_CircularBuffer_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_3->motionDetection_BlockSize);

    MMW_ALLOC_BUF_VS_3(pDopplerWindow, float, pMotionCircularBuffer_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       DOPPLER_WINDOW_SIZE);

    MMW_ALLOC_BUF_VS_3(pVitalSigns_SpectrumCplx, cmplx32ReIm_t,
                       pDopplerWindow_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_3->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_3(pVitalSignsBuffer_Cplx, cmplx32ReIm_t,
                       pVitalSigns_SpectrumCplx_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_3->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_3(pVitalSignsSpectrumTwiddle32x32, cmplx32ReIm_t,
                       pVitalSignsBuffer_Cplx_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_3->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_3(pVitalSigns_Breath_AbsSpectrum, float,
                       pVitalSignsSpectrumTwiddle32x32_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_3->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_3(pVitalSigns_Heart_AbsSpectrum, float,
                       pVitalSigns_Breath_AbsSpectrum_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_3->heartWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_3(
            pFilterCoefsBreath, float, pVitalSigns_Heart_AbsSpectrum_vs3_end,
            MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
            IIR_FILTER_BREATH_NUM_STAGES * IIR_FILTER_COEFS_SECOND_ORDER);

    MMW_ALLOC_BUF_VS_3(pScaleValsBreath, float, pFilterCoefsBreath_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       IIR_FILTER_BREATH_NUM_STAGES + 1);

    MMW_ALLOC_BUF_VS_3(
            pFilterCoefsHeart_4Hz, float, pScaleValsBreath_vs3_end,
            MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
            IIR_FILTER_HEART_NUM_STAGES * IIR_FILTER_COEFS_SECOND_ORDER);

    MMW_ALLOC_BUF_VS_3(pScaleValsHeart_4Hz, float,
                       pFilterCoefsHeart_4Hz_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       IIR_FILTER_HEART_NUM_STAGES + 1);

    MMW_ALLOC_BUF_VS_3(pXcorr, float, pScaleValsHeart_4Hz_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN, XCORR_NUM_LAGS);

    MMW_ALLOC_BUF_VS_3(pTempReal_Prev, float, pXcorr_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj->numRangeBins);

    MMW_ALLOC_BUF_VS_3(pTempImag_Prev, float, pTempReal_Prev_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj->numRangeBins);

    MMW_ALLOC_BUF_VS_3(pRangeProfileClutterRemoved, float,
                       pTempImag_Prev_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj->numRangeBins);

    MMW_ALLOC_BUF_VS_3(pDataOutTemp, float, pRangeProfileClutterRemoved_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_3->heartWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_4(pVitalSigns_Breath_CircularBuffer, float,
                       pDataOutTemp_vs3_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_2->circularBufferSizeBreath);

    MMW_ALLOC_BUF_VS_4(pVitalSigns_Heart_CircularBuffer, float,
                       pVitalSigns_Breath_CircularBuffer_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_4->circularBufferSizeHeart);

    MMW_ALLOC_BUF_VS_4(pMotionCircularBuffer, float,
                       pVitalSigns_Heart_CircularBuffer_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_4->motionDetection_BlockSize);

    MMW_ALLOC_BUF_VS_4(pDopplerWindow, float, pMotionCircularBuffer_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       DOPPLER_WINDOW_SIZE);

    MMW_ALLOC_BUF_VS_4(pVitalSigns_SpectrumCplx, cmplx32ReIm_t,
                       pDopplerWindow_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_4->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_4(pVitalSignsBuffer_Cplx, cmplx32ReIm_t,
                       pVitalSigns_SpectrumCplx_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_4->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_4(pVitalSignsSpectrumTwiddle32x32, cmplx32ReIm_t,
                       pVitalSignsBuffer_Cplx_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_4->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_4(pVitalSigns_Breath_AbsSpectrum, float,
                       pVitalSignsSpectrumTwiddle32x32_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_4->breathingWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_4(pVitalSigns_Heart_AbsSpectrum, float,
                       pVitalSigns_Breath_AbsSpectrum_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj_VS_4->heartWfm_Spectrum_FftSize);

    MMW_ALLOC_BUF_VS_4(
            pFilterCoefsBreath, float, pVitalSigns_Heart_AbsSpectrum_vs4_end,
            MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
            IIR_FILTER_BREATH_NUM_STAGES * IIR_FILTER_COEFS_SECOND_ORDER);

    MMW_ALLOC_BUF_VS_4(pScaleValsBreath, float, pFilterCoefsBreath_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       IIR_FILTER_BREATH_NUM_STAGES + 1);

    MMW_ALLOC_BUF_VS_4(
            pFilterCoefsHeart_4Hz, float, pScaleValsBreath_vs4_end,
            MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
            IIR_FILTER_HEART_NUM_STAGES * IIR_FILTER_COEFS_SECOND_ORDER);

    MMW_ALLOC_BUF_VS_4(pScaleValsHeart_4Hz, float,
                       pFilterCoefsHeart_4Hz_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       IIR_FILTER_HEART_NUM_STAGES + 1);

    MMW_ALLOC_BUF_VS_4(pXcorr, float, pScaleValsHeart_4Hz_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN, XCORR_NUM_LAGS);

    MMW_ALLOC_BUF_VS_4(pTempReal_Prev, float, pXcorr_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj->numRangeBins);

    MMW_ALLOC_BUF_VS_4(pTempImag_Prev, float, pTempReal_Prev_vs4_end,
                       MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
                       obj->numRangeBins);

    MMW_ALLOC_VS4_NO_END(pRangeProfileClutterRemoved, float,
                         pTempImag_Prev_vs4_end,
                         MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN);

//    MMW_ALLOC_BUF_VS_2(pDataOutTemp, float,
//                  pRangeProfileClutterRemoved_vs2_end, MMWDEMO_MEMORY_ALLOC_DOUBLE_WORD_ALIGN,
//                  obj_VS_2->heartWfm_Spectrum_FftSize);

//////////////  Parameters :  IIR-Cascade Bandpass Filter  //////////////////////////////////////

    // 0.1 Hz to 0.5 Hz Bandpass Filter coefficients
    float pFilterCoefsBreath[IIR_FILTER_BREATH_NUM_STAGES
            * IIR_FILTER_COEFS_SECOND_ORDER] = { 1.0000, 0, -1.0000, 1.0000,
                                                 -1.9196, 0.9252, 1.0000, 0,
                                                 -1.0000, 1.0000, -1.6624,
                                                 0.7390 };
    float pScaleValsBreath[IIR_FILTER_BREATH_NUM_STAGES + 1] = { 0.0602, 0.0602,
                                                                 1.0000 };
    memcpy(vitalsigns_dataPathObj_zone_1.pFilterCoefsBreath, pFilterCoefsBreath,
           sizeof(pFilterCoefsBreath));
    memcpy(vitalsigns_dataPathObj_zone_1.pScaleValsBreath, pScaleValsBreath,
           sizeof(pScaleValsBreath));
    memcpy(vitalsigns_dataPathObj_zone_2.pFilterCoefsBreath, pFilterCoefsBreath,
           sizeof(pFilterCoefsBreath));
    memcpy(vitalsigns_dataPathObj_zone_2.pScaleValsBreath, pScaleValsBreath,
           sizeof(pScaleValsBreath));
    memcpy(vitalsigns_dataPathObj_zone_3.pFilterCoefsBreath, pFilterCoefsBreath,
           sizeof(pFilterCoefsBreath));
    memcpy(vitalsigns_dataPathObj_zone_3.pScaleValsBreath, pScaleValsBreath,
           sizeof(pScaleValsBreath));
    memcpy(vitalsigns_dataPathObj_zone_4.pFilterCoefsBreath, pFilterCoefsBreath,
           sizeof(pFilterCoefsBreath));
    memcpy(vitalsigns_dataPathObj_zone_4.pScaleValsBreath, pScaleValsBreath,
           sizeof(pScaleValsBreath));

    // Heart Beat Rate    0.8 - 4.0 Hz
    float pFilterCoefsHeart_4Hz[IIR_FILTER_HEART_NUM_STAGES
            * IIR_FILTER_COEFS_SECOND_ORDER] = { 1.0000, 0, -1.0000, 1.0000,
                                                 -1.4522, 0.6989, 1.0000, 0,
                                                 -1.0000, 1.0000, 1.5573,
                                                 0.7371, 1.0000, 0, -1.0000,
                                                 1.0000, 1.2189, 0.3932, 1.0000,
                                                 0, -1.0000, 1.0000, -1.0947,
                                                 0.3264 };
    float pScaleValsHeart_4Hz[IIR_FILTER_HEART_NUM_STAGES + 1] = { 0.4188,
                                                                   0.4188,
                                                                   0.3611,
                                                                   0.3611,
                                                                   1.0000 };
    memcpy(vitalsigns_dataPathObj_zone_1.pFilterCoefsHeart_4Hz,
           pFilterCoefsHeart_4Hz, sizeof(pFilterCoefsHeart_4Hz));
    memcpy(vitalsigns_dataPathObj_zone_1.pScaleValsHeart_4Hz,
           pScaleValsHeart_4Hz, sizeof(pScaleValsHeart_4Hz));
    memcpy(vitalsigns_dataPathObj_zone_2.pFilterCoefsHeart_4Hz,
           pFilterCoefsHeart_4Hz, sizeof(pFilterCoefsHeart_4Hz));
    memcpy(vitalsigns_dataPathObj_zone_2.pScaleValsHeart_4Hz,
           pScaleValsHeart_4Hz, sizeof(pScaleValsHeart_4Hz));
    memcpy(vitalsigns_dataPathObj_zone_3.pFilterCoefsHeart_4Hz,
           pFilterCoefsHeart_4Hz, sizeof(pFilterCoefsHeart_4Hz));
    memcpy(vitalsigns_dataPathObj_zone_3.pScaleValsHeart_4Hz,
           pScaleValsHeart_4Hz, sizeof(pScaleValsHeart_4Hz));
    memcpy(vitalsigns_dataPathObj_zone_4.pFilterCoefsHeart_4Hz,
           pFilterCoefsHeart_4Hz, sizeof(pFilterCoefsHeart_4Hz));
    memcpy(vitalsigns_dataPathObj_zone_4.pScaleValsHeart_4Hz,
           pScaleValsHeart_4Hz, sizeof(pScaleValsHeart_4Hz));

    memset(vitalsigns_dataPathObj_zone_1.pDelayHeart, 0,
           sizeof(vitalsigns_dataPathObj_zone_1.pDelayHeart));
    memset(vitalsigns_dataPathObj_zone_1.pDelayBreath, 0,
           sizeof(vitalsigns_dataPathObj_zone_1.pDelayBreath));
    memset(vitalsigns_dataPathObj_zone_2.pDelayHeart, 0,
           sizeof(vitalsigns_dataPathObj_zone_2.pDelayHeart));
    memset(vitalsigns_dataPathObj_zone_2.pDelayBreath, 0,
           sizeof(vitalsigns_dataPathObj_zone_2.pDelayBreath));
    memset(vitalsigns_dataPathObj_zone_3.pDelayHeart, 0,
           sizeof(vitalsigns_dataPathObj_zone_3.pDelayHeart));
    memset(vitalsigns_dataPathObj_zone_3.pDelayBreath, 0,
           sizeof(vitalsigns_dataPathObj_zone_3.pDelayBreath));
    memset(vitalsigns_dataPathObj_zone_4.pDelayHeart, 0,
           sizeof(vitalsigns_dataPathObj_zone_4.pDelayHeart));
    memset(vitalsigns_dataPathObj_zone_4.pDelayBreath, 0,
           sizeof(vitalsigns_dataPathObj_zone_4.pDelayBreath));

    float DopplerWindowCoefs[DOPPLER_WINDOW_SIZE] = { 0.0800, 0.0894, 0.1173,
                                                      0.1624, 0.2231, 0.2967,
                                                      0.3802, 0.4703, 0.5633,
                                                      0.6553, 0.7426, 0.8216,
                                                      0.8890, 0.9422, 0.9789,
                                                      0.9976 };
    memcpy(vitalsigns_dataPathObj_zone_1.pDopplerWindow, DopplerWindowCoefs,
           sizeof(DopplerWindowCoefs));
    memcpy(vitalsigns_dataPathObj_zone_2.pDopplerWindow, DopplerWindowCoefs,
           sizeof(DopplerWindowCoefs));
    memcpy(vitalsigns_dataPathObj_zone_3.pDopplerWindow, DopplerWindowCoefs,
           sizeof(DopplerWindowCoefs));
    memcpy(vitalsigns_dataPathObj_zone_4.pDopplerWindow, DopplerWindowCoefs,
           sizeof(DopplerWindowCoefs));

    float pFilterCoefs[FIR_FILTER_SIZE] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
                                            0.1, 0.1, 0.1 };
    memcpy(vitalsigns_dataPathObj_zone_1.pFilterCoefs, pFilterCoefs,
    FIR_FILTER_SIZE * sizeof(float));
    memcpy(vitalsigns_dataPathObj_zone_2.pFilterCoefs, pFilterCoefs,
    FIR_FILTER_SIZE * sizeof(float));
    memcpy(vitalsigns_dataPathObj_zone_3.pFilterCoefs, pFilterCoefs,
    FIR_FILTER_SIZE * sizeof(float));
    memcpy(vitalsigns_dataPathObj_zone_4.pFilterCoefs, pFilterCoefs,
    FIR_FILTER_SIZE * sizeof(float));

    heapUsed = radarCube_end - heapL3start;

    MmwDemo_printHeapStats("L3", heapUsed, sizeof(gMmwL3));
}

void MmwDemo_dataPathConfigFFTs(MmwDemo_DSS_DataPathObj *obj)
{

    MmwDemo_genWindow((void*) obj->window1D,
    FFT_WINDOW_INT16,
                      obj->numAdcSamples, obj->numAdcSamples / 2,
                      ONE_Q15,
                      MMW_WIN_BLACKMAN);

    /* Generate twiddle factors for 1D FFT. This is one time */
    MmwDemo_gen_twiddle_fft16x16_fast((int16_t*) obj->twiddle16x16_1D,
                                      obj->numRangeBins);

    /* Generate twiddle factors for vital signs FFT operation. This is one time */
    //MmwDemo_gen_twiddle_fft32x32_fast((int32_t *)vitalsigns_dataPathObj.pVitalSignsSpectrumTwiddle32x32, vitalsigns_dataPathObj.breathingWfm_Spectrum_FftSize, 2147483647.5);
    MmwDemo_gen_twiddle_fft32x32_fast(
            (int32_t*) vitalsigns_dataPathObj_zone_1.pVitalSignsSpectrumTwiddle32x32,
            vitalsigns_dataPathObj_zone_1.breathingWfm_Spectrum_FftSize,
            2147483647.5);
    MmwDemo_gen_twiddle_fft32x32_fast(
            (int32_t*) vitalsigns_dataPathObj_zone_2.pVitalSignsSpectrumTwiddle32x32,
            vitalsigns_dataPathObj_zone_2.breathingWfm_Spectrum_FftSize,
            2147483647.5);
    MmwDemo_gen_twiddle_fft32x32_fast(
            (int32_t*) vitalsigns_dataPathObj_zone_3.pVitalSignsSpectrumTwiddle32x32,
            vitalsigns_dataPathObj_zone_3.breathingWfm_Spectrum_FftSize,
            2147483647.5);
    MmwDemo_gen_twiddle_fft32x32_fast(
            (int32_t*) vitalsigns_dataPathObj_zone_4.pVitalSignsSpectrumTwiddle32x32,
            vitalsigns_dataPathObj_zone_4.breathingWfm_Spectrum_FftSize,
            2147483647.5);
}

/**
 *  @b Description
 *  @n
 *      Function generates a single FFT window samples. It calls single precision
 *      sine and cosine functions from mathlib library for the first sample only, and then recursively
 *      generates cosine values for other samples.
 *
 *  @param[out] win             Pointer to calculated window samples.
 *  @param[in]  windowDatumType Window samples data format. For windowDatumType = @ref FFT_WINDOW_INT16,
 *              the samples format is int16_t. For windowDatumType = @ref FFT_WINDOW_INT32,
 *              the samples format is int32_t.
 *  @param[in]  winLen          Nominal window length
 *  @param[in]  winGenLen       Number of generated samples
 *  @param[in]  oneQformat      Q format of samples, oneQformat is the value of
 *                              one in the desired format.
 *  @param[in]  winType         Type of window, one of @ref MMW_WIN_BLACKMAN, @ref MMW_WIN_HANNING,
 *              or @ref MMW_WIN_RECT.
 *  @retval none.
 */
void MmwDemo_genWindow(void *win, uint32_t windowDatumType, uint32_t winLen,
                       uint32_t winGenLen, int32_t oneQformat, uint32_t winType)
{
    uint32_t winIndx;
    int32_t winVal;
    int16_t *win16 = (int16_t*) win;
    int32_t *win32 = (int32_t*) win;

    float eR = 1.;
    float eI = 0.;
    float e2R = 1.;
    float e2I = 0.;
    float ephyR, ephyI;
    float e2phyR, e2phyI;
    float tmpR;

    float phi = 2 * PI_ / ((float) winLen - 1);

    ephyR = cossp(phi);
    ephyI = sinsp(phi);

    e2phyR = ephyR * ephyR - ephyI * ephyI;
    e2phyI = 2 * ephyR * ephyI;

    if (winType == MMW_WIN_BLACKMAN)
    {
        //Blackman window
        float a0 = 0.42;
        float a1 = 0.5;
        float a2 = 0.08;
        for (winIndx = 0; winIndx < winGenLen; winIndx++)
        {
            winVal = (int32_t) ((oneQformat * (a0 - a1 * eR + a2 * e2R)) + 0.5);
            if (winVal >= oneQformat)
            {
                winVal = oneQformat - 1;
            }
            if (windowDatumType == FFT_WINDOW_INT16)
            {
                win16[winIndx] = (int16_t) winVal;
            }
            if (windowDatumType == FFT_WINDOW_INT32)
            {
                win32[winIndx] = (int32_t) winVal;
            }
            tmpR = eR;
            eR = eR * ephyR - eI * ephyI;
            eI = tmpR * ephyI + eI * ephyR;

            tmpR = e2R;
            e2R = e2R * e2phyR - e2I * e2phyI;
            e2I = tmpR * e2phyI + e2I * e2phyR;
        }
    }
    else if (winType == MMW_WIN_HANNING)
    {
        //Hanning window
        for (winIndx = 0; winIndx < winGenLen; winIndx++)
        {
            winVal = (int32_t) ((oneQformat * 0.5 * (1 - eR)) + 0.5);
            if (winVal >= oneQformat)
            {
                winVal = oneQformat - 1;
            }
            if (windowDatumType == FFT_WINDOW_INT16)
            {
                win16[winIndx] = (int16_t) winVal;
            }
            if (windowDatumType == FFT_WINDOW_INT32)
            {
                win32[winIndx] = (int32_t) winVal;
            }
            tmpR = eR;
            eR = eR * ephyR - eI * ephyI;
            eI = tmpR * ephyI + eI * ephyR;
        }
    }
    else if (winType == MMW_WIN_RECT)
    {
        //Rectangular window
        for (winIndx = 0; winIndx < winGenLen; winIndx++)
        {
            if (windowDatumType == FFT_WINDOW_INT16)
            {
                win16[winIndx] = (int16_t) (oneQformat - 1);
            }
            if (windowDatumType == FFT_WINDOW_INT32)
            {
                win32[winIndx] = (int32_t) (oneQformat - 1);
            }
        }
    }
}

void MmwDemo_checkDynamicConfigErrors(MmwDemo_DSS_DataPathObj *obj)
{
    MmwDemo_CliCfg_t *cliCfg = obj->cliCfg;

    MmwDemo_dssAssert(
            !((cliCfg->extendedMaxVelocityCfg.enabled == 1)
                    && (cliCfg->multiObjBeamFormingCfg.enabled == 1)));

    MmwDemo_dssAssert(
            !((cliCfg->extendedMaxVelocityCfg.enabled == 1)
                    && (cliCfg->nearFieldCorrectionCfg.enabled == 1)));

    MmwDemo_dssAssert(
            !((cliCfg->extendedMaxVelocityCfg.enabled == 1)
                    && (obj->numTxAntennas == 1)));
}
