/**
 * oddemo_heatmap.h
 *
 * Header file relating to the Occupancy Detection Heatmap
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
#ifndef ODDEMO_HEATMAP
#define ODDEMO_HEATMAP

#include "dss_data_path.h"

#define ODDEMO_HEATMAP_MAX16 32000
#define ODDEMO_HEATMAP_MAX8  252

typedef struct
{
    float min;
    float max;
    float avg;
} ODDemo_Heatmap_Stats;


extern void ODDemo_Heatmap_aoaEstCaponBF_range(ODDemo_DataPathObj *obj);

extern void ODDemo_Heatmap_aoaEstCaponBF_doppler(ODDemo_DataPathObj *obj, uint8_t rangeIndx, uint8_t azimuthIndx);

extern void ODDemo_Heatmap_aoaEstCaponBF_heatmap(uint8_t  bfFlag,
                                                 int32_t  nRxAnt,
                                                 int32_t  numAzimuthBins,
                                                 cplxf_t *steeringVec,
                                                 cplxf_t *invRnMatrices,
                                                 float   *rangeAzimuthHeatMap);


extern void ODDemo_Heatmap_aoaEstCaponBF_covInv(uint8_t   invFlag,
                                                uint8_t   clutterRmFlag,
                                                float     gamma,
                                                int32_t   nRxAnt,
                                                int32_t   nChirps,
                                                int32_t   log2Chirps,
                                                int32_t  *scratch,
                                                cplx16_t *inputAntSamples,
                                                cplxf_t  *invRnMatrices);

extern void ODDEMO_Heatmap_aoaEstCaponBF_dopplerEstInput(uint8_t bfFlag,
                                                        int32_t nRxAnt,
                                                        int32_t nChirps,
                                                        cplx16_t * inputAntSamples,
                                                        cplxf_t * steeringVec,
                                                        cplxf_t  * invRnMatrices,
                                                        int32_t * scratch,
                                                        float rangeAzimuthHeatMap,
                                                        float * bfOutput);

extern void ODDemo_Heatmap_steeringVecGen(ODDemo_DataPathObj *obj);

extern void ODDemo_Heatmap_get_stats(float *heatmap);

extern void ODDemo_Heatmap_arc_removal(float *heatmap);

extern void ODDemo_Heatmap_scale_heatmap16(float *heatin, uint16_t *heatout);

extern void ODDemo_Heatmap_scale_heatmap8(float *heatin, uint8_t *heatout);

void heatmapGetPeaks(float *heatmap);

void addPeakToArray(uint8_t rng_idx, uint8_t az_idx, float power,
                    float *maxPeaks, uint8_t *numPeaks);

#endif //ODDEMO_HEATMAP
