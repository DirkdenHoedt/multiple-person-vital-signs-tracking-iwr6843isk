/**
 * oddemo_common.h
 *
 * Global constant and structure definitions for Occupancy Detection
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

#ifndef ODDEMO_COMMON
#define ODDEMO_COMMON

//#include "swpform.h"
#ifdef _TMS320C6X
#include "c6x.h"
#endif

#define ODDEMO_MAX_RANGE               64
#define ODDEMO_MAX_AZIMUTH             48
#define ODDEMO_MAX_ZONES               4
#define ODDEMO_ZONE_PAIR              (ODDEMO_MAX_ZONES / 2)
#define ODDEMO_MAX_FRAME_HIST          16
#define ODDEMO_MATRIX_ROW_SIZE         6
#define ODDEMO_MATRIX_NUM_ROWS         4
#define ODDEMO_MATRIX_SIZE            (ODDEMO_MATRIX_ROW_SIZE * ODDEMO_MATRIX_NUM_ROWS)
#define ODDEMO_ANGLE_RANGE             60.0
#define ODDEMO_ANGLE_RESOLUTION       ((ODDEMO_ANGLE_RANGE * 2) / (float)ODDEMO_MAX_AZIMUTH)
#define ODDEMO_ROW_NOISE_FRAMES        64

#define ODDEMO_ANGLEHEATMAP_BUF_SIZE   (ODDEMO_MAX_RANGE * ODDEMO_MAX_AZIMUTH)

// SteeringVec Size (bytes): (nRxAnt-1) * steeringVecSize * sizeof(cplxf_t); For Ant0, since steering vector is = 1, it is not saved.
// steeringVecSize =   (uint32_t) ((2.f * estAngleRange) / estAngleResolution) + 0.5f);
// For estAngleRange = 60 (-60, 60) and 64 nAzimuthBins (angles), the resolution = 120/64 = 1.875
// steeringVecSize = 64

#define ODDEMO_STEERINGVEC_L1_BUF_SIZE (ODDEMO_MAX_AZIMUTH * 8)
#define ODDEMO_SCRATCH_L1_BUF_SIZE      360 // for nRxAnt=8 for uint32_t

// size = [(nRxAnt * ( 1 + nRxAnt))/2] * sizeof(cplxf_t)
// upper triangle of nRxAnt x nRxAnt Hermitian matrix
#define ODDEMO_INVRNMATRIX_BUF_SIZE     36  // for nRxAnt=8 for cplxf_t
#define ODDEMO_INVRNMATRIX_BUFFER      (ODDEMO_INVRNMATRIX_BUF_SIZE * ODDEMO_MAX_RANGE)


//This type defines the position and size of an occupancy zone in the heat map.
typedef struct
{
  uint16_t range_start;
  uint16_t range_length;
  uint16_t azimuth_start;
  uint16_t azimuth_width;
} ODDEMO_Zone;

//This type defines the feature vector contents for a pair of zones.
typedef struct
{
  float powerMA[2];
  float powRatio[2];
  float crossCorr;
} ODDEMO_Feature;

//This type defines the feature vector contents for Vital Signs information
typedef struct
{
  float unwrapped_waveform[2];
  float heart_waveform[2];
  float breathing_waveform[2];
  float heart_rate[2];
  float breathing_rate[2];
} VS_Feature;

//This type defines the demo's configurable parameters.
typedef struct
{
    float    gamma;
    uint32_t windowLen;
} ODDEMO_Parms;

#endif //ODDEMO_COMMON
