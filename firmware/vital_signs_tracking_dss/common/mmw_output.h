/**
 * mmw_output.h
 *
 * This is the interface/message header file for the Millimeter Wave Demo
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

#ifndef MMW_OUTPUT_H
#define MMW_OUTPUT_H

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Output packet length is a multiple of this value, must be power of 2*/
#define MMWDEMO_OUTPUT_MSG_SEGMENT_LEN 32

/*!
 * @brief
 *  Message types used in Millimeter Wave Demo for the communication between
 *  target and host, and also for Mailbox communication
 *  between MSS and DSS on the XWR16xx platform. Message types are used to indicate
 *  different type detection information sent out from the target.
 *
 */
typedef enum MmwDemo_output_message_type_e
{
    /*! @brief   List of detected points */
    MMWDEMO_OUTPUT_MSG_DETECTED_POINTS = 1,

    /*! @brief   Range profile */
    MMWDEMO_OUTPUT_MSG_RANGE_PROFILE,

    /*! @brief   Noise floor profile */
    MMWDEMO_OUTPUT_MSG_NOISE_PROFILE,

    /*! @brief   Samples to calculate static azimuth  heatmap */
    MMWDEMO_OUTPUT_MSG_AZIMUT_STATIC_HEAT_MAP,

    /*! @brief   Range/Doppler detection matrix */
    MMWDEMO_OUTPUT_MSG_RANGE_DOPPLER_HEAT_MAP,

    /*! @brief   Stats information */
    MMWDEMO_OUTPUT_MSG_STATS,

    ODDEMO_OUTPUT_MSG_RANGE_AZIMUT_HEAT_MAP = 8,

    ODDEMO_OUTPUT_MSG_DECISION,

    VS_OUTPUT_HEART_BREATHING_RATES,

    ODDEMO_OUTPUT_MSG_ROW_NOISE,

    PEAK_POSITIONS,

    MMWDEMO_OUTPUT_MSG_MAX
} MmwDemo_output_message_type;

/*!
 * @brief
 *  Message header for reporting detection information from data path.
 *
 * @details
 *  The structure defines the message header.
 */
typedef struct MmwDemo_output_message_header_t
{
    /*! @brief   Output buffer magic word (sync word). It is initialized to  {0x0102,0x0304,0x0506,0x0708} */
    uint16_t    magicWord[4];

    /*! @brief   Total packet length including header in Bytes */
    uint32_t    totalPacketLen;

    /*! @brief   platform type */
    uint32_t    platform;

    /*! @brief   Frame number */
    uint32_t    frameNumber;

    /*! @brief   Time in CPU cycles when the message was created. For XWR16xx: DSP CPU cycles, for XWR14xx: R4F CPU cycles */
    uint32_t    timeCpuCycles;

    /*! @brief   Number of detected objects */
    uint32_t    numDetectedObj;

    /*! @brief   Number of TLVs */
    uint32_t    numTLVs;

} MmwDemo_output_message_header;


/**
 * @brief
 *  Message for reporting detected objects from data path.
 *
 * @details
 *  The structure defines the message body for detected objects from from data path.
 */
typedef struct MmwDemo_output_message_tl_t
{
    /*! @brief   TLV type */
    uint32_t    type;

    /*! @brief   Length in bytes */
    uint32_t    length;

} MmwDemo_output_message_tl;

typedef struct VitalSignsDemo_OutputStats_t
{
    uint16_t rangeBinIndexMax;     // 1
    uint16_t rangeBinIndexPhase;   // 1
    float maxVal;                  //2
    uint32_t processingCyclesOut;  //3
    uint16_t rangeBinStartIndex;   //4
    uint16_t rangeBinEndIndex;     //4
    float unwrapPhasePeak_mm;           // 5
    float outputFilterBreathOut;        // 6
    float outputFilterHeartOut;         // 7
    float heartRateEst_FFT;             // 8
    float heartRateEst_FFT_4Hz;         // 9
    float heartRateEst_xCorr;           // 10
    float heartRateEst_peakCount_filtered;  // 11
    float breathingRateEst_FFT;          // 12
    float breathingRateEst_xCorr;        // 13
    float breathingRateEst_peakCount;    // 14
    float confidenceMetricBreathOut;        // 15
    float confidenceMetricBreathOut_xCorr;  // 16
    float confidenceMetricHeartOut;         // 17
    float confidenceMetricHeartOut_4Hz;     // 18
    float confidenceMetricHeartOut_xCorr;   // 19
    float sumEnergyBreathWfm;               // 20
    float sumEnergyHeartWfm;                // 21
    float motionDetectedFlag;               // 22
    float breathingRateEst_harmonicEnergy;  // 23
    float heartRateEst_harmonicEnergy;      // 24
    float reserved7;  //25
    float reserved8;  //26
    float reserved9;  //27
    float reserved10;  //28
    float reserved11;  //29
    float reserved12;  //30
    float reserved13;  //31
    float reserved14;  //32
} VitalSignsDemo_OutputStats;


#ifdef __cplusplus
}
#endif

#endif /* MMW_OUTPUT_H */
