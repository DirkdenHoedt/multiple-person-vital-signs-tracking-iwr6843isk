/**
 * cli.c
 *
 * Mmw (Milli-meter wave) DEMO CLI Implementation
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

/* mmWave SDK Include Files: */
#include <ti/common/sys_common.h>
#include <ti/common/mmwave_sdk_version.h>
#include <ti/drivers/uart/UART.h>
#include <ti/control/mmwavelink/mmwavelink.h>
#include <ti/utils/cli/cli.h>

/* Demo Include Files */
#include "./common/oddemo_common.h"
#include "./common/mmw_messages.h"
#include "./mss_mmw.h"

float       oddemo_coeffMatrix[ODDEMO_ZONE_PAIR][ODDEMO_MATRIX_SIZE];
ODDEMO_Zone oddemo_zone[ODDEMO_MAX_ZONES];
float       oddemo_meanVec[ODDEMO_ZONE_PAIR][ODDEMO_MATRIX_ROW_SIZE-1];
float       oddemo_stdVec[ODDEMO_ZONE_PAIR][ODDEMO_MATRIX_ROW_SIZE-1];
uint16_t    oddemo_num_zones;
float       oddemo_row_noise[ODDEMO_MAX_RANGE];

/**************************************************************************
 *************************** Local Definitions ****************************
 **************************************************************************/

/* CLI Extended Command Functions */
static int32_t MmwDemo_CLICalibDcRangeSig (int32_t argc, char* argv[]);
static int32_t MmwDemo_CLIClutterRemoval (int32_t argc, char* argv[]);
static int32_t MmwDemo_CLISensorStart (int32_t argc, char* argv[]);
static int32_t MmwDemo_CLISensorStop (int32_t argc, char* argv[]);
static int32_t MmwDemo_CLIADCBufCfg (int32_t argc, char* argv[]);

/* ODDEMO Specific Extended Command Functions */
static int32_t ODDemo_CLIGuiMonSel(int32_t argc, char* argv[]);
static int32_t ODDemo_CLIZoneDef(int32_t argc, char *argv[]);
static int32_t ODDemo_CLICoeffMatrix(int32_t argc, char *argv[]);
static int32_t ODDemo_CLIMeanVector(int32_t argc, char *argv[]);
static int32_t ODDemo_CLIStdVector(int32_t argc, char *argv[]);
static int32_t ODDemo_CLIParms(int32_t argc, char *argv[]);
static int32_t ODDemo_CLIRowNoise(int32_t argc, char *argv[]);

/* Vital Signs Specific Extended Command Functions */
static int32_t VitalSignsDemo_CLIGuiMonSel (int32_t argc, char* argv[]);
static int32_t VitalSignsDemo_CLIvitalSignsParamsCfg (int32_t argc, char* argv[]);
static int32_t VitalSignsDemo_CLIMotionDetection (int32_t argc, char* argv[]);

/**************************************************************************
 *************************** Extern Definitions *******************************
 **************************************************************************/

extern MmwDemo_MCB    gMmwMssMCB;
extern int32_t MmwDemo_mboxWrite(MmwDemo_message     * message);

/**************************************************************************
 *************************** CLI  Function Definitions **************************
 **************************************************************************/

/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for the sensor start command
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
static int32_t MmwDemo_CLISensorStart (int32_t argc, char* argv[])
{
    bool doReconfig = true;
    if (argc==2)
    {
        doReconfig = (bool) atoi (argv[1]);
    }
    /* Post sensorSTart event to notify configuration is done */
    MmwDemo_notifySensorStart(doReconfig);
    /* Pend for completion */
    return (MmwDemo_waitSensorStartComplete());
}


/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for the sensor stop command
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
static int32_t MmwDemo_CLISensorStop (int32_t argc, char* argv[])
{
    /* Post sensorSTOP event to notify sensor stop command */
    MmwDemo_notifySensorStop();
    /* Pend for completion */
    MmwDemo_waitSensorStopComplete();
    return 0;
}

static int32_t MmwDemo_CLIGetSubframe (int32_t argc, char* argv[], int32_t expectedArgc, int8_t* subFrameNum)
{
    int8_t subframe;

    /* Sanity Check: Minimum argument check */
    if (argc != expectedArgc)
    {
        CLI_write ("Error: Invalid usage of the CLI command\n");
        return -1;
    }

    /*Subframe info is always in position 1*/
    subframe = (int8_t) atoi(argv[1]);

    if(subframe >= (int8_t)RL_MAX_SUBFRAMES)
    {
        CLI_write ("Error: Subframe number is invalid\n");
        return -1;
    }

    *subFrameNum = (int8_t)subframe;

    return 0;
}

static void MmwDemo_mssCfgUpdate(void *srcPtr, uint32_t offset, uint32_t size, int8_t subFrameNum)
{
    /* if subFrameNum undefined, broadcast to all sub-frames */
    if(subFrameNum == MMWDEMO_SUBFRAME_NUM_FRAME_LEVEL_CONFIG)
    {
        uint8_t  indx;
        for(indx = 0; indx < RL_MAX_SUBFRAMES; indx++)
        {
            memcpy((void *)((uint32_t) &gMmwMssMCB.cliCfg[indx] + offset), srcPtr, size);
        }

    }
    else
    {
        /* Apply configuration to specific subframe (or to position zero for the legacy case
           where there is no advanced frame config) */
        memcpy((void *)((uint32_t) &gMmwMssMCB.cliCfg[subFrameNum] + offset), srcPtr, size);
    }
}


uint32_t log2Approx(uint32_t x)
{
    uint32_t idx, detectFlag = 0;

    if ( x < 2)
    {
        return (0);
    }

    idx = 32U;
    while((detectFlag==0U) || (idx==0U))
    {
        if(x & 0x80000000U)
        {
            detectFlag = 1;
        }
        x <<= 1U;
        idx--;
    }

    if(x != 0)
    {
        idx = idx + 1;
    }

    return(idx);
}

/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for DC range calibration
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
static int32_t MmwDemo_CLICalibDcRangeSig (int32_t argc, char* argv[])
{
    MmwDemo_CalibDcRangeSigCfg cfg;
    MmwDemo_message            message;
    uint32_t                   log2NumAvgChirps;
    int8_t                     subFrameNum;

    if(MmwDemo_CLIGetSubframe(argc, argv, 6, &subFrameNum) < 0)
    {
        return -1;
    }

    /* Initialize configuration for DC range signature calibration */
    memset ((void *)&cfg, 0, sizeof(MmwDemo_CalibDcRangeSigCfg));

    /* Populate configuration: */
    cfg.enabled          = (uint16_t) atoi (argv[2]);
    cfg.negativeBinIdx   = (int16_t)  atoi (argv[3]);
    cfg.positiveBinIdx   = (int16_t)  atoi (argv[4]);
    cfg.numAvgChirps     = (uint16_t)  atoi (argv[5]);

    if (cfg.negativeBinIdx > 0)
    {
        CLI_write ("Error: Invalid negative bin index\n");
        return -1;
    }
    if ((cfg.positiveBinIdx - cfg.negativeBinIdx + 1) > DC_RANGE_SIGNATURE_COMP_MAX_BIN_SIZE)
    {
        CLI_write ("Error: Number of bins exceeds the limit\n");
        return -1;
    }
    log2NumAvgChirps = (uint32_t) log2Approx (cfg.numAvgChirps);
    if (cfg.numAvgChirps != (1 << log2NumAvgChirps))
    {
        CLI_write ("Error: Number of averaged chirps is not power of two\n");
        return -1;
    }

    /* Save Configuration to use later */
    MmwDemo_mssCfgUpdate((void *)&cfg, offsetof(MmwDemo_CliCfg_t, calibDcRangeSigCfg),
        sizeof(MmwDemo_CalibDcRangeSigCfg), subFrameNum);

    /* Send configuration to DSS */
    memset((void *)&message, 0, sizeof(MmwDemo_message));

    message.type = MMWDEMO_MSS2DSS_CALIB_DC_RANGE_SIG;
    message.subFrameNum = subFrameNum;
    memcpy((void *)&message.body.calibDcRangeSigCfg, (void *)&cfg, sizeof(MmwDemo_CalibDcRangeSigCfg));

    if (MmwDemo_mboxWrite(&message) == 0)
        return 0;
    else
        return -1;
}


/**
 *  @b Description
 *  @n
 *      Clutter removal Configuration
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
static int32_t MmwDemo_CLIClutterRemoval (int32_t argc, char* argv[])
{
    MmwDemo_ClutterRemovalCfg cfg;
    MmwDemo_message     message;
    int8_t              subFrameNum;

    if(MmwDemo_CLIGetSubframe(argc, argv, 3, &subFrameNum) < 0)
    {
        return -1;
    }

    /* Initialize configuration for clutter removal */
    memset ((void *)&cfg, 0, sizeof(MmwDemo_ClutterRemovalCfg));

    /* Populate configuration: */
    cfg.enabled          = (uint16_t) atoi (argv[2]);


    /* Save Configuration to use later */
    MmwDemo_mssCfgUpdate((void *)&cfg, offsetof(MmwDemo_CliCfg_t, clutterRemovalCfg),
        sizeof(MmwDemo_ClutterRemovalCfg), subFrameNum);

    /* Send configuration to DSS */
    memset((void *)&message, 0, sizeof(MmwDemo_message));

    message.type = MMWDEMO_MSS2DSS_CLUTTER_REMOVAL;
    message.subFrameNum = subFrameNum;
    memcpy((void *)&message.body.clutterRemovalCfg, (void *)&cfg, sizeof(MmwDemo_ClutterRemovalCfg));

    if (MmwDemo_mboxWrite(&message) == 0)
        return 0;
    else
        return -1;
}


/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for data logger set command
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
static int32_t MmwDemo_CLIADCBufCfg (int32_t argc, char* argv[])
{
    MmwDemo_ADCBufCfg   adcBufCfg;
    MmwDemo_message     message;
    int8_t              subFrameNum;

    if(MmwDemo_CLIGetSubframe(argc, argv, 6, &subFrameNum) < 0)
    {
        return -1;
    }

    /* Initialize the ADC Output configuration: */
    memset ((void *)&adcBufCfg, 0, sizeof(MmwDemo_ADCBufCfg));

    /* Populate configuration: */
    adcBufCfg.adcFmt          = (uint8_t) atoi (argv[2]);
    adcBufCfg.iqSwapSel       = (uint8_t) atoi (argv[3]);
    adcBufCfg.chInterleave    = (uint8_t) atoi (argv[4]);
    adcBufCfg.chirpThreshold  = (uint8_t) atoi (argv[5]);

    /* Save Configuration to use later */
    MmwDemo_mssCfgUpdate((void *)&adcBufCfg, offsetof(MmwDemo_CliCfg_t, adcBufCfg),
        sizeof(MmwDemo_ADCBufCfg), subFrameNum);

    /* Send configuration to DSS */
    memset((void *)&message, 0, sizeof(MmwDemo_message));
    message.type = MMWDEMO_MSS2DSS_ADCBUFCFG;
    message.subFrameNum = subFrameNum;
    memcpy((void *)&message.body.adcBufCfg, (void *)&adcBufCfg, sizeof(MmwDemo_ADCBufCfg));

    if (MmwDemo_mboxWrite(&message) == 0)
        return 0;
    else
        return -1;
}


/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for gui monitoring configuration
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
static int32_t ODDemo_CLIGuiMonSel(int32_t argc, char* argv[])
{
    MmwDemo_GuiMonSel   guiMonSel;
    MmwDemo_message     message;

    if (argc != 4)
    {
        CLI_write ("Error: Invalid usage of the guiMonitor command\n");
        return -1;
    }

    /* Initialize the guiMonSel configuration: */
    memset ((void *)&guiMonSel, 0, sizeof(MmwDemo_GuiMonSel));

    /* Populate configuration: */
    guiMonSel.decision            = atoi (argv[1]);
    guiMonSel.rangeAzimuthHeatMap = atoi (argv[2]);
    guiMonSel.vitalSign           = atoi (argv[3]);

    if ((guiMonSel.rangeAzimuthHeatMap != 0) &&
        (guiMonSel.rangeAzimuthHeatMap != 8) &&
        (guiMonSel.rangeAzimuthHeatMap != 16) &&
        (guiMonSel.rangeAzimuthHeatMap != 32))
    {
        CLI_write ("Invalid heatmap size\n");
        return -1;
    }


    MmwDemo_mssCfgUpdate((void *)&guiMonSel, offsetof(MmwDemo_CliCfg_t, guiMonSel),
        sizeof(MmwDemo_GuiMonSel), -1);

    /* Send configuration to DSS */
    memset((void *)&message, 0, sizeof(MmwDemo_message));

    message.type = MMWDEMO_MSS2DSS_GUIMON_CFG;
    message.subFrameNum = -1;
    memcpy((void *)&message.body.guiMonSel, (void *)&guiMonSel, sizeof(MmwDemo_GuiMonSel));

    if (MmwDemo_mboxWrite(&message) == 0)
        return 0;
    else
        return -1;
}


/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for Occupancy Demo zone definition
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
static int32_t ODDemo_CLIZoneDef(int32_t argc, char *argv[])
{
    uint16_t        idx;
    MmwDemo_message message;


    /* Populate configuration: */
    oddemo_num_zones = (uint16_t) atoi (argv[1]);

    /* Sanity Check: Minimum argument check */
    if (argc != (2 + (4 * oddemo_num_zones)))
    {
        CLI_write ("Error: Invalid usage of the zoneDef command\n");
        return -1;
    }

    if (oddemo_num_zones > ODDEMO_MAX_ZONES)
    {
        CLI_write ("Error: Exceeded max ODDemo zones\n");
        return -1;
    }

    /* Initialize configuration for DC range signature calibration */
    memset ((void *)&oddemo_zone, 0xff, sizeof(oddemo_zone));

    for (idx = 0; idx < oddemo_num_zones; idx ++)
    {
        oddemo_zone[idx].range_start    = (int16_t)atoi(argv[2 + (idx * 4)]);
        oddemo_zone[idx].range_length   = (int16_t)atoi(argv[3 + (idx * 4)]);
        oddemo_zone[idx].azimuth_start  = (int16_t)atoi(argv[4 + (idx * 4)]);
        oddemo_zone[idx].azimuth_width  = (int16_t)atoi(argv[5 + (idx * 4)]);
    }

    /* Sanity Check the zone definitions */
    for (idx = 0; idx < oddemo_num_zones; idx ++)
    {
        if (((oddemo_zone[idx].range_start + oddemo_zone[idx].range_length) > ODDEMO_MAX_RANGE) ||
            ((oddemo_zone[idx].azimuth_start + oddemo_zone[idx].azimuth_width) > ODDEMO_MAX_AZIMUTH))
        {
            CLI_write ("Error: ODDemo zone %d exceeds heat map size\n", idx);
            return -1;
        }
    }

    /* Send configuration to DSS */
    memset((void *)&message, 0, sizeof(MmwDemo_message));

    message.type = MMWDEMO_MSS2DSS_ODDEMO_ZONE_DEF;
    memcpy((void *)&message.body.zone, (void *)oddemo_zone, sizeof(oddemo_zone));

    if (MmwDemo_mboxWrite(&message) == 0)
        return 0;
    else
        return -1;
}


static int32_t ODDemo_CLICoeffMatrix(int32_t argc, char *argv[])
{
    uint16_t idx, pair, row;
    float   *ptr;
    MmwDemo_message message;

    /* Sanity Check: Minimum argument check */
    if (argc != (3 + ODDEMO_MATRIX_ROW_SIZE))
    {
        CLI_write ("Error: Invalid usage of the coeffMatrixRow command\n");
        return -1;
    }

    /* Get the pair number and sanity check */
    pair = (uint16_t) atoi (argv[1]);

    /* Sanity Check: Minimum argument check */
    if (pair >= ODDEMO_ZONE_PAIR)
    {
        CLI_write ("Error: Exceeded max ODDemo coeffMatrix pairs\n");
        return -1;
    }

    /* Get the row number and sanity check */
    row = (uint16_t) atoi (argv[2]);

    /* Sanity Check: Minimum argument check */
    if (row >= ODDEMO_MATRIX_NUM_ROWS)
    {
        CLI_write ("Error: Exceeded max ODDemo coeffMatrix rows\n");
        return -1;
    }

    ptr = &oddemo_coeffMatrix[pair][row * ODDEMO_MATRIX_ROW_SIZE];

    for (idx = 0; idx < ODDEMO_MATRIX_ROW_SIZE; idx ++)
        *ptr++ = (float) atof(argv[3+idx]);

    /* Send configuration to DSS */
    memset((void *)&message, 0, sizeof(MmwDemo_message));

    if (row == (ODDEMO_MATRIX_NUM_ROWS-1)) //send the matrix when it is completed
    {
        message.type = MMWDEMO_MSS2DSS_ODDEMO_COEFF_MATRIX;
        message.subFrameNum = pair;
        memcpy((void *)&message.body.coeffMatrix,
               (void *)&oddemo_coeffMatrix[pair][0], sizeof(float) * ODDEMO_MATRIX_SIZE);

        if (MmwDemo_mboxWrite(&message) == 0)
            return 0;
        else
            return -1;
    }
    else
        return 0;
}


static int32_t ODDemo_CLIMeanVector(int32_t argc, char *argv[])
{
    uint16_t idx, pair;
    MmwDemo_message message;

    /* Sanity Check: Minimum argument check */
    if (argc != (1 + ODDEMO_MATRIX_ROW_SIZE))
    {
        CLI_write ("Error: Invalid usage of the meanVector command\n");
        return -1;
    }

    /* Get the pair number and sanity check */
    pair = (uint16_t) atoi (argv[1]);

    /* Sanity Check: Minimum argument check */
    if (pair >= ODDEMO_ZONE_PAIR)
    {
        CLI_write ("Error: Exceeded max ODDemo meanVec pairs\n");
        return -1;
    }

    for (idx = 0; idx < (ODDEMO_MATRIX_ROW_SIZE-1); idx ++)
        oddemo_meanVec[pair][idx] = (float) atof(argv[2+idx]);

    /* Send configuration to DSS */
    memset((void *)&message, 0, sizeof(MmwDemo_message));

    message.type = MMWDEMO_MSS2DSS_ODDEMO_MEAN_VECTOR;
    message.subFrameNum = pair;
    memcpy((void *)&message.body.meanVec,
           (void *)&oddemo_meanVec[pair][0], sizeof(float) * (ODDEMO_MATRIX_ROW_SIZE-1));

    if (MmwDemo_mboxWrite(&message) == 0)
        return 0;
    else
        return -1;
}


static int32_t ODDemo_CLIStdVector(int32_t argc, char *argv[])
{
    uint16_t idx, pair;
    MmwDemo_message message;

    /* Sanity Check: Minimum argument check */
    if (argc != (1 + ODDEMO_MATRIX_ROW_SIZE))
    {
        CLI_write ("Error: Invalid usage of the stdVector command\n");
        return -1;
    }

    /* Get the pair number and sanity check */
    pair = (uint16_t) atoi (argv[1]);

    /* Sanity Check: Minimum argument check */
    if (pair >= ODDEMO_ZONE_PAIR)
    {
        CLI_write ("Error: Exceeded max ODDemo stdVector pairs\n");
        return -1;
    }

    for (idx = 0; idx < (ODDEMO_MATRIX_ROW_SIZE-1); idx ++)
        oddemo_stdVec[pair][idx] = (float) atof(argv[2+idx]);

    /* Send configuration to DSS */
    memset((void *)&message, 0, sizeof(MmwDemo_message));

    message.type = MMWDEMO_MSS2DSS_ODDEMO_STD_VECTOR;
    message.subFrameNum = pair;
    memcpy((void *)&message.body.stdVec,
           (void *)&oddemo_stdVec[pair][0], sizeof(float) * (ODDEMO_MATRIX_ROW_SIZE-1));

    if (MmwDemo_mboxWrite(&message) == 0)
        return 0;
    else
        return -1;
}

static int32_t ODDemo_CLIParms(int32_t argc, char *argv[])
{
    uint32_t windowLen;
    MmwDemo_message message;

    /* Sanity Check: Minimum argument check */
    if (argc != 3)
    {
        CLI_write ("Error: Invalid usage of the oddemoParms command\n");
        return -1;
    }

    windowLen = (uint32_t)atoi(argv[1]);

    if (windowLen > ODDEMO_MAX_FRAME_HIST)
        windowLen = ODDEMO_MAX_FRAME_HIST;

    /* Send configuration to DSS */
    memset((void *)&message, 0, sizeof(MmwDemo_message));

    message.type = MMWDEMO_MSS2DSS_ODDEMO_PARMS;
    message.body.parms.windowLen = windowLen;
    message.body.parms.gamma     = (float) atof(argv[2]);

    if (MmwDemo_mboxWrite(&message) == 0)
        return 0;
    else
        return -1;
}


static int32_t ODDemo_CLIRowNoise(int32_t argc, char *argv[])
{
    uint16_t idx, row, num;
    MmwDemo_message message;

    /* Get the first row number */
    row = (uint16_t) atoi (argv[1]);

    /* Get the number of rows */
    num = (uint16_t) atoi (argv[2]);

    /* Sanity Check: */
    if (num > 16)
    {
        CLI_write ("Error: Exceeded max row parms\n");
        return -1;
    }

    /* Sanity Check: */
    if ((row + num) > ODDEMO_MAX_RANGE)
    {
        CLI_write ("Error: Exceeded max ODDemo range rows\n");
        return -1;
    }

    /* Sanity Check: Minimum argument check */
    if (argc != (3 + num))
    {
        CLI_write ("Error: Invalid usage of the rowNoise command\n");
        return -1;
    }

    for (idx = 0; idx < num; idx ++)
        oddemo_row_noise[row+idx] = (float) atof(argv[3+idx]);

    //Send to DSS when we have all of them.
    if ((oddemo_row_noise[0] >= 0) && (oddemo_row_noise[ODDEMO_MAX_RANGE-1] >= 0))
    {
         /* Send configuration to DSS */
        memset((void *)&message, 0, sizeof(MmwDemo_message));

        message.type = MMWDEMO_MSS2DSS_ODDEMO_ROW_NOISE;
        memcpy((void *)&message.body.rowNoise,
               (void *)&oddemo_row_noise[0], sizeof(float) * ODDEMO_MAX_RANGE);

        if (MmwDemo_mboxWrite(&message) == 0)
            return 0;
        else
            return -1;
    }

    return 0;
}


/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for gui monitoring configuration
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
static int32_t VitalSignsDemo_CLIGuiMonSel (int32_t argc, char* argv[])
{
    VitalSignsDemo_GuiMonSel   guiMonSel;
    MmwDemo_message     message;

    int8_t              subFrameNum;
    subFrameNum = -1;
   /* if(MmwDemo_CLIGetSubframe(argc, argv, 8, &subFrameNum) < 0)
    {
        return -1;
    }
*/
    /* Initialize the guiMonSel configuration: */
    memset ((void *)&guiMonSel, 0, sizeof(VitalSignsDemo_GuiMonSel));

    /* Populate configuration: */
    guiMonSel.guiFlag_Param1         = atoi(argv[1]);
    guiMonSel.guiFlag_Param2         = atoi(argv[2]);
    guiMonSel.guiFlag_ClutterRemoval = atoi(argv[3]);
    guiMonSel.guiFlag_Reset          = atoi(argv[4]);
    guiMonSel.statsInfo              = atoi(argv[5]);

    MmwDemo_mssCfgUpdate((void *)&guiMonSel, offsetof(MmwDemo_CliCfg_t, vitalSigns_GuiMonSel),
        sizeof(guiMonSel), subFrameNum);

    /* Send configuration to DSS */
    memset((void *)&message, 0, sizeof(MmwDemo_message));

    message.type = MMWDEMO_VITALSIGNS_GUIMON_CFG;
    message.subFrameNum = subFrameNum;
    memcpy((void *)&message.body.vitalSigns_GuiMonSel, (void *)&guiMonSel, sizeof(guiMonSel));

    if (MmwDemo_mboxWrite(&message) == 0)
        return 0;
    else
        return -1;
}

/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for Vital Signs parameters configuration
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
static int32_t VitalSignsDemo_CLIvitalSignsParamsCfg (int32_t argc, char* argv[])
{

    VitalSignsDemo_ParamsCfg   vitalSignsParamsCfg;
    MmwDemo_message     message;

    /* Sanity Check: Minimum argument check */
    if (argc < 8)
    {
        CLI_write ("Error: Invalid usage of the CLI command\n");
        return -1;
    }

    /* Initialize the ADC Output configuration: */
    memset ((void *)&vitalSignsParamsCfg, 0, sizeof(VitalSignsDemo_ParamsCfg));

    /* Populate configuration: */
    vitalSignsParamsCfg.startRange_m         = (float) atof (argv[1]);
    vitalSignsParamsCfg.endRange_m           = (float) atof (argv[2]);
    vitalSignsParamsCfg.winLen_breathing     = (uint16_t) atoi (argv[3]);
    vitalSignsParamsCfg.winLen_heartRate     = (uint16_t) atoi (argv[4]);
    vitalSignsParamsCfg.rxAntennaProcess     = (float)   atof (argv[5]);
    vitalSignsParamsCfg.alpha_breathingWfm   = (float)   atof (argv[6]);
    vitalSignsParamsCfg.alpha_heartWfm       = (float)   atof (argv[7]);
    vitalSignsParamsCfg.scale_breathingWfm   = (float)   atof (argv[8]);
    vitalSignsParamsCfg.scale_heartWfm       = (float)   atof (argv[9]);

    /* Save Configuration to use later */
    MmwDemo_mssCfgUpdate((void *)&vitalSignsParamsCfg, offsetof(MmwDemo_CliCfg_t, vitalSignsParamsCfg),
        sizeof(vitalSignsParamsCfg), -1);

    /* Send configuration to DSS */
    memset((void *)&message, 0, sizeof(MmwDemo_message));
    message.type = MMWDEMO_MSS2DSS_VITALSIGNS_MEASUREMENT_PARAMS;
    memcpy((void *)&message.body.vitalSignsParamsCfg, (void *)&vitalSignsParamsCfg, sizeof(vitalSignsParamsCfg));

    if (MmwDemo_mboxWrite(&message) == 0)
        return 0;
    else
        return -1;

}

/**
 *  @b Description
 *  @n
 *      This is the CLI Handler for Vital Signs Motion Detection
 *
 *  @param[in] argc
 *      Number of arguments
 *  @param[in] argv
 *      Arguments
 *
 *  @retval
 *      Success -   0
 *  @retval
 *      Error   -   <0
 */
static int32_t VitalSignsDemo_CLIMotionDetection (int32_t argc, char* argv[])
{

    //CLI_write ("Error: Motion Detection\n");

    VitalSignsDemo_MotionDetection   motionDetectionParamsCfg;
    MmwDemo_message     message;

    /* Sanity Check: Minimum argument check */
    if (argc < 4)
    {
        CLI_write ("Error: Invalid usage of the CLI command\n");
        return -1;
    }

    /* Initialize the ADC Output configuration: */
    memset ((void *)&motionDetectionParamsCfg, 0, sizeof(motionDetectionParamsCfg));

    /* Populate configuration: */

    motionDetectionParamsCfg.enabled     = (uint16_t) atoi (argv[1]);
    motionDetectionParamsCfg.blockSize   = (uint16_t) atoi (argv[2]);
    motionDetectionParamsCfg.threshold    = (float) atof (argv[3]);
    motionDetectionParamsCfg.gainControl  = (uint16_t) atof (argv[4]);

    /* Save Configuration to use later */
    MmwDemo_mssCfgUpdate((void *)&motionDetectionParamsCfg, offsetof(MmwDemo_CliCfg_t, motionDetectionParamsCfg),
        sizeof(motionDetectionParamsCfg), -1);

    /* Send configuration to DSS */
    memset((void *)&message, 0, sizeof(MmwDemo_message));
    message.type = MMWDEMO_MSS2DSS_VITALSIGNS_MOTION_DETECTION;
    memcpy((void *)&message.body.motionDetectionParamsCfg, (void *)&motionDetectionParamsCfg, sizeof(motionDetectionParamsCfg));


    if (MmwDemo_mboxWrite(&message) == 0)
        return 0;
    else
        return -1;

}

/**
 *  @b Description
 *  @n
 *      This is the CLI Execution Task
 *
 *  @retval
 *      Not Applicable.
 */
void MmwDemo_CLIInit (void)
{
    CLI_Cfg     cliCfg;
    char        demoBanner[256];
    uint32_t    cnt;

    /* Create Demo Banner to be printed out by CLI */
    sprintf(&demoBanner[0],
                       "*********************************************************\n" \
                       "IWR6843ISK Multiple Vital Signs Tracking %02d.%02d.%02d.%02d\n"  \
                       "*********************************************************\n",
                        MMWAVE_SDK_VERSION_MAJOR,
                        MMWAVE_SDK_VERSION_MINOR,
                        MMWAVE_SDK_VERSION_BUGFIX,
                        MMWAVE_SDK_VERSION_BUILD
            );

    /* Initialize the CLI configuration: */
    memset ((void *)&cliCfg, 0, sizeof(CLI_Cfg));

    cliCfg.socHandle = gMmwMssMCB.socHandle;

    /* Populate the CLI configuration: */
    cliCfg.cliPrompt                    = "MVST:/>";
    cliCfg.cliBanner                    = demoBanner;
    cliCfg.cliUartHandle                = gMmwMssMCB.commandUartHandle;
    cliCfg.taskPriority                 = 3;
    cliCfg.mmWaveHandle                 = gMmwMssMCB.ctrlHandle;
    cliCfg.enableMMWaveExtension        = 1U;
    cliCfg.usePolledMode                = true;
    cnt=0;
    cliCfg.tableEntry[cnt].cmd            = "sensorStart";
    cliCfg.tableEntry[cnt].helpString     = "[doReconfig(optional, default:enabled)]";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = MmwDemo_CLISensorStart;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "sensorStop";
    cliCfg.tableEntry[cnt].helpString     = "No arguments";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = MmwDemo_CLISensorStop;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "guiMonitor";
    cliCfg.tableEntry[cnt].helpString     = "<subFrameIdx> <featureVector> <decision> <rangeAzimuthHeatMap> <VitalSign>";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = ODDemo_CLIGuiMonSel;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "dvsGuiMonitor";
    cliCfg.tableEntry[cnt].helpString     = "<guiFlag_Param1> <guiFlag_Param2> <guiFlag_ClutterRemoval> <guiFlag_Reset> <statsInfo>";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = VitalSignsDemo_CLIGuiMonSel;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "calibDcRangeSig";
    cliCfg.tableEntry[cnt].helpString     = "<subFrameIdx> <enabled> <negativeBinIdx> <positiveBinIdx> <numAvgFrames>";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = MmwDemo_CLICalibDcRangeSig;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "clutterRemoval";
    cliCfg.tableEntry[cnt].helpString     = "<enabled>";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = MmwDemo_CLIClutterRemoval;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "adcbufCfg";
    cliCfg.tableEntry[cnt].helpString     = "<subFrameIdx> <adcOutputFmt> <SampleSwap> <ChanInterleave> <ChirpThreshold>";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = MmwDemo_CLIADCBufCfg;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "zoneDef";
    cliCfg.tableEntry[cnt].helpString     = "<zone count> [<range start> <range len> <az start> <az len>]";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = ODDemo_CLIZoneDef;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "coeffMatrixRow";
    cliCfg.tableEntry[cnt].helpString     = "<pair#> <row#> <6 float values>";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = ODDemo_CLICoeffMatrix;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "meanVector";
    cliCfg.tableEntry[cnt].helpString     = "<pair#> <5 float values>";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = ODDemo_CLIMeanVector;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "stdVector";
    cliCfg.tableEntry[cnt].helpString     = "<pair#> <5 float values>";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = ODDemo_CLIStdVector;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "oddemoParms";
    cliCfg.tableEntry[cnt].helpString     = "<window len> <gamma>";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = ODDemo_CLIParms;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "rowNoise";
    cliCfg.tableEntry[cnt].helpString     = "<1st row> <N=num rows> <N float values>";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = ODDemo_CLIRowNoise;
    cnt++;

    oddemo_row_noise[0] = -1.0;
    oddemo_row_noise[ODDEMO_MAX_RANGE-1] = -1.0;

    cliCfg.tableEntry[cnt].cmd            = "vitalSignsCfg";
    cliCfg.tableEntry[cnt].helpString     = "<rangeStart_meters> <rangeEnd_meters> <winLenBreath> <winLenHeart> <RX_recieverProcess> <alpha_breathingWfm> <alpha_heartWfm> <scale_breathingWfm> <scale_heartWfm>";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = VitalSignsDemo_CLIvitalSignsParamsCfg;
    cnt++;

    cliCfg.tableEntry[cnt].cmd            = "motionDetection";
    cliCfg.tableEntry[cnt].helpString     = "<enabled> <threshold> <blockSize> <gainControl>";
    cliCfg.tableEntry[cnt].cmdHandlerFxn  = VitalSignsDemo_CLIMotionDetection;
    cnt++;

    /* Open the CLI: */
    if (CLI_open (&cliCfg) < 0)
    {
        System_printf ("Error: Unable to open the CLI\n");
        return;
    }
    System_printf ("Debug: CLI is operational\n");
    return;
}
