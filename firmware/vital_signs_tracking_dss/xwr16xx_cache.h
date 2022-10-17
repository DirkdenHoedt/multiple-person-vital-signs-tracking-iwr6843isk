/**
 * xwr16xx_cache.h
 *
 * Header file for cache related functions for the C674x DSP
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
 
#ifndef _IWR16XX_CACHE_
#define  _IWR16XX_CACHE_
#ifdef _TMS320C6X
#include "c6x.h"
#endif



#define MAR (0x01848000)
#define Cache_PC (1)
#define Cache_WTE (2)
#define Cache_PCX (4)
#define Cache_PFX (8)

#define CACHE_L1_32KCACHE 0x7
#define CACHE_L1_16KCACHE 0x3
#define CACHE_L1_8KCACHE 0x2
#define CACHE_L1_4KCACHE 0x1
#define CACHE_0KCACHE    0x0



extern void cache_setMar(unsigned int * baseAddr, unsigned int byteSize, unsigned int value);

extern void cache_setL1PSize(int cacheConf);
extern void cache_setL1DSize(int cacheConf);
extern void cache_setL2Size(int cacheConf);
extern void     cache_wbInvAllL2Wait();

#endif //_IWR16XX_CACHE_
