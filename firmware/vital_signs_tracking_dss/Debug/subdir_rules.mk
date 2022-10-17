################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Each subdirectory must supply rules for building sources it contributes
%.oe674: ../%.c $(GEN_OPTS) | $(GEN_FILES) $(GEN_MISC_FILES)
	@echo 'Building file: "$<"'
	@echo 'Invoking: C6000 Compiler'
	"/home/dirk/ti/ti-cgt-c6000_8.3.3/bin/cl6x" -mv6740 --abi=eabi -O3 --opt_for_speed=3 --include_path="/home/dirk/Documents/TI/workspace3/vital_signs_tracking_dss" --include_path="/home/dirk/ti/mmwave_sdk_03_05_00_04/packages" --include_path="/home/dirk/ti/mathlib_c674x_3_1_2_1/packages" --include_path="/home/dirk/ti/dsplib_c64Px_3_4_0_0/packages/ti/dsplib/src/DSP_fft16x16/c64P" --include_path="/home/dirk/ti/dsplib_c64Px_3_4_0_0/packages/ti/dsplib/src/DSP_fft32x32/c64P" --include_path="/home/dirk/Documents/TI/workspace3/vital_signs_tracking_dss/dss" --include_path="/home/dirk/Documents/TI/workspace3/vital_signs_tracking_dss/common" --include_path="/home/dirk/ti/mathlib_c674x_3_1_2_1/packages" --include_path="/home/dirk/ti/dsplib_c64Px_3_4_0_0/packages/" --include_path="/home/dirk/ti/dsplib_c64Px_3_4_0_0/packages/ti/dsplib/src/DSP_fft16x16/c64P" --include_path="/home/dirk/ti/dsplib_c64Px_3_4_0_0/packages/ti/dsplib/src/DSP_fft32x32/c64P" --include_path="/home/dirk/ti/ti-cgt-c6000_8.3.3/include" --define=SOC_XWR68XX --define=SUBSYS_DSS --define=MMWAVE_L3RAM_SIZE=0xC0000 --define=MMWAVE_L3RAM_NUM_BANK=6 --define=MMWAVE_SHMEM_BANK_SIZE=0x20000 --define=DOWNLOAD_FROM_CCS --define=DebugP_ASSERT_ENABLED -g --gcc --diag_warning=225 --diag_wrap=off --display_error_number --gen_func_subsections=on --obj_extension=.oe674 --preproc_with_compile --preproc_dependency="$(basename $(<F)).d_raw" $(GEN_OPTS__FLAG) "$(shell echo $<)"
	@echo 'Finished building: "$<"'
	@echo ' '

build-955202767:
	@$(MAKE) --no-print-directory -Onone -f subdir_rules.mk build-955202767-inproc

build-955202767-inproc: ../dss_mmw.cfg
	@echo 'Building file: "$<"'
	@echo 'Invoking: XDCtools'
	"/home/dirk/ti/xdctools_3_50_08_24_core/xs" --xdcpath="/home/dirk/ti/bios_6_73_01_01/packages;" xdc.tools.configuro -o configPkg -t ti.targets.elf.C674 -p ti.platforms.c6x:IWR68XX:false:600 -r debug -c "/home/dirk/ti/ti-cgt-c6000_8.3.3" "$<"
	@echo 'Finished building: "$<"'
	@echo ' '

configPkg/linker.cmd: build-955202767 ../dss_mmw.cfg
configPkg/compiler.opt: build-955202767
configPkg/: build-955202767


