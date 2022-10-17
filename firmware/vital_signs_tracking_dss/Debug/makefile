################################################################################
# Automatically-generated file. Do not edit!
################################################################################

CG_TOOL_ROOT := /home/dirk/ti/ti-cgt-c6000_8.3.3

GEN_OPTS__FLAG := --cmd_file="configPkg/compiler.opt" 
GEN_CMDS__FLAG := -l"configPkg/linker.cmd" 

ORDERED_OBJS += \
"./dss_config_edma_util.oe674" \
"./dss_data_path.oe674" \
"./dss_main.oe674" \
"./dss_vitalSignsDemo_utilsFunc.oe674" \
"./oddemo_decision.oe674" \
"./oddemo_feature.oe674" \
"./oddemo_heatmap.oe674" \
"./xwr16xx_cache.oe674" \
"../c674x_68xx_linker.cmd" \
"../dss_68xx_linker.cmd" \
$(GEN_CMDS__FLAG) \
-ldsplib.lib \
-lmathlib.ae674 \
-llibsoc_xwr68xx.ae674 \
-llibosal_xwr68xx.ae674 \
-llibcrc_xwr68xx.ae674 \
-llibmailbox_xwr68xx.ae674 \
-llibmmwavelink_xwr68xx.ae674 \
-llibmmwave_xwr68xx.ae674 \
-llibedma_xwr68xx.ae674 \
-llibadcbuf_xwr68xx.ae674 \
-llibmmwavealg_xwr68xx.ae674 \
-ldsplib.ae64P \
-llibcbuff_xwr68xx.ae674 \
-llibhsiheader_xwr68xx.ae674 \
-lrts6740_elf.lib \
-llibc.a \

-include ../makefile.init

RM := rm -rf
RMDIR := rm -rf

# All of the sources participating in the build are defined here
-include sources.mk
-include subdir_vars.mk
-include subdir_rules.mk
-include objects.mk

ifneq ($(MAKECMDGOALS),clean)
ifneq ($(strip $(C55_DEPS)),)
-include $(C55_DEPS)
endif
ifneq ($(strip $(C_UPPER_DEPS)),)
-include $(C_UPPER_DEPS)
endif
ifneq ($(strip $(S67_DEPS)),)
-include $(S67_DEPS)
endif
ifneq ($(strip $(S62_DEPS)),)
-include $(S62_DEPS)
endif
ifneq ($(strip $(S_DEPS)),)
-include $(S_DEPS)
endif
ifneq ($(strip $(OPT_DEPS)),)
-include $(OPT_DEPS)
endif
ifneq ($(strip $(C??_DEPS)),)
-include $(C??_DEPS)
endif
ifneq ($(strip $(ASM_UPPER_DEPS)),)
-include $(ASM_UPPER_DEPS)
endif
ifneq ($(strip $(S??_DEPS)),)
-include $(S??_DEPS)
endif
ifneq ($(strip $(C64_DEPS)),)
-include $(C64_DEPS)
endif
ifneq ($(strip $(CXX_DEPS)),)
-include $(CXX_DEPS)
endif
ifneq ($(strip $(S64_DEPS)),)
-include $(S64_DEPS)
endif
ifneq ($(strip $(INO_DEPS)),)
-include $(INO_DEPS)
endif
ifneq ($(strip $(CLA_DEPS)),)
-include $(CLA_DEPS)
endif
ifneq ($(strip $(S55_DEPS)),)
-include $(S55_DEPS)
endif
ifneq ($(strip $(SV7A_DEPS)),)
-include $(SV7A_DEPS)
endif
ifneq ($(strip $(C62_DEPS)),)
-include $(C62_DEPS)
endif
ifneq ($(strip $(C67_DEPS)),)
-include $(C67_DEPS)
endif
ifneq ($(strip $(PDE_DEPS)),)
-include $(PDE_DEPS)
endif
ifneq ($(strip $(K_DEPS)),)
-include $(K_DEPS)
endif
ifneq ($(strip $(C_DEPS)),)
-include $(C_DEPS)
endif
ifneq ($(strip $(CC_DEPS)),)
-include $(CC_DEPS)
endif
ifneq ($(strip $(C++_DEPS)),)
-include $(C++_DEPS)
endif
ifneq ($(strip $(C43_DEPS)),)
-include $(C43_DEPS)
endif
ifneq ($(strip $(S43_DEPS)),)
-include $(S43_DEPS)
endif
ifneq ($(strip $(ASM_DEPS)),)
-include $(ASM_DEPS)
endif
ifneq ($(strip $(S_UPPER_DEPS)),)
-include $(S_UPPER_DEPS)
endif
ifneq ($(strip $(CPP_DEPS)),)
-include $(CPP_DEPS)
endif
ifneq ($(strip $(SA_DEPS)),)
-include $(SA_DEPS)
endif
endif

-include ../makefile.defs

# Add inputs and outputs from these tool invocations to the build variables 
EXE_OUTPUTS += \
vital_signs_tracking_dss.xe674 \

EXE_OUTPUTS__QUOTED += \
"vital_signs_tracking_dss.xe674" \

BIN_OUTPUTS += \
vital_signs_tracking_dss.hex \

BIN_OUTPUTS__QUOTED += \
"vital_signs_tracking_dss.hex" \


# All Target
all:
	@$(MAKE) --no-print-directory pre-build
	@$(MAKE) --no-print-directory main-build

# Main-build Target
main-build: $(OBJS) $(CMD_SRCS) $(GEN_CMDS)
	@$(MAKE) --no-print-directory -Onone "vital_signs_tracking_dss.xe674"

# Tool invocations
vital_signs_tracking_dss.xe674: $(OBJS) $(CMD_SRCS) $(GEN_CMDS)
	@echo 'Building target: "$@"'
	@echo 'Invoking: C6000 Linker'
	"/home/dirk/ti/ti-cgt-c6000_8.3.3/bin/cl6x" -mv6740 --abi=eabi -O3 --opt_for_speed=3 --define=SOC_XWR68XX --define=SUBSYS_DSS --define=MMWAVE_L3RAM_SIZE=0xC0000 --define=MMWAVE_L3RAM_NUM_BANK=6 --define=MMWAVE_SHMEM_BANK_SIZE=0x20000 --define=DOWNLOAD_FROM_CCS --define=DebugP_ASSERT_ENABLED -g --gcc --diag_warning=225 --diag_wrap=off --display_error_number --gen_func_subsections=on --obj_extension=.oe674 -z -m"vital_signs_tracking_dss.map" --heap_size=0x800 --stack_size=0x800 -i"/home/dirk/ti/ti-cgt-c6000_8.3.3/lib" -i"/home/dirk/Documents/TI/workspace3/vital_signs_tracking_dss" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/control/mmwave/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/control/mmwavelink/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/alg/mmwavelib/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/crc/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/mailbox/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/soc/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/osal/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/edma/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/adcbuf/lib" -i"/home/dirk/ti/dsplib_c64Px_3_4_0_0/packages/ti/dsplib/lib" -i"/home/dirk/ti/mathlib_c674x_3_1_2_1/packages/ti/mathlib/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/cbuff/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/utils/hsiheader/lib" -i"/home/dirk/ti/dsplib_c64Px_3_4_0_0/packages/ti/dsplib/lib" -i"/home/dirk/ti/mathlib_c674x_3_1_2_1/packages/ti/mathlib/lib" -i"/home/dirk/ti/ti-cgt-c6000_8.3.3/include" --reread_libs --disable_auto_rts --define=MMWAVE_L3RAM_NUM_BANK=6 --define=MMWAVE_SHMEM_TCMA_NUM_BANK=0 --define=MMWAVE_SHMEM_TCMB_NUM_BANK=0 --define=MMWAVE_SHMEM_BANK_SIZE=0x20000 --diag_warning=225 --diag_wrap=off --display_error_number --warn_sections --xml_link_info="xwr68xx_vital_signs_tracking_dss_linkInfo.xml" --ram_model --unused_section_elimination=on -o "vital_signs_tracking_dss.xe674" $(ORDERED_OBJS)
	@echo 'Finished building target: "$@"'
	@echo ' '
	@$(MAKE) --no-print-directory post-build

vital_signs_tracking_dss.hex: $(EXE_OUTPUTS)
	@echo 'Building secondary target: "$@"'
	@echo 'Invoking: C6000 Hex Utility'
	"/home/dirk/ti/ti-cgt-c6000_8.3.3/bin/hex6x" --diag_wrap=off -o "vital_signs_tracking_dss.hex" $(EXE_OUTPUTS__QUOTED)
	@echo 'Finished building secondary target: "$@"'
	@echo ' '
	@$(MAKE) --no-print-directory post-build

# Other Targets
clean:
	-$(RM) $(BIN_OUTPUTS__QUOTED)$(GEN_FILES__QUOTED)$(EXE_OUTPUTS__QUOTED)
	-$(RM) "dss_config_edma_util.oe674" "dss_data_path.oe674" "dss_main.oe674" "dss_vitalSignsDemo_utilsFunc.oe674" "oddemo_decision.oe674" "oddemo_feature.oe674" "oddemo_heatmap.oe674" "xwr16xx_cache.oe674" 
	-$(RM) "dss_config_edma_util.d" "dss_data_path.d" "dss_main.d" "dss_vitalSignsDemo_utilsFunc.d" "oddemo_decision.d" "oddemo_feature.d" "oddemo_heatmap.d" "xwr16xx_cache.d" 
	-$(RMDIR) $(GEN_MISC_DIRS__QUOTED)
	-@echo 'Finished clean'
	-@echo ' '

pre-build:
	-rm -f /home/dirk/Documents/TI/workspace3/vital_signs_tracking_dss/vital_signs_tracking_dss.bin
	-@echo ' '

post-build:
	-/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/scripts/ImageCreator/out2rprc/out2rprc.exe vital_signs_tracking_dss.xe674 vital_signs_tracking_dss.bin
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets

