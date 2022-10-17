################################################################################
# Automatically-generated file. Do not edit!
################################################################################

CG_TOOL_ROOT := /home/dirk/ti/ti-cgt-arm_16.9.6.LTS

GEN_OPTS__FLAG := --cmd_file="configPkg/compiler.opt" 
GEN_CMDS__FLAG := -l"configPkg/linker.cmd" 

ORDERED_OBJS += \
"./cli.obj" \
"./mss_main.obj" \
"../mss_mmw_linker.cmd" \
"../r4f_linker.cmd" \
$(GEN_CMDS__FLAG) \
-llibosal_xwr68xx.aer4f \
-llibesm_xwr68xx.aer4f \
-llibgpio_xwr68xx.aer4f \
-llibsoc_xwr68xx.aer4f \
-llibmailbox_xwr68xx.aer4f \
-llibpinmux_xwr68xx.aer4f \
-llibcrc_xwr68xx.aer4f \
-llibdma_xwr68xx.aer4f \
-llibuart_xwr68xx.aer4f \
-llibmmwavelink_xwr68xx.aer4f \
-llibmmwave_xwr68xx.aer4f \
-llibcli_xwr68xx.aer4f \
-llibadcbuf_xwr68xx.aer4f \
-llibedma_xwr68xx.aer4f \
-llibmathutils.aer4f \
-llibqspi_xwr68xx.aer4f \
-llibqspiflash_xwr68xx.aer4f \
-lrtsv7R4_T_le_v3D16_eabi.lib \
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
vital_signs_tracking_mss.xer4f \

EXE_OUTPUTS__QUOTED += \
"vital_signs_tracking_mss.xer4f" \

BIN_OUTPUTS += \
vital_signs_tracking_mss.hex \

BIN_OUTPUTS__QUOTED += \
"vital_signs_tracking_mss.hex" \


# All Target
all:
	@$(MAKE) --no-print-directory pre-build
	@$(MAKE) --no-print-directory main-build

# Main-build Target
main-build: $(OBJS) $(CMD_SRCS) $(GEN_CMDS)
	@$(MAKE) --no-print-directory -Onone "vital_signs_tracking_mss.xer4f"

# Tool invocations
vital_signs_tracking_mss.xer4f: $(OBJS) $(CMD_SRCS) $(GEN_CMDS)
	@echo 'Building target: "$@"'
	@echo 'Invoking: Arm Linker'
	"/home/dirk/ti/ti-cgt-arm_16.9.6.LTS/bin/armcl" -mv7R4 --code_state=16 --float_support=VFPv3D16 -me -O3 --define=_LITTLE_ENDIAN --define=SOC_XWR68XX --define=SUBSYS_MSS --define=DOWNLOAD_FROM_CCS --define=MMWAVE_L3RAM_SIZE=0xC0000 --define=DebugP_ASSERT_ENABLED -g --c99 --diag_warning=225 --diag_wrap=off --display_error_number --gen_func_subsections=on --enum_type=int --abi=eabi -z -m"vital_signs_tracking_mss.map" --heap_size=0x800 --stack_size=0x800 -i"/home/dirk/ti/ti-cgt-arm_16.9.6.LTS/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/control/mmwave/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/control/mmwavelink/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/crc/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/esm/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/gpio/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/osal/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/pinmux/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/mailbox/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/dma/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/soc/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/uart/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/utils/cli/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/adcbuf/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/edma/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/utils/mathutils/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/qspi/lib" -i"/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/ti/drivers/qspiflash/lib" -i"/home/dirk/ti/ti-cgt-arm_16.9.6.LTS/include" --reread_libs --disable_auto_rts --define=MMWAVE_L3RAM_NUM_BANK=6 --define=MMWAVE_SHMEM_TCMA_NUM_BANK=0 --define=MMWAVE_SHMEM_TCMB_NUM_BANK=0 --define=MMWAVE_SHMEM_BANK_SIZE=0x20000 --define=MMWAVE_L3RAM_SIZE=0xC0000 --diag_warning=225 --diag_wrap=off --display_error_number --warn_sections --xml_link_info="vital_signs_tracking_mss.xml" --rom_model --unused_section_elimination=on -o "vital_signs_tracking_mss.xer4f" $(ORDERED_OBJS)
	@echo 'Finished building target: "$@"'
	@echo ' '
	@$(MAKE) --no-print-directory post-build

vital_signs_tracking_mss.hex: $(EXE_OUTPUTS)
	@echo 'Building secondary target: "$@"'
	@echo 'Invoking: Arm Hex Utility'
	"/home/dirk/ti/ti-cgt-arm_16.9.6.LTS/bin/armhex" --diag_wrap=off -o "vital_signs_tracking_mss.hex" $(EXE_OUTPUTS__QUOTED)
	@echo 'Finished building secondary target: "$@"'
	@echo ' '
	@$(MAKE) --no-print-directory post-build

# Other Targets
clean:
	-$(RM) $(BIN_OUTPUTS__QUOTED)$(GEN_FILES__QUOTED)$(EXE_OUTPUTS__QUOTED)
	-$(RM) "cli.obj" "mss_main.obj" 
	-$(RM) "cli.d" "mss_main.d" 
	-$(RMDIR) $(GEN_MISC_DIRS__QUOTED)
	-@echo 'Finished clean'
	-@echo ' '

pre-build:
	-rm -f /home/dirk/Documents/TI/workspace3/vital_signs_tracking_mss/vital_signs_tracking_mss.bin
	-@echo ' '

post-build:
	-/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/scripts/ImageCreator/out2rprc/out2rprc.exe vital_signs_tracking_mss.xer4f vital_signs_tracking_mss.bin
	-/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/scripts/ImageCreator/multicore_image_generator/MulticoreImageGen LE 37 0x00000006 vital_signs_tracking.bin 0x35510000 vital_signs_tracking_mss.bin 0xb5510000 /home/dirk/ti/mmwave_sdk_03_05_00_04/firmware/radarss/xwr6xxx_radarss_rprc.bin 0xd5510000 /home/dirk/Documents/TI/workspace3/vital_signs_tracking_dss/Debug/vital_signs_tracking_dss.bin
	-/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/scripts/ImageCreator/crc_multicore_image/crc_multicore_image vital_signs_tracking.bin vital_signs_tracking.bin.tmp
	-/home/dirk/ti/mmwave_sdk_03_05_00_04/packages/scripts/ImageCreator/append_bin_crc/gen_bincrc32 vital_signs_tracking.bin
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:

-include ../makefile.targets

