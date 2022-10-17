################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CFG_SRCS += \
../mss_mmw.cfg 

CMD_SRCS += \
../mss_mmw_linker.cmd \
../r4f_linker.cmd 

C_SRCS += \
../cli.c \
../mss_main.c 

GEN_CMDS += \
./configPkg/linker.cmd 

GEN_FILES += \
./configPkg/linker.cmd \
./configPkg/compiler.opt 

GEN_MISC_DIRS += \
./configPkg/ 

C_DEPS += \
./cli.d \
./mss_main.d 

GEN_OPTS += \
./configPkg/compiler.opt 

OBJS += \
./cli.obj \
./mss_main.obj 

GEN_MISC_DIRS__QUOTED += \
"configPkg/" 

OBJS__QUOTED += \
"cli.obj" \
"mss_main.obj" 

C_DEPS__QUOTED += \
"cli.d" \
"mss_main.d" 

GEN_FILES__QUOTED += \
"configPkg/linker.cmd" \
"configPkg/compiler.opt" 

C_SRCS__QUOTED += \
"../cli.c" \
"../mss_main.c" 


