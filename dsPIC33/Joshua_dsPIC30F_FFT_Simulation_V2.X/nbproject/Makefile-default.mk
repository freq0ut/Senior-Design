#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Include project Makefile
ifeq "${IGNORE_LOCAL}" "TRUE"
# do not include local makefile. User is passing all local related variables already
else
include Makefile
# Include makefile containing local settings
ifeq "$(wildcard nbproject/Makefile-local-default.mk)" "nbproject/Makefile-local-default.mk"
include nbproject/Makefile-local-default.mk
endif
endif

# Environment
MKDIR=mkdir -p
RM=rm -f 
MV=mv 
CP=cp 

# Macros
CND_CONF=default
ifeq ($(TYPE_IMAGE), DEBUG_RUN)
IMAGE_TYPE=debug
OUTPUT_SUFFIX=elf
DEBUGGABLE_SUFFIX=elf
FINAL_IMAGE=dist/${CND_CONF}/${IMAGE_TYPE}/Joshua_dsPIC30F_FFT_Simulation_V2.X.${IMAGE_TYPE}.${OUTPUT_SUFFIX}
else
IMAGE_TYPE=production
OUTPUT_SUFFIX=hex
DEBUGGABLE_SUFFIX=elf
FINAL_IMAGE=dist/${CND_CONF}/${IMAGE_TYPE}/Joshua_dsPIC30F_FFT_Simulation_V2.X.${IMAGE_TYPE}.${OUTPUT_SUFFIX}
endif

# Object Directory
OBJECTDIR=build/${CND_CONF}/${IMAGE_TYPE}

# Distribution Directory
DISTDIR=dist/${CND_CONF}/${IMAGE_TYPE}

# Source Files Quoted if spaced
SOURCEFILES_QUOTED_IF_SPACED=src/main_FFTExample.c src/twiddleFactors.c /Users/betio32/Documents/myGitHub/Senior-Design/dsPIC33/Joshua_dsPIC30F_FFT_Simulation_V2.X/src/chan1_input_cosine40kHz.c

# Object Files Quoted if spaced
OBJECTFILES_QUOTED_IF_SPACED=${OBJECTDIR}/src/main_FFTExample.o ${OBJECTDIR}/src/twiddleFactors.o ${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o
POSSIBLE_DEPFILES=${OBJECTDIR}/src/main_FFTExample.o.d ${OBJECTDIR}/src/twiddleFactors.o.d ${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o.d

# Object Files
OBJECTFILES=${OBJECTDIR}/src/main_FFTExample.o ${OBJECTDIR}/src/twiddleFactors.o ${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o

# Source Files
SOURCEFILES=src/main_FFTExample.c src/twiddleFactors.c /Users/betio32/Documents/myGitHub/Senior-Design/dsPIC33/Joshua_dsPIC30F_FFT_Simulation_V2.X/src/chan1_input_cosine40kHz.c


CFLAGS=
ASFLAGS=
LDLIBSOPTIONS=

############# Tool locations ##########################################
# If you copy a project from one host to another, the path where the  #
# compiler is installed may be different.                             #
# If you open this project with MPLAB X in the new host, this         #
# makefile will be regenerated and the paths will be corrected.       #
#######################################################################
# fixDeps replaces a bunch of sed/cat/printf statements that slow down the build
FIXDEPS=fixDeps

.build-conf:  ${BUILD_SUBPROJECTS}
ifneq ($(INFORMATION_MESSAGE), )
	@echo $(INFORMATION_MESSAGE)
endif
	${MAKE}  -f nbproject/Makefile-default.mk dist/${CND_CONF}/${IMAGE_TYPE}/Joshua_dsPIC30F_FFT_Simulation_V2.X.${IMAGE_TYPE}.${OUTPUT_SUFFIX}

MP_PROCESSOR_OPTION=30F6012
MP_LINKER_FILE_OPTION=,--script=p30F6012.gld
# ------------------------------------------------------------------------------------
# Rules for buildStep: compile
ifeq ($(TYPE_IMAGE), DEBUG_RUN)
${OBJECTDIR}/src/main_FFTExample.o: src/main_FFTExample.c  nbproject/Makefile-${CND_CONF}.mk
	@${MKDIR} "${OBJECTDIR}/src" 
	@${RM} ${OBJECTDIR}/src/main_FFTExample.o.d 
	@${RM} ${OBJECTDIR}/src/main_FFTExample.o 
	${MP_CC} $(MP_EXTRA_CC_PRE)  src/main_FFTExample.c  -o ${OBJECTDIR}/src/main_FFTExample.o  -c -mcpu=$(MP_PROCESSOR_OPTION)  -MMD -MF "${OBJECTDIR}/src/main_FFTExample.o.d"      -g -D__DEBUG -D__MPLAB_DEBUGGER_SIMULATOR=1    -omf=elf -O0 -msmart-io=1 -Wall -msfr-warn=off
	@${FIXDEPS} "${OBJECTDIR}/src/main_FFTExample.o.d" $(SILENT)  -rsi ${MP_CC_DIR}../ 
	
${OBJECTDIR}/src/twiddleFactors.o: src/twiddleFactors.c  nbproject/Makefile-${CND_CONF}.mk
	@${MKDIR} "${OBJECTDIR}/src" 
	@${RM} ${OBJECTDIR}/src/twiddleFactors.o.d 
	@${RM} ${OBJECTDIR}/src/twiddleFactors.o 
	${MP_CC} $(MP_EXTRA_CC_PRE)  src/twiddleFactors.c  -o ${OBJECTDIR}/src/twiddleFactors.o  -c -mcpu=$(MP_PROCESSOR_OPTION)  -MMD -MF "${OBJECTDIR}/src/twiddleFactors.o.d"      -g -D__DEBUG -D__MPLAB_DEBUGGER_SIMULATOR=1    -omf=elf -O0 -msmart-io=1 -Wall -msfr-warn=off
	@${FIXDEPS} "${OBJECTDIR}/src/twiddleFactors.o.d" $(SILENT)  -rsi ${MP_CC_DIR}../ 
	
${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o: /Users/betio32/Documents/myGitHub/Senior-Design/dsPIC33/Joshua_dsPIC30F_FFT_Simulation_V2.X/src/chan1_input_cosine40kHz.c  nbproject/Makefile-${CND_CONF}.mk
	@${MKDIR} "${OBJECTDIR}/_ext/102761262" 
	@${RM} ${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o.d 
	@${RM} ${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o 
	${MP_CC} $(MP_EXTRA_CC_PRE)  /Users/betio32/Documents/myGitHub/Senior-Design/dsPIC33/Joshua_dsPIC30F_FFT_Simulation_V2.X/src/chan1_input_cosine40kHz.c  -o ${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o  -c -mcpu=$(MP_PROCESSOR_OPTION)  -MMD -MF "${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o.d"      -g -D__DEBUG -D__MPLAB_DEBUGGER_SIMULATOR=1    -omf=elf -O0 -msmart-io=1 -Wall -msfr-warn=off
	@${FIXDEPS} "${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o.d" $(SILENT)  -rsi ${MP_CC_DIR}../ 
	
else
${OBJECTDIR}/src/main_FFTExample.o: src/main_FFTExample.c  nbproject/Makefile-${CND_CONF}.mk
	@${MKDIR} "${OBJECTDIR}/src" 
	@${RM} ${OBJECTDIR}/src/main_FFTExample.o.d 
	@${RM} ${OBJECTDIR}/src/main_FFTExample.o 
	${MP_CC} $(MP_EXTRA_CC_PRE)  src/main_FFTExample.c  -o ${OBJECTDIR}/src/main_FFTExample.o  -c -mcpu=$(MP_PROCESSOR_OPTION)  -MMD -MF "${OBJECTDIR}/src/main_FFTExample.o.d"        -g -omf=elf -O0 -msmart-io=1 -Wall -msfr-warn=off
	@${FIXDEPS} "${OBJECTDIR}/src/main_FFTExample.o.d" $(SILENT)  -rsi ${MP_CC_DIR}../ 
	
${OBJECTDIR}/src/twiddleFactors.o: src/twiddleFactors.c  nbproject/Makefile-${CND_CONF}.mk
	@${MKDIR} "${OBJECTDIR}/src" 
	@${RM} ${OBJECTDIR}/src/twiddleFactors.o.d 
	@${RM} ${OBJECTDIR}/src/twiddleFactors.o 
	${MP_CC} $(MP_EXTRA_CC_PRE)  src/twiddleFactors.c  -o ${OBJECTDIR}/src/twiddleFactors.o  -c -mcpu=$(MP_PROCESSOR_OPTION)  -MMD -MF "${OBJECTDIR}/src/twiddleFactors.o.d"        -g -omf=elf -O0 -msmart-io=1 -Wall -msfr-warn=off
	@${FIXDEPS} "${OBJECTDIR}/src/twiddleFactors.o.d" $(SILENT)  -rsi ${MP_CC_DIR}../ 
	
${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o: /Users/betio32/Documents/myGitHub/Senior-Design/dsPIC33/Joshua_dsPIC30F_FFT_Simulation_V2.X/src/chan1_input_cosine40kHz.c  nbproject/Makefile-${CND_CONF}.mk
	@${MKDIR} "${OBJECTDIR}/_ext/102761262" 
	@${RM} ${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o.d 
	@${RM} ${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o 
	${MP_CC} $(MP_EXTRA_CC_PRE)  /Users/betio32/Documents/myGitHub/Senior-Design/dsPIC33/Joshua_dsPIC30F_FFT_Simulation_V2.X/src/chan1_input_cosine40kHz.c  -o ${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o  -c -mcpu=$(MP_PROCESSOR_OPTION)  -MMD -MF "${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o.d"        -g -omf=elf -O0 -msmart-io=1 -Wall -msfr-warn=off
	@${FIXDEPS} "${OBJECTDIR}/_ext/102761262/chan1_input_cosine40kHz.o.d" $(SILENT)  -rsi ${MP_CC_DIR}../ 
	
endif

# ------------------------------------------------------------------------------------
# Rules for buildStep: assemble
ifeq ($(TYPE_IMAGE), DEBUG_RUN)
else
endif

# ------------------------------------------------------------------------------------
# Rules for buildStep: assemblePreproc
ifeq ($(TYPE_IMAGE), DEBUG_RUN)
else
endif

# ------------------------------------------------------------------------------------
# Rules for buildStep: link
ifeq ($(TYPE_IMAGE), DEBUG_RUN)
dist/${CND_CONF}/${IMAGE_TYPE}/Joshua_dsPIC30F_FFT_Simulation_V2.X.${IMAGE_TYPE}.${OUTPUT_SUFFIX}: ${OBJECTFILES}  nbproject/Makefile-${CND_CONF}.mk  /Applications/microchip/xc16/v1.25/lib/libdsp-elf.a  
	@${MKDIR} dist/${CND_CONF}/${IMAGE_TYPE} 
	${MP_CC} $(MP_EXTRA_LD_PRE)  -o dist/${CND_CONF}/${IMAGE_TYPE}/Joshua_dsPIC30F_FFT_Simulation_V2.X.${IMAGE_TYPE}.${OUTPUT_SUFFIX}  ${OBJECTFILES_QUOTED_IF_SPACED}    /Applications/microchip/xc16/v1.25/lib/libdsp-elf.a  -mcpu=$(MP_PROCESSOR_OPTION)        -D__DEBUG -D__MPLAB_DEBUGGER_SIMULATOR=1  -omf=elf     -Wl,,--defsym=__MPLAB_BUILD=1,--defsym=__MPLAB_DEBUG=1,--defsym=__DEBUG=1,--defsym=__MPLAB_DEBUGGER_SIMULATOR=1,$(MP_LINKER_FILE_OPTION),--stack=16,--check-sections,--data-init,--pack-data,--handles,--isr,--no-gc-sections,--fill-upper=0,--stackguard=16,--no-force-link,--smart-io,-Map="${DISTDIR}/${PROJECTNAME}.${IMAGE_TYPE}.map",--report-mem,--memorysummary,dist/${CND_CONF}/${IMAGE_TYPE}/memoryfile.xml$(MP_EXTRA_LD_POST) 
	
else
dist/${CND_CONF}/${IMAGE_TYPE}/Joshua_dsPIC30F_FFT_Simulation_V2.X.${IMAGE_TYPE}.${OUTPUT_SUFFIX}: ${OBJECTFILES}  nbproject/Makefile-${CND_CONF}.mk  /Applications/microchip/xc16/v1.25/lib/libdsp-elf.a 
	@${MKDIR} dist/${CND_CONF}/${IMAGE_TYPE} 
	${MP_CC} $(MP_EXTRA_LD_PRE)  -o dist/${CND_CONF}/${IMAGE_TYPE}/Joshua_dsPIC30F_FFT_Simulation_V2.X.${IMAGE_TYPE}.${DEBUGGABLE_SUFFIX}  ${OBJECTFILES_QUOTED_IF_SPACED}    /Applications/microchip/xc16/v1.25/lib/libdsp-elf.a  -mcpu=$(MP_PROCESSOR_OPTION)        -omf=elf -Wl,,--defsym=__MPLAB_BUILD=1,$(MP_LINKER_FILE_OPTION),--stack=16,--check-sections,--data-init,--pack-data,--handles,--isr,--no-gc-sections,--fill-upper=0,--stackguard=16,--no-force-link,--smart-io,-Map="${DISTDIR}/${PROJECTNAME}.${IMAGE_TYPE}.map",--report-mem,--memorysummary,dist/${CND_CONF}/${IMAGE_TYPE}/memoryfile.xml$(MP_EXTRA_LD_POST) 
	${MP_CC_DIR}/xc16-bin2hex dist/${CND_CONF}/${IMAGE_TYPE}/Joshua_dsPIC30F_FFT_Simulation_V2.X.${IMAGE_TYPE}.${DEBUGGABLE_SUFFIX} -a  -omf=elf  
	
endif


# Subprojects
.build-subprojects:


# Subprojects
.clean-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r build/default
	${RM} -r dist/default

# Enable dependency checking
.dep.inc: .depcheck-impl

DEPFILES=$(shell "${PATH_TO_IDE_BIN}"mplabwildcard ${POSSIBLE_DEPFILES})
ifneq (${DEPFILES},)
include ${DEPFILES}
endif
