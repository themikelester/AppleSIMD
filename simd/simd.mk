#------------------------------------------------------------------------------
# Notes:  Simd Makefile include
#
# Author: Jeffrey Exterkate
# Date C: 6-25-2019
#------------------------------------------------------------------------------

SIMD_LIB	:= $(TGC_3RDPARTY_DIR)/Lib/libSimd$(CONFIG)$(PLATFORM)$(LIBRARY_SUFFIX)

#----------------------------------------------------------
# Additional flags
#----------------------------------------------------------
CXXFLAGS 	+= -I$(TGC_3RDPARTY_DIR)/clang
LIBS		+= $(SIMD_LIB)

#-------------------------------------------------------------
# Rules
#-------------------------------------------------------------
$(EXE) : $(SIMD_LIB)

pre_build : build_simd
build_simd :
	@$(MAKE) build --no-print-directory -C$(TGC_3RDPARTY_DIR)/clang/simd TGC_DIR=../../../..

pre_clean : clean_simd
clean_simd :
	@$(MAKE) clean --no-print-directory -C$(TGC_3RDPARTY_DIR)/clang/simd TGC_DIR=../../../..
