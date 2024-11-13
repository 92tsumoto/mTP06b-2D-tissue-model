#INCLUDES = syspara.h

bin_PROGRAMS = 2D_TP06b
2D_TP06b_CFLAGS = -O3 -xCORE-AVX2 -qopt-zmm-usage=high -qmkl -parallel -Winline -inline-level=2 -ip -qopt-mem-layout-trans=3 -pg -qopt-report-phase=all -qopt-report-annotate
#2D_TP06b_CFLAGS = -fcommon -fiopenmp -fopenmp-targets=spir64 -O3 -xCORE-AVX2 -qopt-zmm-usage=high -qmkl=parallel -Winline -inline-level=2 -ip -qopt-mem-layout-trans=3 -pg -qopt-report=3 -qopt-report-annotate -flto -Wunknown-pragmas -qopt-report-stdout
#2D_Kurata05_CFLAGS = -O2 -axCORE-AVX2 -parallel -qopt-report-phase=all -qopt-report=2

#2D_Kurata05_LDADD = -L/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
2D_TP06b_LDADD = -L/opt/intel/oneapi/mkl/lib/intel64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core

2D_TP06b_SOURCES = ./src/syspara.h ./src/TP06b.c ./src/exptab.c ./src/mem.c ./src/vars_para.c ./src/system.c

