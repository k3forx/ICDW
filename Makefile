FC=gfortran -O3
FC2=-fbounds-check -ffpe-trap=invalid,zero,overflow -fbacktrace
FC3=-Wall -fbounds-check -O -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace

all : ICDW_CL ICDW_HMC jackbin

ICDW_CL : mt19937.f90 ICDW_CL.F90
#	$(FC) $^ -o $@ -fno-range-check
	$(FC) $^ -o $@ -fno-range-check $(FC2)

ICDW_HMC : mt19937.f90 ICDW_HMC.F90
	$(FC) $^ -o $@ -fno-range-check
#	$(FC) -D_CHECK_REVERSE $^ -o $@ -fno-range-check $(FC2)

jackbin : jackbin.f
	$(FC) $^ -o $@

clean :
	rm -f *.o *.mod *~ ICDW_CL ICDW_HMC jackbin

# -Wall : 全てのコンパイル時警告メッセージを有効にする
# -pedantic -std=f95 : 標準外機能の利用を警告
# -fbounds-check : 配列の領域外参照を検出
# -O -Wuninitialized : 初期化されていない変数を検出
# -ffpe-trap=invalid,zero,overflow : 浮動小数点例外発生時に異常終了
# -fbacktrace : 異常終了時にソースコードの行番号を表示
