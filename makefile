
all:
	mpic++ -DUSEMPI -I /home/yiminliu/boost_1_56_0/ -I /data/cees/yiminliu/global/include/opca -I /data/cees/yiminliu/global/include/simlym -I /data/cees/yiminliu/global/include/nomad -L/data/cees/yiminliu/global/lib -L/home/yiminliu/nomad.3.6.2/lib -L/usr/local/openmpi-1.4.3_intel_14/lib -o "OPCA_HM_NOMAD_SIMLYM" nomad_lib.cpp -lnomadMPI -lsimlym -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lopca -llapack -lblas -lf2c -ltmg -lnlopt -lm -lmpi

