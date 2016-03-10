################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../nomad_lib.cpp 

OBJS += \
./nomad_lib.o 

CPP_DEPS += \
./nomad_lib.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel Intel(R) 64 C++ Compiler '
	icpc -g -I/usr/local/openmpi-1.6.2_intel13/include -I/data/cees/yiminliu/global/include/simlym -I/data/cees/yiminliu/global/include/lapack -I/data/cees/yiminliu/global/include/opca -I/data/cees/yiminliu/global/include/nomad -I/data/cees/yiminliu/global/include/armadillo -I/home/yiminliu/boost_1_56_0 -DUSE_MPI -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


