function(set_gpu)
	execute_process(
		COMMAND nvidia-smi --query-gpu=compute_cap --format=csv
		COMMAND tail -n 1
		COMMAND tr -d .
		OUTPUT_VARIABLE GPU_CC
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	execute_process(
		COMMAND nvcc --version
		COMMAND grep release
		COMMAND cut -d " " -f 5
		COMMAND cut -d "," -f 1
		OUTPUT_VARIABLE CUDA_VERSION
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
	message(STATUS "GPU Compute Capability: ${GPU_CC}")
	message(STATUS "CUDA Version: ${CUDA_VERSION}")
	set(GPU_CC ${GPU_CC} PARENT_SCOPE)
	set(CUDA_VERSION ${CUDA_VERSION} PARENT_SCOPE)
endfunction(set_gpu)
