#!/bin/bash

export XTENSA_SYSTEM=/home/yu-chen/xtensa/XtDevTools/install/tools/RI-2021.6-linux/XtensaTools/config
export XTENSA_BASE=/home/yu-chen/xtensa/XtDevTools/install
export XTENSA_CORE=hifi4_ss_spfpu_7

DATETIME=$(date +%Y%m%d_%H%M%S)
LOG_FILE="test_results_${DATETIME}.txt"

echo "=== FFT/IFFT 波形測試開始 (${DATETIME}) ===" | tee -a "$LOG_FILE"

# 要測試的波形類型（1-9）
WAVEFORMS=(1 2 3 4 5 6 7 8 9)

for waveform in "${WAVEFORMS[@]}"; do
    echo -e "\n\n=== 測試波形類型 $waveform ===" | tee -a "$LOG_FILE"
    
    # 修改波形類型
    sed -i "s/#define WAVEFORM_TYPE .*/#define WAVEFORM_TYPE $waveform/" fftifft_32_32_test_other_signal1.c
    
    # 編譯
    echo "編譯程式..." | tee -a "$LOG_FILE"
    xt-clang fftifft_32_32_test_other_signal1.c -o fftifft_test_program_q15 \
        -I/home/yu-chen/下載/HiFi4_VFPU_Library_v4_1_0/hifi4_library/include \
        -L. -lAll_fft_NatureDSP_Signal-hifi4_ss_spfpu_7_llvm-Xtensa-release -lm 2>&1 | tee -a "$LOG_FILE"
    
    # 運行測試（設置30秒超時）
    echo "運行測試..." | tee -a "$LOG_FILE"
    timeout 30s xt-run fftifft_test_program_q15 2>&1 | tee -a "$LOG_FILE"
    
    # 檢查是否超時
    if [ $? -eq 124 ]; then
        echo "警告：波形類型 $waveform 執行超過30秒，跳過" | tee -a "$LOG_FILE"
    fi
    
    echo "=== 波形類型 $waveform 測試完成 ===" | tee -a "$LOG_FILE"
done

echo -e "\n=== 所有測試完成 ===" | tee -a "$LOG_FILE"
