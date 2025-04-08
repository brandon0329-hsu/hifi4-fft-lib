#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "fft_wrapper.h"

const char* get_waveform_name(signal_type_t type) {
    switch (type) {
        case WAVEFORM_DUAL_SINE: return "雙正弦波";
        case WAVEFORM_SINGLE_COSINE: return "單餘弦波";
        case WAVEFORM_SQUARE: return "方波";
        case WAVEFORM_IMPULSE: return "衝激函數";
        case WAVEFORM_AM: return "AM信號";
        case WAVEFORM_FM: return "窄帶FM";
        case WAVEFORM_CHIRP: return "線性掃頻";
        case WAVEFORM_GAUSSIAN: return "高斯脈衝";
        case WAVEFORM_WHITE_NOISE: return "白噪聲";
        default: return "未知";
    }
}

fft_result_t test_waveform(signal_type_t type, float prescale) {
    fft_result_t result = {0};
    result.prescale_factor = prescale;

    // 分配記憶體
    float* input_signal = (float*)malloc(FFT_N * sizeof(float));
    int32_t* fft_output = (int32_t*)malloc(2 * FFT_N * sizeof(int32_t));
    float* recovered_signal = (float*)malloc(FFT_N * sizeof(float));

    if (!input_signal || !fft_output || !recovered_signal) {
        fprintf(stderr, "記憶體分配失敗!\n");
        result.success = -1;
        goto cleanup;
    }

    // 生成測試信號
    int ret = generate_test_signal(input_signal, type, prescale);
    if (ret != 0) {
        fprintf(stderr, "信號生成失敗: %d\n", ret);
        result.success = -2;
        goto cleanup;
    }

    // 執行 FFT
    ret = perform_fft_32_32(input_signal, fft_output, &result.fft_shifts, 3);
    if (ret != 0) {
        fprintf(stderr, "FFT 執行失敗: %d\n", ret);
        result.success = -3;
        goto cleanup;
    }

    // 執行 IFFT
    ret = perform_ifft_32_32(fft_output, recovered_signal, result.fft_shifts, 
                            &result.ifft_shifts, 3);
    if (ret != 0) {
        fprintf(stderr, "IFFT 執行失敗: %d\n", ret);
        result.success = -4;
        goto cleanup;
    }

    // 計算誤差指標
    float max_error = 0.0f;
    double rmse = 0.0;
    calculate_signal_error_metrics(input_signal, recovered_signal, FFT_N, 
                                 &max_error, &rmse);
    result.max_error = max_error;
    result.rmse = rmse;
    result.success = (max_error <= DEFAULT_TOLERANCE) ? 1 : 0;

cleanup:
    if (input_signal) free(input_signal);
    if (fft_output) free(fft_output);
    if (recovered_signal) free(recovered_signal);

    return result;
}

void save_test_report(const char* filename, fft_result_t* results, int count) {
    FILE* f = fopen(filename, "w");
    if (!f) return;

    fprintf(f, "FFT/IFFT 測試報告 (生成於 %s)\n", __DATE__);
    fprintf(f, "================================================\n");
    fprintf(f, "測試配置:\n");
    fprintf(f, "- FFT點數: %d\n", FFT_N);
    fprintf(f, "- 取樣率: %.1f Hz\n", SAMPLE_RATE);
    fprintf(f, "- 預設縮放因子: %.3f\n", DEFAULT_PRESCALE_FACTOR);
    fprintf(f, "- 容許誤差: %.3f\n", DEFAULT_TOLERANCE);
    fprintf(f, "\n測試結果總結:\n");
    
    int passed = 0;
    for (int i = 0; i < count; i++) {
        if (results[i].success == 1) passed++;
    }
    
    fprintf(f, "總測試數: %d\n", count);
    fprintf(f, "通過測試: %d\n", passed);
    fprintf(f, "失敗測試: %d\n", count - passed);
    fprintf(f, "通過率: %.1f%%\n\n", (float)passed / count * 100);

    fprintf(f, "詳細測試結果:\n");
    fprintf(f, "%-15s | %-8s | %-12s | %-12s | %-12s | %-12s\n",
            "波形類型", "結果", "最大誤差", "RMSE", "FFT移位", "IFFT移位");
    fprintf(f, "----------------|----------|--------------|--------------|--------------|-------------\n");

    for (int i = 0; i < count; i++) {
        fprintf(f, "%-15s | %-8s | %12.6f | %12.6f | %12d | %12d\n",
                get_waveform_name(i + 1),
                results[i].success == 1 ? "通過" : "失敗",
                results[i].max_error,
                results[i].rmse,
                results[i].fft_shifts,
                results[i].ifft_shifts);
    }

    fprintf(f, "\n註: 測試通過標準為最大誤差 <= %.3f\n", DEFAULT_TOLERANCE);
    fclose(f);
}

void run_static_library_tests() {
    printf("\n=== 執行靜態庫測試 ===\n");
    static_lib_test_result_t results[10];  // 預留足夠空間
    int num_tests = 0;
    
    int ret = test_static_library(results, &num_tests);
    if (ret != 0) {
        printf("靜態庫測試執行失敗!\n");
        return;
    }
    
    // 打印測試結果
    printf("\n靜態庫測試結果:\n");
    printf("%-20s | %-40s | %-10s | %-15s\n", "測試名稱", "描述", "結果", "誤差值");
    printf("-------------------|------------------------------------------|------------|----------------\n");
    
    int passed = 0;
    for (int i = 0; i < num_tests; i++) {
        printf("%-20s | %-40s | %-10s | %15.6f\n",
               results[i].test_name,
               results[i].description,
               results[i].test_passed ? "通過" : "失敗",
               results[i].error_value);
        if (results[i].test_passed) passed++;
    }
    
    printf("\n總結:\n");
    printf("總測試數: %d\n", num_tests);
    printf("通過測試: %d\n", passed);
    printf("失敗測試: %d\n", num_tests - passed);
    printf("通過率: %.1f%%\n", (float)passed / num_tests * 100);
}

int main() {
    printf("=== FFT/IFFT 測試開始 ===\n");
    printf("FFT點數: %d, 取樣率: %.1f Hz\n", FFT_N, SAMPLE_RATE);

    // 首先執行靜態庫測試
    run_static_library_tests();
    
    // 準備存儲所有波形的測試結果
    fft_result_t results[9];
    int total_tests = 9;

    // 測試所有波形
    for (signal_type_t type = WAVEFORM_DUAL_SINE; type <= WAVEFORM_WHITE_NOISE; type++) {
        printf("\n=== 測試波形: %s ===\n", get_waveform_name(type));
        results[type - 1] = test_waveform(type, DEFAULT_PRESCALE_FACTOR);
        
        // 輸出即時結果
        if (results[type - 1].success == 1) {
            printf("測試通過 - 最大誤差: %.6f, RMSE: %.6f\n", 
                   results[type - 1].max_error, results[type - 1].rmse);
        } else if (results[type - 1].success == 0) {
            printf("測試失敗 - 最大誤差: %.6f (超過容許值 %.3f)\n", 
                   results[type - 1].max_error, DEFAULT_TOLERANCE);
        } else {
            printf("測試執行錯誤: %d\n", results[type - 1].success);
        }
    }

    // 生成測試報告
    char filename[100];
    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    sprintf(filename, "test_report_%d%02d%02d_%02d%02d%02d.txt",
            t->tm_year + 1900, t->tm_mon + 1, t->tm_mday,
            t->tm_hour, t->tm_min, t->tm_sec);
    
    save_test_report(filename, results, total_tests);
    printf("\n測試報告已保存至: %s\n", filename);

    // 計算總結果
    int passed_tests = 0;
    for (int i = 0; i < total_tests; i++) {
        if (results[i].success == 1) passed_tests++;
    }

    printf("\n=== 測試結果總結 ===\n");
    printf("總測試數: %d\n", total_tests);
    printf("通過測試: %d\n", passed_tests);
    printf("失敗測試: %d\n", total_tests - passed_tests);
    printf("通過率: %.1f%%\n", (float)passed_tests / total_tests * 100);

    return (passed_tests == total_tests) ? 0 : 1;
}