#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include "fft_wrapper.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Q31 精確轉換實現
int32_t float_to_q31(float val) {
    double scaled_val = (double)val * Q31_SCALE;
    if (scaled_val >= 2147483647.0) return 2147483647;
    if (scaled_val <= -2147483648.0) return -2147483648;
    return (int32_t)round(scaled_val);
}

float q31_to_float(int32_t val) {
    return (float)((double)val / Q31_SCALE);
}

// FFT 實現
int perform_fft_32_32(const float* input, int32_t* output, int* shifts, int scaling_option) {
    if (!input || !output || !shifts) return -1;

    // 獲取FFT句柄
    fft_handle_t handle_fft = cfft32_128;
    if (!handle_fft) return -2;

    // 轉換輸入為Q31並檢查飽和
    int saturated_count = 0;
    int32_t min_q31 = 2147483647, max_q31 = -2147483647 - 1;
    
    for (int i = 0; i < FFT_N; i++) {
        int32_t real_q31 = float_to_q31(input[i]);
        output[2 * i] = real_q31;
        output[2 * i + 1] = 0;

        if (input[i] >= 1.0f || input[i] < -1.0f) {
            if (real_q31 == 2147483647 || real_q31 == (-2147483647 - 1)) {
                saturated_count++;
            }
        }
        if (real_q31 < min_q31) min_q31 = real_q31;
        if (real_q31 > max_q31) max_q31 = real_q31;
    }

    // 執行FFT
    *shifts = fft_cplx32x32(output, output, handle_fft, scaling_option);
    return (*shifts >= 0) ? 0 : -3;
}

// IFFT 實現
int perform_ifft_32_32(const int32_t* input, float* output, int fft_shifts, int* shifts, int scaling_option) {
    if (!input || !output || !shifts) return -1;

    // 獲取IFFT句柄
    fft_handle_t handle_ifft = cifft32_128;
    if (!handle_ifft) return -2;

    // 分配IFFT輸出緩衝區
    int32_t* q31_output = (int32_t*)malloc(2 * FFT_N * sizeof(int32_t));
    if (!q31_output) return -4;

    // 複製輸入到工作緩衝區
    memcpy(q31_output, input, 2 * FFT_N * sizeof(int32_t));

    // 執行IFFT
    *shifts = ifft_cplx32x32(q31_output, q31_output, handle_ifft, scaling_option);
    if (*shifts < 0) {
        free(q31_output);
        return -3;
    }

    // 計算縮放並轉換回浮點
    int log2N = 7;  // log2(128)
    int total_shifts = fft_shifts + *shifts;
    double scale = ldexp(1.0, -(total_shifts));

    // 轉換回浮點並應用縮放
    for (int i = 0; i < FFT_N; i++) {
        output[i] = (float)(((double)q31_output[2 * i] / Q31_SCALE) * scale);
    }

    free(q31_output);
    return 0;
}

// 獲取波形推薦的預縮放因子
float get_recommended_prescale(signal_type_t type) {
    switch (type) {
        case WAVEFORM_SINGLE_COSINE:
            return 0.020f;  // 降低以確保不超過閾值
        case WAVEFORM_IMPULSE:
            return 0.015f;  // 顯著降低以處理尖峰
        case WAVEFORM_AM:
            return 0.018f;  // AM信號需要更小的縮放
        case WAVEFORM_GAUSSIAN:
            return 0.015f;  // 高斯脈衝也需要更小的縮放
        default:
            return 0.025f;  // 其他波形使用默認值
    }
}

// 信號生成實現 - 完整的9種波形支持
int generate_test_signal(float* buffer, signal_type_t type, float prescale) {
    if (!buffer) return -1;
    
    // 如果提供的prescale為預設值，使用推薦值
    if (fabsf(prescale - DEFAULT_PRESCALE_FACTOR) < 1e-6f) {
        prescale = get_recommended_prescale(type);
    }

    switch (type) {
        case WAVEFORM_DUAL_SINE: {
            const float FREQ1 = 500.0f, FREQ2 = 1500.0f;
            for (int i = 0; i < FFT_N; i++) {
                double time = (double)i / SAMPLE_RATE;
                buffer[i] = ((float)(0.7 * sin(2.0 * M_PI * FREQ1 * time) + 
                                   0.3 * sin(2.0 * M_PI * FREQ2 * time))) * prescale;
            }
            break;
        }
        case WAVEFORM_SINGLE_COSINE: {
            const float FREQ = 500.0f, AMP = 0.8f;
            for (int i = 0; i < FFT_N; i++) {
                double time = (double)i / SAMPLE_RATE;
                buffer[i] = AMP * cosf(2.0f * M_PI * FREQ * time) * prescale;
            }
            break;
        }
        case WAVEFORM_SQUARE: {
            const float FREQ = 500.0f, AMP = 0.5f;
            for (int i = 0; i < FFT_N; i++) {
                double time = (double)i / SAMPLE_RATE;
                float phase = 2.0f * M_PI * FREQ * time;
                buffer[i] = (fmod(phase, 2.0 * M_PI) < M_PI ? AMP : -AMP) * prescale;
            }
            break;
        }
        case WAVEFORM_IMPULSE: {
            const int IMPULSE_POS = 32;
            memset(buffer, 0, FFT_N * sizeof(float));
            buffer[IMPULSE_POS] = prescale;
            break;
        }
        case WAVEFORM_AM: {
            const float FC = 1000.0f, FM = 200.0f;
            const float M = 0.8f, AMP = 0.6f;
            for (int i = 0; i < FFT_N; i++) {
                double time = (double)i / SAMPLE_RATE;
                float mod = 1.0f + M * cosf(2.0f * M_PI * FM * time);
                float carrier = cosf(2.0f * M_PI * FC * time);
                buffer[i] = AMP * mod * carrier * prescale;
            }
            break;
        }
        case WAVEFORM_FM: {
            const float FC = 1000.0f, FM = 200.0f;
            const float DELTA_F = 50.0f;
            const float BETA = DELTA_F / FM;
            const float AMP = 0.7f;
            for (int i = 0; i < FFT_N; i++) {
                double time = (double)i / SAMPLE_RATE;
                float term1 = AMP * cosf(2.0f * M_PI * FC * time);
                float term2 = -(AMP * BETA / 2.0f) * cosf(2.0f * M_PI * (FC - FM) * time);
                float term3 = (AMP * BETA / 2.0f) * cosf(2.0f * M_PI * (FC + FM) * time);
                buffer[i] = (term1 + term2 + term3) * prescale;
            }
            break;
        }
        case WAVEFORM_CHIRP: {
            const float F0 = 100.0f, F1 = 1000.0f;
            const float AMP = 0.7f;
            const double T = (double)FFT_N / SAMPLE_RATE;
            const double k = (F1 - F0) / T;
            for (int i = 0; i < FFT_N; i++) {
                double time = (double)i / SAMPLE_RATE;
                float phase = 2.0f * M_PI * (F0 * time + k * time * time / 2.0);
                buffer[i] = AMP * cosf(phase) * prescale;
            }
            break;
        }
        case WAVEFORM_GAUSSIAN: {
            const float AMP = 1.0f;
            const double T0 = (double)(FFT_N / 2) / SAMPLE_RATE;
            const double SIGMA = (double)(FFT_N / 16) / SAMPLE_RATE;
            for (int i = 0; i < FFT_N; i++) {
                double time = (double)i / SAMPLE_RATE;
                double exponent = -pow((time - T0) / (2.0 * SIGMA), 2.0);
                buffer[i] = AMP * exp(exponent) * prescale;
            }
            break;
        }
        case WAVEFORM_WHITE_NOISE: {
            const float AMP = 0.8f;
            srand(0);
            for (int i = 0; i < FFT_N; i++) {
                float rand_val = ((float)rand() / (float)RAND_MAX) * 2.0f - 1.0f;
                buffer[i] = AMP * rand_val * prescale;
            }
            break;
        }
        default:
            return -2;
    }
    return 0;
}

// 工具函數實現
void print_signal_preview(const float* signal, int length, int preview_len) {
    printf("\nSignal Preview (first %d points):\n", preview_len);
    float max_val = 0;
    for (int i = 0; i < preview_len; i++) {
        if (fabsf(signal[i]) > max_val) max_val = fabsf(signal[i]);
    }
    for (int i = 0; i < preview_len; i++) {
        printf("%4d |", i);
        float normalized = fabsf(signal[i]) / max_val;
        int pos = (int)(normalized * 30);
        for (int j = 0; j < 30; j++) {
            printf("%s", (j == pos) ? "*" : " ");
        }
        printf(" %.6f\n", signal[i]);
    }
}

double calculate_signal_error_metrics(const float* original, const float* recovered, 
                                    int length, float* max_error, double* rmse) {
    double sum_squared_error = 0.0;
    *max_error = 0.0f;
    
    for (int i = 0; i < length; i++) {
        float error = fabsf(original[i] - recovered[i]);
        if (error > *max_error) *max_error = error;
        sum_squared_error += error * error;
    }
    
    *rmse = sqrt(sum_squared_error / length);
    return *rmse;
}

// 靜態庫測試實現
int test_static_library(static_lib_test_result_t* results, int* num_tests) {
    *num_tests = 0;
    int test_index = 0;
    
    // 測試1: FFT句柄可用性測試
    results[test_index].test_name = "FFT句柄測試";
    results[test_index].description = "檢查FFT靜態庫句柄是否可用";
    fft_handle_t handle_fft = cfft32_128;
    if (handle_fft) {
        results[test_index].test_passed = 1;
        results[test_index].error_message = "成功";
        results[test_index].error_value = 0.0f;
    } else {
        results[test_index].test_passed = 0;
        results[test_index].error_message = "FFT句柄無效";
        results[test_index].error_value = -1.0f;
    }
    test_index++;

    // 測試2: IFFT句柄可用性測試
    results[test_index].test_name = "IFFT句柄測試";
    results[test_index].description = "檢查IFFT靜態庫句柄是否可用";
    fft_handle_t handle_ifft = cifft32_128;
    if (handle_ifft) {
        results[test_index].test_passed = 1;
        results[test_index].error_message = "成功";
        results[test_index].error_value = 0.0f;
    } else {
        results[test_index].test_passed = 0;
        results[test_index].error_message = "IFFT句柄無效";
        results[test_index].error_value = -1.0f;
    }
    test_index++;

    // 測試3: 基本FFT運算測試
    results[test_index].test_name = "基本FFT運算測試";
    results[test_index].description = "執行簡單信號的FFT轉換";
    float input[FFT_N] = {0};
    int32_t output[2 * FFT_N] = {0};
    int shifts = 0;
    
    // 生成簡單的測試信號
    input[0] = 1.0f * DEFAULT_PRESCALE_FACTOR;
    int ret = perform_fft_32_32(input, output, &shifts, 3);
    
    if (ret == 0 && shifts >= 0) {
        results[test_index].test_passed = 1;
        results[test_index].error_message = "成功";
        results[test_index].error_value = 0.0f;
    } else {
        results[test_index].test_passed = 0;
        results[test_index].error_message = "FFT運算失敗";
        results[test_index].error_value = (float)ret;
    }
    test_index++;

    // 測試4: Q31轉換精度測試
    results[test_index].test_name = "Q31轉換精度測試";
    results[test_index].description = "測試浮點數與Q31格式轉換的精度";
    float test_value = 0.5f;
    int32_t q31_value = float_to_q31(test_value);
    float recovered_value = q31_to_float(q31_value);
    float conversion_error = fabsf(test_value - recovered_value);
    
    if (conversion_error < 1e-6f) {
        results[test_index].test_passed = 1;
        results[test_index].error_message = "成功";
        results[test_index].error_value = conversion_error;
    } else {
        results[test_index].test_passed = 0;
        results[test_index].error_message = "Q31轉換誤差過大";
        results[test_index].error_value = conversion_error;
    }
    test_index++;

    *num_tests = test_index;
    return 0;
}