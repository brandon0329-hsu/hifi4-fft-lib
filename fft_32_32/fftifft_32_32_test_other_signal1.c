#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <complex.h>

#include "NatureDSP_types.h"
#include "NatureDSP_Signal_fft.h" // 確保包含 32x32 函數聲明

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- 基本配置 ---
#define FFT_N 128
#define Q31_SCALE 2147483647.0f // 使用 Q31 的縮放因子
#define SAMPLE_RATE 8000.0f

// --- 波形選擇 ---
// 1=雙正弦, 2=單餘弦, 3=方波, 4=衝激
// 5=AM信號, 6=窄帶FM近似, 7=線性掃頻(Chirp), 8=高斯脈衝, 9=白噪聲
#define WAVEFORM_TYPE 9

// --- 輸入預縮放 ---
// Q31 範圍更大，飽和風險降低，但仍建議使用預縮放控制幅度
#define INPUT_PRESCALE_FACTOR 0.05f

// --- 比較容差 (浮點單位) ---
// 32 位精度更高，理論上誤差可能更小，但先保持之前的容差
#define FLOAT_COMPARISON_TOLERANCE 0.01f

// --- 輔助函數：浮點數轉 Q31 ---
int32_t float_to_q31(float val) {
    double scaled_val = (double)val * Q31_SCALE; // 使用 double 提高中間精度
    // 進行飽和處理 (對應 int32_t)
    if (scaled_val >= 2147483647.0) return 2147483647;
    if (scaled_val <= -2147483648.0) return -2147483648;
    return (int32_t)round(scaled_val); // 使用 round 進行四捨五入
}

// --- 輔助函數：Q31 轉浮點數 ---
float q31_to_float(int32_t val) {
    return (float)((double)val / Q31_SCALE); // 使用 double 提高精度
}

// --- 浮點比較輔助函數 (保持不變) ---
double calculate_float_rmse_real(const float complex* arr_complex, const float* arr_real, size_t n) {
    double mse = 0.0; if (n == 0) return 0.0;
    for (size_t i = 0; i < n; ++i) { double diff = crealf(arr_complex[i]) - arr_real[i]; mse += diff * diff; }
    return sqrt(mse / n);
}
float calculate_float_max_abs_imag(const float complex* arr_complex, size_t n) {
    float max_abs_imag = 0.0f;
    for (size_t i = 0; i < n; ++i) { float abs_imag = fabsf(cimagf(arr_complex[i])); if (abs_imag > max_abs_imag) max_abs_imag = abs_imag; }
    return max_abs_imag;
}

// --- 輔助函數：打印訊號預覽 ---
void print_signal_preview(const float* signal, int n, int preview_length) {
    printf("\nSignal Preview (first %d points):\n", preview_length);
    float max_val = 0;
    // 找出最大值用於縮放
    for(int i = 0; i < preview_length; i++) {
        if(fabs(signal[i]) > max_val) max_val = fabs(signal[i]);
    }
    // 打印預覽
    for(int i = 0; i < preview_length; i++) {
        printf("%4d |", i);
        int pos = (int)(fabs(signal[i]) / max_val * 30); // 縮放到30個字符寬度
        for(int j = 0; j < 30; j++) {
            if(j == pos) printf("*");
            else printf(" ");
        }
        printf("\n");
    }
}

int main() {
    printf("--- NatureDSP Fixed-point FFT -> IFFT Demonstration (N=%d, 32x32) ---\n", FFT_N); // 標明 32x32
    printf("--- Verifying in FLOAT domain using scalingOption=3, Non-Inplace ---\n");

    // --- 1. 初始化與準備 ---
    float   *input_signal_float = NULL;      // 原始浮點時域信號 (預縮放後)
    int32_t *q31_input_fft = NULL;        // FFT 的 Q31 輸入緩衝區
    int32_t *q31_output_fft = NULL;       // FFT 的 Q31 輸出緩衝區 (也作為 IFFT 的輸入)
    int32_t *q31_output_ifft = NULL;      // IFFT 的 Q31 輸出緩衝區
    float complex *ifft_output_float_scaled = NULL; // 存放最終縮放後的浮點結果
    fft_handle_t handle_fft = NULL;
    fft_handle_t handle_ifft = NULL;
    int ret = 0;
    int total_shift_fft = 0;
    int total_shift_ifft = 0;
    int scalingOption = 3;

    printf("\n[Step 1] Allocating Memory (Separate Buffers, Q31)...\n");
    input_signal_float = (float *)malloc(FFT_N * sizeof(float));
    q31_input_fft = (int32_t *)malloc(2 * FFT_N * sizeof(int32_t)); // 使用 int32_t
    q31_output_fft = (int32_t *)malloc(2 * FFT_N * sizeof(int32_t)); // 使用 int32_t
    q31_output_ifft = (int32_t *)malloc(2 * FFT_N * sizeof(int32_t)); // 使用 int32_t
    ifft_output_float_scaled = (float complex *)malloc(FFT_N * sizeof(float complex));

    if (!input_signal_float || !q31_input_fft || !q31_output_fft || !q31_output_ifft || !ifft_output_float_scaled) {
        fprintf(stderr, "ERROR: Memory allocation failed!\n"); ret = 1; goto cleanup; }
    printf("Memory allocated successfully.\n");

    // --- Step 2: 生成選擇的浮點信號並預縮放 (邏輯不變) ---
#if WAVEFORM_TYPE == 1 // 雙正弦波
    printf("\n[Step 2] Generating Float Input Signal (Dual Sine) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float FREQ1 = 500.0f; const float FREQ2 = 1500.0f;
    for (unsigned int i = 0; i < FFT_N; i++) {
        double time = (double)i / SAMPLE_RATE;
        float raw_signal = (float)(0.7 * sin(2.0 * M_PI * FREQ1 * time) + 0.3 * sin(2.0 * M_PI * FREQ2 * time));
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#elif WAVEFORM_TYPE == 2 // 單一餘弦波
    printf("\n[Step 2] Generating Float Input Signal (Single Cosine) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float COS_FREQ = 500.0f; const float COS_AMP = 0.8f;
    for (unsigned int i = 0; i < FFT_N; i++) {
        double time = (double)i / SAMPLE_RATE;
        float raw_signal = COS_AMP * cosf(2.0f * M_PI * COS_FREQ * time);
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#elif WAVEFORM_TYPE == 3 // 方波
    printf("\n[Step 2] Generating Float Input Signal (Square Wave) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float SQ_FREQ = 500.0f;
    const float SQ_AMP = 0.5f;
    for (unsigned int i = 0; i < FFT_N; i++) {
        double time = (double)i / SAMPLE_RATE;
        float phase = 2.0f * M_PI * SQ_FREQ * time;
        float raw_signal = (fmod(phase, 2.0 * M_PI) < M_PI) ? SQ_AMP : -SQ_AMP;
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#elif WAVEFORM_TYPE == 4 // 衝激函數
    printf("\n[Step 2] Generating Float Input Signal (Impulse) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const int IMPULSE_POS = 32;
    const float IMPULSE_AMP = 1.0f;
    memset(input_signal_float, 0, FFT_N * sizeof(float));
    if (IMPULSE_POS >= 0 && IMPULSE_POS < FFT_N) {
        input_signal_float[IMPULSE_POS] = IMPULSE_AMP * INPUT_PRESCALE_FACTOR;
    }
#elif WAVEFORM_TYPE == 5 // AM信號
    printf("\n[Step 2] Generating Float Input Signal (AM) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float AM_FC = 1000.0f;
    const float AM_FM = 200.0f;
    const float AM_M = 0.8f;
    const float AM_AMP = 0.6f;
    for (unsigned int i = 0; i < FFT_N; i++) {
        double time = (double)i / SAMPLE_RATE;
        float mod_signal = 1.0f + AM_M * cosf(2.0f * M_PI * AM_FM * time);
        float carrier_signal = AM_AMP * cosf(2.0f * M_PI * AM_FC * time);
        float raw_signal = mod_signal * carrier_signal;
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#elif WAVEFORM_TYPE == 6 // 窄帶FM近似
    printf("\n[Step 2] Generating Float Input Signal (Narrowband FM) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float FM_FC = 1000.0f;
    const float FM_FM = 200.0f;
    const float FM_DELTA_F = 50.0f;
    const float FM_BETA = FM_DELTA_F / FM_FM;
    const float FM_AMP = 0.7f;
    printf("FM Beta: %.3f\n", FM_BETA);
    for (unsigned int i = 0; i < FFT_N; i++) {
        double time = (double)i / SAMPLE_RATE;
        float term1 = FM_AMP * cosf(2.0f * M_PI * FM_FC * time);
        float term2 = -(FM_AMP * FM_BETA / 2.0f) * cosf(2.0f * M_PI * (FM_FC - FM_FM) * time);
        float term3 = (FM_AMP * FM_BETA / 2.0f) * cosf(2.0f * M_PI * (FM_FC + FM_FM) * time);
        float raw_signal = term1 + term2 + term3;
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#elif WAVEFORM_TYPE == 7 // 線性掃頻
    printf("\n[Step 2] Generating Float Input Signal (Linear Chirp) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float CHIRP_F0 = 100.0f;
    const float CHIRP_F1 = 1000.0f;
    const float CHIRP_AMP = 0.7f;
    const double T = (double)FFT_N / SAMPLE_RATE;
    const double k = (CHIRP_F1 - CHIRP_F0) / T;
    for (unsigned int i = 0; i < FFT_N; i++) {
        double time = (double)i / SAMPLE_RATE;
        float phase = 2.0f * M_PI * (CHIRP_F0 * time + k * time * time / 2.0);
        float raw_signal = CHIRP_AMP * cosf(phase);
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#elif WAVEFORM_TYPE == 8 // 高斯脈衝
    printf("\n[Step 2] Generating Float Input Signal (Gaussian Pulse) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float GAUSS_AMP = 1.0f;
    const double GAUSS_T0 = (double)(FFT_N / 2) / SAMPLE_RATE;
    const double GAUSS_SIGMA = (double)(FFT_N / 16) / SAMPLE_RATE;
    for (unsigned int i = 0; i < FFT_N; i++) {
        double time = (double)i / SAMPLE_RATE;
        double exponent = -pow((time - GAUSS_T0) / (2.0 * GAUSS_SIGMA), 2.0);
        float raw_signal = GAUSS_AMP * exp(exponent);
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#elif WAVEFORM_TYPE == 9 // 白噪聲
    printf("\n[Step 2] Generating Float Input Signal (White Noise) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float NOISE_AMP = 0.8f;
    srand(0);
    for (unsigned int i = 0; i < FFT_N; i++) {
        float rand_val = ((float)rand() / (float)RAND_MAX) * 2.0f - 1.0f;
        float raw_signal = NOISE_AMP * rand_val;
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#endif
    printf("Input signal generated and prescaled.\n");
    print_signal_preview(input_signal_float, FFT_N, 20);

    // --- Step 3: 轉換為 Q31 FFT 輸入 ---
    printf("\n[Step 3] Converting Prescaled Float to Q31 FFT Input...\n");
    int32_t min_q31 = 2147483647, max_q31 = -2147483647 - 1; int saturated_count = 0; // 使用 int32_t 範圍
    for (unsigned int i = 0; i < FFT_N; i++) {
        int32_t real_q31 = float_to_q31(input_signal_float[i]); // 調用 Q31 轉換函數
        q31_input_fft[2 * i + 0] = real_q31;
        q31_input_fft[2 * i + 1] = 0;
        // 檢查飽和 (基於浮點值)
        if (input_signal_float[i] >= 1.0f || input_signal_float[i] < -1.0f ) {
           if (real_q31 == 2147483647 || real_q31 == (-2147483647 - 1)){ saturated_count++; }
        }
        if (real_q31 < min_q31) min_q31 = real_q31; if (real_q31 > max_q31) max_q31 = real_q31;
    }
    // 打印 Q31 範圍
    printf("Q31 conversion complete. Q31 Input Range: [%d, %d]. Saturated samples: %d\n", min_q31, max_q31, saturated_count);
    if (saturated_count > 0) printf("WARNING: Saturation occurred during Q31 conversion.\n");

    // --- Step 4: 執行 FFT (32x32, 非原地) ---
    printf("\n[Step 4] Performing Fixed-point FFT (fft_cplx32x32, Non-Inplace)...\n");
    // *** 關鍵：替換為正確的 32x32 FFT 句柄！***
    handle_fft = cfft32_128; // <-- 假設的句柄名，請務必核對！
    if (!handle_fft) { fprintf(stderr, "ERROR: FFT handle cfft32_128 is NULL. Check library.\n"); ret = 1; goto cleanup; }

    total_shift_fft = fft_cplx32x32(q31_output_fft, q31_input_fft, handle_fft, scalingOption); // 調用 32x32 FFT

    if (total_shift_fft < 0) { fprintf(stderr, "ERROR: fft_cplx32x32 failed (%d).\n", total_shift_fft); ret = 1; goto cleanup; }
    printf("FFT finished. Total scaling shifts (FFT): %d\n", total_shift_fft);

    // --- Step 5: 執行 IFFT (32x32, 非原地) ---
    printf("\n[Step 5] Performing Fixed-point IFFT (ifft_cplx32x32, Non-Inplace)...\n");
    // *** 關鍵：替換為正確的 32x32 IFFT 句柄！***
    handle_ifft = cifft32_128; // <-- 假設的句柄名，請務必核對！
     if (!handle_ifft) { fprintf(stderr, "ERROR: IFFT handle cifft32_128 is NULL. Check library.\n"); ret = 1; goto cleanup; }

    total_shift_ifft = ifft_cplx32x32(q31_output_ifft, q31_output_fft, handle_ifft, scalingOption); // 調用 32x32 IFFT

    if (total_shift_ifft < 0) { fprintf(stderr, "ERROR: ifft_cplx32x32 failed (%d).\n", total_shift_ifft); ret = 1; }
    else {
        printf("IFFT finished. Total scaling shifts (IFFT): %d\n", total_shift_ifft);

        // --- Step 6: 轉換 IFFT Q31 輸出到浮點並應用整體縮放 ---
        printf("\n[Step 6] Converting IFFT Q31 Output to Float and Applying *Overall* Scaling...\n");
        int total_shifts = total_shift_fft + total_shift_ifft;
        int log2N = 0; if (FFT_N > 0) { int tempN = FFT_N; while (tempN > 1) { tempN >>= 1; log2N++; } }
        int shift_for_vecToFp64 = log2N - total_shifts;
        printf("Calculated equivalent total shift for vecToFp64: %d (log2N - total_shifts)\n", shift_for_vecToFp64);

        int32_t* ifft_q31_result_ptr = q31_output_ifft; // 使用 Q31 結果

        // 手動模擬 vecToFp64 的核心邏輯 (Q31 版本)
        for (unsigned int i = 0; i < FFT_N; i++) {
            // *** 關鍵：使用 -shift - 31 ***
            double real_fp64 = ldexp((double)ifft_q31_result_ptr[2 * i + 0], -shift_for_vecToFp64 - 31);
            double imag_fp64 = ldexp((double)ifft_q31_result_ptr[2 * i + 1], -shift_for_vecToFp64 - 31);
            ifft_output_float_scaled[i] = (float)real_fp64 + (float)imag_fp64 * I;
        }
        printf("Q31 to Float conversion and overall scaling complete.\n");

        // --- Step 7: 在浮點域比較 (邏輯不變) ---
        printf("\n[Step 7] Comparing Scaled IFFT Float Output with Original Prescaled Float Input...\n");
        // ... (比較代碼、表格打印、成功/失敗判斷 - 與之前相同) ...
        // --- 比較代碼開始 ---
        double rmse_real = calculate_float_rmse_real(ifft_output_float_scaled, input_signal_float, FFT_N);
        float max_abs_imag = calculate_float_max_abs_imag(ifft_output_float_scaled, FFT_N);
        float max_real_error = 0.0f;
        for(unsigned int i = 0; i < FFT_N; ++i) { float current_real_error = fabsf(crealf(ifft_output_float_scaled[i]) - input_signal_float[i]); if (current_real_error > max_real_error) max_real_error = current_real_error; }
        printf("\nOverall Comparison Metrics (Float units):\n");
        printf("  Real Part RMSE vs Expected: %e\n", rmse_real);
        printf("  Max Abs Real Part Error:    %e\n", max_real_error);
        printf("  Max Abs Imaginary Part:     %e\n", max_abs_imag);
        printf("\nFloat Comparison Table (first/last 8 samples):\n");
        printf("Index | Orig Prescaled(F)    | Scaled IFFT(F)       | Abs. Error       | Rel. Error (%%)\n");
        printf("------|----------------------|----------------------|------------------|------------------\n");
        unsigned int print_count = 8;
        for(unsigned int i=0; i < FFT_N; ++i) {
            if (i < print_count || i >= FFT_N - print_count) {
                 float real_exp = input_signal_float[i]; float imag_exp = 0.0f;
                 float real_out = crealf(ifft_output_float_scaled[i]); float imag_out = cimagf(ifft_output_float_scaled[i]);
                 float abs_err_real = fabsf(real_out - real_exp); float abs_err_imag = fabsf(imag_out - imag_exp);
                 float rel_err_real = 0.0f; if (fabsf(real_exp) > 1e-9) { rel_err_real = (abs_err_real / fabsf(real_exp)) * 100.0f; } else if (fabsf(real_out) > 1e-9) { rel_err_real = 100.0f; }
                 float rel_err_imag = 0.0f; if (fabsf(INPUT_PRESCALE_FACTOR) > 1e-9) { float denom = fabsf(INPUT_PRESCALE_FACTOR); rel_err_imag = (abs_err_imag / (denom+1e-9)) * 100.0f; }
                 printf(" %4u | R:%+8.5f I:%+8.5f | R:%+8.5f I:%+8.5f | R:%8.2e I:%8.2e | R:%7.2f I:%7.2f\n", i, real_exp, imag_exp, real_out, imag_out, abs_err_real, abs_err_imag, rel_err_real, rel_err_imag);
                 if (i == print_count - 1 && print_count < FFT_N / 2) { printf("  ... (middle samples omitted) ...\n"); }
            }
        }
        if (max_real_error <= FLOAT_COMPARISON_TOLERANCE && max_abs_imag <= FLOAT_COMPARISON_TOLERANCE) {
            printf("\nSUCCESS: Max errors in float domain are within tolerance (%.4f).\n", FLOAT_COMPARISON_TOLERANCE);
            printf("This demonstrates IFFT(FFT(x)) [32x32] recovers the signal in float domain when overall scaling (log2N - total_shifts) is compensated.\n");
        } else {
            printf("\nFAILURE: Errors in float domain exceed tolerance (%.4f).\n", FLOAT_COMPARISON_TOLERANCE);
             if (max_real_error > FLOAT_COMPARISON_TOLERANCE) printf("  - Real part max error (%.4f) is too large.\n", max_real_error);
             if (max_abs_imag > FLOAT_COMPARISON_TOLERANCE) printf("  - Imaginary part max error (%.4f) is too large (should be near zero).\n", max_abs_imag);
            ret = 1;
            printf("  --> Verify the overall scaling logic (shift = log2N - total_shifts) and float comparison tolerance.\n");
            printf("  --> Check if INPUT_PRESCALE_FACTOR is suitable for this waveform (check Q31 saturation).\n");
        }
        // --- 比較代碼結束 ---

    } // end if (IFFT success)

cleanup:
    // --- Step 8: 清理 ---
    printf("\n[Step 8] Cleaning up memory...\n");
    if(input_signal_float) free(input_signal_float);
    if(q31_input_fft) free(q31_input_fft);     // 釋放 Q31 buffer
    if(q31_output_fft) free(q31_output_fft);   // 釋放 Q31 buffer
    if(q31_output_ifft) free(q31_output_ifft); // 釋放 Q31 buffer
    if(ifft_output_float_scaled) free(ifft_output_float_scaled);

    printf("\nDemonstration finished with exit code %d.\n", ret);
    return ret;
}