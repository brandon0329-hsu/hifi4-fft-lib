#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <complex.h>

#include "NatureDSP_types.h"
#include "NatureDSP_Signal_fft.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- 基本配置 ---
#define FFT_N 128
#define Q15_SCALE 32767.0f
#define SAMPLE_RATE 8000.0f

// --- 波形選擇 ---
// 1=雙正弦, 2=單餘弦, 3=方波, 4=衝激
// 5=AM信號, 6=窄帶FM近似, 7=線性掃頻(Chirp), 8=高斯脈衝, 9=白噪聲
#define WAVEFORM_TYPE 8 // <--- 修改這裡選擇進階波形

// --- 輸入預縮放 ---
// !!! 進階信號可能需要不同的預縮放因子 !!!
// AM/FM/Chirp/Gaussian: 0.05 可能是一個起點
// 白噪聲: 可能需要更小的因子，例如 0.01 或更低，取決於隨機數範圍
#define INPUT_PRESCALE_FACTOR 0.05f

// --- 比較容差 (浮點單位) ---
// 對於更複雜或含噪聲的信號，可能需要適當放寬容差
#define FLOAT_COMPARISON_TOLERANCE 0.02f // 例如放寬到 2%

// --- 確保包含之前的實現 ---
int16_t float_to_q15(float val) {
    float scaled_val = val * Q15_SCALE;
    if (scaled_val >= 32767.0f) return 32767;
    if (scaled_val <= -32768.0f) return -32768;
    return (int16_t)roundf(scaled_val);
}
float q15_to_float(int16_t val) {
    return (float)val / Q15_SCALE;
}
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
// --------------------------


int main() {
    printf("--- NatureDSP Fixed-point FFT -> IFFT Demonstration (N=%d) ---\n", FFT_N);
    printf("--- Verifying in FLOAT domain using scalingOption=3, Non-Inplace ---\n");

    // --- 1. 初始化與準備 (與之前相同) ---
    float   *input_signal_float = NULL;
    int16_t *q15_input_fft = NULL;
    int16_t *q15_output_fft = NULL;
    int16_t *q15_output_ifft = NULL;
    float complex *ifft_output_float_scaled = NULL;
    fft_handle_t handle_fft = NULL;
    fft_handle_t handle_ifft = NULL;
    int ret = 0;
    int total_shift_fft = 0;
    int total_shift_ifft = 0;
    int scalingOption = 3;

    printf("\n[Step 1] Allocating Memory (Separate Buffers)...\n");
    input_signal_float = (float *)malloc(FFT_N * sizeof(float));
    q15_input_fft = (int16_t *)malloc(2 * FFT_N * sizeof(int16_t));
    q15_output_fft = (int16_t *)malloc(2 * FFT_N * sizeof(int16_t));
    q15_output_ifft = (int16_t *)malloc(2 * FFT_N * sizeof(int16_t));
    ifft_output_float_scaled = (float complex *)malloc(FFT_N * sizeof(float complex));
    if (!input_signal_float || !q15_input_fft || !q15_output_fft || !q15_output_ifft || !ifft_output_float_scaled) {
        fprintf(stderr, "ERROR: Memory allocation failed!\n"); ret = 1; goto cleanup; }
    printf("Memory allocated successfully.\n");

    // --- Step 2: 生成選擇的浮點信號並預縮放 ---
#if WAVEFORM_TYPE == 1 // 雙正弦波
    printf("\n[Step 2] Generating Float Input Signal (Dual Sine) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float FREQ1 = 500.0f; const float FREQ2 = 1500.0f;
    for (unsigned int i = 0; i < FFT_N; i++) { double time = (double)i / SAMPLE_RATE; float raw_signal = (float)(0.7 * sin(2.0 * M_PI * FREQ1 * time) + 0.3 * sin(2.0 * M_PI * FREQ2 * time)); input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR; }
#elif WAVEFORM_TYPE == 2 // 單一餘弦波
    printf("\n[Step 2] Generating Float Input Signal (Single Cosine) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float COS_FREQ = 500.0f; const float COS_AMP = 0.8f;
    for (unsigned int i = 0; i < FFT_N; i++) { double time = (double)i / SAMPLE_RATE; float raw_signal = COS_AMP * cosf(2.0f * M_PI * COS_FREQ * time); input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR; }
#elif WAVEFORM_TYPE == 3 // 方波
    printf("\n[Step 2] Generating Float Input Signal (Square Wave) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float SQ_FREQ = 500.0f; const float SQ_AMP = 0.5f;
    for (unsigned int i = 0; i < FFT_N; i++) { double time = (double)i / SAMPLE_RATE; float phase = 2.0f * M_PI * SQ_FREQ * time; float raw_signal = (fmod(phase, 2.0 * M_PI) < M_PI) ? SQ_AMP : -SQ_AMP; input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR; }
#elif WAVEFORM_TYPE == 4 // 衝激函數
    printf("\n[Step 2] Generating Float Input Signal (Impulse) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const int IMPULSE_POS = 32; const float IMPULSE_AMP_RAW = 10.0f;
    memset(input_signal_float, 0, FFT_N * sizeof(float)); if (IMPULSE_POS >= 0 && IMPULSE_POS < FFT_N) { input_signal_float[IMPULSE_POS] = IMPULSE_AMP_RAW * INPUT_PRESCALE_FACTOR; } printf("Impulse at index %d with prescaled value %.5f\n", IMPULSE_POS, input_signal_float[IMPULSE_POS]);
#elif WAVEFORM_TYPE == 5 // AM 信號
    printf("\n[Step 2] Generating Float Input Signal (AM) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float AM_FC = 1000.0f; // 載波頻率
    const float AM_FM = 200.0f;  // 調製頻率
    const float AM_M = 0.8f;     // 調製指數
    const float AM_AMP = 0.6f;   // 載波幅度
    for (unsigned int i = 0; i < FFT_N; i++) {
        double time = (double)i / SAMPLE_RATE;
        float mod_signal = 1.0f + AM_M * cosf(2.0f * M_PI * AM_FM * time);
        float carrier_signal = AM_AMP * cosf(2.0f * M_PI * AM_FC * time);
        float raw_signal = mod_signal * carrier_signal;
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#elif WAVEFORM_TYPE == 6 // 窄帶 FM 近似
    printf("\n[Step 2] Generating Float Input Signal (Narrowband FM Approx) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float FM_FC = 1000.0f; // 載波頻率
    const float FM_FM = 200.0f;  // 調製頻率
    const float FM_DELTA_F = 50.0f; // 最大頻偏 ( << FM_FC )
    const float FM_BETA = FM_DELTA_F / FM_FM; // 調製指數 (應較小)
    const float FM_AMP = 0.7f;   // 幅度
    printf("FM Beta: %.3f\n", FM_BETA);
    for (unsigned int i = 0; i < FFT_N; i++) {
        double time = (double)i / SAMPLE_RATE;
        float term1 = FM_AMP * cosf(2.0f * M_PI * FM_FC * time);
        float term2 = -(FM_AMP * FM_BETA / 2.0f) * cosf(2.0f * M_PI * (FM_FC - FM_FM) * time);
        float term3 = (FM_AMP * FM_BETA / 2.0f) * cosf(2.0f * M_PI * (FM_FC + FM_FM) * time);
        float raw_signal = term1 + term2 + term3;
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#elif WAVEFORM_TYPE == 7 // 線性掃頻 (Chirp)
    printf("\n[Step 2] Generating Float Input Signal (Linear Chirp) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float CHIRP_F0 = 100.0f;    // 起始頻率
    const float CHIRP_F1 = 1000.0f;   // 結束頻率
    const float CHIRP_AMP = 0.7f;    // 幅度
    const double T = (double)FFT_N / SAMPLE_RATE; // 信號總時長
    const double k = (CHIRP_F1 - CHIRP_F0) / T;   // 掃頻速率
    for (unsigned int i = 0; i < FFT_N; i++) {
        double time = (double)i / SAMPLE_RATE;
        float phase = 2.0f * M_PI * (CHIRP_F0 * time + k * time * time / 2.0);
        float raw_signal = CHIRP_AMP * cosf(phase);
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#elif WAVEFORM_TYPE == 8 // 高斯脈衝
    printf("\n[Step 2] Generating Float Input Signal (Gaussian Pulse) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float GAUSS_AMP = 1.0f; // 幅度
    const double GAUSS_T0 = (double)(FFT_N / 2) / SAMPLE_RATE; // 中心時間在中間
    const double GAUSS_SIGMA = (double)(FFT_N / 16) / SAMPLE_RATE; // 控制寬度
    for (unsigned int i = 0; i < FFT_N; i++) {
        double time = (double)i / SAMPLE_RATE;
        double exponent = -pow((time - GAUSS_T0) / (2.0 * GAUSS_SIGMA), 2.0);
        float raw_signal = GAUSS_AMP * exp(exponent);
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#elif WAVEFORM_TYPE == 9 // 白噪聲
    printf("\n[Step 2] Generating Float Input Signal (White Noise) & Prescaling (Factor: %.3f)...\n", INPUT_PRESCALE_FACTOR);
    const float NOISE_AMP = 0.8f; // 控制噪聲幅度
    // 初始化隨機數種子 (可選, 為了可重複性)
    srand(0);
    for (unsigned int i = 0; i < FFT_N; i++) {
        // 生成 [-1, 1] 範圍內的偽隨機浮點數
        float rand_val = ((float)rand() / (float)RAND_MAX) * 2.0f - 1.0f;
        float raw_signal = NOISE_AMP * rand_val;
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
#else
    #error "Invalid WAVEFORM_TYPE selected!"
#endif
    printf("Input signal generated and prescaled.\n");

    // --- Step 3: 轉換為 Q15 FFT 輸入 (與之前相同) ---
    printf("\n[Step 3] Converting Prescaled Float to Q15 FFT Input...\n");
    int16_t min_q15 = 32767, max_q15 = -32768; int saturated_count = 0;
    for (unsigned int i = 0; i < FFT_N; i++) {
        int16_t real_q15 = float_to_q15(input_signal_float[i]);
        q15_input_fft[2 * i + 0] = real_q15; q15_input_fft[2 * i + 1] = 0;
        if (input_signal_float[i] >= 1.0f || input_signal_float[i] < -1.0f ) { if (real_q15 == 32767 || real_q15 == -32768){ saturated_count++; } }
        if (real_q15 < min_q15) min_q15 = real_q15; if (real_q15 > max_q15) max_q15 = real_q15;
    }
    printf("Q15 conversion complete. Q15 Input Range: [%d, %d]. Saturated samples: %d\n", min_q15, max_q15, saturated_count);
     if (saturated_count > 0) printf("WARNING: Saturation occurred during Q15 conversion. Consider adjusting INPUT_PRESCALE_FACTOR.\n");

    // --- Step 4: 執行 FFT (非原地) (不變) ---
    printf("\n[Step 4] Performing Fixed-point FFT (Non-Inplace)...\n");
    handle_fft = cfft16_128; if (!handle_fft) { /*...*/ }
    total_shift_fft = fft_cplx16x16(q15_output_fft, q15_input_fft, handle_fft, scalingOption);
    if (total_shift_fft < 0) { /*...*/ }
    printf("FFT finished. Total scaling shifts (FFT): %d\n", total_shift_fft);

    // --- Step 5: 執行 IFFT (非原地) (不變) ---
    printf("\n[Step 5] Performing Fixed-point IFFT (Non-Inplace)...\n");
    handle_ifft = cifft16_128; if (!handle_ifft) { /*...*/ }
    total_shift_ifft = ifft_cplx16x16(q15_output_ifft, q15_output_fft, handle_ifft, scalingOption);
    if (total_shift_ifft < 0) { /*...*/ ret = 1; }
    else {
        printf("IFFT finished. Total scaling shifts (IFFT): %d\n", total_shift_ifft);

        // --- Step 6: 轉換 IFFT Q15 輸出到浮點並應用整體縮放 (不變) ---
        printf("\n[Step 6] Converting IFFT Q15 Output to Float and Applying *Overall* Scaling...\n");
        int total_shifts = total_shift_fft + total_shift_ifft;
        int log2N = 0; if (FFT_N > 0) { int tempN = FFT_N; while (tempN > 1) { tempN >>= 1; log2N++; } }
        int shift_for_vecToFp64 = log2N - total_shifts;
        printf("Calculated equivalent total shift for vecToFp64: %d (log2N - total_shifts)\n", shift_for_vecToFp64);
        int16_t* ifft_q15_result_ptr = q15_output_ifft;
        for (unsigned int i = 0; i < FFT_N; i++) {
            double real_fp64 = ldexp((double)ifft_q15_result_ptr[2 * i + 0], -shift_for_vecToFp64 - 15);
            double imag_fp64 = ldexp((double)ifft_q15_result_ptr[2 * i + 1], -shift_for_vecToFp64 - 15);
            ifft_output_float_scaled[i] = (float)real_fp64 + (float)imag_fp64 * I;
        }
        printf("Q15 to Float conversion and overall scaling complete.\n");

        // --- Step 7: 在浮點域比較 (不變) ---
        printf("\n[Step 7] Comparing Scaled IFFT Float Output with Original Prescaled Float Input...\n");
        double rmse_real = calculate_float_rmse_real(ifft_output_float_scaled, input_signal_float, FFT_N);
        float max_abs_imag = calculate_float_max_abs_imag(ifft_output_float_scaled, FFT_N);
        float max_real_error = 0.0f;
        for(unsigned int i = 0; i < FFT_N; ++i) { float current_real_error = fabsf(crealf(ifft_output_float_scaled[i]) - input_signal_float[i]); if (current_real_error > max_real_error) max_real_error = current_real_error; }
        printf("\nOverall Comparison Metrics (Float units):\n");
        printf("  Real Part RMSE vs Expected: %e\n", rmse_real);
        printf("  Max Abs Real Part Error:    %e\n", max_real_error);
        printf("  Max Abs Imaginary Part:     %e\n", max_abs_imag);

        // --- 打印浮點比較表格 (可以保持打印部分樣本，或根據需要修改) ---
        printf("\nFloat Comparison Table (first/last 8 samples):\n"); // 可以修改標題
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
        // --------------------------------------

        // --- 判斷成功/失敗 (不變) ---
        if (max_real_error <= FLOAT_COMPARISON_TOLERANCE && max_abs_imag <= FLOAT_COMPARISON_TOLERANCE) {
            printf("\nSUCCESS: Max errors in float domain are within tolerance (%.4f).\n", FLOAT_COMPARISON_TOLERANCE);
            printf("This demonstrates IFFT(FFT(x)) recovers the signal in float domain when overall scaling (log2N - total_shifts) is compensated.\n");
        } else {
            printf("\nFAILURE: Errors in float domain exceed tolerance (%.4f).\n", FLOAT_COMPARISON_TOLERANCE);
             if (max_real_error > FLOAT_COMPARISON_TOLERANCE) printf("  - Real part max error (%.4f) is too large.\n", max_real_error);
             if (max_abs_imag > FLOAT_COMPARISON_TOLERANCE) printf("  - Imaginary part max error (%.4f) is too large (should be near zero).\n", max_abs_imag);
            ret = 1;
            printf("  --> Verify the overall scaling logic (shift = log2N - total_shifts) and float comparison tolerance.\n");
            printf("  --> Check if INPUT_PRESCALE_FACTOR is suitable for this waveform (check Q15 saturation).\n");
        }

    } // end if (IFFT success)

cleanup:
    // --- Step 8: 清理 (不變) ---
    printf("\n[Step 8] Cleaning up memory...\n");
    if(input_signal_float) free(input_signal_float);
    if(q15_input_fft) free(q15_input_fft);
    if(q15_output_fft) free(q15_output_fft);
    if(q15_output_ifft) free(q15_output_ifft);
    if(ifft_output_float_scaled) free(ifft_output_float_scaled);

    printf("\nDemonstration finished with exit code %d.\n", ret);
    return ret;
}