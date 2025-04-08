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

#define FFT_N 128
#define Q15_SCALE 32767.0f
#define SAMPLE_RATE 8000.0f
#define FREQ1 500.0f
#define FREQ2 1500.0f
#define INPUT_PRESCALE_FACTOR 0.05f
#define FLOAT_COMPARISON_TOLERANCE 0.01f // 浮點比較容差

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
    double mse = 0.0;
    if (n == 0) return 0.0;
    for (size_t i = 0; i < n; ++i) {
        double diff = crealf(arr_complex[i]) - arr_real[i];
        mse += diff * diff;
    }
    return sqrt(mse / n);
}
float calculate_float_max_abs_imag(const float complex* arr_complex, size_t n) {
    float max_abs_imag = 0.0f;
     for (size_t i = 0; i < n; ++i) {
         float abs_imag = fabsf(cimagf(arr_complex[i]));
        if (abs_imag > max_abs_imag) max_abs_imag = abs_imag;
     }
    return max_abs_imag;
}


int main() {
    printf("--- NatureDSP Fixed-point FFT -> IFFT Demonstration (N=%d) ---\n", FFT_N);
    printf("--- Verifying in FLOAT domain using scalingOption=3 --- \n");
    printf("--- *** USING SEPARATE BUFFERS (NON-INPLACE) *** ---\n"); // 強調使用獨立緩衝區

    // --- 1. 初始化與準備 ---
    float   *input_signal_float = NULL;      // 原始浮點時域信號 (預縮放後)
    int16_t *q15_input_fft = NULL;        // FFT 的 Q15 輸入緩衝區
    int16_t *q15_output_fft = NULL;       // FFT 的 Q15 輸出緩衝區 (也作為 IFFT 的輸入)
    int16_t *q15_output_ifft = NULL;      // IFFT 的 Q15 輸出緩衝區
    float complex *ifft_output_float_scaled = NULL; // 存放最終縮放後的浮點結果
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
    q15_output_ifft = (int16_t *)malloc(2 * FFT_N * sizeof(int16_t)); // 為 IFFT 輸出分配獨立空間
    ifft_output_float_scaled = (float complex *)malloc(FFT_N * sizeof(float complex));

    if (!input_signal_float || !q15_input_fft || !q15_output_fft || !q15_output_ifft || !ifft_output_float_scaled) {
        fprintf(stderr, "ERROR: Memory allocation failed!\n");
        ret = 1;
        goto cleanup;
    }
    printf("Memory allocated successfully.\n");

    // --- Step 2: 生成浮點信號 (不變) ---
    printf("\n[Step 2] Generating Float Input Signal & Prescaling...\n");
    for (unsigned int i = 0; i < FFT_N; i++) {
        double time = (double)i / SAMPLE_RATE;
        float raw_signal = (float)(0.7 * sin(2.0 * M_PI * FREQ1 * time) +
                                   0.3 * sin(2.0 * M_PI * FREQ2 * time));
        input_signal_float[i] = raw_signal * INPUT_PRESCALE_FACTOR;
    }
    printf("Input signal generated and prescaled.\n");

    // --- Step 3: 轉換為 Q15 FFT 輸入 (不變) ---
    printf("\n[Step 3] Converting Prescaled Float to Q15 FFT Input...\n");
     int16_t min_q15 = 32767, max_q15 = -32768; int saturated_count = 0; // 省略範圍細節打印
    for (unsigned int i = 0; i < FFT_N; i++) {
        int16_t real_q15 = float_to_q15(input_signal_float[i]);
        q15_input_fft[2 * i + 0] = real_q15;
        q15_input_fft[2 * i + 1] = 0;
        // Saturation check omitted for brevity
         if (real_q15 < min_q15) min_q15 = real_q15; if (real_q15 > max_q15) max_q15 = real_q15;
    }
    printf("Q15 conversion for FFT input complete.\n");

    // --- 2. 執行 FFT (非原地) ---
    printf("\n[Step 4] Performing Fixed-point FFT (Non-Inplace)...\n");
    handle_fft = cfft16_128; // (***確認名稱***)
    if (!handle_fft) { /* Error handling */ ret = 1; goto cleanup; }

    total_shift_fft = fft_cplx16x16(q15_output_fft, q15_input_fft, handle_fft, scalingOption); // 輸出到 q15_output_fft

    if (total_shift_fft < 0) { /* Error handling */ ret = 1; goto cleanup; }
    printf("FFT finished. Total scaling shifts (FFT): %d\n", total_shift_fft);

    // --- 3. 執行 IFFT (非原地) ---
    printf("\n[Step 5] Performing Fixed-point IFFT (Non-Inplace)...\n");
    handle_ifft = cifft16_128; // (***確認名稱***)
     if (!handle_ifft) { /* Error handling */ ret = 1; goto cleanup; }

    // 輸入是 FFT 的輸出 q15_output_fft, 輸出到獨立的 q15_output_ifft
    total_shift_ifft = ifft_cplx16x16(q15_output_ifft, q15_output_fft, handle_ifft, scalingOption);

    if (total_shift_ifft < 0) { /* Error handling */ ret = 1; }
    else {
        printf("IFFT finished. Total scaling shifts (IFFT): %d\n", total_shift_ifft);

        // --- 4. 轉換 IFFT 輸出到浮點並應用 *整體* 縮放 (邏輯同上一版) ---
        printf("\n[Step 6] Converting IFFT Q15 Output to Float and Applying *Overall* Scaling...\n");

        int total_shifts = total_shift_fft + total_shift_ifft;
        int log2N = 0;
        if (FFT_N > 0) { int tempN = FFT_N; while (tempN > 1) { tempN >>= 1; log2N++; } }
        int shift_for_vecToFp64 = log2N - total_shifts;
        printf("Calculated equivalent total shift for vecToFp64: %d (log2N - total_shifts)\n", shift_for_vecToFp64);

        int16_t* ifft_q15_result_ptr = q15_output_ifft; // 使用獨立的 IFFT 輸出緩衝區

        for (unsigned int i = 0; i < FFT_N; i++) {
            double real_fp64 = ldexp((double)ifft_q15_result_ptr[2 * i + 0], -shift_for_vecToFp64 - 15);
            double imag_fp64 = ldexp((double)ifft_q15_result_ptr[2 * i + 1], -shift_for_vecToFp64 - 15);
            ifft_output_float_scaled[i] = (float)real_fp64 + (float)imag_fp64 * I;
        }
        printf("Q15 to Float conversion and overall scaling complete.\n");

        // --- 5. 在浮點域比較 (與上一版相同) ---
        printf("\n[Step 7] Comparing Scaled IFFT Float Output with Original Prescaled Float Input...\n");
        // ... (比較代碼、表格打印、成功/失敗判斷 - 與上一版完全相同) ...
        // --- 比較代碼開始 ---
        double rmse_real = calculate_float_rmse_real(ifft_output_float_scaled, input_signal_float, FFT_N);
        float max_abs_imag = calculate_float_max_abs_imag(ifft_output_float_scaled, FFT_N);
        float max_real_error = 0.0f;
        for(unsigned int i = 0; i < FFT_N; ++i) {
             float current_real_error = fabsf(crealf(ifft_output_float_scaled[i]) - input_signal_float[i]);
             if (current_real_error > max_real_error) max_real_error = current_real_error;
        }
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
            printf("This demonstrates IFFT(FFT(x)) recovers the signal in float domain when overall scaling (log2N - total_shifts) is compensated.\n");
        } else {
            printf("\nFAILURE: Errors in float domain exceed tolerance (%.4f).\n", FLOAT_COMPARISON_TOLERANCE);
             if (max_real_error > FLOAT_COMPARISON_TOLERANCE) printf("  - Real part max error (%.4f) is too large.\n", max_real_error);
             if (max_abs_imag > FLOAT_COMPARISON_TOLERANCE) printf("  - Imaginary part max error (%.4f) is too large (should be near zero).\n", max_abs_imag);
            ret = 1;
            printf("  --> Verify the overall scaling logic (shift = log2N - total_shifts) and float comparison tolerance.\n");
            printf("  --> Non-inplace operation used. If errors persist, issue might be scaling logic or library behavior.\n");
        }
        // --- 比較代碼結束 ---

    } // end if (IFFT success)

cleanup:
    // --- 6. 清理 ---
    printf("\n[Step 8] Cleaning up memory...\n");
    if(input_signal_float) free(input_signal_float);
    if(q15_input_fft) free(q15_input_fft);
    if(q15_output_fft) free(q15_output_fft);
    if(q15_output_ifft) free(q15_output_ifft); // 釋放獨立的 IFFT 輸出緩衝區
    if(ifft_output_float_scaled) free(ifft_output_float_scaled);

    printf("\nDemonstration finished with exit code %d.\n", ret);
    return ret;
}