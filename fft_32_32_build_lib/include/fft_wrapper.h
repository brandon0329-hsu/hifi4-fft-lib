#ifndef FFT_WRAPPER_H
#define FFT_WRAPPER_H

#include <stdint.h>
#include "NatureDSP_types.h"
#include "NatureDSP_Signal_fft.h"

// Constants
#define FFT_N 128
#define Q31_SCALE 2147483647.0f
#define SAMPLE_RATE 8000.0f
#define DEFAULT_PRESCALE_FACTOR 0.025f
#define DEFAULT_TOLERANCE 0.02f

// Signal types
typedef enum {
    WAVEFORM_DUAL_SINE = 1,
    WAVEFORM_SINGLE_COSINE = 2,
    WAVEFORM_SQUARE = 3,
    WAVEFORM_IMPULSE = 4,
    WAVEFORM_AM = 5,
    WAVEFORM_FM = 6,
    WAVEFORM_CHIRP = 7,
    WAVEFORM_GAUSSIAN = 8,
    WAVEFORM_WHITE_NOISE = 9
} signal_type_t;

// FFT configuration structure
typedef struct {
    float prescale_factor;
    int scaling_option;
    float tolerance;
} fft_config_t;

// Function prototypes
// Q31轉換函數
int32_t float_to_q31(float val);
float q31_to_float(int32_t val);

// FFT/IFFT核心函數
int perform_fft_32_32(const float* input, 
                     int32_t* output, 
                     int* shifts,
                     int scaling_option);

int perform_ifft_32_32(const int32_t* input, 
                      float* output, 
                      int fft_shifts, 
                      int* shifts,
                      int scaling_option);

// 信號生成函數
int generate_test_signal(float* buffer, 
                        signal_type_t type, 
                        float prescale);

// 工具函數
void print_signal_preview(const float* signal, int length, int preview_len);
double calculate_signal_error_metrics(const float* original, 
                                    const float* recovered, 
                                    int length,
                                    float* max_error,
                                    double* rmse);

// FFT結果結構體
typedef struct {
    int success;
    int fft_shifts;
    int ifft_shifts;
    double rmse;
    float max_error;
    float prescale_factor;
} fft_result_t;

// 靜態庫測試結果結構體
typedef struct {
    int test_passed;
    const char* test_name;
    const char* description;
    float error_value;
    const char* error_message;
} static_lib_test_result_t;

// 靜態庫測試函數
int test_static_library(static_lib_test_result_t* results, int* num_tests);

#endif // FFT_WRAPPER_H