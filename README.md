# FFT/IFFT Implementation for HiFi4

本專案實現了基於Xtensa HiFi4 DSP的FFT（快速傅立葉變換）和IFFT（逆快速傅立葉變換）功能。

## 功能特點

- 支援32x32位元FFT/IFFT運算
- 提供完整的靜態庫實現
- 支援9種不同波形的測試：
  1. 雙正弦波
  2. 單餘弦波
  3. 方波
  4. 衝激函數
  5. AM信號
  6. 窄帶FM信號
  7. 線性掃頻信號
  8. 高斯脈衝
  9. 白噪聲

## 目錄結構

```
fft_32_32_build_lib/
├── include/          # 頭文件
├── src/             # 源代碼
├── lib/             # 編譯後的庫文件
└── test/            # 測試程序
```

## 環境要求

- Xtensa工具鏈
- HiFi4 DSP庫
- Linux作業系統

## 環境設置

```bash
export XTENSA_SYSTEM=/home/yu-chen/xtensa/XtDevTools/install/tools/RI-2021.6-linux/XtensaTools/config
export XTENSA_BASE=/home/yu-chen/xtensa/XtDevTools/install
export XTENSA_CORE=hifi4_ss_spfpu_7
```

## 編譯和運行

1. 編譯靜態庫：
```bash
make all
```

2. 運行測試：
```bash
make test
```

## 測試報告

測試程序會自動生成包含以下內容的報告：
- 靜態庫功能測試結果
- 9種波形的FFT/IFFT轉換測試結果
- 精度分析和誤差統計

## 更新日誌

### 2025-04-08
- 添加靜態庫測試功能
- 實現9種波形的完整支持
- 優化Q31轉換精度