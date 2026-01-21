@echo off
REM 简单的MASCI批量处理脚本（Windows版本）
REM 直接调用MASCI处理RAW目录下的所有.raw文件

setlocal enabledelayedexpansion

REM 配置参数
if not defined MASCI_PATH (
    set MASCI_PATH=C:\Users\MASIC\MASIC_Console.exe
)

REM 支持命令行参数：第一个参数为输入目录，第二个为输出目录，第三个为TMT类型
if not "%~1"=="" (
    set INPUT_DIR=%~1
) else if not defined INPUT_DIR (
    set INPUT_DIR=.\RAW
)

if not "%~2"=="" (
    set OUTPUT_DIR=%~2
) else if not defined OUTPUT_DIR (
    set OUTPUT_DIR=.\OutputData\QuantificationResults
)

if not "%~3"=="" (
    set TMT_TYPE=%~3
) else if not defined TMT_TYPE (
    set TMT_TYPE=TMT10
)

REM 根据TMT类型选择参数文件
if "%TMT_TYPE%"=="TMT10" (
    set PARAM_FILE=.\MASCI_param\TMT10_LTQ-FT_10ppm_ReporterTol0.003Da_2014-08-06.xml
) else if "%TMT_TYPE%"=="TMT11" (
    set PARAM_FILE=.\MASCI_param\TMT11_LTQ-FT_10ppm_ReporterTol0.003Da_2017-03-17.xml
) else (
    echo Error: TMT_TYPE must be TMT10 or TMT11
    exit /b 1
)

REM 检查MASCI是否存在
if not exist "%MASCI_PATH%" (
    echo Error: MASCI not found at %MASCI_PATH%
    echo Please set MASCI_PATH environment variable
    exit /b 1
)

REM 检查参数文件是否存在
if not exist "%PARAM_FILE%" (
    echo Error: Parameter file not found: %PARAM_FILE%
    exit /b 1
)

REM 检查输入目录
if not exist "%INPUT_DIR%" (
    echo Error: Input directory not found: %INPUT_DIR%
    exit /b 1
)

REM 创建输出目录
if not exist "%OUTPUT_DIR%" mkdir "%OUTPUT_DIR%"

REM 处理所有.raw文件
echo Starting MASCI processing...
echo MASCI: %MASCI_PATH%
echo Input: %INPUT_DIR%
echo Output: %OUTPUT_DIR%
echo TMT Type: %TMT_TYPE%
echo.

set count=0
for %%f in ("%INPUT_DIR%\*.raw") do (
    set /a count+=1
    set "filename=%%~nxf"
    echo [!count!] Processing: !filename!
    
    REM 调用MASCI
    "%MASCI_PATH%" /P:"%PARAM_FILE%" /I:"%%f" /O:"%OUTPUT_DIR%"
    
    if !ERRORLEVEL! EQU 0 (
        echo   [OK] Completed: !filename!
    ) else (
        echo   [FAILED] Failed: !filename!
    )
    echo.
)

if %count% EQU 0 (
    echo No .raw files found in %INPUT_DIR%
    exit /b 1
)

echo Processing complete! Processed %count% file(s).
echo Results saved to: %OUTPUT_DIR%

endlocal

