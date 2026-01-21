#!/bin/bash

# 简单的MASCI批量处理脚本
# 直接调用MASCI处理RAW目录下的所有.raw文件

set -e

# 配置参数
# 自动检测WSL环境并转换Windows路径
if [ -z "$MASCI_PATH" ]; then
    # 默认Windows路径
    DEFAULT_WIN_PATH="C:/Users/MASIC/MASIC_Console.exe"
    
    # 检测是否在WSL中
    if [ -d "/mnt/c" ] || [ -d "/mnt/host/c" ]; then
        # 尝试WSL路径
        if [ -f "/mnt/host/c/Users/MASIC/MASIC_Console.exe" ]; then
            MASCI_EXE="/mnt/host/c/Users/MASIC/MASIC_Console.exe"
        elif [ -f "/mnt/c/Users/MASIC/MASIC_Console.exe" ]; then
            MASCI_EXE="/mnt/c/Users/MASIC/MASIC_Console.exe"
        else
            MASCI_EXE="$DEFAULT_WIN_PATH"
        fi
    else
        MASCI_EXE="$DEFAULT_WIN_PATH"
    fi
else
    MASCI_EXE="$MASCI_PATH"
    # 如果提供了Windows路径，尝试转换为WSL路径
    if [[ "$MASCI_EXE" == C:* ]] || [[ "$MASCI_EXE" == c:* ]]; then
        # 转换 C:/path 或 C:\path 为 /mnt/c/path
        WSL_PATH=$(echo "$MASCI_EXE" | sed 's|^[Cc]:|/mnt/c|' | sed 's|\\|/|g')
        if [ -f "$WSL_PATH" ]; then
            MASCI_EXE="$WSL_PATH"
        elif [ -f "/mnt/host${WSL_PATH#/mnt}" ]; then
            MASCI_EXE="/mnt/host${WSL_PATH#/mnt}"
        fi
    fi
fi

# 支持命令行参数：第一个参数为输入目录，第二个为输出目录，第三个为TMT类型
if [ -n "$1" ]; then
    INPUT_DIR="$1"
elif [ -z "$INPUT_DIR" ]; then
    INPUT_DIR="./RAW"
fi

if [ -n "$2" ]; then
    OUTPUT_DIR="$2"
elif [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="./OutputData/QuantificationResults"
fi

if [ -n "$3" ]; then
    TMT_TYPE="$3"
elif [ -z "$TMT_TYPE" ]; then
    TMT_TYPE="TMT10"
fi

# 转换Windows路径到WSL路径（如果输入的是Windows路径）
convert_windows_path() {
    local path="$1"
    # 如果是Windows路径格式 (C:\, E:\, 等)
    if [[ "$path" =~ ^[A-Za-z]: ]]; then
        # 提取盘符（转换为小写）
        drive_letter=$(echo "$path" | cut -c1 | tr '[:upper:]' '[:lower:]')
        # 移除盘符和冒号，保留路径部分
        path_part=$(echo "$path" | sed 's/^[A-Za-z]://')
        # 转换反斜杠为正斜杠
        path_part=$(echo "$path_part" | sed 's/\\/\//g')
        # 移除开头的斜杠（如果有）
        path_part=$(echo "$path_part" | sed 's|^/||')
        # 构建WSL路径，优先尝试/mnt/host（WSL2新格式）
        if [ -d "/mnt/host/$drive_letter" ] 2>/dev/null; then
            wsl_path="/mnt/host/$drive_letter/$path_part"
        else
            wsl_path="/mnt/$drive_letter/$path_part"
        fi
        echo "$wsl_path"
    else
        echo "$path"
    fi
}

# 转换输入和输出路径
if [[ "$INPUT_DIR" =~ ^[A-Za-z]: ]]; then
    INPUT_DIR=$(convert_windows_path "$INPUT_DIR")
fi
if [[ "$OUTPUT_DIR" =~ ^[A-Za-z]: ]]; then
    OUTPUT_DIR=$(convert_windows_path "$OUTPUT_DIR")
fi

# 根据TMT类型选择参数文件
if [ "$TMT_TYPE" = "TMT10" ]; then
    PARAM_FILE="./MASCI_param/TMT10_LTQ-FT_10ppm_ReporterTol0.003Da_2014-08-06.xml"
elif [ "$TMT_TYPE" = "TMT11" ]; then
    PARAM_FILE="./MASCI_param/TMT11_LTQ-FT_10ppm_ReporterTol0.003Da_2017-03-17.xml"
else
    echo "Error: TMT_TYPE must be TMT10 or TMT11"
    exit 1
fi

# 检查MASCI是否存在
if [ ! -f "$MASCI_EXE" ]; then
    echo "Error: MASCI not found at $MASCI_EXE"
    echo ""
    echo "Please set MASCI_PATH environment variable:"
    echo "  In WSL: export MASCI_PATH='/mnt/c/Users/MASIC/MASIC_Console.exe'"
    echo "  Or:     export MASCI_PATH='/mnt/host/c/Users/MASIC/MASIC_Console.exe'"
    echo "  In Windows: set MASCI_PATH=C:\\Users\\MASIC\\MASIC_Console.exe"
    exit 1
fi

# 检查参数文件是否存在
if [ ! -f "$PARAM_FILE" ]; then
    echo "Error: Parameter file not found: $PARAM_FILE"
    exit 1
fi

# 检查输入目录
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory not found: $INPUT_DIR"
    echo ""
    echo "Please create the input directory and place your .raw files there:"
    echo "  mkdir -p $INPUT_DIR"
    echo "  # Then copy your .raw files to $INPUT_DIR/"
    echo ""
    echo "Or set a custom input directory:"
    echo "  export INPUT_DIR=\"./your_input_directory\""
    exit 1
fi

# 创建输出目录
mkdir -p "$OUTPUT_DIR"

# 处理所有.raw文件
echo "Starting MASCI processing..."
echo "MASCI: $MASCI_EXE"
echo "Input: $INPUT_DIR"
echo "Output: $OUTPUT_DIR"
echo "TMT Type: $TMT_TYPE"
echo ""

count=0
for raw_file in "$INPUT_DIR"/*.raw; do
    # 检查文件是否存在（如果没有.raw文件，通配符会返回字面量）
    if [ ! -f "$raw_file" ]; then
        echo "No .raw files found in $INPUT_DIR"
        exit 1
    fi
    
    count=$((count + 1))
    filename=$(basename "$raw_file")
    echo "[$count] Processing: $filename"
    
    # 调用MASCI
    "$MASCI_EXE" /P:"$PARAM_FILE" /I:"$raw_file" /O:"$OUTPUT_DIR"
    
    if [ $? -eq 0 ]; then
        echo "  ✓ Completed: $filename"
    else
        echo "  ✗ Failed: $filename"
    fi
    echo ""
done

echo "Processing complete! Processed $count file(s)."
echo "Results saved to: $OUTPUT_DIR"

