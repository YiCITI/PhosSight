#!/bin/bash
###############################################################################
# PhosDetect 肽段可检测性预测脚本
# 用途: 从FASTA文件读取肽段序列，使用PhosSight模型预测detectability scores
# 
# 使用方法:
#   ./predict_detectability.sh [输入FASTA文件] [输出TXT文件]
#   
# 示例:
#   ./predict_detectability.sh peptides.fasta results.txt
#   ./predict_detectability.sh /path/to/peptides.fasta /path/to/results.txt
###############################################################################

# 设置颜色输出
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# ============================================================================
# 配置参数 - 请根据实际情况修改这些参数
# ============================================================================

# PhosSight 模型路径
PHOSSIGHT_MODEL_PATH="/data0/wangb/cd/duibi0826/0826comparison/model_weights/phosight_v2_best.pth"

# PhosSight 模型代码路径
PHOSSIGHT_CODE_PATH="/data0/wangb/cd/PhosSight-main/PhosSight/sotaimprov2"

# Python 脚本路径
PYTHON_SCRIPT="/data0/wangb/wbscy/czy/0925/fasta_peptide_scoring.py"

# 默认输入输出文件（如果未提供参数）
DEFAULT_INPUT_FASTA="/data0/wangb/wbscy/czy/0925/human_all_peptides_2_7_46.fasta"

# ============================================================================
# 函数定义
# ============================================================================

print_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_header() {
    echo -e "${BLUE}$1${NC}"
}

show_usage() {
    echo "用法: $0 [输入FASTA文件] [输出TXT文件]"
    echo ""
    echo "参数说明:"
    echo "  输入FASTA文件: 包含肽段序列的FASTA格式文件"
    echo "                  (默认: $DEFAULT_INPUT_FASTA)"
    echo "  输出TXT文件:   预测结果的输出文件"
    echo "                 (默认: <输入文件名>_detectability_scores.txt)"
    echo ""
    echo "示例:"
    echo "  $0 peptides.fasta results.txt"
    echo "  $0 /path/to/peptides.fasta"
    echo ""
    echo "注意: 请确保在脚本中配置正确的模型路径"
}

# ============================================================================
# 主程序
# ============================================================================

main() {
    print_header "=========================================="
    print_header "PhosDetect 肽段可检测性预测"
    print_header "=========================================="
    echo ""
    
    # 显示帮助信息
    if [ "$1" = "-h" ] || [ "$1" = "--help" ]; then
        show_usage
        exit 0
    fi
    
    # 获取输入文件
    if [ -z "$1" ]; then
        INPUT_FASTA_FILE="$DEFAULT_INPUT_FASTA"
        print_warning "未指定输入文件，使用默认: $INPUT_FASTA_FILE"
    else
        INPUT_FASTA_FILE="$1"
    fi
    
    # 获取输出文件
    if [ -z "$2" ]; then
        # 从输入文件名生成输出文件名
        OUTPUT_TXT_FILE="${INPUT_FASTA_FILE%.fasta}_detectability_scores.txt"
        OUTPUT_TXT_FILE="${OUTPUT_TXT_FILE%.fa}_detectability_scores.txt"
    else
        OUTPUT_TXT_FILE="$2"
    fi
    
    # 检查输入文件
    print_info "检查输入文件: $INPUT_FASTA_FILE"
    if [ ! -f "$INPUT_FASTA_FILE" ]; then
        print_error "输入文件不存在: $INPUT_FASTA_FILE"
        exit 1
    fi
    
    # 检查模型文件
    print_info "检查模型文件: $PHOSSIGHT_MODEL_PATH"
    if [ ! -f "$PHOSSIGHT_MODEL_PATH" ]; then
        print_error "模型文件不存在: $PHOSSIGHT_MODEL_PATH"
        print_error "请修改脚本中的 PHOSSIGHT_MODEL_PATH 变量"
        exit 1
    fi
    
    # 检查模型代码路径
    print_info "检查模型代码路径: $PHOSSIGHT_CODE_PATH"
    if [ ! -d "$PHOSSIGHT_CODE_PATH" ]; then
        print_error "模型代码路径不存在: $PHOSSIGHT_CODE_PATH"
        print_error "请修改脚本中的 PHOSSIGHT_CODE_PATH 变量"
        exit 1
    fi
    
    # 检查Python脚本
    if [ ! -f "$PYTHON_SCRIPT" ]; then
        print_error "Python脚本不存在: $PYTHON_SCRIPT"
        exit 1
    fi
    
    # 检查Python环境
    if ! command -v python3 &> /dev/null; then
        print_error "未找到python3，请先安装Python"
        exit 1
    fi
    
    # 检查PyTorch
    python3 -c "import torch" 2>/dev/null
    if [ $? -ne 0 ]; then
        print_error "PyTorch未安装，请先安装PyTorch"
        exit 1
    fi
    
    # 检查CUDA可用性
    CUDA_AVAILABLE=$(python3 -c "import torch; print('cuda' if torch.cuda.is_available() else 'cpu')" 2>/dev/null)
    print_info "使用设备: $CUDA_AVAILABLE"
    
    # 创建输出目录（如果不存在）
    OUTPUT_DIR=$(dirname "$OUTPUT_TXT_FILE")
    if [ ! -d "$OUTPUT_DIR" ] && [ "$OUTPUT_DIR" != "." ]; then
        print_info "创建输出目录: $OUTPUT_DIR"
        mkdir -p "$OUTPUT_DIR"
    fi
    
    # 显示配置信息
    echo ""
    print_header "配置信息:"
    echo "  模型路径: $PHOSSIGHT_MODEL_PATH"
    echo "  模型代码: $PHOSSIGHT_CODE_PATH"
    echo "  输入文件: $INPUT_FASTA_FILE"
    echo "  输出文件: $OUTPUT_TXT_FILE"
    echo ""
    
    # 创建临时Python脚本，修改路径参数
    print_info "准备执行预测..."
    TEMP_SCRIPT=$(mktemp /tmp/phosdetect_predict_XXXXXX.py)
    
    # 复制原始脚本并修改路径
    sed "s|model_path = \".*\"|model_path = \"$PHOSSIGHT_MODEL_PATH\"|" "$PYTHON_SCRIPT" | \
    sed "s|fasta_file = \".*\"|fasta_file = \"$INPUT_FASTA_FILE\"|" | \
    sed "s|output_file = \".*\"|output_file = \"$OUTPUT_TXT_FILE\"|" | \
    sed "s|sys.path.append('.*')|sys.path.append('$PHOSSIGHT_CODE_PATH')|" > "$TEMP_SCRIPT"
    
    # 执行预测
    print_header "=========================================="
    print_info "开始预测..."
    print_header "=========================================="
    echo ""
    
    python3 "$TEMP_SCRIPT"
    EXIT_CODE=$?
    
    # 清理临时文件
    rm -f "$TEMP_SCRIPT"
    
    # 检查执行结果
    if [ $EXIT_CODE -eq 0 ] && [ -f "$OUTPUT_TXT_FILE" ]; then
        echo ""
        print_header "=========================================="
        print_info "预测完成！"
        print_info "结果已保存到: $OUTPUT_TXT_FILE"
        
        # 显示结果统计
        TOTAL_LINES=$(wc -l < "$OUTPUT_TXT_FILE" 2>/dev/null || echo "0")
        if [ "$TOTAL_LINES" -gt 1 ]; then
            TOTAL_PEPTIDES=$((TOTAL_LINES - 1))  # 减去标题行
            print_info "预测的肽段数量: $TOTAL_PEPTIDES"
            
            # 显示前几行结果
            echo ""
            print_info "结果预览（前5行）:"
            head -6 "$OUTPUT_TXT_FILE" | column -t
        fi
    else
        print_error "预测失败，请检查错误信息"
        exit 1
    fi
    
    echo ""
    print_info "处理完成！"
}

# 运行主程序
main "$@"