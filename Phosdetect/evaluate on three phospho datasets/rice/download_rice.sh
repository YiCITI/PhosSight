#!/bin/bash
mkdir -p /data0/wangb/cd/duibi0826/0826comparison/data/rice
cd /data0/wangb/cd/duibi0826/0826comparison/data/rice
echo "下载 Rice GSB 数据..."
curl -L -o Rice_phospho_GSB.csv "https://peptideatlas.org/builds/rice/Rice_phospho_GSB_all_proteins.csv"
echo "下载 Rice FASTA 文件..."
curl -L -o Rice.fasta "https://peptideatlas.org/builds/rice/202204/Rice.fasta"
echo "完成"
ls -la
