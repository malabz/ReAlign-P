#ifndef REFINE_P_UTILS_H
#define REFINE_P_UTILS_H

#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include "fasta.h"

inline utils::fasta read_from(std::string file_path) {
    std::ifstream file(file_path);
    if (!file) {
        std::cerr << "Error: cannot open file " << file_path << std::endl;
        std::cerr << "Please check that the file path is correct and make sure the file exists." << std::endl;
        exit(1);
    }
    utils::fasta fasta(file);
    file.close();
    return fasta;
}

// 将比对拆分为左右部分和中间部分（去掉gap）
inline void locate_middle_part(size_t& left_idx, size_t& right_idx, const std::vector<std::string>& alignment, const std::vector<int>& col_status) {
    // 初始值为超出范围，表示未找到第一个1
    left_idx = 0;
    // 初始值为超出范围，表示未找到最后一个1
    right_idx = col_status.size() - 1;

    // 单次遍历查找第一个和最后一个1的索引
    for (size_t i = 0; i < col_status.size(); ++i) {
        if (col_status[i] == 1) {
            if (left_idx == col_status.size()) {
                left_idx = i;  // 记录第一个1的位置
            }
            right_idx = i;  // 每次找到1就更新最后一个1的位置
        }
    }

    // 如果需要，可以在此打印调试信息
    // std::cout << "First 1 found at index: " << left_idx << std::endl;
    // std::cout << "Last 1 found at index: " << right_idx << std::endl;
}

inline std::vector<int> column_status(const std::vector<std::string>& alignment) {
    size_t num_cols = alignment[0].size();  // 获取列数
    std::vector<int> col_status(num_cols, 0);  // 初始化每列的状态为 0

    for (size_t col = 0; col < num_cols; ++col) {
        bool has_gap = false;
        for (const auto& row : alignment) {
            if (row[col] == '-') {
                has_gap = true;
                break;  // 如果这一列有 gap，跳出循环
            }
        }
        if (!has_gap) {
            col_status[col] = 1;  // 如果这一列没有 gap，状态设为 1
        }
    }
    return col_status;  // 返回列状态数组
}

inline void split_alignment(size_t left_idx, size_t right_idx, std::vector<std::string>& alignment,
                     std::vector<std::string>& splited_sequences, std::vector<std::string>& splited_alignment, std::vector<std::string>& tail_alignment)
{
    std::vector<std::string> head_alignment;
    for (const auto& row : alignment) {
        head_alignment.push_back(row.substr(0, left_idx));
        std::string curr_sequences;
        std::string curr_alignment;
        for (size_t col = left_idx; col <= right_idx; ++col) {
            if (row[col] != '-') {  // 跳过 gap 字符
                curr_sequences.push_back(row[col]);
            }
            curr_alignment.push_back(row[col]);
        }
        splited_sequences.push_back(curr_sequences);
        splited_alignment.push_back(curr_alignment);
        tail_alignment.push_back(row.substr(right_idx + 1));
    }
    alignment = head_alignment;
}



#endif //REFINE_P_UTILS_H