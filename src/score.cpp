#include "score.h"

long long aa_gap;  // GAP与氨基酸的得分
long long gap_gap; // GAP与GAP的得分
std::map<char, std::map<char, long long>> used_score_matrix;

// 定义 BLOSUM62 矩阵
std::map<char, std::map<char, long long>> blosum62 = {
        {'A', {{'A', 4}, {'R', -1}, {'N', -2}, {'D', -2}, {'C', 0}, {'Q', -1}, {'E', -1}, {'G', 0}, {'H', -1}, {'I', -1}, {'L', -1}, {'K', -1}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -3}, {'Y', -2}, {'V', 0}}},
        {'R', {{'A', -1}, {'R', 5}, {'N', 0}, {'D', -2}, {'C', -3}, {'Q', 1}, {'E', 0}, {'G', -2}, {'H', 0}, {'I', -3}, {'L', -2}, {'K', 2}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -3}, {'Y', -1}, {'V', -2}}},
        {'N', {{'A', -2}, {'R', 0}, {'N', 6}, {'D', 2}, {'C', -3}, {'Q', 0}, {'E', 0}, {'G', -1}, {'H', 1}, {'I', -3}, {'L', -4}, {'K', 0}, {'M', -2}, {'F', -3}, {'P', -2}, {'S', 1}, {'T', 0}, {'W', -4}, {'Y', -2}, {'V', -3}}},
        {'D', {{'A', -2}, {'R', -2}, {'N', 2}, {'D', 6}, {'C', -3}, {'Q', 0}, {'E', 2}, {'G', 0}, {'H', 0}, {'I', -3}, {'L', -4}, {'K', 0}, {'M', -2}, {'F', -3}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -4}, {'Y', -3}, {'V', -3}}},
        {'C', {{'A', 0}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', 9}, {'Q', -3}, {'E', -3}, {'G', -2}, {'H', -2}, {'I', -1}, {'L', -1}, {'K', -2}, {'M', -1}, {'F', -2}, {'P', -2}, {'S', -1}, {'T', -1}, {'W', -2}, {'Y', -2}, {'V', -1}}},
        {'Q', {{'A', -1}, {'R', 1}, {'N', 0}, {'D', 0}, {'C', -3}, {'Q', 5}, {'E', 2}, {'G', -2}, {'H', 0}, {'I', -3}, {'L', -2}, {'K', 1}, {'M', 0}, {'F', -3}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -2}, {'Y', -1}, {'V', -2}}},
        {'E', {{'A', -1}, {'R', 0}, {'N', 0}, {'D', 2}, {'C', -3}, {'Q', 2}, {'E', 5}, {'G', -2}, {'H', 0}, {'I', -3}, {'L', -2}, {'K', 1}, {'M', 0}, {'F', -3}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -3}, {'Y', -2}, {'V', -2}}},
        {'G', {{'A', 0}, {'R', -2}, {'N', -1}, {'D', 0}, {'C', -2}, {'Q', -2}, {'E', -2}, {'G', 6}, {'H', -2}, {'I', -4}, {'L', -4}, {'K', -2}, {'M', -3}, {'F', -3}, {'P', -2}, {'S', 0}, {'T', -1}, {'W', -2}, {'Y', -3}, {'V', -3}}},
        {'H', {{'A', -1}, {'R', 0}, {'N', 1}, {'D', 0}, {'C', -2}, {'Q', 0}, {'E', 0}, {'G', -2}, {'H', 8}, {'I', -3}, {'L', -3}, {'K', 0}, {'M', -2}, {'F', -1}, {'P', -2}, {'S', 0}, {'T', -1}, {'W', -2}, {'Y', 2}, {'V', -3}}},
        {'I', {{'A', -1}, {'R', -3}, {'N', -3}, {'D', -3}, {'C', -1}, {'Q', -3}, {'E', -3}, {'G', -4}, {'H', -3}, {'I', 4}, {'L', 2}, {'K', -3}, {'M', 1}, {'F', 0}, {'P', -2}, {'S', -2}, {'T', -1}, {'W', -3}, {'Y', -2}, {'V', 3}}},
        {'L', {{'A', -1}, {'R', -2}, {'N', -4}, {'D', -4}, {'C', -1}, {'Q', -2}, {'E', -2}, {'G', -4}, {'H', -3}, {'I', 2}, {'L', 4}, {'K', -2}, {'M', 2}, {'F', 0}, {'P', -2}, {'S', -2}, {'T', -1}, {'W', -2}, {'Y', -2}, {'V', 1}}},
        {'K', {{'A', -1}, {'R', 2}, {'N', 0}, {'D', 0}, {'C', -2}, {'Q', 1}, {'E', 1}, {'G', -2}, {'H', 0}, {'I', -3}, {'L', -2}, {'K', 5}, {'M', 0}, {'F', -3}, {'P', -1}, {'S', 0}, {'T', -1}, {'W', -2}, {'Y', -1}, {'V', -2}}},
        {'M', {{'A', -1}, {'R', -1}, {'N', -2}, {'D', -2}, {'C', -1}, {'Q', 0}, {'E', 0}, {'G', -3}, {'H', -2}, {'I', 1}, {'L', 2}, {'K', 0}, {'M', 6}, {'F', -1}, {'P', -2}, {'S', -1}, {'T', -1}, {'W', -3}, {'Y', -1}, {'V', -1}}},
        {'F', {{'A', -2}, {'R', -2}, {'N', -3}, {'D', -3}, {'C', -2}, {'Q', -3}, {'E', -3}, {'G', -3}, {'H', -1}, {'I', 0}, {'L', 0}, {'K', -3}, {'M', -1}, {'F', 6}, {'P', -3}, {'S', -3}, {'T', -2}, {'W', 1}, {'Y', 3}, {'V', -1}}},
        {'P', {{'A', -1}, {'R', -1}, {'N', -2}, {'D', -1}, {'C', -2}, {'Q', -1}, {'E', -1}, {'G', -2}, {'H', -2}, {'I', -2}, {'L', -2}, {'K', -1}, {'M', -2}, {'F', -3}, {'P', 7}, {'S', -1}, {'T', -1}, {'W', -4}, {'Y', -3}, {'V', -2}}},
        {'S', {{'A', 0}, {'R', 0}, {'N', 1}, {'D', 0}, {'C', -1}, {'Q', 0}, {'E', 0}, {'G', 0}, {'H', 0}, {'I', -2}, {'L', -2}, {'K', 0}, {'M', -1}, {'F', -3}, {'P', -1}, {'S', 4}, {'T', 1}, {'W', -3}, {'Y', -2}, {'V', -2}}},
        {'T', {{'A', -1}, {'R', -1}, {'N', 0}, {'D', -1}, {'C', -1}, {'Q', -1}, {'E', -1}, {'G', -1}, {'H', -1}, {'I', -1}, {'L', -1}, {'K', -1}, {'M', -1}, {'F', -2}, {'P', -1}, {'S', 1}, {'T', 5}, {'W', -3}, {'Y', -2}, {'V', 0}}},
        {'W', {{'A', -3}, {'R', -3}, {'N', -4}, {'D', -4}, {'C', -2}, {'Q', -2}, {'E', -3}, {'G', -2}, {'H', -2}, {'I', -3}, {'L', -2}, {'K', -2}, {'M', -3}, {'F', 1}, {'P', -4}, {'S', -3}, {'T', -3}, {'W', 11}, {'Y', 2}, {'V', -3}}},
        {'Y', {{'A', -2}, {'R', -1}, {'N', -2}, {'D', -3}, {'C', -2}, {'Q', -1}, {'E', -2}, {'G', -3}, {'H', 2}, {'I', -2}, {'L', -2}, {'K', -1}, {'M', -1}, {'F', 3}, {'P', -3}, {'S', -2}, {'T', -2}, {'W', 2}, {'Y', 7}, {'V', -2}}},
        {'V', {{'A', 0}, {'R', -2}, {'N', -3}, {'D', -3}, {'C', -1}, {'Q', -2}, {'E', -2}, {'G', -3}, {'H', -3}, {'I', 3}, {'L', 1}, {'K', -2}, {'M', -1}, {'F', -1}, {'P', -2}, {'S', -2}, {'T', 0}, {'W', -3}, {'Y', -2}, {'V', 4}}},
};

std::map<char, std::map<char, long long>> pam250 = {
        {'A', {{'A', 2}, {'R', -2}, {'N', 0}, {'D', 0}, {'C', -2}, {'Q', 0}, {'E', 0}, {'G', 1}, {'H', -1}, {'I', -1}, {'L', -2}, {'K', -1}, {'M', -1}, {'F', -3}, {'P', 1}, {'S', 1}, {'T', 1}, {'W', -6}, {'Y', -3}, {'V', 0}}},
        {'R', {{'A', -2}, {'R', 6}, {'N', 0}, {'D', -1}, {'C', -4}, {'Q', 1}, {'E', -1}, {'G', -3}, {'H', 2}, {'I', -2}, {'L', -3}, {'K', 3}, {'M', 0}, {'F', -4}, {'P', 0}, {'S', 0}, {'T', -1}, {'W', 2}, {'Y', -4}, {'V', -2}}},
        {'N', {{'A', 0}, {'R', 0}, {'N', 2}, {'D', 2}, {'C', -4}, {'Q', 1}, {'E', 1}, {'G', 0}, {'H', 2}, {'I', -2}, {'L', -3}, {'K', 1}, {'M', -2}, {'F', -3}, {'P', 0}, {'S', 1}, {'T', 0}, {'W', -4}, {'Y', -2}, {'V', -2}}},
        {'D', {{'A', 0}, {'R', -1}, {'N', 2}, {'D', 4}, {'C', -5}, {'Q', 2}, {'E', 3}, {'G', 1}, {'H', 1}, {'I', -2}, {'L', -4}, {'K', 0}, {'M', -3}, {'F', -6}, {'P', -1}, {'S', 0}, {'T', 0}, {'W', -7}, {'Y', -4}, {'V', -2}}},
        {'C', {{'A', -2}, {'R', -4}, {'N', -4}, {'D', -5}, {'C', 12}, {'Q', -5}, {'E', -5}, {'G', -3}, {'H', -3}, {'I', -2}, {'L', -6}, {'K', -5}, {'M', -5}, {'F', -4}, {'P', -3}, {'S', 0}, {'T', -2}, {'W', -8}, {'Y', 0}, {'V', -2}}},
        {'Q', {{'A', 0}, {'R', 1}, {'N', 1}, {'D', 2}, {'C', -5}, {'Q', 4}, {'E', 2}, {'G', -1}, {'H', 3}, {'I', -2}, {'L', -2}, {'K', 1}, {'M', -1}, {'F', -5}, {'P', 0}, {'S', -1}, {'T', -1}, {'W', -5}, {'Y', -4}, {'V', -2}}},
        {'E', {{'A', 0}, {'R', -1}, {'N', 1}, {'D', 3}, {'C', -5}, {'Q', 2}, {'E', 4}, {'G', 0}, {'H', 1}, {'I', -2}, {'L', -3}, {'K', 0}, {'M', -2}, {'F', -5}, {'P', -1}, {'S', 0}, {'T', 0}, {'W', -7}, {'Y', -4}, {'V', -2}}},
        {'G', {{'A', 1}, {'R', -3}, {'N', 0}, {'D', 1}, {'C', -3}, {'Q', -1}, {'E', 0}, {'G', 5}, {'H', -2}, {'I', -3}, {'L', -4}, {'K', -2}, {'M', -3}, {'F', -5}, {'P', 0}, {'S', 1}, {'T', 0}, {'W', -7}, {'Y', -5}, {'V', -1}}},
        {'H', {{'A', -1}, {'R', 2}, {'N', 2}, {'D', 1}, {'C', -3}, {'Q', 3}, {'E', 1}, {'G', -2}, {'H', 6}, {'I', -2}, {'L', -2}, {'K', 0}, {'M', -2}, {'F', -2}, {'P', 0}, {'S', -1}, {'T', -1}, {'W', -3}, {'Y', 0}, {'V', -2}}},
        {'I', {{'A', -1}, {'R', -2}, {'N', -2}, {'D', -2}, {'C', -2}, {'Q', -2}, {'E', -2}, {'G', -3}, {'H', -2}, {'I', 5}, {'L', 2}, {'K', -2}, {'M', 2}, {'F', 1}, {'P', -2}, {'S', -1}, {'T', 0}, {'W', -5}, {'Y', -1}, {'V', 4}}},
        {'L', {{'A', -2}, {'R', -3}, {'N', -3}, {'D', -4}, {'C', -6}, {'Q', -2}, {'E', -3}, {'G', -4}, {'H', -2}, {'I', 2}, {'L', 6}, {'K', -3}, {'M', 4}, {'F', 2}, {'P', -3}, {'S', -3}, {'T', -2}, {'W', -2}, {'Y', -1}, {'V', 2}}},
        {'K', {{'A', -1}, {'R', 3}, {'N', 1}, {'D', 0}, {'C', -5}, {'Q', 1}, {'E', 0}, {'G', -2}, {'H', 0}, {'I', -2}, {'L', -3}, {'K', 5}, {'M', 0}, {'F', -5}, {'P', -1}, {'S', 0}, {'T', 0}, {'W', -3}, {'Y', -4}, {'V', -2}}},
        {'M', {{'A', -1}, {'R', 0}, {'N', -2}, {'D', -3}, {'C', -5}, {'Q', -1}, {'E', -2}, {'G', -3}, {'H', -2}, {'I', 2}, {'L', 4}, {'K', 0}, {'M', 6}, {'F', 0}, {'P', -2}, {'S', -2}, {'T', -1}, {'W', -4}, {'Y', -2}, {'V', 2}}},
        {'F', {{'A', -3}, {'R', -4}, {'N', -3}, {'D', -6}, {'C', -4}, {'Q', -5}, {'E', -5}, {'G', -5}, {'H', -2}, {'I', 1}, {'L', 2}, {'K', -5}, {'M', 0}, {'F', 9}, {'P', -5}, {'S', -3}, {'T', -3}, {'W', 0}, {'Y', 7}, {'V', -1}}},
        {'P', {{'A', 1}, {'R', 0}, {'N', 0}, {'D', -1}, {'C', -3}, {'Q', 0}, {'E', -1}, {'G', 0}, {'H', 0}, {'I', -2}, {'L', -3}, {'K', -1}, {'M', -2}, {'F', -5}, {'P', 6}, {'S', 1}, {'T', 0}, {'W', -6}, {'Y', -5}, {'V', -1}}},
        {'S', {{'A', 1}, {'R', 0}, {'N', 1}, {'D', 0}, {'C', 0}, {'Q', -1}, {'E', 0}, {'G', 1}, {'H', -1}, {'I', -1}, {'L', -3}, {'K', 0}, {'M', -2}, {'F', -3}, {'P', 1}, {'S', 2}, {'T', 1}, {'W', -2}, {'Y', -3}, {'V', -1}}},
        {'T', {{'A', 1}, {'R', -1}, {'N', 0}, {'D', 0}, {'C', -2}, {'Q', -1}, {'E', 0}, {'G', 0}, {'H', -1}, {'I', 0}, {'L', -2}, {'K', 0}, {'M', -1}, {'F', -3}, {'P', 0}, {'S', 1}, {'T', 3}, {'W', -5}, {'Y', -3}, {'V', 0}}},
        {'W', {{'A', -6}, {'R', 2}, {'N', -4}, {'D', -7}, {'C', -8}, {'Q', -5}, {'E', -7}, {'G', -7}, {'H', -3}, {'I', -5}, {'L', -2}, {'K', -3}, {'M', -4}, {'F', 0}, {'P', -6}, {'S', -2}, {'T', -5}, {'W', 17}, {'Y', 0}, {'V', -6}}},
        {'Y', {{'A', -3}, {'R', -4}, {'N', -2}, {'D', -4}, {'C', 0}, {'Q', -4}, {'E', -4}, {'G', -5}, {'H', 0}, {'I', -1}, {'L', -1}, {'K', -4}, {'M', -2}, {'F', 7}, {'P', -5}, {'S', -3}, {'T', -3}, {'W', 0}, {'Y', 10}, {'V', -2}}},
        {'V', {{'A', 0}, {'R', -2}, {'N', -2}, {'D', -2}, {'C', -2}, {'Q', -2}, {'E', -2}, {'G', -1}, {'H', -2}, {'I', 4}, {'L', 2}, {'K', -2}, {'M', 2}, {'F', -1}, {'P', -1}, {'S', -1}, {'T', 0}, {'W', -6}, {'Y', -2}, {'V', 4}}}
 };


long long score_column(const std::vector<std::string> sequences, unsigned j) {
    long long score = 0;

    // 遍历所有的序列，计算该位置的得分
    for (size_t i = 0; i < sequences.size(); ++i) {
        for (size_t k = i + 1; k < sequences.size(); ++k) {
            char a = sequences[i][j];
            char b = sequences[k][j];

            // 处理GAP的情况
            if (a == '-' && b == '-') {
                score += gap_gap;
            } else if (a == '-' || b == '-') {
                score += aa_gap;
            } else {
                // 使用BLOSUM62矩阵计算氨基酸对的得分
                score += used_score_matrix[a][b];
            }
        }
    }

    return score;
}

long long score(const std::vector<std::string> sequences, unsigned l, unsigned r)
{
    long long s = 0;
    for (unsigned i = l; i != r; ++i)
        s += score_column(sequences, i);
    return s;
}