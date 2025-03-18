#ifndef REALIGN_P_SCORE_H
#define REALIGN_P_SCORE_H

#include <string>
#include <numeric>
#include <vector>
#include <map>
#include <iostream>

extern long long aa_gap;  // GAP与氨基酸的得分
extern long long gap_gap; // GAP与GAP的得分
extern std::map<char, std::map<char, long long>> used_score_matrix;
extern std::map<char, std::map<char, long long>> blosum62; // 定义 BLOSUM62 矩阵
extern std::map<char, std::map<char, long long>> pam250; // 定义 PAM250 矩阵

// 定义 PAM250 矩阵
//static const std::unordered_map<char, std::unordered_map<char, int>> pam250;

long long score_column(const std::vector<std::string> sequences, unsigned j);

long long score(const std::vector<std::string> sequences, unsigned l, unsigned r);


#endif //REALIGN_P_SCORE_H
