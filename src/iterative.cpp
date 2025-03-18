#include "iterative.h"

std::string used_msa;

void preprocess(const std::vector<std::string>& identifications, const std::vector<std::string>& sequences,
                utils::fasta& realigned_fasta)
{
    // 确保文件存储在固定位置，如 "realign_p_tmp" 文件夹
    std::filesystem::create_directory("realign_p_tmp");
    std::vector<std::string> empty_idx;
    std::vector<std::string> non_empty_idx;
    std::vector<std::string> non_empty_seq;

    bool contains_empty = false;

    // 遍历 sequences，筛选非空序列
    for (size_t i = 0; i < identifications.size(); ++i) {
        if (sequences[i].empty()) {
            empty_idx.push_back(identifications[i]);
            contains_empty = true;
        } else {
            non_empty_idx.push_back(identifications[i]);
            non_empty_seq.push_back(sequences[i]);
        }
    }

    // 随机生成文件名
    std::srand(std::time(0));
    std::string random_filename = std::to_string(std::rand());

    std::string input_name = "realign_p_tmp/" + random_filename + ".fasta";
    std::string output_name = "realign_p_tmp/" + random_filename + "_aligned.fasta";

    // 创建并写入文件
    utils::fasta non_empty_fasta;
    non_empty_fasta.identifications = non_empty_idx;
    non_empty_fasta.sequences = non_empty_seq;

    std::ofstream ofs(input_name);
    if (!ofs) {
        std::cerr << "Error: cannot open file " << input_name << std::endl;
        exit(1);
    }
    non_empty_fasta.write_to(ofs);
    ofs.close();

    // 使用 MUSCLE 重新比对
    utils::run_msa msa(input_name, output_name);
    if (used_msa == "mafft")
    {
        msa.run_fftnsi();
    }
    if (used_msa == "muscle5")
    {
        msa.run_muscle5();
    }

    // 读取重新比对后的数据
    realigned_fasta = read_from(output_name);

    // 如果有缺失的序列，插入填充
    if (contains_empty) {
        for (size_t i = 0; i < empty_idx.size(); ++i) {
            realigned_fasta.identifications.push_back(empty_idx[i]);
            realigned_fasta.sequences.push_back(std::string(realigned_fasta.sequences[0].size(), '-'));
        }
    }
}

std::vector<std::string> reorder_sequences(const std::vector<std::string>& original_ids,
                                           const utils::fasta& realigned_fasta)
{
    // 根据 original_ids 的顺序重排 aligned_sequences
    std::vector<std::string> reordered_sequences;
    for (const auto& id : original_ids) {
        // 查找 ID 在 realigned_fasta.identifications 中的索引
        auto it = std::find(realigned_fasta.identifications.begin(),
                            realigned_fasta.identifications.end(), id);

        if (it != realigned_fasta.identifications.end()) {
            size_t index = std::distance(realigned_fasta.identifications.begin(), it);
            reordered_sequences.push_back(realigned_fasta.sequences[index]);
        } else {
            reordered_sequences.push_back("");  // 如果找不到对应的 ID，返回空字符串
        }
    }

    return reordered_sequences;
}

std::vector<std::vector<int>> recoded_alignment(const std::vector<std::string>& alignment)
{
    // 获取序列数和列数
    size_t num_sequences = alignment.size();
    size_t num_columns = alignment[0].size();

    // 初始化位置矩阵，默认所有位置为 0
    std::vector<std::vector<int>> position_matrix(num_sequences, std::vector<int>(num_columns, 0));

    // 遍历每行并更新位置矩阵
    for (size_t row_idx = 0; row_idx < num_sequences; ++row_idx) {
        int position = 1;
        for (size_t col_idx = 0; col_idx < num_columns; ++col_idx) {
            if (alignment[row_idx][col_idx] != '-') {  // 跳过gap字符
                position_matrix[row_idx][col_idx] = position;
                ++position;  // 位置自增
            }
        }
    }

    return position_matrix;
}

std::vector<std::pair<int, int>> find_identical_columns(const std::vector<std::vector<int>>& matrix1, const std::vector<std::vector<int>>& matrix2) {
    std::vector<std::pair<int, int>> identical_columns;

    // 获取两个矩阵的列数
    size_t num_cols1 = matrix1[0].size();
    size_t num_cols2 = matrix2[0].size();

    // 遍历第一个矩阵的所有列
    for (size_t col_idx1 = 0; col_idx1 < num_cols1; ++col_idx1) {
        // 遍历第二个矩阵的所有列
        for (size_t col_idx2 = 0; col_idx2 < num_cols2; ++col_idx2) {
            bool identical = true;

            // 检查两列是否相同
            for (size_t row_idx = 0; row_idx < matrix1.size(); ++row_idx) {
                if (matrix1[row_idx][col_idx1] != matrix2[row_idx][col_idx2]) {
                    identical = false;
                    break;
                }
            }

            // 如果列相同，记录列的索引对
            if (identical) {
                identical_columns.emplace_back(col_idx1, col_idx2);
            }
        }
    }

    return identical_columns;
}

// 字符串链接
void concatenate_sequences(std::vector<std::string>& concatenated_result, const std::vector<std::string>& part_to_add) {
    // 如果final_middle为空，直接赋值为part_to_add的内容
    if (concatenated_result.empty())
    {
        concatenated_result = part_to_add;  // 如果final_middle为空，直接赋值为part_to_add
    }
    else
    {
//        // 否则将final_middle和part_to_add中的序列连接
        size_t size = part_to_add.size();
//        final_middle.reserve(final_middle.size() + size); // 预分配内存，避免多次动态分配
        for (size_t i = 0; i < size; ++i) {
            concatenated_result[i] = concatenated_result[i] + part_to_add[i];  // 连接字符串并添加到 final_middle
        }
    }
}

void get_aligned_intervals(const std::vector<std::pair<int, int>>& index_pairs,
                           std::vector<std::pair<int, int>>& raw_intervals,
                           std::vector<std::pair<int, int>>& new_intervals) {
    if (index_pairs.empty()) {
        return;
    }

    int start1 = index_pairs[0].first, start2 = index_pairs[0].second; // 当前区间的起点
    std::vector<int> temp1 = {start1}, temp2 = {start2}; // 临时存储连续区间的索引

    for (size_t i = 1; i < index_pairs.size(); ++i) {
        int cur1 = index_pairs[i].first, cur2 = index_pairs[i].second;
        int prev1 = index_pairs[i - 1].first, prev2 = index_pairs[i - 1].second;

        // 如果两个索引都连续
        if (cur1 == prev1 + 1 && cur2 == prev2 + 1) {
            temp1.push_back(cur1);
            temp2.push_back(cur2);
        } else {
            // 只有当区间长度 >= 5 时才加入
            if (temp1.size() >= 5) {
                raw_intervals.push_back({temp1.front(), temp1.back()});
                new_intervals.push_back({temp2.front(), temp2.back()});
            }
            // 开始新的区间
            temp1 = {cur1};
            temp2 = {cur2};
        }
    }

    // 处理最后一个区间
    if (temp1.size() >= 5) {
        raw_intervals.push_back({temp1.front(), temp1.back()});
        new_intervals.push_back({temp2.front(), temp2.back()});
    }
}

std::vector<std::string> extract_columns(const std::vector<std::string>& alignment, size_t a, size_t b) {
    size_t num_sequences = alignment.size();
    std::vector<std::string> extracted_columns(num_sequences);

    // 遍历所有序列，截取每个序列从列 a 到列 b 的部分
    for (size_t i = 0; i < num_sequences; ++i) {
        extracted_columns[i] = alignment[i].substr(a, b - a + 1); // substr(start, length)
    }

    return extracted_columns;
}

bool is_all_gaps(const std::string& seq) {
    return std::all_of(seq.begin(), seq.end(), [](char c) { return c == '-'; });
}

std::vector<std::string> remove_dashes(const std::vector<std::string>& strings) {
    std::vector<std::string> result;
    for (const auto& s : strings) {
        std::string modified_str = s;
        modified_str.erase(std::remove(modified_str.begin(), modified_str.end(), '-'), modified_str.end());
        result.push_back(modified_str);
    }
    return result;
}

void choose_better_alignment_and_concatenate(const std::vector<std::string>& identifications,
                                             std::vector<std::string>& realigned_alignment,
                                             const std::vector<std::string>& diff_raw,
                                             const std::vector<std::string>& diff_new)
{
    // 计算 raw_sp 和 new_sp
    long long raw_sp = score(diff_raw, 0, diff_raw[0].size());
    long long new_sp = score(diff_new, 0, diff_new[0].size());

    if (new_sp > raw_sp) {
        // 如果 new_sp 更大，处理 diff_new
        std::vector<std::string> diff_seq_new = remove_dashes(diff_new);
        std::vector<std::string> realigned_new = core_iterative(identifications, diff_seq_new, diff_new);
        concatenate_sequences(realigned_alignment, realigned_new);
    }
    else {
        // 否则，处理 diff_raw
        std::vector<std::string> diff_seq_raw = remove_dashes(diff_raw);
        std::vector<std::string> realigned_raw = core_iterative(identifications, diff_seq_raw, diff_raw);
        concatenate_sequences(realigned_alignment, realigned_raw);
    }
}

std::vector<std::string> core_iterative(const std::vector<std::string>& identifications, const std::vector<std::string>& sequences,
                                        const std::vector<std::string>& alignment)
{

    std::vector<std::string> realigned_alignment;

    utils::fasta realigned_fasta;
    preprocess(identifications, sequences, realigned_fasta);
    std::vector<std::string> reordered_alignment = reorder_sequences(identifications, realigned_fasta);


    std::vector<std::vector<int>> recoded_raw_alignment = recoded_alignment(alignment);
    std::vector<std::vector<int>> recoded_new_alignment = recoded_alignment(reordered_alignment);
    auto identical_columns = find_identical_columns(recoded_raw_alignment, recoded_new_alignment);

    if (identical_columns.empty())
    {
        long long raw_sp = score(alignment, 0, alignment[0].size());
        long long new_sp = score(reordered_alignment, 0, reordered_alignment[0].size());
        if (new_sp > raw_sp)
        {
            concatenate_sequences(realigned_alignment, reordered_alignment);
        }
        else
        {
            concatenate_sequences(realigned_alignment, alignment);
        }
        return realigned_alignment;
    }

    std::vector<std::pair<int, int>> raw_intervals, new_intervals;
    get_aligned_intervals(identical_columns, raw_intervals, new_intervals);

    std::vector<std::pair<int, int>> correct_raw_region, correct_new_region;
    for (int i = 0; i < raw_intervals.size(); ++i)
    {
        bool is_correct = true;
        std::vector<std::string> part_seq = extract_columns(alignment, raw_intervals[i].first, raw_intervals[i].second);
        for (auto & j: part_seq)
        {
            if (is_all_gaps(j))
            {
                is_correct = false;
                break;
            }
        }
        if (is_correct)
        {
            correct_raw_region.push_back(raw_intervals[i]);
            correct_new_region.push_back(new_intervals[i]);
        }
    }

    std::vector<std::string> same_part, diff_raw, diff_new;
    if (correct_raw_region.empty())
    {
        long long raw_sp = score(alignment, 0, alignment[0].size());
        long long new_sp = score(reordered_alignment, 0, reordered_alignment[0].size());
        if (new_sp > raw_sp)
        {
            concatenate_sequences(realigned_alignment, reordered_alignment);
        }
        else
        {
            concatenate_sequences(realigned_alignment, alignment);
        }

    }
    else if (correct_raw_region.size() == 1)
    {
        if ((correct_raw_region[0].first == 0) and (correct_raw_region[0].second == alignment[0].size() - 1))
        {
            concatenate_sequences(realigned_alignment, alignment);
        }
        else if (correct_raw_region[0].first == 0)
        {
            same_part = extract_columns(alignment, correct_raw_region[0].first, correct_raw_region[0].second);
            concatenate_sequences(realigned_alignment, same_part);

            diff_raw = extract_columns(alignment, correct_raw_region[0].second + 1, alignment[0].size() - 1);
            diff_new = extract_columns(reordered_alignment, correct_new_region[0].second + 1, reordered_alignment[0].size() - 1);
            choose_better_alignment_and_concatenate(identifications, realigned_alignment, diff_raw, diff_new);
        }
        else if (correct_raw_region[0].second == alignment[0].size() - 1)
        {
            diff_raw = extract_columns(alignment, 0, correct_raw_region[0].first - 1);
            diff_new = extract_columns(reordered_alignment, 0, correct_new_region[0].first - 1);
            choose_better_alignment_and_concatenate(identifications, realigned_alignment, diff_raw, diff_new);

            same_part = extract_columns(alignment, correct_raw_region[0].first, correct_raw_region[0].second);
            concatenate_sequences(realigned_alignment, same_part);
        }
        else
        {
            diff_raw = extract_columns(alignment, 0, correct_raw_region[0].first - 1);
            diff_new = extract_columns(reordered_alignment, 0, correct_new_region[0].first - 1);
            choose_better_alignment_and_concatenate(identifications, realigned_alignment, diff_raw, diff_new);

            same_part = extract_columns(alignment, correct_raw_region[0].first, correct_raw_region[0].second);
            concatenate_sequences(realigned_alignment, same_part);

            diff_raw = extract_columns(alignment, correct_raw_region[0].second + 1, alignment[0].size() - 1);
            diff_new = extract_columns(reordered_alignment, correct_new_region[0].second + 1, reordered_alignment[0].size() - 1);
            choose_better_alignment_and_concatenate(identifications, realigned_alignment, diff_raw, diff_new);
        }
    }
    else
    {
        for (int i = 0; i < correct_raw_region.size(); ++i)
        {
            if (i == 0)
            {
                if (correct_raw_region[i].first == 0)
                {
                    same_part = extract_columns(alignment, correct_raw_region[i].first, correct_raw_region[i].second);
                    concatenate_sequences(realigned_alignment, same_part);

                    diff_raw = extract_columns(alignment, correct_raw_region[i].second + 1, correct_raw_region[i+1].first - 1);
                    diff_new = extract_columns(reordered_alignment, correct_new_region[i].second + 1, correct_new_region[i+1].first - 1);
                    choose_better_alignment_and_concatenate(identifications, realigned_alignment, diff_raw, diff_new);
                }
                else
                {
                    diff_raw = extract_columns(alignment, 0, correct_raw_region[i].first - 1);
                    diff_new = extract_columns(reordered_alignment, 0, correct_new_region[i].first - 1);
                    choose_better_alignment_and_concatenate(identifications, realigned_alignment, diff_raw, diff_new);

                    same_part = extract_columns(alignment, correct_raw_region[i].first, correct_raw_region[i].second);
                    concatenate_sequences(realigned_alignment, same_part);

                    diff_raw = extract_columns(alignment, correct_raw_region[i].second + 1, correct_raw_region[i+1].first - 1);
                    diff_new = extract_columns(reordered_alignment, correct_new_region[i].second + 1, correct_new_region[i+1].first - 1);
                    choose_better_alignment_and_concatenate(identifications, realigned_alignment, diff_raw, diff_new);
                }
            }
            else if (i == correct_raw_region.size() - 1)
            {
                same_part = extract_columns(alignment, correct_raw_region[i].first, correct_raw_region[i].second);
                concatenate_sequences(realigned_alignment, same_part);

                if (correct_raw_region[i].second != alignment[0].size() - 1)
                {
                    diff_raw = extract_columns(alignment, correct_raw_region[i].second + 1, alignment[0].size() - 1);
                    diff_new = extract_columns(reordered_alignment, correct_new_region[i].second + 1, reordered_alignment[0].size() - 1);
                    choose_better_alignment_and_concatenate(identifications, realigned_alignment, diff_raw, diff_new);
                }
            }
            else
            {
                same_part = extract_columns(alignment, correct_raw_region[i].first, correct_raw_region[i].second);
                concatenate_sequences(realigned_alignment, same_part);

                diff_raw = extract_columns(alignment, correct_raw_region[i].second + 1, correct_raw_region[i+1].first - 1);
                diff_new = extract_columns(reordered_alignment, correct_new_region[i].second + 1, correct_new_region[i+1].first - 1);
                choose_better_alignment_and_concatenate(identifications, realigned_alignment, diff_raw, diff_new);
            }
        }
    }

    return realigned_alignment;
}