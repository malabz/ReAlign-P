#ifndef REALIGN_P_ITERATIVE_H
#define REALIGN_P_ITERATIVE_H

#include <filesystem>
#include <ctime>
#include <fstream>
#include "fasta.h"
#include "utils.h"
#include "score.h"
#include "run_msa.h"

bool is_all_gaps(const std::string& seq);

std::vector<std::string> remove_dashes(const std::vector<std::string>& strings);

void preprocess(const std::vector<std::string>& identifications, const std::vector<std::string>& sequences,
                utils::fasta& realigned_fasta);

std::vector<std::string> reorder_sequences(const std::vector<std::string>& original_ids,
                                           const utils::fasta& realigned_fasta);

std::vector<std::vector<int>> recoded_alignment(const std::vector<std::string>& alignment);

std::vector<std::pair<int, int>> find_identical_columns(const std::vector<std::vector<int>>& matrix1, const std::vector<std::vector<int>>& matrix2);

void get_aligned_intervals(const std::vector<std::pair<int, int>>& index_pairs,
                           std::vector<std::pair<int, int>>& raw_intervals,
                           std::vector<std::pair<int, int>>& new_intervals);

std::vector<std::string> extract_columns(const std::vector<std::string>& alignment, size_t a, size_t b);

void concatenate_sequences(std::vector<std::string>& concatenated_result, const std::vector<std::string>& part_to_add);

void choose_better_alignment_and_concatenate(const std::vector<std::string>& identifications, std::vector<std::string>& realigned_alignment,
                                             const std::vector<std::string>& diff_raw, const std::vector<std::string>& diff_new);

std::vector<std::string> core_iterative(const std::vector<std::string>& identifications, const std::vector<std::string>& sequences,
                                        const std::vector<std::string>& alignment);


#endif //REALIGN_P_ITERATIVE_H
