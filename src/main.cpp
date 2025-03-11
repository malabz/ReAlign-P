#include <iostream>
#include "fasta.h"
#include "utils.h"
#include "score.h"
#include "iterative.h"

int main(int argc, char **argv) {

    std::string input_file = argv[1];
    std::string output_file = argv[2];

    std::vector<std::string> final_sequence;
    utils::fasta alignment = read_from(input_file);

//    加载蛋白质矩阵

//  确定列是否包含 gap
    std::vector<int> columns_summary = column_status(alignment.sequences);

    size_t left_idx, right_idx;
    locate_middle_part(left_idx, right_idx, alignment.sequences, columns_summary);

    std::vector<std::string> middle_sequences, middle_alignment, tail_alignment;
    split_alignment(left_idx, right_idx, alignment.sequences, middle_sequences, middle_alignment, tail_alignment);

    std::vector<std::string> realigned_middle_alignment = core_iterative(alignment.identifications, middle_sequences, middle_alignment);

    concatenate_sequences(alignment.sequences, realigned_middle_alignment);
    concatenate_sequences(alignment.sequences, tail_alignment);

    std::ofstream ofs(output_file);
    if (!ofs) {
        std::cerr << "Error: cannot open file " << output_file << std::endl;
        exit(1);
    }
    alignment.write_to(ofs);
    ofs.close();

    return 0;
}
