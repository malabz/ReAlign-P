#include <iostream>
#include "fasta.h"
#include "utils.h"
#include "score.h"
#include "iterative.h"

int main(int argc, char **argv) {

    // Check if the program was called with no arguments or with the "-h" option
    if (argc == 1 || (argc == 2 && std::string(argv[1]) == "-h")) {
        displayHelp();
        return 0;
    }

    std::string input_file, score_matrix, aa_gap_score, gap_gap_score, msa;
    std::string output_file = "realign_p_result.fasta";


    bool have_score_matrix = false;
    bool have_aa_gap_score = false;
    bool have_gap_gap_score = false;
    bool have_msa = false;

    for (int i = 1; i < argc; i += 2) {
        std::string option = argv[i];
        if (i + 1 < argc) {
            std::string value = argv[i + 1];
            if (option == "-i") {
                input_file = value;
            } else if (option == "-o") {
                output_file = value;
            } else if (option == "-s") {
                have_score_matrix = true;
                score_matrix = value;
            } else if (option == "-a") {
                have_aa_gap_score = true;
                aa_gap_score = value;
            } else if (option == "-g") {
                have_gap_gap_score = true;
                gap_gap_score = value;
            } else if (option == "-m") {
                have_msa = true;
                msa = value;
            } else {
                std::cerr << "** Unknown option: " << option << std::endl;
                displayHelp();
                return 1;
            }
        } else {
            std::cerr << "** Missing value for option: " << option << std::endl;
            displayHelp();
            return 1;
        }
    }

    if (! have_score_matrix) {
        score_matrix = "blosum62";
    }

    if (! have_aa_gap_score) {
        aa_gap_score = "-2";
    }

    if (! have_gap_gap_score) {
        gap_gap_score = "0";
    }

    if (! have_msa) {
        msa = "muscle5";
    }

    // Check if required -i option is provided
    if (input_file.empty()) {
        std::cerr << "** Error: -i option is required." << std::endl;
        displayHelp();
        return 1;
    }

    if (msa != "mafft" && msa != "muscle5") {
        std::cerr << "** Error: This MSA tool is not supported." << std::endl;
        displayHelp();
        return 1;
    }

    if (score_matrix != "blosum62" && score_matrix != "pam250") {
        std::cerr << "** Error: This score matrix is not supported." << std::endl;
        displayHelp();
        return 1;
    }

    // 加载蛋白质矩阵和罚分
    if (score_matrix == "blosum62")
    {
        used_score_matrix = blosum62;
    }

    if (score_matrix == "pam250")
    {
        used_score_matrix = pam250;
    }

    aa_gap = atoll(aa_gap_score.c_str());
    gap_gap = atoll(gap_gap_score.c_str());

    used_msa = msa;

//  确定列是否包含 gap

    std::vector<std::string> final_sequence;
    utils::fasta alignment = read_from(input_file);

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
