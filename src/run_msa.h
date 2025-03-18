#ifndef REALIGN_P_RUN_MSA_H
#define REALIGN_P_RUN_MSA_H

#include <iostream>
#include <stdexcept>
#include <string>
#include <cstdlib>  // For system()
#include <cstdio>   // For perror()

namespace utils
{
    class run_msa
    {
    public:
        std::string input_fasta, output_fasta;

        // 构造函数，用于初始化 input_fasta 和 output_fasta
        run_msa(const std::string& input, const std::string& output)
                : input_fasta(input), output_fasta(output)
        {}
        void run_muscle5() const;
        void run_muscle3() const;
        void run_fftnsi() const;
    };
}


#endif //REALIGN_P_RUN_MSA_H
