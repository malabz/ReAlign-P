#include "run_msa.h"

void utils::run_msa::run_muscle5() const {
    try {
        // 构造命令字符串
        std::string command = "muscle -align " + input_fasta + " -output " + output_fasta;

        // 执行命令
        int result = std::system(command.c_str());

        // 检查命令执行结果
        if (result != 0) {
            throw std::runtime_error("Error in MUSCLE5 alignment.");
        }

        // 输出成功消息
        std::cout << "Alignment completed. Output saved to " << output_fasta << "." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
}

void utils::run_msa::run_muscle3() const {
    try {
        // 构造命令字符串
        std::string command = "muscle -in " + input_fasta + " -out " + output_fasta;

        // 执行命令
        int result = std::system(command.c_str());

        // 检查命令执行结果
        if (result != 0) {
            throw std::runtime_error("Error in MUSCLE3 alignment.");
        }

        // 输出成功消息
        std::cout << "Alignment completed. Output saved to " << output_fasta << "." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
}

void utils::run_msa::run_fftnsi() const {
    try {
        // 构造命令字符串
        std::string command = "fftnsi " + input_fasta + " > " + output_fasta;

        // 执行命令
        int result = std::system(command.c_str());

        // 检查命令执行结果
        if (result != 0) {
            throw std::runtime_error("Error in MAFFT alignment.");
        }

        // 输出成功消息
        std::cout << "Alignment completed. Output saved to " << output_fasta << "." << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
    }
}