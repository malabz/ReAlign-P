# 设置编译器
CXX = g++
CXXFLAGS = -std=c++17 -Wall -O2

# 源文件目录
SRC_DIR = src

# 查找src目录下所有cpp文件
SRC = $(wildcard $(SRC_DIR)/*.cpp)

# 将cpp文件对应为o文件
OBJ = $(SRC:$(SRC_DIR)/%.cpp=$(SRC_DIR)/%.o)

# 目标程序名
TARGET = realign_p

# 默认目标
all: $(TARGET)

# 链接所有对象文件，生成可执行文件
$(TARGET): $(OBJ)
	$(CXX) $(OBJ) -o $@

# 编译每个源文件为目标文件
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 清理目标文件和可执行文件
clean:
	rm -f $(OBJ) $(TARGET)

.PHONY: all clean
