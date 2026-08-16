#include <vector>
#include <cstdint>
#include <fstream>
#include <string>
#include <iostream>
#include <cassert>

//Type alias for std::vector<std::vector<uint8_t>>
//which is the 2D-vector array that is being made
using Panel = std::vector<std::vector<uint8_t>>;

using prefixArrs = std::vector<std::vector<uint32_t>>;

prefixArrs buildPrefixArrays(const Panel& panel, size_t M, size_t N){
    prefixArrs result(N+1, std::vector<uint32_t>(M));
    
    //needs comments
    for(size_t i{0}; i < M; i++){
        result[0][i] = i;
    }

    std::vector<uint32_t> zeroes, ones;
    for(size_t k{0}; k < N; k++){

        //Resets sizes to 0 but keeps underlying buffer capacity
        zeroes.clear();
        ones.clear();

        for(size_t x{0}; x < M; x++){
            size_t hap = result[k][x];
            if(panel[k][hap] == 0){
                zeroes.push_back(hap);
            }
            else{
                ones.push_back(hap);
            }
        }
        size_t pos = 0;
        for(size_t i{0}; i < zeroes.size(); i++){
            result[k+1][pos] = zeroes[i];
            pos++;
        }
        for(size_t j{0}; j < ones.size(); j++){
            result[k+1][pos] = ones[j];
            pos++;
        }
    }

    return result;
}

//Loads 2D-array
Panel loadPanel(std::string path, size_t& M, size_t& N){
    
    std::ifstream file(path);
    if(!file){
        std::cerr << "Invalid file.";
        exit(1);
    }

    std::vector<std::string> grid;
    std::string line;
    //Loops through the file and pushes each line to the grid
    while(std::getline(file, line)){
        grid.push_back(line);
    }

    if(grid.empty()){
        std::cerr << "No file contents. Grid is empty.";
        exit(1);
    }
    //Gets # rows
    M = grid.size();
    //Gets # cols and checks that every column has equal rows, if not, fails
    N = grid[0].size();
    for(size_t i{0}; i < M; i++){
        if(grid[i].size() != N){
            std::cerr << "Column lengths vary. Failed.";
            exit(1);
        }
    }

    Panel panel(N, std::vector<uint8_t>(M));
    //Adding to panel vector by transposing
    for(size_t i{0}; i < M; i++){
        for(size_t j{0}; j < N; j++){
            //Subtracting by '0' converts the ASCII value to int
            panel[j][i] = grid[i][j] - '0';
        }
    }

    return panel;
}

int main(int argc, char** argv) {
    //taking inputs from cmd line
    if (argc < 2){
        //tells user how to properly use this (arguments to pass in)
        std::cerr << "usage: " << argv[0] << " <panel.txt>\n";
        return 1;
    }

    // M: Number of Haplotypes (rows)
    // N: Number of sites (columns)
    // Directed manipulated, as it is directly referenced in loadPanel
    size_t M, N;
    std::string path = argv[1];
    Panel panel = loadPanel(path, M, N);
    
    //Test block
    for(size_t i{0}; i < M; i++){
        std::string whole = "";
        for(size_t j{0}; j < N; j++){
            whole += panel[j][i] + '0';
        }
        std::cout << whole << "\n";
    }

    prefixArrs prefixArrays = buildPrefixArrays(panel, M, N);

    //Prefix test block
    // Loops N+1 times becuase it gets one initial array before any sites are processed, then N passes (one per site)
    for (size_t i{0}; i <= N; i++){
        for(size_t j{0}; j < M; j++){
            std::cout << prefixArrays[i][j] << " ";
        }
        std::cout << "\n";
    }
    
    return 0;
}