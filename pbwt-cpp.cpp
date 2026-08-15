#include <vector>
#include <cstdint>
#include <fstream>
#include <string>
#include <iostream>
#include <cassert>

//Type alias for std::vector<std::vector<uint8_t>>
//which is the 2D-vector array that is being made
using Panel = std::vector<std::vector<uint8_t>>;

//Loads 2D-array
Panel loadPanel(std::string path, size_t& M, size_t& N){

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
    
    return 0;
}