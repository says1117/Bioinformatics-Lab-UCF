#include <vector>
#include <cstdint>
#include <fstream>
#include <string>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <utility>

//Type alias for std::vector<std::vector<uint8_t>>
//which is the 2D-vector array that is being made
using Panel = std::vector<std::vector<uint8_t>>;

using prefixArrs = std::vector<std::vector<uint32_t>>;
using divergenceArrs = std::vector<std::vector<uint32_t>>;

struct PrefixAndDivergenceArrays{
    prefixArrs a;
    divergenceArrs div;
};

void ReportLongMatches(const Panel& panel, size_t M, size_t N, size_t L, const PrefixAndDivergenceArrays& output){

    // Each entry is (haplotype id, its index in output.a[k]).
    // The index is needed because the start of a maximal match between two
    // haplotypes is the largest divergence value between their two positions.
    std::vector<std::pair<uint32_t, size_t>> a, b;

    for(size_t k{0}; k <= N; k++){
        a.clear();
        b.clear();

        const std::vector<uint32_t>& div = output.div[k];

        //Reports one pair: start is max(div[lo+1 .. hi]), end is the site before k
        auto emit = [&](std::pair<uint32_t, size_t> u, std::pair<uint32_t, size_t> v){
            size_t lo = std::min(u.second, v.second);
            size_t hi = std::max(u.second, v.second);

            uint32_t start = 0;
            for(size_t t{lo + 1}; t <= hi; t++){
                if(div[t] > start){
                    start = div[t];
                }
            }

            std::cout << std::min(u.first, v.first) << '\t'
                      << std::max(u.first, v.first) << '\t'
                      << start << '\t' << (k - 1) << '\n';
        };

        //Lambda for flush to combat repetitive logic
        auto flush = [&](){
            if(k == N){
                //No site N to split on, so every surviving pair must be reported now
                for(size_t u{0}; u < a.size(); u++){
                    for(size_t v{u + 1}; v < a.size(); v++){
                        emit(a[u], a[v]);
                    }
                }
            }
            else{
                for(size_t u{0}; u < a.size(); u++){
                    for(size_t v{0}; v < b.size(); v++){
                        emit(a[u], b[v]);
                    }
                }
            }
        };

        //At k == N everything lands in a, so b is always empty there
        auto ready = [&](){
            return (k == N) ? (a.size() > 1) : (!a.empty() && !b.empty());
        };

        for(size_t i{0}; i < M; i++){
            if(div[i] + L > k){
                if(ready()){
                    flush();
                }
                a.clear();
                b.clear();
            }

            uint32_t hap = output.a[k][i];
            //k == N short-circuits before panel[k] is touched
            if(k == N || panel[k][hap] == 0){
                a.push_back({hap, i});
            }
            else{
                b.push_back({hap, i});
            }
        }

        if(ready()){
            flush();
        }
    }
}

PrefixAndDivergenceArrays buildPrefixArraysAndDivergenceArrays(const Panel& panel, size_t M, size_t N){
    
    PrefixAndDivergenceArrays output;
    output.a = prefixArrs(N+1, std::vector<uint32_t>(M));
    output.div = divergenceArrs(N+1, std::vector<uint32_t>(M));

    //needs comments
    for(size_t i{0}; i < M; i++){
        output.a[0][i] = i;
        output.div[0][i] = 0;
    }

    std::vector<uint32_t> zeroes, ones, d, e;
    for(size_t k{0}; k < N; k++){
        size_t p = k+1, q = k+1;

        //Resets sizes to 0 but keeps underlying buffer capacity
        zeroes.clear();
        ones.clear();
        d.clear();
        e.clear();

        for(size_t x{0}; x < M; x++){
            size_t hap = output.a[k][x];

            if(output.div[k][x] > p){
                p = output.div[k][x];
            }
            if(output.div[k][x] > q){
                q = output.div[k][x];
            }
            if(panel[k][hap] == 0){
                zeroes.push_back(hap);
                d.push_back(p);
                p = 0;
            }
            else{
                ones.push_back(hap);
                e.push_back(q);
                q = 0;
            }
        }
        size_t pos = 0;
        for(size_t i{0}; i < zeroes.size(); i++){
            output.a[k+1][pos] = zeroes[i];
            pos++;
        }
        for(size_t j{0}; j < ones.size(); j++){
            output.a[k+1][pos] = ones[j];
            pos++;
        }
        pos = 0;
                for(size_t i{0}; i < zeroes.size(); i++){
            output.div[k+1][pos] = d[i];
            pos++;
        }
        for(size_t j{0}; j < ones.size(); j++){
            output.div[k+1][pos] = e[j];
            pos++;
        }
    }

    return output;
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
    N = grid.size();
    //Gets # cols and checks that every column has equal rows, if not, fails
    M = grid[0].size();
    for(size_t k{0}; k < N; k++){
        if(grid[k].size() != M){
            std::cerr << "Row lengths vary. Failed.";
            exit(1);
        }
    }

    Panel panel(N, std::vector<uint8_t>(M));
    //Adding to panel vector by transposing
    //Grid rows are already sites and columns are already haplotypes,
    //so no transpose is needed
    for(size_t k{0}; k < N; k++){
        for(size_t h{0}; h < M; h++){
            //Subtracting by '0' converts the ASCII value to int
            panel[k][h] = grid[k][h] - '0';
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

    // Passing in L (minimum match length used for Algorithm 3) as the 3rd command-line argument for now, can be changed
    if (argc < 3){
        std::cerr << "Invalid Minimum Match Length. Needs to come after .txt file.";
        return 1;
    }

    // M: Number of Haplotypes (rows)
    // N: Number of sites (columns)
    // Directed manipulated, as it is directly referenced in loadPanel
    size_t M, N;
    std::string path = argv[1];
    //Converts Cmd-line argument string to integer
    size_t L = std::stoi(argv[2]);
    Panel panel = loadPanel(path, M, N);
    
    // //Test block
    // for(size_t i{0}; i < M; i++){
    //     std::string whole = "";
    //     for(size_t j{0}; j < N; j++){
    //         whole += panel[j][i] + '0';
    //     }
    //     std::cout << whole << "\n";
    //}

    //Algorithms 1 and 2, output is a variable that has direct access to the Prefix and Divergence arrays
    PrefixAndDivergenceArrays output = buildPrefixArraysAndDivergenceArrays(panel, M, N);

    // //Prefix test block
    // // Loops N+1 times becuase it gets one initial array before any sites are processed, then N passes (one per site)
    // for (size_t i{0}; i <= N; i++){
    //     for(size_t j{0}; j < M; j++){
    //         std::cout << output.a[i][j] << " ";
    //     }
    //     std::cout << "\n";
    // }
    // //Divergence test block
    // for (size_t i{0}; i <= N; i++){
    //     for(size_t j{0}; j < M; j++){
    //         std::cout << output.div[i][j] << " ";
    //     }
    //     std::cout << "\n";
    // }

    //Algorithm 3
    ReportLongMatches(panel, M, N, L, output);
    
    return 0;
}