#include <vector>
#include <cstdint>
#include <fstream>
#include <string>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <utility>
#include <map>
#include <functional>

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

// Algorithm 4: ReportSetMaximalMatches -- greedy Kruskal-style consumption over the divergence-array gaps, see chat for full writeup with examples
void ReportSetMaximalMatches(const Panel& panel, size_t M, size_t N, const PrefixAndDivergenceArrays& output){

    std::vector<size_t> parent(M);
    std::vector<uint8_t> allele(M);   // valid only while a set is still open
    std::vector<bool> closed(M);      // true once a set has been matched away
    std::vector<uint32_t> groupSize(M);

    // iterative path-compressing find
    auto find = [&](size_t x){
        while(parent[x] != x){
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };

    for(size_t k{0}; k <= N; k++){
        const std::vector<uint32_t>& a = output.a[k];
        const std::vector<uint32_t>& div = output.div[k];

        if(k == N){
            // k == N: no site to compare alleles against, so no consumption -- independent per-row tie scan instead
            for(size_t i{0}; i < M; i++){
                uint32_t bestUp = UINT32_MAX, bestDown = UINT32_MAX;
                std::vector<size_t> upTies, downTies;

                if(i > 0){
                    uint32_t running = div[i];
                    size_t r = i - 1;
                    bestUp = running;
                    upTies.push_back(r);
                    while(r > 0 && div[r] <= bestUp){
                        running = std::max(running, div[r]);
                        if(running > bestUp) break;
                        r--;
                        upTies.push_back(r);
                    }
                }
                if(i + 1 < M){
                    uint32_t running = div[i + 1];
                    size_t r = i + 1;
                    bestDown = running;
                    downTies.push_back(r);
                    while(r + 1 < M && div[r + 1] <= bestDown){
                        running = std::max(running, div[r + 1]);
                        if(running > bestDown) break;
                        r++;
                        downTies.push_back(r);
                    }
                }

                uint32_t best = std::min(bestUp, bestDown);
                if(best == UINT32_MAX) continue;
                uint32_t length = static_cast<uint32_t>(k) - best;
                if(length == 0) continue;
                uint32_t end1Based = static_cast<uint32_t>(k) + 1;

                if(bestUp == best) for(size_t r : upTies) std::cout << a[i] << '\t' << a[r] << '\t' << end1Based << '\t' << length << '\n';
                if(bestDown == best) for(size_t r : downTies) std::cout << a[i] << '\t' << a[r] << '\t' << end1Based << '\t' << length << '\n';
            }
            continue;
        }

        for(size_t i{0}; i < M; i++){
            parent[i] = i;
            closed[i] = false;
            groupSize[i] = 1;
            allele[i] = panel[k][a[i]];
        }

        // gaps = div[1..M-1]; sorted and processed tier-by-tier (equal-weight gaps batched together)
        std::vector<size_t> gapOrder(M > 0 ? M - 1 : 0);
        for(size_t i{0}; i < gapOrder.size(); i++) gapOrder[i] = i + 1;
        std::sort(gapOrder.begin(), gapOrder.end(), [&](size_t x, size_t y){ return div[x] < div[y]; });

        size_t batchStart = 0;
        while(batchStart < gapOrder.size()){
            size_t batchEnd = batchStart;
            while(batchEnd < gapOrder.size() && div[gapOrder[batchEnd]] == div[gapOrder[batchStart]]) batchEnd++;

            // A: same-allele merges within this tier
            for(size_t bi{batchStart}; bi < batchEnd; bi++){
                size_t gi = gapOrder[bi];
                size_t L = find(gi - 1), R = find(gi);
                if(L == R) continue;
                if(closed[L] || closed[R]) continue; // wall from a strictly earlier (tighter) tier
                if(allele[L] == allele[R]){
                    parent[R] = L;
                    groupSize[L] += groupSize[R];
                }
            }

            // B: differing-allele events, grouped into chains (see chat for the chain-length rules)
            uint32_t end1Based = static_cast<uint32_t>(k) + 1;
            uint32_t length = static_cast<uint32_t>(k) - div[gapOrder[batchStart]];

            std::vector<std::pair<size_t,size_t>> candidates; // (L, R) roots
            std::vector<size_t> touchCount(M, 0);
            std::vector<size_t> toClose;
            for(size_t bi{batchStart}; bi < batchEnd; bi++){
                size_t gi = gapOrder[bi];
                size_t L = find(gi - 1), R = find(gi);
                if(L == R) continue;                  // merged away in Phase A
                if(closed[L] && closed[R]) continue;   // both already gone, nothing to do
                if(closed[L] || closed[R]){
                    // one side already used up elsewhere -- the open side is spent too, no report
                    toClose.push_back(closed[L] ? R : L);
                    continue;
                }
                candidates.push_back({L, R});
                touchCount[L]++;
                touchCount[R]++;
            }

            // scratch union-find, just to cluster this tier's edges into chains
            std::vector<size_t> chainParent(M);
            for(size_t i{0}; i < M; i++) chainParent[i] = i;
            std::function<size_t(size_t)> chainFind = [&](size_t x){
                while(chainParent[x] != x){ chainParent[x] = chainParent[chainParent[x]]; x = chainParent[x]; }
                return x;
            };
            for(auto [L, R] : candidates) chainParent[chainFind(L)] = chainFind(R);

            std::map<size_t, std::vector<std::pair<size_t,size_t>>> chains;
            for(auto& e : candidates) chains[chainFind(e.first)].push_back(e);

            for(auto& [chainRoot, edges] : chains){
                if(edges.size() == 1){
                    auto [L, R] = edges[0];
                    if(length > 0){
                        if(groupSize[L] == 1) for(size_t y{0}; y < M; y++) if(find(y) == R) std::cout << a[L] << '\t' << a[y] << '\t' << end1Based << '\t' << length << '\n';
                        if(groupSize[R] == 1) for(size_t x{0}; x < M; x++) if(find(x) == L) std::cout << a[R] << '\t' << a[x] << '\t' << end1Based << '\t' << length << '\n';
                    }
                    toClose.push_back(L);
                    toClose.push_back(R);
                }
                else if(edges.size() == 2){
                    // Exactly one of the three groups touches both edges -- that's the hub.
                    size_t hub = (touchCount[edges[0].first] == 2) ? edges[0].first : edges[0].second;
                    if(length > 0 && groupSize[hub] == 1){
                        for(auto [L, R] : edges){
                            size_t other = (L == hub) ? R : L;
                            for(size_t y{0}; y < M; y++) if(find(y) == other) std::cout << a[hub] << '\t' << a[y] << '\t' << end1Based << '\t' << length << '\n';
                        }
                    }
                    for(auto [L, R] : edges){ toClose.push_back(L); toClose.push_back(R); }
                }
                else{
                    // 3+ edges: unresolvable, no report, but still closes (can't retry at a looser tier)
                    for(auto [L, R] : edges){ toClose.push_back(L); toClose.push_back(R); }
                }
            }
            for(size_t r : toClose) closed[r] = true;

            batchStart = batchEnd;
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

// Loads a phased VCF; sample column s contributes haplotype ids 2s, 2s+1
Panel loadVCF(std::string path, size_t& M, size_t& N){

    std::ifstream file(path);
    if(!file){
        std::cerr << "Invalid file.";
        exit(1);
    }

    std::string line;
    size_t numSamples = 0;

    // skip ## metadata, read #CHROM header to count sample columns (9 fixed columns before them)
    while(std::getline(file, line)){
        if(line.rfind("##", 0) == 0) continue;
        if(line.rfind("#CHROM", 0) == 0){
            size_t tabs = 0;
            for(char c : line) if(c == '\t') tabs++;
            numSamples = (tabs >= 8) ? (tabs - 8) : 0;
            break;
        }
        std::cerr << "Malformed VCF: expected ## metadata then a #CHROM header.";
        exit(1);
    }
    if(numSamples == 0){
        std::cerr << "No sample columns found in VCF header.";
        exit(1);
    }

    M = numSamples * 2;

    std::vector<std::vector<uint8_t>> rows; // one entry per site, each of size M

    while(std::getline(file, line)){
        if(line.empty()) continue; // tolerate a trailing blank line

        std::vector<std::string> cols;
        size_t start = 0;
        for(size_t p{0}; p <= line.size(); p++){
            if(p == line.size() || line[p] == '\t'){
                cols.push_back(line.substr(start, p - start));
                start = p + 1;
            }
        }
        if(cols.size() < 9 + numSamples){
            std::cerr << "Malformed VCF data line (too few columns).";
            exit(1);
        }

        std::vector<uint8_t> site(M);
        for(size_t s{0}; s < numSamples; s++){
            const std::string& gt = cols[9 + s];
            // Expect "a|b" (phased); split on the '|' separator.
            size_t bar = gt.find('|');
            if(bar == std::string::npos){
                std::cerr << "Expected phased genotype (a|b), got: " << gt;
                exit(1);
            }
            site[2 * s]     = gt[0] - '0';
            site[2 * s + 1] = gt[bar + 1] - '0';
        }
        rows.push_back(std::move(site));
    }

    N = rows.size();
    Panel panel(N, std::vector<uint8_t>(M));
    for(size_t k{0}; k < N; k++){
        panel[k] = rows[k];
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

    // optional 4th arg: "3" (default) or "4" to pick the algorithm
    std::string algo = (argc >= 4) ? argv[3] : "3";

    // M: Number of Haplotypes (rows)
    // N: Number of sites (columns)
    // Directed manipulated, as it is directly referenced in loadPanel/loadVCF
    size_t M, N;
    std::string path = argv[1];
    //Converts Cmd-line argument string to integer
    size_t L = std::stoi(argv[2]);

    // dispatch on file extension: .vcf vs the plain 0/1-per-line panel format
    bool isVCF = path.size() >= 4 && path.compare(path.size() - 4, 4, ".vcf") == 0;
    Panel panel = isVCF ? loadVCF(path, M, N) : loadPanel(path, M, N);

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

    if(algo == "4"){
        //Algorithm 4
        ReportSetMaximalMatches(panel, M, N, output);
    }
    else{
        //Algorithm 3
        ReportLongMatches(panel, M, N, L, output);
    }

    return 0;
}