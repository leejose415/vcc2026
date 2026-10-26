// Stream a Feng gene x cell UMI-count CSV.gz once and write, per output source, its counts as CSC column blocks
// (one column per gene, rows = the source's cells in a given order). Used by scripts/convert_new_contexts.py.
//
//   feng_csv_to_csc <RNA.csv.gz> <assign.bin> <n_sources> <outdir>
//
// assign.bin: int32 per RNA column (header order): the output source (0..n_sources-1) or -1 to drop the cell, then
// int32 per RNA column: the cell's row within its source. Every count must be a non-negative integer (exits otherwise).
// Writes <outdir>/src<k>.data (uint32), src<k>.rows (uint32), src<k>.colptr (uint64, n_genes + 1), and genes.txt.
#include <zlib.h>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;

struct Reader {
    gzFile f; vector<char> b; size_t p = 0, n = 0;
    explicit Reader(const char* path) : b(1 << 24) {
        f = gzopen(path, "rb"); if (!f) throw runtime_error("cannot open gzip");
        gzbuffer(f, 1 << 24);
    }
    inline int get() {
        if (p == n) { int k = gzread(f, b.data(), (unsigned)b.size()); if (k <= 0) return -1; n = (size_t)k; p = 0; }
        return (unsigned char)b[p++];
    }
};

int main(int argc, char** argv) {
    try {
        if (argc != 5) throw runtime_error("usage: feng_csv_to_csc RNA.csv.gz assign.bin n_sources outdir");
        const string out = argv[4]; const int S = stoi(argv[3]);
        Reader r(argv[1]);
        // header: first field is empty (row-name column), then one field per cell
        size_t N = 0; int c;
        while ((c = r.get()) != -1 && c != '\n') if (c == ',') N++;
        ifstream af(argv[2], ios::binary); vector<int32_t> src(N), row(N);
        af.read((char*)src.data(), N * 4); af.read((char*)row.data(), N * 4);
        if (!af) throw runtime_error("assign.bin shorter than 2 x RNA columns");
        vector<ofstream> fd(S), fr(S); vector<vector<uint64_t>> colptr(S, vector<uint64_t>(1, 0)); vector<uint64_t> nnz(S, 0);
        vector<vector<uint32_t>> bd(S), br(S);
        for (int k = 0; k < S; k++) {
            fd[k].open(out + "/src" + to_string(k) + ".data", ios::binary);
            fr[k].open(out + "/src" + to_string(k) + ".rows", ios::binary);
        }
        ofstream genes(out + "/genes.txt"); size_t g = 0;
        auto flush = [&](int k) {
            fd[k].write((char*)bd[k].data(), bd[k].size() * 4); fr[k].write((char*)br[k].data(), br[k].size() * 4);
            bd[k].clear(); br[k].clear();
        };
        while ((c = r.get()) != -1) {
            string feature;
            while (c != ',' && c != -1 && c != '\n') { if (c != '"' && c != '\r') feature += (char)c; c = r.get(); }
            if (c == -1) break;
            if (c != ',') throw runtime_error("missing separator after gene " + feature);
            genes << feature << '\n';
            for (size_t i = 0; i < N; i++) {
                uint64_t v = 0; bool any = false;
                while ((c = r.get()) != -1 && c != ',' && c != '\n' && c != '\r') {
                    if (c == '"') continue;
                    if (c < '0' || c > '9') throw runtime_error("non-integer count in gene " + feature);
                    v = v * 10 + (uint64_t)(c - '0'); any = true;
                }
                if (!any) throw runtime_error("empty field in gene " + feature);
                if (i + 1 < N && c != ',') throw runtime_error("short row in gene " + feature);
                if (i + 1 == N && c == ',') throw runtime_error("long row in gene " + feature);
                const int k = src[i];
                if (k >= 0 && v > 0) {
                    if (v > 0xFFFFFFFFull) throw runtime_error("count overflow");
                    bd[k].push_back((uint32_t)v); br[k].push_back((uint32_t)row[i]); nnz[k]++;
                    if (bd[k].size() >= (1u << 20)) flush(k);
                }
            }
            if (c == '\r') { c = r.get(); if (c != '\n') throw runtime_error("bad CRLF"); }
            for (int k = 0; k < S; k++) colptr[k].push_back(nnz[k]);
            if (++g % 2000 == 0) cerr << g << " genes\n";
        }
        for (int k = 0; k < S; k++) {
            flush(k); fd[k].close(); fr[k].close();
            ofstream fp(out + "/src" + to_string(k) + ".colptr", ios::binary);
            fp.write((char*)colptr[k].data(), colptr[k].size() * 8);
        }
        ofstream done(out + "/complete.txt"); done << g << " genes " << N << " cells\n";
        cerr << "done: " << g << " genes, " << N << " cells\n";
    } catch (exception& e) { cerr << "error: " << e.what() << '\n'; return 1; }
    return 0;
}
