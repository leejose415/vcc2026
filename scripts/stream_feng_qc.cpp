// Stream gene-by-cell gzip CSV into RNA depth and sparse on-target/control evidence.
#include <zlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <cstdint>
using namespace std;
struct Reader{gzFile f;vector<char>b;size_t p=0,n=0;Reader(char*path):b(1<<22){f=gzopen(path,"rb");if(!f)throw runtime_error("open gzip");gzbuffer(f,1<<22);}int get(){if(p==n){int k=gzread(f,b.data(),b.size());if(k<0)throw runtime_error("gzip read error");n=k;p=0;if(!n)return -1;}return b[p++];}~Reader(){gzclose(f);}};
struct Entry{uint32_t row,col,value;};
template<class T>void save(string p,vector<T>&v){ofstream f(p,ios::binary);f.write((char*)v.data(),v.size()*sizeof(T));}
int main(int argc,char**argv){try{if(argc!=4)throw runtime_error("RNA.gz setup.tsv outdir");string out=argv[3];ifstream setup(argv[2]);vector<int>target,ctrl;unordered_map<string,int>code;vector<string>names;string line;
while(getline(setup,line)){auto k=line.find('\t');string g=line.substr(0,k);if(!code.count(g)){code[g]=names.size();names.push_back(g);}target.push_back(code[g]);ctrl.push_back(stoi(line.substr(k+1)));}
size_t N=target.size();vector<double>tot(N),mito(N);vector<uint32_t>ng(N),measured(names.size());ofstream evidence(out+"/evidence.bin",ios::binary);Reader r(argv[1]);int c;size_t header=0;while((c=r.get())!=-1&&c!='\n')if(c==',')header++;if(header!=N)throw runtime_error("header count mismatch");int nr=0;
while((c=r.get())!=-1){string feature;while(c!=','&&c!=-1&&c!='\n'){if(c!='\"'&&c!='\r')feature+=char(c);c=r.get();}if(c==-1)break;if(c!=',')throw runtime_error("missing gene separator");auto a=feature.find(':');string gene=a==string::npos?feature:feature.substr(a+1,feature.find(':',a+1)-a-1);auto it=code.find(gene);int gc=it==code.end()?-1:it->second;if(gc>=0)measured[gc]=1;bool mt=gene.rfind("MT-",0)==0;
for(uint32_t i=0;i<N;i++){uint64_t v=0;bool digit=false,decimal=false;while((c=r.get())!=-1&&c!=','&&c!='\n'&&c!='\r'){if(c=='\"')continue;if(c=='.'&&!decimal){decimal=true;continue;}if(c<'0'||c>'9')throw runtime_error("invalid/negative RNA count");digit=true;if(decimal){if(c!='0')throw runtime_error("noninteger RNA count");}else v=v*10+(c-'0');if(v>UINT32_MAX)throw runtime_error("count overflow");}if(!digit)throw runtime_error("empty RNA count");if(v){tot[i]+=v;ng[i]++;if(mt)mito[i]+=v;if(gc>=0&&(ctrl[i]||target[i]==gc)){Entry e{i,(uint32_t)gc,(uint32_t)v};evidence.write((char*)&e,sizeof(e));}}
if(i+1<N&&c!=',')throw runtime_error("short row");if(i+1==N&&c==',')throw runtime_error("long row");}
if(c=='\r'){c=r.get();if(c!='\n')throw runtime_error("bad CRLF");}nr++;if(nr%2000==0)cerr<<nr<<" genes, "<<N<<" cells\n";}
save(out+"/totals.bin",tot);save(out+"/ngenes.bin",ng);save(out+"/mito.bin",mito);save(out+"/measured.bin",measured);ofstream g(out+"/genes.txt");for(auto&s:names)g<<s<<'\n';ofstream done(out+"/complete.txt");done<<N<<" "<<nr<<"\n";cerr<<"COMPLETE "<<N<<" cells "<<nr<<" genes\n";
}catch(exception&e){cerr<<e.what()<<'\n';return 1;}}
