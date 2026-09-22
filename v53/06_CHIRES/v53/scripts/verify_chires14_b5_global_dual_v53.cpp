#include <bits/stdc++.h>
using namespace std;
static inline void dec(uint32_t x, int a[9]){for(int i=0;i<9;i++) a[i]=(x>>(3*i))&7u;}
static inline uint32_t enc(const int a[9]){uint32_t x=0;for(int i=0;i<9;i++)x|=(uint32_t)(a[i]&7)<<(3*i);return x;}
static inline uint32_t mul8(uint32_t x,uint32_t y){int A[9],B[9],C[9]={0};dec(x,A);dec(y,B);for(int i=0;i<3;i++)for(int j=0;j<3;j++){int s=0;for(int k=0;k<3;k++)s+=A[3*i+k]*B[3*k+j];C[3*i+j]=s&7;}return enc(C);}
struct Term{uint8_t g;uint32_t p;};
struct Rel{int id;vector<Term> t;};
int main(int argc,char**argv){
 if(argc<5){cerr<<"usage: verify b5.bin fox.txt raw_solution.bin out.json\n";return 2;}
 string fb5=argv[1],fterms=argv[2],fsol=argv[3],fout=argv[4];
 ifstream b(fb5,ios::binary); b.seekg(0,ios::end); size_t bs=b.tellg(); b.seekg(0); if(bs%4){cerr<<"bad b5 size\n";return 3;} vector<uint32_t> H(bs/4); b.read((char*)H.data(),bs); b.close();
 ifstream tf(fterms); vector<Rel> R; string line; size_t termcnt=0; while(getline(tf,line)){if(line.empty())continue; stringstream ss(line); int id,w;ss>>id>>w;Rel r;r.id=id;for(int i=0;i<w;i++){string s;ss>>s;auto pos=s.find(':');if(pos==string::npos)return 4;int g=stoi(s.substr(0,pos));uint32_t p=(uint32_t)stoul(s.substr(pos+1));r.t.push_back({(uint8_t)g,p});} if((int)r.t.size()!=w)return 5;termcnt+=w;R.push_back(move(r));}tf.close();
 ifstream sf(fsol,ios::binary); uint32_t nsol=0;sf.read((char*)&nsol,4);vector<uint8_t>x(nsol);sf.read((char*)x.data(),nsol); if(!sf||sf.peek()!=EOF){cerr<<"bad sol\n";return 6;}sf.close();
 const uint64_t raw_eq=(uint64_t)H.size()*R.size();
 vector<vector<uint32_t>> eqs; eqs.reserve(raw_eq); vector<uint8_t> rhs;rhs.reserve(raw_eq);
 unordered_map<uint64_t,uint32_t> id; id.reserve(800000); id.max_load_factor(0.7);
 uint64_t inc=0, fail=0,target_fail=0,nontarget_fail=0,target_eq=0; uint64_t eidx=0;
 for(auto &r:R){for(uint32_t h:H){vector<uint32_t> q;q.reserve(r.t.size());uint8_t s=0; for(auto &t:r.t){uint32_t hp=mul8(h,t.p);uint64_t key=((uint64_t)t.g<<32)|hp;auto it=id.find(key);uint32_t v;if(it==id.end()){v=id.size();id.emplace(key,v);}else v=it->second;q.push_back(v);if(v>=x.size()){cerr<<"sol too short\n";return 7;}s^=x[v];}uint8_t rr=(r.id==42); if(rr)target_eq++; if(s!=rr){fail++;if(rr)target_fail++;else nontarget_fail++;}inc+=q.size();eqs.push_back(move(q));rhs.push_back(rr);eidx++;}}
 if(id.size()!=nsol){cerr<<"var mismatch "<<id.size()<<" vs "<<nsol<<"\n";return 8;}
 // Incidence and exact degree-one peeling to the 2-core.
 vector<vector<uint32_t>> ve(id.size());for(uint32_t e=0;e<eqs.size();e++)for(uint32_t v:eqs[e])ve[v].push_back(e);
 vector<int> deg(id.size());for(size_t v=0;v<ve.size();v++)deg[v]=ve[v].size(); vector<uint8_t> alive(eqs.size(),1);deque<uint32_t>Q;for(uint32_t v=0;v<deg.size();v++)if(deg[v]<=1)Q.push_back(v);
 uint64_t peeled=0; while(!Q.empty()){uint32_t v=Q.front();Q.pop_front();if(deg[v]>1)continue;if(deg[v]==0)continue;uint32_t e=UINT32_MAX;for(uint32_t a:ve[v])if(alive[a]){e=a;break;}if(e==UINT32_MAX){deg[v]=0;continue;}alive[e]=0;peeled++;for(uint32_t u:eqs[e]){if(deg[u]>0){deg[u]--;if(deg[u]<=1)Q.push_back(u);}}}
 uint64_t ce=0,ci=0,crhs=0,cv=0;int mind=INT_MAX,maxd=0;for(uint32_t e=0;e<alive.size();e++)if(alive[e]){ce++;ci+=eqs[e].size();crhs+=rhs[e];}for(int d:deg)if(d>0){cv++;mind=min(mind,d);maxd=max(maxd,d);}if(mind==INT_MAX)mind=0;
 uint64_t ones=0;for(auto z:x)ones+=z;
 bool ok = H.size()==22240 && R.size()==43 && termcnt==198 && raw_eq==956320 && id.size()==747196 && inc==4403520 && target_eq==22240 && fail==0 && ce==856952 && cv==465651 && ci==3528672 && crhs==44 && mind==2 && maxd==22 && peeled==99368 && ones==121219;
 ofstream o(fout);o<<"{\n"
 <<"  \"status\": \""<<(ok?"PASS":"FAIL")<<"\",\n"
 <<"  \"transport_order\": \"left multiplication h*p; relator-major then B5-class-major variable enumeration\",\n"
 <<"  \"b5_mod8_classes\": "<<H.size()<<",\n"
 <<"  \"relators\": "<<R.size()<<",\n"
 <<"  \"fox_terms_total\": "<<termcnt<<",\n"
 <<"  \"raw_equations\": "<<raw_eq<<",\n"
 <<"  \"raw_variables\": "<<id.size()<<",\n"
 <<"  \"raw_incidence\": "<<inc<<",\n"
 <<"  \"target_equations\": "<<target_eq<<",\n"
 <<"  \"raw_solution_weight\": "<<ones<<",\n"
 <<"  \"raw_verification_failures\": "<<fail<<",\n"
 <<"  \"target_verification_failures\": "<<target_fail<<",\n"
 <<"  \"nontarget_verification_failures\": "<<nontarget_fail<<",\n"
 <<"  \"degree_one_peeling\": {\n"
 <<"    \"peeled_equations\": "<<peeled<<",\n"
 <<"    \"core_equations\": "<<ce<<",\n"
 <<"    \"core_variables\": "<<cv<<",\n"
 <<"    \"core_incidence\": "<<ci<<",\n"
 <<"    \"core_rhs_one\": "<<crhs<<",\n"
 <<"    \"min_active_degree\": "<<mind<<",\n"
 <<"    \"max_active_degree\": "<<maxd<<"\n  },\n"
 <<"  \"decision\": \""<<(ok?"PROVEN-GLOBAL-MOD8-DUAL-OBSTRUCTION / CHIRES14-RADIUS<=5-REFUTED-TYPED":"NONPUBLISHABLE")<<"\"\n}\n";o.close();
 cerr<<"raw eq="<<raw_eq<<" vars="<<id.size()<<" inc="<<inc<<" target="<<target_eq<<" fail="<<fail<<"\n";
 cerr<<"core eq="<<ce<<" vars="<<cv<<" inc="<<ci<<" rhs1="<<crhs<<" mindeg="<<mind<<" maxdeg="<<maxd<<" peeled="<<peeled<<"\n";
 return ok?0:9;
}
