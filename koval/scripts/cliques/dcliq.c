/* count cliques of size W (Delsarte bound) in SRGs given as graph6, one per line.
   Usage: dcliq W E [maxclique_seconds] < file.g6
   W = clique size to count, E = number of neighbors every outside vertex has in a W-clique
   (Delsarte clique property, used as an exact shortcut: once |S| > E, every common neighbor
   of S must lie in the clique). E < 0 disables the shortcut (plain enumeration, for cross-checks).
   Output per graph: "G idx n count nodes ms [maxclique]" and "C idx v1 v2 ..." per clique.     */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#define MAXN 512
#define WD (MAXN/64)
typedef struct { uint64_t w[WD]; } bs;
static int n, nw;
static bs adj[MAXN];
static int W, E;
static long long nodes, count;
static int S[MAXN], dS;
static int gidx;
static FILE *out;
static inline int pc(const bs *a){int c=0;for(int i=0;i<nw;i++)c+=__builtin_popcountll(a->w[i]);return c;}
static inline int empty(const bs *a){for(int i=0;i<nw;i++)if(a->w[i])return 0;return 1;}
static inline void band(bs *r,const bs *a,const bs *b){for(int i=0;i<nw;i++)r->w[i]=a->w[i]&b->w[i];}
static inline int subset(const bs *a,const bs *b){for(int i=0;i<nw;i++)if(a->w[i]&~b->w[i])return 0;return 1;}
static inline int first(const bs *a){for(int i=0;i<nw;i++)if(a->w[i])return i*64+__builtin_ctzll(a->w[i]);return -1;}
static inline void clr(bs *a,int v){a->w[v>>6]&=~(1ULL<<(v&63));}
static inline int has(const bs *a,int v){return (a->w[v>>6]>>(v&63))&1;}
static int isclique(const bs *Q){int v[MAXN],m=0;for(int i=0;i<n;i++)if(has(Q,i))v[m++]=i;
  for(int i=0;i<m;i++)for(int j=i+1;j<m;j++)if(!has(&adj[v[i]],v[j]))return 0;return 1;}
static void emit(const bs *Q){ fprintf(out,"C %d",gidx); for(int i=0;i<dS;i++)fprintf(out," %d",S[i]);
  if(Q)for(int i=0;i<n;i++)if(has(Q,i))fprintf(out," %d",i); fprintf(out,"\n"); }
/* greedy coloring of P: order[] and col[] (col nondecreasing) */
static int colour(const bs *P,int *order,int *col){
  bs U=*P; int m=0,c=0;
  while(!empty(&U)){
    c++; bs Q=U; int v;
    while((v=first(&Q))>=0){
      clr(&U,v); order[m]=v; col[m]=c; m++;
      for(int i=0;i<nw;i++) Q.w[i]&=~adj[v].w[i];
      clr(&Q,v);
    }
  }
  return m;
}
static void expand(bs P, bs Qfull){
  nodes++;
  if(dS==W){ count++; emit(NULL); return; }
  if(E>=0 && dS>E){ /* Delsarte shortcut: C = S u Qfull exactly */
    if(pc(&Qfull)!=W-dS) return;
    if(!subset(&Qfull,&P)) return;
    if(!isclique(&Qfull)) return;
    count++; emit(&Qfull); return; }
  int order[MAXN], col[MAXN];
  int m=colour(&P,order,col);
  for(int i=m-1;i>=0;i--){
    if(dS+col[i]<W) return;
    int v=order[i];
    bs NP, NQ; band(&NP,&P,&adj[v]); band(&NQ,&Qfull,&adj[v]);
    S[dS++]=v; expand(NP,NQ); dS--;
    clr(&P,v);
  }
}
/* max clique (MCQ), with deadline */
static int best; static double deadline; static int timed_out;
static double now(){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
static void mc(bs P,int d){
  if(timed_out) return;
  if((nodes++ & 4095)==0 && now()>deadline){timed_out=1;return;}
  if(empty(&P)){ if(d>best)best=d; return; }
  int order[MAXN], col[MAXN]; int m=colour(&P,order,col);
  for(int i=m-1;i>=0;i--){ if(d+col[i]<=best) return; int v=order[i];
    bs NP; band(&NP,&P,&adj[v]); mc(NP,d+1); clr(&P,v); if(timed_out)return; }
}
static int readg6(const char *s){
  const unsigned char *p=(const unsigned char*)s; if(*p=='>')p+=10;
  if(p[0]==126){ n=((p[1]-63)<<12)|((p[2]-63)<<6)|(p[3]-63); p+=4;} else { n=p[0]-63; p++; }
  if(n>MAXN) return 0; nw=(n+63)/64; memset(adj,0,sizeof(adj));
  int k=0, bit=6, val=0; p--; 
  for(int j=1;j<n;j++)for(int i=0;i<j;i++){
    if(bit==6){ p++; val=*p-63; bit=0; }
    if((val>>(5-bit))&1){ adj[i].w[j>>6]|=1ULL<<(j&63); adj[j].w[i>>6]|=1ULL<<(i&63);} bit++; k++; }
  return 1;
}
/* returns "k,lambda,mu" if strongly regular, else "notSRG" */
static void srgcheck(char *buf){
  int k=pc(&adj[0]), lam=-1, mu=-1;
  for(int i=0;i<n;i++){ if(pc(&adj[i])!=k||has(&adj[i],i)){strcpy(buf,"notSRG");return;}
    for(int j=i+1;j<n;j++){ bs t; band(&t,&adj[i],&adj[j]); int c=pc(&t);
      if(has(&adj[i],j)){ if(lam<0)lam=c; else if(c!=lam){strcpy(buf,"notSRG");return;} }
      else { if(mu<0)mu=c; else if(c!=mu){strcpy(buf,"notSRG");return;} } } }
  sprintf(buf,"srg(%d,%d,%d,%d)",n,k,lam,mu);
}
int main(int argc,char**argv){
  if(argc<3){fprintf(stderr,"usage\n");return 1;}
  W=atoi(argv[1]); E=atoi(argv[2]); double mcs = argc>3? atof(argv[3]) : 0;
  out=stdout; static char line[1<<20]; gidx=0;
  while(fgets(line,sizeof line,stdin)){
    size_t L=strlen(line); while(L&&(line[L-1]=='\n'||line[L-1]=='\r'))line[--L]=0; if(!L)continue;
    if(!readg6(line)){fprintf(stderr,"too big\n");return 1;}
    double t0=now(); nodes=0; count=0; dS=0;
    bs P; memset(&P,0,sizeof P); for(int i=0;i<n;i++)P.w[i>>6]|=1ULL<<(i&63);
    expand(P,P);
    double t1=now(); long long cn=nodes;
    int mcv=-1; if(mcs>0){ best=0; timed_out=0; nodes=0; deadline=now()+mcs; mc(P,0); mcv= timed_out? -best : best; }
    char sb[64]; srgcheck(sb);
    fprintf(out,"G %d %d %lld %lld %.0f %d %.0f %s\n",gidx,n,count,cn,(t1-t0)*1000,mcv,(now()-t1)*1000,sb);
    fflush(out); gidx++;
  }
  return 0;
}
