#ifndef GRBINTERFACE_H
#define GRBINTERFACE_H
#include "gurobi_c++.h"
#include "KGraph.h"

string itos(int i);

vector<bool> boolify(vector<long> &S, long n);

/* Find a large k-club and then do preprocessing, as follows:
1. Find a large distance k-clique using k-th power graph G^k and clique heuristic in G^k
2. Use DROP to find a large k-club within the large k-clique. This gives a lower bound LB.
3. Use the LB to do preprocessing in a k-core "peeling" approach.
Here, k = LB-1 (and not the k from k-club).
*/
vector<long> HeuristicAndPreprocess(KGraph &g, long k);

// DROP heuristic for max k-club, due to Bourjolly et al. (2000).
vector<long> DROP(KGraph &g, long k);
vector<long> DROP(KGraph &g, vector<bool> &W, long k);


// returns the size of the k-hop neighborhood of vertex v in G (or in G[W]).
long KNeighborhoodSize(KGraph &g, long k, long v);
long KNeighborhoodSize(KGraph &g, vector<bool> &W, long k, long v);


// Our cut-like formulation for max k-club. Uses the routine MINIMALIZE to strengthen the cuts.
vector<long> solveMaxKClub_CutLike(KGraph &g, long k, vector<long> &HeuristicSolution, bool &subOpt);
vector<long> MINIMALIZE(KGraph &g, long a, long b, vector<long> length_k_ab_separator, long k);

// Combination of graph decomposition and our cut-like formulation
vector<long>  ICUT(KGraph &g, long k, vector<long> &BestKClub);

vector<long> ICUT_Subproblem(KGraph &g, long k, long v, long LowerBound, long &SubProblemLB, long &SubProblemUB);


// callback functions for cut-like formulation
class Kclub_callback : public GRBCallback
{
public:
	GRBVar *vars;
	KGraph g1;
	long k1;

	Kclub_callback(GRBVar *xvars, KGraph &g, long k)
	{
		vars = xvars;
		g1.Duplicate(g);
		k1 = k;
	}
	void callback();
	static long numCallbacks;
	static double TotalCallbackTime;
	static long numLazyCuts;
};


#endif
