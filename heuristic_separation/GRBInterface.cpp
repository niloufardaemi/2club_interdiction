#include "GRBInterface.h"
#include <sstream>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <queue>
#include <limits.h>
#include <string.h>

using namespace std;

string itos(int i) { stringstream s; s << i; return s.str(); }

vector<bool> boolify(vector<long> &S, long n)
{
	vector<bool> Sbool(n, false);
	for (long i = 0; i < S.size(); i++) Sbool[S[i]] = true;
	return Sbool;
}

vector<long> HeuristicAndPreprocess(KGraph &g, long k)
{
	// find k-th power of graph g and call it gk
	KGraph gk = g.CreatePowerGraph(k);
	vector<long> rd;  // right-degree of degeneracy ordering.
	vector<long> degeneracyordering = gk.FindDegeneracyOrdering(rd);
	vector<long> kclique = gk.FindHeuristicClique(degeneracyordering, rd);
	cerr << "k-clique (heuristic) size is = " << kclique.size() << endl;

	// perform DROP heuristic, with kclique as starting point.
	vector<bool> HeuristicSolution = boolify(kclique, g.n); // make the heuristic kclique into a boolean form
	vector<long> DROP_Solution = DROP(g, HeuristicSolution, k);
	long lb = DROP_Solution.size();
	cerr << "Drop heuristic gives LB = " << lb << endl;

	// perform preprocessing
	vector<long> kcorevertices = gk.FindVerticesOfKCore(degeneracyordering, rd, lb - 1);
	cerr << "Preprocessed instances has this many vertices = " << kcorevertices.size() << endl;
	g.FindInducedGraph(kcorevertices);

	return DROP_Solution;
}

vector<long> DROP(KGraph &g, long k)
{
	vector<bool> W(g.n, true);	// our k-club (eventually)
	return DROP(g, W, k);
}
vector<long> DROP(KGraph &g, vector<bool> &W, long k)
{
	vector<long> near(g.n, 0);	// the number of vertices "nearby" to a vertex (within k-hops)

	long Wsize = count(W.begin(), W.end(), true);
	// while W is not an k-club, delete some "bad" vertex. The number of vertices in W is size.
	for (long size = Wsize; size >= 0; size--)
	{
		// compute how many vertices are "nearby" to each vertex.
		for (long i = 0; i < g.n; i++)
		{
			if (!W[i]) continue;
			near[i] = KNeighborhoodSize(g, W, k, i);
		}

		// pick vertex w (in W) with smallest "nearby" vertices
		long w;	// vertex with smallest number of nearby vertices
		long smallestNearby = size; // number of vertices nearby to w.
		for (long i = 0; i < g.n; i++)
		{
			if (!W[i]) continue;
			if (near[i] < smallestNearby)
			{
				w = i;
				smallestNearby = near[i];
			}
		}

		// check if k-club. If not, then remove vertex w.
		if (smallestNearby == size) break;
		W[w] = false;
	}

	// convert W to vector<long> form, and return
	vector<long> Wvector;
	for (long i = 0; i < g.n; i++)
	{
		if (W[i])
		{
			Wvector.push_back(i);
		}
	}
	return Wvector;
}

long Kclub_callback::numCallbacks = 0;
double Kclub_callback::TotalCallbackTime = 0;
long Kclub_callback::numLazyCuts = 0;



// calback function for our main method - using MINIMALIZE 
// to obtain minimal subset of length-k a,b-separator
void Kclub_callback::callback()
{
	try {
		if (where == GRB_CB_MIPSOL)
		{
			numCallbacks++;
			time_t start = clock();

			// get the ''solution (S)'' from Gurobi
			double *x = new double[g1.n];
			x = getSolution(vars, g1.n);

			// make it boolean and call it D
			vector<bool> D(g1.n, false);

			for (long i = 0; i < g1.n; i++)
			{
				if (x[i] > 0.5) D[i] = true;
			}

			if (count(D.begin(), D.end(), false) != g1.n)
			{
				vector<long> SelectedVertices;
				vector<long> C_Prime;

				for (long i = 0; i < g1.n; i++)
				{
					if (D[i]) SelectedVertices.push_back(i);
					else C_Prime.push_back(i);
				}

				// create G[D], which we call g2
				vector<long> Rmap;
				KGraph g2 = g1.CreateInducedGraph(D, Rmap);
				for (long i = 0; i < g2.n; i++)
				{
					vector<long> dist_from_i_to = g2.ShortestPathsUnweighted(i);
					for (long j = i + 1; j < g2.n; j++)
					{
						if (dist_from_i_to[j] > k1)
						{
							long a = SelectedVertices[i];
							long b = SelectedVertices[j];

							// C_Prime is a length-s a,b-separator. now make it minimal
							vector<long> minimal_length_k_ab_separator = MINIMALIZE(g1, a, b, C_Prime, k1);
							GRBLinExpr expr = vars[a] + vars[b];
							for (long q = 0; q < minimal_length_k_ab_separator.size(); q++)
							{
								long v = minimal_length_k_ab_separator[q];
								expr -= vars[v];
							}

							addLazy(expr <= 1);
							numLazyCuts++;
						}
					}
				}
			}
			TotalCallbackTime += (double)(clock() - start) / CLOCKS_PER_SEC;
			delete[] x;
		}
	}
	catch (GRBException e) {
		cout << "Error number: " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during callback" << endl;
	}
}


vector<long> MINIMALIZE(KGraph &g, long a, long b, vector<long> ab_Separator, long k)
{
	vector<long> MinimalCut; // what we return at end of function
	vector<bool> Cut(g.n, false); // a boolean representation of the cut. 

								  // first do some linear-time preprocessing to remove vertices that are not on length-k a,b-path
								  //		for example, this removes vertices that belong to a different component in G.
	vector<long> dist_from_a = g.ShortestPathsUnweighted(a);
	vector<long> dist_from_b = g.ShortestPathsUnweighted(b);
	for (long i = 0; i < ab_Separator.size(); i++)
	{
		long v = ab_Separator[i];
		if (dist_from_a[v] + dist_from_b[v] <= k)
		{
			Cut[v] = true;  // vertex v was in ab_Separator AND belongs to a length-k a,b-path in graph G
		}
	}

	// initialize VerticesInGc = V \ Cut
	vector<bool> VerticesInGc(g.n, true);
	for (long i = 0; i < g.n; i++) if (Cut[i]) VerticesInGc[i] = false;

	// now run the normal minimalize algorithm.
	for (long c = 0; c < g.n; c++)
	{
		if (!Cut[c]) continue; // only test for cut vertices
		if (dist_from_a[c] == 1 && dist_from_b[c] == 1) continue; // if c \in N(a)\cap N(b), then c belongs to every minimal cut.

																  // put vertex c in G_c
		VerticesInGc[c] = true;

		// compute distances from c in G_c
		vector<long> dist_from_c = g.ShortestPathsUnweighted(c, VerticesInGc);

		// check if vertex c would close the distance from a to b to be at most k.
		if (dist_from_c[a] + dist_from_c[b] <= k) // vertex c remains in cut
		{
			VerticesInGc[c] = false;
		}
		else // vertex c is pulled from the cut.
		{
			Cut[c] = false;
		}
	}
	for (long i = 0; i < g.n; i++) if (Cut[i]) MinimalCut.push_back(i);
	return MinimalCut;
}


long KNeighborhoodSize(KGraph &g, vector<bool> &W, long k, long v)
{
	if (!W[v])
	{
		cerr << "\n ERROR: K-neighborhoodSize. You are calculating distances across W nodes starting from some vertex v which does not belong to W.";
		return 0;
	}
	vector<long> dist = g.ShortestPathsUnweighted(v, W);
	long nsize = 0;		// the number of nodes j with dist(v,j) <= s (including node v).
	for (long i = 0; i < g.n; i++) if (dist[i] <= k) nsize++;
	return nsize;
}
long KNeighborhoodSize(KGraph &g, long k, long v)
{
	vector<bool> W(g.n, true);
	return KNeighborhoodSize(g, W, k, v);
}

vector<long> solve2Club(KGraph &g, long k, vector<long> &HeuristicSolution, bool &subOpt)
{

	vector<long> S;
	subOpt = true;
	vector< vector< long> > components;
	vector<long> degreeZero;
	g.FindConnectedComponents(components, degreeZero);
	long lb = HeuristicSolution.size();

	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_IntParam_Method, 3); //concurrent  
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar *Y = model.addVars(components.size(), GRB_BINARY);
		model.update();


		// Adding objective function
		for (long w = 0; w < g.n; w++)
			X[w].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);


		// Fixing singletons to zero.
		for (long i = 0; i < degreeZero.size(); i++)
		{
			long v = degreeZero[i];
			model.addConstr(X[v] == 0);
		}

		// Fixing y[i]=0 when |V(G_i)| < lb.
		for (long i = 0; i < components.size(); i++)
		{
			if (components[i].size() < lb)
			{
				model.addConstr(Y[i] == 0);
			}
		}

		// Adding X[v] + X[q] - X[N(a) \cap N(b)] <= 1 
		for (long i = 0; i < components.size(); i++)
		{
			// add all constraints within component i
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j]; // v is a vertex in component i.
				vector<bool> neighbors(g.n, false);
				long w;
				for (long q = 0; q < g.degree[v]; q++)
				{
					w = g.adj[v][q];
					neighbors[w] = true;
				}

				for (long q = v + 1; q < g.n; q++)
				{
					if (neighbors[q]) continue;
					vector<long> p = g.CommonNeighborsList(v, q);

					GRBLinExpr expr = 0;
					for (long k = 0; k < p.size(); k++)
					{
						expr += X[p[k]];
					}
					expr = X[v] + X[q] - expr;
					model.addConstr(expr <= 1);
				}
			}
		}


		// Adding \sigma_i y_i <= 1 constraints.
		GRBLinExpr expr = 0;
		for (long i = 0; i < components.size(); i++)
		{
			expr += Y[i];
		}
		model.addConstr(expr <= 1);



		// Adding X[v] <= Y[i] constraints.
		for (long i = 0; i < components.size(); i++)
		{
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j];  // v is a vertex in component i
				model.addConstr(X[v] <= Y[i]);
			}
		}
		model.update();

		// Providing initial solution
		for (long i = 0; i < g.n; i++)
			X[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < HeuristicSolution.size(); i++)
		{
			long v = HeuristicSolution[i];
			X[v].set(GRB_DoubleAttr_Start, 1);
		}

		model.optimize();

		double bestLB = model.get(GRB_DoubleAttr_ObjVal);
		double bestUB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subOpt = false;
			for (long l = 0; l < g.n; l++)
				if (X[l].get(GRB_DoubleAttr_X) > 0.5)
					S.push_back(l);
		}
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return S;
}

vector<long> solveMaxKClub_CutLike(KGraph &g, long k, vector<long> &HeuristicSolution, bool &subOpt)
{
	vector<long> S;
	subOpt = true;
	vector< vector< long> > components;
	vector<long> degreeZero;
	g.FindConnectedComponents(components, degreeZero);
	long lb = HeuristicSolution.size();

	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_Method, 3); //concurrent 
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar *Y = model.addVars(components.size(), GRB_BINARY);
		model.update();


		// Adding objective function
		for (long w = 0; w < g.n; w++)
			X[w].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		model.update();

		// Fixing singletons to zero.
		for (long i = 0; i < degreeZero.size(); i++)
		{
			long v = degreeZero[i];
			model.addConstr(X[v] == 0);
		}

		// Fixing y[i]=0 when |V(G_i)| < lb.
		for (long i = 0; i < components.size(); i++)
		{
			if (components[i].size() < lb)
			{
				model.addConstr(Y[i] == 0);
			}
		}

		// Adding x_u + x_v <= 1 constraints
		// only need to add such constraints when u and v belong to same component.
		for (long i = 0; i < components.size(); i++)
		{
			// add all constraints within component i
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j]; // v is a vertex in component i.
				vector<long> dist = g.ShortestPathsUnweighted(v);
				for (long q = j + 1; q < components[i].size(); q++)
				{
					long u = components[i][q];  // u is another vertex in component i
												// vertices u and v belong to same component i
												// if dist(v,u) > s, then add constraint x_u + x_v <= 1.
					if (dist[u] > k)
					{
						model.addConstr(X[u] + X[v] <= 1);
					}
				}
			}
		}
		model.update();


		// Adding \sigma_i y_i <= 1 constraints.
		GRBLinExpr expr = 0;
		for (long i = 0; i < components.size(); i++)
		{
			expr += Y[i];
		}
		model.addConstr(expr <= 1);
		model.update();


		// Adding X[v] <= Y[i] constraints.
		for (long i = 0; i < components.size(); i++)
		{
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j];  // v is a vertex in component i
				model.addConstr(X[v] <= Y[i]);
			}
		}
		model.update();


		//Adding lazy constraints.
		Kclub_callback cb = Kclub_callback(X, g, k);
		model.setCallback(&cb);

		// Providing initial solution
		for (long i = 0; i < g.n; i++)
			X[i].set(GRB_DoubleAttr_Start, 0);

		for (long i = 0; i < HeuristicSolution.size(); i++)
		{
			long v = HeuristicSolution[i];
			X[v].set(GRB_DoubleAttr_Start, 1);
		}

		model.optimize();

		long NumOfBandBNodes = (long)model.get(GRB_DoubleAttr_NodeCount);
		cerr << "# B&B nodes = " << NumOfBandBNodes << endl;
		cerr << "# callbacks = " << Kclub_callback::numCallbacks << endl;
		cerr << "# lazy cuts = " << Kclub_callback::numLazyCuts << endl;

		double bestLB = model.get(GRB_DoubleAttr_ObjVal);
		double bestUB = model.get(GRB_DoubleAttr_ObjBound);
		cout << bestLB << " " << bestUB << " ";

		int status = model.get(GRB_IntAttr_Status);
		if (status == GRB_OPTIMAL)
		{
			subOpt = false;
			for (long i = 0; i < g.n; i++)
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)
					S.push_back(i);
		}

		delete[] X;
	}

	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return S;
}

vector<long> ICUT_Subproblem(KGraph &g, long k, long v_i, long LowerBound, long &SubProblemLB, long &SubProblemUB, double teta)
{
	vector<long> S;
	vector< vector< long> > components;
	vector<long> degreeZero;
	g.FindConnectedComponents(components, degreeZero);

	try
	{
		GRBEnv env = GRBEnv();
		env.set(GRB_IntParam_Method, 3); //concurrent 
		env.set(GRB_IntParam_OutputFlag, 0);
		env.set(GRB_DoubleParam_TimeLimit, 3600);
		env.set(GRB_IntParam_LazyConstraints, 1);
		env.set(GRB_DoubleParam_MIPGap, 0.0);
		env.set(GRB_DoubleParam_Cutoff, LowerBound);

		// when |s-club| > teta + 1.5, the violation is found and solver stops.
		env.set(GRB_DoubleParam_BestObjStop, teta + 1.5);
	
		GRBModel model = GRBModel(env);
		GRBVar *X = model.addVars(g.n, GRB_BINARY);
		GRBVar *Y = model.addVars(components.size(), GRB_BINARY);
		model.update();


		// Adding objective function
		for (long w = 0; w < g.n; w++)
			X[w].set(GRB_DoubleAttr_Obj, 1);
		model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
		model.update();

		// Fixing singletons to zero.
		for (long i = 0; i < degreeZero.size(); i++)
		{
			long v = degreeZero[i];
			model.addConstr(X[v] == 0);
		}

		// Fixing y[i]=0 when |V(G_i)| < lb.
		for (long i = 0; i < components.size(); i++)
		{
			if (components[i].size() < LowerBound)
			{
				model.addConstr(Y[i] == 0);
			}
		}

		// Adding x_u + x_v <= 1 constraints
		// only need to add such constraints when u and v belong to same component.
		for (long i = 0; i < components.size(); i++)
		{
			// add all constraints within component i
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j]; // v is a vertex in component i.
				vector<long> dist = g.ShortestPathsUnweighted(v);
				for (long q = j + 1; q < components[i].size(); q++)
				{
					long u = components[i][q];  // u is another vertex in component i
												// vertices u and v belong to same component i
												// if dist(v,u) > s, then add constraint x_u + x_v <= 1.
					if (dist[u] > k)
					{
						model.addConstr(X[u] + X[v] <= 1);
					}
				}
			}
		}
		model.update();

		// Adding \sigma_i y_i <= 1 constraints.
		GRBLinExpr expr = 0;
		for (long i = 0; i < components.size(); i++)
		{
			expr += Y[i];
		}
		model.addConstr(expr <= 1);
		model.update();


		// Adding X[v] <= Y[i] constraints.
		for (long i = 0; i < components.size(); i++)
		{
			for (long j = 0; j < components[i].size(); j++)
			{
				long v = components[i][j];  // v is a vertex in component i
				model.addConstr(X[v] <= Y[i]);
			}
		}
		// Fixing X[v_i] to 1.
		model.addConstr(X[v_i] == 1);
		model.update();

		//Adding lazy constraints.
		Kclub_callback cb = Kclub_callback(X, g, k);
		model.setCallback(&cb);

		model.optimize();

		SubProblemLB = model.get(GRB_DoubleAttr_ObjVal);
		SubProblemUB = model.get(GRB_DoubleAttr_ObjBound);

		int status = model.get(GRB_IntAttr_Status);

		if (status == GRB_CUTOFF)
		{
			SubProblemLB = LowerBound;
			SubProblemUB = LowerBound;
			return S;
		}	
		else if (status == GRB_USER_OBJ_LIMIT)
		{
			SubProblemLB = model.get(GRB_DoubleAttr_ObjVal);
			SubProblemUB = model.get(GRB_DoubleAttr_ObjBound);
			for (long i = 0; i < g.n; i++)
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)
					S.push_back(i);
		}
		else if (status == GRB_OPTIMAL)
		{
			SubProblemLB = model.get(GRB_DoubleAttr_ObjVal);
			SubProblemUB = model.get(GRB_DoubleAttr_ObjBound);
			for (long i = 0; i < g.n; i++)
				if (X[i].get(GRB_DoubleAttr_X) > 0.5)
					S.push_back(i);
		}

		delete[] X;
	}

	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Error during optimization" << endl;
	}
	return S;
}

vector <long> ICUT(KGraph &g, long k, vector<long> &BestKClub, double teta, long &counter)
{
	time_t start = clock();
	KGraph gk = g.CreatePowerGraph(k);

	long BestLowerBound = BestKClub.size();
	long BestUpperBound = BestKClub.size();

	vector<long> degeneracyordering = gk.FindDegeneracyOrdering();
	vector<bool> T(g.n, false);

	for (long i = g.n - 1; i >= 0; i--)
	{
		vector<long> SubproblemVertices;

		long v = degeneracyordering[i];
		T[v] = true;
		vector<long> dist_from_v = g.ShortestPathsUnweighted(v, T);

		for (long j = 0; j < g.n; j++)
			if (dist_from_v[j] <= k)
				SubproblemVertices.push_back(j);

		if (SubproblemVertices.size() <= BestLowerBound) continue;

		vector<long> Rmap;
		KGraph g_subproblem = g.CreateInducedGraph(SubproblemVertices, Rmap);

		long SubProblemLB;
		long SubProblemUB;

		vector<long> SubproblemSolution = ICUT_Subproblem(g_subproblem, k, Rmap[v], BestLowerBound, SubProblemLB, SubProblemUB, teta);

		BestLowerBound = max(BestLowerBound, SubProblemLB);
		BestUpperBound = max(BestUpperBound, SubProblemUB);

		if (SubproblemSolution.size() >= BestLowerBound)
		{
			BestKClub = SubproblemSolution;
			for (long q = 0; q < BestLowerBound; q++)
			{
				long v = BestKClub[q];
				long w = SubproblemVertices[v];
				BestKClub[q] = w;
			}
			if (SubproblemSolution.size() > teta + 1.5)    // the violation is found
			{
				counter++;        // count the number of times violation is found without solving ICUT to optimality
				break;            // don't need to call the ICUT_subproblem anymore.
			}
		}
	}

	cout << "bestLB = " << BestLowerBound << ", bestUB = " << BestUpperBound << " ";
	return BestKClub;
}
