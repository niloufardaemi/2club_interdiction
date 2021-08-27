#include "gurobi_c++.h"
#include "GRBInterface.h"
#include "KGraph.h"
#include <sstream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <chrono>
#include <iostream>
#include <algorithm>

using namespace std;
using namespace std::chrono;


bool sortcol(const vector<long>& v1, const vector<long>& v2)
{
	return v1[1] > v2[1];
}

long s = 2;    //the parameter s in s-club

// global variables to count the number of callbacks and lazy cuts
long num_callbacks_interdiction;
long num_Lazycuts_interdiction_1;   // lazy cuts for s-clubs that form a star
long num_Lazycuts_interdiction_2;	// lazy cuts for s-clubs with leavs
long num_Lazycuts_interdiction_3;	// regular lazy cuts

// the total time to solve the maximum s-club problem in all the iterations
double SclubTime = 0;   

// count the number of iterations the solution is found by heuristics
long HS_Counter = 0;
long ICUT_Counter = 0;


class addlazycut_theta : public GRBCallback
{
public:

	GRBVar* Xvar;
	GRBVar Theta;
	long n_Xvar;
	KGraph graph_1;

	addlazycut_theta(GRBVar* xvar, GRBVar Teta, long N_xvar, KGraph &graph_2)
	{
		Xvar = xvar;
		Theta = Teta;
		n_Xvar = N_xvar;
		graph_1 = graph_2;
	}

protected:

	void callback() {
		try
		{
			if (where == GRB_CB_MIPSOL)
			{
				num_callbacks_interdiction++;  // count the number of callbacks

				// get the solution of the master problem
				double *x_master = new double[n_Xvar];
				x_master = getSolution(Xvar, n_Xvar);
				double THETA = getSolution(Theta);

				// if value of variable x is less than or equal to 0.5, it is not interdicted so x = 0
				vector <long> non_interdicted_vertices;
				for (long p2 = 0; p2 < n_Xvar; p2++)
				{
					if (x_master[p2] <= 0.5)
					{
						non_interdicted_vertices.push_back(p2);
					}
				}

				// mapping to find the adjacency list of the graph after interdiction
				vector<long> ReverseMap;
				KGraph induced_g = graph_1.CreateInducedGraph(non_interdicted_vertices, ReverseMap);
				
				// define requrired structures to solve the separation problem (max s-club) and add the lazy cuts
				vector <long> sclb_index;
				vector <long> HS;     // to keep the heuuristic s-club
				GRBLinExpr cut_in_master = 0;
				GRBLinExpr Star = 0;
				GRBLinExpr Leaves = 0;
				GRBLinExpr Critical_vertices = 0;
				bool leaf = false;
				bool lazycut_added = false;

				// solve the separation only if the interdicted graph has at least 2 vertices
				if (induced_g.n >= 2)
				{								
					auto start_HS_sclub = chrono::steady_clock::now();   // begin to compute the time for solving the separtion problem
					HS = HeuristicAndPreprocess(induced_g, s);         // find the heuristic s-club in the interdicted graph
					chrono::duration <double> duration_sclb = chrono::steady_clock::now() - start_HS_sclub;  // duration of finding the heuristic s-club
					SclubTime += duration_sclb.count();

					//  if heuristic approach finds the violation, add a lazy cut using the heuristic s-club:
					if (THETA + 1.5 < HS.size())
					{
						// count number of times the violation is found by heuritic method without calling ICUT 
						HS_Counter++;       
						
						// converting the index of vertices to the indices in the original graph (before interdiction)
						vector<long> HS_original_index;
						for (int i = 0; i < HS.size(); i++)
						{
							HS_original_index.push_back(non_interdicted_vertices[HS[i]]);
						}
						KGraph induced_HS = graph_1.CreateInducedGraph(HS_original_index, ReverseMap);


						// add a lazy cut:
						for (long i = 0; i < HS.size(); i++)
						{
							if (lazycut_added == false)
							{
								cut_in_master += Xvar[non_interdicted_vertices[HS[i]]];

								if (induced_HS.degree[i] == 1)
								{
									Leaves += Xvar[non_interdicted_vertices[HS[i]]];  
									leaf = true;
								}
								if (induced_HS.degree[i] > 1)
								{
									Critical_vertices += Xvar[non_interdicted_vertices[HS[i]]];
								}
								if (induced_HS.degree[i] == HS.size() - 1)
								{
									// in this case, the s-club is a star with center i 
									// add cut where H is all the vertices except i
									Leaves = 0;
									Critical_vertices = 0;
									leaf = false;
									cut_in_master = 0;
									for (int j = 0; j < HS.size(); j++)
									{
										Star += Xvar[non_interdicted_vertices[HS[j]]];
									}
									Star -= Xvar[non_interdicted_vertices[HS[i]]];								
									addLazy(Theta >= (1 - Xvar[non_interdicted_vertices[HS[i]]]) * (HS.size()) - Star);
									num_Lazycuts_interdiction_1++;
									Star = 0;
									lazycut_added = true;
								}
							}
						}
						if (leaf == true && lazycut_added == false)
						{
							// in this case, the sclub contains leavs
							// add cut where H is the set of the leavs in the s-club
							cut_in_master = 0;
							addLazy(Theta >= (1 - Critical_vertices) * (HS.size()) - Leaves);
							num_Lazycuts_interdiction_2++;
							Leaves = 0;
							Critical_vertices = 0;
							leaf = false;
							lazycut_added = true;
						}
						if (lazycut_added == false)
						{
							// no H is found, add a regular cut where H is empty
							Leaves = 0;
							Critical_vertices = 0;
							leaf = false;
							addLazy(Theta >= HS.size() * (1 - cut_in_master));
							num_Lazycuts_interdiction_3++;
							cut_in_master = 0;
							lazycut_added = true;
						}
					}
					
					
					//if heuristic approach does not the violation, call ICUT to find the maximum s-club:
					else  
					{
						auto start_sclub = chrono::steady_clock::now();   
						// call ICUT to solve the maximum the s-club problem 
						// ICUT is not solved to optimality; when a solution greater that THETA + 1.5 is found, solver stops
						sclb_index = ICUT(induced_g, s, HS, THETA, ICUT_Counter);
						chrono::duration <double> duration_sclb = chrono::steady_clock::now() - start_sclub;
						SclubTime += duration_sclb.count();

						// converting the index of vertices to the indices in the original graph (before interdiction)
						vector<long> sclb_original_index;
						for (int i = 0; i < sclb_index.size(); i++)
						{
							sclb_original_index.push_back(non_interdicted_vertices[sclb_index[i]]);
						}
						KGraph induced_kclb = graph_1.CreateInducedGraph(sclb_original_index, ReverseMap);

						// add lazy cut if ICUT has found a violation
						if (THETA < sclb_index.size())
						{
							for (long i = 0; i < sclb_index.size(); i++)
							{
								if (lazycut_added == false)
								{
									cut_in_master += Xvar[non_interdicted_vertices[sclb_index[i]]];

									if (induced_kclb.degree[i] == 1)
									{
										Leaves += Xvar[non_interdicted_vertices[sclb_index[i]]];     // linear expr for cuts for leaves
										leaf = true;
									}

									if (induced_kclb.degree[i] > 1)
									{
										Critical_vertices += Xvar[non_interdicted_vertices[sclb_index[i]]];
									}
									if (induced_kclb.degree[i] == sclb_index.size() - 1)
									{
										// in this case, the s-club is a star with center i 
										// add cut where H is all the vertices except i
										Leaves = 0;
										Critical_vertices = 0;
										leaf = false;
										cut_in_master = 0;

										for (int j = 0; j < sclb_index.size(); j++)
										{
											Star += Xvar[non_interdicted_vertices[sclb_index[j]]];
										}
										Star -= Xvar[non_interdicted_vertices[sclb_index[i]]];										
										addLazy(Theta >= (1 - Xvar[non_interdicted_vertices[sclb_index[i]]]) * (sclb_index.size()) - Star);
										num_Lazycuts_interdiction_1++;
										Star = 0;
										lazycut_added = true;
									}
								}
							}
							if (leaf == true && lazycut_added == false)
							{
								// in this case, the sclub contains leavs
								// add cut where H is the set of the leavs in the s-club
								cut_in_master = 0;					
								addLazy(Theta >= (1 - Critical_vertices) * (sclb_index.size()) - Leaves);
								num_Lazycuts_interdiction_2++;
								Leaves = 0;
								Critical_vertices = 0;
								leaf = false;
								lazycut_added = true;
							}
							if (lazycut_added == false)
							{
								// no H is found, add a regular cut where H is empty
								Leaves = 0;
								Critical_vertices = 0;
								leaf = false;
								addLazy(Theta >= sclb_index.size() * (1 - cut_in_master));
								num_Lazycuts_interdiction_3++;
								cut_in_master = 0;
								lazycut_added = true;
							}
						}
					} 
				} 
				// reset the expressions and variables
				Critical_vertices = 0;
				Star = 0;
				Leaves = 0;
				leaf = false;
				lazycut_added = false;
				cut_in_master = 0;
				non_interdicted_vertices.clear();
				sclb_index.clear();
				HS.clear();
				delete[] x_master;
			} 
		}

		catch (GRBException e)
		{
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...)
		{
			cout << "Error during callback" << endl;
		}

	}     

};



int main(int argc, char *argv[])
{

	auto start = chrono::steady_clock::now();
	if (argc < 2)
		cerr << "ERROR: Not enough arguments.";

	// read the input graph, format of the dataset and the value of the penalty (alpha)
	KGraph grph(argv[2], argv[2], argv[1]);
	double alpha = atof(argv[3]);

	// sort vertices based on the degree 
	vector< vector <long>> sorted_degree(grph.n);
	for (long i = 0; i < grph.n; i++)
	{
		sorted_degree[i].push_back(i);
		sorted_degree[i].push_back(grph.degree[i]);
	}
	sort(sorted_degree.begin(), sorted_degree.end(), sortcol);


	// Master (interdiction) Problem 
	try {
		GRBEnv env_master = GRBEnv();
		GRBModel model_Master = GRBModel(env_master);

		// veriables
		GRBVar* X_Master = model_Master.addVars(grph.n, GRB_BINARY);
		GRBVar theta = model_Master.addVar(0.0, INFINITY, 1.0, GRB_CONTINUOUS);
		model_Master.update();

		//set obj coefficient
		for (long p1 = 0; p1 < grph.n; ++p1)
		{
			X_Master[p1].set(GRB_DoubleAttr_Obj, alpha);
		}
		theta.set(GRB_DoubleAttr_Obj, 1);

		// add constr for the top 20% of vertices by degree	
		long v2;
		long highest_degree = floor(grph.n * 0.2);   // the number of vertices we add a constraint for them
		GRBLinExpr neighbors_of_v = GRBLinExpr();
		for (long v1 = 0; v1 < highest_degree; v1++)
		{
			v2 = sorted_degree[v1][0];
			for (long v3 = 0; v3 < grph.degree[v2]; v3++)
			{
				neighbors_of_v += X_Master[grph.adj[v2][v3]];
			}
			model_Master.addConstr(theta >= ((1 - X_Master[v2]) * (grph.degree[v2] + 1)) - neighbors_of_v);
			neighbors_of_v = 0;
		}
		neighbors_of_v = 0;


		//SET GUROBI PARAMETERS

		//Specify the use of lazy constraints
		model_Master.getEnv().set(GRB_IntParam_LazyConstraints, 1);

		//Set feasibility vs optimality balance
		model_Master.getEnv().set(GRB_IntParam_MIPFocus, 0);
		//1-feasible sols quickly;2-prove optimality;3-focus on MIP bound; default is 0

		//Set threads; review guidance on Gurobi.com; 0-default;
		model_Master.getEnv().set(GRB_IntParam_Threads, 0);

		//Set root node LPR solver
		model_Master.getEnv().set(GRB_IntParam_Method, -1);
		//-1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent

		//Set BC node LPR solver
		model_Master.getEnv().set(GRB_IntParam_NodeMethod, 1);
		//0=primal simplex, 1=dual simplex, 2=barrier

		//Set global cut aggressiveness; over-ridden by individual cut settings
		model_Master.getEnv().set(GRB_IntParam_Cuts, -1);
		//0=no cuts;1=moderate;2=aggressive;3=very aggressive;-1=default

		//Set maximum time limit
		model_Master.getEnv().set(GRB_DoubleParam_TimeLimit, 3600);

		//Set termination gap limit; as needed; default is 1e-4
		model_Master.getEnv().set(GRB_DoubleParam_MIPGap, 1e-4);

		//Set Gurobi screen display flag
		model_Master.getEnv().set(GRB_IntParam_OutputFlag, 1);
		//0=switch off; 1=default

		//Set objective to minimize
		model_Master.set(GRB_IntAttr_ModelSense, 1);

		// Save log file
		model_Master.set(GRB_StringParam_LogFile, "logfile");

		model_Master.update();

		// set callback
		addlazycut_theta cb1 = addlazycut_theta(X_Master, theta, grph.n, grph);
		model_Master.setCallback(&cb1);
		model_Master.optimize();

		long num_BB_Nodes = (long)model_Master.get(GRB_DoubleAttr_NodeCount);  // number of explored nodes

		//update curr_best_kclb when better solution is available
		if (model_Master.get(GRB_IntAttr_SolCount) == 0)
		{
			cout << "No solution found, Gurobi optimization status = " << model_Master.get(GRB_IntAttr_Status) << endl;
		}
		else
		{
			int status = model_Master.get(GRB_IntAttr_Status);
			double obj_master = model_Master.get(GRB_DoubleAttr_ObjVal);  // objective value of the master problem
			cout << "obj = " << obj_master << endl;
			long num_interdicted_vertices = 0;
			for (int i = 0; i < grph.n; i++)
			{
				if (X_Master[i].get(GRB_DoubleAttr_X) > 0.5)
				{
					num_interdicted_vertices++;
				}
			}
			cout << "num_interdicted_vertices : " << num_interdicted_vertices << endl;  // number of interdicted vertices
			cout << "theta : " << theta.get(GRB_DoubleAttr_X) << endl;   // size of the maximum s-club in the interdicted graph
			cout << "# B&B nodes in interdiction = " << num_BB_Nodes << endl;
			cout << "# of callbacks in interdiction  = " << num_callbacks_interdiction << endl;
			cout << "# of lazy cuts in interdiction(Star) = " << num_Lazycuts_interdiction_1 << endl;
			cout << "# of lazy cuts in interdiction (leaves) = " << num_Lazycuts_interdiction_2 << endl;
			cout << "# of lazy cuts in interdiction (regular) = " << num_Lazycuts_interdiction_3 << endl;
			cout << "number of times HS finds cut: " << HS_Counter << endl;
			cout<< "number of times ICUT sub finds cut: " << ICUT_Counter << endl;
		}
	}

	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;

	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}

	// print total time 
	chrono::duration <double> duration = chrono::steady_clock::now() - start;
	printf("Total Time : %.2fs\n", duration.count());

	//print time to solve max sclub problem
	printf("kclb Time : %.2fs\n", SclubTime);

	return 0;
}
