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
#include<iostream>
#include<algorithm>

using namespace std;
using namespace std::chrono;


bool sortcol(const vector<long>& v1, const vector<long>& v2)
{
	return v1[1] > v2[1];
}


int S_in_Sclb;
long num_BB_interdiction;
long num_callbacks_interdiction;
long num_Lazycuts_interdiction_1;
long num_Lazycuts_interdiction_2;
long num_Lazycuts_interdiction_3;
double KclubTime = 0;
double CallbackTime = 0;
double total_sclub_time = 0;
int method;



class addcut_theta : public GRBCallback
{
public:

	GRBVar* var;
	GRBVar* Tta;
	long nvar;
	KGraph graph_1;

	addcut_theta(GRBVar* xvar, GRBVar* Teta, long Nvar, KGraph &graph_2)
	{
		var = xvar;
		Tta = Teta;
		nvar = Nvar;
		graph_1 = graph_2;
	}

protected:

	void callback() {
		try
		{
			if (where == GRB_CB_MIPSOL)
			{

				num_callbacks_interdiction++;
				double *x_master = new double[nvar];
				x_master = getSolution(var, nvar);

				vector <long> non_interdicted_vertices;

				for (int p2 = 0; p2 < nvar; p2++)
				{
					if (x_master[p2] <= 0.5)
					{
						non_interdicted_vertices.push_back(p2);
					}
				}

				double *THETA = new double[1];
				THETA = getSolution(Tta, 1);

				// mapping to find new adj list
				vector<long> ReverseMap;
				KGraph induced_g = graph_1.CreateInducedGraph(non_interdicted_vertices, ReverseMap);

				vector <long> kclb_index;
				GRBLinExpr cut_in_master = 0;
				GRBLinExpr Star = 0;
				GRBLinExpr Leaves = 0;
				GRBLinExpr Critical_vertices = 0;
				bool leaf = false;
				bool lazycut_added = false;

				if (induced_g.n >= 2)
				{
					kclb_index.clear();
					cut_in_master = 0;
					Star = 0;
					Leaves = 0;
					Critical_vertices = 0;
					leaf = false;
					lazycut_added = false;

					auto start_sclub = chrono::steady_clock::now();
					vector <long> HS;
					HS = HeuristicAndPreprocess(induced_g, S_in_Sclb);
					bool Sub_Opt;
					if (method == 0)
					{
						//use cutlike
						kclb_index = solveMaxKClub_CutLike(induced_g, S_in_Sclb, HS, Sub_Opt);
					}
					else if (method == 1)
					{
						//use icut
						kclb_index = ICUT(induced_g, S_in_Sclb, HS);
					}
					chrono::duration <double> duration_sclb = chrono::steady_clock::now() - start_sclub;
					total_sclub_time += duration_sclb.count();

					vector<long> kclb_original_index;
					for (int i = 0; i < kclb_index.size(); i++)
					{
						kclb_original_index.push_back(non_interdicted_vertices[kclb_index[i]]);
					}
					KGraph induced_kclb = graph_1.CreateInducedGraph(kclb_original_index, ReverseMap);


					if (THETA[0] < kclb_index.size())
					{
						for (int i = 0; i < kclb_index.size(); i++)
						{
							if (lazycut_added == false)
							{
								cut_in_master += var[non_interdicted_vertices[kclb_index[i]]];

								if (induced_kclb.degree[i] == 1)
								{
									Leaves += var[non_interdicted_vertices[kclb_index[i]]];     // linear expr for cuts for leaves
									leaf = true;
								}

								if (induced_kclb.degree[i] > 1)
								{
									Critical_vertices += var[non_interdicted_vertices[kclb_index[i]]];
								}
								if (induced_kclb.degree[i] == kclb_index.size() - 1)
								{
									// add cut for a star
									Leaves = 0;
									Critical_vertices = 0;
									leaf = false;
									cut_in_master = 0;

									for (int j = 0; j < kclb_index.size(); j++)
									{
										Star += var[non_interdicted_vertices[kclb_index[j]]];
									}
									Star -= var[non_interdicted_vertices[kclb_index[i]]];
									addLazy(Tta[0] >= (1 - var[non_interdicted_vertices[kclb_index[i]]]) * (kclb_index.size()) - Star);
									num_Lazycuts_interdiction_1++;
									Star = 0;
									lazycut_added = true;
								}
							}
						}
						if (leaf == true && lazycut_added == false)
						{
							cut_in_master = 0;
							addLazy(Tta[0] >= (1 - Critical_vertices) * (kclb_index.size()) - Leaves);
							num_Lazycuts_interdiction_2++;
							Leaves = 0;
							Critical_vertices = 0;
							leaf = false;
							lazycut_added = true;
						}
						if (lazycut_added == false)
						{
							Leaves = 0;
							Critical_vertices = 0;
							leaf = false;
							addLazy(Tta[0] >= kclb_index.size() * (1 - cut_in_master));
							num_Lazycuts_interdiction_3++;
							cut_in_master = 0;
							lazycut_added = true;
						}
					}

				} //end if of (g.size >= 2)
				Critical_vertices = 0;
				Star = 0;
				Leaves = 0;
				leaf = false;
				lazycut_added = false;
				cut_in_master = 0;
				non_interdicted_vertices.clear();
				kclb_index.clear();
				delete[] x_master;
				delete[] THETA;
			} // end if loop

		}  //end try

		catch (GRBException e)
		{
			cout << "Error number: " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...)
		{
			cout << "Error during callback" << endl;
		}

	}      //end callback

};  // end addcut_theta class



int main(int argc, char *argv[])
{

	auto start = chrono::steady_clock::now();
	if (argc < 2)
		cerr << "ERROR: Not enough arguments.";

	else if (strcmp(argv[1], "CutLike") == 0)
	{
		method = 0;
	}
	else if (strcmp(argv[1], "ICUT") == 0)
	{
		method = 1;
	}

	KGraph grph(argv[3], argv[3], argv[2]);
	double alpha = atof(argv[4]);
	S_in_Sclb = atoi(argv[5]);
	long num_interdicted_vertices = 0;


	// sort vertices based on the degree 
	int highest_degree = floor(grph.n * 0.2);

	vector< vector <long>> sorted_degree(grph.n);

	for (int i = 0; i < grph.n; i++)
	{
		sorted_degree[i].push_back(i);
		sorted_degree[i].push_back(grph.degree[i]);
	}

	sort(sorted_degree.begin(), sorted_degree.end(), sortcol);

	// Master (interdiction) Problem 
	try {
		GRBEnv env_master = GRBEnv();
		GRBModel model_Master = GRBModel(env_master);

		//Specify the use of lazy constraints
		model_Master.getEnv().set(GRB_IntParam_LazyConstraints, 1);

		// veriables
		GRBVar* X_Master = model_Master.addVars(grph.n, GRB_BINARY);
		GRBVar* theta = model_Master.addVars(1, GRB_CONTINUOUS);

		model_Master.addConstr(theta[0] >= 0);
		model_Master.update();

		//set obj coefficient
		for (int p1 = 0; p1 < grph.n; ++p1)
		{
			X_Master[p1].set(GRB_DoubleAttr_Obj, alpha);
		}
		theta[0].set(GRB_DoubleAttr_Obj, 1);

		// add constr for 20% of the vertices
		GRBLinExpr neighbors_of_v = GRBLinExpr();
		long v2;

		for (int v1 = 0; v1 < highest_degree; v1++)
		{
			v2 = sorted_degree[v1][0];
			for (int v3 = 0; v3 < grph.degree[v2]; v3++)
			{
				neighbors_of_v += X_Master[grph.adj[v2][v3]];
			}

			model_Master.addConstr(theta[0] >= ((1 - X_Master[v2]) * (grph.degree[v2] + 1)) - neighbors_of_v);
			neighbors_of_v = 0;
		}
		neighbors_of_v = 0;


		//SET GUROBI PARAMETERS
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

		model_Master.update();

		// set callback
		addcut_theta cb1 = addcut_theta(X_Master, theta, grph.n, grph);

		model_Master.setCallback(&cb1);
		model_Master.optimize();

		long num_BB_Nodes = (long)model_Master.get(GRB_DoubleAttr_NodeCount);

		//update curr_best_kclb when better solution is available
		if (model_Master.get(GRB_IntAttr_SolCount) == 0)
		{
			cout << "No solution found, Gurobi optimization status = " << model_Master.get(GRB_IntAttr_Status) << endl;
		}
		else
		{

			double obj_master = model_Master.get(GRB_DoubleAttr_ObjVal);
			cout << "obj = " << obj_master << endl;
			for (int i = 0; i < grph.n; i++)
			{
				if (X_Master[i].get(GRB_DoubleAttr_X) > 0.5)
				{
					num_interdicted_vertices++;
				}
			}
			cout << "num_interdicted_vertices : " << num_interdicted_vertices << endl;
			cout << "theta : " << theta[0].get(GRB_DoubleAttr_X) << endl;
			cout << endl;
			cout << "*****************************" << endl;
			cout << "# B&B nodes in interdiction = " << num_BB_Nodes << endl;
			cout << "# of callbacks in interdiction  = " << num_callbacks_interdiction << endl;
			cout << "# of lazy cuts in interdiction(Star) = " << num_Lazycuts_interdiction_1 << endl;
			cout << "# of lazy cuts in interdiction (leaves) = " << num_Lazycuts_interdiction_2 << endl;
			cout << "# of lazy cuts in interdiction (regular) = " << num_Lazycuts_interdiction_3 << endl;

		}

	}

	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;

	}
	catch (...) {
		cout << "Exception during optimization" << endl;

	}

	// Total Time 
	chrono::duration <double> duration = chrono::steady_clock::now() - start;
	cout << "Total Time : " << duration.count() << "s " << endl;

	//print time for max kclub problem
	printf("kclb Time : %.2fs\n", total_sclub_time);

	return 0;
}