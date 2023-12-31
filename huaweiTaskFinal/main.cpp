#include <iostream>               
//#include <utility>                   
//#include <algorithm>
#include <boost/config.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <random>

//#include <typeinfo>

double delay(std::normal_distribution<double>& d, std::mt19937& gen) {
	double temp = d(gen);
	return temp < 0 ? 0 : temp;          // No sense in negative delay edges
}

double loss(std::vector<double>& losses, std::discrete_distribution<>& d, std::mt19937& gen)
{
	int temp = d(gen);

	if (temp >= losses.size() || temp < 0) {
		return 0;
	}
	return losses[temp];
}



typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
//typedef boost::property<boost::vertex_distance_t, double> Distance;
//typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, Distance, EdgeWeightProperty> Graph;
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS, boost::no_property, EdgeWeightProperty> Graph;

typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;
typedef boost::graph_traits<Graph>::vertex_descriptor vertex_descriptor;



class Solution
{
private:
	Graph* G_delay;
	Graph* G_loss;

public:
	Solution(Graph& d, Graph& l) {
		G_delay = &d;
		G_loss = &l;
	}

	/*
	* Firstly, I will calculate the shortest paths with delay and with loss as a metric,
	* 
	* then I will construct SPT of delay metric and fill the edges weight with losses, and calculate loss distances
	* 
	* in the delay tree. 
	* 
	* To compare the cost of paths in losses to every vertex I will compare final loss value.
	* 
	* It is easier than path comparison and it prevents a situation, when another path, but with the same cost value 
	* 
	* is considered as wrong path.
	* 
	* To solve this task I use boost library.
	*/
	
	int optimal_routes(int root)
	{
		if (root >= boost::num_vertices(*G_delay))
		{
			std::cout << "root vertex out of range" << std::endl;
			return 0;
		}
		
		// Dijkstra for the graph with delay as an edge weight 

		boost::property_map< Graph, boost::edge_weight_t >::type weightmap = get(boost::edge_weight, *G_delay);
		
		std::vector< vertex_descriptor > parents_delay(num_vertices(*G_delay));  // vector of parent vertices in delay graph SPT
		std::vector< double > distance_delay(num_vertices(*G_delay));            // vector of distances to the vertices of the delay graph
		
		vertex_descriptor s = vertex(root, *G_delay);  // start vertex
		

		auto predecessor = boost::predecessor_map(boost::make_iterator_property_map(parents_delay.begin(), get(boost::vertex_index, *G_delay)))
			.distance_map(boost::make_iterator_property_map(
			distance_delay.begin(), get(boost::vertex_index, *G_delay)));

		boost::dijkstra_shortest_paths(*G_delay, s, predecessor);


		// Dijkstra for the graph with loss as an edge weight
		weightmap = get(boost::edge_weight, *G_loss);
		
		std::vector< vertex_descriptor > parents_loss(num_vertices(*G_loss));    // vector of parent vertices in loss graph SPT
		std::vector< double > distance_loss(num_vertices(*G_loss));              // vector of distances to the vertices of the loss graph
		
		s = vertex(root, *G_loss);

		predecessor = boost::predecessor_map(boost::make_iterator_property_map(parents_loss.begin(), get(boost::vertex_index, *G_loss))).distance_map(boost::make_iterator_property_map(
			distance_loss.begin(), get(boost::vertex_index, *G_loss)));

		boost::dijkstra_shortest_paths(*G_loss, s, predecessor);


		// Construction of SPT for delay graph with losses as edge weight
		Graph delay_SPT(boost::num_vertices(*G_loss));   //SPT-graph with losses as edge weight
		for (int i = 0; i < parents_delay.size(); i++)
		{
			if (i == parents_delay[i]) { //parent vertex of the root is root 
				continue;
			}
			add_edge(i, parents_delay[i], get(boost::edge_weight_t(), *G_loss, boost::edge(i, parents_delay[i], *G_loss).first), delay_SPT);
		}

		weightmap = get(boost::edge_weight, delay_SPT);

		std::vector< vertex_descriptor > parents_spt(num_vertices(delay_SPT));
		std::vector< double > distance_spt(num_vertices(delay_SPT));
		
		s = vertex(root, delay_SPT);
		predecessor = boost::predecessor_map(boost::make_iterator_property_map(parents_loss.begin(), get(boost::vertex_index, delay_SPT))).distance_map(boost::make_iterator_property_map(
			distance_loss.begin(), get(boost::vertex_index, delay_SPT)));
		
		boost::dijkstra_shortest_paths(delay_SPT, s, predecessor);

		int optimal_routes_num = 0;
		for (int i = 0; i < distance_spt.size(); i++)
		{
			if (distance_spt[i] == distance_loss[i]) {
				optimal_routes_num++;
			}
		}
		//return double(optimal_routes_num) / boost::num_vertices(*G_loss); 

		return optimal_routes_num;    // for optimisation I will divide it once in the end 
	}
};


void two_grids_constructor(Graph& G_delay, Graph& G_loss, int N, std::vector<double>& losses, std::normal_distribution<double>& norm_d, std::discrete_distribution<>& discr_d)
{
	std::random_device rd;
	std::mt19937 gen(rd());

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N - 1; j++)
		{
			add_edge(i * N + j, i * N + j + 1, delay(norm_d, gen), G_delay);
			add_edge(i * N + j, i * N + j + 1, loss(losses, discr_d, gen), G_loss);
			if (i < N - 1) {
				add_edge(i * N + j, (i + 1) * N + j, delay(norm_d, gen), G_delay);
				add_edge(i * N + j, (i + 1) * N + j, loss(losses, discr_d, gen), G_loss);
			}
		}
		if (i < 31) {
			add_edge(i * N + (N - 1), (i + 1) * N + (N - 1), delay(norm_d, gen), G_delay);
			add_edge(i * N + (N - 1), (i + 1) * N + (N - 1), loss(losses, discr_d, gen), G_loss);
		}
	}
}

int main()
{
	const int N = 32;
	int NUM_EXPERIMENTS = 100;

	/* To determine the number of experiments that we need I used the Chebyshev theorem.
	* 
	* It requires an estimate of the variance, which I roughly estimated as (10/1000)^2.
	* 
	* With 100 experiments, this should be 0.002 precision of E estimation with 0.75 probability.
	* 
	* But after launching this program many times and compare results, I think, the precision is much lower,
	* 
	* so variance should be also much lower.
	*/
	
	std::vector<double> losses = { 0, 0.01, 0.02, 0.03, 0.04, 0.05 };
	std::normal_distribution<double> norm_d {7.5, 1.25};
	std::discrete_distribution<> discr_d({ 0.5, 0.1, 0.1, 0.1, 0.1, 0.1 });


	std::vector<double> optimal_routes_average(N); // vector with number of optimal routes for first N vertices (boundary vertices)

	time_t start, end;
	time(&start);


	/*
	* Optimal routes for bound vertices of the same graph are not independent,
	* 
	* but it does not matter, because to i-vertex results I add only new i-vertex result computed
	* 
	* on the new graph.
	*/


	for(int i = 0; i < NUM_EXPERIMENTS; i++)
	{

		Graph G_delay(N * N);
		Graph G_loss(N * N);
		
		two_grids_constructor(G_delay, G_loss, N, losses, norm_d, discr_d);
		
		Solution S(G_delay, G_loss);
		

		for (int j = 0; j < N; j++)
		{	
			optimal_routes_average[j] += S.optimal_routes(j);
		}
	}

	time(&end);
	
	for (int i = 0; i < N; i++)
	{
		std::cout << "average for vertex " << i << " : " << optimal_routes_average[i] / (NUM_EXPERIMENTS * N * N) << std::endl;
	}
	
	/*double time_taken = double(end - start);
	std::cout << "Time  : "
		<< time_taken << " ";
	std::cout << " sec " << std::endl;
	*/
	
	return 0;
}