#include <vector>

#include <valgrind/callgrind.h>

#include "finite_hypergraph_model.h"
#include "finite_network.h"
#include "infinite_hypergraph_model.h"
#include "infinite_network.h"
#include "neighborhood.h"

void test_finite_hypergraph();
void test_infinite_hypergraph();

int main()
{
    test_finite_hypergraph();
    return 0;
}

void test_finite_hypergraph()
{
    const auto seed{0U};

    const std::vector<double> model_params{
        2U,    // max_dimension
        1000., // network_size
        1000., // interaction_intensity
        0.001, // beta
        0.7,   // gamma
        0.2,   // gamma_prime
        true   // weighted
    };

    std::cout << "\rCreating finite hypergraph model" << std::endl;
    const FiniteHypergraphModel model{model_params, seed};
    std::cout << "\rGenerating finite network" << std::endl;
    auto network_interface{model.generate_network()};
    auto network{std::get<0>(network_interface)};

    // std::cout << "\rCalculating degree sequence (0, 1)" << std::endl;
    // network.calc_degree_sequence(0, 1);
    // std::cout << "\rCalculating degree sequence (1, 2)" << std::endl;
    // network.calc_degree_sequence(1, 2);
    // std::cout << "\rCalculating facet dimension distribution" << std::endl;
    // network.calc_facet_dimension_distribution();
    // std::cout << "\rCalculating interaction dimension distribution" << std::endl;
    // network.calc_interaction_dimension_distribution();
    // std::cout << "\rCalculating simplex dimension distribution" << std::endl;
    // network.calc_simplex_dimension_distribution();
    // std::cout << "\rCalculating vertex interaction degree distribution" << std::endl;
    // network.calc_vertex_interaction_degree_distribution();
    std::cout << "\rCalculating Betti numbers" << std::endl;
    network.calc_betti_numbers();
    std::cout << "\rCalculating number of triangles" << std::endl;
    // network.num_simplices(2);
    // std::cout << "\rCalculating persistence pairs" << std::endl;
    // network.calc_persistence_pairs();
    // std::cout << "\rDone" << std::endl;
    std::cout << "\rCalculating persistence intervals" << std::endl;
    network.calc_persistence_intervals();
    std::cout << "\rDone" << std::endl;
}

void test_infinite_hypergraph()
{
    const auto seed{0U};

    const std::vector<double> model_params{
        2U,   // max_dimension
        1.,   // network_size
        1.,   // interaction_intensity
        1,    // beta
        0.7,  // gamma
        0.2}; // gamma_prime

    std::cout << "\rCreating infinite hypergraph model" << std::endl;
    const InfiniteHypergraphModel model{model_params, seed};
    std::cout << "\rGenerating infinite network" << std::endl;
    auto network_interface{model.generate_networks(1000)};
    auto network{std::get<0>(network_interface[0])};

    // std::cout << "\rCalculating degree sequence (0, 1)" << std::endl;
    // network.calc_degree_sequence(0, 1);
    std::cout << "\rCalculating degree sequence (1, 2)" << std::endl;
    network.calc_degree_sequence(1, 2);
    std::cout << "\rCalculating facet dimension distribution" << std::endl;
    network.calc_facet_dimension_distribution();
    std::cout << "\rCalculating interaction dimension distribution" << std::endl;
    network.calc_interaction_dimension_distribution();
    std::cout << "\rCalculating simplex dimension distribution" << std::endl;
    network.calc_simplex_dimension_distribution();
    std::cout << "\rCalculating vertex interaction degree distribution" << std::endl;
    network.calc_vertex_interaction_degree_distribution();
    std::cout << "\rCalculating number of triangles" << std::endl;
    network.num_simplices(2);
    std::cout << "\rDone" << std::endl;
}
