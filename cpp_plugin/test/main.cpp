#include <vector>

#include <valgrind/callgrind.h>

#include "finite_hypergraph_model.h"
#include "finite_network.h"
#include "infinite_hypergraph_model.h"
#include "infinite_network.h"
#include "neighborhood.h"

void test_finite_hypergraph();
void test_infinite_hypergraph();
void test_neighborhoods();
void save_points(const PointList &points, const std::string &filename);

int main()
{
    test_neighborhoods();
    return 0;
}

void test_finite_hypergraph()
{
    const auto seed{0U};

    const std::vector<double> model_params{
        2U,   // max_dimension
        10.,  // network_size
        10.,  // interaction_intensity
        0.1,  // beta
        0.7,  // gamma
        0.2}; // gamma_prime

    std::cout << "\rCreating finite hypergraph model" << std::endl;
    const FiniteHypergraphModel model{model_params, seed};
    std::cout << "\rGenerating finite network" << std::endl;
    auto network_interface{model.generate_network()};
    auto network{std::get<0>(network_interface)};

    // std::cout << "\rCalculating degree sequence (0, 1)" << std::endl;
    // network.calc_degree_sequence(0, 1);
    std::cout << "\rCalculating degree sequence (1, 2)" << std::endl;
    network.calc_degree_sequence(1, 2);
    // std::cout << "\rCalculating facet dimension distribution" << std::endl;
    // network.calc_facet_dimension_distribution();
    // std::cout << "\rCalculating interaction dimension distribution" << std::endl;
    // network.calc_interaction_dimension_distribution();
    // std::cout << "\rCalculating simplex dimension distribution" << std::endl;
    // network.calc_simplex_dimension_distribution();
    // std::cout << "\rCalculating vertex interaction degree distribution" << std::endl;
    // network.calc_vertex_interaction_degree_distribution();
    // std::cout << "\rCalculating Betti numbers" << std::endl;
    // network.calc_betti_numbers();
    // std::cout << "\rCalculating number of triangles" << std::endl;
    // network.num_simplices(2);
    // std::cout << "\rCalculating persistence pairs" << std::endl;
    // network.calc_persistence_pairs();
    // std::cout << "\rDone" << std::endl;
}

void test_infinite_hypergraph()
{
    const auto seed{0U};

    const std::vector<double> model_params{
        2U,   // max_dimension
        10.,  // network_size
        10.,  // interaction_intensity
        0.1,  // beta
        0.7,  // gamma
        0.2}; // gamma_prime

    std::cout << "\rCreating infinite hypergraph model" << std::endl;
    const InfiniteHypergraphModel model{model_params, seed};
    std::cout << "\rGenerating infinite network" << std::endl;
    auto network_interface{model.generate_networks(100)};
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

void test_neighborhoods()
{
    const auto seed{0U};

    const std::vector<double> model_params{
        2U,    // max_dimension
        1000., // network_size
        10.,   // interaction_intensity
        0.1,   // beta
        0.7,   // gamma
        0.2};  // gamma_prime

    HypergraphModel::Parameters parameters{model_params};

    std::cout << "\rCreating infinite hypergraph model" << std::endl;
    const InfiniteHypergraphModel model{model_params, seed};

    const Mark mark{0.5};
    const Position position{2.};
    const Mark transformed_mark{parameters.beta * std::pow(mark, -parameters.gamma_prime)};
    const Point interaction{transformed_mark, position};

    const auto left_tail{model.get_neighborhood_left_tail(interaction)};
    const auto center{model.get_neighborhood_center(interaction)};
    const auto right_tail{model.get_neighborhood_right_tail(interaction)};

    std::cout << "\rGenerating points in neighborhood" << std::endl;
    std::mt19937 random_number_generator{seed};
    const PointList vertices_left{left_tail.create_points(parameters, random_number_generator)};
    const PointList vertices_center{center.create_points(parameters, random_number_generator)};
    const PointList vertices_right{right_tail.create_points(parameters, random_number_generator)};

    save_points(vertices_left, "vertices_left.txt");
    save_points(vertices_center, "vertices_center.txt");
    save_points(vertices_right, "vertices_right.txt");

    std::cout << "\rDone" << std::endl;
}

void save_points(const PointList &points, const std::string &filename)
{
    std::ofstream file{filename};
    for (const auto &point : points)
    {
        file << point.position() << " " << point.mark() << std::endl;
    }
    file.close();
}