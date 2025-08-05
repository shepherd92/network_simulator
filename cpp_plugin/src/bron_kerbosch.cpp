#include "bron_kerbosch.h"

void bron_kerbosch(
    std::unordered_set<PointId> R,
    std::unordered_set<PointId> P,
    std::unordered_set<PointId> X,
    const std::unordered_map<PointId, std::unordered_set<PointId>> &graph,
    ISimplexList &cliques);

ISimplexList find_maximal_cliques(const PointIdList &vertices, const ConnectionList &edges)
{
    std::unordered_map<PointId, std::unordered_set<PointId>> graph;

    // Build undirected adjacency list
    for (const auto &vertex : vertices)
    {
        graph[vertex] = {};
    }
    for (const auto &edge : edges)
    {
        std::cout << "Adding edge: " << edge.first << " - " << edge.second << std::endl;
        graph[edge.first].insert(edge.second);
        graph[edge.second].insert(edge.first);
    }

    ISimplexList cliques;
    std::unordered_set<PointId> P(vertices.begin(), vertices.end());
    std::unordered_set<PointId> R, X;

    bron_kerbosch(R, P, X, graph, cliques);
    return cliques;
}

void bron_kerbosch(
    std::unordered_set<PointId> R, // current clique being built
    std::unordered_set<PointId> P, // potential candidates to extend the clique
    std::unordered_set<PointId> X, // vertices already considered
    const std::unordered_map<PointId, std::unordered_set<PointId>> &graph,
    ISimplexList &cliques)
{
    if (P.empty() && X.empty())
    {
        cliques.emplace_back(R.begin(), R.end());
        return;
    }

    std::unordered_set<PointId> Pcopy = P; // Copy so we can modify P in loop
    for (const int v : Pcopy)
    {
        std::unordered_set<PointId> newR = R;
        newR.insert(v);

        std::unordered_set<PointId> newP;
        std::unordered_set<PointId> newX;

        for (int u : graph.at(v))
        {
            if (P.count(u))
            {
                newP.insert(u);
            }
            if (X.count(u))
            {
                newX.insert(u);
            }
        }

        bron_kerbosch(newR, newP, newX, graph, cliques);

        P.erase(v);
        X.insert(v);
    }
}