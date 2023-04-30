#include <algorithm>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

#include "facet_finder.h"

using namespace std;

vector<vector<int>> sort_simplices(const vector<vector<int>> &simplices)
{
  vector<vector<int>> unique_simplices = simplices;
  unique_simplices.erase(unique(unique_simplices.begin(), unique_simplices.end()), unique_simplices.end());

  // sort list of simplices
  sort(unique_simplices.begin(), unique_simplices.end(), [](const vector<int> &a, const vector<int> &b)
       { return a.size() < b.size(); });

  // sort simplices themselves
  for (auto &simplex : unique_simplices)
  {
    sort(simplex.begin(), simplex.end());
  }

  return unique_simplices;
}

vector<vector<int>> extract_facets(const vector<vector<int>> &simplices)
{
  vector<vector<int>> sorted_simplices = sort_simplices(simplices);
  vector<vector<int>> facets;
  int counter = 0;

  for (auto first = sorted_simplices.begin(); first != sorted_simplices.end(); ++first)
  {
    if (++counter % 1000 == 0)
    {
      cout << '\r'
           << setprecision(3)
           << static_cast<float>(counter) / static_cast<float>(sorted_simplices.size()) * 100
           << "%    ";
    }

    bool first_is_subset = false;
    for (auto second = std::next(first); second != sorted_simplices.end(); ++second)
    {
      first_is_subset = false;
      if (includes(second->begin(), second->end(), first->begin(), first->end()))
      {
        first_is_subset = true;
        break;
      }
    }
    if (!first_is_subset)
    {
      facets.push_back(*first);
    }
  }

  cout << endl;
  return facets;
}
