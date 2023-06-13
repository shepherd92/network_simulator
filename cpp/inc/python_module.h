#ifndef _PYTHON_MODULE_H_
#define _PYTHON_MODULE_H_

#include "adrcm_model.h"
#include "network.h"
#include "typedefs.h"

class AdrcmModelInterface
{
public:
private:
    AdrcmModel model;
};

class NetworkInterface
{
public:
    NetworkInterface(const dimension max_dimension);
    py::array_t<vertex_id> get_simplices_by_dimension(const dimension dimension) const;
    void add_simplices(const py::list &simplices);
    py::list facets();

private:
    Network network;
};

#endif