#ifndef _FINITE_NETWORK_H_
#define _FINITE_NETWORK_H_

#include "network.h"
#include "typedefs.h"

class FiniteNetwork : public Network
{
public:
    FiniteNetwork(const Dimension max_dimension);
    SimplexHandleList get_simplices() override;
    uint32_t num_simplices() override;
};

#endif