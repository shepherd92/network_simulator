#ifndef _BRON_KERBOSCH_H_
#define _BRON_KERBOSCH_H_

#include "typedefs.h"

ISimplexList find_maximal_cliques(
    const PointIdList &vertices,
    const ConnectionList &edges);

#endif // _BRON_KERBOSCH_H_
