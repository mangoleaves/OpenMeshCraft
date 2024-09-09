// Triangle soup for input and output
#include "OpenMeshCraft/Mesh/TriSoup.h"
// Kernel
#include "OpenMeshCraft/Geometry/ExactIndirectPredicatesApproxConstructions.h"

#include "MeshBoolean.h"

namespace OMC {

template class MeshBoolean<EIAC, TriSoupTraits>;

} // namespace OMC