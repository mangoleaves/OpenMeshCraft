#include "OBJReader.h"
#include "OBJWriter.h"
#include "STLReader.h"
#include "STLWriter.h"

namespace OMC {

template class OBJReader<TriSoupTraits>;
template class OBJWriter<TriSoupTraits>;
template class STLReader<TriSoupTraits>;
template class STLWriter<TriSoupTraits>;

} // namespace OMC