// The purpose of this file is to instantiate templates
// This is a workaround template instantiation in g++. It might
// need some editing if using a different compiler version.
// This file should work well with versions 2.8.1 to 2.95.1

#include "CircuitManager.h"

// Vectors
// template class vector<string>;

template class std::vector<unsigned>;

template class std::vector<GraphNode*>;
template class std::vector<Terminal*>;
template class std::vector<Element*>;
template class std::vector<Circuit*>;

// Maps
#if __GNUC_MINOR__ == 95

// only for gcc 2.95.x ?
template class _Rb_tree<unsigned int, pair<unsigned int const, Element *>,
_Select1st<pair<unsigned int const, Element *> >, ltunsigned, allocator<Element *> >;
template class _Rb_tree<unsigned int, pair<unsigned int const, Terminal *>,
_Select1st<pair<unsigned int const, Terminal *> >, ltunsigned, allocator<Terminal *> >;

#else

template class std::map<unsigned int, Terminal*, ltunsigned>;
template class std::map<unsigned int, Element*, ltunsigned>;

#endif

