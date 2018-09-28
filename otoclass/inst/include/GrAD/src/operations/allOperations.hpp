// This file is part of GrAD, a C++ template library for
// graph based reverse mode automatic differentiation in C++
//
// Copyright (C), Christoffer Moesgaard Albertsen <cmoe@aqua.dtu.dk>

#include "ADparameter.hpp"
#include "ADconstant.hpp"

// Arithmetic operations
#include "arithmetic/ADprod.hpp"
#include "arithmetic/ADdiv.hpp"
#include "arithmetic/ADsum.hpp"
#include "arithmetic/ADsubtract.hpp"
#include "arithmetic/ADusubtract.hpp"

// Unary math
#include "unary_math/ADsign.hpp"
#include "unary_math/ADabs.hpp" // Uses sign
#include "unary_math/ADfabs.hpp" // Uses sign
#include "unary_math/ADsqrt.hpp"
#include "unary_math/ADacos.hpp" // Uses sqrt
#include "unary_math/ADasin.hpp"  // Uses sqrt
#include "unary_math/ADatan.hpp"
#include "unary_math/ADcos.hpp"
#include "unary_math/ADcosh.hpp"
#include "unary_math/ADexp.hpp"
#include "unary_math/ADlog.hpp"
#include "unary_math/ADlog10.hpp"
#include "unary_math/ADsin.hpp"
#include "unary_math/ADsinh.hpp"
#include "unary_math/ADtan.hpp"
#include "unary_math/ADtanh.hpp"

#ifdef __cplusplus
#if __cplusplus >= 201103L
#include "unary_math/ADacosh.hpp" // Uses sqrt
#include "unary_math/ADasinh.hpp" // Uses sqrt
#include "unary_math/ADatanh.hpp"
#include "unary_math/ADexpm1.hpp"
#include "unary_math/ADerf.hpp"	// Uses exp
#include "unary_math/ADlog1p.hpp"
#include "unary_math/ADexp2.hpp"
#include "unary_math/ADlog2.hpp"
#endif
#endif


// Binary math
#include "binary_math/ADatan2.hpp"
#include "binary_math/ADpow.hpp"

// // R math
