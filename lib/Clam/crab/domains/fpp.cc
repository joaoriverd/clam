#include <clam/config.h>
#include <clam/CrabDomain.hh>
#include <clam/RegisterAnalysis.hh>
#include <crab/config.h>
#include "fpp.hh"

namespace clam {
#ifdef INCLUDE_ALL_DOMAINS
#if defined(HAVE_APRON) || defined(HAVE_ELINA)
REGISTER_DOMAIN(clam::CrabDomain::FPP, fpp_domain)
#else
UNREGISTER_DOMAIN(fpp_domain)
#endif
#else
UNREGISTER_DOMAIN(fpp_domain)
#endif  
} // end namespace clam

