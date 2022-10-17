#pragma once

#include <crab/domains/apron_domains.hpp>
#include <crab/domains/elina_domains.hpp>
#include "crab_defs.hh"

namespace clam {
#ifdef HAVE_APRON
using BASE(tvpi_domain_t) =
  crab::domains::apron_domain<number_t, region_subdom_varname_t,
			      crab::domains::apron_domain_id_t::APRON_TVPI>;
#else
using BASE(tvpi_domain_t) =
  crab::domains::elina_domain<number_t, region_subdom_varname_t,
			      crab::domains::elina_domain_id_t::ELINA_TVPI>;
#endif
using tvpi_domain_t = RGN_FUN(ARRAY_FUN(BOOL_NUM(BASE(tvpi_domain_t))));
} // end namespace clam
