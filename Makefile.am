AUTOMAKE_OPTIONS = check-news
ACLOCAL_AMFLAGS = -I m4

bin_PROGRAMS = mincmorph

EXTRA_DIST = $(m4_files) \
	epm-header.in

m4_files = m4/mni_REQUIRE_LIB.m4 m4/mni_REQUIRE_MNILIBS.m4	 \
	m4/mni_REQUIRE_OPENINVENTOR.m4				 \
	m4/mni_cxx_have_koenig_lookup.m4 m4/smr_CGAL_MAKEFILE.m4 \
	m4/smr_OPTIONAL_LIB.m4 m4/smr_REQUIRED_LIB.m4		 \
	m4/smr_WITH_BUILD_PATH.m4

mincmorph_SOURCES = kernel_io.c kernel_ops.c mincmorph.c kernel_io.h kernel_ops.h
