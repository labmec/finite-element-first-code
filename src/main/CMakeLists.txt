

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

########### next target ###############

set(femsc_STAT_SRCS
	${headers}
	${sources}
	${femsc_STAT_SRCS}
	${tinyxml_STAT_SRCS}
#	PARENT_SCOPE
	)


add_executable(Fem_sc ${femsc_STAT_SRCS})

target_link_libraries(Fem_sc pz pthread)



#original Makefile.am contents follow:

##
## NeoPZ StrMatrix makefile
##
## It will be processed by GNU automake to produce Makefile.in
##
## $Id#
##
#
#DISTCLEANFILES = @PZ_DISTCLEANFILES@
#CLEANFILES     = @PZ_CLEANFILES@
#
#-include $(top_builddir)/@PZ_MK@
#
#strmatrixdir = $(prefix)/StrMatrix
#
#strmatrixincdir = $(includedir)/StrMatrix
#strmatrixinc_HEADERS= $(top_srcdir)/StrMatrix/*.h
#
#noinst_LIBRARIES = libpzstrmatrix.a
#
#libpzstrmatrix_a_SOURCES = TPBSpStructMatrix.cpp TPZFrontStructMatrix.cpp \
#	TPZMTAssemble.cpp TPZParFrontStructMatrix.cpp TPZParSkylineStructMatrix.cpp \
#	TPZSpStructMatrix.cpp pzbdstrmatrix.cpp pzbstrmatrix.cpp pzfstrmatrix.cpp pzmaterialids \
#	pzsbstrmatrix.cpp pzskylstrmatrix.cpp pzstrmatrix.cpp \
#	tpzsparseblockdiagonalstructmatrix.cpp TPZMTAssemble.cpp
#
### --| NeoPZ |------------------------------------------------------------------
#noinst_HEADERS = TPZMTAssemble.h
