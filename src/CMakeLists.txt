

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)

## change another change and another going on

set(femsc_STAT_SRCS
	${headers}
	${sources}
	)

add_executable(Fem_sc ${femsc_STAT_SRCS})

