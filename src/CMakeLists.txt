

file(GLOB headers *.h)
install(FILES ${headers} DESTINATION include)

file(GLOB sources *.cpp)
file(GLOB input input.txt)

## change another change and another going on

set(femsc_STAT_SRCS
	${headers}
	${sources}
	${input}
	)

add_subdirectory(test)

set(testefem
${testefem}
PARENT_SCOPE
)

add_subdirectory(main)

