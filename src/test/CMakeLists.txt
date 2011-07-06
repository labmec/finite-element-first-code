
file(GLOB headers *.h)
file(GLOB sources *.cpp)

########### next target ###############
	set(testefem
		${femsc_STAT_SRCS}
		${testefem}
		${headers}
		${sources}
		PARENT_SCOPE
		)


