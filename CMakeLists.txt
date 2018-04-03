cmake_minimum_required (VERSION 3.1)

project (Neo_Fem_sc)

include_directories(
	"${PROJECT_SOURCE_DIR}/src"
 )

set (PTHREAD_INCLUDE "/usr/include" CACHE PATH "directory where pthread.h can be found")
set (PTHREAD_LIB_INCLUDE "/usr/lib" CACHE PATH "directory where the pthread library can be found")



# ############### ############### ############### ############### ############### ############### ############### ############### ###############
# Begin Pointer to Numerical_Recipes_source directory 

option(USING_Numerical_Recipes_source "Whether the Numerical_Recipes_source package will be linked with the library" ON)

if(USING_Numerical_Recipes_source)
	
	set(Numerical_Recipes_source_INCLUDE "/usr/local/include" CACHE PATH "directory to the source of GMM")
	add_definitions(-D_AUTODIFF)
	add_definitions(-D_USESTDVECTOR_)
    include_directories(${Numerical_Recipes_source_INCLUDE},${Numerical_Recipes_source_INCLUDE}/Numerical_Recipes_source)
endif()

# End Pointer to Numerical_Recipes_source directory 
# ############### ############### ############### ############### ############### ############### ############### ############### ###############


# ############### ############### ############### ############### ############### ############### ############### ############### ###############
# Begin Pointer to PZ directory 

option(USING_PZ "Whether the PZ package will be linked with the library" ON)

if(USING_PZ)
	set(CMAKE_CXX_STANDARD 11)
	set(CMAKE_CXX_STANDARD_REQUIRED ON)

	if (WIN32)
		find_package( PZ PATHS "C:/Arquivos de Programas/PZ" REQUIRED NO_DEFAULT_PATH)
	else()
		find_package( PZ REQUIRED )
	endif(WIN32)

	include_directories(${PZ_INCLUDE_DIRS})
    link_libraries(pz)
endif()

# End Pointer to PZ directory 
# ############### ############### ############### ############### ############### ############### ############### ############### ###############


# ############### ############### ############### ############### ############### ############### ############### ############### ###############
# Begin Pointer to Boost directory 


# ############### ############### ############### ############### ############### ############### ############### ############### ###############
# Begin Pointer to Armadillo library

option(USING_ARMADILLO "Whether the Armadillo library will be linked in" OFF)
if(USING_ARMADILLO)
	find_package( Armadillo REQUIRED NO_MODULE)
	include_directories(${ARMADILLO_INCLUDE_DIRS})
	link_libraries(${ARMADILLO_LIBRARIES})
endif()

# End Pointer to Armadillo directory 
# ############### ############### ############### ############### ############### ############### ############### ############### ###############


# ############### ############### ############### ############### ############### ############### ############### ############### ###############
# Begin Pointer to Eigen3 library

option(USING_EIGEN3 "Whether the Eigen3 library will be linked in" OFF)
if(USING_EIGEN3)
	find_package( Eigen3 REQUIRED )
	include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# End Pointer to Eigen3 directory 
# ############### ############### ############### ############### ############### ############### ############### ############### ###############



add_subdirectory(src)
add_subdirectory(Test)

enable_testing()

add_test (TestShapeRuns Test/femtest)