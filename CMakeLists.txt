cmake_minimum_required (VERSION 2.6)

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
# Begin Pointer to gmm directory 

option(USING_GMM "Whether the gmm package will be linked with the library" ON)

if(USING_GMM)
	
	set(GMM_INCLUDE "/usr/local/include" CACHE PATH "directory to the source of GMM")
	add_definitions(-D_AUTODIFF)
    include_directories(${GMM_INCLUDE},${GMM_INCLUDE}/gmm,${GMM_INCLUDE}/gmm/gmm-4.1)
endif()

# End Pointer to gmm directory 
# ############### ############### ############### ############### ############### ############### ############### ############### ###############



# ############### ############### ############### ############### ############### ############### ############### ############### ###############
# Begin Pointer to tinyxml directory 

option(USING_tinyxml "Whether the tinyxml package will be linked with the library" ON)

if(USING_tinyxml)
	
	set(tinyxml_INCLUDE "/usr/local/include" CACHE PATH "directory to the source of tinyxml")
	add_definitions(-D_AUTODIFF)
    include_directories(${tinyxml_INCLUDE},${tinyxml_INCLUDE}/tinyxml)
endif()

# End Pointer to tinyxml directory 
# ############### ############### ############### ############### ############### ############### ############### ############### ###############



# ############### ############### ############### ############### ############### ############### ############### ############### ###############
# Begin Pointer to PZ directory 

option(USING_PZ "Whether the PZ package will be linked with the library" ON)

if(USING_PZ)
    
    set(PZ_INCLUDE "/usr/local/include" CACHE PATH "directory to the source of PZ")
    add_definitions()
    link_directories(${PZ_INCLUDE}/lib)
    include_directories(${PZ_INCLUDE}/include/pz)
    link_libraries(pzmatrix)

endif()

# End Pointer to PZ directory 
# ############### ############### ############### ############### ############### ############### ############### ############### ###############



# ############### ############### ############### ############### ############### ############### ############### ############### ###############
# Begin Pointer to Boost directory 

option(USING_BOOST "Whether the Boost package will be linked with the library" ON)

if(USING_BOOST)
	
	find_package( Boost 1.36.0 )
   		if(Boost_FOUND)
			add_definitions(-DUSING_BOOST)
      		include_directories(${Boost_INCLUDE_DIRS})
   		endif()
endif()

# End Pointer to Boost directory 
# ############### ############### ############### ############### ############### ############### ############### ############### ###############



# ############### ############### ############### ############### ############### ############### ############### ############### ###############
# Begin Pointer to Fad directory

option(USING_FAD "Whether the Boost package will be linked with the library" ON)

if(USING_FAD)
	
	set(FAD_INCLUDE "/usr/local/include" CACHE PATH "directory to the source of FAD")
	add_definitions(-D_AUTODIFF)
    include_directories(${FAD_INCLUDE},${FAD_INCLUDE}/Fad,${FAD_INCLUDE}/TinyFad,${FAD_INCLUDE}/TinyFadET)
endif()

# End Pointer to Fad directory 
# ############### ############### ############### ############### ############### ############### ############### ############### ###############



# ############### ############### ############### ############### ############### ############### ############### ############### ###############
# Begin Pointer to Log4cxx directory

option(USING_LOG4CXX "Whether the LOG4CXX library will be linked in" OFF)
if(USING_LOG4CXX)
	set (LOG4CXX_INCLUDE "/usr/local/include" CACHE PATH "directory where log4cxx header files can be found")
	set (LOG4CXX_LIB_INCLUDE "/usr/local/lib" CACHE PATH "directory where the log4cxx library can be found")
	include_directories(${LOG4CXX_INCLUDE})
	link_directories(${LOG4CXX_LIB_INCLUDE})
	#add_definitions(-DLOG4CXX -DPZSOURCEDIR=\"${PROJECT_SOURCE_DIR}/NeoPZ\")
	link_libraries(log4cxx)
endif()

# End Pointer to Log4cxx directory 
# ############### ############### ############### ############### ############### ############### ############### ############### ###############


add_subdirectory(src)
add_subdirectory(Test)

enable_testing()

add_test (TestShapeRuns Test/femtest)