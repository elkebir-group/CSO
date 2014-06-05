cmake_minimum_required( VERSION 2.8 )

project( cso )

set( Analysis_Hdr
  src/analysis/linkageanalysis.h
  src/analysis/lowerbounds.h
)

set( Analysis_Src
  #  src/analysis/analysemain.cpp
  src/analysis/linkageanalysis.cpp
  src/analysis/lowerbounds.cpp
)

set( Common_Hdr
  src/common/cso.h
  src/common/data.h
  src/common/genotype.h
  src/common/genotypegamete.h
  src/common/stopwatch.h
)

set( Common_Src
  src/common/cso.cpp
  src/common/data.cpp
  src/common/genotype.cpp
  src/common/stopwatch.cpp
)

set( Generator_Src
  src/generator/inputgen.cpp
)

set( Tinyxml_Hdr
  src/externlib/tinyxml/tinystr.h
  src/externlib/tinyxml/tinyxml.h
)

set( Tinyxml_Src
  src/externlib/tinyxml/tinystr.cpp
  src/externlib/tinyxml/tinyxml.cpp
  src/externlib/tinyxml/tinyxmlerror.cpp
  src/externlib/tinyxml/tinyxmlparser.cpp
)

set( Cso_Hdr
  src/ilpsolverheuristic.h
  src/crossingschedule.h
  src/ilpsolver.h
  src/ilpsolverext.h
)

set( Cso_Src
  src/ilpsolverheuristic.cpp
  src/crossingschedule.cpp
  src/ilpsolver.cpp
  src/ilpsolverext.cpp
  src/main.cpp
)

set( GitCommand1 "symbolic-ref HEAD 2> /dev/null | cut -b 12-" )
execute_process( COMMAND git symbolic-ref HEAD  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                 COMMAND tail -c -7 COMMAND tr -d "\\n" OUTPUT_VARIABLE GitBranch )
set( GitCommand2 "log --pretty=format:\\\"%h\\\" -1" )
execute_process( COMMAND git log --pretty=format:%h WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
                 COMMAND head -n 1  COMMAND tr -d "\\n" OUTPUT_VARIABLE GitHashtag )
set( CSO_VERSION \"${GitBranch}-${GitHashtag}\" )

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
	set( CPLEX_AUTODETECTION_PREFIXES "$ENV{HOME}/ILOG;/opt/ibm/ILOG" CACHE STRING "CPLEX automatic lookup directories" )
	foreach( CplexInstallPrefix ${CPLEX_AUTODETECTION_PREFIXES} )
		file( GLOB_RECURSE CplexLibraryPaths ${CplexInstallPrefix}/*/cplex/lib/*/libcplex.a )
		if( CplexLibraryPaths )
			list( LENGTH CplexLibraryPaths NumPath )
			if( NumPath GREATER 1 )
				message( WARNING "Found multiple CPLEX installations in ${CplexInstallPrefix}, automatically selecting the first one, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually if this is not the right selection" )
			else()
				message( "Found CPLEX's installation in ${CPlexInstallPrefix}" )
			endif()
			list( GET CplexLibraryPaths 0 CplexLibraryPath )

			set( MyRegex "${CplexInstallPrefix}/(.*)/cplex/?/lib/(.*)/libcplex.a" )
			string( REGEX MATCH "${MyRegex}" ValidatedCplexLibPath ${CplexLibraryPath} )
			if( NOT ValidatedCplexLibPath )
				message( WARNING "Unable to recognize ${CplexLibraryPath} as a valid CPLEX path" )
			else()
				string( REGEX REPLACE "${MyRegex}" "\\1" CplexHome ${ValidatedCplexLibPath} )
				string( REGEX REPLACE "${MyRegex}" "\\2" CplexArch ${ValidatedCplexLibPath} )

				set( CPLEX_INC_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/include"            CACHE PATH "CPLEX's include dir" )
				set( CPLEX_LIB_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/lib/${CplexArch}"   CACHE PATH "CPLEX's library dir" )
				set( CONCERT_INC_DIR "${CplexInstallPrefix}/${CplexHome}/concert/include"          CACHE PATH "Concert's include dir" )
        set( CONCERT_LIB_DIR "${CplexInstallPrefix}/${CplexHome}/concert/lib/${CplexArch}" CACHE PATH "Concert's library dir" )

				break()
			endif()
		else()
			message( "Unable to find CPLEX and Concert installation in ${CplexInstallPrefix}" )
		endif()
	endforeach()

	if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
		message( "Unable to autodetect CPLEX's and Concert's installation directories, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually" )
	else()
		message( "Autodetected CPLEX's include dir   : ${CPLEX_INC_DIR}" )
		message( "Autodetected CPLEX's library dir   : ${CPLEX_LIB_DIR}" )
		message( "Autodetected Concert's include dir : ${CONCERT_INC_DIR}" )
		message( "Autodetected Concert's library dir : ${CONCERT_LIB_DIR}" )
	endif()
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -g -ggdb -O0 -DDEBUG -DVERBOSE")
set(CMAKE_CXX_FLAGS_RELEASE "-Wall -Wextra -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_PROFILE "-Wall -Wextra -Wno-long-long -Wno-unused-parameter -Wno-unknown-pragmas -O3 -g -ggdb -DNDEBUG")

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings " )
add_definitions( -DIL_STD )
configure_file( src/config.h.in ${PROJECT_SOURCE_DIR}/src/config.h )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++ " )
endif()

SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libstdc++")

include_directories( ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} "${LIBLEMON_ROOT}/include" "src" "src/externlib" )
link_directories( ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} "${LIBLEMON_ROOT}/lib" )

set( CommonLibs
	emon
)

set( CplexLibs
  m
  ilocplex
  cplex
  concert
  pthread
)

add_executable( csoilp ${Analysis_Hdr} ${Analysis_Src} ${Common_Hdr} ${Common_Src} ${CommonObjs} ${Tinyxml_Hdr} ${Tinyxml_Src} ${Cso_Hdr} ${Cso_Src} )
target_link_libraries( csoilp ${CommonLibs} ${CplexLibs} )
