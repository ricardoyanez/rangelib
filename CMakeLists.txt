cmake_minimum_required(VERSION 3.10)

# project name
project(range VERSION 0.2.0)

set(CMAKE_C_FLAGS "-g -Wall -O2")

include(GNUInstallDirs)

# shared library
add_library(${PROJECT_NAME}-lib SHARED src/rangelib.c src/ranges.c src/nr.c)
set_target_properties(${PROJECT_NAME}-lib PROPERTIES
	VERSION ${PROJECT_VERSION}
	SOVERSION 1
	PUBLIC_HEADER src/range.h
	OUTPUT_NAME ${PROJECT_NAME})

# install library
install(TARGETS ${PROJECT_NAME}-lib
	LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# executable
add_executable(${PROJECT_NAME}-bin src/range.c)
target_link_libraries(${PROJECT_NAME}-bin ${PROJECT_NAME}-lib m)
set_target_properties(${PROJECT_NAME}-bin PROPERTIES
	OUTPUT_NAME ${PROJECT_NAME})

# install executable
install(TARGETS ${PROJECT_NAME}-bin)

# install man pages
install(FILES man/range.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
install(FILES man/rangelib.3 DESTINATION ${CMAKE_INSTALL_MANDIR}/man3)

# install examples
install(DIRECTORY examples DESTINATION ${CMAKE_INSTALL_DOCDIR})

# install docs
install(FILES changelog DESTINATION ${CMAKE_INSTALL_DOCDIR})
