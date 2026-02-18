set (SOURCES
  main.cpp 
  Stopwatch.h
  Hydroinfo_MUSIC.cpp
  ParameterReader.cpp
  arsenal.cpp
  cornelius.cpp
  SurfaceFinder.cpp
  )

add_library(${libname} SHARED ${SOURCES})
add_executable (${exename} main.cpp ${SOURCES})

install(TARGETS ${exename} DESTINATION ${CMAKE_HOME_DIRECTORY})
