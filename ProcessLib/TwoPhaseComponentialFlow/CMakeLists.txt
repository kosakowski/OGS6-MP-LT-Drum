append_source_files(SOURCES)

add_library(TwoPhaseComponentialFlow ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS TwoPhaseComponentialFlow
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(TwoPhaseComponentialFlow PUBLIC ProcessLib PRIVATE ParameterLib)

include(Tests.cmake)
