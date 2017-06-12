set(dynmetidcmd_SRC main.cpp
                    lcmsannotate.cpp
                    stroperation.cpp)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(dynmetid ${dynmetidcmd_SRC})
target_link_libraries(dynmetid m ${MYSQL_LIBRARIES})
install(TARGETS dynmetid DESTINATION bin COMPONENT ApplicationData)
