set(dynmetidcmd_SRC main.cpp
                    lcmsannotate.cpp
                    stroperation.cpp)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_executable(mysqldynmetid ${dynmetidcmd_SRC})
target_link_libraries(mysqldynmetid m ${MYSQL_LIBRARIES})
install(TARGETS mysqldynmetid DESTINATION bin COMPONENT ApplicationData)
