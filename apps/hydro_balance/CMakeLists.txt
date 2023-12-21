add_executable(hydro_balance hydro_balance.cpp)
target_link_libraries(hydro_balance PRIVATE nebulasem)
install(TARGETS hydro_balance)
