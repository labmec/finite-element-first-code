########### next target ###############

set(main_SRCS Testmain.cpp)

add_executable(femtest ${testefem} ${main_SRCS})

target_link_libraries(femtest pthread boost_unit_test_framework)


