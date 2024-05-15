#add_executable(h2o_ovlp h2o_ovlp.c)
add_executable(h2_ovlp h2_ovlp.c)
add_executable(h2_ovlp_diff h2_ovlp_diff.c)
#add_executable(diff diff.c)
#add_executable(hello_world hello_world.c)

#target_compile_definitions(hello_world PUBLIC -DNOVERSION)

#target_link_libraries(h2o_ovlp cint m)
target_link_libraries(h2_ovlp cint m)
#target_link_libraries(h2_ovlp_diff cint m)
#target_link_libraries(diff cint m)
target_link_libraries(h2_ovlp_diff PUBLIC cint m LLDEnzymeFlags)
#target_link_libraries(hello_world PUBLIC cint m LLDEnzymeFlags)
