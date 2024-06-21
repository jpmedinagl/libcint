# add_executable(eig eig.c)
# target_link_libraries(eig diag cint m)

# add_executable(h2_stv h2_stv.c)
# add_executable(h2_rep h2_rep.c)
# add_executable(h2_scf h2_scf.c)

# target_link_libraries(h2_stv cint m)
# target_link_libraries(h2_rep cint m)
# target_link_libraries(h2_scf diag cint m)

add_executable(h2o_scf h2o_scf.c)
target_link_libraries(h2o_scf diag cint m)


# add_executable(h2_stv_grad h2_stv_grad.c)
# add_executable(h2_rep_grad h2_rep_grad.c)
# add_executable(h2_scf_grad h2_scf_grad.c)

# target_link_libraries(h2_stv_grad PUBLIC cint m LLDEnzymeFlags LLDEnzymeLooseTypeFlags)
# target_link_libraries(h2_rep_grad PUBLIC cint m LLDEnzymeFlags LLDEnzymeLooseTypeFlags)
# target_link_libraries(h2_scf_grad PUBLIC cint m LLDEnzymeFlags diag LLDEnzymeLooseTypeFlags)


# add_executable(h2o_stv h2o_stv.c)
# target_link_libraries(h2o_stv cint m)