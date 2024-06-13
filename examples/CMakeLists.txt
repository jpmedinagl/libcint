#add_executable(h2_ovlp h2_ovlp.c)
#add_executable(h2_kin h2_kin.c)
#add_executable(h2_nuc h2_nuc.c)

#add_executable(h2_ovlp_diff h2_ovlp_diff.c)
#add_executable(h2_kin_diff h2_kin_diff.c)

#add_executable(h2_stv h2_stv.c)
#add_executable(h2_rep h2_rep.c)
# add_executable(h2o_stv h2o_stv.c)
#add_executable(h2_stv_grad h2_stv_grad.c)
#add_executable(h2_rep_grad h2_rep_grad.c)

#target_link_libraries(h2_ovlp cint m)
#target_link_libraries(h2_kin cint m)
#target_link_libraries(h2_nuc cint m)

#target_link_libraries(h2_ovlp_diff PUBLIC cint m LLDEnzymeFlags)
#target_link_libraries(h2_kin_diff PUBLIC cint m LLDEnzymeFlags)
#target_link_libraries(h2_nuc_diff PUBLIC cint m LLDEnzymeFlags)

#target_link_libraries(h2_stv cint m)
#target_link_libraries(h2_rep cint m)
#target_link_libraries(h2o_stv cint m)
#target_link_libraries(h2_stv_grad PUBLIC cint m LLDEnzymeFlags)
#target_link_libraries(h2_rep_grad PUBLIC cint m LLDEnzymeFlags)

# add_subdirectory(diag)

# add_executable(scf scf.c)
# add_executable(scf_grad scf_grad.c)
# target_link_libraries(scf diag cint m)
# target_link_libraries(scf_grad PUBLIC cint m LLDEnzymeFlags diag)
