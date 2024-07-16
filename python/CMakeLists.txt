# add_executable(stv stvgrad.c)
# target_link_libraries(stv PUBLIC cint m LLDEnzymeFlags diag LLDEnzymeLooseTypeFlags)

add_executable(scfcp scf_cp.c)
target_link_libraries(scfcp PUBLIC cint m LLDEnzymeFlags diag LLDEnzymeLooseTypeFlags)

# add_executable(scfnew scf_new.c)
# target_link_libraries(scfnew PUBLIC cint m LLDEnzymeFlags diag LLDEnzymeLooseTypeFlags)

# add_executable(scf scf.c)
# target_link_libraries(scf PUBLIC cint m LLDEnzymeFlags diag LLDEnzymeLooseTypeFlags)
