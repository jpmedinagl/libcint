add_executable(buf buf.c)
target_link_libraries(buf PUBLIC cint m LLDEnzymePrintFlags) # LLDEnzymeFlags)

# add_executable(stv stv.c)
# target_link_libraries(stv PUBLIC cint m LLDEnzymeFlags LLDEnzymeLooseTypeFlags)

# add_executable(scfnew scf_new.c)
# target_link_libraries(scfnew PUBLIC cint m LLDEnzymeFlags diag LLDEnzymeLooseTypeFlags)

# add_executable(scfcp scf_cp.c)
# target_link_libraries(scfcp PUBLIC cint m LLDEnzymeFlags diag LLDEnzymeLooseTypeFlags)

# add_executable(scf scf.c)
# target_link_libraries(scf PUBLIC cint m LLDEnzymeFlags diag LLDEnzymeLooseTypeFlags)
