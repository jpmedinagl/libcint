add_executable(scf scf.c)
target_link_libraries(scf PUBLIC cint m LLDEnzymeFlags diag LLDEnzymeLooseTypeFlags)