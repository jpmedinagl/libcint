# add_executable(buf buf.c)
# target_link_libraries(buf PUBLIC cint m LLDEnzymePrintFlags) # LLDEnzymeFlags)

# add_executable(stv stv.c)
# target_link_libraries(stv PUBLIC cint m LLDEnzymeFlags LLDEnzymeLooseTypeFlags)

add_executable(scfp scf.c)
target_link_libraries(scfp PUBLIC cint m LLDEnzymeFlags diag LLDEnzymeLooseTypeFlags)
