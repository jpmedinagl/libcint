add_executable(buf buf.c)
target_link_libraries(buf PUBLIC cint m) # LLDEnzymeFlags)

add_executable(stv stv.c)
target_link_libraries(stv PUBLIC cint m LLDEnzymeFlags LLDEnzymeLooseTypeFlags)

add_executable(scf scf.c)
target_link_libraries(scf PUBLIC cint m diag LLDEnzymeFlags LLDEnzymeLooseTypeFlags)
