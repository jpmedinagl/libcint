add_executable(bug bug.c)
target_link_libraries(bug PUBLIC cint m LLDEnzymeFlags LLDEnzymePrintFlags LLDEnzymePrintTypeFlags)

# add_executable(buf ${CMAKE_SOURCE_DIR}/scf/buf.c)
# target_link_libraries(buf PUBLIC cint m LLDEnzymeFlags)

# add_executable(stv stv.c)
# target_link_libraries(stv PUBLIC cint m LLDEnzymeFlags LLDEnzymeLooseTypeFlags)

# add_executable(scf scf.c)
# target_link_libraries(scf PUBLIC cint m diag LLDEnzymeFlags LLDEnzymeLooseTypeFlags)
