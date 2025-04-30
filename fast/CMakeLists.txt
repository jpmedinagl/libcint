add_executable(tcache tcache.c)
target_link_libraries(tcache PUBLIC cint m quadmath)# diag LLDEnzymeFlags LLDEnzymeLooseTypeFlags)

add_executable(nocache nocache.c)
target_link_libraries(nocache PUBLIC cint m quadmath)# diag LLDEnzymeFlags LLDEnzymeLooseTypeFlags)