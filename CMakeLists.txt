set(LibraryVersion "1.0")
add_definitions(-DLIBRARY_VERSION="${LibraryVersion}")

COMPILELIB("")

file(GLOB_RECURSE MAC "${CMAKE_CURRENT_SOURCE_DIR}/macros/*")
INSTALL(FILES ${MAC} DESTINATION ./macros/wimp)

file(GLOB_RECURSE MAC "${CMAKE_CURRENT_SOURCE_DIR}/examples/*")
INSTALL(FILES ${MAC} DESTINATION ./examples/wimp)

ADD_LIBRARY_TEST()
