set(LibraryVersion "1.0")
add_definitions(-DLIBRARY_VERSION="${LibraryVersion}")

compilelib("")

file(GLOB_RECURSE MAC "${CMAKE_CURRENT_SOURCE_DIR}/macros/*")
install(FILES ${MAC} DESTINATION ./macros/wimp)

file(GLOB_RECURSE MAC "${CMAKE_CURRENT_SOURCE_DIR}/examples/*")
install(FILES ${MAC} DESTINATION ./examples/wimp)

add_library_test()
