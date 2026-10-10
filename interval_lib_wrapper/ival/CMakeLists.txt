# The ival interval library (crnumerics' ival): tight IEEE 1788.1 intervals, CORE-MATH's correctly rounded functions,
# on x86-64 and aarch64; version 0.2.0 or later (ival1_scale). Found through pkg-config (ival.pc); IVAL_ROOT may
# name its install prefix.
set (IVAL_NAME "ival" PARENT_SCOPE)
set (IVAL_TARGET "ival_wrapper" PARENT_SCOPE)

if (IVAL_ROOT)
  set (ENV{PKG_CONFIG_PATH} "${IVAL_ROOT}/lib/pkgconfig:${IVAL_ROOT}/lib64/pkgconfig:$ENV{PKG_CONFIG_PATH}")
endif ()
find_package (PkgConfig REQUIRED)
pkg_check_modules (LIBIVAL REQUIRED ival>=0.2.0)
set (IVAL_VERSION ${LIBIVAL_VERSION} PARENT_SCOPE)
message (STATUS "ival ${LIBIVAL_VERSION}: ${LIBIVAL_LINK_LIBRARIES}")

add_library (ival_wrapper INTERFACE)
target_include_directories (ival_wrapper INTERFACE ${LIBIVAL_INCLUDE_DIRS})
target_link_libraries (ival_wrapper INTERFACE ${LIBIVAL_LINK_LIBRARIES})
install (TARGETS ival_wrapper EXPORT ibexExport COMPONENT devtools)
