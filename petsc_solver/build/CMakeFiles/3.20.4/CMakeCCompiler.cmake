set(CMAKE_C_COMPILER "/es01/paratera/parasoft/oneAPI/2022.1/mpi/2021.5.0/bin/mpiicc")
set(CMAKE_C_COMPILER_ARG1 "")
set(CMAKE_C_COMPILER_ID "Intel")
set(CMAKE_C_COMPILER_VERSION "2021.5.0.20211109")
set(CMAKE_C_COMPILER_VERSION_INTERNAL "")
set(CMAKE_C_COMPILER_WRAPPER "")
set(CMAKE_C_STANDARD_COMPUTED_DEFAULT "11")
set(CMAKE_C_COMPILE_FEATURES "c_std_90;c_function_prototypes;c_std_99;c_restrict;c_variadic_macros;c_std_11;c_static_assert")
set(CMAKE_C90_COMPILE_FEATURES "c_std_90;c_function_prototypes")
set(CMAKE_C99_COMPILE_FEATURES "c_std_99;c_restrict;c_variadic_macros")
set(CMAKE_C11_COMPILE_FEATURES "c_std_11;c_static_assert")

set(CMAKE_C_PLATFORM_ID "Linux")
set(CMAKE_C_SIMULATE_ID "GNU")
set(CMAKE_C_COMPILER_FRONTEND_VARIANT "")
set(CMAKE_C_SIMULATE_VERSION "7.3.0")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_C_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_C_COMPILER_RANLIB "")
set(CMAKE_LINKER "/usr/bin/ld")
set(CMAKE_MT "")
set(CMAKE_COMPILER_IS_GNUCC )
set(CMAKE_C_COMPILER_LOADED 1)
set(CMAKE_C_COMPILER_WORKS TRUE)
set(CMAKE_C_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_C_COMPILER_ENV_VAR "CC")

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_C_COMPILER_ID_RUN 1)
set(CMAKE_C_SOURCE_FILE_EXTENSIONS c;m)
set(CMAKE_C_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_C_LINKER_PREFERENCE 10)

# Save compiler ABI information.
set(CMAKE_C_SIZEOF_DATA_PTR "8")
set(CMAKE_C_COMPILER_ABI "ELF")
set(CMAKE_C_BYTE_ORDER "LITTLE_ENDIAN")
set(CMAKE_C_LIBRARY_ARCHITECTURE "")

if(CMAKE_C_SIZEOF_DATA_PTR)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_C_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_C_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_C_COMPILER_ABI}")
endif()

if(CMAKE_C_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()

set(CMAKE_C_CL_SHOWINCLUDES_PREFIX "")
if(CMAKE_C_CL_SHOWINCLUDES_PREFIX)
  set(CMAKE_CL_SHOWINCLUDES_PREFIX "${CMAKE_C_CL_SHOWINCLUDES_PREFIX}")
endif()





set(CMAKE_C_IMPLICIT_INCLUDE_DIRECTORIES "/es01/paratera/parasoft/oneAPI/2022.1/mpi/2021.5.0/include;/es01/paratera/sce0588/zl/ASC/lapack-3.11/CBLAS/include;/es01/paratera/sce0588/zl/ASC/lapack-3.11/LAPACKE/include;/es01/paratera/parasoft/soft/7.3.0-para/include;/es01/paratera/parasoft/oneAPI/2022.1/vpl/2022.0.0/include;/es01/paratera/parasoft/oneAPI/2022.1/tbb/2021.5.0/include;/es01/paratera/parasoft/oneAPI/2022.1/mkl/2022.0.0/include;/es01/paratera/parasoft/oneAPI/2022.1/ipp/2021.5.0/include;/es01/paratera/parasoft/oneAPI/2022.1/ippcp/2021.5.0/include;/es01/paratera/parasoft/oneAPI/2022.1/dpl/2021.6.0/linux/include;/es01/paratera/parasoft/oneAPI/2022.1/dpcpp-ct/2022.0.0/include;/es01/paratera/parasoft/oneAPI/2022.1/dnnl/2022.0.0/cpu_dpcpp_gpu_dpcpp/lib;/es01/paratera/parasoft/oneAPI/2022.1/dev-utilities/2021.5.0/include;/es01/paratera/parasoft/oneAPI/2022.1/dal/2021.5.0/include;/es01/paratera/parasoft/oneAPI/2022.1/ccl/2021.5.0/include/cpu_gpu_dpcpp;/es01/paratera/parasoft/oneAPI/2022.1/compiler/2022.0.0/linux/compiler/include/intel64;/es01/paratera/parasoft/oneAPI/2022.1/compiler/2022.0.0/linux/compiler/include/icc;/es01/paratera/parasoft/oneAPI/2022.1/compiler/2022.0.0/linux/compiler/include;/usr/local/include;/es01/paratera/parasoft/soft/7.3.0-para/lib/gcc/x86_64-pc-linux-gnu/7.3.0/include;/es01/paratera/parasoft/soft/7.3.0-para/lib/gcc/x86_64-pc-linux-gnu/7.3.0/include-fixed;/usr/include")
set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "mpifort;mpi;dl;rt;pthread;imf;svml;irng;m;ipgo;decimal;cilkrts;stdc++;gcc;gcc_s;irc;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_C_IMPLICIT_LINK_DIRECTORIES "/es01/paratera/parasoft/oneAPI/2022.1/mpi/2021.5.0/lib/release;/es01/paratera/parasoft/oneAPI/2022.1/mpi/2021.5.0/lib;/es01/paratera/parasoft/soft/7.3.0-para/lib64;/es01/paratera/parasoft/soft/7.3.0-para/lib;/es01/paratera/parasoft/oneAPI/2022.1/vpl/2022.0.0/lib;/es01/paratera/parasoft/oneAPI/2022.1/tbb/2021.5.0/lib/intel64/gcc4.8;/es01/paratera/parasoft/oneAPI/2022.1/mpi/2021.5.0/libfabric/lib;/es01/paratera/parasoft/oneAPI/2022.1/mkl/2022.0.0/lib/intel64;/es01/paratera/parasoft/oneAPI/2022.1/ipp/2021.5.0/lib/intel64;/es01/paratera/parasoft/oneAPI/2022.1/ippcp/2021.5.0/lib/intel64;/es01/paratera/parasoft/oneAPI/2022.1/dnnl/2022.0.0/cpu_dpcpp_gpu_dpcpp/lib;/es01/paratera/parasoft/oneAPI/2022.1/dal/2021.5.0/lib/intel64;/es01/paratera/parasoft/oneAPI/2022.1/compiler/2022.0.0/linux/compiler/lib/intel64_lin;/es01/paratera/parasoft/oneAPI/2022.1/compiler/2022.0.0/linux/lib;/es01/paratera/parasoft/oneAPI/2022.1/clck/2021.5.0/lib/intel64;/es01/paratera/parasoft/oneAPI/2022.1/ccl/2021.5.0/lib/cpu_gpu_dpcpp;/es01/paratera/parasoft/soft/7.3.0-para/lib/gcc/x86_64-pc-linux-gnu/7.3.0;/lib64;/usr/lib64;/lib;/usr/lib")
set(CMAKE_C_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
