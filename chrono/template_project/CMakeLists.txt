#--------------------------------------------------------------
# CMake configuration script for linking an external project to
# Chrono libraries.
# For a particular user project, modify sections 1 - 4 below,
# as appropriate.
#--------------------------------------------------------------
 
cmake_minimum_required(VERSION 3.18)
cmake_policy(SET CMP0091 NEW)

#--------------------------------------------------------------
# 1. Set the project name
#--------------------------------------------------------------

project(my_project)

#--------------------------------------------------------------
# 2. Find REQUIRED and OPTIONAL Chrono modules
#
# Invoke the find_package function in CONFIG mode:
#    find_package(Chrono
#                 COMPONENTS req_module1 req_module1 ...
#                 OPTIONAL_COMPONENTS opt_module1 opt_module2 ...
#                 CONFIG)
# This hard-trimmed fork exports the core Chrono multibody framework only.
# 
# Note that you will have to set the variable Chrono_DIR to 
# specify the location of the chrono-config.cmake script, if
# it is not in a default install location. Chrono_DIR can point
# to either a Chrono build tree or a Chrono install tree.
# 
# The following variables are set and can be used further down:
# Chrono_FOUND
#       set to true if Chrono and all required components were found
# CHRONO_STATIC
#       set to ON if Chrono has been built as static library, OFF if shared library
# CHRONO_TARGETS
#       list of exported Chrono targets
# CHRONO_DLL_NAMES
#       list of Chrono DLLs (without path)
# CHRONO_DATA_DIR
#       path to the Chrono data directory
#
# In addition, for each requested component [COMPONENT], the following variable
# is set to true (ON) or false (OFF):
# CHRONO_[COMPONENT]_FOUND
# 
# In this trimmed fork, request only the core Chrono package.
#--------------------------------------------------------------

find_package(Chrono CONFIG)

# Return now if Chrono or a required component was not found.
if (NOT Chrono_FOUND)
  message("Could not find Chrono or one of its required modules")
  return()
endif()

#--------------------------------------------------------------
# 3. Specify project sources and add executable
#--------------------------------------------------------------

set(MY_FILES my_example.cpp)

add_executable(my_demo ${MY_FILES})

#--------------------------------------------------------------
# Set properties for the executable target
#--------------------------------------------------------------

# Here, we define a macro CHRONO_DATA_DIR which will contain the
# path to the Chrono data directory, either in its source tree
# (if using a build version of Chrono), or in its install tree
# (if using an installed version of Chrono).
target_compile_definitions(my_demo PRIVATE "CHRONO_DATA_DIR=\"${CHRONO_DATA_DIR}\"") 

if(MSVC)
    set_target_properties(my_demo PROPERTIES MSVC_RUNTIME_LIBRARY ${CHRONO_MSVC_RUNTIME_LIBRARY})
endif()

#--------------------------------------------------------------
# Link to Chrono targets for the requested modules
#--------------------------------------------------------------

target_link_libraries(my_demo PRIVATE ${CHRONO_TARGETS})

#--------------------------------------------------------------
# 4. OPTIONAL (Windows only)
# 
# Add a custom command for copying all Chrono and
# dependency DLLs to the appropriate binary output folder.
# This function has effect only on Windows.
#--------------------------------------------------------------

add_CHRONO_DLLS_copy_command()
