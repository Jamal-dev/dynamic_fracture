#
#  CMake script for the step-3 tutorial program:
##

# Set the name of the project and target:
SET(TARGET "dynamic_fracture")
#SET(TARGET "test")


SET(TARGET_SRC 
${TARGET}.cc 
${TARGET_SRC}  
${TARGET_INC}
)



# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 9.2 QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

PROJECT(${TARGET})


DEAL_II_INVOKE_AUTOPILOT()



#set_property(TARGET cell PROPERTY COMPILE_DEFINITIONS deal_II_dimension=2
add_definitions(-Ddeal_II_dimension=2)
