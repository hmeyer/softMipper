cmake_minimum_required(VERSION 2.6)
PROJECT(softMipper)

FIND_PATH(MYLIBS itkbasics.h ../myLibs ../../myLibs)

SET(SRC ${MYLIBS}/itkbasics.cpp ${MYLIBS}/itkbasics.h ${MYLIBS}/toString.h  ${MYLIBS}/itkVariableProjectImageFilter.txx ${MYLIBS}/itkVariableProjectImageFilter.h ${MYLIBS}/itkFilterFunctions.h)

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)

SET(Libraries
 ITKCommon
 ITKAlgorithms
 ITKIO
)

SET(CurrentExe "softMipCmd")
ADD_EXECUTABLE(${CurrentExe} softmipcmd.cpp ${SRC})
INCLUDE_DIRECTORIES(${MYLIBS})
TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})


