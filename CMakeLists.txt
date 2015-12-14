cmake_minimum_required(VERSION 2.8)
 
project(ImageRegistration)
 
find_package(ITK REQUIRED)
include(${ITK_USE_FILE})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

if (ITKVtkGlue_LOADED)
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
else()
  find_package(ItkVtkGlue REQUIRED)
  include(${ItkVtkGlue_USE_FILE})
  set(Glue ItkVtkGlue)
endif()

find_package(OpenIGTLink REQUIRED)
include(${OpenIGTLink_USE_FILE})

set( USR_SOURCES
	usrClientIGT.cxx
	usrVolume.cxx
	usrImage.cxx
	usrImageRegistration.cxx
	usrVolumeReslice.cxx
	usrTransformMatrix.cxx
	)

set( USR_HEADERS
	usrClientIGT.h
	usrVolume.h
	usrImage.h
	usrImageRegistration.h
	usrVolumeReslice.h
	usrTransformMatrix.h
	)

add_library(usrLib ${USR_SOURCES} ${USR_HEADERS})
#TARGET_LINK_LIBRARIES(usrLib OpenIGTLink)

add_executable(ImageRegistrationTest MainTestRegistration.cxx)
target_link_libraries(ImageRegistrationTest usrLib)
TARGET_LINK_LIBRARIES(ImageRegistrationTest OpenIGTLink)
if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(ImageRegistrationTest ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(ImageRegistrationTest ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )
target_link_Libraries(ImageRegistrationTest ${VTK_LIBRARIES})
target_link_libraries(ImageRegistrationTest ${Glue})

add_executable(TestRegistration testRegistration.cxx)
target_link_libraries(TestRegistration usrLib)
TARGET_LINK_LIBRARIES(TestRegistration OpenIGTLink)
if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(TestRegistration ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(TestRegistration ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )
target_link_Libraries(TestRegistration ${VTK_LIBRARIES})
target_link_libraries(TestRegistration ${Glue})
