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

find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

set( USR_SOURCES
	usrClientIGT.cxx
	usrVolume.cxx
	usrImage.cxx
	usrImageRegistration.cxx
	usrVolumeRegistration.cxx
	usrVolumeReslice.cxx
	usrTransformMatrix.cxx
	usrVolumeCropping.cxx
	usrImageCropping.cxx
	usrRotationMatrix.cxx
	usrMoveImage.cxx
	)

set( USR_HEADERS
	usrClientIGT.h
	usrVolume.h
	usrImage.h
	usrImageRegistration.h
	usrVolumeRegistration.h
	usrVolumeReslice.h
	usrTransformMatrix.h
	usrVolumeCropping.h
	usrImageCropping.h
	usrRotationMatrix.h
	usrMoveImage.h
	)

add_library(usrLib ${USR_SOURCES} ${USR_HEADERS})
#TARGET_LINK_LIBRARIES(usrLib OpenIGTLink)

#add_executable(TestIGTImageRegistration TestIGTImageRegistration.cxx)
#target_link_libraries(TestIGTImageRegistration usrLib)
#TARGET_LINK_LIBRARIES(TestIGTImageRegistration OpenIGTLink)
#if( "${ITK_VERSION_MAJOR}" LESS 4 )
#  target_link_libraries(TestIGTImageRegistration ITKReview ${ITK_LIBRARIES})
#else( "${ITK_VERSION_MAJOR}" LESS 4 )
#  target_link_libraries(TestIGTImageRegistration ${ITK_LIBRARIES})
#endif( "${ITK_VERSION_MAJOR}" LESS 4 )
#target_link_Libraries(TestIGTImageRegistration ${VTK_LIBRARIES})
#target_link_libraries(TestIGTImageRegistration ${Glue})
#target_link_libraries(TestIGTImageRegistration ${ARMADILLO_LIBRARIES})

#add_executable(TestCroppedImageRegistration TestCroppedImageRegistration.cxx)
#target_link_libraries(TestCroppedImageRegistration usrLib)
#TARGET_LINK_LIBRARIES(TestCroppedImageRegistration OpenIGTLink)
#if( "${ITK_VERSION_MAJOR}" LESS 4 )
#  target_link_libraries(TestCroppedImageRegistration ITKReview ${ITK_LIBRARIES})
#else( "${ITK_VERSION_MAJOR}" LESS 4 )
#  target_link_libraries(TestCroppedImageRegistration ${ITK_LIBRARIES})
#endif( "${ITK_VERSION_MAJOR}" LESS 4 )
#target_link_Libraries(TestCroppedImageRegistration ${VTK_LIBRARIES})
#target_link_libraries(TestCroppedImageRegistration ${Glue})
#target_link_libraries(TestCroppedImageRegistration ${ARMADILLO_LIBRARIES})

add_executable(TestReslicedImageRegistration TestReslicedImageRegistration.cxx)
target_link_libraries(TestReslicedImageRegistration usrLib)
TARGET_LINK_LIBRARIES(TestReslicedImageRegistration OpenIGTLink)
if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(TestReslicedImageRegistration ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(TestReslicedImageRegistration ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )
target_link_Libraries(TestReslicedImageRegistration ${VTK_LIBRARIES})
target_link_libraries(TestReslicedImageRegistration ${Glue})
target_link_libraries(TestReslicedImageRegistration ${ARMADILLO_LIBRARIES})

#add_executable(TestCroppedVolumeRegistration TestCroppedVolumeRegistration.cxx)
#target_link_libraries(TestCroppedVolumeRegistration usrLib)
#TARGET_LINK_LIBRARIES(TestCroppedVolumeRegistration OpenIGTLink)
#if( "${ITK_VERSION_MAJOR}" LESS 4 )
#  target_link_libraries(TestCroppedVolumeRegistration ITKReview ${ITK_LIBRARIES})
#else( "${ITK_VERSION_MAJOR}" LESS 4 )
#  target_link_libraries(TestCroppedVolumeRegistration ${ITK_LIBRARIES})
#endif( "${ITK_VERSION_MAJOR}" LESS 4 )
#target_link_Libraries(TestCroppedVolumeRegistration ${VTK_LIBRARIES})
#target_link_libraries(TestCroppedVolumeRegistration ${Glue})
#target_link_libraries(TestCroppedVolumeRegistration ${ARMADILLO_LIBRARIES})

add_executable(SimulationTest SimulationTest.cxx)
target_link_libraries(SimulationTest usrLib)
TARGET_LINK_LIBRARIES(SimulationTest OpenIGTLink)
if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(SimulationTest ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(SimulationTest ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )
target_link_Libraries(SimulationTest ${VTK_LIBRARIES})
target_link_libraries(SimulationTest ${Glue})
target_link_libraries(SimulationTest ${ARMADILLO_LIBRARIES})

