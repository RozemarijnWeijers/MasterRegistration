cmake_minimum_required(VERSION 2.8)
 
project(ImageRegistrationMethodAffine)
 
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

set(OpenIGTLink_USE_FILE_DIR "/home/snr/Support/OpenIGTLink-build")
find_package(OpenIGTLink REQUIRED)
include(${OpenIGTLink_USE_FILE})

add_executable(ImageRegistrationMethodAffine MACOSX_BUNDLE ImageRegistrationMethodAffine.cxx)
TARGET_LINK_LIBRARIES(ImageRegistrationMethodAffine OpenIGTLink)
if( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(ImageRegistrationMethodAffine ITKReview ${ITK_LIBRARIES})
else( "${ITK_VERSION_MAJOR}" LESS 4 )
  target_link_libraries(ImageRegistrationMethodAffine ${ITK_LIBRARIES})
endif( "${ITK_VERSION_MAJOR}" LESS 4 )
target_link_Libraries(ImageRegistrationMethodAffine ${VTK_LIBRARIES})
target_link_libraries(ImageRegistrationMethodAffine ${Glue})

