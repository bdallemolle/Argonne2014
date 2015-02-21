cmake_minimum_required(VERSION 2.8)

find_package( OpenCV REQUIRED )

SET(OPENCV_CPP_SAMPLES_REQUIRED_DEPS opencv_core opencv_flann opencv_imgproc
    opencv_highgui opencv_ml opencv_video opencv_objdetect opencv_photo opencv_nonfree opencv_softcascade
    opencv_features2d opencv_calib3d opencv_legacy opencv_contrib opencv_stitching opencv_videostab opencv_shape)

ocv_check_dependencies(${OPENCV_CPP_SAMPLES_REQUIRED_DEPS})

include_directories( ${OPENCV_INCLUDE_DIR} ${PROJECT_BINARY_DIR} )

#Executables
add_executable( find_rotrect find_rotrect.cpp )
add_executable( find_quad find_quad.cpp )
add_executable( polyhist polyhist.cpp )
add_executable( autorealign autorealign.cpp )
add_executable( results sum_results.cpp )

#Link libs
target_link_libraries( find_rotrect ${OpenCV_LIBS} )
target_link_libraries( find_quad ${OpenCV_LIBS} )
target_link_libraries( polyhist ${OpenCV_LIBS} )
target_link_libraries( autorealign ${OpenCV_LIBS} )
target_link_libraries( results ${OpenCV_LIBS} )

