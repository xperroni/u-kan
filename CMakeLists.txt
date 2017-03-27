cmake_minimum_required (VERSION 2.8)

project(ExtendedKF)

add_definitions(-std=c++0x)

add_executable(ukan
    "src/ukf.cpp"
    "src/main.cpp"
    "src/model.cpp"
    "src/model_ctrv.cpp"
    "src/sensors.cpp"
    "src/state.cpp"
)
