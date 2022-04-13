cmake_minimum_required(VERSION 2.8)

project(Filtration)			# Название проекта

set(SOURCE_EXE Main.cpp)		# Установка переменной со списком исходников для исполняемого файла

set(SOURCE_LIB ArrayFunctions.cpp ChemicalKinetics.cpp EnergyTransfer.cpp Functions.cpp MassTransfer.cpp PhaseSolver.cpp Visualizer.cpp Well.cpp)			# Тоже самое, но для библиотеки

add_library(FunctionsLib STATIC ${SOURCE_LIB})	# Создание статической библиотеки с именем foo

add_executable(exe ${SOURCE_EXE})	# Создает исполняемый файл с именем main

add_custom_target(cleanOutput rm -f -r ../Output/*.vtk)

target_link_libraries(exe FunctionsLib)	
