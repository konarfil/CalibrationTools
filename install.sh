rm -rf build lib
mkdir build lib

cd build
cmake ..
make
cd ..
cp -f build/libCalibrationModule.so lib/