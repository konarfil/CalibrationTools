rm -rf build lib
mkdir build lib

cd build
cmake ..
make
cd ..
cp -f build/libCharge2EnergyModule.so lib/