# Execute the sequence of commands to build the library
echo "compiling stochastic TACS [STACS]"
cd src/cpp
make clean
rm lib/*.so
rm lib/*.a
make
cd -
mkdir -p lib
mv src/cpp/*.so lib/.
mv src/cpp/*.a lib/.
