#!/bin/bash
# REQUIRED LIBRARIES: libcairo libseqan2-dev
#To avoid permission issue: chmod +x install.sh
#chmod +x install.sh

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"
CURRENT_DIR=$(pwd)

cd $DIR
# Check if tmp exists or destroys it
if [ ! -d tmp ]; then
	mkdir tmp
else
	rm -rf ./tmp
	mkdir tmp
fi

# Check if bin exists or destroys it
if [ ! -d bin ]; then
	mkdir bin
else
	rm -rf ./bin
	mkdir bin
fi

# Makes into temp directory
cd ./tmp
cmake .. || { # Catch
	echo "ERROR: \"cmake\" command failed!"
	echo "WARNING: Compilation failed!"
	echo "Please check that cairo and python3-dev libraries are installed"
	exit 1
}
make || { # Catch
	echo "ERROR: \"make\" command failed!"
	echo "WARNING: Compilation failed!"
	echo "Please check that cairo and python3-dev libraries are installed"
	exit 1
}
# Moves compiled to bin
mv COCOA ../bin || { # Catch
	echo "ERROR: Could not find COCOA binary!"
	echo "WARNING: Compilation failed!"
	exit 1
}
mv papaya ../bin || { # Catch
	echo "ERROR: Could not find papaya binary!"
	echo "WARNING: Compilation failed!"
	exit 1
}
cd ../
rm -rf ./tmp
# SUCCESS
echo ""
echo "All done!"
echo ""
echo "To run from any folder, add the following line to ~/.bashrc (or other shell file):"
echo "	export PATH=$(pwd)/bin:\$PATH"
echo ""
exit 0
