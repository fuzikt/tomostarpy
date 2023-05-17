#!/bin/bash

if [[ ! -f lib/requirements.txt ]]; then
echo "Please run this script from inside the tomostarpy directory!"
exit
fi

if ! grep -q tomostarpy lib/requirements.txt; then
echo "Please run this script from inside the tomostarpy directory!"
exit
fi

if [[ -d venv ]]; then
  echo ">>Old venv drectory found. Removing it..."
  rm -r venv
fi

echo ">>Creating virtual environment for tomostarpy..."
if ! python3 -m venv ./venv; then
  echo "!!! There was an error during the creation of the virtual environment."
  echo "Check that you have python venv installed!"
  echo "On Ubuntu systems use: apt install python3.8-venv"
  echo "Installation failed !!!"
  exit
fi

source venv/bin/activate

pip3 install --upgrade pip

echo ">>Installing required libs into the new environment..."
pip3 install -r lib/requirements.txt

echo ">>Compiling the required libraries..."
./build_libs.bsh

echo "export CUDA_PATH=$(pwd)/venv/lib64/python3.8/site-packages/cuda" >> venv/bin/activate
echo "export LD_LIBRARY_PATH=/soft/tomostarpy/venv/lib64/python3.7/site-packages/nvidia/cuda_nvrtc/lib:$LD_LIBRARY_PATH" >> venv/bin/activate
echo "---------------------------------------------------------------"
echo "!!! Before using the tomostarpy scripts don't forget to source:"
echo "source <tomostarpydir>/venv/bin/activate"
echo "All done! Have fun!"