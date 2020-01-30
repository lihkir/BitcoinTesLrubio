#!/bin/bash

echo "Copying .cpp files"
find . -name '*.cpp' | cpio -pdm /mnt/d/MEGA/SharedUbuntu/pvm_sed_ptrs_git
hecho "Done"

echo "Copying .h files"
find . -name '*.h' | cpio -pdm /mnt/d/MEGA/SharedUbuntu/pvm_sed_ptrs_git
echo "Done"

echo "Copying Makefile"
find . -name 'Makefile' | cpio -pdm /mnt/d/MEGA/SharedUbuntu/pvm_sed_ptrs_git
echo "Done"

echo "Copying .sh files"
find . -name '*sh' | cpio -pdm /mnt/d/MEGA/SharedUbuntu/pvm_sed_ptrs_git
echo "Done"