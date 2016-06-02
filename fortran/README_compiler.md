In order to use "Makeifort" Makefile first rename it to Makefile
then execute the following sed command:

sed -i 's/\!USE IFPORT/USE IFPORT/g' *.f90

This adds a module with contains rand() and getpid()

If you want to return to the old version, just restor ethe old Makefile and execute

sed -i 's/USE IFPORT/\!USE IFPORT/g' *.f90
