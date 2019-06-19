#! /bin/csh

# NR is a predefined awk variable that indicates the number of line being read
# the () indicates conditional propositions, {} indicates the action to do in the case it's true.
# $0 is the whole line, $1 is the first column, and so on. 


set file=$argv

if ($#argv != 1)  then
        echo
        echo Usage:    ./program your_file
        echo
	exit
endif

# important to set it in english, otherwise the decimal separator screws the calculations
setenv LC_NUMERIC en

cat $file |  awk '(NR==1){x0=$1;  sum=0;}; (NR != 1) {print x0+($1-x0)/2,sum} { sum+= $2*($1-x0); x0=$1; }'

