#!/bin/csh 
setenv ROOT `pwd`

set input = "list" 

foreach line ( "`cat $input`" )
    cd $line
    `sed -i -e 's/mem=100/mem=200/g' optstep.sh`
    sbatch optstep.sh
    cd $ROOT
end
