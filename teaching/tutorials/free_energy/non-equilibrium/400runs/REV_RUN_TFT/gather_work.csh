#! /bin/csh

set count=0


while ( $count < 400 )


cd run_$count

../turbo_integration.csh dhdl.xvg | tail -n1 | awk '{print $2/5}'

@ count ++

cd ..

end
