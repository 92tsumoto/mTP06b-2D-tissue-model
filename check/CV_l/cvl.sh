#!/bin/sh

cat ../act_time.out | awk '{if(NR%300 == 50)print}' | cat -n > cvl1
cat ../act_time.out | awk '{if(NR%300 == 100)print}' | cat -n > cvl2
cat ../act_time.out | awk '{if(NR%300 == 150)print}' | cat -n > cvl3
cat ../act_time.out | awk '{if(NR%300 == 200)print}' | cat -n > cvl4
cat ../act_time.out | awk '{if(NR%300 == 250)print}' | cat -n > cvl5
awk -f ./cvl.awk cvl1 | awk '{print "Line  50;",1.0/(0.001*($2-$1)) " (cm/s)"}'
awk -f ./cvl.awk cvl2 | awk '{print "Line 100:",1/(0.001*($2-$1)) " (cm/s)"}'
awk -f ./cvl.awk cvl3 | awk '{print "Line 150:",1/(0.001*($2-$1)) " (cm/s)"}'
awk -f ./cvl.awk cvl4 | awk '{print "Line 200:",1/(0.001*($2-$1)) " (cm/s)"}'
awk -f ./cvl.awk cvl5 | awk '{print "Line 250:",1/(0.001*($2-$1)) " (cm/s)"}'
