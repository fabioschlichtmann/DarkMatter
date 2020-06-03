#!/bin/bash

for item in `cat $1`
do 
 echo FILE name = $item
 cp ./job.submit_template ./sub_job${item}.submit
 echo  "arguments       =" $item >> sub_job${item}.submit
 echo  "output          =/home/ceres/schmah/ALICE/DarkMatter/submit/log/out"$item".log" >> sub_job${item}.submit
 echo  "error           =/home/ceres/schmah/ALICE/DarkMatter/submit/log/err"$item".log" >> sub_job${item}.submit
 echo  "log             =/home/ceres/schmah/ALICE/DarkMatter/submit/log/log"$item".log" >> sub_job${item}.submit
 echo  'requirements    =((UtsnameNodename =!= "alice-serv13"))' >> sub_job${item}.submit
 echo  'requirements    =!stringListMember(UtsnameNodename, "alice-serv13")' >> sub_job${item}.submit
 echo  "queue" >> sub_job${item}.submit
 echo "submit: condor_submit" sub_job${item}.submit
 condor_submit sub_job${item}.submit
done
