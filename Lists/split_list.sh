#!/bin/ksh
# split_list.sh file_list.txt  --> loops over file_list.txt textfile which contains
# a list of files. file_list.txt is then splittet with split level nsplit
if [ $# -eq 1 ] 
then
        counter=0
        nsplit=5
        nlow=1
        nup=$nsplit
        filename=Split_DM_V6_
        fileb=-
        fileend=.txt
        fullfilename="$filename$nlow$fileb$nup$fileend"
        echo fullfilename = $fullfilename
        rm ${fullfilename}
        for item in `cat $1`
	do 
          echo FILE name = $item
          echo ${item}>>${fullfilename}
          let "counter=counter+1"
          if [ $counter -ge $nup ]
          then
          let "nlow=nlow+nsplit"
          let "nup=nup+nsplit"
          fullfilename="$filename$nlow$fileb$nup$fileend"
          echo fullfilename = $fullfilename
          rm ${fullfilename}
          echo N Counter: ${n}
          fi
          
        #echo "**************************************************************************************"
        #echo   FILE COUNTER: ${counter}
        #echo "**************************************************************************************"
     
        done
else
	echo "**************************************************************************************"
 	echo "\t" usage: mycompare text-file1 text-file2 
	echo 
	echo "\t\t" returns lines of text-file1 which are not contained in text-file2
	echo "**************************************************************************************" 
fi