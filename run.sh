#!/bin/sh
# named the folder as you want its better in this way: folder="case_year_month_day-execution". Modify lines that say FOLDER, FILES, RUNNING.

# main foulder name
foldercase="shakenumM"

# subfoulder name
folder="ShakeM_2024_8_15_N"


echo "*"
echo "*************CLEANING PREVIOUS COMPILATION*************"
echo "*"
#cmake .
#make clean
echo "*"
echo "******************COMPILING THE CODE******************"
echo "*"


cd ../$foldercase
               #  0.6 1.0 2.5 6.3 15.8 25.0 39.0 63.0 REO
               #  1.5 2.4 6.0 15.0 38.0 60.0 90.0
               #   1 2 3 4 5  6 7 8 9 10   12 15 18 21 23   28 33 38 45 55 65 75
#---------------SIMULATION CASE------------------------
for vari in    28 33 38 45 55 65 75   
do
          

     echo " "
     echo "------------------Checking folders------------------"
     echo " "
     
     if ! [[ -d $folder$vari ]]; then # if there is NOT a folder with this name, make the folder
          mkdir $folder$vari
     fi
     
     cd /home/gezhuan/projects/shaking/${foldercase}/${folder}${vari}
     
     if ls -1qA | grep -q .; then # if folder is NOT EMPTY remove everything and copy the new files
          mv /home/gezhuan/projects/${foldercase}/${folder}${vari} ~/bak
     fi
     
     echo " "
     echo "--------------------Copying files-------------------"
     echo " "
     cp /home/gezhuan/projects/shaking/test/shakecpu .
     cp /home/gezhuan/projects/shaking/test/shake .
     cp /home/gezhuan/projects/shaking/test/input-shake.inp .
     cp /home/gezhuan/projects/shaking/test/post.pvsm .
     cp /home/gezhuan/projects/shaking/test/postcu.pvsm .

     sed -i "s/vari/$vari/g" input-shake.inp

     echo "*"
     echo "******************EXECUTING THE CODE******************"
     echo "*"
     ./shakecpu input-shake 5
     cd ..
     
     echo "Saving in folder: "$folder$vari
done

echo " "
squeue -u gezhuan
