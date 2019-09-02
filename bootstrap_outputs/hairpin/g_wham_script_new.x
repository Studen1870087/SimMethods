#!/bin/bash
#
let k=1	
let nst=21		# Nr of struct		
let nBootstrap=20	# Nr of bootstap iterations
#
g_wham=g_wham_s
#
rm tpr-files.dat
rm pullf-files.dat
#
echo ' '
#
while [ $k -le $nst ] 
do
    let k1=k
    #
    echo ' Processing structure n...' $k1
    echo $k1-struct/umbrella.tpr >> tpr-files.dat
    echo $k1-struct/pullf-umbrella.xvg >> pullf-files.dat
    #
    let k=k+1
done
#
echo ' '
echo ' g_wham is running...'
echo ' '
#
$g_wham -it tpr-files.dat -if pullf-files.dat -o pmf.xvg -hist -ac -bs-method traj-gauss -unit kCal -v -nBootstrap $nBootstrap >& g_wham.log
#
echo ' '
echo ' Done'
echo ' '
