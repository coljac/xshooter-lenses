#!/bin/bash

input=regions
output=regions/fixed
base=.
band=g
scale="mode 99"

if [ "$2" != "" ]; then
    band=$2
fi
for f in `cat $1|cut -d, -f2|grep -v name`; do
    opends9=1
    if [ -e $output/$f.reg ] 
    then
        continue
    fi
    echo 
    echo $f
    echo "#! /bin/bash\nxpaset -p ds9 regions system wcs
\nxpaset -p ds9 regions save regions/$f.reg" > savereg.sh
    targetdone=0
    while [ $targetdone == 0 ]
    do
        # ./ds9.sh $f $band "mode 99.5" &
        ra=`grep $f $1|cut -d, -f3`
        dec=`grep $f $1|cut -d, -f4`
        if [ "$opends9" == "1" ]
        then
            ds9 fits/1024/${f}_${band}.fits -regions regions/${f}.reg -zoom 4 -scale `echo $scale` -geometry 1100x1100+300+100 -frame center -pan to $ra $dec wcs -mode region &
            pid=$!
            sleep 1
        fi
        opends9=1

        # xpaset -p ds9 scale zscal
        # eog -f fits/100/$f.png &
        # read -k1 -s "task?q to redo, p to change pa and redo, g or r to change band "
        read -n 1 -r -s -p "q to redo, p to change pa and redo, g or r to change band " task
        echo
        if [ "$task" == "q" ]
        then
            xpaset -p ds9 quit
            # kill $pid
            # Don't move on
        elif [ "$task" == "p" ]
        then
            newpa=""
            read -p "New position angle? " newpa
            bash changepa.sh $f $newpa
            xpaset -p ds9 quit
            # Then go again
        elif [ "$task" == "r" ]
        then
            band=r
        elif [ "$task" == "g" ]
        then
            band=g
        elif [ "$task" == "i" ]
        then
            feh fits/100/$f.png &
            opends9=0
        elif [ "$task" == "s" ]
        then
            targetdone=1
            xpaset -p ds9 regions system wcs
            xpaset -p ds9 regions save $output/$f.reg
            xpaset -p ds9 quit
            # killall eog
        else
            opends9=0
        fi
    done
done

