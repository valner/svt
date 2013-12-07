#!/bin/bash -e 
# configure working with supercomputers
readonly CONFIG_NAME=svt.conf
CONFIG=$PWD/$CONFIG_NAME

echo "Please answer for some questions for configure connection"
echo -n "Root dir [default: current dir]? "
read root
if [ $root ];
then
    if [ -d $root ];
    then
        CONFIG=$root/svt.conf
        echo "ROOT=$root" > $CONFIG
    else
        echo "Bad path. Setting to default value"
    fi
else
    echo "ROOT=$PWD" > $CONFIG
fi

echo ""

echo -n "Login for regatta and BG [default: edu-vmk-stud-525-003]? "
read rg_bg_login
if [ $rg_bg_login ];
then
    echo "RG_BG_LOGIN=$rg_bg_login" >> $CONFIG
else
    echo "RG_BG_LOGIN=edu-vmk-stud-525-003" >> $CONFIG
fi

echo ""

echo -n "Login for lomonosov [default: edu-vmk13-studvalner47_251659]? "
read lom_login
if [ $lom_login ];
then
    echo "LOM_LOGIN=$lom_login" >> $CONFIG
else
    echo "LOM_LOGIN=edu-vmk13-studvalner47_251659" >> $CONFIG
fi

echo ""

echo -n "GUI Terminal binary [default: xterm]? "
read terminal
if [ $terminal ];
then
    if command -v $terminal > /dev/null 2>&1; 
    then
        echo "TERMINAL=$terminal" >> $CONFIG
    else
        echo "Bad terminal command: $terminal. Setting default"
        echo "TERMINAL=xterm" >> $CONFIG
    fi
else
    echo "TERMINAL=xterm" >> $CONFIG
fi


