#!/bin/bash -e

# initialize work with regatta
source ../../svt.conf

readonly REGATTA=regatta.cs.msu.ru
readonly HOST=$RG_BG_LOGIN@$REGATTA
readonly KEY=$ROOT/keys/$RG_BG_LOGIN
scp  -i $KEY $ROOT/script/regatta/* $ROOT/src/*  $HOST:/home/$RG_BG_LOGIN  

$TERMINAL -e "ssh $HOST -i $KEY" >/dev/null 2>&1 &

