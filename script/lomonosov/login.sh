#!/bin/bash -e

# initialize work with regatta
source ../../svt.conf

readonly LOMONOSOV=lomonosov.parallel.ru
readonly HOST=$LOM_LOGIN@$LOMONOSOV
readonly KEY=$ROOT/keys/$RG_BG_LOGIN
scp  -i $KEY $ROOT/script/lomonosov/* $ROOT/src/*  $HOST:/home/$LOM_LOGIN/_scratch  

$TERMINAL -e "ssh $HOST -i $KEY" >/dev/null 2>&1 &

