#!/bin/bash -e

# initialize work with regatta
source ../../svt.conf

readonly BLUEGENE=10.6.7.50
readonly HOST=$RG_BG_LOGIN@$BLUEGENE
readonly KEY=$ROOT/keys/$RG_BG_LOGIN
scp  -i $KEY $ROOT/script/regatta/* $ROOT/src/*  $HOST:/home/$RG_BG_LOGIN  

$TERMINAL -e "openvpn --config $ROOT/keys/cmc.ovpn"
$TERMINAL -e "ssh $HOST -i $KEY" >/dev/null 2>&1 &

