#!/bin/bash -e

# initialize work with regatta
source ../../svt.conf

readonly BLUEGENE=10.6.7.50
readonly HOST=$RG_BG_LOGIN@$BLUEGENE
readonly KEY=$ROOT/keys/$RG_BG_LOGIN
if [ $USER != "root" ];
then
    echo "Please execute this script with root rights"
    exit 1
fi
$TERMINAL -e "cd $ROOT/vpn/; openvpn --config cmc.ovpn" >/dev/null 2>&1 &
echo "Waiting 10 second for vpn connection established..."
sleep 10
scp  -i $KEY $ROOT/script/regatta/* $ROOT/src/*  $HOST:/gpfs/data/vmk-edu-2013/$RG_BG_LOGIN  

$TERMINAL -e "ssh $HOST -i $KEY" >/dev/null 2>&1 &

