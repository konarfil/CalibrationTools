#! /bin/bash

rm -rf CMakeFiles
rm -f CMakeChache.txt

chmod 755 load_enviroment.sh
chmod 755 sub.sh
chmod 755 run.sh

TKEVENT_P='./TKEvent/TKEvent'

cmake -DTKEvent_PATH=$TKEVENT_P ../

make
