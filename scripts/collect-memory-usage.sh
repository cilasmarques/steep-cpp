#!/bin/bash

export LC_NUMERIC="C"
TIME_BETWEEN_COMMANDS=1
echo TIMESTAMP, PID, USED, COMMAND
while [ -e /proc/$1 ]; do
  /usr/bin/pidstat -p $1 -r $TIME_BETWEEN_COMMANDS 1 | grep 'Average:' | awk -v date="$( date +"%s" )" '{print date", "$3", "$8", "$9;}' 2> /dev/null
done

exit