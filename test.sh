#!/usr/bin/sh

executable="${1:-./main}"
status=0

if $executable -n 7 -i 300  | cmp reference/7.txt; then
    echo "n= 7 i=300 correct"
else
    echo "n= 7 i=300 failed!"
    status=1
fi

if $executable -n 15 -i 300  | cmp reference/15.txt; then
    echo "n=15 i=300 correct"
else
    echo "n=15 i=300 failed!"
    status=1
fi

if $executable -n 51 -i 300  | cmp reference/51.txt; then
    echo "n=51 i=300 correct"
else
    echo "n=51 i=300 failed!"
    status=1
fi

exit $status
