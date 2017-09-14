# list node names, in range from $1 to $2
# use example: ./station_name_list.sh 101 116
for a in $(seq $1 $2) ; do ssh stud$a -q -o "ConnectTimeout 1" echo stud$a ; done
