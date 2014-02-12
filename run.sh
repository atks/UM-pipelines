cat $1 | grep -v "#"| tr -d '\n\\' | sh
