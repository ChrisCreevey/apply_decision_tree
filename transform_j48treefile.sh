grep "^N[0-9]" test_tree.txt | sed 's/\[label="//g' | sed 's/" *\]//g'| sed 's/\".*\]$//g' | sed "s/\(N[0-9]*\) /\1   /g" | sed '/->/s/^/EDGE /g' | sed '/^N[0-9]/s/^/NODE    /g' | awk '{print $1"\t"$2"\t"$3"\t"$4"#"}' > formatted_tree_file.txt


