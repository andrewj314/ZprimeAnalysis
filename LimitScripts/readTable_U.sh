counter=1
limitCard=zprime_mutau
while read col1 col2 col3 col4 col5 col6 col7 col8 col9 col10 col11 col12
do
    counter=`expr $counter + 1`

    for sig in 2 3 4 5 6 7
    do

     eval signalcol=\$col$sig
       
#     cp $limitCard".txt" $limitCard"_signal_"$sig"_"$counter".txt"
     sed -i .bak "s:US:$signalcol:g"   $limitCard"_signal_"$sig"_"$counter".txt"
     sed -i .bak "s:UQ:$col8:g"    $limitCard"_signal_"$sig"_"$counter".txt"
     sed -i .bak "s:UDY:$col9:g"   $limitCard"_signal_"$sig"_"$counter".txt"
     sed -i .bak "s:UTT:$col10:g"    $limitCard"_signal_"$sig"_"$counter".txt"
     sed -i .bak "s:UW:$col11:g"    $limitCard"_signal_"$sig"_"$counter".txt"
     sed -i .bak "s:UDI:$col12:g" $limitCard"_signal_"$sig"_"$counter".txt"
    done
done < output_zprime_U.txt
rm *.bak
