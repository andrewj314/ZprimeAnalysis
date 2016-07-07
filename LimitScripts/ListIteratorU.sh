for file in `cat LIST.txt`; do
  cp runCodesU_COPY.sh runCodes.sh
  sed -i .bak "s:LIST:$file:g" runCodes.sh 
  sh runCodes.sh $file
done