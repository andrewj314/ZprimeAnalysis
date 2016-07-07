inputDir=$1
for file in `cat LIST.txt`; do
  echo $file
  cp dumpBinContentU_COPY.C dumpBinContent.C 
  sed -i .bak "s:DIRECTORY:$inputDir:g" dumpBinContent.C
  sed -i .bak "s:FILE:$file:g" dumpBinContent.C 
  root -l -b -q dumpBinContent.C
done