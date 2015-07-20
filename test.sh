#!/bin/sh

#set case file
case=examples/cavity/ 
grid=cavity

#prepare directory
rm -rf run$1
cp -r $case run$1
cp run.sh run$1

#change controls
cd run$1
sed -i '' -e "s/3 {2 1 1}/3 {"$2" "$3" 1}/g" controls
./run.sh $grid $1
cd ..
