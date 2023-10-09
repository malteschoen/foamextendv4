#!/bin/bash
dp=`find ./ -iname eta`
echo $dp
dp=`echo $dp| cut -d/ -f2`
echo $dp
cd ./$dp
sed -i '/#{/,/#}/s@^@//@' U
