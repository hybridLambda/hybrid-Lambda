#!/bin/bash

#if [ ! -z "$src_dir+x" ]; then
src_dir=$PWD
#fi

if [ ! -f "${src_dir}/main.cpp" ]; then
  "Error: please execute \".checkfiles.sh\" from hybrid-Lambda src/ directory"
  exit 1
fi


OS=`uname -s`
#echo "operating system is ${OS}"
VERSION=$(cat ../version | tr -d '\n')
echo -n "." # echo "ok not git clone, now check plot and figure" 

#rm master.tar.gz  # just making sure download the one we actually want to 

extract_tmp(){
	if [ ! -d ${tmp}_dir ]; then mkdir ${tmp}_dir; fi
	echo -n "." # echo "ok, do something" # 

	if [[ "${OS}" == "Linux" ]]; then
		wget --no-check-certificate https://github.com/hybridLambda/${tmp}/archive/${VERSION}.tar.gz -o /dev/null
	elif [[ "${OS}" == "Darwin"* ]]; then
		curl -LOk https://github.com/hybridLambda/${tmp}/archive/${VERSION}.tar.gz
	else
		echo "Unknown OS, fail to download package from https://github.com/hybridLambda/${tmp}/archive/${VERSION}.tar.gz" 
		echo "Please contact Joe at sha.joe.zhu@gmail.com if assistance is needed"
		exit 1
	fi

	if [ ! -f ${VERSION}.tar.gz ]; then
		echo "Error: Download package from https://github.com/hybridLambda/${tmp}/archive/${VERSION}.tar.gz failed"
		echo "Please contact Joe at sha.joe.zhu@gmail.com if assistance is needed"
		exit 1
	fi
	tar -xf ${VERSION}.tar.gz -C ${tmp}_dir
	if [ ! -d ${tmp} ]; then mkdir ${tmp}; fi
	cp ${tmp}_dir/*/* ${tmp}/
	rm -r ${tmp}_dir ${VERSION}.tar.gz	
}

tmp=plot
tmp_dir=${src_dir}/${tmp}
if [ -f ${tmp_dir}/figure.cpp  -a -f ${tmp_dir}/figure.hpp  -a -f ${tmp_dir}/test_figure.cpp ];then
	echo -n "." # echo "ok, do nothing" 
else
	extract_tmp
fi

tmp=freq
tmp_dir=${src_dir}/${tmp}
if [ -f ${tmp_dir}/freq.cpp  -a -f ${tmp_dir}/freq.hpp  -a -f ${tmp_dir}/freq_extra.cpp ];then
	echo -n "." # echo "ok, do nothing" 
else
	extract_tmp
fi

