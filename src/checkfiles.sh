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
	mkdir ${tmp}_dir
	echo -n "." # echo "ok, do something" # 

	if [[ "${OS}" == "Linux" ]]; then
		wget --no-check-certificate https://github.com/hybridLambda/${tmp}/archive/${VERSION}.tar.gz -o /dev/null
	elif [[ "${OS}" == "Darwin"* ]]; then
		curl -LOk https://github.com/hybridLambda/${tmp}/archive/master.zip
	else
		echo "Unknown OS, fail to download package from https://github.com/hybridLambda/${tmp}/archive/master.tar.gz" 
		exit 1
	fi

	if [ ! -f master.tar.gz ]; then
		echo "Error: Download package from https://github.com/hybridLambda/${tmp}/archive/master.tar.gz failed"
		exit 1
	fi
	tar -xf master.tar.gz -C ${tmp}_dir
	cp ${tmp}_dir/*/* ${tmp}/
	rm -r ${tmp}_dir ${VERSION}.tar.gz	
}

tmp=plot
tmp_dir=${src_dir}/${tmp}
if [ -f ${tmp_dir}/figure.cpp  -a -f ${tmp_dir}/figure.hpp  -a -f ${tmp_dir}/test_figure.cpp ];then
	echo -n "." # echo "ok, do nothing" 
else
	extract_tmp()
fi

tmp=freq
tmp_dir=${src_dir}/${tmp}
if [ -f ${tmp_dir}/freq.cpp  -a -f ${tmp_dir}/freq.hpp  -a -f ${tmp_dir}/freq_extra.cpp ];then
	echo -n "." # echo "ok, do nothing" 
else
	extract_tmp()
	#mkdir freq_tmp_dir
	#echo -n "." # echo "ok, do something" # 

	#if [[ "${OS}" == "Linux" ]]; then
		#wget --no-check-certificate https://github.com/hybridLambda/freq/archive/master.tar.gz -o /dev/null
	#elif [[ "${OS}" == "Darwin"* ]]; then
		#curl -LOk https://github.com/hybridLambda/freq/archive/master.zip
	#else
		#echo "Unknown OS, fail to download package from https://github.com/hybridLambda/freq/archive/master.tar.gz" 
		#exit 1
	#fi
	
	#if [ ! -f master.tar.gz ]; then
		#echo "Error: Download package from https://github.com/hybridLambda/freq/archive/master.tar.gz failed"
		#exit 1
	#fi
	#tar -xf master.tar.gz -C freq_tmp_dir
	#cp freq_tmp_dir/*/* freq/
	#rm -r freq_tmp_dir master.tar.gz
fi



#if [[ "$OSTYPE" == "linux-gnu" ]]; then
        ## ...
#elif [[ "$OSTYPE" == "darwin"* ]]; then
        ## Mac OSX
#elif [[ "$OSTYPE" == "cygwin" ]]; then
        ## POSIX compatibility layer and Linux environment emulation for Windows
#elif [[ "$OSTYPE" == "msys" ]]; then
        ## Lightweight shell and GNU utilities compiled for Windows (part of MinGW)
#elif [[ "$OSTYPE" == "win32" ]]; then
        ## I'm not sure this can happen.
#elif [[ "$OSTYPE" == "freebsd"* ]]; then
        ## ...
#else
        ## Unknown.
#fi
