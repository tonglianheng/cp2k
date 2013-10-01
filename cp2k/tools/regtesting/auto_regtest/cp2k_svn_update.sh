#!/bin/bash
# cp2k_svn_update
# CP2K repository update script
# For those remote systems that do not support SVN in compute nodes
# Script accepts three arguments

if [[ $# -ne 3 ]]
then
  echo "Command line: cp2k_svn_update conf_file cp2kdir wwwtestdir"
  exit 500
fi

# First argument is configuration file  
conf_file=$1

# Second and third arguments are cp2k and wwwtest directories
# cp2kdir is equivalent to dir_base in do_regtest script
cp2kdir=$2
wwwtestdir=$3

# Extracting necessary info to form LAST directory
dir_triplet=`cat ${conf_file} | grep "dir_triplet" | sed 's/.*=\(.*\)/\1/'`
cp2k_version=`cat ${conf_file} | grep "cp2k_version" | sed 's/.*=\(.*\)/\1/'`
emptycheck=`cat ${conf_file} | grep "emptycheck" | sed 's/.*="\(.*\)"/\1/'`
dir_last=${cp2kdir}/LAST-${dir_triplet}-${cp2k_version}

# Products of the script
changelog_diff=${wwwtestdir}/ChangeLog.diff
changelog_diff_tests=${wwwtestdir}/ChangeLog-tests.diff
error_description_file=${wwwtestdir}/svn_error_summary

# *******************************************************************
# From this point, SVN repository in cp2kdir is checked for updates  
# *******************************************************************

# *** svn update src
cd ${cp2kdir}/cp2k/src

# If ChangeLog does not exist for the first time, it is created 
if [[ ! -s ${dir_last}/ChangeLog ]]  
then
${cp2kdir}/cp2k/tools/svn2cl/svn2cl.sh -i -o ${dir_last}/ChangeLog 
fi

svn update &> out
if (( $? )); then
echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" >> ${error_description_file}
tail -20 out >> ${error_description_file}
echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" >> ${error_description_file}
echo "error happened : no svn update ... bailing out"
cat ${error_description_file}
exit 1
fi

echo "--- svn update src ---"
cat out
echo "svn src update went fine"
cp2k_lines=`wc *.F | tail -1 |  awk  '{print $1}'`
echo "cp2k is now ${cp2k_lines} lines .F"

# *** using svn2cl.pl to generate GNU like changlelog
${cp2kdir}/cp2k/tools/svn2cl/svn2cl.sh --limit 100 -i -o ChangeLog.new &> out

line1="$(head -n 1 ${dir_last}/ChangeLog)"
nline=$(grep -n "${line1}" ChangeLog.new | head -n 1 | cut -f 1 -d:)
head -n $((nline - 1)) ChangeLog.new >ChangeLog
cat ${dir_last}/ChangeLog >>ChangeLog

diff ChangeLog ${dir_last}/ChangeLog > ${changelog_diff}
echo "------- differences --------" >> ${changelog_diff}

rm ChangeLog.new 

echo "---  changelog diff src  ---"
cat ${changelog_diff} 
echo "----------------------------"

cp ChangeLog ${dir_last}/ChangeLog

rm out

# *** svn update test 
cd ${cp2kdir}/cp2k/tests

if [[ ! -s ${dir_last}/ChangeLog-tests ]]  
then
${cp2kdir}/cp2k/tools/svn2cl/svn2cl.sh -i -o ${dir_last}/ChangeLog-tests
fi

svn update &> out
if (( $? )); then
echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" >> ${error_description_file}
tail -20 out >> ${error_description_file}
echo "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx" >> ${error_description_file}
echo "error happened : no svn update ... bailing out"
cat ${error_description_file}
exit 1
fi

echo "--- svn update test ---"
cat out
echo "svn tests update went fine"

# *** using svn2cl.pl to generate GNU like changlelog
${cp2kdir}/cp2k/tools/svn2cl/svn2cl.sh --limit 100 -i -o ChangeLog-tests.new &> out

line1="$(head -n 1 ${dir_last}/ChangeLog-tests)"
nline=$(grep -n "${line1}" ChangeLog-tests.new | head -n 1 | cut -f 1 -d:)
head -n $((nline - 1)) ChangeLog-tests.new >ChangeLog-tests
cat ${dir_last}/ChangeLog-tests >>ChangeLog-tests

diff ChangeLog-tests ${dir_last}/ChangeLog-tests > ${changelog_diff_tests}
echo "------- differences --------" >> ${changelog_diff_tests}

rm ChangeLog-tests.new 

echo "---  changelog diff tests  ---"
cat ${changelog_diff_tests} 
echo "----------------------------"

cp ChangeLog-tests ${dir_last}/ChangeLog-tests

rm out

# In the end, script checks if there are any changes
# -emptycheck equivalent from do_regtest script
if [[ ${emptycheck} == "YES" ]] 
then
  isempty_1=`nl ${changelog_diff} | awk '{print $1}'`
  isempty_2=`nl ${changelog_diff_tests} | awk '{print $1}'`

  if [[ ${isempty_1} == "1" && ${isempty_2} == "1" ]]; then
    echo "No changes since last run -- clean exit without testing"
    exit 2
  else
    echo "Code has changed since last run -- continue regtest"
    exit 0
  fi
fi

exit 0
