#!/bin/bash

# --------------------------------------------------------------------
# Script for recrusively updating TAGS file used by emacs for a source
# tree.
# --------------------------------------------------------------------
# Lianheng Tong (2014), tonglianheng@gmail.com
# --------------------------------------------------------------------

######################################################################
# You should modify the following to suit your source code
######################################################################

ETAGS=etags
SRC_ROOT=.
INCLUDE_REGEXP='.*\.(F|f90|h)'
IGNORE_REGEXP='[^a-zA-Z]*svn[^a-zA-Z]*'

# If $ETAGS_REGEXP_FILE is found, then -r option of etags is used.
ETAGS_REGEXP_FILE=./etags_regexp

######################################################################
# You do not need to modify the code below
######################################################################

if [ "x$ETAGS_REGEXP_FILE" != "x" ] && \
   [ -f "$ETAGS_REGEXP_FILE" ] ; then
    echo "$ETAGS_REGEXP_FILE file found, using -r option of etags"
    etags_command="$ETAGS -r @etags_regexp"
else
    etags_command="$ETAGS"
fi

find -E . -regex $INCLUDE_REGEXP -and \
     -not -regex $IGNORE_REGEXP \
     -print | xargs $etags_command
