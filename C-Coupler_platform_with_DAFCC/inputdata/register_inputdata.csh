#!/bin/csh

# == Get the path of this script ==
set called=($_)

if ("$called" != "") then
        set me="$called[2]"
else
        set me="$0"
endif
set MYPATH=`readlink -f "$me"`
set MYPATH=`dirname "$MYPATH"`
# =================================

setenv DATAROOT $MYPATH
