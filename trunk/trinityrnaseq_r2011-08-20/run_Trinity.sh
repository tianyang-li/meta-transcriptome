#!/bin/bash -l

reuse Java-1.6
reuse GCC-4.3
reuse Perl-5.8

ulimit -s unlimited


ALLPATHSLG_BASEDIR=/seq/annotation/bio_tools/ALLPATHSLG/current/
export ALLPATHSLG_BASEDIR

`dirname $0`/Trinity.pl $*

