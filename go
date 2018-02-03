#! /bin/sh
com=test_lift
com=example
in_img='in.pgm'
out_img='out.pgm'

make lib

cd $com
make -W $com.cpp

$com $in_img $out_img
#gdb --arg $com $in_img $out_img
#kuickshow $out_img
