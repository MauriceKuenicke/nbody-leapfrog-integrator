#!/bin/bash
#
# make_movie.sh: script to generate an avi from all PNG images in current
#                directory; to use type the following commands
#
#                $ ssh celeborn
#                $ cd dir/where/PNGs/are
#                $ ./make_movie.sh
#                $ exit
#                $ mplayer out.avi
#
#                the first steps is needed because the mencoder programm
#                is not yet installed on the lab machines; you probably
#                need to add execute permission the shell-script
#
#                you can modify this script to read other files than PNGs
#                or to rename the output

# setting some options
opt=vbitrate=2160000:mbd=2:keyint=132:vqblur=1.0:cmp=2:subcmp=2:dia=2:o=mpv_flags=+mv0:last_pred=3

# encoding the movie using Microsofts MPEG4 encoder for compatibility
mencoder mf://*.png -mf fps=24 -ovc lavc -lavcopts vcodec=msmpeg4v2:$opt -o out.avi