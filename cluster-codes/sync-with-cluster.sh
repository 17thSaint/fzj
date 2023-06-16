#!/bin/bash

which_node="$1"

if [ "$2" == "all" ]
then
	unison -ignorearchives -prefer ~/fzj/main-git/other-funcs "${which_node}"-otherfuncs.prf
	unison -ignorearchives -prefer ~/fzj/main-git/review-practice-codes "${which_node}"-review.prf
	unison -ignorearchives -prefer ~/fzj/main-git/synth-dims "${which_node}"-synthdims.prf
	unison -ignorearchives -prefer newer "${which_node}"-clustercodes.prf
fi

