#!/bin/bash

which_node="$1"

unison -ignorearchives -prefer newer "${which_node}"-quasielectron.prf
unison -ignorearchives -prefer newer "${which_node}"-clustercodes.prf
unison -ignorearchives -prefer ~/fzj/main-git/synth-dims "${which_node}"-synthdims.prf
unison -ignorearchives -prefer ~/fzj/main-git/review-practice-codes "${which_node}"-review.prf
unison -ignorearchives -prefer ~/fzj/main-git/other-funcs "${which_node}"-otherfuncs.prf
unison -ignorearchives -prefer ~/fzj/main-git/fci-transition "${which_node}"-fci-transition.prf
unison -ignorearchives -prefer ~/fzj/main-git/chemical-potential "${which_node}"-chemical-potential.prf
unison -ignorearchives -prefer ~/fzj/main-git/exact-diag "${which_node}"-exact-diag.prf
unison -ignorearchives -prefer ~/fzj/main-git/j1j2 "${which_node}"-j1j2.prf

if [ "$2" == "all" ]
then
	sh data-sync.sh
fi

