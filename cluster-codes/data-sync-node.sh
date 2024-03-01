#!/bin/bash

which_node="$1"

unison -prefer newer "${which_node}"-clusterdata.prf
