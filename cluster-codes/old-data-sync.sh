#!/bin/bash


# triangle of data storage with nodes at Desktop, PGI8, and MyLaptop
# 	MyLaptop will never be preferred
# This keeps two sets of backup
# First sync between Desktop and PGI8 preferring newest
# Both have full list of data and thus won't overwrite duplicates
# 	the only way a duplicate is made is if codes on each cluster
#	are running simultaneously and create files of the same name
#	before a sync takes place. Don't run same program on two clusters
# Then sync one with MyLaptop where Cluster node is preferred

ssh iff1500 "unison -prefer newer desktop-to-pgi8-data" 
echo "Synced Desktop with PGI8, now doing PGI8 to MyLaptop" 
unison -prefer newer pgi8-clusterdata 
exit


