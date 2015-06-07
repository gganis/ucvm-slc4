# File to be sourced to setup the ALEPH env

# Main ALEPH software directory
export ALEPH=/cvmfs/aleph.cern.ch/i386_redhat42

# ALEPH Banks and database
export BANKALFMT=$ALEPH/dbase/bankal.fmt
export DBASBANK=$ALEPH/dbase/dbas.bank
export ADBSCONS=$ALEPH/dbase/adbscons.daf

# Job cards
export ALPHACARDS=MIT.cards
# Not needed: required perhaps by ALPHARUN
# export STAGE_HOST=local
