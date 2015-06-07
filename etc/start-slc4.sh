# Run the relevant docker container 
docker run -i -t -v /cvmfs/cernvm-slc4.cern.ch/cvm3:/cernvm -v /cvmfs/aleph.cern.ch:/cvmfs/aleph.cern.ch -v /home/ganis/local/MIT:/MIT -v /cvmfs/sft.cern.ch:/cvmfs/sft.cern.ch ucvm-slc4-docker /init
 