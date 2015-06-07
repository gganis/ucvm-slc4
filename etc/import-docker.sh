# Import docker container
cat $HOME/local/docker/ucvm-slc4-docker.tar | docker import - ucvm-slc4-docker
# Check
docker images
