# Import docker container
cat ucvm-slc4-docker.tar | docker import - ucvm-slc4-docker
# Check
docker images
