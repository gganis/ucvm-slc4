#
# Configuring cernvm-slc4.cern.ch
cwddir=`pwd`
# Fix /etc/cvmfs/config.d/cernvm-slc4.cern.ch.conf
cd /etc/cvmfs/
if test -f config.d/cernvm-slc4.cern.ch.conf ; then
   mv config.d/cernvm-slc4.cern.ch.conf config.d/cernvm-slc4.cern.ch.conf.backup
fi
cat > config.d/cernvm-slc4.cern.ch.conf <<EOF
CVMFS_SERVER_URL=http://hepvm.cern.ch/cvmfs/cernvm-slc4.cern.ch
CVMFS_KEYS_DIR=/etc/cvmfs/keys
EOF

# Get the key
cd keys
wget https://raw.githubusercontent.com/cernvm/cernvm-micro/master/packages.d/cvmfs-config/etc/cvmfs/keys/cernvm-slc4.cern.ch/cernvm-slc4.cern.ch.pub
# Display it
cat cernvm-slc4.cern.ch.pub 

# Access the repository
cd $cwddir
ls /cvmfs/cernvm-slc4.cern.ch
cvm3  new_repository  update-packs
