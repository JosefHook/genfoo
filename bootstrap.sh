
echo " THIS SCRIPT ASSUMES YOU HAVE CHECKED OUT GenFoo in ~/GenFoo "
echo " It should look like '~/GenFoo/genfoo '"
echo "Bootstrap installation of GenFoo"
export GENFOO_HOME=`pwd`

mkdir Deps
cd Deps
mkdir src
mkdir include

wget "http://launchpad.net/dorsal/trunk/0.8.0/+download/dorsal-0.8.0.tar.bz2"
bunzip2 dorsal-0.8.0.tar.bz2
tar -xvf  dorsal-0.8.0.tar
pwd
cp -R ../Dorsal/GenFoo dorsal-0.8.0 
cp    ../Dorsal/dorsal.cfg dorsal-0.8.0/dorsal.cfg
cd dorsal-0.8.0

echo "Bootstrapping GenFoo "





./dorsal.sh GenFoo/platforms/sid.platform

echo "Done!"

