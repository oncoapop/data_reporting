#~/bin/sh

echo "Script to install R.."

#echo "Installing key..."
#sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

#echo "Add in repo..."
#sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

echo "Updating..."
sudo apt update

#echo "Installing dependencies..."
#sudo apt install r-base

cd /opt
# install prereq
echo "Installing readline..."
sudo wget -c -N https://ftp.gnu.org/gnu/readline/readline-7.0.tar.gz
sudo tar -xzf readline-7.0.tar.gz
cd readline-7.0/
sudo ./configure --prefix=`pwd`
sudo make
sudo make install

sudo wget https://cran.r-project.org/src/base/R-4/R-4.0.0.tar.gz
cd /opt/R-4.0.0
#sudo nano config.site
#CPPFLAGS='-I/opt/readline-7.0/include/'
#LDFLAGS='-L/opt/readline-7.0/lib/'

export LD_LIBRARY_PATH=/opt/readline-7.0/lib/

sudo ./configure --prefix=/opt/R-4.0.0 --with-x=yes --enable-R-shlib=yes \
	--with-cairo=yes --with-libpng=yes

make
sudo make install


echo "Testing R..."
sudo -i R

