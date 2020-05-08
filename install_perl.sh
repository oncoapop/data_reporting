#!/bin/bash

#sudo apt update
#sudo apt-get install perl
#sudo cpan App::cpanminus
#sudo cpan install YAML
sudo apt-get install libdbd-pg-perl
sudo apt-get install libpq-dev
# Install older local version
# https://www.cpan.org/src/
src="https://www.cpan.org/src/5.0/perl-5.14.4.tar.bz2"
file="perl-5.14.4.tar.bz2"
file2="perl-5.14.4.tar"
ver="perl-5.14.4"
echo "file="$file
echo "ver="$ver
rm -f $file
cd $HOME
     wget $src
	bzip2 -dk $file
     	tar -xf $file2

     cd $ver
     sudo ./Configure -des -Dprefix=$HOME/localperl
     sudo make
     sudo make test
     sudo make install

myperl="$HOME/perl-5.14.4/perl"

sudo $myperl -MCPAN -e 'install DBD::Pg'
sudo $myperl -MCPAN -e 'install DBI'

exit

