#!/bin/bash

#sudo apt update
#sudo apt-get install perl
#sudo cpan App::cpanminus
#sudo cpan install YAML
#sudo apt-get install libdbd-pg-perl
#sudo apt-get install libpq-dev


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
     ./Configure -des -Dprefix=$HOME/localperl
     make
     make test
     make install

myperl="$HOME/perl-5.14.4/perl"


echo "export PERL5LIB=\'$HOME/perl5/lib/perl5\'">>~/.bashrc && \
echo "export PERL_MB_OPT=\"--install_base '$HOME/perl5'\">>~/.bashrc && \
echo "export PERL_LOCAL_LIB_ROOT=$HOME/perl5">>~/.bashrc

#perl Makefile.PL --bootstrap && make test && make install
#$myperl -MCPAN -e 'install DBD::Pg'
#$myperl -MCPAN -e 'install DBI'

# Perl withuot root
 curl -L http://cpanmin.us | perl - App::cpanminus
~/perl5/bin/cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
~/perl5/bin/cpanm  Test::More
~/perl5/bin/cpanm  Expect
~/perl5/bin/cpanm  DBD::Pg
~/perl5/bin/cpanm  XML::Parser
~/perl5/bin/cpanm  XML::LibXML
~/perl5/bin/cpanm  LWP::Protocol::https



exit
