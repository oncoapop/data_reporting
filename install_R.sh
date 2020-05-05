#~/bin/sh

echo "Script to install R.."

echo "Installing key..."
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

echo "Add in repo..."
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

echo "Updating..."
sudo apt update

echo "Installing R..."
sudo apt install r-base

echo "Testing R..."
sudo -i R

