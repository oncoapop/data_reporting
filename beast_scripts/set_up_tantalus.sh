#!/bin/bash

pip install pandas
pip install numpy
pip install click
pip install requests

echo " " >> ~/.bash_profile
echo "export TANTALUS_API_URL='https://tantalus.canadacentral.cloudapp.azure.com/api/'" >> ~/.bash_profile
echo "export TANTALUS_API_USERNAME='dyap'"  >> ~/.bash_profile
echo "export TANTALUS_API_PASSWORD='dyap'"  >> ~/.bash_profile

source ~/.bash_profile

