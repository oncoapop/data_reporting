#/bin/sh

clear
# Check operating system and process architecture

uname --all

# Check Processors
cat /proc/cpuinfo | grep processor

# Check disk space
df -h

# Check Memory
free -m

exit;
