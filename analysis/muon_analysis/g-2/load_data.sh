#!/bin/bash

# check if the username is present, otherwise stop the script

if [[ $1 == "" ]]; then
  echo "ERROR --> you must provide your username for the gate as first argument."
  echo "          nothing provided, abort."
  exit
fi

SSH_PID=""

# Function to clean up on script exit
cleanup() {
  if [ -n "$SSH_PID" ]; then
    echo -e "\nStopping SSH tunnel..."
    umount $DATA_FOLDER
    if [ -z "$(ls -A "$DATA_FOLDER")" ]; then
      rm -r $DATA_FOLDER
    fi
    kill "$SSH_PID"
    echo "SSH tunnel stopped"
  fi
}

# Trap EXIT signal and call the cleanup function
trap cleanup EXIT

DATA_FOLDER="data"
USERNAME=$1

# Check if there is the folder data

if [ -d $DATA_FOLDER ]; then
  # if there is data check if it is empty or not
  if [ -n "$(ls $DATA_FOLDER)"  ]; then
   # if it is NOT empty stop the execution
   echo "FATAL ERROR --> the folder $DATA_FOLDER is already present and NOT empty"
   echo "                It can't be used as mount point, abort."
   exit
  fi
else
  # if there isn't a folder called $DATA_FOLDER we create it
  mkdir $DATA_FOLDER
fi



# Port forwarding (SSH tunneling)
echo "LOGIN TO GET ACCESS TO THE DATA USING SSH"
echo "you have 10 seconds to login"
echo ""
ssh -L 1111:147.162.55.96:22 $USERNAME@gate.pd.infn.it -N &
SSH_PID=$!

# Wait for a moment to ensure the tunnel is established
sleep 10

# Esegui il comando desiderato
# Sostituisci 'your_command' con il comando effettivo che desideri eseguire
echo "ACCESSING DATA ON REMOTE MACHINE AND MOUNTING FILES ON LOCAL FOLDER"
echo ""
sshfs -o ro -p 1111 lhcb@localhost:/home/lhcb/g-2/workdir/natale/ ./$DATA_FOLDER/ & 

wait
