#!/bin/bash

# Port forwarding (SSH tunneling)
echo "LOGIN TO GET ACCESS TO THE DATA USING SSH"
echo "you have 10 seconds to login"
echo ""
ssh -L 4000:147.162.55.96:22 salimben@gate.pd.infn.it -N &

# Wait for a moment to ensure the tunnel is established
sleep 10

# Esegui il comando desiderato
# Sostituisci 'your_command' con il comando effettivo che desideri eseguire
echo "ACCESSING DATA ON REMOTE MACHINE AND MOUNTING FILES ON LOCAL FOLDER"
echo ""
sshfs -o ro -p 4000 lhcb@localhost:/home/lhcb/g-2/workdir/natale/ ./data/
