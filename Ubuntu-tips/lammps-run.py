
1) To run lammps in background

Step A) create screen env
screen
hit enter

Step B) run the program
./submit.sh
hit enter

Step C) detach and run it in bckgrnd
press Ctrl + A then type D


2) To kill processes in screen

Step A)
# To view current list of screens
screen -ls

There are screens on:
	12346.pts-0.sumith-lab-1-l	(05/05/2024 05:11:32 PM)	(Detached)
	12242.pts-0.sumith-lab-1-l	(05/05/2024 05:09:28 PM)	(Detached)
2 Sockets in /run/screen/S-ydsumith.

Step B)
Get attached to the detached screen session
screen -r 12346.pts-0.sumith-lab-1-l

Step C)
Once connected to the session press Ctrl + A then type :quit and hit enter
