# Tout est compilé sous la forme d'une librairie: VFCLibrary.
# Cette librairie utilise la librairie glu (glu32 dans MinGW).
# Avec la commande suivante il est possible de compiler avec MinGW, et de 
# linker le tout.

mingw32-g++.exe main.cpp -L./ -lVFCLibrary -lglu32 -o main.exe
