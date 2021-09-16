# massive_superpoly_recovery
The codes and recovered superpolies are for the paper "Massive Superpoly Recovery with Nested Monomial Predictions" accepted by Asiacrypt 2021. 

## Usage 
The "code" fold contains the source codes for Trivium, Grain128AEAD and Kreyvium.
For each stream cipher, you can edit main function in main.cpp to modify the cube and the attacked round. 
Then you need to edit makefile according to the path of your Gurobi installation. 
Finally use make to compile the source codes. 
