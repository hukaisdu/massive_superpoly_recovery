# Reference codes for the paper "Massive Superpoly Recovery with Nested Monomial Predictions"

## 1 Structure of the reference codes

The codes for Trivium, Grain128AEAD and Kreyvium are in the folders named by the algorithm names. The structures in every folder is similar. We take Trivium as an example.

main.cpp : the main function to evaluate the superpoly. 

trivium.cpp, trivium.h : the MILP models for trivium including the algorithms for expanding and evaluating $\mathsf{Coe}(\boldsymbol{s}^{\boldsymbol{t}}, \boldsymbol{x}^{\boldsymbol{u}})$ .

deg.cpp, deg.h : implementation of the numeric mapping to filter our some terms in solved-0 sets beforehand. 

Log and Term: the two folders are for the log file and the parts of the superpoly, respectively.

## 2 Usage of the codes

1. `cd Trivium(Kreyvium or Grain128AEAD)`
2. If there are no folders Log and Term, use `mkdir Log` and `mkdir Term`  to generate them. 
3. Edit the main.cpp to change the *cubeIndex* and *ROUND* to test the cube you are interested in.
4. type `make` to compile the executive file. Note you should install the Gurobi and edit the *makefile* to fit into the edition of the Gurobi
5. `./main` to run the program.



 

