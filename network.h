#ifndef NETWORK_H
#define NETWORK_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/src/IterativeSolvers/GMRES.h>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>

// custom includes, node is included in element.h
#include "element.h"

using namespace Eigen;
using namespace std;

class Network
{
public:
    // type definitions
    typedef Triplet<double> T;

    // constructor and destructor
    Network();
    Network(vector<Node>,vector<Element>,int,int,int);
    ~Network();

    // general methods
    //void printEquations(); // DEBUG!
    //void printSolution(); // DEBUG!
    //void printForces(); // DEBUG!
    //double returnJacobi(int, int); // DEBUG
    //int returnNumberOfVariables();

    // public variables (network components)
    // REPLACE WITH POINTERS *ragaaaaaaahhh*!
    vector<Element> elements;
    vector<Node> nodes;

private:
    // private methods
    void _assembleLinearSystem(int);
    bool _solveLinearSystem();
    void _updateAll(bool);
    void _solveIncrement(int);
    void _solveNetwork();

    // helper methods
    double _getNorm(VectorXd);
    VectorXd _reverseSign(VectorXd,int);

    // network properties
    int _Nel;
    int _Np;
    int _Nbc;
    int _Nvar;
    int _Ninc;

    // boundary point indices
    vector<int> _boundaryIndices; // of length _Nbc !

    // equations of motion as F(x,f)=0
    VectorXd _F;

    // point loads (LP)
    VectorXd _f;

    // jacobian matrix
    SparseMatrix<double> _Jac;

    // triplet list for jacobian assembly
    vector<T> _tripletJac;

    // temporary solution vector: use to update nodes(current)
    // as well as point loads _f !
    VectorXd _sol;
    VectorXd _guess;

    // solver
    BiCGSTAB<SparseMatrix<double> > _solver;
    //LLC<<SparseMatrix<double> > _solver;

    // temporary variables

};

#endif // NETWORK_H
