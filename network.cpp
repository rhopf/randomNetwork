#include "network.h"

// constructor and destructor
Network::Network()
{
}
Network::Network(vector<Node> tempNodes, vector<Element> tempElements,int Np, int Nel, int Ninc)
{
    // network properties
    _Nel = Nel;
    _Np = Np;
    _Ninc = Ninc;

    bool existKBC;
    _Nbc = 0; // number of points (=corresponds to two var each) with kinematic constraints imposed

    // add nodes -------------------------------------
    for (int i=0; i<_Np; i++)
    {
        nodes.push_back(tempNodes[i]);
        existKBC = tempNodes[i].getBoundaryFlag();
        if (existKBC)
        {
            _Nbc++;
            _boundaryIndices.push_back(i);
        }
    }

    // add elements -----------------------------------
    int idx1,idx2;
    for (int i=0; i<_Nel; i++)
    {
        elements.push_back(tempElements[i]);

        // reset node pointers
        idx1 = elements[i].getFirstIndex();
        idx2 = elements[i].getSecondIndex();
        elements[i].resetNodePointers(nodes[idx1],nodes[idx2]);
    }

    // get total number of variables
    _Nvar = 2*_Np + 2*_Nbc;

    // ======= initialize equations of motion and point loads to zero =======
    _F.resize(_Nvar);
    _F.setZero(_Nvar);
    _f.resize(2*_Nbc);
    _f.setZero(2*_Nbc);
    _sol.resize(_Nvar);
    _sol.setZero(_Nvar);
    _guess.resize(_Nvar);
    _guess.setZero(_Nvar);

    // initialize triplet list
    _tripletJac.reserve(10*_Np);

    // initialize jacobian matrix
    _Jac.resize(_Nvar,_Nvar);
    _Jac.reserve(10*_Np);

    // solve the network for all increments
    _solveNetwork();
}

Network::~Network()
{
}

// general methods --------------------------------------------------------------------


// solve equations: ONLY 1 INCREMENT IMPLEMENTED YET!
// TODO: update all quantities, TOL limit based on norm(n(n+1)-x(n)) etc.


// PRIVATE functions ------------------------------------------------------------
void Network::_assembleLinearSystem(int currentIncrement)
{
    // reset equations of motion and jacobian matrix and triplets
    _F.setZero(_Nvar);
    _Jac.setZero(); // UGLY, FIND ANOTHER WAY!?
    _tripletJac.clear();

    // current element-nodes (node index "Left" and "Right")
    int idxL;
    int idxR;

    // current coordinates
    vector<double> xL;
    vector<double> xR;

    // alpha and beta coefficients
    double alpha;
    double alphaD;

    // temporary Jacobian matrix contributions
    // "11" / +22" refers to derivates of the same coordinate dimension
    // "sym" refers to derivatives of the same node index
    double J11sym;
    double J12sym;
    double J22sym;
    double J11asym;
    double J12asym;
    double J22asym;

    // PART I - ELEMENTS: loop through all elements.
    // This will yield the first F(1:2Np) components. Adding all 16 jacobian contributions also.
    for (int elIdx=0; elIdx<_Nel; elIdx++)
    {

        // get nodal indices
        idxL = elements[elIdx].getFirstIndex();
        idxR = elements[elIdx].getSecondIndex();

        // get alpha and beta coefficients
        alpha = elements[elIdx].getAlpha();
        alphaD = elements[elIdx].getAlphaD();

        // get current coordinates
        xL = elements[elIdx].getFirstCoordinatesCurrent();
        xR = elements[elIdx].getSecondCoordinatesCurrent();

        // update guess
        _guess[2*idxL]   = xL[0];
        _guess[2*idxL+1] = xL[1];
        _guess[2*idxR]   = xR[0];
        _guess[2*idxR+1] = xR[1];

        // ------------------ update equilibrium equations -------------------
        _F[2*idxL]   += -alpha*(xL[0]-xR[0]);
        _F[2*idxL+1] += -alpha*(xL[1]-xR[1]);
        _F[2*idxR]   += alpha*(xL[0]-xR[0]);
        _F[2*idxR+1] += alpha*(xL[1]-xR[1]);

        // --------------------- update jacobian matrix ----------------------
        // contributions based on symmetry and coordinate dimension
        J11sym = -alphaD*(xL[0]-xR[0])*(xL[0]-xR[0]) - alpha;
        J22sym = -alphaD*(xL[1]-xR[1])*(xL[1]-xR[1]) - alpha;
        J12sym = -alphaD*(xL[0]-xR[0])*(xL[1]-xR[1]);
        J11asym = -J11sym;
        J22asym = -J22sym;
        J12asym = -J12sym;

        // update triplets
        _tripletJac.push_back(T(2*idxL,     2*idxL,     J11sym));
        _tripletJac.push_back(T(2*idxL,     2*idxL+1,   J12sym));
        _tripletJac.push_back(T(2*idxL,     2*idxR,     J11asym));
        _tripletJac.push_back(T(2*idxL,     2*idxR+1,   J12asym));
        _tripletJac.push_back(T(2*idxL+1,   2*idxL,     J12sym));
        _tripletJac.push_back(T(2*idxL+1,   2*idxL+1,   J22sym));
        _tripletJac.push_back(T(2*idxL+1,   2*idxR,     J22asym));
        _tripletJac.push_back(T(2*idxL+1,   2*idxR+1,   J12asym));
        _tripletJac.push_back(T(2*idxR,     2*idxL,     J11asym));
        _tripletJac.push_back(T(2*idxR,     2*idxL+1,   J12asym));
        _tripletJac.push_back(T(2*idxR,     2*idxR,     J11sym));
        _tripletJac.push_back(T(2*idxR,     2*idxR+1,   J12sym));
        _tripletJac.push_back(T(2*idxR+1,   2*idxL,     J12asym));
        _tripletJac.push_back(T(2*idxR+1,   2*idxL+1,   J22asym));
        _tripletJac.push_back(T(2*idxR+1,   2*idxR,     J12sym));
        _tripletJac.push_back(T(2*idxR+1,   2*idxR+1,   J22sym));

    }

    // PART II - NODES: Complete F(x,f) with the all thats related to d()/df, i.e. all boundary
    // clamping constraints f*(x-xPrescribed). This will yield the last F((2Np+1):4Np) components
    // The procedure checks every node for a boundary constraint flag. If a flag is found:
    // F(2i) = x1Current-x1Prescribed, F(2i+1) = x2Current-x2Prescribed, i in ((2Np+1):4Np)

    vector<double> currentPos;
    vector<double> referencePos;
    vector<double> displacement(2);
    vector<double> prescribedDispl;

    // boundary index in global node numbering
    int globalIdx = 0;

    // loop through all boundary nodes
    for (int bcIdx=0; bcIdx<_Nbc; bcIdx++)
    {
        // get current global node index
        globalIdx = _boundaryIndices[bcIdx];

        // get current position and fixed constraint of node
        currentPos = nodes[globalIdx].getCurrent();
        referencePos = nodes[globalIdx].getReference();
        prescribedDispl = nodes[globalIdx].getPrescribedDisplacement();
        displacement[0] = prescribedDispl[0] - referencePos[0];
        displacement[1] = prescribedDispl[1] - referencePos[1];
        _F[2*_Np + 2*bcIdx] = currentPos[0] - (referencePos[0] + displacement[0]*currentIncrement/_Ninc);
        _F[2*_Np + 2*bcIdx+1] = currentPos[1] - (referencePos[1] + displacement[1]*currentIncrement/_Ninc);

        // make triplets with dFi/dfj for Jacobian
        _tripletJac.push_back(T(2*globalIdx,2*_Np + 2*bcIdx,1));
        _tripletJac.push_back(T(2*globalIdx+1,2*_Np + 2*bcIdx+1,1));

        // make triplets with dF_n+i/dui for Jacobian
        _tripletJac.push_back(T(2*_Np + 2*bcIdx,2*globalIdx,1));
        _tripletJac.push_back(T(2*_Np + 2*bcIdx + 1,2*globalIdx + 1,1));

        // update eqn with constraint forces
        _F[2*globalIdx] += _f[2*bcIdx];
        _F[2*globalIdx+1] += _f[2*bcIdx+1];
    }

    // assemble jacobian matrix from triplets
    _Jac.setFromTriplets(_tripletJac.begin(),_tripletJac.end());
}

bool Network::_solveLinearSystem()
{
    // solver settings
    _solver.setMaxIterations(100000);
    _solver.setTolerance(1E-2);
    _solver.compute(_Jac);
    _solver.preconditioner();

    // DEBUG
    MatrixXd mat(_Nvar,_Nvar);
    mat = _Jac.toDense();

    // decomp ok?
    if (_solver.info() == Success)
    {
        //qDebug() << "yay decomp OK :)";
    }
    else
    {
        qDebug() << "shit decomp fail :(";
        return false;
    }

    // reverse the sign of _F
    _F = _reverseSign(_F,_Nvar);

    // SOLVE SYSTEM
    //_sol=_solver.solveWithGuess(_F,_guess);
    _sol = mat.colPivHouseholderQr().solve(_F);

    // converged?
    if (_solver.info() == Success)
    {
        qDebug() << "Yay, converged. :) ...after " << _solver.iterations() << " iterations.";
        return true;
    }
    else
    {
        qDebug() << "Shit, diverged. :(";
        return false;
    }
}

void Network::_updateAll(bool converged)
{
    if (converged)
    {
        vector<double> tempNew, tempCurrent;

        // update nodal coordinates
        for (int i=0; i<_Np; i++)
        {
            tempCurrent.clear();
            tempNew.clear();
            tempCurrent = nodes[i].getCurrent();
            tempNew.push_back(tempCurrent[0] + _sol[2*i]);
            tempNew.push_back(tempCurrent[1] + _sol[2*i+1]);
            nodes[i].setCurrent(tempNew);
        }

        // update nodal forces
        for (int i=0; i<2*_Nbc; i++)
        {
            _f[i] += _sol[2*_Np + i];
        }
    }
}

void Network::_solveIncrement(int currentIncrement)
{
    // temporary variables
    bool converged;

    // assemble equilibrium equations
    _assembleLinearSystem(currentIncrement);

    // solve equations
    converged = _solveLinearSystem();

    // update everything if converged
    _updateAll(converged);

    // norm of eqn of mot
    double fNorm = _getNorm(_F);
    //qDebug() << "norm is now: " << fNorm;

    // tolerance
    double TOL = 1E-4;

    if (fNorm>TOL && converged)
    {
        _solveIncrement(currentIncrement);
    }
    if (!converged || !(_sol[0]==_sol[0]))
    {
        qDebug() << "error solving eqn";
    }
    else
    {
        //qDebug() << "tolerance reached :)";
    }
}

void Network::_solveNetwork()
{
    // file output
    vector<double> coords(2);
    vector<double> forces(2);
    int gIdx;
    ofstream outNodes, outForces, jacobi;
    /*
    outNodes.open("/Users/rhopf/Dropbox/div_work/network_model/outputNodes.dat");
    outNodes.close();
    outNodes.open("/Users/rhopf/Dropbox/div_work/network_model/outputNodes.dat",ios::app);
    outForces.open("/Users/rhopf/Dropbox/div_work/network_model/outputForces.dat");
    outForces.close();
    outForces.open("/Users/rhopf/Dropbox/div_work/network_model/outputForces.dat",ios::app);
    */
    outNodes.open("D:/dropbox/Dropbox/div_work/network_model/outputNodes.dat");
    outNodes.close();
    jacobi.open("D:/dropbox/Dropbox/div_work/network_model/jacobi.dat");
    jacobi.close();
    jacobi.open("D:/dropbox/Dropbox/div_work/network_model/jacobi.dat");
    outNodes.open("D:/dropbox/Dropbox/div_work/network_model/outputNodes.dat",ios::app);
    outForces.open("D:/dropbox/Dropbox/div_work/network_model/outputForces.dat");
    outForces.close();
    outForces.open("D:/dropbox/Dropbox/div_work/network_model/outputForces.dat",ios::app);

    for (int inc=1; inc<=_Ninc; inc++)
    {
        _solveIncrement(inc);

        if (inc==1)
        {
            for (int i=0; i<_Nvar; i++)
            {
                for (int j=0; j <_Nvar; j++)
                {
                    jacobi << _Jac.coeff(i,j) << "\t";
                }

                jacobi << endl;
            }
        }

        // write file output: nodes
        for (int i=0; i<_Np; i++)
        {
            coords = nodes[i].getCurrent();
            outNodes << i << "," << coords[0] << "," << coords[1] << endl;
        }

        // write file output: forces
        for (int i=0; i<_Nbc; i++)
        {
            gIdx = _boundaryIndices[i];
            forces[0] = _f[2*i];
            forces[1] = _f[2*i+1];
            outForces << gIdx << "," << forces[0] << "," << forces[1] << endl;
        }
    }

    outNodes.close();
    outForces.close();
    jacobi.close();
}

double Network::_getNorm(VectorXd vec)
{
    double vecNorm = 0;
    int length = vec.rows();

    for (int i=0; i<length; i++)
    {
        vecNorm += vec(i)*vec(i);
    }

    vecNorm = sqrt(vecNorm);

    return vecNorm;
}

VectorXd Network::_reverseSign(VectorXd V, int size)
{

    VectorXd temp(size);

    for (int i=0; i<size; i++)
    {
        temp[i] = -V[i];
    }

    return temp;
}

// DEBUG methods ----------------------------------------------------------------

/*
 *void Network::printEquations()
{
    for (int i=0; i<_Nvar; i++)
    {
        qDebug() << _F[i];
    }
}

void Network::printSolution()
{
    for (int i=0; i<_Nvar; i++)
    {
        qDebug() << _sol[i];
    }
}

void Network::printForces()
{
    int gIdx;

    for (int i=0; i<_Nbc; i++)
    {
        gIdx = _boundaryIndices[i];
        qDebug() << "Force f" << gIdx << "x = " << _f[2*i];
        qDebug() << "Force f" << gIdx << "y = " << _f[2*i+1];
    }
}

double Network::returnJacobi(int i, int j)
{
    return _Jac.coeff(i,j);
}

int Network::returnNumberOfVariables()
{
    return _Nvar;
}
*/
