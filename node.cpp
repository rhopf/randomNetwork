#include "node.h"

// constructor and destructor
Node::Node()
{
}

Node::~Node()
{
}

// initialize node in reference state (ref=def)
void Node::initializeNode(vector<double> coordinates,int nIndex, int dimension, bool boundary)
{
    // set dimension
    _DIM = dimension;

    // loop over dimension, set coordinates
    for (int i=0; i<_DIM; i++)
    {
        // initial coordinates, ref=def
        _coordinateReference.push_back(coordinates[i]);
        _coordinateCurrent.push_back(coordinates[i]);

        // set initial displacement to zero
        _displacement.push_back(0);
    }

    // boundary flag
    _boundaryFlag = boundary;

    // node index
    _nodeIndex = nIndex;
}
// initialize node (ref=def)
void Node::initializeNode(vector<double> coordinatesReference, vector<double> coordinatesCurrent, int nIndex, int dimension, bool boundary)
{
    // set dimension
    _DIM = dimension;

    // loop over dimension, set coordinates and displacement
    for (int i=0; i<_DIM; i++)
    {
        _coordinateReference.push_back(coordinatesReference[i]);
        _coordinateCurrent.push_back(coordinatesCurrent[i]);

        _displacement.push_back(coordinatesCurrent[i]-coordinatesReference[i]);
    }

    // boundary flag
    _boundaryFlag = boundary;

    // node index
    _nodeIndex = nIndex;
}

// SET methods ---------------------------------------------------------------------------

// set reference (debug)
void Node::setReference(vector<double> coordinates)
{
    // loop over dimension, set coordinates
    for (int i=0; i<_DIM; i++)
    {
        _coordinateReference[i] = coordinates[i];
    }
}

// set (update) current coordinates
void Node::setCurrent(vector<double> coordinates)
{
    // loop over dimension, set coordinates
    for (int i=0; i<_DIM; i++)
    {
        _coordinateCurrent[i] = coordinates[i];
    }
}

// set node index
void Node::setNodeIndex(int nIndex)
{
    _nodeIndex = nIndex;
}

// set displacement
void Node::setDisplacement()
{
    // loop over dimension, set coordinates
    for (int i=0; i<_DIM; i++)
    {
        _displacement[i] = _coordinateCurrent[i]-_coordinateReference[i];
    }
}

// set boundary flag
void Node::setBoundaryFlag(bool boundaryFlag)
{
    _boundaryFlag = boundaryFlag;
}

void Node::setPrescribedDisplacement(vector<double> tempDisplacement)
{
    if (_boundaryFlag)
    {
        _prescribedDisplacement = tempDisplacement;
    }
}

// GET methods ---------------------------------------------------------------------------

// get reference coordinates
vector<double> Node::getReference()
{
    return _coordinateReference;
}

// get current coordinates
vector<double> Node::getCurrent()
{
    return _coordinateCurrent;
}

// get displacements
vector<double> Node::getDisplacement()
{
    return _displacement;
}

// get node index
int Node::getNodeIndex()
{
    return _nodeIndex;
}

// get boundary flag
bool Node::getBoundaryFlag()
{
    return _boundaryFlag;
}

// get dimension
int Node::getDimension()
{
    return _DIM;
}

vector<double> Node::getPrescribedDisplacement()
{
    return _prescribedDisplacement;
}

// PRIVATE methods ---------------------------------------------------------------------------
