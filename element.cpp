#include "element.h"

// constructor and destructor
Element::Element(Node &node1,Node &node2,double cLow, double cHigh, double cutoff)
{
    // set nodes
    _firstNode = &node1;
    _secondNode = &node2;

    // stiffnesses
    _cLow = cLow;
    _cHigh = cHigh;

    // cutoff
    _cutoff = cutoff;
}

Element::~Element()
{
}

// SET methods -----------------------------------------------------------------
void Element::resetNodePointers(Node &node1,Node &node2)
{
    // set nodes
    _firstNode = &node1;
    _secondNode = &node2;

}

// GET methods -----------------------------------------------------------------

// element identification
int Element::getFirstIndex()
{
    return _firstNode->getNodeIndex();
}

int Element::getSecondIndex()
{
    return _secondNode->getNodeIndex();
}

bool Element::getFirstBoundaryFlag()
{
    return _firstNode->getBoundaryFlag();
}

bool Element::getSecondBoundaryFlag()
{
    return _secondNode->getBoundaryFlag();
}

// get methods - coordinates
vector<double> Element::getFirstCoordinatesReference()
{
    vector<double> temp = _firstNode->getReference();
    return temp;
}

vector<double> Element::getFirstCoordinatesCurrent()
{
    vector<double> temp = _firstNode->getCurrent();
    return temp;
}

vector<double> Element::getSecondCoordinatesReference()
{
    vector<double> temp = _secondNode->getReference();
    return temp;
}

vector<double> Element::getSecondCoordinatesCurrent()
{
    vector<double> temp = _secondNode->getCurrent();
    return temp;
}

// get methods - geometry
double Element::getReferenceLength()
{
    vector<double> a = _firstNode->getReference();
    vector<double> b = _secondNode->getReference();
    int DIM = _firstNode->getDimension();

    double sumOfSquares = 0;

    for (int i=0; i<DIM; i++)
    {
        sumOfSquares += (a[i]-b[i])*(a[i]-b[i]);
    }

    return sqrt(sumOfSquares);
}

double Element::getCurrentLength()
{
    vector<double> a = _firstNode->getCurrent();
    vector<double> b = _secondNode->getCurrent();
    int DIM = _firstNode->getDimension();

    double sumOfSquares = 0;

    for (int i=0; i<DIM; i++)
    {
        sumOfSquares += (a[i]-b[i])*(a[i]-b[i]);
    }

    return sqrt(sumOfSquares);
}

double Element::getDeltaD()
{
    double tempRef = getReferenceLength();
    double tempCur = getCurrentLength();

    return (tempCur - tempRef);
}

// get methods - coefficients
double Element::getAlpha()
{
    double lRef = getReferenceLength();
    double lCur = getCurrentLength();
    double deltaD = getDeltaD();

    double alpha;

    if (deltaD < _cutoff)
    {
        alpha = _cLow*(1 - lRef/lCur);
        return alpha;
    }

    else
    {
        alpha = _cHigh*(1 - lRef/lCur);
        alpha += _cutoff*(_cLow - _cHigh)/lCur;
        return alpha;
    }
}

double Element::getAlphaD()
{
    double lRef = getReferenceLength();
    double lCur = getCurrentLength();
    double deltaD = getDeltaD();

    double alphaD;

    if (deltaD < _cutoff)
    {
        alphaD = _cLow*lRef/pow(lCur,3);
        return alphaD;
    }

    else
    {
        alphaD = _cHigh*lRef/pow(lCur,3);
        alphaD += -_cutoff*(_cLow-_cHigh)/pow(lCur,3);
        return alphaD;
    }
}
