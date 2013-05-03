#ifndef ELEMENT_H
#define ELEMENT_H

// includes
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <QCustomPlot/qcustomplot.h>
#include "node.h"

class Element
{
public:
    // constructor and destructor
    Element(Node&,Node&,double,double,double);
    ~Element();

    // set methods
    void resetNodePointers(Node&,Node&);

    // get methods - element identification
    int getFirstIndex();
    int getSecondIndex();
    bool getFirstBoundaryFlag();
    bool getSecondBoundaryFlag();

    // get methods - coordinates
    vector<double> getFirstCoordinatesReference();
    vector<double> getFirstCoordinatesCurrent();
    vector<double> getSecondCoordinatesReference();
    vector<double> getSecondCoordinatesCurrent();

    // get methods - geometry
    double getReferenceLength();
    double getCurrentLength();
    double getDeltaD();

    // get methods - coefficients
    double getAlpha();
    double getAlphaD();

private:

    // nodes
    Node* _firstNode;
    Node* _secondNode;

    // stiffnesses
    double _cLow;
    double _cHigh;

    // cutoff displacement
    double _cutoff;
};

#endif // ELEMENT_H
