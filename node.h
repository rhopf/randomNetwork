#ifndef NODE_H
#define NODE_H

// includes
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <cmath>

// namespaces
using namespace std;
using namespace Eigen;

class Node
{
public:
    // constructor and destructor
    Node();
    ~Node();

    // initialize methods
    void initializeNode(vector<double>, int, int, bool); // ref=def
    void initializeNode(vector<double>, vector<double>, int, int, bool);

    // set methods
    void setReference(vector<double>);
    void setCurrent(vector<double>);
    void setPrescribedDisplacement(vector<double>);
    void setNodeIndex(int);
    void setBoundaryFlag(bool);
    void setDisplacement();

    // get methods
    vector<double> getReference();
    vector<double> getCurrent();
    vector<double> getDisplacement();
    vector<double> getPrescribedDisplacement();
    int getNodeIndex();
    bool getBoundaryFlag();
    int getDimension();

private:
    // spatial dimension
    int _DIM;
    int _nodeIndex;
    vector<double> _coordinateReference;
    vector<double> _coordinateCurrent;
    vector<double> _prescribedDisplacement;
    vector<double> _displacement;
    bool _boundaryFlag;
};

#endif // NODE_H
