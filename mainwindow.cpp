#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // temp vars
    QStringList pieces;
    QString line;
    QString delimiterPattern(","), tempData;

    int dimension = 2;
    int nInc = 100;

    // ========================= NODES FILE (open and read nodes file) =========================
    QString nodesFile = QFileDialog::getOpenFileName(this, QObject::tr("Open Nodes File"),"",QObject::tr("Data (*.dat*);; Text (*.txt*)"));

    QFile fileNodes(nodesFile);
    fileNodes.open(QIODevice::ReadOnly | QIODevice::Text);

    // get number of lines
    QTextStream inNodes(&fileNodes);

    NodeList nodes;
    DoubleList coords(2);
    DoubleList coordsPrescribed(2);
    double x1,x2,x1t,x2t;
    bool bflag;

    int idxNode = 0;
    while (!inNodes.atEnd())
    {
        // read node
        line = inNodes.readLine();
        pieces = line.split(delimiterPattern);

        tempData = pieces.value(0);
        x1 = tempData.toDouble();
        tempData = pieces.value(1);
        x2 = tempData.toDouble();
        tempData = pieces.value(2);
        bflag = tempData.toDouble();
        bflag = (bool) bflag;
        tempData = pieces.value(3);
        x1t = tempData.toDouble();
        tempData = pieces.value(4);
        x2t = tempData.toDouble();

        coords[0] = x1;
        coords[1] = x2;

        coordsPrescribed[0] = x1t;
        coordsPrescribed[1] = x2t;

        // write nodes array
        if (bflag)
        {
            nodes.push_back(Node());
            nodes[idxNode].initializeNode(coords,idxNode,dimension,bflag);
            nodes[idxNode].setPrescribedDisplacement(coordsPrescribed);
        }
        else
        {
            nodes.push_back(Node());
            nodes[idxNode].initializeNode(coords,idxNode,dimension,bflag);
        }

        // update counter
        idxNode++;

    }

    // total number of nodes
    int nNodes = idxNode;

    fileNodes.close();


    // ========================= ELEMENTS FILE (open and read nodes file) =========================

    QString elementsFile = QFileDialog::getOpenFileName(this, QObject::tr("Open Elements File"),"",QObject::tr("Data (*.dat*);; Text (*.txt*)"));

    QFile fileElements(elementsFile);
    fileElements.open(QIODevice::ReadOnly | QIODevice::Text);

    // get number of lines
    QTextStream inElements(&fileElements);

    ElementList elements;
    int idxN1,idxN2;
    double cl,ch,cutoff;

    int idxElement = 0;
    while (!inElements.atEnd())
    {
        // read node
        line = inElements.readLine();
        pieces = line.split(delimiterPattern);

        tempData = pieces.value(0);
        idxN1 = tempData.toDouble();
        idxN1 = (int) idxN1;
        tempData = pieces.value(1);
        idxN2 = tempData.toDouble();
        idxN2 = (int) idxN2;
        tempData = pieces.value(2);
        cl = tempData.toDouble();
        tempData = pieces.value(3);
        ch = tempData.toDouble();
        tempData = pieces.value(4);
        cutoff = tempData.toDouble();

        // write elements array
        elements.push_back(Element(nodes[idxN1],nodes[idxN2],cl,ch,cutoff));

        // update counter
        idxElement++;

    }

    // total number of nodes
    int nElements = idxElement;

    fileElements.close();


    // ========================= INITIALIZE NETWORK AND SOLVE =========================

    Network network(nodes,elements,nNodes,nElements,nInc);

}

MainWindow::~MainWindow()
{
    delete ui;
}

