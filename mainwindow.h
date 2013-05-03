#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QVector>
#include <QDebug>
#include <QtCore>
#include <QtGui>
#include <QFileDialog>
#include <QString>
#include <QFile>
#include <QDir>
#include <QLabel>
#include <QMessageBox>
#include <QTextStream>

// custom classes
#include "network.h"

// typedefs
typedef vector<double> DoubleList;
typedef vector<Node> NodeList;
typedef vector<Element> ElementList;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private slots:

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
