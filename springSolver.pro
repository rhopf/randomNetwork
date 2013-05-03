#-------------------------------------------------
#
# Project created by QtCreator 2013-04-20T03:28:24
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

greaterThan(QT_MAJOR_VERSION, 4) {
    QT += printsupport
}

TARGET = springSolver
TEMPLATE = app

#INCLUDEPATH += /Library/LibrariesCPP/
#DEPENDPATH += /Library/LibrariesCPP/

INCLUDEPATH += D:/library/
DEPENDPATH += D:/library/

SOURCES += main.cpp\
        mainwindow.cpp \
    network.cpp \
    node.cpp \
    element.cpp

HEADERS  += mainwindow.h \
    network.h \
    node.h \
    element.h

FORMS    += mainwindow.ui
