#-------------------------------------------------
#
# Project created by QtCreator 2012-06-01T12:23:40
#
#-------------------------------------------------

QT       += core gui opengl

TARGET = VViewer

TEMPLATE = app


SOURCES += main.cpp \
    mainwindow.cpp \
    glwidget.cpp \
    camera.cpp \
    translator.cpp \
    zoomer.cpp \
    rotator.cpp \
    yimage.cpp \
    mesh.cpp \
    controller.cpp

HEADERS += \
    mainwindow.h \
    glwidget.h \
    camera.h \
    translator.h \
    zoomer.h \
    rotator.h \
    yimage.h \
    mesh.h \
    controller.h

LIBS    += -lGLU -lpng -L/home/evouga/OpenMesh-2.0.1/build/Build/lib/OpenMesh/ -lOpenMeshCore

INCLUDEPATH    += /home/evouga/eigen/ /home/evouga/OpenMesh-2.0.1/src/

QMAKE_CXXFLAGS += -g

FORMS += \
    mainwindow.ui
