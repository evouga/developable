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
    controller.cpp \
    developablemesh.cpp \
    schwarzdialog.cpp \
    glwidget2d.cpp \
    periodicmesh.cpp \
    developablemesh-energies.cpp \
    developablemesh-topology.cpp \
    mathutil.cpp \
    materialmesh.cpp \
    developablemesh-shapematch.cpp

HEADERS += \
    mainwindow.h \
    glwidget.h \
    camera.h \
    translator.h \
    zoomer.h \
    rotator.h \
    yimage.h \
    mesh.h \
    controller.h \
    developablemesh.h \
    schwarzdialog.h \
    autodiffTemplates.h \
    glwidget2d.h \
    periodicmesh.h \
    mathutil.h \
    materialmesh.h

LIBS    += -lGLU -lpng -L$${PWD}/ext/OpenMeshLib -lOpenMeshCore

INCLUDEPATH    += $${PWD}/ext/eigen/ $${PWD}/ext/OpenMesh/src $${PWD}/FADBAD++

macx {
    ## png from macports (X11 png didn't work) and GLU from OpenGL.framework
    LIBS += -L/opt/local/lib -L/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries
    INCLUDEPATH += /opt/local/include
}

DEFINES += HAVE_CSTDDEF

QMAKE_CXXFLAGS += -g

FORMS += \
    mainwindow.ui \
    schwarzdialog.ui




































