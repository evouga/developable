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
    meshcurvature.cpp \
    meshcontours.cpp

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
    meshcurvature.h \
    meshcontours.h

LIBS    += -lGLU -lpng -L$${PWD}/ext/OpenMesh/build/Build/lib/OpenMesh/ -lOpenMeshCore

INCLUDEPATH    += $${PWD}/ext/eigen/ $${PWD}/ext/OpenMesh/src

macx {
    ## png from macports (X11 png didn't work) and GLU from OpenGL.framework
    LIBS += -L/opt/local/lib -L/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries
    INCLUDEPATH += /opt/local/include
}

QMAKE_CXXFLAGS += -g

FORMS += \
    mainwindow.ui
