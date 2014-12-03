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
    projectionnlp.cpp \
    springmesh.cpp

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
    materialmesh.h \
    projectionnlp.h \
    springmesh.h

<<<<<<< HEAD
LIBS    += -lGLU -lpng -L$${PWD}/ext/OpenMesh/lib -lOpenMeshCored -L/Users/work/libs/ipopt//lib -L/usr/local/lib/gcc/x86_64-apple-darwin13.1.0/4.9.0 -L/usr/local/lib/gcc/x86_64-apple-darwin13.1.0/4.9.0/../../.. -L/Users/work/libs/ipopt//lib -lipopt -framework vecLib -lm -ldl -lcoinmumps -framework vecLib -lgfortran -lSystem -lquadmath -lm -lcoinmetis

INCLUDEPATH    += $${PWD}/ext/ $${PWD}/ext/OpenMesh/include $${PWD}/FADBAD++
=======
LIBS    += -lGLU -lpng -L/Users/work/libs/newopenmesh/lib/ -lOpenMeshCore -lOpenMeshCored -lgfortran -lblas -llapack -L/Users/work/libs/ipopt/lib -L/usr/local/lib/gcc/x86_64-apple-darwin13.1.0/4.9.0 -L/usr/local/lib/gcc/x86_64-apple-darwin13.1.0/4.9.0/../../.. -L/Users/work/libs/ipopt/lib -lipopt -framework vecLib -lm -ldl -lcoinmumps -framework vecLib -lgfortran -lSystem -lquadmath -lm -lcoinmetis
# -lOpenMeshCored -lOpenMeshToolsd
INCLUDEPATH    += $${PWD}/ext/ $${PWD}/ext/OpenMesh/include $${PWD}/FADBAD
>>>>>>> b9ebcec3dbf38739862d535980a9d29f5a82a734

macx {
    ## png from macports (X11 png didn't work) and GLU from OpenGL.framework
    LIBS += -L/opt/local/lib -L/System/Library/Frameworks/OpenGL.framework/Versions/A/Libraries
    INCLUDEPATH += /opt/local/include/
}

DEFINES += HAVE_CSTDDEF OM_STATIC_BUILD

QMAKE_CXXFLAGS += -g
QMAKE_MACOSX_DEPLOYMENT_TARGET = 10.9

FORMS += \
    mainwindow.ui \
    schwarzdialog.ui




































