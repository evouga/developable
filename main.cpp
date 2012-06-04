#include <QApplication>
#include <QDesktopWidget>

#include "mainwindow.h"
#include "controller.h"

 int main(int argc, char *argv[])
 {
     QApplication app(argc, argv);
     MainWindow window;
     Controller cont(window);
     window.setController(cont);

     int desktopArea = QApplication::desktop()->width() *
                      QApplication::desktop()->height();
     int widgetArea = window.width() * window.height();
     if (((float)widgetArea / (float)desktopArea) < 0.75f)
         window.show();
     else
         window.showMaximized();
     return app.exec();
 }
