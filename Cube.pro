QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

INCLUDEPATH += ./include

QMAKE_RPATHDIR += ../lib/win32 #可执行程序在部署环境运行时依赖的库文件的路径

#LIBS += -lglut32 -lopengl32 -lopencv_world400d

CONFIG(debug, debug|release): {
LIBS += -LE:/Courses/SoftwareTech/Cube/lib/win32\ -lglut -lglut32 -lopencv_world400d
} else:CONFIG(release, debug|release): {
LIBS += -LE:/Courses/SoftwareTech/Cube/lib/win32\ -lopencv_world400 -lglut32 -lglut
}

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    main.cpp \
    mainwindow.cpp \
    openglwidget.cpp \
    src/Deformer.cpp \
    src/LinearSystemSolver.cpp \
    src/Mesh.cpp \
    src/Parameterization.cpp

HEADERS += \
    mainwindow.h \
    openglwidget.h \
    src/Deformer.h \
    src/GLProjector.h \
    src/LinearSystemSolver.h \
    src/Mesh.h \
    src/Parameterization.h

FORMS += \
    mainwindow.ui \
    openglwidget.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

RESOURCES += \
    res.qrc \
