QT += core gui
QT += openglwidgets
QT += opengl
QT += gui
QT += core
QT += widgets

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++20

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    #arcball.cpp \
    cbox.cpp \
    foldline.cpp \
    glwidget_2d.cpp \
    glwidget_3d.cpp \
    gradationwidget.cpp \
    gtoolwnd.cpp \
    main.cpp \
    make3d.cpp \
    mathtool.cpp \
    originalbutton.cpp \
    setrulings.cpp \
    widget.cpp

HEADERS += \
    #arcball.h \
    cbox.hpp \
    foldline.hpp \
    glwidget_2d.hpp \
    glwidget_3d.hpp \
    gradationwidget.hpp \
    gtoolwnd.hpp \
    make3d.hpp \
    mathtool.hpp \
    originalbutton.hpp \
    setrulings.hpp \
    widget.hpp

FORMS += \
    gtoolwnd.ui \
    widget.ui

INCLUDEPATH += $$PWD/includes \
                $$PWD/includes/nlopt/api

QMAKE_RPATHDIR += $$PWD/includes/dlls
LIBS += $$PWD/includes/libs/nlopt.lib
LIBS += $$PWD/includes/dlls/nlopt.dll

#dllファイルの追加：横メニューバーのプロジェクト→ビルドと実行メニューのうちメニューを選択→環境→PATHを編集してdllが入ったフォルダを追記するとできる
#dllファイルの追加：横メニューバーのプロジェクト→BUild Environmentのシステム環境変数を使用をクリック→Pathの中にdllのフォルダを追加
# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
