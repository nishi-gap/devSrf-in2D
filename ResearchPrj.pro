QT += core gui
QT += openglwidgets
QT += opengl
QT += gui
QT += core
QT += widgets

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++20
CONFIG(release, debug|release) {
    CONFIG += optimize_full
 }
# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    cbox.cpp \
    fitting_plane.cpp \
    foldline.cpp \
    glwidget_2d.cpp \
    glwidget_3d.cpp \
    main.cpp \
    make3d.cpp \
    mathtool.cpp \
    originalbutton.cpp \
    setrulings.cpp \
    transform.cpp \
    widget.cpp

HEADERS += \
    cbox.h \
    fitting_plane.h \
    foldline.h \
    glwidget_2d.h \
    glwidget_3d.h \
    make3d.h \
    mathtool.h \
    originalbutton.h \
    setrulings.h \
    transform.h \
    widget.h \

FORMS += \
    widget.ui

INCLUDEPATH += $$PWD/includes
INCLUDEPATH += $$PWD/includes/nlopt/api


LIBS += $$PWD/libs/nlopt.lib



#QMAKE_POST_LINK += cp $$PWD/includes/dlls/nlopt.dll $$OUT_PWD/nlopt.dll
#dllファイルの追加：横メニューバーのプロジェクト→ビルドと実行メニューのうちメニューを選択→環境→PATHを編集してdllが入ったフォルダを追記するとできる
#dllファイルの追加：横メニューバーのプロジェクト→BUild Environmentのシステム環境変数を使用をクリック→Pathの中にdllのフォルダを追加
# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

DISTFILES += \
    Key_Assign.md
