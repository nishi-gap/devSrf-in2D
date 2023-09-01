/********************************************************************************
** Form generated from reading UI file 'gtoolwnd.ui'
**
** Created by: Qt User Interface Compiler version 6.3.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_GTOOLWND_H
#define UI_GTOOLWND_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QWidget>
#include "gradationwidget.h"

QT_BEGIN_NAMESPACE

class Ui_GToolWnd
{
public:
    GradationWidget *Gview;

    void setupUi(QWidget *GToolWnd)
    {
        if (GToolWnd->objectName().isEmpty())
            GToolWnd->setObjectName(QString::fromUtf8("GToolWnd"));
        GToolWnd->resize(623, 542);
        Gview = new GradationWidget(GToolWnd);
        Gview->setObjectName(QString::fromUtf8("Gview"));
        Gview->setGeometry(QRect(20, 20, 581, 511));

        retranslateUi(GToolWnd);

        QMetaObject::connectSlotsByName(GToolWnd);
    } // setupUi

    void retranslateUi(QWidget *GToolWnd)
    {
        GToolWnd->setWindowTitle(QCoreApplication::translate("GToolWnd", "Form", nullptr));
    } // retranslateUi

};

namespace Ui {
    class GToolWnd: public Ui_GToolWnd {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_GTOOLWND_H
