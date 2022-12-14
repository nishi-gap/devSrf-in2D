#ifndef ORIGINALBUTTON_H
#define ORIGINALBUTTON_H

#include <QPushButton>
#include <QMouseEvent>
#include "setrulings.h"

class OriginalButton : public QPushButton
{
    Q_OBJECT
public:
    explicit  OriginalButton(QString &text, QWidget *parent = nullptr);
protected slots:

public slots:
signals:

};

class Btn4Crv : public OriginalButton
{
    Q_OBJECT
public:
    Btn4Crv(CRV *_crv, QString &text, QWidget *parent = nullptr);
protected:
    CRV *crv;
    void mousePressEvent(QMouseEvent *e);
    void mouseReleaseEvent(QMouseEvent *e);
signals:
    void clicked(Btn4Crv *click, QMouseEvent *e);
private:

};

#endif // ORIGINALBUTTON_H
