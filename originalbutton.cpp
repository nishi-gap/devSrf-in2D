#include "originalbutton.h"

OriginalButton::OriginalButton(QString &text, QWidget *parent) : QPushButton(text, parent)
{
}


Btn4Crv::Btn4Crv(std::shared_ptr<CRV>& _crv, QString &text, QWidget *parent) : OriginalButton(text, parent)
{
    crv = _crv;
}

void Btn4Crv::mousePressEvent(QMouseEvent *e){
    emit clicked(this, e);
}

void Btn4Crv::mouseReleaseEvent(QMouseEvent *e){
    emit clicked(nullptr, e);
}

