/****************************************************************************
** Meta object code from reading C++ file 'originalbutton.h'
**
** Created by: The Qt Meta Object Compiler version 68 (Qt 6.3.2)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include <memory>
#include "../../originalbutton.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'originalbutton.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 68
#error "This file was generated using the moc from 6.3.2. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
QT_WARNING_PUSH
QT_WARNING_DISABLE_DEPRECATED
struct qt_meta_stringdata_OriginalButton_t {
    uint offsetsAndSizes[2];
    char stringdata0[15];
};
#define QT_MOC_LITERAL(ofs, len) \
    uint(sizeof(qt_meta_stringdata_OriginalButton_t::offsetsAndSizes) + ofs), len 
static const qt_meta_stringdata_OriginalButton_t qt_meta_stringdata_OriginalButton = {
    {
        QT_MOC_LITERAL(0, 14)   // "OriginalButton"
    },
    "OriginalButton"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_OriginalButton[] = {

 // content:
      10,       // revision
       0,       // classname
       0,    0, // classinfo
       0,    0, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

       0        // eod
};

void OriginalButton::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    (void)_o;
    (void)_id;
    (void)_c;
    (void)_a;
}

const QMetaObject OriginalButton::staticMetaObject = { {
    QMetaObject::SuperData::link<QPushButton::staticMetaObject>(),
    qt_meta_stringdata_OriginalButton.offsetsAndSizes,
    qt_meta_data_OriginalButton,
    qt_static_metacall,
    nullptr,
qt_incomplete_metaTypeArray<qt_meta_stringdata_OriginalButton_t
, QtPrivate::TypeAndForceComplete<OriginalButton, std::true_type>



>,
    nullptr
} };


const QMetaObject *OriginalButton::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *OriginalButton::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_OriginalButton.stringdata0))
        return static_cast<void*>(this);
    return QPushButton::qt_metacast(_clname);
}

int OriginalButton::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QPushButton::qt_metacall(_c, _id, _a);
    return _id;
}
struct qt_meta_stringdata_Btn4Crv_t {
    uint offsetsAndSizes[14];
    char stringdata0[8];
    char stringdata1[8];
    char stringdata2[1];
    char stringdata3[9];
    char stringdata4[6];
    char stringdata5[13];
    char stringdata6[2];
};
#define QT_MOC_LITERAL(ofs, len) \
    uint(sizeof(qt_meta_stringdata_Btn4Crv_t::offsetsAndSizes) + ofs), len 
static const qt_meta_stringdata_Btn4Crv_t qt_meta_stringdata_Btn4Crv = {
    {
        QT_MOC_LITERAL(0, 7),  // "Btn4Crv"
        QT_MOC_LITERAL(8, 7),  // "clicked"
        QT_MOC_LITERAL(16, 0),  // ""
        QT_MOC_LITERAL(17, 8),  // "Btn4Crv*"
        QT_MOC_LITERAL(26, 5),  // "click"
        QT_MOC_LITERAL(32, 12),  // "QMouseEvent*"
        QT_MOC_LITERAL(45, 1)   // "e"
    },
    "Btn4Crv",
    "clicked",
    "",
    "Btn4Crv*",
    "click",
    "QMouseEvent*",
    "e"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_Btn4Crv[] = {

 // content:
      10,       // revision
       0,       // classname
       0,    0, // classinfo
       1,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags, initial metatype offsets
       1,    2,   20,    2, 0x06,    1 /* Public */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3, 0x80000000 | 5,    4,    6,

       0        // eod
};

void Btn4Crv::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        auto *_t = static_cast<Btn4Crv *>(_o);
        (void)_t;
        switch (_id) {
        case 0: _t->clicked((*reinterpret_cast< std::add_pointer_t<Btn4Crv*>>(_a[1])),(*reinterpret_cast< std::add_pointer_t<QMouseEvent*>>(_a[2]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<QMetaType *>(_a[0]) = QMetaType(); break;
        case 0:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<QMetaType *>(_a[0]) = QMetaType(); break;
            case 0:
                *reinterpret_cast<QMetaType *>(_a[0]) = QMetaType::fromType< Btn4Crv* >(); break;
            }
            break;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        {
            using _t = void (Btn4Crv::*)(Btn4Crv * , QMouseEvent * );
            if (_t _q_method = &Btn4Crv::clicked; *reinterpret_cast<_t *>(_a[1]) == _q_method) {
                *result = 0;
                return;
            }
        }
    }
}

const QMetaObject Btn4Crv::staticMetaObject = { {
    QMetaObject::SuperData::link<OriginalButton::staticMetaObject>(),
    qt_meta_stringdata_Btn4Crv.offsetsAndSizes,
    qt_meta_data_Btn4Crv,
    qt_static_metacall,
    nullptr,
qt_incomplete_metaTypeArray<qt_meta_stringdata_Btn4Crv_t
, QtPrivate::TypeAndForceComplete<Btn4Crv, std::true_type>, QtPrivate::TypeAndForceComplete<void, std::false_type>, QtPrivate::TypeAndForceComplete<Btn4Crv *, std::false_type>, QtPrivate::TypeAndForceComplete<QMouseEvent *, std::false_type>



>,
    nullptr
} };


const QMetaObject *Btn4Crv::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *Btn4Crv::qt_metacast(const char *_clname)
{
    if (!_clname) return nullptr;
    if (!strcmp(_clname, qt_meta_stringdata_Btn4Crv.stringdata0))
        return static_cast<void*>(this);
    return OriginalButton::qt_metacast(_clname);
}

int Btn4Crv::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = OriginalButton::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 1)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 1;
    }
    return _id;
}

// SIGNAL 0
void Btn4Crv::clicked(Btn4Crv * _t1, QMouseEvent * _t2)
{
    void *_a[] = { nullptr, const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t1))), const_cast<void*>(reinterpret_cast<const void*>(std::addressof(_t2))) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_WARNING_POP
QT_END_MOC_NAMESPACE
