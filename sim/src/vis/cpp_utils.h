/** This file is part of the HYB simulation platform.
 *
 *  Copyright 2014- Finnish Meteorological Institute
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CPP_UTILS_H
#define CPP_UTILS_H

/** \brief Reference counting smart pointer that works like shared_ptr
 *         in C++0x and TR1 (but is simpler).
 *
 * Counts how many SharedPtrs point to object and automatically deletes
 * that object when this count goes to zero.
 */
template<class T>
class SharedPtr
{
public:
    SharedPtr(const SharedPtr<T>& other) {
        m_ref = other.m_ref;
        m_refCount = other.m_refCount;
        ++(*m_refCount);
    }
    explicit SharedPtr(T* ref) : m_ref(ref), m_refCount(new unsigned int(1)) { }
    SharedPtr<T>& operator=(const SharedPtr<T>& rhs) {
        if (&rhs != this) {
            abandon();
            m_ref = rhs.m_ref;
            m_refCount = rhs.m_refCount;
            ++(*m_refCount);
        }
        return *this;
    }
    bool operator==(const SharedPtr& other) const {
        return m_ref == other.m_ref;
    }
    T& operator*() const {
        return *m_ref;
    }
    T* operator->() const {
        return m_ref;
    }
    ~SharedPtr() {
        abandon();
    }
private:
    T* m_ref;
    unsigned int* m_refCount;

    void abandon() {
        --(*m_refCount);
        if (*m_refCount == 0) {
            delete m_ref;
            delete m_refCount;
        }
    }
};

/** \brief Smart pointer that has copy semantics.
 *
 * When pointer is copied, the object it points to is copied also.
 * So this works like "polymorphic variable". Object is also copied
 * when it is given to constructor of CopyPtr, and deleted when
 * destructor of CopyPtr is called.
 *
 * Warning: slicing occurs if real type of object is subclass of type
 * of reference given to constructor.
 */
template <class T>
class CopyPtr
{
public:
    CopyPtr(const CopyPtr<T>& other) : m_holder(other.m_holder->clone()) { }
    template <class ConcType>
    explicit CopyPtr(ConcType value)
        : m_holder(new ConcreteHolder<ConcType>(value)) { }
    ~CopyPtr() {
        delete m_holder;
    }
    CopyPtr<T>& operator=(const CopyPtr<T>& other) {
        if (&other != this) {
            delete m_holder;
            m_holder = other.m_holder->clone();
        }
        return *this;
    }
    T& operator*() const {
        return *(m_holder->getRef());
    }
    T* operator->() const {
        return m_holder->getRef();
    }
private:
    class HolderInterface;
    HolderInterface* m_holder;
    class HolderInterface
    {
    public:
        virtual ~HolderInterface() { };
        T* getRef() const {
            return m_ref;
        }
        virtual HolderInterface* clone() const = 0;
    protected:
        HolderInterface(T* ref) : m_ref(ref) { }
        T* const m_ref;
    };

    template <class ConcType>
    class ConcreteHolder : public HolderInterface
    {
    public:
        ConcreteHolder(const ConcType& value) : HolderInterface(new ConcType(value)) { }
        virtual ~ConcreteHolder() {
            delete this->m_ref;
        }
        virtual ConcreteHolder<ConcType>* clone() const {
            return new ConcreteHolder<ConcType>(*(static_cast<ConcType*>(this->m_ref)));
        }
    };
};

/** \brief SCopyPtr is like CopyPtr, but it can store only type "T"
 *         (not its subclass), and it works with incomplete
 *         ("forward declared") types.
 */
template <class T>
class SCopyPtr
{
    typedef T* (*copyf)(const T*);
    typedef void (*delf)(T*);
public:
    explicit SCopyPtr(T* ref) : m_copy(cpy), m_ref(ref), m_del(pdel) { }
    SCopyPtr(const SCopyPtr& other);
    ~SCopyPtr() {
        m_del(m_ref);
    }
    SCopyPtr operator=(const SCopyPtr<T>& other) {
        if (this != &other) {
            delete m_ref;
            m_ref = m_copy(other.m_ref);
        }
        return *this;
    }
    T& operator*() const {
        return *m_ref;
    }
    T* operator->() const {
        return m_ref;
    }
private:
    /* Copying and deleting is done through function pointers, which allows
       use with incomplete types. */
    copyf m_copy;
    T* m_ref;
    delf m_del;
    static void pdel(T* ref) {
        delete ref;
    }
    static T* cpy(const T* ref) {
        return new T(*ref);
    }
};

template <class T>
inline SCopyPtr<T>::SCopyPtr(const SCopyPtr<T>& other)
    : m_copy(other.m_copy), m_ref(m_copy(other.m_ref)), m_del(other.m_del) { }

/** \brief ClonePtr is like CopyPtr, but is copies by calling their
 *         "clone()" method (so they have to implement it), and
 *         slicing doesn't occur.
 *
 * clone() method has to return a clone (identical copy) of object
 * and real type of clone has to be same as real type of cloned objects
 * (or slicing occurs).
 */
template <class T>
class ClonePtr
{
public:
    explicit ClonePtr(T& ref) : m_ref(ref.clone()) { }
    ClonePtr(const ClonePtr& other) : m_ref(other.m_ref.clone()) { }
    ~ClonePtr() {
        delete &m_ref;
    }
    ClonePtr& operator=(const ClonePtr& other) {
        if (this != &other) {
            delete &m_ref;
            m_ref = other.clone();
        }
        return *this;
    }
    T& operator*() const {
        return m_ref;
    }
    T* operator->() const {
        return &m_ref;
    }
private:
    T& m_ref;
};

/** \brief Empty struct that is used by Function as "empty" template arguments.
 */
struct EmptyArg { };

/** \brief Class that works like "funtion pointer". Works like
 *         type "function" in C++Ox and TR1, but is simpler.
 *
 * Class Function can store function pointers and function objects.
 * Signature of function is defined by template parameters, for
 * example "Function<double (int, double>".
 *
 * Function call through this class contains virtual method call,
 * so it has some performance overhead.
 */
template <class Ret = EmptyArg, class Param1 = EmptyArg,
         class Param2 = EmptyArg, class Param3 = EmptyArg>
class Function;

template <class Ret>
class Function<Ret ()>
{
public:
    template <class Func> Function(Func f) : m_func(FuncWrapper<Func>(f)) { }
    Ret operator()() {
        return (*m_func)();
    }
private:
    struct FuncInterface;
    CopyPtr<FuncInterface> m_func;
    struct FuncInterface {
        virtual Ret operator()() = 0;
        virtual ~FuncInterface() { }
    };
    template <class Func>
    struct FuncWrapper : FuncInterface {
        FuncWrapper(Func f) : m_f(f) { }
        virtual Ret operator()() {
            return m_f();
        }
    private:
        Func m_f;
    };
};

template <class Ret, class Param1>
class Function<Ret (Param1)>
{
public:
    template <class Func> Function(Func f) : m_func(FuncWrapper<Func>(f)) { }
    Ret operator()(Param1 p1) const {
        return (*m_func)(p1);
    }
private:
    struct FuncInterface;
    CopyPtr<FuncInterface> m_func;
    struct FuncInterface {
        virtual Ret operator()(Param1 p1) = 0;
        virtual ~FuncInterface() { }
    };
    template <class Func>
    struct FuncWrapper : FuncInterface {
        FuncWrapper(Func f) : m_f(f) { }
        virtual Ret operator()(Param1 p1) {
            return m_f(p1);
        }
    private:
        Func m_f;
    };
};

template <class Ret, class Param1, class Param2>
class Function<Ret (Param1, Param2)>
{
public:
    template <class Func> Function(Func f) : m_func(FuncWrapper<Func>(f)) { }
    Ret operator()(Param1 p1, Param2 p2) {
        return (*m_func)(p1, p2);
    }
private:
    struct FuncInterface;
    CopyPtr<FuncInterface> m_func;
    struct FuncInterface {
        virtual Ret operator()(Param1 p1, Param2 p2) = 0;
        virtual ~FuncInterface() { }
    };
    template <class Func>
    struct FuncWrapper : FuncInterface {
        FuncWrapper(Func f) : m_f(f) { }
        virtual Ret operator()(Param1 p1, Param2 p2) {
            return m_f(p1, p2);
        }
    private:
        Func m_f;
    };
};

template <class Ret, class Param1, class Param2, class Param3>
class Function<Ret (Param1, Param2, Param3)>
{
public:
    template <class Func> Function(Func f) : m_func(FuncWrapper<Func>(f)) { }
    Ret operator()(Param1 p1, Param2 p2, Param3 p3) {
        return (*m_func)(p1, p2, p3);
    }
private:
    struct FuncInterface;
    CopyPtr<FuncInterface> m_func;
    struct FuncInterface {
        virtual Ret operator()(Param1 p1, Param2 p2, Param3 p3) = 0;
        virtual ~FuncInterface() { }
    };
    template <class Func>
    struct FuncWrapper : FuncInterface {
        FuncWrapper(Func f) : m_f(f) { }
        virtual Ret operator()(Param1 p1, Param2 p2, Param3 p3) {
            return m_f(p1, p2, p3);
        }
    private:
        Func m_f;
    };
};

#endif

