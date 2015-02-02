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

#ifndef CONTAINER_H
#define CONTAINER_H

#include <typeinfo>
#include <cstddef>
#include <cassert>
#include <iterator>
#include <vector>
#include <ostream>
#include "cpp_utils.h"

typedef unsigned int Dimension;

//! Which position in container an iterator points
enum IteratorPosition { BEGIN, END };

//! Common interface for iterators
template <class T>
class AbstractIterator : public std::iterator<std::forward_iterator_tag, T>
{
public:
    virtual ~AbstractIterator() { }
    virtual AbstractIterator& operator++() = 0;
    virtual void operator++(int) = 0;
    virtual bool operator==(const AbstractIterator& other) const = 0;
    virtual bool operator!=(const AbstractIterator& other) const {
        return !(*this == other);
    }
    virtual T& operator*() const = 0;
    virtual T* operator->() const = 0;
    virtual AbstractIterator& clone() = 0;
};

//! Common interface for const iterators
template <class T>
class AbstractConstIterator
    : public std::iterator<std::input_iterator_tag, const T>
{
public:
    virtual ~AbstractConstIterator() { }
    virtual AbstractConstIterator& operator++() = 0;
    virtual void operator++(int) = 0;
    virtual bool operator==(const AbstractConstIterator& other) const = 0;
    virtual bool operator!=(const AbstractConstIterator& other) const {
        return !(*this == other);
    }
    virtual const T& operator*() const = 0;
    virtual const T* operator->() const = 0;
    virtual AbstractConstIterator& clone() = 0;
};


/** \brief Iterator that uses iterator given to its constructor as implementation
 *
 * This class can be used to "convert" any iterator to subclass
 * of AbstractIterator. Template parameter Impl is type of
 * implementation iterator.
 */
template <class T, class Impl>
class GeneralIteratorWrapper : public AbstractIterator<T>
{
public:
    GeneralIteratorWrapper(Impl implementationBegin, const Impl implementationEnd,
                           IteratorPosition startPosition)
        : m_impl(startPosition == BEGIN ? implementationBegin : implementationEnd),
          m_implEnd(implementationEnd) {
    }
    GeneralIteratorWrapper(const GeneralIteratorWrapper& other)
        : m_impl(other.m_impl), m_implEnd(other.m_implEnd) {
    }
    virtual GeneralIteratorWrapper& operator++() {
        ++m_impl;
        return *this;
    }
    virtual ~GeneralIteratorWrapper() { }
    virtual void operator++(int) {
        ++(*this);
    }
    virtual bool operator==(const AbstractIterator<T>& other) const {
        try {
            return *this ==
                   dynamic_cast<const GeneralIteratorWrapper<T, Impl>&>(other);
        } catch (std::bad_cast& e) {
            return false;
        }
    }
    virtual bool operator==(const GeneralIteratorWrapper<T, Impl>& other) const {
        return m_impl == other.m_impl;
    }
    virtual T& operator*() const {
        return *m_impl;
    }
    virtual T* operator->() const {
        return m_impl;
    }
    virtual GeneralIteratorWrapper& clone() {
        return *new GeneralIteratorWrapper(*this);
    }
private:
    Impl m_impl;
    const Impl m_implEnd;
};

//! Like GeneralIteratorWrapper, but for const iterators
template <class T, class Impl>
class GeneralConstIteratorWrapper : public AbstractConstIterator<T>
{
public:
    GeneralConstIteratorWrapper(Impl implementationBegin,
                                const Impl implementationEnd,
                                IteratorPosition startPosition)
        : m_impl(startPosition == BEGIN ? implementationBegin : implementationEnd),
          m_implEnd(implementationEnd) { }
    GeneralConstIteratorWrapper(const GeneralConstIteratorWrapper& other)
        : m_impl(other.m_impl), m_implEnd(other.m_implEnd) { }
    virtual GeneralConstIteratorWrapper& operator++() {
        ++m_impl;
        return *this;
    }
    virtual ~GeneralConstIteratorWrapper() { }

    virtual void operator++(int) {
        ++(*this);
    }
    virtual bool operator==(const AbstractConstIterator<T>& other) const {
        try {
            return *this ==
                   dynamic_cast<const GeneralConstIteratorWrapper<T, Impl>&>(other);
        } catch (std::bad_cast) {
            return false;
        }
    }
    virtual bool operator==
    (const GeneralConstIteratorWrapper<T, Impl>& other) const {
        return m_impl == other.m_impl;
    }
    virtual const T& operator*() const {
        return *m_impl;
    }
    virtual const T* operator->() const {
        return m_impl;
    }
    virtual GeneralConstIteratorWrapper& clone() {
        return *new GeneralConstIteratorWrapper(*this);
    }
private:
    Impl m_impl;
    const Impl m_implEnd;
};

/** \brief This class is used as iterator for container given to its constructor
 *
 * This class can be used to create an iterator for container
 * that doesn't have one. Iterator uses operator[] to get values
 * from container.
 */
template <class T, class Cont>
class ContainerWrapperIterator : public AbstractIterator<T>
{
public:
    ContainerWrapperIterator(Cont& implContainer, IteratorPosition startPosition)
        : m_position(startPosition == BEGIN ? 0 : implContainer.size()),
          m_container(implContainer) { }
    ContainerWrapperIterator(const ContainerWrapperIterator& other)
        : m_position(other.m_position), m_container(other.m_container) { }
    ContainerWrapperIterator& operator++() {
        ++m_position;
        return *this;
    }
    virtual ~ContainerWrapperIterator() { }
    void operator++(int) {
        ++(*this);
    }
    bool operator==(const AbstractIterator<T>& other) const {
        try {
            return *this == dynamic_cast<ContainerWrapperIterator<T, Cont>&>(other);
        } catch (std::bad_cast) {
            return false;
        }
    }
    bool operator==(const ContainerWrapperIterator<T, Cont>&(other)) {
        return this.m_position == other.m_position;
    }
    const T& operator*() const {
        return m_container[m_position];
    }
    const T* operator->() const {
        return &m_container[m_position];
    }
    virtual ContainerWrapperIterator& clone() {
        return *new ContainerWrapperIterator(this);
    }
private:
    typename Cont::size_type m_position;
    Cont& m_container;
};


//! Like ContainerWrapperIterator, but this is const iterator
template <class T, class Cont>
class ContainerWrapperConstIterator : public AbstractConstIterator<T>
{
public:
    ContainerWrapperConstIterator(const Cont& implContainer,
                                  IteratorPosition startPosition)
        : m_position(startPosition == BEGIN ? 0 : implContainer.size()),
          m_container(implContainer) { }
    ContainerWrapperConstIterator(const ContainerWrapperConstIterator& other)
        : m_position(other.m_position), m_container(other.m_container) { }
    ContainerWrapperConstIterator& operator++() {
        ++m_position;
        return *this;
    }
    virtual ~ContainerWrapperConstIterator() { }
    void operator++(int) {
        ++(*this);
    }
    bool operator==(const AbstractConstIterator<T>& other) const {
        try {
            return *this ==
                   dynamic_cast<const ContainerWrapperConstIterator<T, Cont>&>(other);
        } catch (std::bad_cast) {
            return false;
        }
    }
    bool operator==(const ContainerWrapperConstIterator<T, Cont>& other) const {
        return m_position == other.m_position;
    }
    const T& operator*() const {
        return m_container[m_position];
    }
    const T* operator->() const {
        return &m_container[m_position];
    }
    virtual ContainerWrapperConstIterator& clone() {
        return *new ContainerWrapperConstIterator(*this);
    }
private:
    typename Cont::size_type m_position;
    const Cont& m_container;
};

/** \brief This is iterator that uses any iterator (that has to be subclass
 *         of AbstractIterator) to its constructor as implementation
 *
 * This class is not abstract, so it can be used to "convert" subclass
 * of AbstractIterator of which type is not known to known, "concrete"
 * type. So with this they can be used as values in addition to pointers and
 * references.
 */
template <class T>
class GeneralIterator
    : public std::iterator<std::forward_iterator_tag, T>
{
public:
    GeneralIterator(const GeneralIterator& other) : m_impl(other.m_impl) { }
    template <class IterImpl>
    explicit GeneralIterator(IterImpl implementation)
        : m_impl(implementation) { }
    GeneralIterator& operator++() {
        ++(*m_impl);
        return *this;
    }
    GeneralIterator operator++(int) {
        GeneralIterator prev = *this;
        (*m_impl)++;
        return prev;
    }
    bool operator==(const GeneralIterator& other) const {
        return *m_impl == *(other.m_impl);
    }
    bool operator!=(const GeneralIterator& other) const {
        return *m_impl != *(other.m_impl);
    }
    T& operator*() const {
        return *(*m_impl);
    }
    T* operator->() const {
        return &*(*m_impl);
    }
private:
    ClonePtr<AbstractIterator<T> > m_impl;
};

//! Like GeneralIterator, but this is const iterator
template <class T>
class GeneralConstIterator
    : public std::iterator<std::input_iterator_tag, const T>
{
public:
    GeneralConstIterator(const GeneralConstIterator& other) : m_impl(other.m_impl) { }
    template <class IterImpl>
    explicit GeneralConstIterator(IterImpl implementation)
        : m_impl(implementation) { }
    GeneralConstIterator& operator++() {
        ++(*m_impl);
        return *this;
    }
    GeneralConstIterator operator++(int) {
        GeneralConstIterator prev = *this;
        (*m_impl)++;
        return prev;
    }
    bool operator==(const GeneralConstIterator& other) const {
        return *m_impl == *(other.m_impl);
    }
    bool operator!=(const GeneralConstIterator& other) const {
        return *m_impl != *(other.m_impl);
    }
    const T& operator*() const {
        return *(*m_impl);
    }
    const T* operator->() const {
        return &*(*m_impl);
    }
private:
    ClonePtr<AbstractConstIterator<T> > m_impl;
};


//! Common interface for all sequence containers that have (at least) const operations.
template <class T>
class ConstAbstractSequence
{
public:
    typedef std::size_t size_type;
    typedef GeneralIterator<T> iterator;
    typedef GeneralConstIterator<T> const_iterator;
    virtual ~ConstAbstractSequence() { }
    virtual const T& operator[] (size_type idx) const = 0;
    virtual size_type size() const = 0;
    virtual const_iterator begin() const = 0;
    virtual const_iterator end() const = 0;
};

//! Common interface for all sequence containers. Unlike ConstAbstractSequence, this has both const and non-const operations.
template <class T>
class AbstractSequence : public ConstAbstractSequence<T>
{
public:
    typedef typename ConstAbstractSequence<T>::size_type size_type;
    typedef typename ConstAbstractSequence<T>::iterator iterator;
    typedef typename ConstAbstractSequence<T>::const_iterator const_iterator;
    virtual ~AbstractSequence() { }
    using ConstAbstractSequence<T>::begin;
    using ConstAbstractSequence<T>::end;
    using ConstAbstractSequence<T>::operator[];
    virtual T& operator[] (size_type idx) = 0;
    virtual iterator begin() = 0;
    virtual iterator end() = 0;
};

/** \brief Sequence container that gets values of elements from function that is given to constructor
 *
 * Non-const member functions are commented out, because
 * function container is not mutable.
 */
template <class T, class Func>
class FunctionSequence : public ConstAbstractSequence<T>
{
public:
    typedef typename ConstAbstractSequence<T>::size_type size_type;
    typedef typename ConstAbstractSequence<T>::iterator iterator;
    typedef typename ConstAbstractSequence<T>::const_iterator const_iterator;
    FunctionSequence(Func function, size_type size)
        : m_function(function), m_size(size), m_hasStoredValue(false),
          m_storedValueIdx(0) { }
    virtual ~FunctionSequence() { }
    // virtual T operator[](size_type idx) { return m_function(idx); }
    virtual const T& operator[](size_type idx) const {
        if (m_storedValueIdx != idx || (!m_hasStoredValue)) {
            tmp = m_function(idx);
            m_hasStoredValue = true;
            m_storedValueIdx = idx;
        }
        return tmp;
    }
    virtual size_type size() const {
        return m_size;
    }
    /*
    virtual iterator begin()
    { return GeneralIterator<T>
        (ContainerWrapperIterator<T, FunctionSequence> (*this, BEGIN)); }
    virtual iterator end()
    { return GeneralIterator<T>
        (ContainerWrapperIterator<T, FunctionSequence> (*this, END)); }
    */
    virtual const_iterator begin() const {
        return GeneralConstIterator<T>
               (ContainerWrapperConstIterator<T, FunctionSequence> (*this, BEGIN));
    }
    virtual const_iterator end() const {
        return GeneralConstIterator<T>
               (ContainerWrapperConstIterator<T, FunctionSequence> (*this, END));
    }
protected:
    Func m_function;
    size_type m_size;
    mutable bool m_hasStoredValue;
    mutable size_type m_storedValueIdx;
    mutable T tmp;
};

/** \brief Like FunctionSequence, but caches values of elements, so that
 *         function is not called multiple times if value of same element is
 *         accessed multiple times.
 */
template <class T, class Func>
class CachedFunctionSequence : public FunctionSequence<T, Func>
{
public:
    typedef typename FunctionSequence<T, Func>::size_type size_type;
    CachedFunctionSequence(Func function, size_type size)
        : FunctionSequence<T, Func>(function, size), cachedValues(size, 0) { }
    virtual ~CachedFunctionSequence() {
        for (size_type n = 0; n < cachedValues.size(); ++n)
            delete cachedValues[n];
    }
    virtual const T& operator[](size_type idx) const {
        if (cachedValues[idx] == 0)
            cachedValues[idx] = new T(m_function(idx));
        return *cachedValues[idx];
    }
protected:
    using FunctionSequence<T, Func>::m_function;
    mutable std::vector<T*> cachedValues;
};

//! Sequence container that is implemented with array.
template <class T>
class ArraySequence : public AbstractSequence<T>
{
public:
    typedef typename AbstractSequence<T>::size_type size_type;
    typedef typename AbstractSequence<T>::iterator iterator;
    typedef typename AbstractSequence<T>::const_iterator const_iterator;
    ArraySequence(const ArraySequence& other)
        : m_size(other.m_size) {
        m_data = makeArray(other);
    }
    ArraySequence(T* values, std::size_t size)
        : m_data(values), m_size(size) { }
    explicit ArraySequence(const std::vector<T>& values)
        : m_data(new T[values.size()]), m_size(values.size()) {
        for (std::size_t idx = 0; idx < values.size(); ++idx)
            m_data[idx] = values[idx];
    }
    ArraySequence(const T& value0, const T& value1, const T& value2)
        : m_data(new T[3]), m_size(3) {
        m_data[0] = value0;
        m_data[1] = value1;
        m_data[2] = value2;
    }
    virtual ~ArraySequence() {
        ownDeleteData();
    }
    ArraySequence& operator=(const ArraySequence& other) {
        if (this != &other) {
            deleteData();
            m_data = makeArray(other);
        }
        return *this;
    }
    virtual T& operator[] (size_type idx) {
        return m_data[idx];
    }
    virtual const T& operator[] (size_type idx) const {
        return m_data[idx];
    }
    virtual size_type size() const {
        return m_size;
    }
    virtual iterator begin() {
        return GeneralIterator<T>(GeneralIteratorWrapper<T, T*>
                                  (&m_data[0], &m_data[m_size], BEGIN));
    }
    virtual iterator end() {
        return GeneralIterator<T>(GeneralIteratorWrapper<T, T*>
                                  (&m_data[0], &m_data[m_size], END));
    }
    virtual const_iterator begin() const {
        return GeneralConstIterator<T>(GeneralConstIteratorWrapper<T, const T*>
                                       (&m_data[0], &m_data[m_size], BEGIN));
    }
    virtual const_iterator end() const {
        return GeneralConstIterator<T>(GeneralConstIteratorWrapper<T, const T*>
                                       (&m_data[0], &m_data[m_size], END));
    }
protected:
    T* m_data;
    const size_t m_size;
    T* makeArray(ArraySequence<T> seq) {
        T* newAr = new T[seq.size()];
        for (std::size_t idx = 0; idx < seq.size(); ++idx)
            newAr[idx] = seq[idx];
        return newAr;
    }
    virtual void deleteData() {
        ownDeleteData();
    }
private:
    void ownDeleteData() {
        delete[] m_data;
    }
};

/** \brief Like ArraySequence, but when destructor is called,
 *         calls delete[] statement on all elements of sequence.
 *
 * Can be used to automatically delete arrays in sequence.
 */
template <class T>
class AutoADelArraySequence : public ArraySequence<T>
{
public:
    AutoADelArraySequence(const AutoADelArraySequence& other)
        : ArraySequence<T>(other) { }
    AutoADelArraySequence(T* values, std::size_t size)
        : ArraySequence<T>(values, size) { }
    explicit AutoADelArraySequence(const std::vector<T>& values)
        : ArraySequence<T>(values) { }
    virtual ~AutoADelArraySequence() {
        ownDeleteData();
    }
protected:
    using ArraySequence<T>::m_size;
    using ArraySequence<T>::m_data;
    virtual void deleteData() {
        ownDeleteData();
        ArraySequence<T>::deleteData();
    }
private:
    void ownDeleteData() {
        for (size_t idx = 0; idx < m_size; ++idx)
            delete[] m_data[idx];
    }
};

/** \brief Generic class to implement handle classes by subclassing.
 *
 * Constructor takes pointer implementation as argument, and this
 * implementation automatically deleted when destructor of handle
 * is called. Handle follows reference semantics: when handle
 * is copied, implementation is not copied, but reference to
 * implementation is shared.
 */
template<class Impl>
class GenericHandle
{
protected:
    GenericHandle() : m_implementation(0) { }
    GenericHandle(Impl* implementation) : m_implementation(implementation) { }
    Impl& getImpl() {
        return *m_implementation;
    }
    const Impl& getImpl() const {
        return *m_implementation;
    }
private:
    SharedPtr<Impl> m_implementation;
};

/** \brief Handle to sequence container that has const operations.
 *
 * Follows reference semantics like GenericHandle. Forwards
 * operations to implementation, that constuctor takes as argument.
 */
template <class T, class Impl = ConstAbstractSequence<T> >
class ConstSequenceHandle : public GenericHandle<Impl>
{
public:
    typedef typename Impl::size_type size_type;
    typedef typename Impl::iterator iterator;
    typedef typename Impl::const_iterator const_iterator;
    using GenericHandle<Impl>::getImpl;
    ConstSequenceHandle() : GenericHandle<Impl>() { }
    ConstSequenceHandle(Impl* implementation)
        : GenericHandle<Impl>(implementation) { }
    const T& operator[] (size_type idx) const {
        return getImpl()[idx];
    }
    size_type size() const {
        return getImpl().size();
    }
    const_iterator begin() const {
        return const_iterator(getImpl().begin());
    }
    const_iterator end() const {
        return const_iterator(getImpl().end());
    }
};

//! Like ConstSequenceHandle, but has also non-const operations
template <class T, class Impl = AbstractSequence<T> >
class SequenceHandle : public ConstSequenceHandle<T, Impl>
{
public:
    typedef typename Impl::size_type size_type;
    typedef typename Impl::iterator iterator;
    typedef typename Impl::const_iterator const_iterator;
    using ConstSequenceHandle<T, Impl>::getImpl;
    using ConstSequenceHandle<T, Impl>::operator[];
    using ConstSequenceHandle<T, Impl>::begin;
    using ConstSequenceHandle<T, Impl>::end;
    SequenceHandle() : ConstSequenceHandle<T, Impl>() { }
    SequenceHandle(Impl* implementation)
        : ConstSequenceHandle<T, Impl>(implementation) { }
    T& operator[] (size_type idx) {
        return getImpl()[idx];
    }
    iterator begin() {
        return iterator(getImpl().begin());
    }
    iterator end() {
        return iterator(getImpl().end());
    }
};

/** \brief Common operations for Coordinate-templates, so that
 *         there's no need to duplicate them in template specializations.
 *         This class is for implementation, not for users.
 */
template<Dimension nDims>
class CoordinatesCommon
{
public:
    CoordinatesCommon() { }
    CoordinatesCommon(const CoordinatesCommon<nDims>& coords) {
        for (unsigned int i = 0; i < nDims; ++i)
            m_coords[i] = coords[i];
    }
    std::size_t& operator[] (Dimension coord) {
        return m_coords[coord];
    }
    const std::size_t& operator[] (Dimension coord) const {
        return m_coords[coord];
    }
protected:
    std::size_t m_coords[nDims];
};

//! Represents coordinates with nDims dimensions
template<Dimension nDims>
class Coordinates : public CoordinatesCommon<nDims>
{
public:
    Coordinates() : CoordinatesCommon<nDims>() { }
    Coordinates(const Coordinates<nDims>& other)
        : CoordinatesCommon<nDims>(other) { }
};

//! Represents coordinates with 3 dimensions
template<>
class Coordinates<3> : public CoordinatesCommon<3>
{
public:
    Coordinates() : CoordinatesCommon<3>() { }
    Coordinates(const Coordinates<3>& other) : CoordinatesCommon<3>(other) { }
    Coordinates(std::size_t x, std::size_t y, std::size_t z) {
        m_coords[0] = x;
        m_coords[1] = y;
        m_coords[2] = z;
    }
};

template <Dimension nDims>
bool operator==(const Coordinates<nDims>& coords1,
                const Coordinates<nDims>& coords2)
{
    for (Dimension dim = 0; dim < nDims; ++dim)
        if (coords1[dim] != coords2[dim])
            return false;
    return true;
}

template <Dimension nDims>
bool operator!=(const Coordinates<nDims>& coords1,
                const Coordinates<nDims>& coords2)
{
    return !(coords1 == coords2);
}

template <Dimension nDims>
bool operator<(const Coordinates<nDims>& coords1,
               const Coordinates<nDims>& coords2)
{
    for (Dimension dim = nDims - 1; dim < nDims && dim >= 0; --dim) {
        if (coords1[dim] < coords2[dim])
            return true;
        else if (coords1[dim] > coords2[dim])
            return false;
    }
    return false;
}

template <Dimension nDims>
bool operator>(const Coordinates<nDims>& coords1,
               const Coordinates<nDims>& coords2)
{
    return coords2 < coords1;
}

template <Dimension nDims>
bool operator<=(const Coordinates<nDims>& coords1,
                const Coordinates<nDims>& coords2)
{
    for (Dimension dim = nDims - 1; dim < nDims && dim >= 0; --dim) {
        if (coords1[dim] < coords2[dim])
            return true;
        else if (coords1[dim] > coords2[dim])
            return false;
    }
    return true;
}

template <Dimension nDims>
bool operator>=(const Coordinates<nDims>& coords1,
                const Coordinates<nDims>& coords2)
{
    return coords2 <= coords1;
}

template <Dimension nDims>
Coordinates<nDims>& operator+=(Coordinates<nDims>& coords1,
                               const Coordinates<nDims>& coords2)
{
    for (Dimension dim = 0; dim < nDims; ++dim)
        coords1[dim] += coords2[dim];
    return coords1;
}

template <Dimension nDims>
Coordinates<nDims>& operator-=(Coordinates<nDims>& coords1,
                               const Coordinates<nDims>& coords2)
{
    for (Dimension dim = 0; dim < nDims; ++dim)
        coords1[dim] -= coords2[dim];
    return coords1;
}

template <Dimension nDims>
const Coordinates<nDims> operator+(const Coordinates<nDims>& coords1,
                                   const Coordinates<nDims>& coords2)
{
    Coordinates<nDims> sum = coords1;
    sum += coords2;
    return sum;
}

template <Dimension nDims>
const Coordinates<nDims> operator-(const Coordinates<nDims>& coords1,
                                   const Coordinates<nDims>& coords2)
{
    Coordinates<nDims> diff = coords1;
    diff -= coords2;
    return diff;
}

template <Dimension nDims, class ScalarType>
Coordinates<nDims>& operator*=(Coordinates<nDims>& coords, ScalarType mult)
{
    for (Dimension dim = 0; dim < nDims; ++dim)
        coords[dim] *= mult;
    return coords;
}

template <Dimension nDims, class ScalarType>
Coordinates<nDims>& operator/=(Coordinates<nDims>& coords, ScalarType div)
{
    for (Dimension dim = 0; dim < nDims; ++dim)
        coords[dim] /= div;
    return coords;
}

template <Dimension nDims, class ScalarType>
const Coordinates<nDims> operator*(Coordinates<nDims> coords,
                                   ScalarType mult)
{
    Coordinates<nDims> product = coords;
    product *= mult;
    return product;
}

template <Dimension nDims, class ScalarType>
const Coordinates<nDims> operator/(Coordinates<nDims> coords,
                                   ScalarType div)
{
    Coordinates<nDims> quotient = coords;
    quotient /= div;
    return quotient;
}

//! Prints Coordinates to stream in textual form
template <Dimension nDims>
std::ostream& operator<<(std::ostream& output, const Coordinates<nDims>& coords)
{
    output << "(";
    if (nDims > 0) {
        for (Dimension dim = 0; dim < nDims - 1; ++dim)
            output << coords[dim] << ", ";
        output << coords[nDims - 1] << ")";
    }
    return output;
}


/** \brief Represents interval with nDims dimensions.
 *
 * Begin coordinates of interval are inclusive, end coordinates are exclusive.
 */
template<Dimension nDims>
class Interval
{
public:
    class const_iterator
        : public std::iterator<std::input_iterator_tag, const Coordinates<nDims> >
    {
    public:
        const_iterator(const Interval& i, IteratorPosition startPosition)
            : m_currentPosition(i.getBegin()), m_atEnd(startPosition == END),
              m_interval(i) { }
        const_iterator& operator++() {
            if (!m_atEnd)
                m_atEnd = !incrementCoordinatesInInterval
                          (m_currentPosition, m_interval);
            return *this;
        }
        const_iterator& operator++(int) {
            iterator prev = *this;
            ++(*this);
            return prev;
        }
        bool operator==(const const_iterator& other) const {
            return m_atEnd == other.m_atEnd
                   && m_currentPosition == other.m_currentPosition;
        }
        bool operator!=(const const_iterator& other) const {
            return !((*this) == other);
        }
        const Coordinates<nDims>& operator*() const {
            return m_currentPosition;
        }
        const Coordinates<nDims>* operator->() const {
            return &m_currentPosition;
        }
    protected:
        Coordinates<nDims> m_currentPosition;
        bool m_atEnd;
        const Interval& m_interval;
    };
    class iterator
        : public std::iterator<std::forward_iterator_tag, Coordinates<nDims> >,
      public const_iterator
    {
    public:
        iterator(const Interval& i, IteratorPosition startPosition)
            : const_iterator(i, startPosition) { }
        Coordinates<nDims> operator*() const {
            return this->m_currentPosition;
        }
        Coordinates<nDims>* operator->() const {
            return &(this->m_currentPosition);
        }
    };
    Interval() { }
    Interval(const Interval<nDims>& other)
        : m_begin(other.m_begin), m_end(other.m_end) { }
    Interval(const Coordinates<nDims>& begin, const Coordinates<nDims>& end)
        : m_begin(begin), m_end(end) { }
    Coordinates<nDims>& getBegin() {
        return m_begin;
    }
    const Coordinates<nDims>& getBegin() const {
        return m_begin;
    }
    Coordinates<nDims>& getEnd() {
        return m_end;
    }
    const Coordinates<nDims>& getEnd() const {
        return m_end;
    }
    iterator begin() {
        return iterator(*this, BEGIN);
    }
    const_iterator begin() const {
        return const_iterator(*this, BEGIN);
    }
    iterator end() {
        return iterator(*this, END);
    }
    const_iterator end() const {
        return const_iterator(*this, END);
    }
private:
    Coordinates<nDims> m_begin;
    Coordinates<nDims> m_end;
};

//! Returns volume of interval
template<Dimension nDims>
std::size_t volume(const Interval<nDims>& interval)
{
    if (nDims == 0)
        return 0;
    std::size_t size = 1;
    for (Dimension dim = 0; dim < nDims; ++dim)
        size *= interval.getEnd()[dim] - interval.getBegin()[dim];
    return size;
}

template <Dimension nDims, class ScalarType>
Interval<nDims> operator*(Interval<nDims> interval, ScalarType factor)
{
    return Interval<nDims>(interval.getBegin() * factor,
                           interval.getEnd() * factor);
}

//! Prints Interval to stream in textual form
template <Dimension nDims>
std::ostream& operator<<(std::ostream& output, const Interval<nDims>& i)
{
    output << "(" << i.getBegin() << ", " << i.getEnd() << ")";
    return output;
}


//! Return true if "interval1" is subset of "interval2", otherwise returns false.
template <Dimension nDims>
bool isSubSet(const Interval<nDims>& interval1,
              const Interval<nDims>& interval2)
{
    for (Dimension dim = 0; dim < nDims; ++dim)
        if (interval1.getBegin()[dim] < interval2.getBegin()[dim] ||
            interval1.getEnd()[dim] > interval2.getEnd()[dim])
            return false;
    return true;
}


/** \brief Returns true if "interval1" (that is given to constructor)
 *         intersects with "interval2" (that is given to operator()),
 *         otherwise returns false.
 */
template <Dimension nDims>
struct IntervalsIntersect {
    IntervalsIntersect(const Interval<nDims>& interval1)
        : m_interval1(interval1) { }
    bool operator()(const Interval<nDims>& interval2) {
        for (Dimension dim = 0; dim < nDims; ++dim)
            if (m_interval1.getEnd()[dim] <= interval2.getBegin()[dim] ||
                m_interval1.getBegin()[dim] >= interval2.getEnd()[dim])
                return false;
        return true;
    }
private:
    const Interval<nDims>& m_interval1;
};


/** \brief Return true if coordinates "coords" are in interval "interval",
 *         otherwise returns false.
 */
template <Dimension nDims>
bool coordinatesInInterval(const Coordinates<nDims>& coords,
                           const Interval<nDims>& interval)
{
    for (Dimension dim = 0; dim < nDims; ++dim)
        if (coords[dim] < interval.getBegin()[dim] ||
            coords[dim] >= interval.getEnd()[dim])
            return false;
    return true;
}

/** \brief Increments coordinates "coords" by one in interval "interval".
 *
 * Goes through interval so that first ("x") component of coordinates
 * changes most quickly, second component secondly quickly and so on.
 * Return true if incrementing succeeded (coordinates were not at end
 * of interval), otherwise returns false.
 */
template <Dimension nDims>
bool incrementCoordinatesInInterval(Coordinates<nDims>& coords,
                                    const Interval<nDims>& interval)
{
    Dimension dim;
    for (dim = 0; dim < nDims && coords[dim] >= interval.getEnd()[dim] - 1;
         ++dim)
        coords[dim] = interval.getBegin()[dim];
    if (dim == nDims)
        return false;
    else {
        ++coords[dim];
        return true;
    }
}


//! Retuns true if given coordinates are in some of given intervals, otherwise returns false.
template <Dimension nDims>
bool coordinatesInInterval(const Coordinates<nDims>& coords,
                           const std::vector<Interval<nDims> >& intervals)
{
    for (typename std::vector<Interval<nDims> >::const_iterator interval
         = intervals.begin(); interval != intervals.end(); ++interval)
        if (coordinatesInInterval(coords, *interval))
            return true;
    return false;
}

#endif

