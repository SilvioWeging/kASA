/***************************************************************************
 *  include/stxxl/bits/common/is_sorted.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_IS_SORTED_HEADER
#define STXXL_COMMON_IS_SORTED_HEADER

#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

template <class ForwardIterator>
bool is_sorted_helper(ForwardIterator first, ForwardIterator last)
{
    if (first == last)
        return true;

    ForwardIterator next = first;
    for (++next; next != last; first = next, ++next) {
        if (*next < *first)
            return false;
    }

    return true;
}

template <class ForwardIterator, class StrictWeakOrdering>
bool is_sorted_helper(ForwardIterator first, ForwardIterator last,
                      StrictWeakOrdering comp)
{
    if (first == last)
        return true;

    ForwardIterator next = first;
    for (++next; next != last; first = next, ++next) {
        if (comp(*next, *first))
            return false;
    }

    return true;
}

template <class ForwardIterator>
bool is_sorted(ForwardIterator first, ForwardIterator last)
{
    return is_sorted_helper(first, last);
}

template <class ForwardIterator, class StrictWeakOrdering>
bool is_sorted(ForwardIterator first, ForwardIterator last,
               StrictWeakOrdering comp)
{
    return is_sorted_helper(first, last, comp);
}

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_IS_SORTED_HEADER
