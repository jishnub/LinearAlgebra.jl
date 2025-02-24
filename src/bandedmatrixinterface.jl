module BandedMatrixInterface
# The following code is borrowed from BandedMatrices.jl, which is distributed under the MIT license
#=
Copyright (c) 2016 ApproxFun

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
=#

export colsupport, rowsupport

public bandwidths, rowrange, colrange

#=
We define and use `rowsupport` and `colsupport` in general.
For `BandedMatrix`es, we define `bandwidths`. This automatically defines
`bandwidths` -> `colrange -> colsupport` and`bandwidths -> rowrange` -> `rowsupport`.
For these matrices, `colsupport` and `colrange` are equivalent, as the support is
literally a range. Similarly, for `rowsupport` and `rowrange`.

In general, the loops should be over `rowsupport` and `colsupport`, and custom matrix
types may extend these functions.
=#

"""
    bandwidths(A)

Returns a tuple containing the lower and upper bandwidth of `A`, in order.

!!! note
    The lower bandwidth represents the subdiagonal index, so it would
    be positive for a matrix containing non-zero subdiagonals.
"""
bandwidths(A::AbstractMatrix) = (size(A,1)-1 , size(A,2)-1)

"""
    bandwidth(A,i)

Returns the lower bandwidth (`i==1`) or the upper bandwidth (`i==2`).
"""
bandwidth(A, k::Integer) = bandwidths(A)[k]

"""
    colsupport(A)

Return an iterator containing the row indices of the possible non-zero entries in `A`.
"""
colsupport(A) = colsupport(A, axes(A,2))
"""
    colsupport(A, j)

Return an iterator containing the row indices of the possible non-zero entries in the `j`-th column of `A`.
"""
colsupport(A, j) = axes(A,1)

"""
    colstart(A, i::Integer)

Return the starting row index of the filled bands in the `i`-th column,
bounded by the actual matrix size.
"""
colstart(A, i::Integer) = max(i-bandwidth(A,2), 1) + firstindex(A,1)-1

"""
    colstop(A, i::Integer)

Return the stopping row index of the filled bands in the `i`-th column,
bounded by the actual matrix size.
"""
colstop(A, i::Integer) = clamp(i+bandwidth(A,1), 0:size(A, 1)) + firstindex(A,1)-1

"""
    rowsupport(A, k)

Return an iterator containing the column indices of the possible non-zero entries in the `k`-th row of `A`.
"""
rowsupport(A, k) = axes(A,2)
"""
    rowsupport(A)

Return an iterator containing the column indices of the possible non-zero entries in `A`.
"""
rowsupport(A) = rowsupport(A, axes(A,1))

"""
    rowstart(A, i::Integer)

Return the starting column index of the filled bands in the `i`-th row,
bounded by the actual matrix size.
"""
rowstart(A, i::Integer) = max(i-bandwidth(A,1), 1) + firstindex(A,2)-1

"""
    colrange(A, i::Integer)

Return the range of rows in the `i`-th column that correspond to filled bands.
"""
colrange(A, i::Integer) = colstart(A,i):colstop(A,i)

"""
    rowrange(A, i::Integer)

Return the range of columns in the `i`-th row that correspond to filled bands.
"""
rowrange(A, i::Integer) = rowstart(A,i):rowstop(A,i)

"""
    rowstop(A, i::Integer)

Return the stopping column index of the filled bands in the `i`-th row,
bounded by the actual matrix size.
"""
rowstop(A, i::Integer) = clamp(i+bandwidth(A,2), 0:size(A, 2)) + firstindex(A,2)-1

end
