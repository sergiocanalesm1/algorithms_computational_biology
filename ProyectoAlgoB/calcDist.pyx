from libc.math cimport pow
from libc.math cimport sqrt
cpdef double calcDist_cython (double x1, double x2, double x3, double y1, double y2, double y3):
    """Compute x^2 + x as double.

    This is a cdef function that can be called from within
    a Cython program, but not from Python.
    """
    return sqrt(pow((x1-y1), 2.0) + pow((x2-y2), 2.0) + pow((x3-y3), 2.0))

cpdef double calcLJ_cython (double eps, double sigma, double r):
    return 4*eps*( pow( (sigma/r) , 12) - pow((sigma/r) ,6) )
#cpdef print_result (double x1, double x2, double x3, double y1, double y2, double y3):
#    """This is a cpdef function that can be called from Python."""
#    print("{}".format(square_and_add( double x1, double x2, double x3, double y1, double y2, double y3)))
