from functools import cached_property
from mathematics.vec3 import Vec3
from typing import Self

ComplexNumber = int | float | complex
TupleElsMat33 = tuple[ComplexNumber, ComplexNumber, ComplexNumber, ComplexNumber, ComplexNumber, ComplexNumber, ComplexNumber, ComplexNumber, ComplexNumber]

class Mat33:
    els: TupleElsMat33

    def __init__(self, *els: ComplexNumber):
        """
        [[
            a, b, c
            d, e, f,
            g, h, i
        ]]
        """
        self.els = tuple(els)

    @staticmethod
    def from_col_vectors(c1: Vec3, c2: Vec3, c3: Vec3) -> "Mat33":
        """Vectors should be dimensionles"""
        return Mat33(
            c1._x, c2._x, c3._x,
            c1._y, c2._y, c3._y,
            c1._z, c2._z, c3._z,
        )
    
    def __getitem__(self, indexes: tuple[int, int]) -> ComplexNumber:
        i, j = indexes
        assert 0 <= i <= 2
        assert 0 <= j <= 2

        return self.els[3*i + j]
    
    def prod_scalar(self, scalar: ComplexNumber) -> Self:
        return Mat33([scalar*el for el in self.els])
    
    def prod_mat(self, mat_right: Self) -> Self:
        a,  b,  c,  d,  e,  f,  g,  h,  i  = self.els
        ra, rb, rc, rd, re, rf, rg, rh, ri = mat_right.els

        return Mat33(
            a*ra + b*rd + c*rg, a*rb + b*re + c*rh, a*rc + b*rf + c*ri,
            d*ra + e*rd + f*rg, d*rb + e*re + f*rh, d*rc + e*rf + f*ri,
            g*ra + h*rd + i*rg, g*rb + h*re + i*rh, g*rc + h*rf + i*ri,
        )

    def prod_vec(self, vec: Vec3) -> Vec3:
        a, b, c, d, e, f, g, h, i = self.els
        x, y, z = vec.xyz
        return Vec3(
            a*x + b*y + c*z,
            d*x + e*y + f*z,
            g*x + h*y + i*z,
            vec.unit
        )
    
    def __mul__(self, scalar: ComplexNumber) -> Self:
        if not isinstance(scalar, ComplexNumber):
            raise TypeError(f"Unsupported operation with {scalar} of type `{type(scalar)}`")
        
        return self.prod_scalar(scalar)
    
    def __matmul__(self, o: Self | Vec3 | ComplexNumber) -> Self:
        if isinstance(o, Mat33):
            return self.prod_mat(o)
        elif isinstance(o, Vec3):
            return self.prod_vec(o)
        
        raise TypeError(f"Unsupported operation with {o} of type `{type(o)}`")
    
    def __add__(self, mat_right: Self) -> "Mat33":
        a,  b,  c,  d,  e,  f,  g,  h,  i  = self.els
        ra, rb, rc, rd, re, rf, rg, rh, ri = mat_right.els

        return Mat33(
            a+ra, b+rb, c+rc,
            d+rd, e+re, f+rf,
            g+rg, h+rh, i+ri
        )
    
    def __sub__(self, mat_right: Self) -> Self:
        a,  b,  c,  d,  e,  f,  g,  h,  i  = self.els
        ra, rb, rc, rd, re, rf, rg, rh, ri = mat_right.els
        
        return Mat33(
            a-ra, b-rb, c-rc,
            d-rd, e-re, f-rf,
            g-rg, h-rh, i-ri
        )
    
    def __neg__(self) -> Self:
        return self.prod_scalar(-1)
    
    @cached_property
    def inverse(self) -> Self:
        raise NotImplementedError()

MAT_I = Mat33(
    1, 0, 0, 
    0, 1, 0, 
    0, 0, 1
)

MAT_0 = Mat33(*[0]*9)