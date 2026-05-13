from abc import ABC, abstractmethod
from typing import Self, cast

from leorbit.mathematics import U, Tensor, TensorBound, TensorKind, matrix33, cos, scalar, sin



class Transform(ABC):
    """Abstract reversible transform between two tensor-like types."""

    @abstractmethod
    def _do(self, tensor: Tensor) -> Tensor: ...

    def do(self, tensor: Tensor) -> Tensor:
        output = self._do(
            self._input_bound.secure(tensor)
        )

        return self._output_bound.secure(output)

    @abstractmethod
    def _undo(self, tensor: Tensor) -> Tensor: ...

    def undo(self, tensor: Tensor) -> Tensor:
        output = self._undo(
            self._output_bound.secure(tensor)
        )

        return self._input_bound.secure(output)

    @abstractmethod
    def copy(self) -> Self: ...

    @property
    @abstractmethod
    def _input_bound(self) -> TensorBound: ...

    @property
    @abstractmethod
    def _output_bound(self) -> TensorBound: ...

    def reverse(self) -> Transform:
        """Return a transform with ``do`` and ``undo`` swapped."""
        t = self.copy()

        do = t.do
        undo = t.undo

        tt = cast(Transform, t)

        tt.do = undo  # type: ignore
        tt.undo = do  # type: ignore

        return tt


class TransformIdentity(Transform):
    """Identity transform that leaves values unchanged."""

    @property
    def _input_bound(self) -> TensorBound:
        return TensorBound()

    @property
    def _output_bound(self) -> TensorBound:
        return TensorBound()

    def _do(self, tensor: Tensor) -> Tensor:
        """Return the input value unchanged."""
        return tensor

    def _undo(self, tensor: Tensor) -> Tensor:
        """Return the input value unchanged."""
        return tensor

    def copy(self) -> Self:
        """Return a new identity transform instance."""
        t = TransformIdentity()
        return cast(Self, t)


class TransformVector3Linear(Transform):
    """Linear transform for vectors using a dimensionless 3×3 matrix."""

    def __init__(self, matrix: Tensor):
        """Initialize with the transformation matrix."""
        matrix.secure(units=U.dimensionless, kind=TensorKind.MATRIX33)
        self.matrix = matrix

    @property
    def _input_bound(self) -> TensorBound:
        return TensorBound(kind=TensorKind.VECTOR3)

    @property
    def _output_bound(self) -> TensorBound:
        return TensorBound(kind=TensorKind.VECTOR3)

    def _do(self, tensor: Tensor) -> Tensor:
        """Apply the linear transformation."""
        return self.matrix.matrix33 @ tensor.vector3

    def _undo(self, tensor: Tensor) -> Tensor:
        """Apply the inverse linear transformation."""
        return self.matrix.matrix33.inverse().matrix33 @ tensor.vector3

    def copy(self) -> Self:
        """Return a deep copy of this transform."""
        t = TransformVector3Linear(
            matrix=self.matrix.copy()
        )
        return cast(Self, t)


class TransformVector3Affine(Transform):
    """Affine transform combining linear map and translation."""

    def __init__(self, matrix: Tensor, translation: Tensor):
        """Initialize with matrix and translation components."""
        matrix.secure(units=U.dimensionless, kind=TensorKind.MATRIX33)
        translation.secure(kind=TensorKind.VECTOR3)
        self.matrix = matrix
        self.translation = translation

    @property
    def _input_bound(self) -> TensorBound:
        return TensorBound(kind=TensorKind.VECTOR3)

    @property
    def _output_bound(self) -> TensorBound:
        return TensorBound(kind=TensorKind.VECTOR3)

    def _do(self, tensor: Tensor) -> Tensor:
        """Apply affine transform ``M @ v + t``."""
        return self.matrix.matrix33 @ tensor.vector3 + self.translation

    def _undo(self, tensor: Tensor) -> Tensor:
        """Apply inverse affine transform ``M⁻¹ @ (v - t)``."""
        return self.matrix.matrix33.inverse().matrix33 @ (tensor - self.translation).vector3

    def copy(self) -> Self:
        """Return a deep copy of this transform."""
        t = TransformVector3Affine(
            matrix=self.matrix.copy(),
            translation=self.translation.copy()
        )
        return cast(Self, t)
    
    def as_linear_transform(self) -> TransformVector3Linear:
        """Returns this transform as linear, getting rid of the translation"""
        return TransformVector3Linear(self.matrix)


class TransformVector3RotationZ(TransformVector3Linear):
    """Rotation around the Z axis by a given angle."""

    def __init__(self, angle: Tensor):
        angle.secure(units=U.radian, kind=TensorKind.SCALAR)
        c = cos(angle).scalar.value("dimensionless")
        s = sin(angle).scalar.value("dimensionless")

        super().__init__(
            matrix33(
                c, -s, 0,
                s, c, 0,
                0, 0, 1,
            )
        )

        self.angle = angle

    def copy(self) -> Self:
        t = TransformVector3RotationZ(self.angle.copy())
        return cast(Self, t)


class TransformChain(Transform):
    def __init__(self, *transforms: Transform):
        """Create a composite transform executed in the given order."""
        self.transforms = list(transforms)

    @property
    def _input_bound(self) -> TensorBound:
        if not self.transforms:
            return TensorBound()
        return self.transforms[0]._input_bound

    @property
    def _output_bound(self) -> TensorBound:
        if not self.transforms:
            return TensorBound()
        return self.transforms[-1]._output_bound

    def _do(self, tensor: Tensor) -> Tensor:
        """Apply all transforms in forward order."""
        result = tensor
        for t in self.transforms:
            result = t.do(result)
        return result

    def _undo(self, tensor: Tensor) -> Tensor:
        """Apply all inverse transforms in reverse order."""
        result = tensor
        for t in reversed(self.transforms):
            result = t.undo(result)
        return result

    def copy(self) -> Self:
        """Return a deep copy of the transform chain."""
        t = TransformChain(
            *(t.copy() for t in self.transforms)
        )

        return cast(Self, t)
    
################################

def _secure_position_transform(transform: Transform) -> Transform:
    """Reshapes a transformation so that it fits the requirements to
    be a 'position' transform"""
    
    if isinstance(transform, TransformVector3Linear):
        return transform
    if isinstance(transform, TransformVector3Affine):
        return transform
    elif isinstance(transform, TransformChain):
        return TransformChain(
            *map(_secure_position_transform, transform.transforms)
        )
    
    raise ValueError("Uncompatible transform type for position transformation")

def _secure_velocity_transform(transform: Transform) -> Transform:
    """Reshapes a transformation so that it fits the requirements to
    be a 'velocity' transform"""

    if isinstance(transform, TransformVector3Linear):
        return transform
    if isinstance(transform, TransformVector3Affine):
        return transform.as_linear_transform()
    elif isinstance(transform, TransformChain):
        return TransformChain(
            *map(_secure_velocity_transform, transform.transforms)
        )
    
    raise ValueError("Uncompatible transform type for velocity transformation")