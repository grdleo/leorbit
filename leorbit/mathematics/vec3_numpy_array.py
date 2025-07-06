from dataclasses import dataclass
from functools import cached_property

from numpy.typing import NDArray


@dataclass
class Vec3NumpyArray:
	x: NDArray
	y: NDArray
	z: NDArray

	@property
	def shape(self) -> tuple[int, int]:
		return self.x.shape

	def __post_init__(self):
		assert self.x.shape == self.y.shape == self.z.shape
		sx, sy = self.x.shape
		assert sx == 1

	def __add__(self, other: "Vec3NumpyArray") -> "Vec3NumpyArray":
		return Vec3NumpyArray(
			self.x + other.x,
			self.y + other.y,
			self.z + other.z
		)
	
	def __sub__(self, other: "Vec3NumpyArray") -> "Vec3NumpyArray":
		return Vec3NumpyArray(
			self.x - other.x,
			self.y - other.y,
			self.z - other.z
		)

	def dot(self, other: "Vec3NumpyArray") -> NDArray:
		return self.x * other.x + self.y * other.y + self.z * other.z

	def cross(self, other: "Vec3NumpyArray") -> "Vec3NumpyArray":
		return Vec3NumpyArray(
			self.y * other.z - self.z * other.y,
			self.z * other.x - self.x * other.z,
			self.x * other.y - self.y * other.x
		)
	
	@cached_property
	def rho(self) -> NDArray:
		return (self.x**2 + self.y**2 + self.z**2)**0.5
	
	@cached_property
	def normalized(self) -> "Vec3NumpyArray":
		rho = self.rho
		if rho == 0:
			raise ValueError("Cannot normalize a zero vector")
		return Vec3NumpyArray(
			self.x / rho,
			self.y / rho,
			self.z / rho
		)