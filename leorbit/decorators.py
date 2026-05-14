from __future__ import annotations

import functools
from typing import Callable, ParamSpec, TypeVar, overload, cast

P = ParamSpec("P")
R = TypeVar("R")
_DEFAULT = object()


@overload
def lru_cache(func: Callable[P, R], /) -> Callable[P, R]:
    ...


@overload
def lru_cache(maxsize: int | None = 128, /, typed: bool = False) -> Callable[[Callable[P, R]], Callable[P, R]]:
    ...


@overload
def lru_cache(*, maxsize: int | None = 128, typed: bool = False) -> Callable[[Callable[P, R]], Callable[P, R]]:
    ...


def lru_cache( # type: ignore
    func: Callable[P, R] | int | None | object = _DEFAULT,
    /,
    *,
    maxsize: int | None = 128,
    typed: bool = False,
) -> Callable[P, R] | Callable[[Callable[P, R]], Callable[P, R]]:
    """Type-preserving wrapper around functools.lru_cache.

    Supports both forms:
    - ``@lru_cache``
    - ``@lru_cache(4096)``
    - ``@lru_cache(maxsize=4096, typed=True)``
    """
    if callable(func):
        wrapped = functools.lru_cache(maxsize=maxsize, typed=typed)(func)
        return cast(Callable[P, R], wrapped)

    if func is _DEFAULT:
        resolved_maxsize = maxsize
    else:
        resolved_maxsize = cast(int | None, func)

    def _decorator(func: Callable[P, R]) -> Callable[P, R]:
        wrapped = functools.lru_cache(maxsize=resolved_maxsize, typed=typed)(func)
        return cast(Callable[P, R], wrapped)

    return _decorator
