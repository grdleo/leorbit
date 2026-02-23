import numpy as np
import pytest

from leorbit.utils import truth_array_to_indices_intervals


@pytest.mark.parametrize(
    "truth_array, expected",
    [
        ([False, True, True, True, False, False, True, True, False], [(1, 3), (6, 7)]),
        ([False, False, False], []),
        ([True, True, True], [(0, 2)]),
        ([True], [(0, 0)]),
        ([False], []),
        ([False, True, False], [(1, 1)]),
        ([True, False, True], [(0, 0), (2, 2)]),
        ([False, True, True, False, True, False, True, True, True, False], [(1, 2), (4, 4), (6, 8)]),
    ],
)
def test_truth_array_to_indices_intervals_various_patterns(truth_array, expected):
    assert truth_array_to_indices_intervals(truth_array) == expected


@pytest.mark.parametrize(
    "truth_array, min_size, expected",
    [
        ([False, True, True, True, False, False, True, True, False], 1, [(1, 3), (6, 7)]),
        ([False, True, True, True, False, False, True, True, False], 2, [(1, 3), (6, 7)]),
        ([False, True, True, True, False, False, True, True, False], 3, [(1, 3)]),
        ([False, True, True, True, False, False, True, True, False], 4, []),
        ([True, True, False, True, True, True, False, True], 2, [(0, 1), (3, 5)]),
        ([True, True, False, True, True, True, False, True], 3, [(3, 5)]),
        ([True, True, False, True, True, True, False, True], 10, []),
    ],
)
def test_truth_array_to_indices_intervals_min_size_filtering(truth_array, min_size, expected):
    assert truth_array_to_indices_intervals(truth_array, min_size_intervals=min_size) == expected


@pytest.mark.parametrize(
    "arr",
    [
        np.array([0, 1, 1, 1, 0, 1, 1, 0, 1], dtype=bool),
        np.array([0, 1, 1, 1, 0, 1, 1, 0, 1], dtype=np.int8),
        np.array([0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0], dtype=np.float64),
        np.array([[False, True, True], [False, True, False]], dtype=bool),
    ],
)
def test_truth_array_to_indices_intervals_numpy_inputs(arr):
    got = truth_array_to_indices_intervals(arr, min_size_intervals=2)
    ref = truth_array_to_indices_intervals(np.asarray(arr, dtype=bool).reshape(-1), min_size_intervals=2)
    assert got == ref


@pytest.mark.parametrize(
    "iterable_input, expected",
    [
        ((x for x in [False, True, True, False]), [(1, 2)]),
        ((x for x in [True, False, True, True]), [(0, 0), (2, 3)]),
        (tuple([False, False, True, True, True]), [(2, 4)]),
    ],
)
def test_truth_array_to_indices_intervals_iterable_inputs(iterable_input, expected):
    assert truth_array_to_indices_intervals(iterable_input) == expected


def test_truth_array_to_indices_intervals_empty_and_invalid_min_size():
    assert truth_array_to_indices_intervals([]) == []

    assert truth_array_to_indices_intervals(np.array([], dtype=bool)) == []

    with pytest.raises(ValueError):
        truth_array_to_indices_intervals([True, False], min_size_intervals=0)

    with pytest.raises(ValueError):
        truth_array_to_indices_intervals([True, False], min_size_intervals=-3)
