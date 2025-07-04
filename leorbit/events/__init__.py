from abc import ABC, abstractmethod
from enum import Enum
from typing import Literal, Self
from itertools import batched, combinations

from pint import Quantity

from leorbit.mathematics.units import Q_
from leorbit.physics.time import Time
from leorbit.physics.time_interval import MIN_DURATION, TimeInterval

EventInfo = dict[str, float | int | str] | None

class StartStop(Enum):
	START = 0
	STOP = 1

def _get_type_step(before: bool, after: bool) -> StartStop | None:
	if before == after:
		return None
	elif not before and after:
		return StartStop.START
	elif before and not after:
		return StartStop.STOP

class Event(ABC):
	def __init__(self, *args, **kwargs):
		self.__pieces: dict[TimeInterval, EventInfo] = {}
          
	def _check_for_overlaps(self) -> bool:
		"""Returns `False` if computed time intervals do overlap. Returns `True` if everything ok."""
		for (ti1, ti2) in combinations(self.__pieces.keys(), 2):
			if ti1.intersection(ti2) is not None:
				return False
		return True
	
	def _build_by_predicate(self, predicate_values: list[bool], start: Time, dt: Quantity):
		steps = len(predicate_values)

		predicate_values = list(predicate_values)
		predicate_values_next = list(predicate_values)
		predicate_values_next.pop(0)
		predicate_values.pop(-1)

		flags_start_or_stop: list[tuple[StartStop, Time]] = [
			(type_step, start + i * dt)
			for i, before, after in zip(
				range(steps - 1), 
				predicate_values, 
				predicate_values_next
			)
			if (type_step := _get_type_step(before, after)) is not None
		]

		if not flags_start_or_stop:
			return

		if flags_start_or_stop[0][0] == StartStop.STOP:
			flags_start_or_stop.insert(0, (StartStop.START, start))
		if flags_start_or_stop[-1][0] == StartStop.START:
			flags_start_or_stop.append((StartStop.STOP, start + steps * dt))

		assert len(flags_start_or_stop) % 2 == 0

		for batch in batched(flags_start_or_stop, 2):
			(start_str, start_time), (stop_str, stop_time) = batch
			assert start_str == StartStop.START and stop_str == StartStop.STOP

			duration = stop_time.delta(start_time)
			interval: TimeInterval

			if duration < dt:
				interval = TimeInterval.make_ponctual(start_time)
			else:
				interval = TimeInterval(start_time, stop_str, dt)

			self.__pieces.setdefault(interval, None)
			
	@abstractmethod
	def compute(self, on: TimeInterval):
		...