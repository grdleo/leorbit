# walkthrough

# LEOrbit walkthrough

Welcome!

## The need for tensors

When you make (literally) astronomical calculation, well, as you can expect, there are **a lot** of maths involved. So, **of course** we are using `numpy`, which makes all the calculations absurdly fast. But `numpy` itself is not enough if you want to create a robust library for this purpose. Why?

Sometimes you have a velocity (km/s) 3D vector; Sometimes you have an array of scalars which are supposed to be angles; Sometimes you have dimensionless 3x3 matrices...

To solve this issue, and so that the developper/user always know what they are dealing with, we must create a `Tensor` wrapper:

- It must encode the *kind* of the tensor (is it a *scalar?* a *vector3?* a *matrix33?*)
- It must encode the *size* of the tensor (aka. is it a *unique* vector3? or does it holds multiple vector3?)
- It must encode the *dimension* (and units) of the tensor (aka. it can have the dimension of a *length*, a *time*, a *velocity*, *acceleration, force,* etc.)

But basically it wraps the numpy array representing the tensor, and exposes mathematical operators and functions.


```python
from leorbit.api import Qty, scalar, vector3, matrix33
```


```python
a = scalar(42.3) # a simple dimensionless scalar
a
```




    <Scalar [42.3] [dimensionless]>




```python
2 * a + 12/5 # arithmetic operations are supported, and the result is still a scalar
```




    <Scalar [87.] [dimensionless]>




```python
b = scalar("12 meter") # a scalar with a length unit

# a scalar with a velocity unit (using pint to generate the units)
c = scalar(2**.5) * (Qty.m / Qty.s)
b, c
```




    (<Scalar [12.] [meter]>, <Scalar [1.41421356] [meter / second]>)




```python
32/3 * b + c * (5.3 * Qty.hour)
```




    <Scalar [135.49533188] [meter]>




```python
b + (1.2 * Qty.minute) # Careful! This is not a valid operation, because the unit of b is not compatible with the unit of 1.2 (which is dimensionless). This will raise an error.
```


    ---------------------------------------------------------------------------

    ValueError                                Traceback (most recent call last)

    Cell In[6], line 1
    ----> 1 b + (1.2 * Qty.minute) # Careful! This is not a valid operation, because the unit of b is not compatible with the unit of 1.2 (which is dimensionless). This will raise an error.


    File ~/Code/leorbit/leorbit/leorbit/mathematics.py:340, in Tensor.__add__(self, right)
        339 def __add__(self, right: Tensor | RealNumber | pint.Quantity | pint.Unit) -> Tensor:
    --> 340     return self.perform_binary_operation(right, TensorBinaryOperator.ADD)


    File ~/Code/leorbit/leorbit/leorbit/mathematics.py:387, in Tensor.perform_binary_operation(self, other, op)
        385 units = op.dimension_result(self._units, other._units)
        386 if units is None:
    --> 387     raise ValueError("Dimensions incompatible for given operator")
        389 return Tensor(
        390     data=op.operator(self._data, other._data), 
        391     units=units
        392 )


    ValueError: Dimensions incompatible for given operator



```python
v = vector3(1, 2, 3)
w = vector3(-4, 5, 2)
v + w, v.vector3.dot(w.vector3), v.vector3.cross(w.vector3)

m = matrix33(
    0, 0, 1,
    1, 0, 0,
    0, 1, 0,
)
m.matrix33 @ v.vector3

m_units = matrix33(
    1, 2, 3,
    4, 5, 6,
    7, 8, 9,
).with_units("meter")
a = vector3(11, 22, 33).with_units("s**2")
m_units.matrix33 @ a.vector3
```




    <Vector3 x=154.0 y=352.0 z=550.0 [meter * second ** 2]>



**For experimented Pint users:** if you intend to use quantities, do not forget to use the unit registry provided in `leorbit.api`!

## Time handling in LEOrbit

Time is everywhere in orbital mechanics, and tiny timing errors can quickly become kilometer-level position errors. That is why LEOrbit does not rely on raw floats or ad-hoc datetime arithmetic for core computations.

The `Timestamp` object represents a precise instant and exposes conversion helpers and formatting tools. It is designed to stay explicit and predictable when you move between human-readable time and computational time scales.

`TimeInterval` represents a sampled time range with a start, a stop, and a fixed step `dt`. This is the natural container for trajectory propagation and event search: instead of manually building loops over datetimes, you define the interval once and LEOrbit gives you a coherent timeline.

In practice, this gives you three major advantages:

- Consistency: every sampled value in a trajectory is tied to a well-defined timestamp.
- Safety: operations use typed quantities (`Quantity.second`, `Quantity.minute`, etc.), reducing unit mistakes.
- Interoperability: the same time primitives are reused across propagation, coordinate transformations, and event computation.


```python
from leorbit.api import Timestamp, TimeInterval
from datetime import timedelta

now = Timestamp.now()
next_15_minutes = TimeInterval(
    start=now,
    stop=now + timedelta(seconds=15 * 60),
    dt=timedelta(seconds=5),
)

next_15_minutes.start.isoformat, next_15_minutes.stop.isoformat, next_15_minutes.steps
```




    ('2026-05-14T20:11:17.457318+00:00', '2026-05-14T20:26:17.457318+00:00', 180)



## Let's talk about frames

Thanks to Newton, we can very easily compute the trajectory of an object affected by gravitationnal forces.

$$ \vec{F} = m\vec{a} \quad \text{and for gravity} \quad \vec{F}_{12} = -G\frac{m_1m_2}{r^2}\hat{r} $$

But in order to actually **make** calculations, well, we need to project the position and velocity in a frame. But... Which frame? I mean, everything is moving: the Earth revolves around itself, and around the Sun... Help!

This walkthrough's purpose is not to rewrite a full course on astrodynamics, so you should definitely take 5min to check a few links:

- [Wikipedia article on astronomical coordinates system](https://en.wikipedia.org/wiki/Astronomical_coordinate_systems)
- ...

And now we list the useful frames used in LEOrbit:

### ITRS frame
- Is a $(x,y,z)$ frame of origin the center of Earth. 
- Its $z$ axis is along the Earth's rotation axis, and points toward north.
- Its $x$ axis points toward the origin of GPS coordinates (Guinea golf, Atlantic ocean).
- And $y = z \times x$

### GCRS frame
- Is a $(x,y,z)$ frame of origin the center of Earth. 
- Its $z$ axis is along the Earth's rotation axis, and points toward north.
- Its $x$ axis points toward the [Vernal point](https://en.wikipedia.org/wiki/First_point_of_Aries), aka a fixed point in the deep sky.
- And $y = z \times x$

So we can say that ITRS = GRFS with a rotation along the $z$ axis, of a angle depending on the time of the day

Those two are considered as "absolute frames" in LEOrbit, which means two things:
- Their definition is not ambiguous, and the same for everyone.
- There is also a well defined transformation, that only depends on the epoch (the time), to convert a coordinates from/to.
- They are used as reference to construct "relative frames"

### Earth relative frames
ITRS and GCRS are great, but the actual observations are made by humans, and, as you may know, humans live *on the surface of the Earth.*

It means that, for every point on the surface of Earth, you can define a local frame:

- $(x,y,z)$ orthonormal frame
- $z$ points to the [zenith](https://en.wikipedia.org/wiki/Zenith)
- $x$ points to the north
- $y = x \times z$ (which makes a [non right-handed basis](https://en.wikipedia.org/wiki/Right-hand_rule) but works better with the definition of the horizontal coordinates)

A Earth local frame is simply defined by an affine transformation from the ITRS absolute reference frame.
LEOrbit users can easilly create such relative frames, from a simple GPS coordinates!




```python
from leorbit import GPS, EarthLocalFrame, N, E

gps_paris = 2.333333 * E + 48.866667 * N
local_frame_paris = gps_paris.earth_local_frame 
```

## Coordinates

Now that you are familiar with frames, you should understand more precisely what *coordinates* are: coordinates are simply a position (or a position-velocity couple) at a given time. A coordinates *physically* exists without the need to define any frame.
But if we need to make actual *calculations* with it (spoiler alert, we definitely need to), we have no choice but to make **projections** of the coordinates in a given frame.
If we need multiple frames, and not just a single one, it is only because for a given calculation, there is always a relevant frame where the calculation is easy and makes perfect sense.

But... everything we just talked about shoud be, in most of the cases, completely transparent for the end user!!

*Why should I care about ITRS, GCRS and stuff if my goal is just to compute passes of my satellite of choice??*

Well this is why LEOrbit handles the concept of *coordinates* in a very user-friendly way: a `Coordinate` is just an object representing the position (and velocity) of an object at a given time. All the frames and projections are handled internally, without any intervention of user. Yet if you need the actual **projection** of the coordinates in your frame of choice, well, of course you can still get it, and very easily!

Also LEOrbit handles `Trajectory` objects, which are simply a collection of coordinates in a given time interval.

## What is an orbit?

If you model a system with one dominant body (for example Earth) and one much smaller object (for example a satellite), Newtonian gravity leads to a conic trajectory. In the bounded case, that conic is an ellipse: this is exactly the classical Keplerian orbit model.

In that ideal two-body model, the shape and orientation of the orbit are fixed by five orbital elements:
- semi-major axis $a$
- eccentricity $e$
- inclination $i$
- right ascension of ascending node $\Omega$
- argument of pericenter $\omega$

The sixth element, mean anomaly $M$, evolves almost linearly with time and tells you where the object is on that fixed ellipse.

Real trajectories are never perfectly two-body. Atmosphere, Earth oblateness, third-body gravity (Sun/Moon), solar radiation pressure, and maneuvers all perturb the motion. In practice, this means the orbital elements can drift with time, sometimes slowly, sometimes very quickly depending on altitude and mission profile.


```python
from leorbit.api import get_satellite, Timestamp, TimeInterval

iss = get_satellite(25544, log=False)
t0 = Timestamp.now()
state_now = iss.coordinates(t0)

short_arc = TimeInterval(
    start=t0,
    stop=t0 + 10 * Qty.min,
    dt=30 * Qty.s,
)
traj = iss.trajectory(short_arc)

state_now.gcrf().human_repr("km"), traj
```

## The two types of bodies

All the moving objects in the universe have their trajectory completely defined by the law of Newton and gravity. 
These objects, in LEOrbit, inherit from `SkyObject`.

In our solar system, all the important bodies (the Sun, planets, moons, ...) are now stable since for very VERY long time. 
It means their coordinates in our sky can be computed with extreme precision for thousand of years, with just simple formulaes.
Such objects in LEOrbit are considered as `Body` objects. 
Currently only `Sun` and `Moon` have been hardcoded. But all the known ones should be in the end retrievable from ephemerids.

But for man-made satellites, aka small objects that orbit planets with an atmosphere, well, computing their coordinates is much harder, mostly because of friction with the air (yes, even when you are in space!). 
The difficulty is not really in the computation itself, but you cannot make predictions that are accurate for more than a few weeks at most. 
Why? 
Because drag depends on atmospheric density, and atmospheric density changes with altitude, local solar time, geomagnetic conditions, and solar activity. On top of that, small modeling errors accumulate over time, and operational satellites can perform maneuvers that are not always known in advance.

In LEOrbit satellites are `Satellite` objects. And we need to use **propagators** to compute their trajectories.

## The propagators in LEOrbit

A propagator is the algorithm that advances a satellite state from one epoch to another. In other words, it is the physical and numerical model that answers: "given what I know now, where will the object be later?"

Different contexts require different propagators. Sometimes you just want a lightweight, deterministic trajectory from a fixed analytical model; sometimes you want the best practical short-term prediction for Earth satellites from public orbital data.

LEOrbit currently exposes two main families:

- `NoPropagator`: it simply propagates the satellite *along* the fixed ellipsis (aka. updates $ M $ only). It implies that the satellite moves only due to Earth gravity, without any perturbations. Therefore it has no real practical use, but handy to simulate perfectly Keplerian trajectories...
- `SGP4`: the standard practical propagator for TLE/GP data. It is the right choice for most operational LEO tracking use cases, with good short-term accuracy at low computational cost.

Choosing a propagator is therefore a trade-off between physical fidelity, required horizon, and data availability.
The goal is to implement [all the main propagators,](https://en.wikipedia.org/wiki/Simplified_perturbations_models) so for now `SGP4` is really the only option, and keep in mind it is only relevant for LEO satellites!!

## Verbose example

A complete low-level flow for ISS typically looks like this:

1. Fetch the latest GP/TLE set from Celestrak for NORAD 25544.
2. Build an `SGP4` propagator from that element set.
3. Create a `Satellite` object with this propagator.
4. Define a `TimeInterval` and ask for coordinates/trajectory over that interval.
5. Convert to the frame you need (ITRS, GCRS, or local horizontal) and inspect or export results.

This is very explicit and great to understand what happens under the hood.

Then, in most real workflows, you can do the same in one line with `get_satellite(25544)`: it wraps data retrieval + propagator setup for you. For common LEO operations this convenience path is usually the best starting point, and you can still switch to the verbose approach whenever you need more control.


```python
from leorbit.api import get_satellite, TimeInterval, Timestamp, GPS

iss = get_satellite(25544, log=False)
now = Timestamp.now()
timeline = TimeInterval(
    start=now,
    stop=now + 1 * Qty.hour,
    dt=20 * Qty.s,
)

gps_paris = GPS(
    longitude=2.333333 * Qty.deg, 
    latitude=48.866667 * Qty.deg, 
    altitude=0 * Qty.m
)

traj = iss.trajectory(timeline)
horizontal_traj = traj.horizontal(gps_paris.earth_local_frame)

traj, horizontal_traj
```

## Computing events in LEOrbit

Trajectory samples are useful, but many practical questions are event-based: "when does the satellite rise?", "when does visibility start?", "when is max elevation?". That is why LEOrbit separates continuous trajectories from event detection tools.

`TimeMap` is the generic sampled timeline container used to evaluate time-dependent quantities over an interval. `Event` is the abstraction that defines a boolean condition over time (for example "visible from this observer").

Event computation works by evaluating the event condition on the timeline, finding transitions (`False -> True` and `True -> False`), then extracting the corresponding intervals and key instants. Internally, LEOrbit can refine boundaries depending on step size and event logic, but conceptually it is a robust transition search over a typed time grid.

To create your own event, define a class that implements the event predicate for a given object/time context, then run it on a `TimeInterval`. This lets you encode domain-specific criteria (illumination constraints, elevation masks, custom station rules, etc.) while reusing the same timeline and interval extraction machinery as built-in events.


```python
from leorbit.api import (
    get_satellite,
    TimeInterval,
    Timestamp,
    GPS,
    VisibleFromEarthLocationEvent,
 )

iss = get_satellite(25544, log=False)
now = Timestamp.now()
timeline = TimeInterval(
    start=now,
    stop=now + scalar("2 day"),
    dt=scalar("10s"),
 )

gps_paris = GPS(
    longitude=scalar("2.333333 deg"),
    latitude=scalar("48.866667 deg"),
    altitude=scalar("0m"),
 )

visibility_event = VisibleFromEarthLocationEvent(
    trajectory=iss.trajectory(timeline),
    gps_observer=gps_paris,
    altitude_angle_min=scalar("0°"),
)

visibility_event.visible_intervals[:3], len(visibility_event.visible_intervals)
```




    ([<TimeInterval from: '2026-05-08 at 22:02:53' to: '2026-05-08 at 22:13:33' dt: 10s>,
      <TimeInterval from: '2026-05-08 at 23:39:53' to: '2026-05-08 at 23:50:33' dt: 10s>,
      <TimeInterval from: '2026-05-09 at 01:16:43' to: '2026-05-09 at 01:27:33' dt: 10s>],
     13)


