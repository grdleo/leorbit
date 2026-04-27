from leorbit.api import get_satellite, get_passes, Timestamp, Quantity, GPS, VisibleFromEarthLocationEvent, TimeInterval

# Let's compute the position of a LEO satellite!

iss = get_satellite(25544) # 25544: ISS NORAD Cat ID
now = Timestamp.now()
c = iss.coordinates(now) # We compute the ISS position at given time
c.gps()
'<GPS:  003° 59′ 29″E,  042° 29′ 24″S>'

# Now that we have computed its position for the given time, we can project
# it into any frame we want!

# Coordinates in GCRF!
c.gcrf().human_repr("km")
'<Vector3 x=-6381.3991361611925 y=146.7020812274849 z=-2333.735440753873 [km]>'

# Coordinates in ITRF!
c.itrf().human_repr("km")
'<Vector3 x=-2639.11482315943 y=5811.957448727188 z=-2333.735440753873 [km]>'

# Horizontal coordinates in any Earth local frame!
# For example, let's try in Paris.
gps_paris = GPS(
    longitude=2.333333 * Quantity.degree, 
    latitude=48.866667 * Quantity.degree, 
    altitude=0 * Quantity.meter
)
c.horizontal(gps_paris.earth_local_frame)
'<Horizontal: Azimuth:  087° 21′ 36″, Altitude: - 058° 36′ 26″>'

# Ugh, negative altitude, it means it is not visible currently...
# Want to get the ISS passes for the next 7 days?

timeline = TimeInterval(
    start=now,
    stop=now + 7 * Quantity.day,
    dt=5 * Quantity.second
)
first_pass, *others = get_passes(iss, timeline, gps_paris)
"""<TimeInterval from: \'2026-04-28 at 00:19:19\' to: \'2026-04-28 at 00:28:54\' dt: 5s>"""

# And finally let's export the horizontal coordinates of the first pass to CSV!
iss.trajectory(first_pass).horizontal(gps_paris.earth_local_frame).to_csv(first_pass)
""""
timestamp,azimuth,elevation,range
2026-04-28T00:19:21.197208+00:00,-159.7,0.3,6790.4
2026-04-28T00:19:26.197208+00:00,-160.1,0.5,6790.4
2026-04-28T00:19:31.197208+00:00,-160.6,0.8,6790.4
2026-04-28T00:19:36.197208+00:00,-161.0,1.1,6790.3
...
"""