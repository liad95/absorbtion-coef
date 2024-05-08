import sys
import numpy as np
from kwave.data import Vector
from kwave.utils.kwave_array import kWaveArray
from kwave.kgrid import kWaveGrid
from kwave.kmedium import kWaveMedium
from kwave.ksource import kSource
from kwave.ksensor import kSensor
from kwave.utils.data import scale_SI
from kwave.utils.mapgen import make_disc, make_ball
from kwave.kspaceFirstOrder2D import kspaceFirstOrder2D
from kwave.reconstruction import
from kwave.options.simulation_execution_options import SimulationExecutionOptions
from kwave.options.simulation_options import SimulationOptions

Nx = 64      # number of grid points in the x (row) direction

x = 1e-3     # size of the domain in the x direction [m]

dx = x / Nx  # grid point spacing in the x direction [m]

sound_speed = 1500

# size of the initial pressure distribution
source_radius = 2              # [grid points]

# distance between the centre of the source and the sensor
source_sensor_distance = 10    # [grid points]

# time array
dt = 2e-9                      # [s]
t_end = 300e-9                 # [s]

medium2 = kWaveMedium(sound_speed=sound_speed)

# create the k-space grid
kgrid2 = kWaveGrid([Nx, Nx], [dx, dx])

# create the time array using an integer number of points per period
kgrid2.setTime(int(np.round(t_end / dt)), dt)

# create instance of a sensor
sensor2 = kSensor()

# set sensor mask: the mask says at which points data should be recorded
sensor2.mask = np.zeros((Nx, Nx), dtype=bool)

# define a single sensor point
sensor2.mask[Nx // 2 + source_sensor_distance, Nx // 2] = True

# set the record type: record the pressure waveform
sensor2.record = ['p']

# make a source object
source2 = kSource()
source2.p0 = make_disc(Vector([Nx, Nx]), Vector([Nx // 2, Nx // 2]), source_radius, plot_disc=False)

simulation_options = SimulationOptions(
    data_cast='single',
    save_to_disk=True)

execution_options = SimulationExecutionOptions(
    is_gpu_simulation=True,
    delete_data=False,
    verbose_level=2)

# run the simulation
sensor_data_2D = kspaceFirstOrder2D(
    medium=medium2,
    kgrid=kgrid2,
    source=source2,
    sensor=sensor2,
    simulation_options=simulation_options,
    execution_options=execution_options)

t_sc, t_scale, t_prefix, _ = scale_SI(t_end)
import matplotlib.pyplot as plt
_, ax1 = plt.subplots()
ax1.plot(np.squeeze(kgrid2.t_array * t_scale), sensor_data_2D['p'] / np.max(np.abs(sensor_data_2D['p'])), 'r-')
ax1.set(xlabel= f"Time [{t_prefix}s]", ylabel='Recorded Pressure [au]')
ax1.grid(True)
plt.show()