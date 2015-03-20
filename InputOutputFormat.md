# Introduction #

To start the simulation, the initial locations and velocities of all the stars in the system must be read in from a file.   At any point during the simulation, the current state of the galaxy can also be outputted to a file.  Both of these files follow the same format so it is simple to begin a simulation where an old one left off.


# Details #

Each line in the file is for a specific star.  The data included in the file are:
  * x position
  * y position
  * z position
  * x velocity
  * y velocity
  * z velocity

These are stored as comma separated values per line. A single line is defined as follows:

```
x_position, y_position, z_position, x_velocity, y_velocity, z_velocity
```


## Future Considerations ##
Mass is currently not included in the star's information since we have made the assumption that all stars are the same mass.  If we choose to add a mass, this file format should be modified accordingly.