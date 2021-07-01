# RWA

Coursework for AE4135 Rotor/Wake Aerodynamics 2020/21 at TU Delft. Written in Python with occasional notes in Jupyter notebooks.

## Assignment 1: BEM

Blade element momentum code. Investigating effect of tip-speed ratio on non-yawed and yawed rotor. Geometry optimisation based on Aerodynamic Design of Wind Turbine Rotors (2013) by Bak.

- Mesh sensitivity analysis
- Non-yawed and yawed rotors rotor
  - Angle of attack and inflow angle
  - Induction factors
  - Loads
  - Thrust and torque
  - Circulation
- Non-yawed
  - Stagnation enthalpy
- Geometry optimisation

## Assignment 2: Lifting Line

Wind turbine vortex lifting line model in Python. Single rotor and double rotor configurations investigating, amongst other things, streamtube independence.

- Sensitivity study
  - Wake convection speed
  - Radial discretisation
  - Azimuthal discretisation
  - Wake length
- Single rotor case
  - Angle of attack and inflow angle
  - Circulation
  - Induction factors
  - Loads
  - Thrust and power coefficients
- Double rotor case
  - Effect of distance
  - Effect of phase

## Assignment 3: Unsteady Aerodynamics

Unsteady, incompressible, potential flow, thin-airfoil, vortex panel code based on Section 13.10 of "Low-Speed Aerodynamics" by Joseph Katz & Allen Plotkin.

- Sensitivity study
  - Number of panels
  - Time step
- Steady aerofoil
- Unsteady aerofoil
- Steady aerofoil with flap
- Unsteady aerofoil experiencing gust
