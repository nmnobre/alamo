.. role:: cpp(code)
   :language: c++
.. role:: bash(code)
   :language: bash

Voronoi Microstructure with Elasticity
======================================


This test problem illustrates the use of the phase field method when initialized with a Voronoi tesselation.
The following is the result in 2D; the problem can be run in 3D as well.

.. image:: voronoielastic.gif
   :scale: 75%
   :align: center

(Notice that there are 100 initial grains in the tesselation, but only 10 actual grains in the simulation.
Therefore, some grains will eventually merge together when the boundary between them disappears.)

.. literalinclude:: ../../../tests/VoronoiElastic/input
   :caption: VoronoiElastic input file (tests/VoronoiElastic/input)
   :language: makefile

Run the test problem (2D or 3D) via

.. code-block:: bash

      ./bin/alamo... tests/VoronoiElastic/input

This will generate an output directory :code:`tests/VoronoiElastic/output/`, and will re-name the old directory if it exists already.
Use VisIt to open :code:`tests/VoronoiElastic/output/output.visit`.

See also:

- :bash:`./src/alamo.cc`: Entry point for the solver
- :ref:`API-Integrator-PhaseFieldMicrostructure`: Integrator that explicitly evolves the order parameter
- :ref:`API-Operator-Elastic`: Operator to do elastic solves
- :ref:`API-IC-PerturbedInterface`: Parameterization of initial perturbation
