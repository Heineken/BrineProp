# Repo has been moved to [GFZ-Git](https://gitext.gfz-potsdam.de/francke/BrineProp/)!

BrineProp
=========
BrineProp is a modelica package that calculates the thermodynamic properties of a specified brine, i.e., an aqueous solution of salts and gases, with a potential gas phase, including degassing/evaporation and solution/condensation.

An Excel version is available in the [download section](https://github.com/Heineken/BrineProp/releases). Its VBA code can also be found in `/VBA`.

Compatibility
-------------
* Works in Dymola 2014 FD01/2015 with MSL 3.2.1
* Works in Dymola 2013/2014 FD01 with MSL 3.2 (changes needed, see "Installation" below)

* Works in OpenModelica 1.9.1Beta2 and JModelica 1.14 via FMU.
* Works partly in OpenModelica, for details see Library documentation (in code).

Installation
------------
Download, unzip and open BrineProp/package.mo

For use with MSL 3.2 you need to make changes marked with "MSL 3.2.1" in PartialMixtureTwoPhaseMedium.

Getting started
---------------
Run models from `BrineProp.Examples.*`.

Documentation
-------------
* in the package (info annotation)
* in PhD thesis *[Thermo-hydraulic model of the two-phase flow in the brine circuit of a geothermal power plant](http://nbn-resolving.de/urn:nbn:de:kobv:83-opus4-47126)*

## Development and contribution
Feedback/contributions welcome to [@Heineken](https://github.com/Heineken) or as [pull request](https://github.com/Heineken/BrineProp/pulls).

## License
Licensed by Helmholtz Centre Potsdam, GFZ German Research Centre for Geosciences under the [Modelica License 2](https://www.modelica.org/licenses/ModelicaLicense2) or newer.

Copyright &copy; 2009-2014 Henning Francke

This Modelica package is free software and the use is completely at your own risk;
it can be redistributed and/or modified under the terms of the [Modelica License 2](https://www.modelica.org/licenses/ModelicaLicense2) or newer.
For license conditions (including the disclaimer of warranty) visit [https://www.modelica.org/licenses/](https://www.modelica.org/licenses).
