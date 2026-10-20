# OPAL — Orbitals, Plotting, Analysis and Layout

OPAL is a small suite of browser-based tools for working with OpenMolcas:
building molecules, looking at orbitals and densities, inspecting RASSI results,
comparing geometries and planning RASSCF calculations. Each tool is a single
HTML file. There is nothing to install, nothing to build, and no server: you open
the file in a web browser and it works, with or without an internet connection.
Files you open stay on your computer; the browser reads them locally and nothing
is uploaded.

## The tools

| File | Tool | What it does |
|---|---|---|
| `molbuild.html` | MolBuild | Draw a molecule in 2D, clean it up, and get 3D coordinates (UFF-optimised); reads and writes XYZ files. |
| `orbviz.html` | OrbViz | Render molecular orbitals, electron densities and density differences from OpenMolcas HDF5 and Molden files. |
| `rassiviz.html` | RassiViz | Explore RASSI states, transition dipole moments and transition densities from HDF5 files; includes an absorption spectrum viewer. |
| `molviz.html` | MolViz | Ball-and-stick viewer for XYZ, Molden and OpenMolcas HDF5 geometries, with bond editing, atom labels, transition dipole vectors and image export. |
| `geodiff.html` | GeoDiff | Compare two structures (XYZ, Molden or HDF5): RMSD and best-fit superposition, ranked bond, angle and dihedral differences. |
| `csfscalc.html` | CSFsCalc | Count configuration state functions and determinants for a RASSCF active space, with an estimate of the memory needed. |
| `index.html` | Launcher | A start page that links to the tools above. |

## How to use

Open `index.html` in a browser (double-click it, or drag it into a browser
window) and pick a tool. Each tool can also be opened directly on its own. To
load a calculation, use the tool's file button or drop the file onto the window.

## Browsers

OPAL is written to run on Safari 14.1 and later (macOS, iPadOS and iOS) and on
current versions of Chrome, Edge and Firefox. The 3D tools need WebGL, which all
of these provide.

The launcher also contains a link to a Google Colab notebook for running
OpenMolcas without installing it. That one link needs an internet connection;
the tools themselves do not.

## Third-party code

The 3D tools include a copy of [three.js](https://threejs.org) r128, which is
distributed under the MIT licence, so that each file is complete on its own.

## Licence

OPAL is free software, distributed under the GNU Lesser General Public License,
version 2.1, or (at your option) any later version. See
<https://www.gnu.org/licenses/>.

## Author

Luca De Vico, Università degli Studi di Siena.
