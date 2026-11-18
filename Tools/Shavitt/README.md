# Shavitt (version 0.9)

Draws the Shavitt graph (the split-GUGA distinct-row table) of a RASSCF/CASSCF
wave function, with the walks of its leading CSFs.

RASSCF writes the graph to the group `SHAVITT` of `$Project.rasscf.h5`
(`src/rasscf/shavitt2h5.F90`); this tool reads that group and the
`CI_VECTORS` of the same file.

## Requirements

Python 3.8+, `numpy`, `h5py`, `matplotlib`; for the graphical interface also a
Qt binding (PySide6, PyQt6, PySide2 or PyQt5, e.g. `pip install PySide6`).
Optional, for the orbital pictures: [Pegamoid](https://gitlab.com/Jellby/Pegamoid)
(`pip install pegamoid`, which brings VTK and qtpy). Pegamoid is a separate
program (GPL v3); Shavitt only loads it when it is installed. It is looked for
in `$PEGAMOID`, next to the Python that runs Shavitt, then on the `PATH`. Use
Pegamoid 2.13 with NumPy 2 (2.9 uses `np.math`, which NumPy 2 no longer has).

## Graphical interface

    shavitt.py h2o.rasscf.h5        # or just shavitt.py, then File > Open

* **Wave function**: root, energy, threshold |c| for the CSF table.
* **Leading CSFs**: select rows (Ctrl/Shift-click) to highlight their walks;
  *Select top N* picks the N largest.
* **Highlight a walk**: by CSF number, or by step vector (`0 u d 2`, level 1
  first, blanks ignored). A step vector of the wave function is located as its
  CSF; any other walk of the graph (another irrep) is drawn as well and marked.
  At most 8 walks; each keeps its colour while it stays highlighted.
* **Display**: number the vertices, Paldus (a,b,c) labels, midlevel line,
  midvertices (squares), step d on each arc, active orbital of each level,
  fading of the graph behind walks, line width by |c|, subgraph of the
  selected vertex. *Reset view* undoes zoom/pan.
* **Orbitals** (Pegamoid): resting the mouse on a level label on the left
  (e.g. `4  b1`) opens a pop-up with the isosurfaces (+-0.05) of the active
  orbital of that level, with its MO number in its irrep, RAS space,
  occupation and energy. Clicking the label opens the orbital in Pegamoid's
  own window. Switch off with *Orbital pop-ups* (View tab). The orbitals are
  the (average natural) orbitals of the rasscf.h5 file.
* **Hover**: resting the mouse on a vertex shows its (a,b,c), number, level,
  N and S in a tooltip.
* **Vertex**: click a vertex for L, N, (a,b,c), S, the number of walks above,
  below and through it, its midvertex number and its arcs; the walks through it
  are darkened.
* **Loops** (tab): while this tab is shown the CSF walks are hidden (they
  return on the other tabs). The loops of the one-body coupling coefficients
  <L|E_ij|R> in the split graph. Click a midvertex in the graph (*Pick
  midvertices with the mouse*, on by default; click a partner of it for the
  crossing loops, the selected one again to clear) or choose it in the list:
  * one midvertex: its **upper loops** (both walks between the graph head and
    the midvertex, levels i > midlevel) and **lower loops** (between the
    midvertex and the tail, j <= midlevel);
  * a second, partner midvertex (same level, N one less, b differing by 1;
    only the valid partners are offered): the **upper (U\*) and lower (L\*)
    partial loops** of the loops crossing the midlevel there. Every U\* x L\*
    pair is a loop; its coupling coefficient is the product of the two values.

  Filters: head level j, tail level i, upper/lower, irrep i = irrep j. With
  nothing selected in the list all its loops are drawn (their bra and ket
  arcs); select loops to see them in detail, with the walk parts common to bra
  and ket darkened and the value of <L|E_ij|R> (for one U\* and one L\*: the
  product). *Common walk* steps through the walks bra and ket share outside
  the loop (*all*: draw them together): each one embeds the loop in another
  pair of complete walks, with the same coupling coefficient; the note under
  the graph gives the two step vectors and their CSF numbers (if of the target
  irrep). Bra (L) is solid violet, ket (R) dashed orange. The segment values
  are those of SGUGA (MKSEG, ASTA); `tests/check_loops.py` verifies them
  against DENSITY_MATRIX and the enumeration against all pairs of walks.
* **Segments** (tab): click a vertex (or give its number) for all the
  segments of the singlet generators that start at it, as SGUGA tabulates
  them (MKSEG, ISGM/VSGM): the vertex as top-left (bra) vertex, and
  optionally as top-right (ket) vertex of a partner pair. Each row: SGUGA
  segment number (1-34), type (walk, head, mid db=-1/+1, tail, E_ii weight),
  role L/R, bra/ket step, value, and the vertices top L,R -> bottom L,R; the
  tooltip gives the classes of the vertex pairs and the value formula in b of
  the top-left vertex. Nothing selected: all listed segments are drawn;
  selected segments are drawn in detail (bra solid violet, ket dashed orange,
  corners ringed, number and value). CSF walks and loops are hidden on this
  tab.
* **Export**: PDF, PNG or SVG of the graph as shown (zoom included), with the
  dpi for PNG. The toolbar above the graph zooms and pans.

## Command line (no interface)

    shavitt.py h2o.rasscf.h5 --text                  # summary and leading CSFs
    shavitt.py h2o.rasscf.h5 -o graph.pdf            # draw it (pdf, png, svg)
    shavitt.py h2o.rasscf.h5 -r 2 -t 0.1 -n 8 -o graph.png --numbers --midvertices
    shavitt.py h2o.rasscf.h5 --csf 18 27 -o g.pdf    # given CSFs (1-based)
    shavitt.py h2o.rasscf.h5 --walk 2ud2ud0 -o g.pdf # given step vectors
    shavitt.py h2o.rasscf.h5 --vertices              # the distinct-row table

Display flags (also the initial settings of the interface): `--numbers`,
`--paldus`, `--midvertices`, `--no-midlevel`, `--steps`, `--orbitals`,
`--width-by-weight`, `--no-dim`, `--no-logo`, `--vertex V`; `--dpi` for bitmaps.
Loops: `--loops V` (midvertex V) or `--loops V W` (partner pair), with
`--loop-head J`, `--loop-tail I`; the loops are listed and drawn, without the
leading CSFs unless `--csf`/`--walk` is given.

## Layout of the drawing

Rows are the levels (head at the top, level 0 at the bottom), labelled with
the irrep of the orbital of the level. Columns are a regular Paldus (a,b)
lattice, x = (amax-a)(bmax+1) + (bmax-b), empty columns included: every step
has a fixed displacement (d=0 vertical, d=1 one column, d=2 bmax columns,
d=3 bmax+1 columns), so a run of equal steps is a straight line. Only the
columns in use are labelled. The dashed line marks the midlevel of the split
graph; the OpenMolcas logo sits in the lower right corner of the figure.

## The SHAVITT group of rasscf.h5

All integers; arrays are written in Fortran order (h5py shows the dimensions
reversed). Indices in the file are 1-based.

| Name | Dimensions (Fortran) | Content |
|---|---|---|
| `DRT` | `[NVERT,5]` | L, N, a, b, c of each vertex; vertex 1 is the head |
| `DOWN` | `[NVERT,4]` | vertex below for the steps d=0..3 (0: no arc) |
| `LEVEL_SYM` | `[NLEV]` | irrep of each level |
| `LEVEL_TO_ACT` | `[NLEV]` | active orbital of each level |
| `NRAS`, `NRASEL`, `NRASMX` | `[NSYM,NRSPRT]`, `[NRSPRT]` | RAS/GAS partitions and cumulative electron bounds |
| `NCSF` | `[NSYM]` | CSFs per total irrep |
| `NOW`, `IOW` | `[2,NSYM,NMIDV]` | upper(1)/lower(2) half walks per irrep and midvertex, offsets in `ICASE` |
| `NOCSF`, `IOCSF` | `[NSYM,NMIDV,NSYM]` | CSFs and offsets per (upper irrep, midvertex, total irrep) |
| `ICASE` | `[NWALK*NIPWLK]` | packed half-walk step vectors, `NPACK` 2-bit steps per integer, lowest level in the lowest bits |

Attributes: `NSYM NACTEL ISPIN NLEV NVERT MIDLEV MVSTA MVEND NMIDV IFRAS NRSPRT NPACK NIPWLK NWALK NUW`.

CSF order of `CI_VECTORS` (symmetry `LSYM`): blocks (midvertex outermost,
upper-walk irrep inside), within a block the upper half walk runs fastest.
`shavittlib.graph.ShavittGraph.step_vector` implements this.

## Checking

    python tests/check_loops.py rasscf.h5 [--no-pairs]

rebuilds the active density matrix from the loop values and compares it with
`DENSITY_MATRIX`, checks that the loops enumerated at the midvertices are
exactly those between all pairs of walks of the graph, and that the common
walks of all the loops give exactly the pairs of walks connected by an E_ij,
each once and with the loop value.

    python tests/check_log.py rasscf.log rasscf.h5

decodes every CSF printed by RASSCF (`PRWF`) and compares occupation string
and coefficient, and checks that the walks counted on the graph equal `NCSF`.
