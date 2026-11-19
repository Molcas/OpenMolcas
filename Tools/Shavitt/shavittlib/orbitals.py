#***********************************************************************
# This file is part of OpenMolcas.                                     *
#                                                                      *
# OpenMolcas is free software; you can redistribute it and/or modify   *
# it under the terms of the GNU Lesser General Public License, v. 2.1. *
# OpenMolcas is distributed in the hope that it will be useful, but it *
# is provided "as is" and without any express or implied warranties.   *
# For more details see the full text of the license in the file        *
# LICENSE or in <http://www.gnu.org/licenses/>.                        *
#                                                                      *
# Copyright (C) 2026, Roland Lindh                                     *
#***********************************************************************
"""The active orbitals of the levels, shown with Pegamoid.

Pegamoid (https://gitlab.com/Jellby/Pegamoid, pip install pegamoid; GPL v3)
is an optional, separately installed program; it is not part of this tool.
Its script has no import guard, so its definitions are executed here without
the final lines that start its application. Pegamoid evaluates an orbital on a
grid (its class Orbitals reads the basis set and the MOs of rasscf.h5); VTK
draws the isosurfaces off screen for the pop-up, and Pegamoid's own window
opens on an orbital on request.
"""

import os
import shutil
import sys
import types

import numpy as np

_PEG = None  # the loaded Pegamoid definitions (module), or the exception why not
_START = ('\n# VTK & Qt not supported on Wayland', '\napp = QApplication(sys.argv)')


def find_pegamoid():
    """Path of pegamoid.py: $PEGAMOID, the one installed with this Python, or the PATH.

    The copy installed with this Python comes before the PATH: it matches its packages (an old
    Pegamoid, e.g. 2.9, uses np.math, which NumPy 2 no longer has).
    """
    cand = [os.environ.get('PEGAMOID'), os.path.join(os.path.dirname(sys.executable), 'pegamoid.py'),
            os.path.join(sys.prefix, 'bin', 'pegamoid.py'), shutil.which('pegamoid.py'), shutil.which('pegamoid')]
    for c in cand:
        if c and os.path.isfile(c):
            return c
    return None


def pegamoid():
    """The Pegamoid definitions as a module; raises RuntimeError if they cannot be loaded."""
    global _PEG
    if isinstance(_PEG, types.ModuleType):
        return _PEG
    if _PEG is not None:
        raise RuntimeError(_PEG)
    path = find_pegamoid()
    if path is None:
        _PEG = 'Pegamoid not found: pip install pegamoid, or set PEGAMOID to the path of pegamoid.py'
        raise RuntimeError(_PEG)
    with open(path, encoding='utf-8') as f:
        src = f.read()
    cut = min([src.find(s) for s in _START if s in src] or [-1])
    if cut < 0:
        _PEG = f'{path}: unknown version of Pegamoid (no start-up section found)'
        raise RuntimeError(_PEG)
    # same Qt binding as the interface
    try:
        from matplotlib.backends.qt_compat import QT_API
        os.environ.setdefault('QT_API', QT_API.lower())
    except ImportError:
        pass
    mod = types.ModuleType('pegamoid_lib')
    mod.__file__ = path
    try:
        exec(compile(src[:cut], path, 'exec'), mod.__dict__)
    except (Exception, SystemExit) as e:  # Pegamoid exits when it finds no Qt binding
        _PEG = f'Pegamoid could not be loaded from {path}: {e!r}'
        raise RuntimeError(_PEG) from None
    _PEG = mod
    return mod


# element colours (CPK-like) for the atoms of the pop-up
_COLOR = {1: (0.9, 0.9, 0.9), 6: (0.35, 0.35, 0.35), 7: (0.2, 0.3, 0.95), 8: (0.95, 0.15, 0.1),
          9: (0.6, 0.9, 0.3), 15: (1.0, 0.5, 0.0), 16: (0.95, 0.85, 0.2), 17: (0.2, 0.85, 0.2)}
_RADIUS = {1: 0.31, 6: 0.76, 7: 0.71, 8: 0.66, 9: 0.57, 15: 1.07, 16: 1.05, 17: 1.02}  # covalent, Angstrom
BOHR = 0.52917721


class ActiveOrbitals:
    """The active MOs of the levels of a rasscf.h5 file, their data and images."""

    def __init__(self, h5file, graph):
        import h5py
        self.file = h5file
        self.g = graph
        with h5py.File(h5file, 'r') as f:
            types_ = [t.decode().strip() for t in f['MO_TYPEINDICES'][()]]
            nbas = [int(n) for n in f.attrs['NBAS']]
            occ = np.asarray(f['MO_OCCUPATIONS'][()])
            ene = np.asarray(f['MO_ENERGIES'][()]) if 'MO_ENERGIES' in f else np.full(len(occ), np.nan)
        # the active MOs in the order of the active orbital indices (irrep by irrep, HDF5 order)
        self.act = []
        o = 0
        for s, n in enumerate(nbas):
            for k in range(n):
                t = types_[o + k]
                if t in ('1', '2', '3'):
                    self.act.append(dict(mo=o + k, sym=s, num=k + 1, ras=int(t), occ=float(occ[o + k]),
                                         ene=float(ene[o + k])))
            o += n
        self._orb = None
        self._images = {}

    def info(self, level):
        """Data of the active orbital of a level (1-based)."""
        a = int(self.g.level_to_act[level - 1])
        d = dict(self.act[a - 1])
        d['act'] = a
        d['level'] = level
        return d

    def caption(self, level):
        d = self.info(level)
        irr = self.g.irrep_name(d['sym'])
        e = '' if np.isnan(d['ene']) or d['ene'] == 0.0 else f', energy {d["ene"]:.4f}'
        return (f'level {level}: active orbital {d["act"]} = MO {d["num"]} {irr} (RAS{d["ras"]}), '
                f'occupation {d["occ"]:.4f}{e}')

    # ------------------------------------------------------------------
    def orbitals(self):
        if self._orb is None:
            peg = pegamoid()
            self._orb = peg.Orbitals(self.file, 'hdf5')
        return self._orb

    def image(self, level, size=320, npts=48, iso=0.05):
        """PNG bytes of the isosurfaces (+iso blue, -iso red) of the active orbital of a level, cached."""
        key = (level, size, npts, iso)
        if key not in self._images:
            self._images[key] = self._render(self.info(level)['mo'], size, npts, iso)
        return self._images[key]

    def _render(self, n, size, npts, iso):
        import vtk
        from vtk.util import numpy_support
        orb = self.orbitals()
        xyz = np.array([c['xyz'] for c in orb.centers], dtype=float)
        lo, hi = xyz.min(axis=0) - 5.0, xyz.max(axis=0) + 5.0
        ax = [np.linspace(lo[k], hi[k], npts) for k in range(3)]
        z, y, x = np.meshgrid(ax[2], ax[1], ax[0], indexing='ij')  # x fastest, as vtkImageData
        val = orb.mo(n, x.ravel(), y.ravel(), z.ravel())
        img = vtk.vtkImageData()
        img.SetDimensions(npts, npts, npts)
        img.SetOrigin(*lo)
        img.SetSpacing(*[(hi[k] - lo[k]) / (npts - 1) for k in range(3)])
        arr = numpy_support.numpy_to_vtk(np.ascontiguousarray(val, dtype=float), deep=True)
        img.GetPointData().SetScalars(arr)

        ren = vtk.vtkRenderer()
        ren.SetBackground(1.0, 1.0, 1.0)
        for v, col in ((iso, (0.16, 0.47, 0.84)), (-iso, (0.89, 0.29, 0.28))):
            cf = vtk.vtkContourFilter()
            cf.SetInputData(img)
            cf.SetValue(0, v)
            nm = vtk.vtkPolyDataNormals()
            nm.SetInputConnection(cf.GetOutputPort())
            mp = vtk.vtkPolyDataMapper()
            mp.SetInputConnection(nm.GetOutputPort())
            mp.ScalarVisibilityOff()
            a = vtk.vtkActor()
            a.SetMapper(mp)
            a.GetProperty().SetColor(*col)
            a.GetProperty().SetOpacity(0.75)
            a.GetProperty().SetSpecular(0.3)
            ren.AddActor(a)
        # atoms and bonds
        Z = [int(round(c['Z'])) for c in orb.centers]
        for p, q in zip(xyz, Z):
            s = vtk.vtkSphereSource()
            s.SetCenter(*p)
            s.SetRadius(0.45 * _RADIUS.get(q, 0.8) / BOHR)
            s.SetThetaResolution(24)
            s.SetPhiResolution(24)
            mp = vtk.vtkPolyDataMapper()
            mp.SetInputConnection(s.GetOutputPort())
            a = vtk.vtkActor()
            a.SetMapper(mp)
            a.GetProperty().SetColor(*_COLOR.get(q, (0.8, 0.5, 0.8)))
            ren.AddActor(a)
        for i in range(len(xyz)):
            for j in range(i):
                d = np.linalg.norm(xyz[i] - xyz[j])
                if d < 1.2 * (_RADIUS.get(Z[i], 0.8) + _RADIUS.get(Z[j], 0.8)) / BOHR:
                    ln = vtk.vtkLineSource()
                    ln.SetPoint1(*xyz[i])
                    ln.SetPoint2(*xyz[j])
                    tb = vtk.vtkTubeFilter()
                    tb.SetInputConnection(ln.GetOutputPort())
                    tb.SetRadius(0.12)
                    tb.SetNumberOfSides(16)
                    mp = vtk.vtkPolyDataMapper()
                    mp.SetInputConnection(tb.GetOutputPort())
                    a = vtk.vtkActor()
                    a.SetMapper(mp)
                    a.GetProperty().SetColor(0.6, 0.6, 0.6)
                    ren.AddActor(a)
        ren.ResetCamera()
        cam = ren.GetActiveCamera()
        cam.Azimuth(35)
        cam.Elevation(25)
        ren.ResetCamera()
        cam.Zoom(1.25)

        win = vtk.vtkRenderWindow()
        win.SetOffScreenRendering(1)
        win.SetMultiSamples(8)
        win.AddRenderer(ren)
        win.SetSize(size, size)
        win.Render()
        w2i = vtk.vtkWindowToImageFilter()
        w2i.SetInput(win)
        w2i.Update()
        png = vtk.vtkPNGWriter()
        png.SetWriteToMemory(1)
        png.SetInputConnection(w2i.GetOutputPort())
        png.Write()
        data = bytes(numpy_support.vtk_to_numpy(png.GetResult()))
        win.Finalize()
        return data

    # ------------------------------------------------------------------
    def open_viewer(self, level, parent=None):
        """Open Pegamoid's own window on the file, showing the active orbital of a level. Returns the window."""
        from matplotlib.backends.qt_compat import QtCore
        peg = pegamoid()
        win = peg.MainWindow()
        win.filename = os.path.abspath(self.file)
        n = self.info(level)['mo'] + 1
        timer = QtCore.QTimer(win)

        tries = [0]

        def select():  # the file is read in a thread; select the orbital once its list is filled
            tries[0] += 1
            box = getattr(win, 'orbitalButton', None)
            k = box.findData(n) if box is not None else -1
            if k >= 0 and getattr(win, 'ready', True):
                timer.stop()
                box.setCurrentIndex(k)  # through the list, so that it shows the orbital (sets win.orbital)
            elif tries[0] > 400:  # about a minute: give up, leave Pegamoid as it is
                timer.stop()
        timer.timeout.connect(select)
        timer.start(150)
        win._shavitt_timer = timer
        win.show()
        return win
