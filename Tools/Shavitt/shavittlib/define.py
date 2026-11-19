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
"""Dialog that defines a CI space as the input of RASSCF does (electrons, spin,
symmetry, CAS/RAS/GAS orbitals and restrictions) and builds its split graph,
without a file. Qt through matplotlib's qt_compat.
"""

from matplotlib.backends.qt_compat import QtCore, QtWidgets

from .graph import ShavittGraph

Qt = QtCore.Qt
KINDS = ('CAS', 'RAS', 'GAS')
RAS_ROWS = ('RAS1', 'RAS2', 'RAS3')


class DefineGraphDialog(QtWidgets.QDialog):
    """Define a CI space; `graph` is the graph built from it after accept()."""

    def __init__(self, parent=None, graph=None):
        super().__init__(parent)
        self.setWindowTitle('New graph')
        self.graph = None
        v = QtWidgets.QVBoxLayout(self)

        form = QtWidgets.QFormLayout()
        self.electrons = QtWidgets.QSpinBox()
        self.electrons.setRange(0, 200)
        self.electrons.setToolTip('Number of active electrons')
        form.addRow('Active electrons', self.electrons)
        self.spin = QtWidgets.QSpinBox()
        self.spin.setRange(1, 101)
        self.spin.setToolTip('Spin multiplicity 2S+1')
        form.addRow('Spin multiplicity 2S+1', self.spin)
        self.nsym = QtWidgets.QComboBox()
        self.nsym.addItems(['1', '2', '4', '8'])
        self.nsym.setToolTip('Number of irreps of the point group (D2h and its subgroups)')
        form.addRow('Irreps', self.nsym)
        self.labels = QtWidgets.QLineEdit()
        self.labels.setPlaceholderText('optional, e.g. a1 b1 b2 a2')
        self.labels.setToolTip('Names of the irreps, separated by spaces')
        form.addRow('Irrep labels', self.labels)
        self.irrep = QtWidgets.QComboBox()
        self.irrep.setToolTip('Irrep of the wave function')
        form.addRow('Wave function irrep', self.irrep)
        self.kind = QtWidgets.QComboBox()
        self.kind.addItems(KINDS)
        self.kind.setToolTip('CAS: one space; RAS: RAS1/2/3 with holes and particles; '
                             'GAS: spaces with cumulative electron bounds')
        form.addRow('CI space', self.kind)
        v.addLayout(form)

        box = QtWidgets.QGroupBox('Active orbitals of each irrep')
        bl = QtWidgets.QVBoxLayout(box)
        self.table = QtWidgets.QTableWidget()
        self.table.setToolTip('Orbitals of each irrep in each space. GAS: Min and Max bound the electrons in the '
                              'spaces 1..k together; the last space holds all electrons.')
        bl.addWidget(self.table)
        h = QtWidgets.QHBoxLayout()
        self.add_space = QtWidgets.QPushButton('Add space')
        self.del_space = QtWidgets.QPushButton('Remove space')
        h.addWidget(self.add_space)
        h.addWidget(self.del_space)
        h.addStretch(1)
        bl.addLayout(h)
        self.ras_box = QtWidgets.QWidget()
        rf = QtWidgets.QHBoxLayout(self.ras_box)
        rf.setContentsMargins(0, 0, 0, 0)
        self.holes = QtWidgets.QSpinBox()
        self.holes.setRange(0, 400)
        self.holes.setToolTip('Maximum number of holes in RAS1')
        self.particles = QtWidgets.QSpinBox()
        self.particles.setRange(0, 400)
        self.particles.setToolTip('Maximum number of electrons in RAS3')
        rf.addWidget(QtWidgets.QLabel('Holes in RAS1'))
        rf.addWidget(self.holes)
        rf.addWidget(QtWidgets.QLabel('Electrons in RAS3'))
        rf.addWidget(self.particles)
        rf.addStretch(1)
        bl.addWidget(self.ras_box)
        v.addWidget(box)

        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        v.addWidget(self.status)
        self.buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        self.buttons.button(QtWidgets.QDialogButtonBox.Ok).setText('Build')
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        v.addWidget(self.buttons)

        self._fill(graph)
        for w in (self.electrons, self.spin, self.holes, self.particles):
            w.valueChanged.connect(self.update_status)
        self.labels.textChanged.connect(self._labels_changed)
        self.irrep.currentIndexChanged.connect(self.update_status)
        self.nsym.currentIndexChanged.connect(self._shape_changed)
        self.kind.currentIndexChanged.connect(self._shape_changed)
        self.table.itemChanged.connect(self.update_status)
        self.add_space.clicked.connect(lambda: self._resize_gas(+1))
        self.del_space.clicked.connect(lambda: self._resize_gas(-1))
        self.resize(560, 520)
        self.update_status()

    # ------------------------------------------------------------------ set up
    def _fill(self, g):
        """Start from the definition of graph g, or from a CAS(8,7) of H2O in C2v."""
        if g is None:
            n, spin, nsym, lsym, labels, kind = 8, 1, 4, 0, ['a1', 'b1', 'b2', 'a2'], 'CAS'
            rows, mins, maxs, holes, particles = [[3, 1, 2, 1]], None, None, 0, 0
        else:
            n, spin, nsym, lsym, labels, kind = g.nactel, g.ispin, g.nsym, g.lsym, list(g.irrep_labels), g.ci_type
            nras = [list(map(int, r)) for r in g.nras.T] if g.nras is not None else [[0] * g.nsym]
            mins = [int(x) for x in g.nrasel] if g.nrasel is not None else None
            maxs = [int(x) for x in g.nrasmx] if g.nrasmx is not None else None
            holes = particles = 0
            if kind == 'RAS' and len(nras) == 3:
                rows = nras
                holes, particles = 2 * sum(nras[0]) - mins[0], n - mins[1]
            elif kind == 'GAS':
                rows = nras
                maxs = [n if m < 0 else m for m in maxs]
            else:
                kind, rows = 'CAS', [[sum(c) for c in zip(*nras)]]
        self.electrons.setValue(n)
        self.spin.setValue(spin)
        self.nsym.setCurrentText(str(nsym))
        self.labels.setText(' '.join(labels))
        self.kind.setCurrentText(kind)
        self.holes.setValue(holes)
        self.particles.setValue(particles)
        self._layout_table(rows, mins, maxs)
        self._fill_irreps(lsym)

    def _irrep_names(self):
        nsym = int(self.nsym.currentText())
        labels = self.labels.text().split()
        return labels if len(labels) == nsym else [str(i + 1) for i in range(nsym)]

    def _fill_irreps(self, select=None):
        cur = self.irrep.currentIndex() if select is None else select
        self.irrep.blockSignals(True)
        self.irrep.clear()
        self.irrep.addItems(self._irrep_names())
        self.irrep.setCurrentIndex(min(max(cur, 0), self.irrep.count() - 1))
        self.irrep.blockSignals(False)

    def _labels_changed(self):
        self._fill_irreps()
        self.table.setHorizontalHeaderLabels(self._columns())
        self.update_status()

    def _columns(self):
        cols = self._irrep_names()
        return cols + ['Min', 'Max'] if self.kind.currentText() == 'GAS' else cols

    def _layout_table(self, rows, mins=None, maxs=None):
        """Rows of orbitals per irrep (and, for GAS, the cumulative bounds) for the current kind and irreps."""
        kind = self.kind.currentText()
        nsym = int(self.nsym.currentText())
        n = self.electrons.value()
        rows = [list(r[:nsym]) + [0] * (nsym - len(r[:nsym])) for r in (rows or [[0] * nsym])]
        if kind == 'CAS':
            rows, names = [[sum(c) for c in zip(*rows)]], ['CAS']
        elif kind == 'RAS':
            rows = (rows + [[0] * nsym] * 3)[:3] if len(rows) == 3 else [[0] * nsym, rows[0], [0] * nsym]
            names = list(RAS_ROWS)
        else:
            names = [f'GAS{k + 1}' for k in range(len(rows))]
            if mins is None or len(mins) != len(rows):
                mins = [0] * (len(rows) - 1) + [n]
                maxs = [n] * len(rows)
        self.table.blockSignals(True)
        self.table.clear()
        cols = self._columns()
        self.table.setColumnCount(len(cols))
        self.table.setRowCount(len(rows))
        self.table.setHorizontalHeaderLabels(cols)
        self.table.setVerticalHeaderLabels(names)
        for i, r in enumerate(rows):
            vals = list(r) + ([mins[i], maxs[i]] if kind == 'GAS' else [])
            for j, x in enumerate(vals):
                it = QtWidgets.QTableWidgetItem(str(x))
                it.setTextAlignment(Qt.AlignCenter)
                self.table.setItem(i, j, it)
        self.table.resizeColumnsToContents()
        self.table.blockSignals(False)
        self.ras_box.setVisible(kind == 'RAS')
        self.add_space.setVisible(kind == 'GAS')
        self.del_space.setVisible(kind == 'GAS')

    def _cells(self):
        """The table as integers: rows of orbitals, and for GAS the bounds (ValueError for a bad entry)."""
        nsym = int(self.nsym.currentText())
        rows, mins, maxs = [], [], []
        for i in range(self.table.rowCount()):
            vals = []
            for j in range(self.table.columnCount()):
                it = self.table.item(i, j)
                t = it.text().strip() if it else ''
                try:
                    vals.append(int(t or 0))
                except ValueError:
                    raise ValueError(f'{self.table.verticalHeaderItem(i).text()}, '
                                     f'{self.table.horizontalHeaderItem(j).text()}: "{t}" is not a whole number') from None
            rows.append(vals[:nsym])
            if len(vals) > nsym:
                mins.append(vals[nsym])
                maxs.append(vals[nsym + 1])
        return rows, mins, maxs

    def _shape_changed(self):
        try:
            rows, mins, maxs = self._cells()
        except ValueError:
            rows, mins, maxs = None, None, None
        self._layout_table(rows, mins or None, maxs or None)
        self._fill_irreps()
        self.update_status()

    def _resize_gas(self, delta):
        try:
            rows, mins, maxs = self._cells()
        except ValueError as e:
            self.status.setText(str(e))
            return
        nsym = int(self.nsym.currentText())
        n = self.electrons.value()
        if delta > 0:
            rows.append([0] * nsym)
            mins.append(n)
            maxs.append(n)
        elif len(rows) > 1:
            rows.pop()
            mins.pop()
            maxs.pop()
        self._layout_table(rows, mins, maxs)
        self.update_status()

    # ------------------------------------------------------------------ build
    def definition(self):
        """Keyword arguments of ShavittGraph.from_definition for the current entries."""
        rows, mins, maxs = self._cells()
        kind = self.kind.currentText()
        d = dict(electrons=self.electrons.value(), spin=self.spin.value(), irrep=self.irrep.currentIndex() + 1,
                 irrep_labels=self.labels.text().split() or None)
        if kind == 'CAS':
            d['cas'] = rows[0]
        elif kind == 'RAS':
            d.update(ras1=rows[0], ras2=rows[1], ras3=rows[2], holes=self.holes.value(),
                     particles=self.particles.value())
        else:
            d['gas'] = list(zip(rows, mins, maxs))
        return d

    def build(self):
        d = self.definition()
        if not any(any(r) for k in ('cas', 'ras1', 'ras2', 'ras3') for r in [d.get(k) or []]) and \
                not any(any(o) for o, _, _ in d.get('gas') or []):
            raise ValueError('no active orbitals')
        return ShavittGraph.from_definition(**d)

    def update_status(self, *_):
        ok = self.buttons.button(QtWidgets.QDialogButtonBox.Ok)
        try:
            g = self.build()
        except ValueError as e:
            self.status.setText(f'<span style="color:#c0392b">{e}</span>')
            ok.setEnabled(False)
            return
        ncsf = int(g.ncsf[g.lsym])
        self.status.setText(f'{g.ci_type}: {g.nactel} electrons in {g.nlev} orbitals, {g.nvert} vertices, '
                            f'{ncsf} CSF{"s" if ncsf != 1 else ""} of irrep {g.irrep_name(g.lsym)}<br>'
                            f'CSFs per irrep: {g.csf_counts()}')
        ok.setEnabled(ncsf > 0)

    def accept(self):
        try:
            self.graph = self.build()
        except ValueError as e:
            self.status.setText(f'<span style="color:#c0392b">{e}</span>')
            return
        super().accept()
