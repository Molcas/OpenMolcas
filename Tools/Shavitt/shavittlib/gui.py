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
"""Graphical interface to the Shavitt graph of a RASSCF wave function.

Qt through matplotlib's qt_compat: any of PySide6, PyQt6, PySide2, PyQt5.
"""

import os

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.backends.qt_compat import QtCore, QtGui, QtWidgets
from matplotlib.figure import Figure

from .draw import MAX_WALKS, SERIES, SURFACE, DrawOptions, LoopOverlay, SegmentOverlay, Walk, draw
from . import segments as sgm
from . import loops as lps
from .orbitals import ActiveOrbitals
from . import __version__
from .graph import ShavittGraph

Qt = QtCore.Qt
QAction = getattr(QtGui, 'QAction', None) or QtWidgets.QAction  # Qt6 / Qt5

# the display options: attribute of DrawOptions, label, tooltip
DISPLAY = [
    ('vertex_numbers', 'Number the vertices', 'Vertex numbers as in the distinct-row table (1-based)'),
    ('paldus_labels', 'Paldus (a,b,c) labels', 'The Paldus numbers (a,b,c) of each vertex'),
    ('midlevel_line', 'Show the midlevel', 'Dashed line at the midlevel of the split graph'),
    ('mark_midvertices', 'Mark the midvertices', 'Draw the vertices at the midlevel as squares'),
    ('step_labels', 'Label arcs with step d', 'The step d = 0, 1, 2, 3 at the middle of each arc'),
    ('level_orbitals', 'Active orbital of each level', 'Add the active orbital index to the level labels'),
    ('dim_graph', 'Fade the graph behind walks', 'Draw the unselected arcs lighter when walks are highlighted'),
    ('width_by_weight', 'Line width by weight', 'Line width of each walk proportional to |c|'),
    ('show_subgraph', 'Subgraph of the selected vertex', 'Darken all walks through the vertex selected by a click'),
    ('show_logo', 'OpenMolcas logo', 'The OpenMolcas logo in the lower right corner'),
]


def _swatch(color):
    pm = QtGui.QPixmap(14, 14)
    pm.fill(QtGui.QColor(color))
    return QtGui.QIcon(pm)


class OrbitalPopup(QtWidgets.QFrame):
    """Frameless pop-up with the picture and the data of an active orbital."""

    def __init__(self, parent=None):
        super().__init__(parent, Qt.ToolTip | Qt.FramelessWindowHint)
        self.setFrameShape(QtWidgets.QFrame.Box)
        self.setStyleSheet('OrbitalPopup { background: #ffffff; border: 1px solid #5f5e58; }')
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        self.picture = QtWidgets.QLabel()
        self.picture.setAlignment(Qt.AlignCenter)
        self.text = QtWidgets.QLabel()
        self.text.setWordWrap(True)
        self.text.setMaximumWidth(340)
        lay.addWidget(self.picture)
        lay.addWidget(self.text)

    def show_orbital(self, pos, caption, png=None, note=''):
        if png:
            pm = QtGui.QPixmap()
            pm.loadFromData(png, 'PNG')
            self.picture.setPixmap(pm)
            self.picture.show()
        else:
            self.picture.clear()
            self.picture.hide()
        self.text.setText(f'<qt><b>{caption}</b>' + (f'<br>{note}' if note else '') + '</qt>')
        self.adjustSize()
        # keep it on the screen
        scr = QtGui.QGuiApplication.screenAt(pos) or QtGui.QGuiApplication.primaryScreen()
        r = scr.availableGeometry()
        x = min(pos.x() + 18, r.right() - self.width())
        y = min(pos.y() + 18, r.bottom() - self.height())
        self.move(max(r.left(), x), max(r.top(), y))
        self.show()
        self.raise_()


class ShavittWindow(QtWidgets.QMainWindow):

    def __init__(self, filename=None, root=1, threshold=0.05, opts=None, graph=None):
        super().__init__()
        self.setWindowTitle(f'Shavitt {__version__}')
        self.g = None
        self.filename = None
        self.opts = opts if opts is not None else DrawOptions()
        # highlighted walks: key (CSF index, or step string for walks of other irreps) -> Walk
        self.walks = {}
        self.from_table = set()  # keys that follow the selection of the CSF table
        self._sync = False
        self._reset_view = True

        self.fig = Figure(facecolor=SURFACE)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('motion_notify_event', self.on_hover)
        self.canvas.mpl_connect('figure_leave_event', lambda ev: self._leave_figure())
        self.canvas.mpl_connect('resize_event', lambda ev: self.refresh())
        self._hover = None  # annotation of the vertex under the mouse
        self._hover_v = None
        # orbital pop-ups (Pegamoid): the level under the mouse, a short delay before drawing
        self.orbs = None
        self.orb_popup = OrbitalPopup(self)
        self._orb_level = None
        self._orb_timer = QtCore.QTimer(self)
        self._orb_timer.setSingleShot(True)
        self._orb_timer.timeout.connect(self.show_orbital_popup)
        self._peg_windows = []
        self.toolbar = NavigationToolbar2QT(self.canvas, self)

        right = QtWidgets.QWidget()
        rl = QtWidgets.QVBoxLayout(right)
        rl.setContentsMargins(0, 0, 0, 0)
        rl.addWidget(self.toolbar)
        rl.addWidget(self.canvas)

        panel = self._build_panel(root, threshold)
        panel.setMinimumWidth(400)

        split = QtWidgets.QSplitter()
        split.addWidget(panel)
        split.addWidget(right)
        split.setStretchFactor(1, 1)
        self.setCentralWidget(split)
        self._build_menu()
        self.statusBar().showMessage('Open a rasscf.h5 file (File > Open)')
        self.resize(1500, 900)

        if filename or graph is not None:
            self.load(filename, graph)

    # ------------------------------------------------------------------ widgets
    def _build_menu(self):
        m = self.menuBar().addMenu('&File')
        for text, key, slot in [('&New graph...', QtGui.QKeySequence.New, self.new_dialog),
                                ('&Open...', QtGui.QKeySequence.Open, self.open_dialog),
                                ('&Save figure...', QtGui.QKeySequence.Save, self.save_dialog),
                                ('&Quit', QtGui.QKeySequence.Quit, self.close)]:
            a = QAction(text, self)
            a.setShortcut(key)
            a.triggered.connect(slot)
            m.addAction(a)

    def _build_panel(self, root, threshold):
        # three tabs, each scrollable: wave function and walks, loops, view
        tabs = QtWidgets.QTabWidget()
        pages = []
        for name in ('Wave function', 'Loops', 'Segments', 'View'):
            page = QtWidgets.QWidget()
            pl = QtWidgets.QVBoxLayout(page)
            sc = QtWidgets.QScrollArea()
            sc.setWidget(page)
            sc.setWidgetResizable(True)
            sc.setFrameShape(QtWidgets.QFrame.NoFrame)
            tabs.addTab(sc, name)
            pages.append(pl)
        lay = pages[0]

        # wave function
        box = QtWidgets.QGroupBox('Wave function')
        f = QtWidgets.QFormLayout(box)
        self.file_label = QtWidgets.QLabel('-')
        self.file_label.setWordWrap(True)
        openb = QtWidgets.QPushButton('Open...')
        openb.setToolTip('Open a rasscf.h5 file')
        openb.clicked.connect(self.open_dialog)
        newb = QtWidgets.QPushButton('New...')
        newb.setToolTip('Build a graph from the definition of a CI space (electrons, spin, symmetry, CAS/RAS/GAS), '
                        'without a file')
        newb.clicked.connect(self.new_dialog)
        h = QtWidgets.QHBoxLayout()
        h.addWidget(openb)
        h.addWidget(newb)
        f.addRow(h)
        f.addRow(self.file_label)
        self.info_label = QtWidgets.QLabel('')
        self.info_label.setWordWrap(True)
        f.addRow(self.info_label)
        self.root_spin = QtWidgets.QSpinBox()
        self.root_spin.setRange(1, 1)
        self.root_spin.setValue(root)
        self.root_spin.valueChanged.connect(self.fill_table)
        f.addRow('Root', self.root_spin)
        self.energy_label = QtWidgets.QLabel('')
        f.addRow('Energy', self.energy_label)
        self.thr_spin = QtWidgets.QDoubleSpinBox()
        self.thr_spin.setDecimals(3)
        self.thr_spin.setRange(0.0, 1.0)
        self.thr_spin.setSingleStep(0.01)
        self.thr_spin.setValue(threshold)
        self.thr_spin.valueChanged.connect(self.fill_table)
        f.addRow('Threshold |c|', self.thr_spin)
        lay.addWidget(box)

        # CSF table
        box = QtWidgets.QGroupBox('Leading CSFs (select to highlight)')
        v = QtWidgets.QVBoxLayout(box)
        self.table = QtWidgets.QTableWidget(0, 4)
        self.table.setHorizontalHeaderLabels(['CSF', 'step vector', 'coeff', 'weight'])
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont))
        self.table.setMinimumHeight(180)
        self.table.itemSelectionChanged.connect(self.on_table_selection)
        v.addWidget(self.table)
        h = QtWidgets.QHBoxLayout()
        self.top_spin = QtWidgets.QSpinBox()
        self.top_spin.setRange(1, MAX_WALKS)
        self.top_spin.setValue(3)
        topb = QtWidgets.QPushButton('Select top')
        topb.clicked.connect(self.select_top)
        h.addWidget(topb)
        h.addWidget(self.top_spin)
        h.addStretch(1)
        v.addLayout(h)
        lay.addWidget(box)

        # highlight a walk
        box = QtWidgets.QGroupBox('Highlight a walk')
        f = QtWidgets.QFormLayout(box)
        self.csf_edit = QtWidgets.QLineEdit()
        self.csf_edit.setPlaceholderText('CSF number, e.g. 18')
        self.csf_edit.returnPressed.connect(self.add_csf)
        b = QtWidgets.QPushButton('Add')
        b.clicked.connect(self.add_csf)
        h = QtWidgets.QHBoxLayout()
        h.addWidget(self.csf_edit)
        h.addWidget(b)
        f.addRow('CSF', h)
        self.step_edit = QtWidgets.QLineEdit()
        self.step_edit.setPlaceholderText('0/u/d/2 for levels 1..n')
        self.step_edit.setToolTip('Step vector, level 1 first, e.g. 2202200; blanks are ignored.\n'
                                  'Any walk of the graph is accepted, also of another irrep.')
        self.step_edit.returnPressed.connect(self.add_steps)
        b = QtWidgets.QPushButton('Add')
        b.clicked.connect(self.add_steps)
        h = QtWidgets.QHBoxLayout()
        h.addWidget(self.step_edit)
        h.addWidget(b)
        f.addRow('Step vector', h)
        self.walk_list = QtWidgets.QListWidget()
        self.walk_list.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont))
        self.walk_list.setMinimumHeight(110)
        f.addRow(self.walk_list)
        h = QtWidgets.QHBoxLayout()
        b = QtWidgets.QPushButton('Remove')
        b.clicked.connect(self.remove_walk)
        h.addWidget(b)
        b = QtWidgets.QPushButton('Clear all')
        b.clicked.connect(self.clear_walks)
        h.addWidget(b)
        f.addRow(h)
        lay.addWidget(box)

        # coupling-coefficient loops
        lay = pages[1]
        box = QtWidgets.QGroupBox('Coupling-coefficient loops')
        box.setToolTip('One midvertex: the loops of E_ij in the upper half (head to midvertex) and in the lower\n'
                       'half (midvertex to tail). Two partner midvertices: the upper and lower partial loops\n'
                       'of the loops crossing the midlevel there; their product is the coupling coefficient.')
        f = QtWidgets.QFormLayout(box)
        self.mv1_combo = QtWidgets.QComboBox()
        self.mv1_combo.currentIndexChanged.connect(self.on_mv1_changed)
        f.addRow('Midvertex', self.mv1_combo)
        self.mv2_combo = QtWidgets.QComboBox()
        self.mv2_combo.currentIndexChanged.connect(self.recompute_loops)
        f.addRow('Partner', self.mv2_combo)
        self.pick_check = QtWidgets.QCheckBox('Pick midvertices with the mouse')
        self.pick_check.setToolTip('Click a midvertex to select it; click a partner of it for the crossing loops;\n'
                                   'click the selected midvertex again to clear (on this tab only)')
        self.pick_check.setChecked(True)
        f.addRow(self.pick_check)
        h = QtWidgets.QHBoxLayout()
        self.head_combo = QtWidgets.QComboBox()
        self.tail_combo = QtWidgets.QComboBox()
        self.head_combo.currentIndexChanged.connect(self.recompute_loops)
        self.tail_combo.currentIndexChanged.connect(self.recompute_loops)
        h.addWidget(QtWidgets.QLabel('j'))
        h.addWidget(self.head_combo)
        h.addWidget(QtWidgets.QLabel('i'))
        h.addWidget(self.tail_combo)
        f.addRow('Levels', h)
        h = QtWidgets.QHBoxLayout()
        self.upper_check = QtWidgets.QCheckBox('upper')
        self.lower_check = QtWidgets.QCheckBox('lower')
        self.sym_check = QtWidgets.QCheckBox('irrep i = irrep j')
        self.sym_check.setToolTip('Only the loops of operators E_ij that conserve the symmetry (complete loops)')
        for c in (self.upper_check, self.lower_check):
            c.setChecked(True)
        for c in (self.upper_check, self.lower_check, self.sym_check):
            c.toggled.connect(self.recompute_loops)
            h.addWidget(c)
        f.addRow('Show', h)
        self.loop_info = QtWidgets.QLabel('click a midvertex in the graph (or choose one above)')
        self.loop_info.setWordWrap(True)
        f.addRow(self.loop_info)
        self.loop_list = QtWidgets.QListWidget()
        self.loop_list.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont))
        self.loop_list.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.loop_list.setMinimumHeight(380)
        self.loop_list.setToolTip(
            '<qt><b>Loop list</b> (hover a row for its explanation)<br>'
            '<b>U</b>, <b>L</b>: complete loop of E<sub>ij</sub> in the upper / lower half of the split graph.<br>'
            '<b>U*</b> E(?,j): upper partial loop of a loop crossing the midlevel; its tail level i is in the '
            'lower half, hence <b>?</b>.<br>'
            '<b>L*</b> E(i,?): lower partial loop; its head level j is in the upper half, hence <b>?</b>.<br>'
            'A U* and an L* of the same midvertex pair form the loop of E<sub>ij</sub>; the coupling '
            'coefficient is the product of their values.<br>'
            'Columns: tag, operator, value, bra L and ket R steps (lowest level first).<br>'
            'Nothing selected: all listed loops are drawn. Select loops to see them in detail.</qt>')
        self.loop_list.itemSelectionChanged.connect(self.on_loop_selection)
        f.addRow(self.loop_list)
        h = QtWidgets.QHBoxLayout()
        self.cw_spin = QtWidgets.QSpinBox()
        self.cw_spin.setRange(0, 0)
        self.cw_spin.setSpecialValueText('all')
        self.cw_spin.setEnabled(False)
        self.cw_spin.setToolTip(
            '<qt>Outside the loop bra and ket share their walk (dark grey). Usually several such common walks '
            'exist: each one embeds the loop in another pair of CSFs, all with the same coupling coefficient.<br>'
            'Step through them here (<b>all</b>: draw them all at once). Available for one selected loop, or '
            'one U* together with one L*.</qt>')
        self.cw_spin.valueChanged.connect(self.on_common_walk)
        h.addWidget(self.cw_spin)
        self.cw_label = QtWidgets.QLabel('')
        h.addWidget(self.cw_label, 1)
        f.addRow('Common walk', h)
        b = QtWidgets.QPushButton('Clear loops')
        b.clicked.connect(lambda: self.mv1_combo.setCurrentIndex(0))
        f.addRow(b)
        lay.addWidget(box)
        self.loop_items = []
        self.loop_note = ''
        self.loop_mvs = []
        self.loop_mv = None
        self.common = None  # CommonWalks of the selected loop(s), when they can be stepped through

        # segments of a vertex
        lay = pages[2]
        box = QtWidgets.QGroupBox('Segments of a vertex')
        box.setToolTip('The one-level pieces of the loops (SGUGA MKSEG): a pair of top vertices (bra L, ket R),\n'
                       'the bra and ket steps, the pair of bottom vertices, and the segment value.')
        f = QtWidgets.QFormLayout(box)
        h = QtWidgets.QHBoxLayout()
        self.seg_spin = QtWidgets.QSpinBox()
        self.seg_spin.setRange(0, 0)
        self.seg_spin.setSpecialValueText('none')
        self.seg_spin.setToolTip('Vertex number (1-based, as in the DRT); or click a vertex in the graph')
        self.seg_spin.valueChanged.connect(lambda n: self.set_segment_vertex(n - 1 if n > 0 else None))
        h.addWidget(self.seg_spin)
        h.addWidget(QtWidgets.QLabel('(click a vertex in the graph)'))
        h.addStretch(1)
        f.addRow('Vertex', h)
        self.seg_vinfo = QtWidgets.QLabel('')
        self.seg_vinfo.setTextFormat(Qt.RichText)
        self.seg_vinfo.setWordWrap(True)
        f.addRow(self.seg_vinfo)
        h = QtWidgets.QHBoxLayout()
        self.seg_kind = {}
        for key, text, tip in [('walk', 'walk', 'walk segments above/below a loop (1-4, 23-26)'),
                               ('head', 'head', 'loop heads (5-8)'), ('mid', 'mid', 'intermediate segments (9-18)'),
                               ('tail', 'tail', 'loop tails (19-22)'), ('weight', 'weight', 'E_ii weights (27-34)')]:
            c = QtWidgets.QCheckBox(text)
            c.setToolTip(tip)
            c.setChecked(key != 'walk')
            c.toggled.connect(self.recompute_segments)
            h.addWidget(c)
            self.seg_kind[key] = c
        f.addRow('Types', h)
        h = QtWidgets.QHBoxLayout()
        self.seg_left = QtWidgets.QCheckBox('as bra (L) vertex')
        self.seg_right = QtWidgets.QCheckBox('as ket (R) vertex')
        self.seg_left.setToolTip('Segments whose top-left (bra) vertex is the vertex: the SGUGA tables ISGM/VSGM')
        self.seg_right.setToolTip('Segments whose top-right (ket) vertex is the vertex, with a partner as bra')
        for c in (self.seg_left, self.seg_right):
            c.setChecked(True)
            c.toggled.connect(self.recompute_segments)
            h.addWidget(c)
        f.addRow('Role', h)
        self.seg_info = QtWidgets.QLabel('no vertex selected')
        self.seg_info.setWordWrap(True)
        f.addRow(self.seg_info)
        self.seg_list = QtWidgets.QListWidget()
        self.seg_list.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont))
        self.seg_list.setSelectionMode(QtWidgets.QAbstractItemView.ExtendedSelection)
        self.seg_list.setMinimumHeight(380)
        self.seg_list.setToolTip(
            '<qt><b>Segment list</b> (hover a row for its explanation)<br>'
            'Columns: SGUGA segment number, type, role of the vertex (L: bra / R: ket top vertex), '
            'bra and ket step (0 u d 2), value, and the vertices top L,R &rarr; bottom L,R.<br>'
            'Nothing selected: all listed segments are drawn. Select segments to see them in detail.</qt>')
        self.seg_list.itemSelectionChanged.connect(self.on_segment_selection)
        f.addRow(self.seg_list)
        lay.addWidget(box)
        self.seg_vertex = None
        self.seg_items = []
        self.seg_note = ''

        # display options
        lay = pages[3]
        box = QtWidgets.QGroupBox('Display')
        v = QtWidgets.QVBoxLayout(box)
        self.checks = {}
        for name, text, tip in DISPLAY:
            c = QtWidgets.QCheckBox(text)
            c.setToolTip(tip)
            c.setChecked(bool(getattr(self.opts, name)))
            c.toggled.connect(lambda on, n=name: self.set_option(n, on))
            v.addWidget(c)
            self.checks[name] = c
        self.orb_check = QtWidgets.QCheckBox('Orbital pop-ups (Pegamoid)')
        self.orb_check.setChecked(True)
        self.orb_check.setToolTip('Hover a level label for a picture of its active orbital; click the label to\n'
                                  'open the orbital in Pegamoid (pip install pegamoid)')
        v.addWidget(self.orb_check)
        self.logo_button = QtWidgets.QPushButton('Logo color...')
        self.logo_button.setToolTip('Color of the OpenMolcas logo (default black)')
        self.logo_button.clicked.connect(self.choose_logo_color)
        self.logo_button.setIcon(_swatch(self.opts.logo_color))
        v.addWidget(self.logo_button)
        b = QtWidgets.QPushButton('Reset view')
        b.clicked.connect(self.reset_view)
        v.addWidget(b)
        lay.addWidget(box)

        # vertex
        box = QtWidgets.QGroupBox('Vertex (click one in the graph)')
        v = QtWidgets.QVBoxLayout(box)
        self.vertex_label = QtWidgets.QLabel('none selected')
        self.vertex_label.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont))
        self.vertex_label.setTextInteractionFlags(Qt.TextSelectableByMouse)
        v.addWidget(self.vertex_label)
        lay.addWidget(box)

        # export
        box = QtWidgets.QGroupBox('Export')
        h = QtWidgets.QHBoxLayout(box)
        b = QtWidgets.QPushButton('Save PNG / PDF / SVG...')
        b.clicked.connect(self.save_dialog)
        h.addWidget(b)
        h.addWidget(QtWidgets.QLabel('dpi'))
        self.dpi_spin = QtWidgets.QSpinBox()
        self.dpi_spin.setRange(72, 1200)
        self.dpi_spin.setValue(300)
        h.addWidget(self.dpi_spin)
        lay.addWidget(box)

        for pl in pages:
            pl.addStretch(1)
        self.tabs = tabs
        tabs.currentChanged.connect(lambda i: self.refresh())
        return tabs

    # ------------------------------------------------------------------ file
    def open_dialog(self):
        fn, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open RASSCF wave function', '',
                                                      'HDF5 files (*.h5 *.hdf5);;All files (*)')
        if fn:
            self.load(fn)

    def new_dialog(self):
        """Define a CI space (electrons, spin, symmetry, CAS/RAS/GAS) and show its graph, without a file."""
        from .define import DefineGraphDialog
        dlg = DefineGraphDialog(self, self.g)
        if dlg.exec() and dlg.graph is not None:
            self.load(None, dlg.graph)

    def load(self, filename, graph=None):
        """Show the graph of a rasscf.h5 file, or a graph built from scratch (filename None)."""
        if graph is None:
            try:
                g = ShavittGraph.from_h5(filename)
            except Exception as e:  # show any reading problem to the user
                QtWidgets.QMessageBox.critical(self, 'Shavitt', f'Could not read {filename}:\n{e}')
                return
        else:
            g = graph
        self.g = g
        self.filename = filename
        self.orb_popup.hide()
        self._orb_level = None
        self.orbs = None
        if filename is not None:
            try:
                self.orbs = ActiveOrbitals(filename, g)
            except Exception:  # no MO data in the file
                self.orbs = None
        name = os.path.basename(filename) if filename else 'Graph from its definition (no wave function)'
        self.walks.clear()
        self.from_table.clear()
        self.opts.selected_vertex = None
        self._reset_view = True
        self.setWindowTitle(f'Shavitt {__version__} - {name}')
        self.file_label.setText(name)
        self.info_label.setText(
            f'{g.ci_type} {g.nactel} electrons in {g.nlev} orbitals, 2S+1 = {g.ispin}, irrep {g.irrep_name(g.lsym)}\n'
            f'{g.nvert} vertices, {len(g.arcs())} arcs, midlevel {g.midlev} ({g.nmidv} midvertices), '
            f'{int(g.ncsf[g.lsym])} CSFs')
        nroots = 0 if g.ci_vectors is None else g.ci_vectors.shape[0]
        self.root_spin.blockSignals(True)
        self.root_spin.setRange(1, max(1, nroots))
        self.root_spin.blockSignals(False)
        self.vertex_label.setText('none selected')
        self.fill_loop_controls()
        self.seg_spin.blockSignals(True)
        self.seg_spin.setRange(0, g.nvert)
        self.seg_spin.setValue(0)
        self.seg_spin.blockSignals(False)
        self.set_segment_vertex(None)
        self.fill_table()
        self.select_top()
        self.statusBar().showMessage(f'Read {filename}' if filename else 'Graph built from its definition', 5000)
        if g.warnings:
            QtWidgets.QMessageBox.warning(self, 'Shavitt', f'{os.path.basename(filename)}:\n\n' + '\n\n'.join(g.warnings))

    # ------------------------------------------------------------------ CSFs
    def fill_table(self):
        g = self.g
        if g is None:
            return
        root = self.root_spin.value() - 1
        e = ''
        if g.root_energies is not None and root < len(g.root_energies):
            e = f'{g.root_energies[root]:.8f}'
        self.energy_label.setText(e)
        # the coefficients of the highlighted CSFs follow the root
        for key, w in self.walks.items():
            if isinstance(key, int):
                w.weight, w.label = self._csf_label(key)
        rows = g.leading_csfs(root, self.thr_spin.value()) if g.ci_vectors is not None else []
        self._sync = True
        self.table.setRowCount(len(rows))
        for r, (i, c, s) in enumerate(rows):
            for col, text in enumerate([f'{i + 1}', g.csf_string(s), f'{c:+.5f}', f'{c * c:.5f}']):
                it = QtWidgets.QTableWidgetItem(text)
                it.setData(Qt.UserRole, i)
                if col != 1:
                    it.setTextAlignment(int(Qt.AlignRight | Qt.AlignVCenter))
                self.table.setItem(r, col, it)
        self.table.resizeColumnsToContents()
        self._restore_table_selection()
        self._sync = False
        self.refresh()

    def _restore_table_selection(self):
        sm = self.table.selectionModel()
        sm.clearSelection()
        for r in range(self.table.rowCount()):
            if self.table.item(r, 0).data(Qt.UserRole) in self.walks:
                sm.select(self.table.model().index(r, 0),
                          QtCore.QItemSelectionModel.Select | QtCore.QItemSelectionModel.Rows)

    def select_top(self):
        if self.g is None:
            return
        n = min(self.top_spin.value(), self.table.rowCount())
        for k in list(self.from_table):
            self.walks.pop(k, None)
        self.from_table.clear()
        for r in range(n):
            self._add_csf_key(self.table.item(r, 0).data(Qt.UserRole), from_table=True)
        self._sync = True
        self._restore_table_selection()
        self._sync = False
        self.refresh()

    def on_table_selection(self):
        if self._sync or self.g is None:
            return
        sel = {self.table.item(ix.row(), 0).data(Qt.UserRole) for ix in self.table.selectionModel().selectedRows()}
        shown = {self.table.item(r, 0).data(Qt.UserRole) for r in range(self.table.rowCount())}
        for k in list(self.walks):
            if k in shown and k not in sel:
                self.walks.pop(k)
                self.from_table.discard(k)
        for k in sorted(sel, key=lambda i: -abs(self._coef(i))):
            if k not in self.walks:
                if not self._add_csf_key(k, from_table=True):
                    break
        self._sync = True
        self._restore_table_selection()
        self._sync = False
        self.refresh()

    def _coef(self, i):
        g = self.g
        if g.ci_vectors is None:
            return 0.0
        return float(g.ci_vectors[self.root_spin.value() - 1, i])

    def _csf_label(self, i):
        g = self.g
        s = g.csf_string(g.step_vector(i))
        if g.ci_vectors is None:
            return None, f'{i + 1:>6d}  {s}'
        c = self._coef(i)
        return c * c, f'{i + 1:>6d}  {s}  {c:+.3f}'

    def _free_slot(self):
        used = {w.slot for w in self.walks.values()}
        for k in range(MAX_WALKS):
            if k not in used:
                return k
        return None

    def _add_csf_key(self, i, from_table=False):
        if i in self.walks:
            return True
        slot = self._free_slot()
        if slot is None:
            self.statusBar().showMessage(f'At most {MAX_WALKS} walks can be highlighted', 5000)
            return False
        weight, label = self._csf_label(i)
        self.walks[i] = Walk(label, self.g.step_vector(i), weight, slot)
        if from_table:
            self.from_table.add(i)
        return True

    # ------------------------------------------------------------------ walks
    def add_csf(self):
        g = self.g
        if g is None:
            return
        t = self.csf_edit.text().strip()
        try:
            i = int(t) - 1
            g.step_vector(i)
        except (ValueError, IndexError):
            self.statusBar().showMessage(f"'{t}' is not a CSF number 1..{int(g.ncsf[g.lsym])}", 6000)
            return
        if self._add_csf_key(i):
            self.csf_edit.clear()
            self._sync = True
            self._restore_table_selection()
            self._sync = False
            self.refresh()

    def add_steps(self):
        g = self.g
        if g is None:
            return
        t = self.step_edit.text()
        try:
            s = g.parse_steps(t)
            if len(s) != g.nlev:
                raise ValueError(f'{len(s)} steps given, the graph has {g.nlev} levels')
            g.walk(s)
        except ValueError as e:
            self.statusBar().showMessage(str(e), 8000)
            return
        i = g.csf_index(s)
        if i is not None:
            self.csf_edit.setText(str(i + 1))
            self.step_edit.clear()
            self.add_csf()
            self.statusBar().showMessage(f'{g.csf_string(s)} is CSF {i + 1}', 6000)
            return
        key = g.csf_string(s)
        if key in self.walks:
            return
        slot = self._free_slot()
        if slot is None:
            self.statusBar().showMessage(f'At most {MAX_WALKS} walks can be highlighted', 5000)
            return
        irr = g.irrep_name(g.walk_irrep(s))
        self.walks[key] = Walk(f'     -  {key}  ({irr})', s, None, slot)
        self.step_edit.clear()
        self.statusBar().showMessage(f'{key} is a walk of irrep {irr}, not a CSF of the wave function', 8000)
        self.refresh()

    def remove_walk(self):
        for it in self.walk_list.selectedItems():
            k = it.data(Qt.UserRole)
            self.walks.pop(k, None)
            self.from_table.discard(k)
        self._sync = True
        self._restore_table_selection()
        self._sync = False
        self.refresh()

    def clear_walks(self):
        self.walks.clear()
        self.from_table.clear()
        self._sync = True
        self.table.clearSelection()
        self._sync = False
        self.refresh()

    # ------------------------------------------------------------------ vertices
    def _vertex_at(self, ev, tol=10.0):
        """The vertex within tol pixels of the mouse event, or None."""
        if self.g is None or ev.inaxes is not self.ax or ev.x is None:
            return None
        xy = self.ax.transData.transform(np.column_stack([self._x, self._y]))
        d = np.hypot(xy[:, 0] - ev.x, xy[:, 1] - ev.y)
        v = int(np.argmin(d))
        return v if d[v] < tol else None

    def _new_hover(self):
        self._hover = self.ax.annotate('', (0, 0), xytext=(12, 12), textcoords='offset points', fontsize=8,
                                       color='#1f1f1e', zorder=40, annotation_clip=False,
                                       bbox=dict(boxstyle='round,pad=0.35', fc='#ffffff', ec='#5f5e58', lw=0.8))
        self._hover.set_in_layout(False)
        self._hover.set_visible(False)
        self._hover_v = None

    def _level_label_at(self, ev):
        """The level (1..nlev) whose label on the left axis is under the mouse, or None."""
        if self.g is None or ev.x is None:
            return None
        r = self.canvas.get_renderer()
        for L, t in enumerate(self.ax.get_yticklabels()):
            if L == 0 or not t.get_text():
                continue
            bb = t.get_window_extent(r).expanded(1.15, 1.4)
            if bb.contains(ev.x, ev.y):
                return L
        return None

    def _leave_figure(self):
        self._orb_timer.stop()
        self._orb_level = None
        self.orb_popup.hide()

    def on_level_hover(self, ev):
        L = self._level_label_at(ev) if (self.orbs is not None and self.orb_check.isChecked()) else None
        if L == self._orb_level:
            return
        self._orb_level = L
        self.orb_popup.hide()
        self._orb_timer.stop()
        if L is not None:
            self._orb_timer.start(250)
            self.statusBar().showMessage(self.orbs.caption(L) + '  -  click to open it in Pegamoid', 6000)

    def show_orbital_popup(self):
        L = self._orb_level
        if L is None or self.orbs is None:
            return
        pos = QtGui.QCursor.pos()
        caption = self.orbs.caption(L)
        try:
            QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
            png = self.orbs.image(L)
            note = 'isosurfaces &plusmn;0.05 (blue +, red &minus;) &middot; click the label for Pegamoid'
        except Exception as e:  # Pegamoid or VTK missing or failing: show the data only
            png = None
            note = f'<i>No picture: {e}</i>'
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        if self._orb_level == L:
            self.orb_popup.show_orbital(pos, caption, png, note)

    def open_pegamoid(self, L):
        self.orb_popup.hide()
        try:
            QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
            win = self.orbs.open_viewer(L)
        except Exception as e:
            QtWidgets.QMessageBox.warning(self, 'Shavitt', f'Pegamoid could not be opened:\n{e}')
            return
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self._peg_windows.append(win)
        self.statusBar().showMessage(f'Pegamoid: {self.orbs.caption(L)}', 8000)

    def on_hover(self, ev):
        self.on_level_hover(ev)
        if self._hover is None:
            return
        v = self._vertex_at(ev)
        if v == self._hover_v:
            return
        self._hover_v = v
        if v is None:
            self._hover.set_visible(False)
        else:
            L, N, a, b, c = (int(t) for t in self.g.drt[v])
            self._hover.xy = (self._x[v], self._y[v])
            # open to the left near the right edge so that the box stays in view
            right = self._x[v] > 0.7 * self._x.max()
            self._hover.set_position((-12, 12) if right else (12, 12))
            self._hover.set_horizontalalignment('right' if right else 'left')
            self._hover.set_text(f'(a,b,c) = ({a},{b},{c})\nvertex {v + 1}, level {L}, N = {N}, S = {b / 2:g}')
            self._hover.set_visible(True)
        self.canvas.draw_idle()

    def on_click(self, ev):
        if self.g is not None and ev.button == 1 and self.orbs is not None and self.orb_check.isChecked():
            L = self._level_label_at(ev)
            if L is not None:
                self.open_pegamoid(L)
                return
        if self.g is None or ev.inaxes is not self.ax or ev.button != 1:
            return
        if self.toolbar.mode:  # zoom or pan active
            return
        v = self._vertex_at(ev, tol=12.0)
        if self.segments_mode():
            self.seg_spin.setValue(0 if v is None else v + 1)  # sets the segment vertex
            return
        if (self.loops_mode() and self.pick_check.isChecked() and v is not None
                and v in self.g.midvertices()):
            self.pick_midvertex(v)
            return
        self.select_vertex(v)

    def select_vertex(self, v):
        self.opts.selected_vertex = v
        g = self.g
        if v is None:
            self.vertex_label.setText('none selected')
        else:
            L, N, a, b, c = (int(t) for t in g.drt[v])
            up, dn = g.path_counts()
            irr = g.irrep_name(g.level_sym[L - 1]) if L > 0 else '-'
            lines = [f'vertex {v + 1}: level {L} (orbital irrep {irr})',
                     f'N = {N}, (a,b,c) = ({a},{b},{c}), S = {b / 2:g}',
                     f'walks head->v: {up[v]}, v->tail: {dn[v]}',
                     f'walks through v: {up[v] * dn[v]} (all irreps)']
            if L == g.midlev:
                lines.append(f'midvertex {v - g.mvsta + 1} of {g.nmidv}')
            arcs = [f'd={d}->{int(g.down[v, d]) + 1}' for d in range(4) if g.down[v, d] >= 0]
            if arcs:
                lines.append('down: ' + ', '.join(arcs))
            self.vertex_label.setText('\n'.join(lines))
        self.refresh()

    # ------------------------------------------------------------------ loops
    def _mv_text(self, v):
        a, b, c = (int(t) for t in self.g.drt[v, 2:5])
        return f'MV {v - self.g.mvsta + 1}: vertex {v + 1} ({a},{b},{c})'

    def fill_loop_controls(self):
        g = self.g
        for cb in (self.mv1_combo, self.mv2_combo, self.head_combo, self.tail_combo):
            cb.blockSignals(True)
            cb.clear()
        self.mv1_combo.addItem('none', None)
        for v in g.midvertices():
            self.mv1_combo.addItem(self._mv_text(v), v)
        self.mv2_combo.addItem('none', None)
        for cb in (self.head_combo, self.tail_combo):
            cb.addItem('any', None)
            for L in range(g.nlev, 0, -1):
                cb.addItem(f'{L} {g.irrep_name(g.level_sym[L - 1])}', L)
        for cb in (self.mv1_combo, self.mv2_combo, self.head_combo, self.tail_combo):
            cb.blockSignals(False)
        self.recompute_loops()

    def on_mv1_changed(self):
        mv = self.mv1_combo.currentData()
        self.mv2_combo.blockSignals(True)
        self.mv2_combo.clear()
        self.mv2_combo.addItem('none', None)
        if mv is not None:
            for w in lps.midvertex_partners(self.g, mv):
                self.mv2_combo.addItem(self._mv_text(w), w)
        self.mv2_combo.blockSignals(False)
        self.recompute_loops()

    def pick_midvertex(self, v):
        mv1 = self.mv1_combo.currentData()
        if mv1 == v:
            self.mv1_combo.setCurrentIndex(0)
            return
        if mv1 is not None and self.mv2_combo.currentData() is None and self.mv2_combo.findData(v) > 0:
            self.mv2_combo.setCurrentIndex(self.mv2_combo.findData(v))
            return
        self.mv1_combo.setCurrentIndex(self.mv1_combo.findData(v))

    def _loop_text(self, lp):
        g = self.g
        if lp.kind in ('upper', 'lower'):
            lo, hi = lp.i, lp.j
        elif lp.kind == 'upper-partial':
            lo, hi = g.midlev + 1, lp.j
        else:
            lo, hi = lp.i, g.midlev
        tag = {'upper': 'U ', 'lower': 'L ', 'upper-partial': 'U*', 'lower-partial': 'L*'}[lp.kind]
        # value first so that it stays visible; the steps of levels lo..hi, lowest level first
        text = f'{tag} {lp.operator():8s} {lp.value:+.5f}  L {lp.steps_string("bra")}  R {lp.steps_string("ket")}'
        tip = self._loop_tip(lp, lo, hi)
        return text, tip

    def _loop_tip(self, lp, lo, hi):
        """Rich-text tooltip that explains a row of the loop list."""
        g = self.g
        m = g.midlev
        bra, ket = lp.steps_string('bra'), lp.steps_string('ket')
        if lp.kind in ('upper', 'lower'):
            half = ('upper half: between the graph head and the midvertex' if lp.kind == 'upper'
                    else 'lower half: between the midvertex and the graph tail')
            what = (f'<b>{"U" if lp.kind == "upper" else "L"} &nbsp;{lp.operator()}</b> &mdash; a complete loop '
                    f'in the {half}.<br>'
                    f'Head level j = {lp.j}, tail level i = {lp.i}: the coupling coefficient '
                    f'&lang;L|E<sub>{lp.i},{lp.j}</sub>|R&rang; = <b>{lp.value:+.6f}</b>.<br>'
                    f'The walk on the other side of the midvertex is a common spectator of bra and ket.')
        elif lp.kind == 'upper-partial':
            what = (f'<b>U* &nbsp;E(?,{lp.j})</b> &mdash; an <i>upper partial loop</i>: the part above the '
                    f'midlevel of a loop that crosses it.<br>'
                    f'Its head level is j = {lp.j}; the tail level i (&le; {m}) is in the lower half, so it is '
                    f'still open: <b>?</b>.<br>'
                    f'Combined with any lower partial loop L* &nbsp;E(i,?) of the same midvertex pair it gives '
                    f'the loop of E<sub>i,{lp.j}</sub>, with<br>'
                    f'&nbsp;&nbsp;&lang;L|E<sub>i,{lp.j}</sub>|R&rang; = (U* value) &times; (L* value) '
                    f'= {lp.value:+.6f} &times; (L* value).<br>'
                    f'Select one U* and one L* to see the full loop and its value.')
        else:
            what = (f'<b>L* &nbsp;E({lp.i},?)</b> &mdash; a <i>lower partial loop</i>: the part below the '
                    f'midlevel of a loop that crosses it.<br>'
                    f'Its tail level is i = {lp.i}; the head level j (&gt; {m}) is in the upper half, so it is '
                    f'still open: <b>?</b>.<br>'
                    f'Combined with any upper partial loop U* &nbsp;E(?,j) of the same midvertex pair it gives '
                    f'the loop of E<sub>{lp.i},j</sub>, with<br>'
                    f'&nbsp;&nbsp;&lang;L|E<sub>{lp.i},j</sub>|R&rang; = (U* value) &times; (L* value) '
                    f'= (U* value) &times; {lp.value:+.6f}.<br>'
                    f'Select one U* and one L* to see the full loop and its value.')
        steps = (f'<br><br>Steps of levels {lo}..{hi} (lowest level first; 0 u d 2):<br>'
                 f'<tt>&nbsp;bra L: {bra}<br>&nbsp;ket R: {ket}</tt><br>'
                 f'Bra L (solid violet) has one electron more than ket R (dashed orange) inside the loop: '
                 f'L = E<sub>ij</sub> R up to the coupling coefficient.<br>'
                 f'Value: product of the segment values (SGUGA, ASTA) = {lp.value:+.6f}.')
        return f'<qt>{what}{steps}</qt>'

    def recompute_loops(self):
        g = self.g
        self.loop_list.blockSignals(True)
        self.loop_list.clear()
        self.loop_items = []
        self.loop_mvs = []
        self.loop_mv = None
        self.loop_note = self.loop_value_note = ''
        self._set_common([])
        limit = 2000
        mv1 = self.mv1_combo.currentData() if g is not None else None
        if mv1 is None:
            self.loop_info.setText('click a midvertex in the graph (or choose one above)')
        else:
            mv2 = self.mv2_combo.currentData()
            head, tail = self.head_combo.currentData(), self.tail_combo.currentData()
            if mv2 is None:
                items = lps.midvertex_loops(g, mv1, head, tail, self.upper_check.isChecked(),
                                            self.lower_check.isChecked(), limit=limit)
                if self.sym_check.isChecked():
                    items = [lp for lp in items if g.level_sym[lp.i - 1] == g.level_sym[lp.j - 1]]
                nu = sum(lp.kind == 'upper' for lp in items)
                self.loop_info.setText(f'{nu} upper and {len(items) - nu} lower loops at '
                                       f'{self._mv_text(mv1)}' + (' (list truncated)' if len(items) >= limit else ''))
                self.loop_mvs = [mv1]
                self.loop_mv = mv1
            else:
                L, R, ups, lows = lps.crossing_loops(g, mv1, mv2, head, tail, limit=limit)
                nu, nl = len(ups), len(lows)
                # the upper/lower switches apply to the partial loops as well
                items = (ups if self.upper_check.isChecked() else []) + (lows if self.lower_check.isChecked() else [])
                hidden = [t for t, c in (('upper', self.upper_check), ('lower', self.lower_check)) if not c.isChecked()]
                self.loop_info.setText(
                    f'crossing loops at L = MV {L - g.mvsta + 1}, R = MV {R - g.mvsta + 1}: {nu} upper x '
                    f'{nl} lower partial loops = {nu * nl} loops. '
                    + (f'Hidden: the {" and ".join(hidden)} partial loops. ' if hidden else '')
                    + 'Select one U* and one L* for a coupling coefficient.')
                self.loop_mvs = [L, R]
            self.loop_items = items
            for k, lp in enumerate(items):
                text, tip = self._loop_text(lp)
                it = QtWidgets.QListWidgetItem(text)
                it.setToolTip(tip)
                it.setData(Qt.UserRole, k)
                self.loop_list.addItem(it)
        self.loop_list.blockSignals(False)
        self.refresh()

    def _selected_loops(self):
        return [self.loop_items[it.data(Qt.UserRole)] for it in self.loop_list.selectedItems()]

    def on_loop_selection(self):
        sel = self._selected_loops()
        self.loop_note = ''
        ups = [lp for lp in sel if lp.kind == 'upper-partial']
        lows = [lp for lp in sel if lp.kind == 'lower-partial']
        if len(ups) == 1 and len(lows) == 1:
            u, lo = ups[0], lows[0]
            self.loop_note = (f'<L|E({lo.i},{u.j})|R> = {u.value:+.5f} x {lo.value:+.5f} = '
                              f'{u.value * lo.value:+.5f}  (upper x lower partial loop)')
        elif len(sel) == 1 and sel[0].kind in ('upper', 'lower'):
            lp = sel[0]
            self.loop_note = f'<L|{lp.operator()}|R> = {lp.value:+.5f}  ({lp.kind} loop)'
        self.loop_value_note = self.loop_note
        self._set_common(sel if (len(sel) == 1 or (len(ups) == 1 and len(lows) == 1 and len(sel) == 2)) else [])
        self.refresh()

    def _set_common(self, sel):
        """The common walks of the selected loop(s); the spin box steps through them."""
        self.common = None
        if sel and self.g is not None:
            self.common = lps.CommonWalks(self.g, sel, self.loop_mv, upper=self.upper_check.isChecked(),
                                          lower=self.lower_check.isChecked())
        n = self.common.count if self.common is not None else 0
        self.cw_spin.blockSignals(True)
        self.cw_spin.setRange(0, min(n, 2**31 - 1))
        self.cw_spin.setValue(1 if n >= 1 else 0)
        self.cw_spin.setEnabled(n > 1)
        self.cw_spin.blockSignals(False)
        self._common_note()

    def on_common_walk(self):
        self._common_note()
        self.refresh()

    def _common_note(self):
        """Label and note for the chosen common walk: the pair of CSFs it makes of the loop."""
        cw, k = self.common, self.cw_spin.value()
        if cw is None:
            self.cw_label.setText('')
            self.loop_note = getattr(self, 'loop_value_note', '')
            return
        n = cw.count
        if not cw.segments:
            text = 'none (bra and ket share no walk outside the loop)'
        elif n == 1:
            text = '1 of 1'
        elif k == 0:
            text = f'all {n} drawn'
        else:
            text = f'{k} of {n}' + (' (first 2^31-1 reachable)' if n > 2**31 - 1 else '')
        self.cw_label.setText(text)
        note = self.loop_value_note
        if n == 1 or k > 0:
            g = self.g
            bra, ket = cw.step_vectors(max(k, 1) - 1)
            if (bra >= 0).all():
                def csf(steps):
                    i = g.csf_index(steps)
                    return (f'CSF {i + 1}' if i is not None
                            else f'irrep {g.irrep_name(g.walk_irrep(steps))}, not the target irrep')
                walks = (f'L = {g.csf_string(bra)} ({csf(bra)}), R = {g.csf_string(ket)} ({csf(ket)})')
            elif any(lp.kind.endswith('partial') for lp in cw.loops) and len(cw.loops) == 1:
                walks = 'select a partner partial loop too for the complete walks'
            else:
                walks = 'switch on upper and lower for the complete walks'
            note = (note + '\n' if note else '') + (f'common walk {text}: ' if cw.segments else '') + walks
        self.loop_note = note
        if note:
            self.statusBar().showMessage(note.replace('\n', ';  '))

    def overlay(self):
        if self.segments_mode():
            return None
        if not self.loop_items:
            if self.loop_mvs:
                return LoopOverlay(midvertices=self.loop_mvs)
            return None
        sel = self._selected_loops()
        if not sel:
            return LoopOverlay(loops=self.loop_items, midvertices=self.loop_mvs, detail=False)
        k = self.cw_spin.value()
        if self.common is not None and (self.common.count == 1 or k > 0):
            ctx = {(u, c) for u, c, _ in self.common.arcs(max(k, 1) - 1)}
        else:
            ctx = set()
            for lp in sel:
                ctx.update(lps.common_parts(self.g, lp, self.loop_mv, upper=self.upper_check.isChecked(),
                                            lower=self.lower_check.isChecked()))
        return LoopOverlay(loops=sel, context=sorted(ctx), midvertices=self.loop_mvs, detail=True, note=self.loop_note)

    # ------------------------------------------------------------------ drawing
    def set_option(self, name, on):
        setattr(self.opts, name, on)
        self.refresh()

    def choose_logo_color(self):
        c = QtWidgets.QColorDialog.getColor(QtGui.QColor(self.opts.logo_color), self, 'Logo color')
        if c.isValid():
            self.set_logo_color(c.name())

    def set_logo_color(self, color):
        self.opts.logo_color = color
        self.logo_button.setIcon(_swatch(color))
        self.refresh()

    def reset_view(self):
        self._reset_view = True
        self.refresh()

    def _title(self):
        g = self.g
        return (f'{g.ci_type} {g.nactel}e/{g.nlev}o, 2S+1={g.ispin}, irrep {g.irrep_name(g.lsym)}'
                + (f', root {self.root_spin.value()}' if g.ci_vectors is not None else ''))

    LOOPS_TAB = 1
    SEGMENTS_TAB = 2

    def loops_mode(self):
        """True while the Loops tab is shown: the CSF walks are then hidden (kept for the other tabs)."""
        return getattr(self, 'tabs', None) is not None and self.tabs.currentIndex() == self.LOOPS_TAB

    def segments_mode(self):
        """True while the Segments tab is shown: only the segments are drawn, a click selects their vertex."""
        return getattr(self, 'tabs', None) is not None and self.tabs.currentIndex() == self.SEGMENTS_TAB

    # ------------------------------------------------------------------ segments
    def set_segment_vertex(self, v):
        self.seg_vertex = v
        g = self.g
        if v is None or g is None:
            self.seg_vinfo.setText('')
        else:
            L, N, a, b, c = (int(t) for t in g.drt[v])
            irr = g.irrep_name(g.level_sym[L - 1]) if L > 0 else '-'
            self.seg_vinfo.setText(f'vertex {v + 1}: level {L} (irrep {irr}), N = {N}, (a,b,c) = ({a},{b},{c}); '
                                   f'the segments span level {L} &rarr; {L - 1}' if L > 0 else
                                   f'vertex {v + 1}: the tail of the graph (level 0), no segments below it')
        self.recompute_segments()

    def _segment_text(self, s):
        g = self.g
        bl, kt = s.steps
        return (f'#{s.number:<2d} {s.category:12s} {s.role}  {bl}/{kt}  {s.value:+.5f}  '
                f'{s.vlt + 1},{s.vrt + 1}->{s.vlb + 1},{s.vrb + 1}')

    def _segment_tip(self, s):
        g = self.g
        bl, kt = s.steps
        L = int(g.drt[s.vlt, 0])
        t, bt = sgm.ITVPT[s.k], sgm.IBVPT[s.k]
        role = ('the vertex is the top-left (bra) vertex, as SGUGA tabulates it (ISGM/VSGM)' if s.role == 'L'
                else 'the vertex is the top-right (ket) vertex; the bra vertex is its partner')
        abc = {u: tuple(int(q) for q in g.drt[u, 2:5]) for u in (s.vlt, s.vrt, s.vlb, s.vrb)}
        return (f'<qt><b>Segment {s.number}: {s.category}</b> &mdash; {sgm.DESCRIPTION[s.category]}.<br>'
                f'Level {L} &rarr; {L - 1}; {role}.<br>'
                f'Top: L = vertex {s.vlt + 1} {abc[s.vlt]}, R = vertex {s.vrt + 1} {abc[s.vrt]} '
                f'(class {t}: {sgm.CLASS[t]}).<br>'
                f'Steps: bra d = {sgm.IC1[s.k]} ({bl}), ket d = {sgm.IC2[s.k]} ({kt}).<br>'
                f'Bottom: L = vertex {s.vlb + 1} {abc[s.vlb]}, R = vertex {s.vrb + 1} {abc[s.vrb]} '
                f'(class {bt}: {sgm.CLASS[bt]}).<br>'
                f'Value: {s.formula} with b = {s.b} (b of the top-left vertex) = <b>{s.value:+.6f}</b> '
                f'(SGUGA formula ISVC = {sgm.ISVC[s.k]}).</qt>')

    def recompute_segments(self):
        self.seg_list.blockSignals(True)
        self.seg_list.clear()
        self.seg_items = []
        self.seg_note = ''
        v = self.seg_vertex
        if self.g is None or v is None:
            self.seg_info.setText('no vertex selected')
        else:
            kinds = {k for k, c in self.seg_kind.items() if c.isChecked()}
            items = sgm.segments_at(self.g, v, self.seg_left.isChecked(), self.seg_right.isChecked(), kinds)
            self.seg_items = items
            nl = sum(s.role == 'L' for s in items)
            self.seg_info.setText(f'{len(items)} segments at vertex {v + 1}: {nl} as bra (L), '
                                  f'{len(items) - nl} as ket (R) vertex')
            for k, s in enumerate(items):
                it = QtWidgets.QListWidgetItem(self._segment_text(s))
                it.setToolTip(self._segment_tip(s))
                it.setData(Qt.UserRole, k)
                self.seg_list.addItem(it)
        self.seg_list.blockSignals(False)
        self.refresh()

    def on_segment_selection(self):
        sel = [self.seg_items[it.data(Qt.UserRole)] for it in self.seg_list.selectedItems()]
        self.seg_note = ''
        if len(sel) == 1:
            s = sel[0]
            self.seg_note = (f'segment {s.number} ({s.category}): bra {s.steps[0]}, ket {s.steps[1]}, '
                             f'value {s.formula} = {s.value:+.6f} (b = {s.b})')
            self.statusBar().showMessage(self.seg_note)
        self.refresh()

    def seg_overlay(self):
        if not self.segments_mode() or self.seg_vertex is None:
            return None
        sel = [self.seg_items[it.data(Qt.UserRole)] for it in self.seg_list.selectedItems()]
        if sel:
            return SegmentOverlay(segments=sel, vertex=self.seg_vertex, detail=True, note=self.seg_note)
        return SegmentOverlay(segments=self.seg_items, vertex=self.seg_vertex, detail=False)

    def _walk_list(self):
        if self.loops_mode() or self.segments_mode():
            return []
        return sorted(self.walks.values(), key=lambda w: w.slot)

    def refresh(self):
        g = self.g
        if g is None:
            return
        self.walk_list.clear()
        for key, w in sorted(self.walks.items(), key=lambda kw: kw[1].slot):
            it = QtWidgets.QListWidgetItem(_swatch(SERIES[w.slot]), w.label.strip())
            it.setData(Qt.UserRole, key)
            self.walk_list.addItem(it)
        keep = None if self._reset_view else (self.ax.get_xlim(), self.ax.get_ylim())
        self.opts.title = self._title()
        self._x, self._y = draw(g, self.ax, self._walk_list(), self.opts, self.overlay(), self.seg_overlay())
        self._new_hover()
        if keep:
            self.ax.set_xlim(keep[0])
            self.ax.set_ylim(keep[1])
        else:
            self.toolbar.update()  # new home view
        self._reset_view = False
        self.fig.tight_layout()
        self.canvas.draw_idle()

    def save_dialog(self):
        if self.g is None:
            return
        base = os.path.splitext(os.path.basename(self.filename))[0] + '.shavitt.pdf'
        fn, filt = QtWidgets.QFileDialog.getSaveFileName(self, 'Save the graph', base,
                                                         'PDF (*.pdf);;PNG (*.png);;SVG (*.svg)')
        if not fn:
            return
        if not os.path.splitext(fn)[1]:
            fn += '.' + filt.split('*.')[1].rstrip(')')
        self.save(fn)

    def save(self, fn):
        """Write the graph as shown (zoom, options, walks) to fn."""
        from .draw import save as save_figure
        g = self.g
        self.opts.title = self._title()
        lim = (self.ax.get_xlim(), self.ax.get_ylim())
        save_figure(g, fn, self._walk_list(), self.opts, dpi=self.dpi_spin.value(), limits=lim, loops=self.overlay(),
                    segments=self.seg_overlay())
        self.statusBar().showMessage(f'Saved {fn}', 6000)


def run(filename=None, root=1, threshold=0.05, opts=None, graph=None):
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    win = ShavittWindow(filename, root, threshold, opts, graph)
    win.show()
    return app.exec()
