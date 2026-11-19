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
# Copyright (C) 2026, Benjamin W. Kaufold                              *
#               2026, Meredith E. Wanta                                *
#               2026, Sijia S. Dong                                    *
#***********************************************************************

class FakeElectronicState:
    def __init__(self, state, dipole, vector=None):
        self.state = state
        self.dipole = dipole
        self.vector = vector


class FakeMultiRefCalc:
    def __init__(self, states, path='mr.log'):
        self.states = states
        self.path = path


class FakeSingleRefCalc:
    def __init__(self, state, path='ref.orca'):
        self.state = state
        self.path = path

class FakeVector:
    def __init__(self, x, y, z):
        self.x, self.y, self.z = x, y, z