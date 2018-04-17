************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2010, Grigory A. Shamov                                *
************************************************************************

                        Additions to dft_util.

Working on a hydrocarbon's DFT benchmarking project, I felt that
I miss some of the density functionals; so I added them in MOLCAS
and would like to share it.

The following routines for functional cores were added:

  KealTozer, CTCA, CWIG, XG96, XB86, XOPT, XRGE2, XSSBSW

Derivation was done using Maxima, following spirit of Pawel Salek's Funclib.
I did not use his scripts nor Funclib Maxima code, for I wanted Fortran77 conversion
and different calling conwentions from his.


Based on these added functionals and the ones that existed in MOLCAS already,
the following combinations were created.

  B2PLYP B86LYP B86PBE BPBE PBESOL GLYP GPBE KT2 KT3 OLYP
  OPBE O2PLYP O3LYP PTCA RGE2 SSBSW

I've added references and minimal info/notes in the source files.
For double-hybrid functionals, MBPT2 code can be modified to check the RunFile
if KS-DFT is in use and whether functionals are double-hybrids (B2PLYP or O2PLYP).
I did a hack that works for me, but probably it is not general enough.


Grigory A Shamov,
University of Manitoba
April 2010

