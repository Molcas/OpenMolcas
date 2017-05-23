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
* Copyright (C) 2004, Alexander Wolf                                   *
*               2004,2006, Markus Reiher                               *
************************************************************************
      subroutine output1 (dkhorder,xorder,paramtype,dkhscfflg)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.2
c
c   last modified: 11.10.2006 (MR, ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer dkhorder,xorder,ramused
      character*(3) paramtype
      logical dkhscfflg
c
      If (DKH_Verbose) Then
      write (stdout,1001)
1001  format (//2X,'DKH PARSER : Symbolic construction and evaluation',
     *        ' of the Douglas--Kroll--Hess Hamiltonian',/2X,
     *        12('-'),/15X,'written by:  M. Reiher and A. Wolf  ',
     *        '(Univ. Bonn, 2004; Univ. Jena 2005, ETH Zurich 2006)',
     *        //15X,'References:  M.Reiher, A.Wolf, J. Chem. Phys. 121',
     *        ' (2004) 10945  (Hamiltonian),',
     *        /28X,'A.Wolf, M.Reiher, J. Chem. Phys. 124 (2006) ',
     *        '064103 (properties)',//2X)
c
      write (stdout,1011) dkhorder,xorder,paramtype
1011  format (5X,'DKH Hamiltonian order:',5X,I2,/5X,
     *        'DKH property order:',8X,I2,/5X,'Chosen parametrization',
     *        ' of unitary transformations U_i:',5X,A3,/)
c
      if (dkhscfflg) then
        write (stdout,1012)
1012    format (5X,'Treat property X variationally, i.e., U_i = ',
     *          'U_i(X).')
      else
        write (stdout,1013)
1013    format (5X,'Treat property X perturbatively, i.e., U_i is ',
     *          'independent of X.')
      endif
c
CMR      write (stdout,1015) maxoperators,maxuops,maxlength
CMR1015  format (/5X,'Max. number of operators:',6X,
CMR     *        'maxoperators =',I9,/5X,'Max. number of terms for U:',4X,
CMR     *        'maxuops',6X,'=',I9,/5X,'Max. length of each operator:',
CMR     *        2X,'maxlength',4X,'=',5X,I4)

      write (stdout,1023) dkhzero
1023  format (5X,'Zero threshold:',16X,'dkhzero',6X,'=',D13.5)
      if (dynthrsh.eq.0) write (stdout,1025) dynthrsh
      if (dynthrsh.eq.1) write (stdout,1026) dynthrsh
1025  format (5X,'No dynamical thresholding:     dynthrsh = ',I1,/)
1026  format (5X,'Use dynamical thresholding:    dynthrsh = ',I1,/)
      End If
c
c   Calculate approximate RAM requirements in MByte
c
      ramused = 4*maxoperators*(maxlength+20)
      ramused = ramused + (maxuops*(maxlength+24))
      ramused = ramused/8388608
c
CMR      write (stdout,1031) ramused
CMR1031  format (5X,'-->  Approximately ',I4,' MB RAM needed for symbolic',
CMR     *        ' DKH PARSER subroutines.'//)
c
      return
      end
