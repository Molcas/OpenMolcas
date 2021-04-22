************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
        subroutine SaveChckData (LunAux)
c
c        this routine do:
c        Save Chck Data to ChkDat file
c
        implicit none
#include "chcc1.fh"
        integer LunAux
c
C       open (unit=LunAux,file='ChKDat',form='unformatted')
        Call Molcas_BinaryOpen_Vanilla(LunAux,'ChKDat')
          write(LunAux) T1c,T2c,OEo,OEv,Q0,Q1,Q21,Q22,Q3,Q4
     c                 ,L0k,L1k,L2k
        close (LunAux)
c
        return
        end
