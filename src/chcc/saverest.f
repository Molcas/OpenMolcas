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
        subroutine SaveRest (wrk,wrksize,LunAux,niter,E1old,E2old)
c
c        this file save 1) T1o,OE
c                      2) E1old,E2old,niter
c        into RstFil file
c
        implicit none
#include "chcc1.fh"
#include "wrk.fh"
        integer LunAux,niter
        real*8 E1old,E2old
c
c        help variables
        integer len
c
*       open (unit=LunAux,File='RstFil',form='unformatted')
        Call MOLCAS_BinaryOpen_Vanilla(LunAux,'RstFil')
        len=no*nv
        call wri_chcc (LunAux,len,wrk(PossT1o))
        write (LunAux) E1old,E2old,niter
        close (LunAux)
c
c
        return
        end
