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
        subroutine GetRest_t3 (t1,t1_tmp,E2old)
c
c        this file read 1) T1o
c                      2) E1old,E2old,niter
c        from RstFil file
c
        implicit none
#include "cht3_ccsd1.fh"
#include "ccsd_t3compat.fh"
        real*8 E1old,E2old
        real*8 t1(*),t1_tmp(*)
c
c        help variables
        integer length,i
c
*       open (unit=LunAux,File='RstFil',form='unformatted')
        Call MOLCAS_BinaryOpen_Vanilla(LunAux,'RstFil')
        length=nv*no
cmp        write (*,*) 'no, nv, length = ',no,nv,length
        call cht3_rea (LunAux,length,t1)
c
        call transp (t1,t1_tmp,nv,no)
c
        do i=1,length
        t1(i+length)=t1_tmp(i)
        t1(i)=t1_tmp(i)
        end do
c
c
        read (LunAux) E1old,E2old,i

        if (printkey.gt.1) then
        write (6,'(A,2(f15.12,1x))') 'Results from CCSD : E1, E2 ',
     & E1old,E2old
        end if

        close (LunAux)
c
c
        return
        end
