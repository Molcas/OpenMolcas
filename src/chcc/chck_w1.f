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
        subroutine Chck_W1 (W1,aSGrp,beSGrp,bSGrp,gaSGrp)
c
c        cek W1
        implicit none
#include "chcc1.fh"
        integer aSGrp,beSGrp,bSGrp,gaSGrp
        real*8 W1(1:16*31,1:16*33)
c
        integer a,b,be,ga,bad,ap,bp,bep,gap,ab,bega
        real*8 s
c
        if (aSGrp.eq.2) then
          ap=nv/2
        else
          ap=0
        end if
c
        if (bSGrp.eq.2) then
          bp=nv/2
        else
          bp=0
        end if
c
        if (gaSGrp.eq.2) then
          gap=nv/2
        else
          gap=0
        end if
c
        if (beSGrp.eq.2) then
          bep=nv/2
        else
          bep=0
        end if
c
        bad=0
        bega=0
        do be=1,nv/2
        do ga=1,be
        bega=bega+1
        ab=0
        do a=2,nv/2
        do b=1,a-1
        ab=ab+1
          s=(Q4(ap+a,bep+be,bp+b,gap+ga)
     c      +Q4(ap+a,gap+ga,bp+b,bep+be))/1
          if (abs(W1(ab,bega)-s).gt.1.0d-10) then
          bad=bad+1
c          write (6,99) a,b,be,ga,ab,bega,s,W1(a,be,ga)
c99        format(4(i2,1x),2(i6,1x),2(f15.10))
          end if
        W1(ab,bega)=s
        end do
        end do
        end do
        end do
c
        if (bad.eq.0) then
        write (6,*) ' Chck W OK ', bad
        else
        write (6,*) ' Chck W Bug !!!!!!! ', bad
        end if
c
        return
        end
