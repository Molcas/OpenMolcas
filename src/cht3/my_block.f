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
        subroutine my_block (vblock,vblock_my)
c
c this subroutine calculates maximum overlap between juzek's
c vblock segmentation and palo's dimgrp
c
        implicit none
c
#include "cht3_ccsd1.fh"
#include "ccsd_t3compat.fh"
#include "cht3_reord.fh"
c
        integer i,j,i_tmp,i_f,i_l,poss
        integer vblock,vblock_my,vblock_my_tmp
        logical found
c
        vblock_my=0
        i_l=0
        i_f=0
c
        do i=1,nv,vblock
c
c - find initial possition of the i-th juzek's block
c
          poss=0
          found=.false.
        do j=1,NvGrp
          poss=poss+DimGrpaR(j)
          if ((i.le.poss).and.(.not.found)) then
            i_f=j
            found=.true.
cmp        write (6,'(A,3(i5,x))') 'i,i_f,poss     = ',
cmp     & i,i_f,poss
          end if
        end do
c
        if ((i+vblock-1).le.nv) then
          i_tmp=i+vblock-1
        else
          i_tmp=nv
        end if
cmp        write (6,'(A,2(i5,x))') 'i,i_tmp        = ',
cmp     & i,i_tmp
c
c - find terminal possition of the i-th juzek's block
c
          poss=0
          found=.false.
        do j=1,NvGrp
          poss=poss+DimGrpaR(j)
          if ((i_tmp.le.poss).and.(.not.found)) then
            i_l=j
            found=.true.
cmp        write (6,'(A,3(i5,x))') 'i_tmp,i_l,poss = ',
cmp     & i_tmp,i_l,poss
          end if
        end do
c
        vblock_my_tmp=0
        do j=i_f,i_l
        vblock_my_tmp=vblock_my_tmp+DimGrpaR(j)
        end do
c
        if (vblock_my_tmp.gt.vblock_my)
     & vblock_my=vblock_my_tmp
c
cmp        write (6,'(A,2(i5,x))') 'vblock_my_tmp, vblock_my',
cmp     & vblock_my_tmp, vblock_my
cmp        write (6,*)
        end do
c
        return
        end
