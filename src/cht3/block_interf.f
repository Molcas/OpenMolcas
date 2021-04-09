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
        subroutine block_interf(
     &   ind1f,ind1l,
     &   ind2f,ind2l,
     &   b1f,b1l,nind_b1f,nind_b1l,
     &   b2f,b2l,nind_b2f,nind_b2l)
c
c this routine do :
c
c Interface between Palo and Juzek virtual orbitals blocking structure
c
c ind1f,ind1l,ind2f,ind2l - first and last absolute of the VO indexes of interess
c
c b1f,b1l   - first and last Palo's block which contain ind1
c b2f,b2l   - first and last Palo's block which contain ind2
c
c nind_b1f  - sum of VOs in blocks 1,2, ..., b1f-1
c nind_b2f  - sum of VOs in blocks 1,2, ..., b2f-1
c
c nind_b1l  - # of VOs in b1f before ind1
c nind_b2l  - # of VOs in b2f before ind2
c
        implicit none
#include "cht3_ccsd1.fh"
#include "cht3_reord.fh"
#include "ccsd_t3compat.fh"
c
        integer i,sum
        integer ind1f,ind1l,ind2f,ind2l
        integer b1f,b2f,b1l,b2l
        integer nind_b1f,nind_b2f,nind_b1l,nind_b2l
        logical found1,found2,found3,found4
c
c set b1f, b2f, b1l, b2l
c
          sum=0
          found1=.false.
          found2=.false.
          found3=.false.
          found4=.false.
c
        do i=1,NvGrp
          sum=sum+DimGrpaR(i)
c
          if ((ind1f.le.sum).and.(.not.found1)) then
            b1f=i
            found1=.true.
          end if
c
          if ((ind1l.le.sum).and.(.not.found2)) then
            b1l=i
            found2=.true.
          end if
c
          if ((ind2f.le.sum).and.(.not.found3)) then
            b2f=i
            found3=.true.
          end if
c
          if ((ind2l.le.sum).and.(.not.found4)) then
            b2l=i
            found4=.true.
          end if
c
        end do
c
cmp        write (*,*) 'b1f, b1l, b2f, b2l ',b1f,b1l,b2f,b2l
c
c set nind_b1f, nind_b1l
c
        if (b1f.gt.1) then
            sum=0
          do i=1,b1f-1
            sum=sum+DimGrpaR(i)
          end do
            nind_b1f=sum
        else
            nind_b1f=0
        end if
            nind_b1l=ind1f-nind_b1f-1
c
c set nind_b1f, nind_b1l
c
        if (b2f.gt.1) then
            sum=0
          do i=1,b2f-1
            sum=sum+DimGrpaR(i)
          end do
            nind_b2f=sum
        else
            nind_b2f=0
        end if
            nind_b2l=ind2f-nind_b2f-1
c
        return
        end
