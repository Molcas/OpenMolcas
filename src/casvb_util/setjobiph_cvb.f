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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      subroutine setjobiph_cvb(iorcore_c,iorclos_c,iorocc_c,mxirrep,
     >  nstsym_c,weight_c,istnel_c,istsy_c,istms2_c,nstats_c,
     >  mxstt_ci,mxstsy_ci,nel_c,norb_c,i2s_c,isym_c,mcore_c,neltot_c)
      implicit real*8 (a-h,o-z)
#include "rasdim.fh"
#include "jobiph_j.fh"
      dimension iorcore_c(mxirrep),iorclos_c(mxirrep),iorocc_c(mxirrep)
      dimension weight_c(mxstt_ci,mxstsy_ci)
      dimension istnel_c(mxstsy_ci),istsy_c(mxstsy_ci)
      dimension istms2_c(mxstsy_ci),nstats_c(mxstsy_ci)

c  Orbitals  --  OCC, CLOSED and CORE cards
      call imove_cvb(nfro_j,iorcore_c,mxirrep)
      call imove_cvb(nish_j,iorclos_c,mxirrep)
      call imove_cvb(nrs2_j,iorocc_c,mxirrep)
c  States  --  WF, STATE and WEIGHT cards
      nstsym_c=1
      call fzero(weight_c,mxstt_ci*mxstsy_ci)
      do i=1,lroots_j
      wgt=0.0d0
      do j=1,nroots_j
      if (iroot_j(j).eq.i) wgt=weight_j(j)
      end do
      if(wgt.ne.0d0)then
        if(i.gt.mxstt_ci)then
        write(6,*)' Root number too large in casrecov_cvb :',i,mxstt_ci
          call abend_cvb()
        endif
      endif
      weight_c(i,1)=wgt
      enddo

      istnel_c(1)=nactel_j
      istsy_c(1)=lsym_j
      istms2_c(1)=ispin_j-1
      nstats_c(1)=lroots_j
c  Set derived info
      nel_c=nactel_j
      i2s_c=ispin_j-1
      isym_c=lsym_j
      norb_c=0
      mcore_c=0
      do i=1,mxstsy_ci
      norb_c=norb_c+nrs2_j(i)
      mcore_c=mcore_c+nfro_j(i)+nish_j(i)
      enddo
      neltot_c=nel_c+2*mcore_c
c  MO common block :
      call setmocom_cvb()
      return
      end
