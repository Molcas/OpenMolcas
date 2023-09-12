!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************
      subroutine ddinit_cvb(method,nparm1,nfrdim1,maxd1,mxit1,          &
     &  ifollow1,isaddle1,ip1,corenrg1,n_div1)
      implicit real*8 (a-h,o-z)
      character*(*) method
#include "direct_cvb.fh"

!  Input parameters:
      nparm=nparm1
      nfrdim=nfrdim1
      maxd=maxd1
      mxit=mxit1
      ifollow=ifollow1
      isaddle=isaddle1
      nroot=max(1,isaddle1+1)
      ip=ip1
      corenrg=corenrg1
      n_div=n_div1
!  Defaults:
      resthr=1d-5
!  Local DIRDIAG parameters:
      orththr=1d-10
      nortiter=50

      if(method.eq.'AxEx')then
        mxrhs=0
        ivrhs=0
!  Arrays:  C       AxC     RES
!           AP      SOLP    SOLP_RES
        idd(1) = mstackr_cvb(nparm*maxd)
        idd(2) = mstackr_cvb(nparm*maxd)
        idd(3) = mstackr_cvb(nparm)
        idd(4) = mstackr_cvb(maxd*maxd)
        idd(5) = mstackr_cvb(maxd)
        idd(6) = mstackr_cvb(maxd)
      elseif(method.eq.'AxESx')then
        mxrhs=0
        ivrhs=0
!  Arrays:  C       AxC     SxC     RES
!           AP      SOLP    SOLP_RES
        idd(1) = mstackr_cvb(nparm*maxd)
        idd(2) = mstackr_cvb(nparm*maxd)
        idd(3) = mstackr_cvb(nparm*maxd)
        idd(4) = mstackr_cvb(nparm)
        idd(5) = mstackr_cvb(maxd*maxd)
        idd(6) = mstackr_cvb(maxd)
        idd(7) = mstackr_cvb(maxd)
      elseif(method.eq.'Axb')then
        mxrhs=1
        ivrhs=4
!  Arrays:  C       SxC     RES     RHS
!           RHSP    SOLP    SOLP_RES
        idd(1) = mstackr_cvb(nparm*maxd)
        idd(2) = mstackr_cvb(nparm*maxd)
        idd(3) = mstackr_cvb(nparm)
        idd(4) = mstackr_cvb(nparm)
        idd(5) = mstackr_cvb(maxd)
        idd(6) = mstackr_cvb(maxd)
        idd(7) = mstackr_cvb(maxd)
      elseif(method.eq.'AxExb')then
        mxrhs=1
        ivrhs=4
!  Arrays:  C       AxC     RES     RHS
!           AP      RHSP    SOLP    SOLP_RES
        idd(1) = mstackr_cvb(nparm*maxd)
        idd(2) = mstackr_cvb(nparm*maxd)
        idd(3) = mstackr_cvb(nparm)
        idd(4) = mstackr_cvb(nparm)
        idd(5) = mstackr_cvb(maxd*maxd)
        idd(6) = mstackr_cvb(maxd)
        idd(7) = mstackr_cvb(maxd)
        idd(8) = mstackr_cvb(maxd)
      endif
      return
      end
