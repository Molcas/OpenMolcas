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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 07, 2022, created this file.               *
!*****************************************************************
      Subroutine UnzipTUVX(TUVX,gtuvx,nTUVX)
      use rasscf_global, only: NACPR2, NAC
      Implicit None


#include "warnings.h"
      INTEGER nTUVX
      Real*8 gtuvx(nTUVX),TUVX(NACPR2)
      INTEGER it,iu,iv,ix,ituvx,ixmax,                                  &
     &        jtuvx,jtuxv,jutvx,jutxv,                                  &
     &        jvxtu,jvxut,jxvtu,jxvut,                                  &
     &        NAC3,NAC2

!      CALL FZero(gtuvx,nTUVX)

      NAC2=NAC**2
      NAC3=NAC2*NAC
      ituvx=0
      DO it=1,NAC
       Do iu=1,it
        dO iv=1,it
         ixmax=iv
         if (it==iv) ixmax=iu
         do ix=1,ixmax
          ituvx=ituvx+1
          jtuvx=(it-1)*NAC3+(iu-1)*NAC2+(iv-1)*NAC+ix
          jtuxv=(it-1)*NAC3+(iu-1)*NAC2+(ix-1)*NAC+iv
          jutvx=(iu-1)*NAC3+(it-1)*NAC2+(iv-1)*NAC+ix
          jutxv=(iu-1)*NAC3+(it-1)*NAC2+(ix-1)*NAC+iv
          jvxtu=(iv-1)*NAC3+(ix-1)*NAC2+(it-1)*NAC+iu
          jvxut=(iv-1)*NAC3+(ix-1)*NAC2+(iu-1)*NAC+it
          jxvtu=(ix-1)*NAC3+(iv-1)*NAC2+(it-1)*NAC+iu
          jxvut=(ix-1)*NAC3+(iv-1)*NAC2+(iu-1)*NAC+it
          Gtuvx(jtuvx)=TUVX(ituvx)
          Gtuvx(jtuxv)=TUVX(ituvx)
          Gtuvx(jutvx)=TUVX(ituvx)
          Gtuvx(jutxv)=TUVX(ituvx)
          Gtuvx(jvxtu)=TUVX(ituvx)
          Gtuvx(jvxut)=TUVX(ituvx)
          Gtuvx(jxvtu)=TUVX(ituvx)
          Gtuvx(jxvut)=TUVX(ituvx)
         end do
        eND dO
       End Do
      END DO
      RETURN
      End Subroutine
