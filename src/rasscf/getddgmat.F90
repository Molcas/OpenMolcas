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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************
      Subroutine GetDDgMat(DDg,GDMat,Gtuvx)
      use rasscf_global, only: lRoots, NAC
      Implicit None


#include "warnings.h"
      Real*8,DIMENSION(lRoots,lRoots,lRoots,lRoots)::DDG
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::Gtuvx
      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GDMat

      INTEGER iI,iJ,iK,iL,it,iu,iv,ix,iII,iJJ,iKK,iLL
      DO iI=1,lRoots
       DO iJ=1,lRoots
        IF(iJ.gt.iI) THEN
         iJJ=iI
         iII=iJ
        ELSE
         iII=iI
         iJJ=iJ
        END IF
        DO iK=1,lRoots
         DO iL=1,lRoots
          IF(iL.gt.iK) THEN
           iLL=iK
           iKK=iL
          ELSE
           iLL=iL
           iKK=iK
          END IF
          DDG(iI,iJ,iK,iL)=0.0d0
          do it=1,NAC
           do iu=1,NAC
            do iv=1,NAC
             do ix=1,NAC
              DDG(iI,iJ,iK,iL)=DDG(iI,iJ,iK,iL)                         &
     & +GDMat(iII*(iII-1)/2+iJJ,it,iu)*GDMat(iKK*(iKK-1)/2+iLL,iv,ix)   &
     & *Gtuvx(it,iu,iv,ix)
             end do
            end do
           end do
          end do
         END DO
        END DO
       END DO
      END DO
      RETURN
      End Subroutine GetDDgMat
