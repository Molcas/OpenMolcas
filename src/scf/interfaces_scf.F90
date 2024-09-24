!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

module Interfaces_SCF

implicit none
private

public :: dOne_SCF, MinDns, OccDef, PMat_SCF, TraClc_i, vOO2OV

interface

  subroutine dOne_SCF(nSym,nBas,nOrb,nFro,CMO,nCMO,Occ,Dlt,alpha_density)
    integer :: nSym, nCMO, nBas(nSym), nOrb(nSym), nFro(nSym)
    real*8, target :: CMO(nCMO), Occ(*), Dlt(*)
    logical :: alpha_density
  end subroutine dOne_SCF

  subroutine MinDns(Dens,mBT,NumD,XCff,ltXCff,nD)
    integer :: mBT, NumD, ltXCff, nD
    real*8 :: XCff(ltXCff,nD)
    real*8, target :: Dens(mBT,nD,NumD)
  end subroutine MinDns

  subroutine OccDef(Occ,mmB,nD,CMO,mBB)
    integer :: mmB, nD, mBB
    real*8 :: Occ(mmB,nD)
    real*8, target :: CMO(mBB,nD)
  end subroutine OccDef

  subroutine PMat_SCF(FstItr,XCf,nXCf,nD)
    logical :: FstItr
    integer :: nXCf, nD
    real*8 :: XCf(nXCf,nD)
  end subroutine PMat_SCF

  subroutine TraClc_i(iterLw,nD)
    integer :: iterLw, nD
  end subroutine TraClc_i

  subroutine vOO2OV(v1,n1,v2,n2,nD,n3)
    integer :: n1, n2, nD
    integer :: n3(nD)
    real*8, target :: v1(n1,nD), v2(n2)
  end subroutine vOO2OV

end interface

end module Interfaces_SCF
