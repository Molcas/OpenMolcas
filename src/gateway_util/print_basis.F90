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
! Copyright (C) 2006, Roland Lindh                                     *
!***********************************************************************

subroutine Print_Basis(lOPTO)
!***********************************************************************
!                                                                      *
!     Object: to print the basis set                                   *
!                                                                      *
!     Author: Roland Lindh, Dept Chem. Phys., Lund University, Sweden  *
!             September 2006                                           *
!***********************************************************************

use Basis_Info
use Center_Info
use RICD_Info, only: Thrshld_CD

implicit real*8(A-H,O-Z)
#include "angtp.fh"
#include "relmp.fh"
#include "real.fh"
#include "print.fh"
character DBas*4
character ChCo*1, ChCa*1, ChSph*1
logical Output, type(0:7), lOPTO

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)
if (iPrint == 0) return
LuWr = 6
!                                                                      *
!***********************************************************************
!                                                                      *
if (Show) then
  write(LuWr,*)
  call CollapseOutput(1,'   Basis set information:')
  write(LuWr,'(3X,A)') '   ----------------------'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Write out basis set information
!                                                                      *
!***********************************************************************
!                                                                      *
! Valence basis set

do iCnttp=1,nCnttp
  mdc = dbsc(iCnttp)%mdci
  lSh = 0
  output = Show
  if (dbsc(iCnttp)%Aux .or. dbsc(iCnttp)%Frag) output = output .and. (iPrint >= 10) .and. (iCnttp /= iCnttp_Dummy)
  if (output) then
    write(LuWr,*)
    write(LuWr,*)
    write(LuWr,'(6X,A,1X,A)') 'Basis set label:',trim(dbsc(iCnttp)%Bsl)
    if (lOPTO) goto 100
    write(LuWr,*)
    dbas = dc(mdc+1)%LblCnt(1:4)
    call Upcase(dbas)
    if (dbas == 'DBAS') then
      write(LuWr,'(6X,A)') 'Diffuse basis set for R-matrix:'
      write(LuWr,'(6X,A)') '==============================='
      if (dbsc(iCnttp)%nCntr /= 1) then
        call WarningMessage(2,'Too many centers, should only be one!')
        call Quit_OnUserError()
      end if
    else
      if (dbsc(iCnttp)%Aux) then
        write(LuWr,'(6X,A)') 'Auxiliary basis set:'
        write(LuWr,'(6X,A)') '=================='
        if (dbsc(iCnttp)%aCD_Thr /= One) then
          write(LuWr,'(6X,A,G9.2)') 'Threshold in the auxiliary basis set generation is modified to ', &
                                    dbsc(iCnttp)%aCD_Thr*Thrshld_CD
        end if
      else if (dbsc(iCnttp)%Frag) then
        write(LuWr,'(6X,A)') 'Fragment basis set:'
        write(LuWr,'(6X,A)') '=================='
      else
        if (dbsc(iCnttp)%fMass == One) then
          write(LuWr,'(6X,A)') 'Electronic valence basis set:'
          write(LuWr,'(6X,A)') '------------------'
        else
          write(LuWr,'(6X,A)') 'Muonic valence basis set:'
          write(LuWr,'(6X,A)') '------------------'
        end if
      end if
    end if
    if (dbsc(iCnttp)%Fixed) write(LuWr,'(6X,A)') 'Centers of this basis set are frozen!'
    if (dbsc(iCnttp)%IsMM == 1) then
      write(LuWr,'(6X,A)') 'This is a MM atom: no basis set'
    else
      if (dbsc(iCnttp)%pChrg) then
        write(LuWr,'(6X,A,F10.6,A)') 'Associated Effective Charge ',dbsc(iCnttp)%Charge,' au (this is a pseudo charge)'
      else
        write(LuWr,'(6X,A,F10.6,A)') 'Associated Effective Charge ',dbsc(iCnttp)%Charge,' au'
      end if
      write(LuWr,'(6X,A,F10.6,A)') 'Associated Actual Charge    ',max(Zero,dble(dbsc(iCnttp)%AtmNr)),' au'

      if (Nuclear_Model == Point_Charge) then
        write(LuWr,'(6X,A)') 'Nuclear Model: Point charge'
      else if (Nuclear_Model == Gaussian_type) then
        write(LuWr,'(6X,A)') 'Nuclear Model: Finite nucleus - Gaussian distribution'
        write(LuWr,'(6X,A,E12.5)') '  Gaussian exponent, Xi/bohr**(-2): ',dbsc(iCnttp)%ExpNuc
      else if (Nuclear_Model == mGaussian_type) then
        write(LuWr,'(6X,A)') 'Nuclear Model: Finite nucleus - Modified Gaussian distribution'
        write(LuWr,'(6X,A,E12.5,A,E12.5)') '  Parameters, Xi/bohr**(-2), w/bohr**(-2): ',dbsc(iCnttp)%ExpNuc,', ', &
                                           dbsc(iCnttp)%w_mGauss
      else
        call WarningMessage(2,'Illegal Nuclear Model!')
        call Abend()
      end if
    end if
    write(LuWr,*)
100 continue

  end if
  kShStr = dbsc(iCnttp)%iVal
  kShEnd = dbsc(iCnttp)%iVal+dbsc(iCnttp)%nVal-1
  type(0) = .false.
  do kSh=kShStr,kShEnd
    nExpk = Shells(kSh)%nExp
    type(0) = type(0) .or. (nExpk*Shells(kSh)%nBasis /= 0)
  end do
  if (output .and. type(0) .and. (.not. lOPTO)) then
    write(LuWr,'(6X,A)') 'Shell  nPrim  nBasis  Cartesian Spherical Contaminant'
  end if
  do kSh=kShStr,kShEnd
    nExpk = Shells(kSh)%nExp
    ChCa = ' '
    ChSph = 'X'
    ChCo = ' '
    if (.not. Shells(kSh)%Transf) then
      ChCa = 'X'
      ChSph = ' '
    end if
    if (Shells(kSh)%Transf .and. (.not. Shells(kSh)%Prjct)) ChCo = 'X'
    if (output .and. (nExpk*Shells(kSh)%nBasis /= 0) .and. (.not. lOPTO)) &
      write(LuWr,'(9X,A,5X,I3,5X,I3,8X,A,8X,A,8X,A)') AngTp(lSh),nExpk,Shells(kSh)%nBasis,ChCa,ChSph,ChCo
    lSh = lSh+1
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Process PP part, if any.

  kShStr = dbsc(iCnttp)%iPP
  kShEnd = kShStr+dbsc(iCnttp)%nPP-1
  if (output .and. (dbsc(iCnttp)%nPP /= 0) .and. (.not. lOPTO)) then
    write(LuWr,*)
    write(LuWr,'(6X,A)') 'Pseudo Potential specification:'
    write(LuWr,'(6X,A)') '======================================='
    type(0) = .false.
    do kSh=kShStr,kShEnd
      nExpk = Shells(kSh)%nExp/3
      type(0) = type(0) .or. (nExpk /= 0)
    end do
    if (type(0)) then
      write(LuWr,*)
      write(LuWr,'(6X,A)') 'Potential  nTerms    '
    end if
  end if
  lSh = 0
  do kSh=kShStr,kShEnd
    nExpk = Shells(kSh)%nExp/3
    !write(6,*) 'kSh,lSh=',kSh,lSh
    if (output .and. (nExpk /= 0) .and. (.not. lOPTO)) then
      if (lSh == 0) then
        write(LuWr,'(9X,A,6X,I2)') '  H',nExpk
      else
        write(LuWr,'(9X,A,6X,I2)') AngTp(lSh-1)//'-H',nExpk
      end if
    end if
    lSh = lSh+1
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Process ECP part, if any.

  kShStr = dbsc(iCnttp)%iPrj
  kShEnd = kShStr+dbsc(iCnttp)%nPrj-1
  if (output .and. dbsc(iCnttp)%ECP .and. (.not. lOPTO)) then
    write(LuWr,*)
    write(LuWr,'(6X,A)') 'Effective Core Potential specification:'
    write(LuWr,'(6X,A)') '======================================='
    if (dbsc(iCnttp)%nM1 > 0) then
      write(LuWr,*)
      write(LuWr,'(6X,A,I5)') ' Number of M1 terms:',dbsc(iCnttp)%nM1
    end if
    if (dbsc(iCnttp)%nM2 > 0) then
      write(LuWr,*)
      write(LuWr,'(6X,A,I5)') ' Number of M2 terms:',dbsc(iCnttp)%nM2
    end if
    type(0) = .false.
    do kSh=kShStr,kShEnd
      nExpk = Shells(kSh)%nExp
      nBasisk = Shells(kSh)%nBasis
      type(0) = type(0) .or. (nExpk*nBasisk /= 0)
    end do
    if (type(0)) then
      write(LuWr,*)
      write(LuWr,'(6X,A)') 'Projection basis set '
      write(LuWr,'(6X,A)') 'Shell  nPrim  nBasis '
    end if
  end if
  lSh = 0
  do kSh=kShStr,kShEnd
    nExpk = Shells(kSh)%nExp
    nBasisk = Shells(kSh)%nBasis
    if (output .and. (nExpk*nBasisk /= 0) .and. (.not. lOPTO)) write(LuWr,'(9X,A,6X,I2,6X,I2)') AngTp(lSh),nExpk,nBasisk
    lSh = lSh+1
  end do
  ! Spectral resolution basis set
  kShStr = dbsc(iCnttp)%iSRO
  kShEnd = kShStr+dbsc(iCnttp)%nSRO-1
  if (output .and. dbsc(iCnttp)%ECP .and. (.not. lOPTO)) then
    type(0) = .false.
    do kSh=kShStr,kShEnd
      nExpk = Shells(kSh)%nExp
      type(0) = type(0) .or. (nExpk /= 0)
    end do
    if (type(0)) then
      if (dbsc(iCnttp)%nOpt /= 0) then
        write(LuWr,*)
        write(LuWr,'(6X,A)') 'Spectral Resolvent Operators :'
        if (iand(2**0,dbsc(iCnttp)%nOpt) /= 0) write(LuWr,'(8X,A)') ' Exchange'
        if (iand(2**1,dbsc(iCnttp)%nOpt) /= 0) write(LuWr,'(8X,A)') ' Mass-Velocity'
        if (iand(2**2,dbsc(iCnttp)%nOpt) /= 0) write(LuWr,'(8X,A)') ' Darwin 1-electron contact term'
        if (iand(2**3,dbsc(iCnttp)%nOpt) /= 0) then
          if (IRELMP == 0) then
            write(LuWr,'(8X,A)') ' No-Pair approximation'
          else if (IRELMP == 1) then
            write(LuWr,'(8X,A)') ' No-Pair approximation (DK1)'
          else if (IRELMP == 2) then
            write(LuWr,'(8X,A)') ' No-Pair approximation (DK2)'
          else if (IRELMP == 3) then
            write(LuWr,'(8X,A)') ' No-Pair approximation (DK3)'
          else if (IRELMP == 4) then
            write(LuWr,'(8X,A)') ' No-Pair approximation (DK3)'
          else if (IRELMP == 11) then
            write(LuWr,'(8X,A)') ' RESC approximation'
          else if (IRELMP == 21) then
            write(LuWr,'(8X,A)') ' ZORA approximation'
          else if (IRELMP == 22) then
            write(LuWr,'(8X,A)') ' ZORA-FP approximation'
          else if (IRELMP == 23) then
            write(LuWr,'(8X,A)') ' IORA approximation'
          end if
        end if
      end if
      write(LuWr,*)
      write(LuWr,'(6X,A)') 'Spectral Resolvent basis set '
      write(LuWr,'(6X,A)') 'Shell  nPrim '
    end if
  end if
  lSh = 0
  do kSh=kShStr,kShEnd
    nExpk = Shells(kSh)%nExp
    if (output .and. (nExpk /= 0) .and. (.not. lOPTO)) write(LuWr,'(9X,A,6X,I2)') AngTp(lSh),nExpk
    lSh = lSh+1
  end do

  ! Auxilliary SO core

  kShStr = dbsc(iCnttp)%iSOC
  kShEnd = kShStr+dbsc(iCnttp)%nSOC-1
  if (output .and. dbsc(iCnttp)%ECP .and. (.not. lOPTO)) then
    type(0) = .false.
    do kSh=kShStr,kShEnd
      nExpk = Shells(kSh)%nExp
      type(0) = type(0) .or. (nExpk /= 0)
    end do
    if (type(0)) then
      if (dbsc(iCnttp)%nOpt /= 0) then
        write(LuWr,*)
        write(LuWr,'(6X,A)') 'Auxilliary core basis'
      end if
      write(LuWr,*)
      write(LuWr,'(6X,A)') 'SOC basis set '
      write(LuWr,'(6X,A)') 'Shell  nPrim '
    end if
  end if
  lSh = 0
  do kSh=kShStr,kShEnd
    nExpk = Shells(kSh)%nExp
    if (output .and. (nExpk /= 0)) write(LuWr,'(9X,A,6X,I2)') AngTp(lSh),nExpk
    lSh = lSh+1
  end do

  if (output .and. (iPrint >= 6)) then
    write(LuWr,*)
    write(LuWr,'(6X,A)') ' Label   Cartesian Coordinates / Bohr'
    write(LuWr,*)
    do iCnt=1,dbsc(iCnttp)%nCntr
      write(LuWr,'(1X,A,1X,3F20.10)') dc(mdc+iCnt)%LblCnt,dbsc(iCnttp)%Coor(1:3,iCnt)
    end do
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (Show) then
  call CollapseOutput(0,'   Basis set information:')
  write(LuWr,*)
end if

return

end subroutine Print_Basis
