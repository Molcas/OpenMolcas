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

use Basis_Info, only: dbsc, Gaussian_type, iCnttp_Dummy, mGaussian_type, nCnttp, Nuclear_Model, Point_Charge, Shells
use Center_Info, only: dc
use RICD_Info, only: Thrshld_CD
use DKH_Info, only: iRELMP
use define_af, only: AngTp
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: lOPTO
#include "print.fh"
integer(kind=iwp) :: iCnt, iCnttp, iPrint, iRout, kSh, kShEnd, kShStr, lSh, mdc, nBasisk, nExpk
character(len=4) :: DBas
character :: ChCa, ChCo, ChSph
logical(kind=iwp) :: Output, type(0:7)

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)
if (iPrint == 0) return
!                                                                      *
!***********************************************************************
!                                                                      *
if (Show) then
  write(u6,*)
  call CollapseOutput(1,'   Basis set information:')
  write(u6,'(3X,A)') '   ----------------------'
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
    write(u6,*)
    write(u6,*)
    write(u6,'(6X,A,1X,A)') 'Basis set label:',trim(dbsc(iCnttp)%Bsl)
    if (.not. lOPTO) then
      write(u6,*)
      dbas = dc(mdc+1)%LblCnt(1:4)
      call Upcase(dbas)
      if (dbas == 'DBAS') then
        write(u6,'(6X,A)') 'Diffuse basis set for R-matrix:'
        write(u6,'(6X,A)') '==============================='
        if (dbsc(iCnttp)%nCntr /= 1) then
          call WarningMessage(2,'Too many centers, should only be one!')
          call Quit_OnUserError()
        end if
      else
        if (dbsc(iCnttp)%Aux) then
          write(u6,'(6X,A)') 'Auxiliary basis set:'
          write(u6,'(6X,A)') '=================='
          if (dbsc(iCnttp)%aCD_Thr /= One) then
            write(u6,'(6X,A,G9.2)') 'Threshold in the auxiliary basis set generation is modified to ', &
                                    dbsc(iCnttp)%aCD_Thr*Thrshld_CD
          end if
        else if (dbsc(iCnttp)%Frag) then
          write(u6,'(6X,A)') 'Fragment basis set:'
          write(u6,'(6X,A)') '=================='
        else
          if (dbsc(iCnttp)%fMass == One) then
            write(u6,'(6X,A)') 'Electronic valence basis set:'
            write(u6,'(6X,A)') '------------------'
          else
            write(u6,'(6X,A)') 'Muonic valence basis set:'
            write(u6,'(6X,A)') '------------------'
          end if
        end if
      end if
      if (dbsc(iCnttp)%Fixed) write(u6,'(6X,A)') 'Centers of this basis set are frozen!'
      if (dbsc(iCnttp)%IsMM == 1) then
        write(u6,'(6X,A)') 'This is a MM atom: no basis set'
      else
        if (dbsc(iCnttp)%pChrg) then
          write(u6,'(6X,A,F10.6,A)') 'Associated Effective Charge ',dbsc(iCnttp)%Charge,' au (this is a pseudo charge)'
        else
          write(u6,'(6X,A,F10.6,A)') 'Associated Effective Charge ',dbsc(iCnttp)%Charge,' au'
        end if
        write(u6,'(6X,A,F10.6,A)') 'Associated Actual Charge    ',max(Zero,real(dbsc(iCnttp)%AtmNr,kind=wp)),' au'

        if (Nuclear_Model == Point_Charge) then
          write(u6,'(6X,A)') 'Nuclear Model: Point charge'
        else if (Nuclear_Model == Gaussian_type) then
          write(u6,'(6X,A)') 'Nuclear Model: Finite nucleus - Gaussian distribution'
          write(u6,'(6X,A,ES12.5)') '  Gaussian exponent, Xi/bohr**(-2): ',dbsc(iCnttp)%ExpNuc
        else if (Nuclear_Model == mGaussian_type) then
          write(u6,'(6X,A)') 'Nuclear Model: Finite nucleus - Modified Gaussian distribution'
          write(u6,'(6X,A,ES12.5,A,ES12.5)') '  Parameters, Xi/bohr**(-2), w/bohr**(-2): ',dbsc(iCnttp)%ExpNuc,', ', &
                                             dbsc(iCnttp)%w_mGauss
        else
          call WarningMessage(2,'Illegal Nuclear Model!')
          call Abend()
        end if
      end if
      write(u6,*)
    end if

  end if
  kShStr = dbsc(iCnttp)%iVal
  kShEnd = dbsc(iCnttp)%iVal+dbsc(iCnttp)%nVal-1
  type(0) = .false.
  do kSh=kShStr,kShEnd
    nExpk = Shells(kSh)%nExp
    type(0) = type(0) .or. (nExpk*Shells(kSh)%nBasis /= 0)
  end do
  if (output .and. type(0) .and. (.not. lOPTO)) then
    write(u6,'(6X,A)') 'Shell  nPrim  nBasis  Cartesian Spherical Contaminant'
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
      write(u6,'(9X,A,5X,I3,5X,I3,8X,A,8X,A,8X,A)') AngTp(lSh),nExpk,Shells(kSh)%nBasis,ChCa,ChSph,ChCo
    lSh = lSh+1
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Process PP part, if any.

  kShStr = dbsc(iCnttp)%iPP
  kShEnd = kShStr+dbsc(iCnttp)%nPP-1
  if (output .and. (dbsc(iCnttp)%nPP /= 0) .and. (.not. lOPTO)) then
    write(u6,*)
    write(u6,'(6X,A)') 'Pseudo Potential specification:'
    write(u6,'(6X,A)') '======================================='
    type(0) = .false.
    do kSh=kShStr,kShEnd
      nExpk = Shells(kSh)%nExp/3
      type(0) = type(0) .or. (nExpk /= 0)
    end do
    if (type(0)) then
      write(u6,*)
      write(u6,'(6X,A)') 'Potential  nTerms    '
    end if
  end if
  lSh = 0
  do kSh=kShStr,kShEnd
    nExpk = Shells(kSh)%nExp/3
    !write(u6,*) 'kSh,lSh=',kSh,lSh
    if (output .and. (nExpk /= 0) .and. (.not. lOPTO)) then
      if (lSh == 0) then
        write(u6,'(9X,A,6X,I2)') '  H',nExpk
      else
        write(u6,'(9X,A,6X,I2)') AngTp(lSh-1)//'-H',nExpk
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
    write(u6,*)
    write(u6,'(6X,A)') 'Effective Core Potential specification:'
    write(u6,'(6X,A)') '======================================='
    if (dbsc(iCnttp)%nM1 > 0) then
      write(u6,*)
      write(u6,'(6X,A,I5)') ' Number of M1 terms:',dbsc(iCnttp)%nM1
    end if
    if (dbsc(iCnttp)%nM2 > 0) then
      write(u6,*)
      write(u6,'(6X,A,I5)') ' Number of M2 terms:',dbsc(iCnttp)%nM2
    end if
    type(0) = .false.
    do kSh=kShStr,kShEnd
      nExpk = Shells(kSh)%nExp
      nBasisk = Shells(kSh)%nBasis
      type(0) = type(0) .or. (nExpk*nBasisk /= 0)
    end do
    if (type(0)) then
      write(u6,*)
      write(u6,'(6X,A)') 'Projection basis set '
      write(u6,'(6X,A)') 'Shell  nPrim  nBasis '
    end if
  end if
  lSh = 0
  do kSh=kShStr,kShEnd
    nExpk = Shells(kSh)%nExp
    nBasisk = Shells(kSh)%nBasis
    if (output .and. (nExpk*nBasisk /= 0) .and. (.not. lOPTO)) write(u6,'(9X,A,6X,I2,6X,I2)') AngTp(lSh),nExpk,nBasisk
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
        write(u6,*)
        write(u6,'(6X,A)') 'Spectral Resolvent Operators :'
        if (btest(dbsc(iCnttp)%nOpt,0)) write(u6,'(8X,A)') ' Exchange'
        if (btest(dbsc(iCnttp)%nOpt,1)) write(u6,'(8X,A)') ' Mass-Velocity'
        if (btest(dbsc(iCnttp)%nOpt,2)) write(u6,'(8X,A)') ' Darwin 1-electron contact term'
        if (btest(dbsc(iCnttp)%nOpt,3)) then
          select case (IRELMP)
            case (0)
              write(u6,'(8X,A)') ' No-Pair approximation'
            case (1)
              write(u6,'(8X,A)') ' No-Pair approximation (DK1)'
            case (2)
              write(u6,'(8X,A)') ' No-Pair approximation (DK2)'
            case (3)
              write(u6,'(8X,A)') ' No-Pair approximation (DK3)'
            case (4)
              write(u6,'(8X,A)') ' No-Pair approximation (DK3)'
            case (11)
              write(u6,'(8X,A)') ' RESC approximation'
            case (21)
              write(u6,'(8X,A)') ' ZORA approximation'
            case (22)
              write(u6,'(8X,A)') ' ZORA-FP approximation'
            case (23)
              write(u6,'(8X,A)') ' IORA approximation'
          end select
        end if
      end if
      write(u6,*)
      write(u6,'(6X,A)') 'Spectral Resolvent basis set '
      write(u6,'(6X,A)') 'Shell  nPrim '
    end if
  end if
  lSh = 0
  do kSh=kShStr,kShEnd
    nExpk = Shells(kSh)%nExp
    if (output .and. (nExpk /= 0) .and. (.not. lOPTO)) write(u6,'(9X,A,6X,I2)') AngTp(lSh),nExpk
    lSh = lSh+1
  end do

  ! Auxiliary SO core

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
        write(u6,*)
        write(u6,'(6X,A)') 'Auxiliary core basis'
      end if
      write(u6,*)
      write(u6,'(6X,A)') 'SOC basis set '
      write(u6,'(6X,A)') 'Shell  nPrim '
    end if
  end if
  lSh = 0
  do kSh=kShStr,kShEnd
    nExpk = Shells(kSh)%nExp
    if (output .and. (nExpk /= 0)) write(u6,'(9X,A,6X,I2)') AngTp(lSh),nExpk
    lSh = lSh+1
  end do

  if (output .and. (iPrint >= 6)) then
    write(u6,*)
    write(u6,'(6X,A)') ' Label   Cartesian Coordinates / Bohr'
    write(u6,*)
    do iCnt=1,dbsc(iCnttp)%nCntr
      write(u6,'(1X,A,1X,3F20.10)') dc(mdc+iCnt)%LblCnt,dbsc(iCnttp)%Coor(1:3,iCnt)
    end do
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
if (Show) then
  call CollapseOutput(0,'   Basis set information:')
  write(u6,*)
end if

return

end subroutine Print_Basis
