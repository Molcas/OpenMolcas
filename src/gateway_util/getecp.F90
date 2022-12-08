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
! Copyright (C) 1990,2020,  Roland Lindh                               *
!               1990, IBM                                              *
!               1993, Per Boussard                                     *
!***********************************************************************
!***********************************************************************
!                                                                      *
!    Objective: To read ECP information, excluding the valence basis-  *
!               set. This means that we read (from input stream) the   *
!               effective charge, the model potential terms (M1 and    *
!               M2), projection parameters and projection orbital      *
!               basis-set (exponents and coefficients)                 *
!                                                                      *
! Called from: Input                                                   *
!                                                                      *
! Calling    : RecPrt                                                  *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!                                                                      *
!     Modified: Per Boussard -93.                                      *
!***********************************************************************

subroutine GetECP(lUnit,iShll,nProj,UnNorm)

use Basis_Info, only: dbsc, nCnttp, Shells
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: lUnit
integer(kind=iwp), intent(inout) :: iShll
integer(kind=iwp), intent(out) :: nProj
logical(kind=iwp), intent(in) :: UnNorm
#include "Molcas.fh"
integer(kind=iwp) :: i, iAng, iEE, ierr, iPP, iPrim, iSS, iStrt, jcr, kcr, mPP(2), n_Elec, n_Occ, nCntrc, ncr, nExpi, nM1, nM2, &
                     nPrim
real(kind=wp) :: ccr, zcr
character(len=180) :: Line
real(kind=wp), allocatable :: Scrt1(:), Scrt2(:)
character(len=180), external :: Get_Ln

!                                                                      *
!***********************************************************************
!                                                                      *
! Process Pseudo Potentials

Line = Get_Ln(lUnit)

if (index(Line,'PP') /= 0) then

  dbsc(nCnttp)%iPP = iShll+1

  !write(u6,*) 'Line=',Line
  if (index(Line,'PPSO') /= 0) then
    call Get_i(4,mPP,2)
  else
    call Get_i(4,mPP,1)
    mPP(2) = 0
  end if
  !write(u6,*) 'mPP=',mPP
  dbsc(nCnttp)%nPP = 1+mPP(1)+mPP(2)

  do iPP=0,mPP(1)
    iShll = iShll+1
    !write(u6,*) 'iPP,dbsc(nCnttp)%nPP=',iPP,dbsc(nCnttp)%nPP
    if (iShll > MxShll) then
      call WarningMessage(2,'Abend in GetECP: Increase MxShll')
      call Abend()
    end if

    ! Pick up the number of terms in the shell
    Line = Get_Ln(lUnit)
    !write(u6,*) 'Line=',Line
    call Get_i1(1,kcr)
    !write(u6,*) 'kcr,iShll=',kcr,iShll
    Shells(iShll)%nExp = 3*kcr
    call mma_allocate(Shells(iShll)%Exp,3*kcr,Label='Exp')

    iStrt = 1
    do jcr=1,kcr
      Line = Get_Ln(lUnit)

      call Get_I1(1,ncr)
      !write(u6,*) 'ncr=',ncr
      Shells(iShll)%Exp(iStrt) = real(ncr,kind=wp)
      iStrt = iStrt+1
      call Get_F1(2,zcr)
      !write(u6,*) 'zcr=',zcr
      Shells(iShll)%Exp(iStrt) = zcr
      iStrt = iStrt+1
      call Get_F1(3,ccr)
      !write(u6,*) 'ccr=',ccr
      Shells(iShll)%Exp(iStrt) = ccr
      iStrt = iStrt+1
    end do

  end do
  do iPP=1,mPP(2)
    iShll = iShll+1
    !write(u6,*) 'iPP,dbsc(nCnttp)%nPP=',iPP,dbsc(nCnttp)%nPP
    if (iShll > MxShll) then
      write(u6,*) 'Abend in GetECP: Increase MxShll'
      call Abend()
    end if

    ! Pick up the number of terms in the shell
    Line = Get_Ln(lUnit)
    !write(u6,*) 'Line=',Line
    call Get_i1(1,kcr)
    !write(u6,*) 'kcr,iShll=',kcr,iShll
    Shells(iShll)%nExp = 3*kcr
    call mma_allocate(Shells(iShll)%Exp,3*kcr,Label='Exp')

    iStrt = 1
    do jcr=1,kcr
      Line = Get_Ln(lUnit)

      call Get_I1(1,ncr)
      !write(u6,*) 'ncr=',ncr
      ncr = ncr+1000
      Shells(iShll)%Exp(iStrt) = real(ncr,kind=wp)
      iStrt = iStrt+1
      call Get_F1(2,zcr)
      !write(u6,*) 'zcr=',zcr
      Shells(iShll)%Exp(iStrt) = zcr
      iStrt = iStrt+1
      call Get_F1(3,ccr)
      !write(u6,*) 'ccr=',ccr
      Shells(iShll)%Exp(iStrt) = ccr
      iStrt = iStrt+1
    end do

  end do
  !write(u6,*) 'Done'

  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! M1 section

!write(u6,*) ' Reading M1'
if (index(Line,'M1') == 0) then
  call WarningMessage(2,'ERROR: Keyword M1 expected, offending line : '//Line)
  call Quit_OnUserError()
end if
Line = Get_Ln(lUnit)
call Get_i1(1,nM1)
dbsc(nCnttp)%nM1 = nM1
if (nM1 > 0) then
  call mma_allocate(dbsc(nCnttp)%M1xp,nM1,Label='dbsc:M1xp')
  call mma_allocate(dbsc(nCnttp)%M1cf,nM1,Label='dbsc:M1cf')
  call Read_v(lUnit,dbsc(nCnttp)%M1xp,1,nM1,1,ierr)
  call Read_v(lUnit,dbsc(nCnttp)%M1cf,1,nM1,1,ierr)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! M2 section

!write(u6,*) ' Reading M2'
Line = Get_Ln(lUnit)
if (index(Line,'M2') == 0) then
  call WarningMessage(2,'ERROR: Keyword M2 expected, offending line : '//Line)
  call Quit_OnUserError()
end if
Line = Get_Ln(lUnit)
call Get_i1(1,nM2)
dbsc(nCnttp)%nM2 = nM2
if (nM2 > 0) then
  call mma_allocate(dbsc(nCnttp)%M2xp,nM2,Label='dbsc:M2xp')
  call mma_allocate(dbsc(nCnttp)%M2cf,nM2,Label='dbsc:M2cf')
  call Read_v(lUnit,dbsc(nCnttp)%M2xp,1,nM2,1,ierr)
  call Read_v(lUnit,dbsc(nCnttp)%M2cf,1,nM2,1,ierr)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Read core-repulsion parameter

!write(u6,*) ' Reading Core-repulsion'
Line = Get_Ln(lUnit)
if (index(Line,'COREREP') == 0) then
  call WarningMessage(2,'ERROR: Keyword COREREP expected, offending line : '//Line)
  call Quit_OnUserError()
end if
Line = Get_Ln(lUnit)
call Get_F1(1,dbsc(nCnttp)%CrRep)
!                                                                      *
!***********************************************************************
!                                                                      *
! Now read Projection parameters

!write(u6,*) ' Reading projection operator'
Line = Get_Ln(lUnit)
if (index(Line,'PROJOP') == 0) then
  call WarningMessage(2,'ERROR: Keyword PROJOP expected, offending line : '//Line)
  call Quit_OnUserError()
end if

! Now read projection basis set
! Loop over each shell type (s,p,d,etc....)
! The ECP is assumed to contain s,p,d and f functions,
! Read(Line,*,Err=993) nProj
Line = Get_Ln(lUnit)
call Get_I1(1,nProj)

if (nProj < 0) return
do iAng=0,nProj
  !write(u6,*) ' iAng=',iAng
  n_Elec = 2*(2*iAng+1)
  iShll = iShll+1
  if (iShll > MxShll) then
    call WarningMessage(2,'Abend in GetECP: Increase MxShll')
    call Quit_OnUserError()
  end if
  !read(Line,*,Err=993) nPrim, nCntrc
  Line = Get_Ln(lUnit)
  call Get_I1(1,nPrim)
  call Get_i1(2,nCntrc)

  Shells(iShll)%nExp = nPrim
  Shells(iShll)%nBasis = nCntrc

  ! Check if occupation number is included on the line

  iSS = 1
  call NxtWrd(Line,iSS,iEE)
  iSS = iEE+1
  call NxtWrd(Line,iSS,iEE)
  iSS = iEE+1
  call NxtWrd(Line,iSS,iEE)
  call mma_allocate(Shells(iShll)%Occ,nCntrc,Label='Occ')
  if (iEE > 0) then
    do i=1,nCntrc
      call Get_i1(2+i,n_Occ)
      Shells(iShll)%Occ(i) = real(n_Occ,kind=wp)/real(n_Elec,kind=wp)
      !write(u6,*) 'n_Occ=',n_Occ
    end do
  else
    Shells(iShll)%Occ(:) = One
  end if

  ! Read "orbital energies"

  !write(u6,*) ' Reading Bk'
  call mma_allocate(Shells(iShll)%Bk,nCntrc,Label='Bk')
  Shells(iShll)%nBk = nCntrc
  if (nCntrc > 0) call Read_v(lUnit,Shells(iShll)%Bk,1,nCntrc,1,ierr)
  if (ierr /= 0) then
    call WarningMessage(2,'Abend in GetBS: Error while reading the exponents')
    call Quit_OnUserError()
  end if

  ! Read gaussian EXPonents

  !write(u6,*) ' Reading Exponents'
  call mma_allocate(Shells(iShll)%Exp,nPrim,Label='Exp')
  Shells(iShll)%nExp = nPrim
  if (nPrim > 0) call Read_v(lUnit,Shells(iShll)%Exp,1,nPrim,1,ierr)
  if (ierr /= 0) then
    call WarningMessage(2,'Abend in GetBS: Error while reading the exponents')
    call Quit_OnUserError()
  end if
  !call RecPrt(' Exponents',Shells(iShll)%nExp,nPrim,1)
  call mma_allocate(Shells(iShll)%Cff_c,nPrim,nCntrc,2,Label='Cff_c')
  call mma_allocate(Shells(iShll)%pCff,nPrim,nCntrc,Label='pCff')
  call mma_allocate(Shells(iShll)%Cff_p,nPrim,nPrim,2,Label='Cff_p')
  Shells(iShll)%Cff_p(:,:,:) = Zero  ! dummy assign
  Shells(iShll)%nBasis = nCntrc

  ! Read contraction coefficients
  ! Observe that the matrix will have nPrim rows and nCntrc columns

  !write(u6,*) ' Reading coefficients'
  do iPrim=1,nPrim
    call Read_v(lUnit,Shells(iShll)%Cff_c,iPrim,nPrim*nCntrc,nPrim,ierr)
    if (ierr /= 0) then
      call WarningMessage(2,'Abend in GetBS: Error while reading the coefficients')
      call Quit_OnUserError()
    end if
  end do

  ! Renormalize

  call mma_allocate(Scrt1,nPrim**2)
  call mma_allocate(Scrt2,nPrim*nCntrc)
  call Nrmlx(Shells(iShll)%Exp,nPrim,Shells(iShll)%Cff_c(1,1,1),nCntrc,Scrt1,nPrim**2,Scrt2,nPrim*nCntrc,iAng)
  call mma_deallocate(Scrt1)
  call mma_deallocate(Scrt2)

  ! Duplicate!

  Shells(iShll)%Cff_c(:,:,2) = Shells(iShll)%Cff_c(:,:,1)

  if (.not. UnNorm) then
    nExpi = Shells(iShll)%nExp
    if (nExpi*Shells(iShll)%nBasis >= 1) then
      call Nrmlz(Shells(iShll)%Exp,nExpi,Shells(iShll)%Cff_c(1,1,1),Shells(iShll)%nBasis,iAng)
    end if
  end if
  Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)

end do

return

end subroutine GetECP
