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
! Calling    : RecPrt, Rdbsl                                           *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!                                                                      *
!     Modified: Per Boussard -93.                                      *
!***********************************************************************

#include "compiler_features.h"
#ifdef _IN_MODULE_

subroutine GetECP(lUnit,iShll,nProj,UnNorm)

use Basis_Info

implicit real*8(A-H,O-Z)
#include "itmax.fh"
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
integer lUnit
integer iShll
integer nProj
logical UnNorm
! Local variables
character(LEN=180) Line, Get_Ln
! External Get_Ln
integer mPP(2)
real*8, dimension(:), allocatable :: Scrt1, Scrt2
!                                                                      *
!***********************************************************************
!                                                                      *
interface
  subroutine Read_v(lunit,Work,istrt,iend,inc,ierr)
    integer lUnit, iStrt, Inc, iErr
    real*8 Work(iend)
  end subroutine Read_v
end interface

!                                                                      *
!***********************************************************************
!                                                                      *
! Process Pseudo Potentials

Line = Get_Ln(lUnit)

if (index(Line,'PP') /= 0) then

  dbsc(nCnttp)%iPP = iShll+1

  !write(6,*) 'Line=',Line
  if (index(Line,'PPSO') /= 0) then
    call Get_i(4,mPP,2)
  else
    call Get_i(4,mPP,1)
    mPP(2) = 0
  end if
  !write(6,*) 'mPP=',mPP
  dbsc(nCnttp)%nPP = 1+mPP(1)+mPP(2)

  do iPP=0,mPP(1)
    iShll = iShll+1
    !write(6,*) 'iPP,dbsc(nCnttp)%nPP=',iPP,dbsc(nCnttp)%nPP
    if (iShll > MxShll) then
      call WarningMessage(2,'Abend in GetECP: Increase MxShll')
      call Abend()
    end if

    ! Pick up the number of terms in the shell
    Line = Get_Ln(lUnit)
    !write(6,*) 'Line=',Line
    call Get_i1(1,kcr)
    !write(6,*) 'kcr,iShll=',kcr,iShll
    Shells(iShll)%nExp = 3*kcr
    call mma_allocate(Shells(iShll)%Exp,3*kcr,Label='Exp')

    iStrt = 1
    do jcr=1,kcr
      Line = Get_Ln(lUnit)

      call Get_I1(1,ncr)
      !write(6,*) 'ncr=',ncr
      Shells(iShll)%exp(iStrt) = dble(ncr)
      iStrt = iStrt+1
      call Get_F1(2,zcr)
      !write(6,*) 'zcr=',zcr
      Shells(iShll)%exp(iStrt) = zcr
      iStrt = iStrt+1
      call Get_F1(3,ccr)
      !write(6,*) 'ccr=',ccr
      Shells(iShll)%exp(iStrt) = ccr
      iStrt = iStrt+1
    end do

  end do
  do iPP=1,mPP(2)
    iShll = iShll+1
    !write(6,*) 'iPP,dbsc(nCnttp)%nPP=',iPP,dbsc(nCnttp)%nPP
    if (iShll > MxShll) then
      call ErrTra()
      write(6,*) 'Abend in GetECP: Increase MxShll'
      call Abend()
    end if

    ! Pick up the number of terms in the shell
    Line = Get_Ln(lUnit)
    !write(6,*) 'Line=',Line
    call Get_i1(1,kcr)
    !write(6,*) 'kcr,iShll=',kcr,iShll
    Shells(iShll)%nExp = 3*kcr
    call mma_allocate(Shells(iShll)%Exp,3*kcr,Label='Exp')

    iStrt = 1
    do jcr=1,kcr
      Line = Get_Ln(lUnit)

      call Get_I1(1,ncr)
      !write(6,*) 'ncr=',ncr
      ncr = ncr+1000
      Shells(iShll)%exp(iStrt) = dble(ncr)
      iStrt = iStrt+1
      call Get_F1(2,zcr)
      !write(6,*) 'zcr=',zcr
      Shells(iShll)%exp(iStrt) = zcr
      iStrt = iStrt+1
      call Get_F1(3,ccr)
      !write(6,*) 'ccr=',ccr
      Shells(iShll)%exp(iStrt) = ccr
      iStrt = iStrt+1
    end do

  end do
  !write(6,*) 'Done'

  return
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! M1 section

!write(6,*) ' Reading M1'
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
  call Read_v(lUnit,dbsc(nCnttp)%M1xp(1),1,nM1,1,ierr)
  call Read_v(lUnit,dbsc(nCnttp)%M1cf(1),1,nM1,1,ierr)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
! M2 section

!write(6,*) ' Reading M2'
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

!write(6,*) ' Reading Core-repulsion'
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

!write(6,*) ' Reading projection operator'
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
  !write(6,*) ' iAng=',iAng
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
      Shells(iShll)%Occ(i) = dble(n_Occ)/dble(n_Elec)
      !write(6,*) 'n_Occ=',n_Occ
    end do
  else
    Shells(iShll)%Occ(:) = One
  end if

  ! Read "orbital energies"

  !write(6,*) ' Reading Bk'
  call mma_allocate(Shells(iShll)%Bk,nCntrc,Label='Bk')
  Shells(iShll)%nBk = nCntrc
  if (nCntrc > 0) call Read_v(lUnit,Shells(iShll)%Bk,1,nCntrc,1,ierr)
  if (ierr /= 0) goto 992

  ! Read gaussian EXPonents

  !write(6,*) ' Reading Exponents'
  call mma_allocate(Shells(iShll)%Exp,nPrim,Label='Exp')
  Shells(iShll)%nExp = nPrim
  if (nPrim > 0) call Read_v(lUnit,Shells(iShll)%Exp,1,nPrim,1,ierr)
  if (ierr /= 0) goto 992
  !call RecPrt(' Exponents',Shells(iShll)%nExp,nPrim,1)
  call mma_allocate(Shells(iShll)%Cff_c,nPrim,nCntrc,2,Label='Cff_c')
  call mma_allocate(Shells(iShll)%pCff,nPrim,nCntrc,Label='pCff')
  call mma_allocate(Shells(iShll)%Cff_p,nPrim,nPrim,2,Label='Cff_p')
  Shells(iShll)%Cff_p(:,:,:) = Zero  ! dummy assign
  Shells(iShll)%nBasis = nCntrc

  ! Read contraction coefficients
  ! Observe that the matrix will have nPrim rows and nCntrc columns

  !write(6,*) ' Reading coefficients'
  do iPrim=1,nPrim
    call Read_v(lUnit,Shells(iShll)%Cff_c(1,1,1),iPrim,nPrim*nCntrc,nPrim,ierr)
    if (ierr /= 0) goto 991
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
!                                                                      *
!***********************************************************************
!                                                                      *
992 continue
call WarningMessage(2,'Abend in GetBS: Error while reading the exponents')
call Quit_OnUserError()
991 continue
call WarningMessage(2,'Abend in GetBS: Error while reading the coefficients')
call Quit_OnUserError()

end subroutine GetECP

#elif !defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(GetECP)

#endif
