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
! Copyright (C) 1991, Roland Lindh                                     *
!***********************************************************************

subroutine Prpt_old(nirrep,nbas,ndim,n2dim,vec,occ,maxscr,scr)
!***********************************************************************
!                                                                      *
! Purpose: calculation of expectation values of different              *
!          operators as available on the 'ONEINT' file                 *
!                                                                      *
! Caution: before calling this subroutine one needs to                 *
!          open the ONEINT file                                        *
!                                                                      *
! Calling parameters:                                                  *
!                                                                      *
!   nIrRep            number of irreducible representations            *
!   nBas(0:nIrRep)    number of basis functions in each repre-         *
!                     sentation                                        *
!   ndim              total number of basis functions                  *
!   n2dim             sum(i,i=0,nirrep-1)(nbas(i)**2): elements        *
!                     of all vectors in all representations            *
!   vec(1:n2dim)      eigenvectors, all for each representation        *
!   occ(1:ndim)       occupation number for all eigenvectors           *
!   maxscr            the maximum available size of the                *
!                     scratch area                                     *
!   scr(1:maxscr)     a scratch area whose size, nscr, can be          *
!                     calculated in the following way:                 *
!                                                                      *
!                     nscr =sum(i,i=0,nirrep-1)(nbas(i)*(nbas(i)+1)/2) *
!                           +3                                         *
!                           +3                                         *
!                           +nComp                                     *
!                           +nComp                                     *
!                           +2*nComp                                   *
!                           +ndim*(ndim+1)/2+4                         *
!                           +2*ntComp*(ntComp+1)                       *
!                                                                      *
!                     and should not exceed maxscr.                    *
!                     nComp is the number of cartresian components     *
!                     for the given operator                           *
!                     ntComp is the number of components of            *
!                     the l-th moment opartors which are transformed   *
!                     into l-pole moment. Currently,                   *
!                     ntComp=15 (hexadecpole moments)                  *
!                                                                      *
! 1991 R. Lindh, Dept. of Theor. Chem. Univ. of Lund, Sweden.          *
!***********************************************************************

use iso_c_binding, only: c_f_pointer, c_loc
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nirrep, nbas(0:nirrep), ndim, n2dim, maxscr
real(kind=wp), intent(in) :: vec(n2dim), occ(ndim)
real(kind=wp), intent(out) :: scr(maxscr)
integer(kind=iwp) :: iadC1, iadC2, iadDen, iadEl, iadLab, iadNuc, iadopr, iadtmp, iadtmt, iCount, idum(1), iOcc, iopt, irc, iscr, &
                     iSmLbl, iVec, jCount, jRC, maxCen, maxGG, mDim, nblock, nComp, nfblock, n_Int, nscr
real(kind=wp) :: dummy
logical(kind=iwp) :: short, NxtOpr, ifallorb
character(len=8) :: label

call Prpt_old_Internal(Scr)

return

! This is to allow type punning without an explicit interface
contains

subroutine Prpt_old_Internal(Scr)

  real(kind=wp), target :: Scr(*)
  character, pointer :: cScr(:)
  integer(kind=iwp) :: i, iComp, iEF, ii, il, ir, j, k

  write(u6,*)
  call CollapseOutput(1,'   Molecular properties:')
  write(u6,'(3X,A)') '   ---------------------'
  write(u6,*)
  short = .true.
  ifallorb = .false.
  mDim = 1
  iadopr = -1 ! dummy initialize
  iadtmt = -1 ! dummy initialize
  iadtmp = -1 ! dummy initialize

  nblock = 0
  do i=0,nirrep-1
    nblock = nblock+nbas(i)*(nbas(i)+1)/2
  end do
  nfblock = ndim*(ndim+1)/2

  ! estimate the minimum scratch area needed for calculatig
  ! the average values of a 3 component operator

  nscr = nblock+nfblock
  iscr = nscr+76

  if (iscr > maxscr) Go To 999

  iadDen = 1

  ! calculate the density matrix with all off-diagonal elements
  ! multipled by 2

  call dcopy_(nblock,[Zero],0,Scr(iadDen),1)
  iCount = iadDen
  iVec = 0
  iOcc = 0
  do i=0,nIrrep-1
    do ii=1,nBas(i)
      iOcc = iOcc+1
      jCount = iCount
      do il=1,nBas(i)
        do ir=1,il-1
          Scr(jCount) = Scr(jCount)+Two*Vec(iVec+il)*Vec(iVec+ir)*Occ(iOcc)
          jCount = jCount+1
        end do
        Scr(jCount) = Scr(jCount)+Vec(iVec+il)*Vec(iVec+il)*Occ(iOcc)
        jCount = jCount+1
      end do
      iVec = iVec+nBas(i)
    end do
    iCount = iCount+nbas(i)*(nBas(i)+1)/2
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Scan the ONEINT file for multipole moment operators

  iadC1 = iadDen+nblock
  iadC2 = iadC1+3
  iadNuc = iadC2+3

  !write(u6,*) ' Starting scan of ONEINT for multipole moments'
  do i=1,99
    NxtOpr = .false.
    nComp = (i+1)*(i+2)/2

    iadEl = iadNuc+nComp
    iadLab = iadEl+nComp
    iscr = nscr+10+4*nComp
    if (iscr > maxscr) Go To 999

    call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
    call dcopy_(nComp,[Zero],0,Scr(iadEl),1)
    call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
    write(label,'(a,i2)') 'MLTPL ',i
    do iComp=1,nComp
      irc = -1
      iopt = 1
      call iRdOne(irc,iopt,label,iComp,idum,iSmLbl)
      if (irc == 0) n_Int = idum(1)
      if (irc /= 0) go to 101
      NxtOpr = .true.
      irc = -1
      iopt = 0
      iadOpr = iadLab+2*nComp
      call RdOne(irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
      if (irc /= 0) go to 101
      if (n_Int /= 0) call CmpInt(scr(iadOpr),n_Int,nBas,nIrrep,iSmLbl)
      scr(iadNuc+icomp-1) = scr(iadOpr+n_Int+3)
      if (iComp == 1) then
        do k=0,2
          scr(iadC1+k) = scr(iadOpr+n_Int+k)
          scr(iadC2+k) = scr(iadOpr+n_Int+k)
        end do
      end if
      if (n_Int == 0) Go To 101
      call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen),nDim,Occ,nblock,Scr(iadOpr),Scr(iadEl+iComp-1))
  101 continue
    end do
    if (.not. NxtOpr) Go To 199
    iadTmt = iadOpr+nblock
    if (i <= 4) then
      iadTmp = iadTmt+nComp**2
      iscr = iadTmp+nComp
      if (iscr > maxscr) Go To 999
    else
      iadTmp = iadTmt
    end if

    call c_f_pointer(c_loc(Scr(iadLab)),cScr,[1])
    call prop(short,label,scr(iadC1),scr(iadC2),nirrep,nBas,mDim,occ,dummy,scr(iadEl),scr(iadNuc),i,cScr,scr(iadTmt),scr(iadTmp), &
              ifallorb)
    nullify(cScr)
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Scan 'ONEINT' for electric field integrals

  199 continue

  !write(u6,*) ' Starting scan of ONEINT for various elec. field integrals'

  do iEF=0,2
    nComp = (iEF+1)*(iEF+2)/2

    iadEl = iadNuc+nComp
    iadLab = iadEl+nComp

    ! loop over differnt operator origins (max.9999)

    maxCen = 9999
    do i=1,maxCen
      call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
      call dcopy_(nComp,[Zero],0,Scr(iadEl),1)
      call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
      write(label,'(a,i1,i5)') 'EF',iEF,i
      NxtOpr = .false.
      do iComp=1,nComp
        irc = -1
        iopt = 1
        call iRdOne(irc,iopt,label,iComp,idum,iSmLbl)
        if (irc == 0) n_Int = idum(1)
        if (irc /= 0) go to 201
        NxtOpr = .true.
        irc = -1
        iopt = 0
        iadOpr = iadLab+2*nComp
        call RdOne(irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
        if (irc /= 0) go to 201
        if (n_Int /= 0) call CmpInt(scr(iadOpr),n_Int,nBas,nIrrep,iSmLbl)
        scr(iadNuc+icomp-1) = scr(iadOpr+n_Int+3)
        if (iComp == 1) then
          do k=0,2
            scr(iadC1+k) = scr(iadOpr+n_Int+k)
            scr(iadC2+k) = scr(iadOpr+n_Int+k)
          end do
        end if
        if (n_Int == 0) Go To 201
        call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen),nDim,Occ,nblock,Scr(iadOpr),Scr(iadEl+iComp-1))
  201   continue
      end do
      if (.not. NxtOpr) go to 299

      call c_f_pointer(c_loc(Scr(iadLab)),cScr,[1])
      call prop(short,label,scr(iadC1),scr(iadC2),nirrep,nBas,mDim,occ,dummy,scr(iadEl),scr(iadNuc),i,cScr,scr(iadTmt), &
                scr(iadTmp),ifallorb)
      nullify(cScr)
    end do

  299 continue
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !write(u6,*) ' Starting scan of ONEINT diamagnetic shielding'

  nComp = 9

  iadEl = iadNuc+nComp
  iadLab = iadEl+nComp

  maxGG = 99
  maxCen = 99
  ! loop over differnt gauge origins (max.99)
  do j=1,maxGG
    call dcopy_(nComp,[Zero],0,Scr(iadNuc),1)
    call dcopy_(nComp,[Zero],0,Scr(iadEl),1)
    call dcopy_(2*nComp,[Zero],0,Scr(iadLab),1)
    jRC = 0
    ! loop over differnt operator origins (max.99)
    do i=1,maxCen
      write(label,'(a,i2,i2)') 'DMS ',j,i
      NxtOpr = .false.
      do iComp=1,nComp
        irc = -1
        iopt = 1
        call iRdOne(irc,iopt,label,iComp,idum,iSmLbl)
        if (irc == 0) n_Int = idum(1)
        if (irc /= 0) go to 402
        NxtOpr = .true.
        irc = -1
        iopt = 0
        iadOpr = iadLab+2*nComp
        call RdOne(irc,iopt,label,iComp,scr(iadOpr),iSmLbl)
        if (irc /= 0) go to 402
        if (n_Int /= 0) call CmpInt(scr(iadOpr),n_Int,nBas,nIrrep,iSmLbl)
        scr(iadNuc+icomp-1) = scr(iadOpr+n_Int+3)
        if (iComp == 1) then
          do k=0,2
            scr(iadC1+k) = scr(iadOpr+n_Int+k)
          end do
        end if
        if (iComp == 2) then
          do k=0,2
            scr(iadC2+k) = scr(iadOpr+n_Int+k)
          end do
        end if
        if (n_Int == 0) Go To 402
        call Xprop(short,ifallorb,nIrrep,nBas,nBlock,Scr(iadDen),nDim,Occ,nblock,Scr(iadOpr),Scr(iadEl+iComp-1))
  402   continue
      end do
      if (.not. NxtOpr) Go To 4000

      call c_f_pointer(c_loc(Scr(iadLab)),cScr,[1])
      call prop(short,label,scr(iadC1),scr(iadC2),nirrep,nBas,mDim,occ,dummy,scr(iadEl),scr(iadNuc),i,cScr,scr(iadTmt), &
                scr(iadTmp),ifallorb)
      nullify(cScr)
      jRC = 1
    end do
  4000 continue
    if (jRC == 0) Go To 499
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  499 continue

  call CollapseOutput(0,'   Molecular properties:')
  write(u6,*)
  return

  999 continue
  write(u6,'(//1x,a/1x,a/1x,a,i8,a,i8)') ' Warrning:',' Not enough scratch area to perform calculations',' Needed at least:',iscr, &
                                        '   available:',maxscr

  call CollapseOutput(0,'   Molecular properties:')
  write(u6,*)

  return

end subroutine Prpt_old_Internal

end subroutine Prpt_old
