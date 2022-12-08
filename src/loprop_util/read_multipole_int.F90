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

subroutine Read_Multipole_Int(lMax,sq_mu,nBas,imu,Ttot,Temp,Origin,rMPq,nElem,nBas1,nBas2,nBasMax,nTemp,nSym,P,Restart,Utility)

use Symmetry_Info, only: Mul
use OneDat, only: sOpSiz
use Data_Structures, only: Alloc1DArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: lMax, nSym, nBas(0:nSym-1), nElem, nBas1, nBas2, nBasMax, nTemp
real(kind=wp), intent(out) :: sq_mu(nBas1**2,0:nElem-1), Temp(nTemp), Origin(3,0:lMax), rMPq(0:nElem-1)
type(Alloc1DArray_Type), intent(out) :: imu(0:nElem-1)
real(kind=wp), intent(in) :: Ttot(nBas2), P(*)
logical(kind=iwp), intent(in) :: Restart, Utility
integer(kind=iwp) :: iCmp, iComp, idum(1), ijSym, iOff, iOffs, iOfft, iOpt0, iOpt1, iRc, iSyLbl, iSym, jSym, l, mu, nComp, nInts, &
                     nInts_Tot, nScr
character(len=8) :: Label
logical(kind=iwp) :: Found
integer(kind=iwp), allocatable :: Comp(:), SyLbl(:)
real(kind=wp), allocatable :: all_ints(:), Scr(:), Tmp(:)
character(len=16), parameter :: RunFile_dLabel = 'LoProp Integrals', &
                                RunFile_iLabel = 'LoProp nInts', &
                                RunFile_iLabel2 = 'LoProp iSyLbl'

!                                                                      *
!***********************************************************************
!                                                                      *
! Note that we will store two sets of the integrals.
!
! Set 1: used to compute the localized properties.
! Set 2: used in the finite difference calculation of the polarizability.
!
! Set 1 will be desymmetrized in case of symmetry.
!
! Pointers to Set 1 are stored in sq_mu.
! Pointers to Set 2 are stored in imu.

nInts_Tot = 0
call mma_allocate(Comp,[0,nElem],label='nComp')
call mma_allocate(SyLbl,[0,nElem],label='SyLbl')
if (Restart) then
  call Qpg_dArray(RunFile_dLabel,Found,nInts_tot)
  if (.not. Found) then
    write(u6,*) 'LoProp Integrals not available on the RunFile.'
    call Abend()
  end if
  call mma_allocate(all_ints,nInts_Tot,label='all_ints')
  call Get_dArray(RunFile_dLabel,all_ints,nInts_tot)
  call Get_iArray(RunFile_iLabel,Comp,nElem)
  call Get_iArray(RunFile_iLabel2,SyLbl,nElem)
end if
nInts = 0
iOpt0 = 0
iOpt1 = ibset(0,sOpSiz)
Label = 'Mltpl  X'
mu = -1
iOff = 1
do l=0,lMax
  nComp = (l+1)*(l+2)/2
  write(Label(8:8),'(I1)') l
  do iComp=1,nComp
#   ifdef _DEBUGPRINT_
    write(u6,*) 'l,iComp=',l,iComp
#   endif
    mu = mu+1
    if (Restart) then
      call mma_allocate(imu(mu)%A,Comp(mu),label='imu')
      nInts = Comp(mu)-4
      call dCopy_(Comp(mu),all_ints(iOff),1,imu(mu)%A,1)
      iSyLbl = SyLbl(mu)
      iOff = iOff+Comp(mu)
    else
      iRc = -1
      iSyLbl = 0
      iCmp = iComp
      call iRdOne(iRc,iOpt1,Label,iCmp,idum,iSyLbl)
      if (iRc /= 0) then
        write(u6,*) 'Polar: error reading length of mu!'
        write(u6,*) 'Mu=',mu
        call Abend()
      end if
      nInts = idum(1)
      call mma_allocate(imu(mu)%A,nInts+4,label='imu')
      call RdOne(iRc,iOpt0,Label,iCmp,imu(mu)%A,iSyLbl)
      if (iRc /= 0) then
        write(u6,*) 'Polar: error reading mu!'
        write(u6,*) 'Mu=',mu
        call Abend()
      end if
      SyLbl(mu) = iSyLbl
      Comp(mu) = nInts+4
      nInts_Tot = nInts_Tot+Comp(mu)
    end if

    ! Transform multipole moment integrals to new basis

    sq_mu(:,mu) = Zero
    if (nSym /= 1) then
      call mma_allocate(Tmp,nBas2,label='Tmp')
      Tmp(:) = Zero
    end if

    ! Now I square the dipole moment integral matrix because it
    ! has been stored as a low triangle

    iOfft = 1
    iOffs = 1
    do iSym=0,nSym-1
      do jSym=0,iSym
        ijSym = Mul(iSym+1,jSym+1)-1
        if (.not. btest(iSyLbl,ijSym)) cycle
        if (nBas(iSym)*nBas(jSym) == 0) cycle
        if (iSym == jSym) then

          ! Now I square the overlap matrix because it has been
          ! stored as a lower triangle

          if (nSym == 1) then
            call Square(imu(mu)%A(iOfft),sq_mu(iOffs,mu),1,nBas(iSym),nBas(iSym))
          else
            call Square(imu(mu)%A(iOfft),Tmp(iOffs),1,nBas(iSym),nBas(iSym))
          end if

          iOfft = iOfft+nBas(iSym)*(nBas(iSym)+1)/2
        else
          if (nSym == 1) then
            call dcopy_(nBas(iSym)*nBas(jSym),imu(mu)%A(iOfft),1,sq_mu(iOffs,mu),1)
          else
            call dcopy_(nBas(iSym)*nBas(jSym),imu(mu)%A(iOfft),1,Tmp(iOffs),1)
          end if
          iOfft = iOfft+nBas(iSym)*nBas(jSym)
        end if
        iOffs = iOffs+nBas(iSym)*nBas(jSym)
      end do
    end do

    if (nSym /= 1) then

      ! Desymmetrize

      nScr = nBasMax*nBas1
      call mma_allocate(Scr,nScr,label='Scr')
      sq_mu(:,mu) = Zero
      call Desymmetrize(Tmp,nBas2,Scr,nScr,sq_mu(:,mu),nBas,nBas1,P,nSym,iSyLbl)
      call mma_deallocate(Scr)
      call mma_deallocate(Tmp)

    end if

    ! Now I transform with my Transformation Matrix Ttot the
    ! multipole moment integral matrices

#   ifdef _DEBUGPRINT_
    call RecPrt('Multipole Integrals in AO Basis',' ',sq_mu(:,mu),nBas1,nBas1)
#   endif
    call TransMu(sq_mu(:,mu),nBas1,Ttot,Temp)
#   ifdef _DEBUGPRINT_
    call RecPrt('Multipole Integrals in LoProp Basis',' ',sq_mu(:,mu),nBas1,nBas1)
#   endif

    ! Pick up the nuclear value for this operator

    rMPq(mu) = imu(mu)%A(nInts+4)
  end do

  ! Pick up the origin of this operator

  Origin(:,l) = imu(mu)%A(nInts+1:nInts+3)

end do
if ((.not. Restart) .and. (.not. Utility)) then
  call mma_allocate(all_ints,nInts_Tot,label='all_ints')
  mu = -1
  iOff = 1
  do l=0,lMax
    nComp = (l+1)*(l+2)/2
    do iComp=1,nComp
      mu = mu+1
      call dCopy_(Comp(mu),imu(mu)%A,1,all_ints(iOff),1)
      iOff = iOff+Comp(mu)
    end do
  end do
  call Put_dArray(RunFile_dLabel,all_ints,nInts_Tot)
  call Put_iArray(RunFile_iLabel,Comp,nElem)
  call Put_iArray(RunFile_iLabel2,SyLbl,nElem)
end if
if (.not. Utility) then
  call mma_deallocate(all_ints)
end if
call mma_deallocate(Comp)
call mma_deallocate(SyLbl)
#ifdef _DEBUGPRINT_
call RecPrt('Origin',' ',Origin,3,lMax+1)
call RecPrt('rMPq',' ',rMPq,1,nElem)
write(u6,*)
write(u6,*) 'Exit  Read_Multipole_Int'
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Read_Multipole_Int
