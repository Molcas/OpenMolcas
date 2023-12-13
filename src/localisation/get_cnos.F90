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

subroutine Get_CNOs(irc,nIF,nRASO,xNrm)
!***********************************************************************
!                                                                      *
!     purpose: Generate constrained orbitals from INPORB               *
!              Analysis of the spin configurations                     *
!                                                                      *
!***********************************************************************

use Localisation_globals, only: CMO, MxConstr, nBas, nConstr, nSym, Occ
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nIF(nSym), nRasO(nSym)
real(kind=wp), intent(out) :: xNrm
integer(kind=iwp) :: i, ic1, ic2, iCount, iDab, iOcc, iOff, iOffS(0:8), indxC(16,2,8), ipDab, iSym, j, jc, jCount, ji, jOcc, jOff, &
                     jSym, k, kbit, kc, kc1, kc2, kOff, l, lc, lc1, lc2, lConstr, lCount, lOcc_, mAdCMOO, MaxBas, nBB, nBLT, nBT, &
                     nSconf
real(kind=wp) :: Etwo, xnorm, xOkk, yOkk
character(len=62) :: Line
logical(kind=iwp) :: DoneCholesky
integer(kind=iwp), allocatable :: Match(:,:)
real(kind=wp), allocatable :: CMO_(:), CMO_ab(:), Corb(:), Da(:), Db(:), Occ_ab(:)
integer(kind=iwp), external :: Cho_irange
real(kind=wp), external :: ddot_

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

call Untested('Get_CNOs')

call DecideonCholesky(DoneCholesky)
if (.not. DoneCholesky) then
  write(u6,*) '*** Constrained NOs implemented only with CD or RI.'
  write(u6,*) '*** Use Cholesky or RICD in Seward and rerun! *****'
  call Abend()
end if

irc = 0
nSconf = 1
iOffS(0) = 0
nBB = 0
nBT = 0
nBLT = 0
MaxBas = 0
lConstr = 0
do iSym=1,nSym
  iOffS(iSym) = iOffS(iSym-1)+nConstr(iSym)
  do j=1,nConstr(iSym)
    nSconf = 2*nSconf
  end do
  lConstr = lConstr+nConstr(iSym)
  nBB = nBB+nBas(iSym)**2
  nBT = nBT+nBas(iSym)
  nBLT = nBLT+nBas(iSym)*(nBas(iSym)+1)/2
  MaxBas = max(MaxBas,nBas(iSym))
end do
write(u6,'(A,I6)') ' Total number of spin configurations: ',nSconf
write(u6,*)

call mma_allocate(Occ_ab,nBT,label='Occb')
call mma_allocate(CMO_,nBB,label='CMOb')
call mma_allocate(CMO_ab,nBB,label='CMOb')
call mma_allocate(Da,nBLT,label='DLT')
call mma_allocate(Db,nBLT,label='DLT')

call mma_allocate(Match,2,MxConstr,label='Match')
call mma_allocate(Corb,MaxBas,label='Corb')

do iCount=0,nSconf-1
  do jCount=0,lConstr-1
    jSym = Cho_Irange(jCount+1,iOffS,nSym,.true.)
    lCount = jCount-iOffS(jSym-1)+1
    kbit = ibits(iCount,jCount,1)
    if (kbit == 0) then
      indxC(lCount,1,jSym) = 1
      indxC(lCount,2,jSym) = 2
    else
      indxC(lCount,1,jSym) = 2
      indxC(lCount,2,jSym) = 1
    end if
  end do
  write(u6,*) ' -------------------------------------------------'
  write(u6,*) ' Configuration of the constrained spins (up/down) '
  write(u6,*) ' -------------------------------------------------'
  write(u6,'(1X,A,I3)') ' nr ',iCount+1
  do iSym=1,nSym
    write(u6,'(1X,A,I1)') ' sym: ',iSym
    Line(1:14) = '         (+) '
    k = 15
    do j=1,nConstr(iSym)
      if (indxC(j,1,iSym) == 1) then
        Line(k:k+2) = ' u '
      else if (indxC(j,1,iSym) == 2) then
        Line(k:k+2) = ' d '
      else
        Line(k:k+2) = '   '
      end if
      k = k+3
    end do
    write(u6,*) Line(1:k-1)
    Line(1:14) = '         (-) '
    k = 15
    do j=1,nConstr(iSym)
      if (indxC(j,2,iSym) == 1) then
        Line(k:k+2) = ' u '
      else if (indxC(j,2,iSym) == 2) then
        Line(k:k+2) = ' d '
      else
        Line(k:k+2) = '   '
      end if
      k = k+3
    end do
    write(u6,*) Line(1:k-1)
  end do
  write(u6,*) ' -------------------------------------------------'

  xNrm = Zero
  iOff = 1
  jOff = 0
  do iSym=1,nSym
    call dcopy_(nBas(iSym)**2,CMO(iOff),1,CMO_(iOff),1)
    call dcopy_(nBas(iSym)**2,CMO(iOff),1,CMO_ab(iOff),1)
    lOcc_ = jOff+nIF(iSym)+1
    call dcopy_(nRASO(iSym),Occ(lOcc_),1,Occ_ab,1)
    call BestMatch(nConstr(iSym),nRASO(iSym),Occ_ab,Match,MxConstr)
    do i=1,nConstr(iSym)
      k = Match(1,i)
      jOcc = jOff+nIF(iSym)+k
      xOkk = Half*Occ(jOcc)
      kc = iOff+nBas(iSym)*(nIF(iSym)+k-1)
      xNrm = xNrm+ddot_(nBas(iSym),CMO_(kc),1,CMO_(kc),1)
      l = Match(2,i)
      iOcc = jOff+nIF(iSym)+l
      yOkk = Half*Occ(iOcc)
      xnorm = sqrt(abs(xOkk)+abs(yOkk)) !ensures correct normaliz
      lc = iOff+nBas(iSym)*(nIF(iSym)+l-1)
      xOkk = sqrt(abs(xOkk))/xnorm
      yOkk = sqrt(abs(yOkk))/xnorm
      call dscal_(nBas(iSym),xOkk,CMO_(kc),1)
      call dscal_(nBas(iSym),yOkk,CMO_(lc),1)
      call dcopy_(nBas(iSym),CMO_(lc),1,Corb,1)
      call daxpy_(nBas(iSym),One,CMO_(kc),1,Corb,1)
      call daxpy_(nBas(iSym),-One,CMO_(kc),1,CMO_(lc),1)
      call dscal_(nBas(iSym),-One,CMO_(lc),1)
      call dcopy_(nBas(iSym),Corb,1,CMO_(kc),1)
    end do
    jc = 1
    kc = nConstr(iSym)+1
    do i=1,nConstr(iSym)
      l = Match(indxC(i,2,iSym),i)
      lc1 = iOff+nBas(iSym)*(nIF(iSym)+l-1)
      lc2 = iOff+nBas(iSym)*(nIF(iSym)+jc-1)
      call dcopy_(nBas(iSym),CMO_(lc1),1,CMO_ab(lc2),1)
      k = Match(indxC(i,1,iSym),i)
      kc1 = iOff+nBas(iSym)*(nIF(iSym)+k-1)
      kc2 = iOff+nBas(iSym)*(nIF(iSym)+kc-1)
      call dcopy_(nBas(iSym),CMO_(kc1),1,CMO_ab(kc2),1)
      jc = jc+1
      kc = kc+1
    end do
    kc = nConstr(iSym)+1
    do i=1,nConstr(iSym)
      ic1 = iOff+nBas(iSym)*(nIF(iSym)+i-1)
      ic2 = iOff+nBas(iSym)*(nIF(iSym)+kc-1)
      call dcopy_(nBas(iSym),CMO_ab(ic1),1,CMO_(ic2),1)
      kc1 = iOff+nBas(iSym)*(nIF(iSym)+kc-1)
      kc2 = iOff+nBas(iSym)*(nIF(iSym)+i-1)
      call dcopy_(nBas(iSym),CMO_ab(kc1),1,CMO_(kc2),1)
      kc = kc+1
    end do
    iOff = iOff+nBas(iSym)**2
    jOff = jOff+nBas(iSym)
  end do

  iOff = 1
  kOff = 1
  do iSym=1,nSym
    ipDab = kOff
    mAdCMOO = iOff+nBas(iSym)*nIF(iSym)
    call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nConstr(iSym),One,CMO_(mAdCMOO),nBas(iSym),CMO_(mAdCMOO),nBas(iSym),Zero, &
                   Da(ipDab),nBas(iSym))
    call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nConstr(iSym),One,CMO_ab(mAdCMOO),nBas(iSym),CMO_ab(mAdCMOO),nBas(iSym),Zero, &
                   Db(ipDab),nBas(iSym))
    do j=1,nBas(iSym)
      do i=1,j-1
        ji = j*(j-1)/2+i
        iDab = ipDab-1+ji
        Da(iDab) = Two*Da(iDab)
        Db(iDab) = Two*Db(iDab)
      end do
    end do
    iOff = iOff+nBas(iSym)**2
    kOff = kOff+nBas(iSym)*(nBas(iSym)+1)/2
  end do

  call Get_Etwo_act(Da,Db,nBLT,nBas,nSym,Etwo)

  write(u6,'(1X,A,F12.7,A)') ' Active-Active repulsion : ',Etwo,'  a.u.'
  write(u6,*) ' -------------------------------------------------'
  write(u6,*)
  xNrm = sqrt(xNrm)

end do

call mma_deallocate(Occ_ab)
call mma_deallocate(CMO_)
call mma_deallocate(CMO_ab)
call mma_deallocate(Da)
call mma_deallocate(Db)
call mma_deallocate(Match)
call mma_deallocate(Corb)

return

end subroutine Get_CNOs
