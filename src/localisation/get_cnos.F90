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

use Localisation_globals, only: ipCMO, ipOcc, MxConstr, nBas, nConstr, nSym
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6, r8

implicit none
integer(kind=iwp), intent(out) :: irc
integer(kind=iwp), intent(in) :: nIF(nSym), nRasO(nSym)
real(kind=wp), intent(out) :: xNrm
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ic1, ic2, iCount, iDaa, iDbb, iOcc, iOff, iOffS(0:8), indxC(16,2,8), ip_Da, ip_Db, ipCorb, ipDaa, ipDbb, &
                     ipMatch, iSym, j, jc, jCount, ji, jOcc, jOff, jSym, k, kbit, kc, kc1, kc2, kOff, l, lc, lc1, lc2, lConstr, &
                     lCount, lOcc_, mAdCMO, mAdCMO_ab, mAdCMOO, mAdOcc_ab, MaxBas, nBB, nBLT, nBT, nSconf
real(kind=wp) :: Etwo, xnorm, xOkk, yOkk
character(len=62) :: Line
logical(kind=iwp) :: DoneCholesky
integer(kind=iwp), external :: Cho_irange
real(kind=r8), external :: ddot_
!***********************************************************************
integer(kind=iwp) :: Match
Match(k,i) = iWork(ipMatch-1+2*(i-1)+k)
!***********************************************************************

!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*

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

call GetMem('Occb','Allo','Real',mAdOcc_ab,nBT)
call GetMem('CMOb','Allo','Real',mAdCMO,2*nBB)
mAdCMO_ab = mAdCMO+nBB
call GetMem('DLT','ALLO','Real',ip_Da,2*nBLT)
ip_Db = ip_Da+nBLT

call GetMem('Match','Allo','Inte',ipMatch,2*MxConstr)
call GetMem('Corb','Allo','Real',ipCorb,MaxBas)

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
      elseif (indxC(j,1,iSym) == 2) then
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
      elseif (indxC(j,2,iSym) == 2) then
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
  iOff = 0
  jOff = 0
  do iSym=1,nSym
    call dcopy_(nBas(iSym)**2,Work(ipCMO+iOff),1,Work(mAdCMO+iOff),1)
    call dcopy_(nBas(iSym)**2,Work(mAdCMO+iOff),1,Work(mAdCMO_ab+iOff),1)
    lOcc_ = ipOcc+jOff+nIF(iSym)
    call dcopy_(nRASO(iSym),Work(lOcc_),1,Work(mAdOcc_ab),1)
    call BestMatch(nConstr(iSym),nRASO(iSym),Work(mAdOcc_ab),iWork(ipMatch),MxConstr)
    do i=1,nConstr(iSym)
      k = Match(1,i)
      jOcc = ipOcc-1+jOff+nIF(iSym)+k
      xOkk = Half*Work(jOcc)
      kc = mAdCMO+iOff+nBas(iSym)*(nIF(iSym)+k-1)
      xNrm = xNrm+ddot_(nBas(iSym),Work(kc),1,Work(kc),1)
      l = Match(2,i)
      iOcc = ipOcc-1+jOff+nIF(iSym)+l
      yOkk = Half*Work(iOcc)
      xnorm = sqrt(abs(xOkk)+abs(yOkk)) !ensures correct normaliz
      lc = mAdCMO+iOff+nBas(iSym)*(nIF(iSym)+l-1)
      xOkk = sqrt(abs(xOkk))/xnorm
      yOkk = sqrt(abs(yOkk))/xnorm
      call dscal_(nBas(iSym),xOkk,Work(kc),1)
      call dscal_(nBas(iSym),yOkk,Work(lc),1)
      call dcopy_(nBas(iSym),Work(lc),1,Work(ipCorb),1)
      call daxpy_(nBas(iSym),One,Work(kc),1,Work(ipCorb),1)
      call daxpy_(nBas(iSym),-One,Work(kc),1,Work(lc),1)
      call dscal_(nBas(iSym),-One,Work(lc),1)
      call dcopy_(nBas(iSym),Work(ipCorb),1,Work(kc),1)
    end do
    jc = 1
    kc = nConstr(iSym)+1
    do i=1,nConstr(iSym)
      l = Match(indxC(i,2,iSym),i)
      lc1 = mAdCMO+iOff+nBas(iSym)*(nIF(iSym)+l-1)
      lc2 = mAdCMO_ab+iOff+nBas(iSym)*(nIF(iSym)+jc-1)
      call dcopy_(nBas(iSym),Work(lc1),1,Work(lc2),1)
      k = Match(indxC(i,1,iSym),i)
      kc1 = mAdCMO+iOff+nBas(iSym)*(nIF(iSym)+k-1)
      kc2 = mAdCMO_ab+iOff+nBas(iSym)*(nIF(iSym)+kc-1)
      call dcopy_(nBas(iSym),Work(kc1),1,Work(kc2),1)
      jc = jc+1
      kc = kc+1
    end do
    kc = nConstr(iSym)+1
    do i=1,nConstr(iSym)
      ic1 = mAdCMO_ab+iOff+nBas(iSym)*(nIF(iSym)+i-1)
      ic2 = mAdCMO+iOff+nBas(iSym)*(nIF(iSym)+kc-1)
      call dcopy_(nBas(iSym),Work(ic1),1,Work(ic2),1)
      kc1 = mAdCMO_ab+iOff+nBas(iSym)*(nIF(iSym)+kc-1)
      kc2 = mAdCMO+iOff+nBas(iSym)*(nIF(iSym)+i-1)
      call dcopy_(nBas(iSym),Work(kc1),1,Work(kc2),1)
      kc = kc+1
    end do
    iOff = iOff+nBas(iSym)**2
    jOff = jOff+nBas(iSym)
  end do

  iOff = 0
  kOff = 0
  do iSym=1,nSym
    ipDaa = ip_Da+kOff
    mAdCMOO = mAdCMO+iOff+nBas(iSym)*nIF(iSym)
    call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nConstr(iSym),One,Work(mAdCMOO),nBas(iSym),Work(mAdCMOO),nBas(iSym),Zero, &
                   Work(ipDaa),nBas(iSym))
    ipDbb = ip_Db+kOff
    mAdCMOO = mAdCMO_ab+iOff+nBas(iSym)*nIF(iSym)
    call DGEMM_tri('N','T',nBas(iSym),nBas(iSym),nConstr(iSym),One,Work(mAdCMOO),nBas(iSym),Work(mAdCMOO),nBas(iSym),Zero, &
                   Work(ipDbb),nBas(iSym))
    do j=1,nBas(iSym)
      do i=1,j-1
        ji = j*(j-1)/2+i
        iDaa = ipDaa-1+ji
        Work(iDaa) = Two*Work(iDaa)
        iDbb = ipDbb-1+ji
        Work(iDbb) = Two*Work(iDbb)
      end do
    end do
    iOff = iOff+nBas(iSym)**2
    kOff = kOff+nBas(iSym)*(nBas(iSym)+1)/2
  end do

  call Get_Etwo_act(Work(ip_Da),Work(ip_Db),nBLT,nBas,nSym,Etwo)

  write(u6,'(1X,A,F12.7,A)') ' Active-Active repulsion : ',Etwo,'  a.u.'
  write(u6,*) ' -------------------------------------------------'
  write(u6,*)
  xNrm = sqrt(xNrm)

end do

call GetMem('Corb','Free','Real',ipCorb,MaxBas)
call GetMem('Match','Free','Inte',ipMatch,2*MxConstr)
call GetMem('DLT','Free','Real',ip_Da,2*nBLT)
call GetMem('CMOb','Free','Real',mAdCMO,2*nBB)
call GetMem('Occb','Free','Real',mAdOcc_ab,nBT)

return

end subroutine Get_CNOs
