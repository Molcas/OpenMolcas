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

subroutine WrDisk(rIn,nrIn,jDisp,iIrrep)
! Sorry about this litle subroutine
! The reason is just that right now I want to use AO for SCF
! and MO for RASSCF, this will hopefully be changed, but if you
! see this mess before that I apologize

use McKinley_global, only: ipDisp, ipDisp2, ipDisp3, ipMO, nMethod, RASSCF
use Index_Functions, only: iTri, nTri_Elem
use Basis_Info, only: nBas
use pso_stuff, only: CMO, G1
use Symmetry_Info, only: iOper, nIrrep
use Etwas, only: nAsh, nIsh
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nrIn, jDisp, iIrrep
real(kind=wp), intent(in) :: rIn(nrIn)
#include "print.fh"
integer(kind=iwp) :: iii, iopt, ip(0:7), ip2(0:7), ipCC, ipCM(0:7), ipIn1, ipOut, irc, jAsh, jIrrep, kAsh, kIrrep, nA(0:7), nin, &
                     nIn2, nna
real(kind=wp) :: rDe
integer(kind=iwp) :: n
character(len=8) :: Label
real(kind=wp), allocatable :: Act(:), InAct(:), rOut(:), TempX(:), TempY(:)
integer(kind=iwp), external :: NrOpr
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
if (Show) then
  write(u6,*)
  write(u6,*) '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
  write(u6,*)
  write(u6,*) 'jDisp=',jDisp
end if
nin = 0
nIn2 = 0
nna = 0
ipCC = 1
do jIrrep=0,nIrrep-1
  kIrrep = NrOpr(ieor(iOper(jIrrep),iOper(iIrrep)))
  if (kIrrep < jIrrep) then
    ip(jIrrep) = nIn
    nIn = nIN+nBas(kIrrep)*nBas(jIrrep)
  else if (kIrrep == jIrrep) then
    ip(jIrrep) = nIn
    nIn = nIN+nTri_Elem(nBas(kIrrep))
  end if
  ip2(kIrrep) = nIn2
  nIn2 = nIn2+nBas(kIrrep)*nBas(jIrrep)
  nA(jIrrep) = nnA
  nnA = nnA+nAsh(jIrrep)
  ipCM(jIrrep) = ipCC
  ipCC = ipCC+nBas(jIrrep)**2
# ifdef __INTEL_COMPILER
  ! To avoid error in intel optimization -O3
  if (.false.) write(u6,*) ip(jIrrep)
# endif
end do
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(Act,nIn2,Label='Act')
Act(:) = Zero
call mma_allocate(InAct,nIn2,Label='InAct')
InAct(:) = Zero
call mma_allocate(rOut,nIn2,Label='rOut')
rOut(:) = Zero
call mma_allocate(TempX,nIn2,Label='TempX')
call mma_allocate(TempY,nIn2,Label='TempY')
!                                                                      *
!***********************************************************************
!                                                                      *
! Fock1

if (Show) then
  write(u6,*)
  write(u6,*) 'Fock1'
  write(u6,*)
end if
do jIrrep=0,nIrrep-1
  kIrrep = NrOpr(ieor(iOper(jIrrep),iOper(iIrrep)))
  if ((nBAs(jIrrep) > 0) .and. (nbas(kIrrep) > 0)) then
    if (kIrrep < jIrrep) then
      call DGEMM_('N','N',nBas(jIrrep),nBas(kIrrep),nBas(kIrrep),One,rIn(ipDisp(jDisp)+ip(jIrrep)),nBas(jIrrep), &
                  CMO(ipCM(kIrrep),1),nBas(kIrrep),Zero,TempY,nBas(jIrrep))
      call DGEMM_('T','N',nBas(jIrrep),nBas(kirrep),nBas(jIrrep),One,CMO(ipCM(jIrrep),1),nBas(jIrrep),TempY,nBas(jIrrep),Zero, &
                  Act(1+ip2(kIrrep)),nBas(jIrrep))
      if (Show) then
        write(u6,*)
        write(u6,*) 'ipDisp(jDisp),ip(jIrrep)=',ipDisp(jDisp),ip(jIrrep)
        call RecPrt('ipDisp',' ',rIn(ipDisp(jDisp)+ip(jIrrep)),nBas(jIrrep),nBas(kIrrep))
        write(u6,'(A,G20.10)') 'ipDisp:', &
                               DDot_(nBas(jIrrep)*nBas(kIrrep),rIn(ipDisp(jDisp)+ip(jIrrep)),1,rIn(ipDisp(jDisp)+ip(jIrrep)),1)
        write(u6,'(A,G20.10)') 'ipCM(kIrrep):',DDot_(nBas(kIrrep)*nBas(kIrrep),CMO(ipCM(kIrrep),1),1,CMO(ipCM(kIrrep),1),1)
        write(u6,'(A,G20.10)') 'ipCM(jIrrep):',DDot_(nBas(jIrrep)*nBas(jIrrep),CMO(ipCM(jIrrep),1),1,CMO(ipCM(jIrrep),1),1)
      end if
      call DGetMO(Act(1+ip2(kIrrep)),Nbas(jIrrep),nbas(jIrrep),nBas(kIrrep),Act(1+ip2(jIrrep)),nBas(kIrrep))
    else if (kIrrep == jIrrep) then
      call Square(rIn(ipDisp(jDisp)+ip(jIrrep)),TempX,1,nBas(kirrep),nBas(kirrep))
      call DGEMM_('N','N',nBas(jIrrep),nBas(kIrrep),nBas(kIrrep),One,TempX,nBas(jIrrep),CMO(ipCM(kIrrep),1),nBas(kIrrep),Zero, &
                  TempY,nBas(jIrrep))
      call DGEMM_('T','N',nBas(jIrrep),nBas(kirrep),nBas(jIrrep),One,CMO(ipCM(jIrrep),1),nBas(jIrrep),TempY,nBas(jIrrep),Zero, &
                  Act(1+ip2(jIrrep)),nBas(jIrrep))
    end if
    if (Show) then
      write(u6,*) 'jIrrep,kIrrep=',jIrrep,kIrrep
      write(u6,'(A,G20.10)') 'Act:',DDot_(nIn2,Act,1,Act,1)
    end if
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *

if (nMethod == RASSCF) then

  ! Fock2

  if (Show) then
    write(u6,*)
    write(u6,*) 'Fock2'
    write(u6,*)
  end if
  do jIrrep=0,nIrrep-1
    kIrrep = NrOpr(ieor(iOper(jIrrep),iOper(iIrrep)))
    if ((nBas(jIrrep) > 0) .and. (nBas(kIrrep) > 0)) then
      if (kIrrep == jIrrep) then
        call Square(rIn(ipDisp2(jDisp)+ip(jIrrep)),TempX,1,nBas(kirrep),nBas(kirrep))
        call DGEMM_('N','N',nBas(jIrrep),nBas(kIrrep),nBas(kIrrep),One,TempX,nBas(jIrrep),CMO(ipCM(kIrrep),1),nBas(kIrrep),Zero, &
                    TempY,nBas(jIrrep))
        call DGEMM_('T','N',nBas(jIrrep),nBas(kirrep),nBas(jIrrep),One,CMO(ipCM(jIrrep),1),nBas(jIrrep),TempY,nBas(jIrrep),Zero, &
                    InAct(1+ip2(jIrrep)),nBas(jIrrep))
      else if (kirrep < jirrep) then
        call DGEMM_('N','N',nBas(jIrrep),nBas(kIrrep),nBas(kIrrep),One,rIn(ipDisp2(jDisp)+ip(jIrrep)),nBas(jIrrep), &
                    CMO(ipCM(kIrrep),1),nBas(kIrrep),Zero,TempY,nBas(jIrrep))
        call DGEMM_('T','N',nBas(jIrrep),nBas(kirrep),nBas(jIrrep),One,CMO(ipCM(jIrrep),1),nBas(jIrrep),TempY,nBas(jIrrep),Zero, &
                    InAct(1+ip2(kIrrep)),nBas(jIrrep))
        call DGetMO(InAct(1+ip2(kIrrep)),Nbas(jIrrep),nBas(jIrrep),nBas(kIrrep),InAct(1+ip2(jIrrep)),nBas(kIrrep))
      end if
      if (Show) then
        write(u6,*) 'jIrrep,kIrrep=',jIrrep,kIrrep
        write(u6,'(A,G20.10)') 'InAct:',DDot_(nIn2,InAct,1,InAct,1)
      end if
    end if
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Fock Tot

  if (Show) then
    write(u6,*)
    write(u6,*) 'Fock Tot'
    write(u6,*)
  end if
  iii = 0
  do jIrrep=0,nIrrep-1
    kIrrep = NrOpr(ieor(iOper(jIrrep),iOper(iIrrep)))

    if (nBas(jIrrep)*nIsh(kIrrep) > 0) then
      n = nIsh(kIrrep)*nBas(jIrrep)
      rOut(ip2(kIrrep)+1:ip2(kIrrep)+n) = rOut(ip2(kIrrep)+1:ip2(kIrrep)+n)+ &
                                          Two*(Act(ip2(kIrrep)+1:ip2(kIrrep)+n)+InAct(ip2(kIrrep)+1:ip2(kIrrep)+n))
    end if

    if (nBas(jIrrep) > 0) then
      do jAsh=1,nAsh(kIrrep)
        do kAsh=1,nAsh(kIrrep)
          rDe = G1(iTri(nA(kIrrep)+jAsh,nA(kIrrep)+kAsh),1)
          ipOut = 1+ip2(kIrrep)+nIsh(kIrrep)*nBas(jIrrep)+nBas(jIrrep)*(kAsh-1)
          ipIn1 = 1+ip2(kIrrep)+nBas(jIrrep)*(jAsh-1+nIsh(kIrrep))
          rOut(ipOut:ipOut+nBas(jIrrep)-1) = rOut(ipOut:ipOut+nBas(jIrrep)-1)+rde*InAct(ipIn1:ipIn1+nBas(jIrrep)-1)
        end do
      end do
    end if

    if (nBas(jIrrep)*nAsh(kIrrep) > 0) then
      ipOut = 1+ip2(kIrrep)+nIsh(kIrrep)*nBas(jIrrep)
      if (Show) then
        write(u6,*) 'jIrrep,kIrrep=',jIrrep,kIrrep
        write(u6,'(A,G20.10)') 'ipDisp3:',DDot_(nBas(jIrrep)*nAsh(kIrrep),rIn(ipDisp3(jDisp)+iii),1,rIn(ipDisp3(jDisp)+iii),1)
      end if
      call DGEMM_('T','N',nBas(jIrrep),nAsh(kIrrep),nBas(jIrrep),One,CMO(ipCM(jIrrep),1),nBas(jIrrep),rIn(ipDisp3(jDisp)+iii), &
                  nBas(jIrrep),Zero,TempY,nBas(jIrrep))
      n = nAsh(kIrrep)*nBas(jIrrep)
      rOut(ipOut:ipOut+n-1) = rOut(ipOut:ipOut+n-1)+TempY(1:n)
      iii = iii+nBas(jIrrep)*nAsh(kIrrep)
    end if
#   ifdef __INTEL_COMPILER
    if (.false.) write(u6,*) kIrrep,iii
#   endif
    if (Show) then
      write(u6,*) 'jIrrep,kIrrep=',jIrrep,kIrrep
      write(u6,'(A,G20.10)') 'rOut:',DDot_(nIn2,rOut,1,rOut,1)
    end if
  end do
  if (Show) write(u6,*)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  irc = -1
  iopt = 0
  Label = 'TOTAL'
  call dWrMck(irc,iopt,Label,jdisp,rOut,2**iIrrep)
  if (iRc /= 0) then
    write(u6,*) 'WrDisk: Error writing to MCKINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if
  if (Show) then
    write(u6,'(A,G20.10)') 'TOTAL:',DDot_(nIn2,rOut,1,rOut,1)
  end if

  irc = -1
  iopt = 0
  Label = 'INACTIVE'
  call dWrMck(irc,iopt,Label,jdisp,InAct,2**iIrrep)
  if (iRc /= 0) then
    write(u6,*) 'WrDisk: Error writing to MCKINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if
  if (Show) then
    write(u6,'(A,G20.10)') 'INACTIVE:',DDot_(nIn2,InAct,1,InAct,1)
    write(u6,*)
  end if

  irc = -1
  iopt = 0
  Label = 'MOPERT'
  call dWrMck(irc,iopt,Label,jdisp,rIn(ipMO(jdisp)),2**iIrrep)
  if (iRc /= 0) then
    write(u6,*) 'WrDisk: Error writing to MCKINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! SCF case
  irc = -1
  iopt = 0
  Label = 'TOTAL'
  call dWrMck(irc,iopt,Label,jdisp,Act,2**iIrrep)
  if (iRc /= 0) then
    write(u6,*) 'WrDisk: Error writing to MCKINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if
  if (Show) then
    write(u6,'(A,G20.10)') 'TOTAL:',DDot_(nIn2,Act,1,Act,1)
  end if

end if
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(TempY)
call mma_deallocate(TempX)
call mma_deallocate(rOut)
call mma_deallocate(InAct)
call mma_deallocate(Act)
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine WrDisk
