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

subroutine INTX(FockI,Temp1,Temp2,Temp3,Temp4,Fock,rMo,loper,idisp)

use Index_Functions, only: iTri
use Symmetry_Info, only: Mul
use MCLR_Data, only: CMO, DspVec, G1t, ipCM, ipMat, ipMatLT, nA, nB, nDens, SWLbl
use input_mclr, only: iMethod, nAsh, nBas, nIsh, nOrb, nSym, nTPert
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
real(kind=wp), intent(out) :: FockI(nDens), Temp1(nDens), Temp2(nDens), Temp3(nDens), Temp4(nDens), Fock(nDens)
real(kind=wp), intent(_OUT_) :: rMO(*)
integer(kind=iwp), intent(in) :: lOper, iDisp
integer(kind=iwp) :: i, iOff, iOp, iOpt, iRC, iS, j, jDisp, jOff, jS
real(kind=wp) :: rde
character(len=8) :: Label

!***********************************************************************
if (nDens == 0) return
jDisp = DspVec(idisp)

!       x
! Read P
!
! Remember two-electron contributions are saved in MO base
! but one electron contributions is saved in AO base.

if (btest(ntpert(idisp),2)) then   ! 2 el contribution
  if (btest(ntpert(idisp),4)) then  ! from mckinley
    if (iMethod == 2) then

      !----------------------------------------------------------------*
      !
      ! RASSCF
      !
      !----------------------------------------------------------------*

      Label = 'TOTAL'
      iop = ibset(0,loper)
      irc = -1
      iopt = 0
      call dRdMck(iRC,iOpt,Label,jDisp,Fock,iop)
      if (iRc /= 0) then
        write(u6,*) 'IntX: Error reading MCKINT'
        write(u6,'(A,A)') 'Label=',Label
        call Abend()
      end if
      call ReLoad(Fock,loper+1,nbas,norb)

      Label = 'INACTIVE'
      iop = ibset(0,loper)
      irc = -1
      iopt = 0
      call dRdMck(iRC,iOpt,Label,jDisp,Focki,iop)
      if (iRc /= 0) then
        write(u6,*) 'IntX: Error reading MCKINT'
        write(u6,'(A,A)') 'Label=',Label
        call Abend()
      end if
      call ReLoad(Focki,loper+1,nbas,norb)

      Label = 'MOPERT'
      iop = ibset(0,loper)
      irc = -1
      iopt = 0
      call dRdMck(iRC,iOpt,Label,jDisp,rMO,iop)
      if (iRc /= 0) then
        write(u6,*) 'IntX: Error reading MCKINT'
        write(u6,'(A,A)') 'Label=',Label
        call Abend()
      end if
    else

      !----------------------------------------------------------------*
      !
      ! SCF
      !
      !----------------------------------------------------------------*

      Label = 'TOTAL'
      iop = ibset(0,loper)
      irc = -1
      iopt = 0
      call dRdMck(iRC,iOpt,Label,jDisp,Focki,iop)
      if (iRc /= 0) then
        write(u6,*) 'IntX: Error reading MCKINT'
        write(u6,'(A,A)') 'Label=',Label
        call Abend()
      end if
      Fock(:) = Zero
      do iS=1,nSym
        js = Mul(is,loper+1)
        Fock(ipMat(is,js):ipMat(is,js)+nOrb(is)*nIsh(js)-1) = Two*Focki(ipMat(is,js):ipMat(is,js)+nOrb(is)*nIsh(js)-1)
      end do
    end if
  end if
end if

!----------------------------------------------------------------------*
!
! Two electron fock matrix done!
! Lets fix  the one electron matrixes
!
!----------------------------------------------------------------------*

if (btest(nTPert(iDisp),1)) then ! 1 el contribution
  iop = ibset(0,loper)
  if (btest(ntpert(idisp),4)) then  ! from mckinley
    Label = 'ONEGRD'
    irc = -1
    iopt = 0
    call dRdMck(iRC,iOpt,Label,jDisp,Temp1,iop)
    if (iRc /= 0) then
      write(u6,*) 'IntX: Error reading MCKINT'
      write(u6,'(A,A)') 'Label=',Label
      call Abend()
    end if
  else                                         ! or seward
    Label = SwLbl(idisp)
    iopt = 0
    irc = -1
    call RdOne(irc,iopt,Label,jDisp,temp1,iop)
    if (iRc /= 0) then
      write(u6,*) 'IntX: Error reading MCKINT'
      write(u6,'(A,A)') 'Label=',Label
      call Abend()
    end if
  end if
end if

Temp2(:) = Zero

do iS=1,nSym
  do jS=1,is
    if (nBas(is)*nBas(js) /= 0) then
      if (Mul(iS,jS) == loper+1) then
        if (is == js) then
          call Square(Temp1(ipMatLt(iS,jS)),Temp4,1,nBas(is),nBas(is))
        else
          Temp4(1:nBas(iS)*nBas(jS)) = Temp1(ipMatLT(iS,Js):ipMatLT(iS,Js)+nBas(iS)*nBas(jS)-1)
        end if
        call DGEMM_('T','N',nOrb(iS),nBas(jS),nBAs(iS),One,CMO(ipCM(iS)),nBas(is),Temp4,nBas(iS),Zero,Temp3,nOrb(iS))
        call DGEMM_('N','N',nOrb(is),nB(jS),nBas(jS),One,Temp3,nOrb(iS),CMO(ipCM(jS)),nBas(jS),Zero,Temp2(ipMat(iS,jS)),nOrb(iS))
        if (is /= js) then
          call DGEMM_('T','T',nOrb(jS),nOrb(iS),nBAs(jS),One,CMO(ipCM(jS)),nBas(js),Temp4,nBas(iS),Zero,Temp3,nOrb(jS))
          call DGEMM_('N','N',nOrb(js),nB(iS),nBas(iS),One,Temp3,nOrb(jS),CMO(ipCM(iS)),nBas(iS),Zero,Temp2(ipMat(jS,iS)),nOrb(jS))
        end if

      end if
    end if
  end do
end do

Temp3(:) = Zero
do iS=1,nSym
  js = Mul(is,loper+1)
  if (nOrb(js) < 1) cycle
  do j=1,nAsh(is)+nish(is)
    do i=1,nAsh(is)+nIsh(is)
      if ((i == j) .and. (i <= nish(is)) .and. (j <= nish(is))) then
        rde = Two
      else if ((i > nish(is)) .and. (j > nish(is))) then
        rde = G1t(iTri(i-nIsh(is)+nA(is),j-nIsh(is)+nA(is)))
      else
        cycle
      end if
      iOff = ipMat(js,is)+(i-1)*nOrb(js)
      jOff = ipMat(js,is)+(j-1)*nOrb(js)
      Temp3(iOff:iOff+nOrb(js)-1) = Temp3(iOff:iOff+nOrb(js)-1)+rDe*Temp2(jOff:jOff+nOrb(js)-1)
    end do
  end do
end do
!***********************************************************************

!----------------------------------------------------------------------*
!
! One electron fock matrix done!
!
!----------------------------------------------------------------------*

if (btest(ntpert(idisp),2)) then
  Focki(:) = Focki(:)+Temp2(:)
  Fock(:) = Fock(:)+Temp3(:)
else
  FockI(:) = Temp2(:)
  Fock(:) = Temp3(:)
end if

end subroutine INTX
