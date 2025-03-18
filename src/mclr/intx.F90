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

subroutine INTX(FockI,Temp1,Temp2,Temp3,Temp4,Fock,rMo,loper,idisp,r)

use Arrays, only: G1t, CMO
use MCLR_Data, only: nDens2, ipCM, ipMat, ipMatLT, nA, nB, nDens
use MCLR_Data, only: DspVec, SWLbl
use input_mclr, only: iMethod, nSym, nAsh, nBas, nIsh, nOrb, nTPert

implicit none
real*8 FockI(nDens2), Temp2(ndens2), Temp3(nDens2), Temp4(ndens2), Temp1(nDens2), Fock(nDens2), rMO(*)
integer lOper, iDisp
real*8 r
character(len=8) Label
integer jDisp, iOp, iRC, iOpt, iS, jS
real*8 rde
integer i, j, iTri
! Statement function
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)

!***********************************************************************
if (nDens2 == 0) return
jDisp = DspVec(idisp)

!       x
! Read P
!
! Remember two-electron contributions are saved in MO base
! but one electron contributions is saved in AO base.

if (iand(ntpert(idisp),2**2) == 4) then   ! 2 el contribution
  if (iand(ntpert(idisp),2**4) == 2**4) then  ! from mckinley
    if (iMethod == 2) then

      !----------------------------------------------------------------*
      !
      ! RASSCF
      !
      !----------------------------------------------------------------*

      Label = 'TOTAL'
      iop = 2**loper
      irc = -1
      iopt = 0
      call dRdMck(iRC,iOpt,Label,jDisp,Fock,iop)
      if (iRc /= 0) then
        write(6,*) 'IntX: Error reading MCKINT'
        write(6,'(A,A)') 'Label=',Label
        call Abend()
      end if
      call ReLoad(Fock,loper+1,nbas,norb)

      Label = 'INACTIVE'
      iop = 2**loper
      irc = -1
      iopt = 0
      call dRdMck(iRC,iOpt,Label,jDisp,Focki,iop)
      if (iRc /= 0) then
        write(6,*) 'IntX: Error reading MCKINT'
        write(6,'(A,A)') 'Label=',Label
        call Abend()
      end if
      call ReLoad(Focki,loper+1,nbas,norb)

      Label = 'MOPERT'
      iop = 2**loper
      irc = -1
      iopt = 0
      call dRdMck(iRC,iOpt,Label,jDisp,rMO,iop)
      if (iRc /= 0) then
        write(6,*) 'IntX: Error reading MCKINT'
        write(6,'(A,A)') 'Label=',Label
        call Abend()
      end if
    else

      !----------------------------------------------------------------*
      !
      ! SCF
      !
      !----------------------------------------------------------------*

      Label = 'TOTAL'
      iop = 2**loper
      irc = -1
      iopt = 0
      call dRdMck(iRC,iOpt,Label,jDisp,Focki,iop)
      if (iRc /= 0) then
        write(6,*) 'IntX: Error reading MCKINT'
        write(6,'(A,A)') 'Label=',Label
        call Abend()
      end if
      call dcopy_(ndens2,[0.0d0],0,fock,1)
      do iS=1,nSym
        js = ieor(is-1,loper)+1
        call Dyax(nOrb(is)*nIsh(js),2.0d0,Focki(ipMat(is,js)),1,Fock(ipMat(is,js)),1)
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

if (iand(nTPert(iDisp),2**1) == 2) then ! 1 el contribution
  iop = 2**loper
  if (iand(ntpert(idisp),2**4) == 2**4) then  ! from mckinley
    Label = 'ONEGRD'
    irc = -1
    iopt = 0
    call dRdMck(iRC,iOpt,Label,jDisp,Temp1,iop)
    if (iRc /= 0) then
      write(6,*) 'IntX: Error reading MCKINT'
      write(6,'(A,A)') 'Label=',Label
      call Abend()
    end if
  else                                         ! or seward
    Label = SwLbl(idisp)
    iopt = 0
    irc = -1
    call RdOne(irc,iopt,Label,jDisp,temp1,iop)
    if (iRc /= 0) then
      write(6,*) 'IntX: Error reading MCKINT'
      write(6,'(A,A)') 'Label=',Label
      call Abend()
    end if
    call DSCAL_(ndens2,1.0d0,Temp1,1)
  end if
end if

call dcopy_(nDens,[0.0d0],0,Temp2,1)
do iS=1,nSym
  do jS=1,is
    if (nBas(is)*nBas(js) /= 0) then
      if (ieor(iS-1,jS-1) == loper) then
        if (is == js) then
          call Square(Temp1(ipMatLt(iS,jS)),Temp4,1,nBas(is),nBas(is))
        else
          call dcopy_(nBas(iS)*nBas(jS),temp1(ipMatLT(iS,Js)),1,Temp4,1)
        end if
        call DGEMM_('T','N',nOrb(iS),nBas(jS),nBAs(iS),1.0d0,CMO(ipCM(iS)),nBas(is),Temp4,nBas(iS),0.0d0,Temp3,nOrb(iS))
        call DGEMM_('N','N',nOrb(is),nB(jS),nBas(jS),1.0d0,Temp3,nOrb(iS),CMO(ipCM(jS)),nBas(jS),0.0d0,Temp2(ipMat(iS,jS)),nOrb(iS))
        if (is /= js) then
          call DGEMM_('T','T',nOrb(jS),nOrb(iS),nBAs(jS),1.0d0,CMO(ipCM(jS)),nBas(js),Temp4,nBas(iS),0.0d0,Temp3,nOrb(jS))
          call DGEMM_('N','N',nOrb(js),nB(iS),nBas(iS),1.0d0,Temp3,nOrb(jS),CMO(ipCM(iS)),nBas(iS),0.0d0,Temp2(ipMat(jS,iS)), &
                      nOrb(jS))
        end if

      end if
    end if
  end do
end do

call dcopy_(ndens2,[0.0d0],0,Temp3,1)
do iS=1,nSym
  js = ieor(is-1,loper)+1
  if (nOrb(js) < 1) cycle
  do j=1,nAsh(is)+nish(is)
    do i=1,nAsh(is)+nIsh(is)
      if ((i == j) .and. (i <= nish(is)) .and. (j <= nish(is))) then
        rde = 2.0d0
      else if ((i > nish(is)) .and. (j > nish(is))) then
        rde = G1t(itri(i-nish(is)+nA(is),j-nIsh(is)+nA(is)))
      else
        rde = 0.0d0
      end if
      if (rde /= 0.0d0) call DaXpY_(nOrb(js),rDe,Temp2(ipMat(js,is)+(j-1)*nOrb(js)),1,Temp3(ipMat(js,is)+(i-1)*nOrb(js)),1)
    end do
  end do
end do
!***********************************************************************

!----------------------------------------------------------------------*
!
! One electron fock matrix done!
!
!----------------------------------------------------------------------*

if (iand(ntpert(idisp),2**2) == 4) then
  call daxpy_(nDens2,1.0d0,Temp2,1,Focki,1)
  call daxpy_(nDens2,1.0d0,Temp3,1,Fock,1)
else
  call dcopy_(ndens2,Temp2,1,FockI,1)
  call dcopy_(ndens2,Temp3,1,Fock,1)
end if

! Avoid unused argument warnings
if (.false.) call Unused_real(r)

end subroutine INTX
