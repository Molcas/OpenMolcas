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

subroutine RelEne(ErelMV,ErelDC,nSym,nBas,CMO,OCC,D,OP)
!***********************************************************************
!                                                                      *
!     Purpose: Compute relativistic correction to energy for a given   *
!              set of natural orbitals and occupation numbers          *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
dimension nBas(*), CMO(*), OCC(*), D(*), OP(*)

!----------------------------------------------------------------------*
! Compute 1-particle density matrix in AO basis                        *
!----------------------------------------------------------------------*
ij = 0
iOff1 = 0
iOff2 = 0
do iSym=1,nSym
  nBs = nBas(iSym)
  if (nBs == 0) goto 10
  do iBas=1,nBs
    do jBas=1,iBas
      ij = ij+1
      D(ij) = 0.0d0
      do iOrb=1,nBs
        iOff3 = iOff1+(iOrb-1)*nBs
        D(ij) = D(ij)+OCC(iOff2+iOrb)*CMO(iBas+iOff3)*CMO(jBas+iOff3)
      end do
      if (iBas /= jBas) D(ij) = 2.0d0*D(ij)
    end do
  end do
  iOff1 = iOff1+nBs*nBs
  iOff2 = iOff2+nBs
10 continue
end do
!----------------------------------------------------------------------*
! Compute energy contributions                                         *
!----------------------------------------------------------------------*
lOp = 0
do iSym=1,nSym
  nBs = nBas(iSym)
  lOp = lOp+(nBs**2+nBs)/2
end do
ErelMV = 0.0
iRc = -1
iOpt = 1
iComp = 1
call RdOne(iRc,iOpt,'MassVel ',iComp,OP,iSyLbl)
if (iRc == 0) then
  iRc = -1
  iOpt = 6
  iComp = 1
  call RdOne(iRc,iOpt,'MassVel ',iComp,OP,iSyLbl)
  ErelMV = DDOT_(lOP,D,1,OP,1)
end if
ErelDC = 0.0
iRc = -1
iOpt = 1
iComp = 1
call RdOne(iRc,iOpt,'Darwin  ',iComp,OP,iSyLbl)
if (iRc == 0) then
  iRc = -1
  iOpt = 6
  iComp = 1
  call RdOne(iRc,iOpt,'Darwin  ',iComp,OP,iSyLbl)
  ErelDC = DDOT_(lOP,D,1,OP,1)
end if

!----------------------------------------------------------------------*
! Exit                                                                 *
!----------------------------------------------------------------------*
return

end subroutine RelEne
