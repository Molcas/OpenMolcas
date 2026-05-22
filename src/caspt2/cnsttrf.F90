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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine CnstTrf(nTrf,Trf0,Trf)

use caspt2_global, only: TraFro
use caspt2_module, only: IfChol, NASH, NBAS, NDEL, NFRO, NISH, NRAS1, NRAS2, NRAS3, NSSH, NSYM
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTrf
real(kind=wp), intent(in) :: Trf0(nTrf)
real(kind=wp), intent(inout) :: Trf(nTrf)
integer(kind=iwp) :: I, iAsh, iIsh, IJ, ipTrfL, iSQ, iSsh, iSym, iTOrb, J, jAsh, jIsh, jSsh, nAshI, nBasI, nCor, nDelI, nFroI, &
                     nIshI, NR1, NR2, NR3, nSshI, nVir

iSQ = 0
iTOrb = 1 ! LTOrb
ipTrfL = 0
do iSym=1,nSym
  nBasI = nBas(iSym)
  nFroI = nFro(iSym)
  nIshI = nIsh(iSym)
  nAshI = nAsh(iSym)
  nSshI = nSsh(iSym)
  nDelI = nDel(iSym)
  NR1 = nRAS1(iSym)
  NR2 = nRAS2(iSym)
  NR3 = nRAS3(iSym)
  nCor = nFroI+nIshI
  nVir = nSshI+nDelI
  ipTrfL = ipTrfL+iSQ
  !! frozen + inactive
  !do iIsh=1,nFroI+nIshI
  !  Trf(ipTrfL+iIsh+nBasI*(iIsh-1)) = One
  !end do
  !! frozen
  if (IfChol) then
    do J=1,nFroI
      Trf(ipTrfL+nBasI*(J-1)+1:ipTrfL+nBasI*(J-1)+nFroI) = TraFro(nFroI*(J-1)+1:nFroI*J)
    end do
  else
    do iIsh=1,nFroI
      Trf(ipTrfL+iIsh+nBasI*(iIsh-1)) = One
    end do
  end if
  !! inactive
  do I=1,nIshI
    iIsh = nFroI+I
    do J=1,nIshI
      jIsh = nFroI+J
      IJ = I-1+nIshI*(J-1)
      Trf(ipTrfL+iIsh+nBasI*(jIsh-1)) = Trf0(iTOrb+IJ)
    end do
  end do
  iTOrb = iTOrb+nIshI*nIshI
  !! RAS1 space
  do I=1,NR1
    iAsh = nCor+I
    do J=1,NR1
      jAsh = nCor+J
      IJ = I-1+NR1*(J-1)
      Trf(ipTrfL+iAsh+nBasI*(jAsh-1)) = Trf0(iTOrb+IJ)
    end do
  end do
  iTOrb = iTOrb+NR1*NR1
  !! RAS2 space
  do I=1,NR2
    iAsh = nCor+NR1+I
    do J=1,NR2
      jAsh = nCor+NR1+J
      IJ = I-1+NR2*(J-1)
      Trf(ipTrfL+iAsh+nBasI*(jAsh-1)) = Trf0(iTOrb+IJ)
    end do
  end do
  iTOrb = iTOrb+NR2*NR2
  !! RAS3 space
  do I=1,NR3
    iAsh = nCor+NR1+NR2+I
    do J=1,NR3
      jAsh = nCor+NR1+NR2+J
      IJ = I-1+NR3*(J-1)
      Trf(ipTrfL+iAsh+nBasI*(jAsh-1)) = Trf0(iTOrb+IJ)
    end do
  end do
  iTOrb = iTOrb+NR3*NR3
  !call sqprt(trf,12)
  !! Active
  !do iAsh0=1,nAshI
  !  iAsh = nCor+iAsh0
  !  do jAsh0=1,nAshI
  !    jAsh = nCor+jAsh0
  !    Work(ipTrfL+iAsh-1+nBasI*(jAsh-1)) = Work(iTOrb+nIshI*nIshI+iAsh0-1+nAshI*(jAsh0-1))
  !    Trf(ipTrfL+iAsh+nBasI*(jAsh-1)) = Trf0(iTOrb+iAsh0-1+nAshI*(jAsh0-1))
  !  end do
  !end do
  !call sqprt(trf,12)
  !! virtual + deleted (deleted is not needed, though)
  !do iSsh=nOcc+1,nOcc+nVir
  !  Trf(ipTrfL+iSsh+nBasI*(iSsh-1)) = One
  !end do
  do I=1,nVir
    iSsh = nCor+nAshI+I
    do J=1,nVir
      jSsh = nCor+nAshI+J
      IJ = I-1+nVir*(J-1)
      Trf(ipTrfL+iSsh+nBasI*(jSsh-1)) = Trf0(iTOrb+IJ)
    end do
  end do
  iTOrb = iTOrb+nSshI*nSshI
  !call sqprt(trf,12)
  iSQ = iSQ+nBasI*nBasI

  !n123 = nAshI*nAshI !! just for CAS at present
  !iTOrb = iTOrb+n123+nSshI*nSshI
  !write(u6,*) 'transformation matrix'
  !call sqprt(trf,nbasi)
end do

return

end subroutine CnstTrf
