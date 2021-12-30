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
! Copyright (C) 2005,2006, Thomas Bondo Pedersen                       *
!***********************************************************************

subroutine GetRawPAOs(R,C,nBas,nOrb,nFro,nOrb2Loc,nSym,Normalize)
! Thomas Bondo Pedersen, December 2005.
! - revised january 2006 (Thomas Bondo Pedersen).
!
! Purpose: compute projected AOs spanning the primary space from the
!          formula
!
!          R = 1 - Do*S = D*S
!
!          where Do is the "AO density" matrix of the orthogonal
!          complement orbitals and D that of the orbital space to be
!          localised (primary space).
!          S is the AO overlap matrix (which is read from disk).
!
!          Which formula is used depends on the dimensions of the
!          two spaces (such that the most economical is computed).
!
!          If (Normalize): normalize each PAO (recommended).
!                          Note that if the norm of the PAO is
!                          smaller than 1.0d-6, it will not be
!                          normalized (it is left unchanged).
!
!-------------------------------------------------------------
!-TODO/FIXME: it is in most cases faster to use Do*S=C*(C^T*S)
!-------------------------------------------------------------

implicit real*8(a-h,o-z)
integer nBas(nSym), nOrb(nSym), nFro(nSym), nOrb2Loc(nSym)
real*8 R(*), C(*)
logical Normalize
#include "WrkSpc.fh"
character*10 SecNam
parameter(SecNam='GetRawPAOs')
character*80 Txt
external ddot_

! Read the overlap matrix from disk.
! ----------------------------------

lOvlp = nBas(1)**2
do iSym=2,nSym
  lOvlp = lOvlp+nBas(iSym)**2
end do
call GetMem('Ovlp','Allo','Real',ipOvlp,lOvlp)
call GetOvlp_Localisation(Work(ipOvlp),'Sqr',nBas,nSym)

! Compute R.
! ----------

lDo = nBas(1)**2
do iSym=2,nSym
  lDo = max(lDo,nBas(iSym)**2)
end do
call GetMem('Do','Allo','Real',ipDo,lDo)

kOff = 1
kOffS = ipOvlp
do iSym=1,nSym

  nB = nBas(iSym)
  if (nB > 0) then

    nF = nFro(iSym)
    nO2L = nOrb2Loc(iSym)
    nRest = nOrb(iSym)-nF-nO2L
    nOrth = nF+nRest ! dim. of orthogonal complement

    if (nO2L < 1) then ! R = 0
      call fZero(R(kOff),nB**2)
    else if (nOrth < 0) then ! error
      call SysAbendMsg(SecNam,'Dimension of orthogonal complement space < 0',' ')
    else if (nOrth == 0) then ! R = 1
      call fZero(R(kOff),nB**2)
      do i=1,nB
        R(kOff-1+nB*(i-1)+i) = 1.0d0
      end do
    else if (nOrth < nO2L) then ! R = 1 - Do*S
      if (nRest > 0) then
        lOff = kOff+nB*(nF+nO2L)
        call GetDens_Localisation(Work(ipDo),C(lOff),nB,nRest)
      else
        call fZero(Work(ipDo),nB**2)
      end if
      if (nF > 0) then
        call GetDens_Localisation(R(kOff),C(kOff),nB,nF)
        call dAXPY_(nB**2,1.0d0,R(kOff),1,Work(ipDo),1)
      end if
      call DGEMM_('N','N',nB,nB,nB,-1.0d0,Work(ipDo),nB,Work(kOffS),nB,0.0d0,R(kOff),nB)
      do i=1,nB
        R(kOff-1+nB*(i-1)+i) = R(kOff-1+nB*(i-1)+i)+1.0d0
      end do
    else ! R = D*S
      lOff = kOff+nB*nF
      call GetDens_Localisation(Work(ipDo),C(lOff),nB,nO2L)
      call DGEMM_('N','N',nB,nB,nB,1.0d0,Work(ipDo),nB,Work(kOffS),nB,0.0d0,R(kOff),nB)
    end if

    kOff = kOff+nB**2
    kOffS = kOffS+nB**2

  end if

end do

! If requested, normalize the PAOs.
! ---------------------------------

if (Normalize) then
  kOff = 1
  kOffS = ipOvlp
  do iSym=1,nSym
    nB = nBas(iSym)
    if (nB > 0) then
      call DGEMM_('N','N',nB,nB,nB,1.0d0,Work(kOffS),nB,R(kOff),nB,0.0d0,Work(ipDo),nB)
      do mu=0,nB-1
        kR = kOff+nB*mu
        kSR = ipDo+nB*mu
        Ovlp = dDot_(nB,R(kR),1,Work(kSR),1)
        if (Ovlp > 1.0d-6) then
          Fac = 1.0d0/sqrt(Ovlp)
          call dScal_(nB,Fac,R(kR),1)
        else if (Ovlp < 0.0d0) then
          write(Txt,'(A,1P,D15.5)') 'Overlap = ',Ovlp
          call SysAbendMsg(SecNam,'Negative raw PAO overlap!',Txt)
        end if
      end do
      kOff = kOff+nB**2
      kOffS = kOffS+nB**2
    end if
  end do
end if

call GetMem('Do','Free','Real',ipDo,lDo)
call GetMem('Ovlp','Free','Real',ipOvlp,lOvlp)

end subroutine GetRawPAOs
