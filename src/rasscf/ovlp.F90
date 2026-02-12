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
! Copyright (C) 1999, Markus P. Fuelscher                              *
!***********************************************************************

subroutine Ovlp(iWay,C1,C2,Smat)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Compute the orbital overlap matrix between the MO sets C1 and C2 *
!                                                                      *
!     calling arguments:                                               *
!     iWay    : integer                                                *
!               =0 : S-matrix for all orbitals                         *
!                    symmetry blocked, triangular                      *
!               =1 : S-matrix for active orbitals only                 *
!                    no symmetry, triangular                           *
!     C1      : real*8                                                 *
!               MO-basis                                               *
!     C2      : real*8                                                 *
!               MO-basis                                               *
!     S       : real*8                                                 *
!               orbital overlap matrix                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1999                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

use OneDat, only: sNoNuc, sNoOri
use rasscf_global, only: NAC
use output_ras, only: LF
use general_data, only: NSYM, NASH, NBAS, NFRO, NISH, NTOT1
use Constants, only: Zero
use stdalloc, only: mma_allocate, mma_deallocate

implicit none
integer iWay
real*8 C1(*), C2(*), Smat(*)
character(len=8) Label
real*8, allocatable :: OAO(:), Scr1(:), Scr2(:)
integer iRC, iOpt, iComp, iSyLbl, ipC, ipO, ipSMat, nAcO, iSym, nBs, nIs, nAs, iiOrb, ij, iOrb, jjOrb, jOrb
#include "warnings.h"

! prologue

call dCopy_(nAc*nAc,[zero],0,Smat,1)
call mma_allocate(OAO,nTot1,Label='OAO')

! read the overlap integrals

iRc = -1
iOpt = ibset(ibset(0,sNoOri),sNoNuc)
iComp = 1
iSyLbl = 1
Label = 'Mltpl  0'
call RdOne(iRc,iOpt,Label,iComp,OAO,iSyLbl)
if (iRc /= 0) then
  write(LF,*)
  write(LF,*) ' *** Error in subroutine Ovlp ***'
  write(LF,*) ' premature abort in subroutine RdOne'
  write(LF,*) ' reading label: ',Label
  write(LF,*) ' RASSCF is trying to orthonormalize orbitals but'
  write(LF,*) ' could not read overlaps from ONEINT. Something'
  write(LF,*) ' is wrong with the file, or possibly with the'
  write(LF,*) ' program. Please check.'
  write(LF,*)
  call Quit(_RC_IO_ERROR_READ_)
end if

! compute the S-matrix for all orbitals and select elements

ipC = 1
ipO = 1
ipSmat = 1
nAcO = 0
do iSym=1,nSym
  nBs = nBas(iSym)
  nIs = nFro(iSym)+nIsh(iSym)
  nAs = nAsh(iSym)
  if (nBs > 0) then
    call mma_allocate(Scr1,nBs*nBs,Label='Scr1')
    call mma_allocate(Scr2,nBs*nBs,Label='Scr2')
    call Square(OAO(ipO),Scr1,1,nBs,nBs)
    call DGEMM_('N','N',nBs,nBs,nBs,1.0d0,Scr1,nBs,C1(ipC),nBs,0.0d0,Scr2,nBs)
    call DGEMM_('T','N',nBs,nBs,nBs,1.0d0,C2(ipC),nBs,Scr2,nBs,0.0d0,Scr1,nBs)
    if (iWay == 0) then
      ij = 1
      do iOrb=1,nBs
        do jOrb=1,nBs
          if (jOrb <= iOrb) then
            Smat(ipSmat) = Scr1(ij)
            ipSmat = ipSmat+1
          end if
          ij = ij+1
        end do
      end do
    else
      ij = 1
      do iOrb=1,nBs
        do jOrb=1,nBs
          !if (jOrb <= iOrb) then
          if ((iOrb > nIs) .and. (iOrb <= (nIs+nAs))) then
            if ((jOrb > nIs) .and. (jOrb <= (nIs+nAs))) then
              iiOrb = iOrb-nIs+nAcO
              jjOrb = jOrb-nIs+nAcO
              ipSmat = jjOrb+(iiOrb*iiOrb-iiOrb)/2
              ipSmat = jjOrb+(iiOrb-1)*nAc
              Smat(ipSmat) = Scr1(ij)
            end if
          end if
          !end if
          ij = ij+1
        end do
      end do
    end if
    call mma_deallocate(Scr2)
    call mma_deallocate(Scr1)
  end if
  ipC = ipC+nBs*nBs
  ipO = ipO+(nBs*nBs+nBs)/2
  nAcO = nAcO+nAs
end do

! epilogue

call mma_deallocate(OAO)

end subroutine Ovlp
