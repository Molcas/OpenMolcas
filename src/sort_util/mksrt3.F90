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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!               1993, Per-Olof Widmark                                 *
!***********************************************************************

subroutine MkSrt3(iRc,iSquar,nIrrep,nBas,nSkip)
!***********************************************************************
!                                                                      *
!    Purpose: Prepare the address table for the virtual disk and       *
!             initialize the memory                                    *
!                                                                      *
!    Global data declarations (Include files) :                        *
!    TwoDat : table of contents and auxiliary information              *
!    TowRc  : Table of return code                                     *
!                                                                      *
!    calling arguments:                                                *
!    iRc    : return code                                              *
!             A value of 0 (zero) is returned upon successful          *
!             completion of the request. A nonzero value indi-         *
!             cates an error.                                          *
!    iSquar  : ordering flag                                           *
!    nIrrep  : number of irreducible representations                   *
!    nBas    : number of basis functions per irred. rep.               *
!    nSkip   : flag to exlude symmetry combinations                    *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M. P. Fuelscher and P.O. Widmark                                 *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

implicit integer(A-Z)

dimension nBas(nIrrep), nSkip(nIrrep)

#include "srt1.fh"
#include "Molcas.fh"
#include "TwoDat.fh"
#include "print.fh"
integer nSln_Temp(mSyBlk), lSll_Temp(mSyBlk)
logical Square

iRout = 80
iPrint = nPrint(iRout)
if (iPrint > 5) write(6,*) ' >>> Enter MKSRT3 <<<'

Square = .true.
if (iSquar == 0) Square = .false.

call iCopy(mSyBlk,nSln,1,nSln_Temp,1)
call iCopy(mSyBlk,lSll,1,lSll_Temp,1)
do iSyBlk=1,mSyBlk
  nSln(iSyBlk) = 0
  lSll(iSyBlk) = 0
end do

! loop over all symmetry blocks
!
!write(6,*)
!write(6,'(2X,A,I3.3,A)') '*** (I)-level message MKSRT3 ***'
!write(6,'(2X,A,I8    )') 'RAMD_size   =',RAMD_size
!write(6,'(2X,A,I8    )') 'RAMD_anchor =',RAMD_anchor
!write(6,*)
!write(6,'(2X,84A1)') ('*',i=1,84)

iRc = 0
iOff = RAMD_anchor
nInts = 0
nSym = TocTwo(isSym)
mxSyP = nSym*(nSym+1)/2
do iSymi=1,nSym
  ib = nBas(iSymi)
  iSkip = nSkip(iSymi)
  do jSymj=1,iSymi
    iSymj = 1+ieor(iSymi-1,jSymj-1)
    iSyblj = jSymj+iSymi*(iSymi-1)/2
    jb = nBas(jSymj)
    ibj = ib*jb
    if (iSymi == jSymj) ibj = ib*(ib+1)/2
    jSkip = nSkip(jSymj)
    kSymMx = iSymi
    if (Square) kSymMx = nSym
    do kSymk=1,kSymMx
      kb = nBas(kSymk)
      kSkip = nSkip(kSymk)
      lSymMx = jSymj
      if ((kSymk /= iSymi) .or. (Square)) lSymMx = kSymk
      do lSyml=1,lSymMx
        kSyml = 1+ieor(kSymk-1,lSyml-1)
        kSybll = lSyml+kSymk*(kSymk-1)/2
        if (ieor(iSymj-1,kSyml-1) == 0) then
          lb = nBas(lSyml)
          kbl = kb*lb
          if (kSymk == lSyml) kbl = kb*(kb+1)/2
          lSkip = nSkip(lSyml)
          iSyBlk = kSybll+mxSyP*(iSyblj-1)
          if ((iSkip+jSkip+kSkip+lSkip == 0) .and. (ibj*kbl /= 0)) then
            !write(6,*) 'iSymi,jSymj,kSymk,lSyml=',iSymi,jSymj,kSymk,lSyml
            !write(6,*) 'iSyBlk=',iSyBlk
            !write(6,*) 'nBatch(iSyBlk)=',nBatch(iSyBlk)
            ! check if there is enough space available

            lBuf = ibj*kbl
            nInts = nInts+lBuf
            if (nInts >= RAMD_size) then
              iRc = 001
              write(6,*)
              write(6,'(2X,A,I3.3,A)') '*** (W)-level message MKSRT3',iRc,' ***'
              write(6,'(2X,A)') 'There is not enough space on the RAM disk'
              write(6,'(2X,A)') 'The program will resume normal activity'
              write(6,*)
              call iCopy(mSyBlk,nSln_Temp,1,nSln,1)
              call iCopy(mSyBlk,lSll_Temp,1,lSll,1)
              return
            end if

            ! save start address of this symmetry block

            iBatch = nBatch(iSyBlk)
            RAMD_adr(iBatch) = iOff
            nSln(iSyBlk) = 1
            lSll(iSyBlk) = lBuf
            !write(6,'(2X,A,4I2,A,2I8,A,2I8)') ' iSym,jSym,kSym,lSym',iSymi,jSymj,kSymk,lSyml,' lBuf,nInts',lBuf,nInts, &
            !                                  ' iBatch,iOff',iBatch,iOff

            ! Init integrals

            call dCopy_(lBuf,[0.0d0],0,RAMD_ints(iOff),1)

            ! update pointers

            iOff = iOff+lBuf

          end if
        end if
      end do
    end do
  end do
end do
!write (6,'(2X,84A1)') ('*',i=1,84)

!----------------------------------------------------------------------*
!     compute offsets                                                  *
!----------------------------------------------------------------------*
!
! Note that we fake the information so that the iSyBlk index is
! really the iSyBlk index as it will be needed in Sort1c!

do iSyBlk=1,mSyBlk
  iStBin(iSyBlk) = iSyBlk
end do

! define initial load point

!RAMD_next = RAMD_adr(1)

if (iPrint > 5) write(6,*) ' >>> Exit MKSRT3 <<<'

return

end subroutine MkSrt3
