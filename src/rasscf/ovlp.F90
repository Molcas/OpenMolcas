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
      Subroutine Ovlp(iWay,C1,C2,Smat)

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
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero
      use rasscf_global, only: NAC
      use output_ras, only: LF
      use general_data, only: NSYM,NASH,NBAS,NFRO,NISH,NTOT1

      Implicit None

#include "warnings.h"

      Integer iWay
      Real*8 C1(*),C2(*),Smat(*)

      Character(LEN=8) Label
      Real*8, Allocatable:: OAO(:), Scr1(:), Scr2(:)
      Integer iRC, iOpt, iComp, iSyLbl, ipC, ipO, ipSMat, nAcO, iSym,   &
     &        nBs, nIs, nAs, iiOrb, ij, iOrb, jjOrb, jOrb

! prologue


      Call dCopy_(nAc*nAc,[zero],0,Smat,1)
      Call mma_allocate(OAO,nTot1,Label='OAO')

! read the overlap integrals

      iRc=-1
      iOpt=ibset(ibset(0,sNoOri),sNoNuc)
      iComp=1
      iSyLbl=1
      Label='Mltpl  0'
      Call RdOne(iRc,iOpt,Label,iComp,OAO,iSyLbl)
      If ( iRc.ne.0 ) Then
         Write(LF,*)
         Write(LF,*) ' *** Error in subroutine Ovlp ***'
         Write(LF,*) ' premature abort in subroutine RdOne'
         Write(LF,*) ' reading label: ',Label
         Write(LF,*)' RASSCF is trying to orthonormalize orbitals but'
         Write(LF,*)' could not read overlaps from ONEINT. Something'
         Write(LF,*)' is wrong with the file, or possibly with the'
         Write(LF,*)' program. Please check.'
         Write(LF,*)
         Call Quit(_RC_IO_ERROR_READ_)
      End If

! compute the S-matrix for all orbitals and select elements

      ipC=1
      ipO=1
      ipSmat=1
      nAcO=0
      Do iSym = 1,nSym
        nBs = nBas(iSym)
        nIs = nFro(iSym)+nIsh(iSym)
        nAs = nAsh(iSym)
        If ( nBs.gt.0 ) then
          Call mma_allocate(Scr1,nBs*nBs,Label='Scr1')
          Call mma_allocate(Scr2,nBs*nBs,Label='Scr2')
          Call Square(OAO(ipO),Scr1,1,nBs,nBs)
          Call DGEMM_('N','N',                                          &
     &                nBs,nBs,nBs,                                      &
     &                1.0d0,Scr1,nBs,                                   &
     &                C1(ipC),nBs,                                      &
     &                0.0d0,Scr2,nBs)
          Call DGEMM_('T','N',                                          &
     &                nBs,nBs,nBs,                                      &
     &                1.0d0,C2(ipC),nBs,                                &
     &                Scr2,nBs,                                         &
     &                0.0d0,Scr1,nBs)
          If ( iWay.eq.0 ) then
            ij = 1
            Do iOrb = 1,nBs
              Do jOrb = 1,nBs
                If ( jOrb.le.iOrb ) then
                  Smat(ipSmat) = Scr1(ij)
                  ipSmat = ipSmat+1
                End If
                ij = ij+1
              End Do
            End Do
          Else
            ij = 1
            Do iOrb = 1,nBs
              Do jOrb = 1,nBs
!               If ( jOrb.le.iOrb ) then
                  If ( (iOrb.gt.nIs) .and. (iOrb.le.(nIs+nAs)) ) then
                    If ( (jOrb.gt.nIs) .and. (jOrb.le.(nIs+nAs)) ) then
                      iiOrb = iOrb-nIs+nAcO
                      jjOrb = jOrb-nIs+nAcO
                      ipSmat = jjOrb+(iiOrb*iiOrb-iiOrb)/2
                      ipSmat = jjOrb+(iiOrb-1)*nAc
                      Smat(ipSmat) = Scr1(ij)
                    End If
                  End If
!               End If
                ij = ij+1
              End Do
            End Do
          End If
          Call mma_deallocate(Scr2)
          Call mma_deallocate(Scr1)
        End If
        ipC=ipC+nBs*nBs
        ipO=ipO+(nBs*nBs+nBs)/2
        nAcO = nAcO+nAs
      End Do

! epilogue

      Call mma_deallocate(OAO)

      End Subroutine Ovlp
