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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               1996, Martin Schuetz                                   *
!               2017, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine TraClc_i(iterLw,nD)
!***********************************************************************
!                                                                      *
! purpose: compute traces                                              *
!                                                                      *
! input:                                                               *
!   Dens    : a few last density matrix differences (nDT,NumDT)        *
!   TwoHam  : a few last two-electron hamiltonians (nDT,NumDT)         *
!   OneHam  : one-electron hamiltonian of length nDT                   *
!   iterLw  : lowest iteration count ...                               *
!             traces are computed from iterLw...iter                   *
!                                                                      *
! output:                                                              *
!   TrDh    : Traces of dD(i)*h of size (nTr)                          *
!   TrDP    : Traces of dD(i)*dP(j) of size (nTr,nTr)                  *
!   TrDD    : Traces of dD(i)*dD(j) of size (nTr,nTr)                  *
!                                                                      *
! called from: OptClc                                                  *
!                                                                      *
! calls to: RWDTG                                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Written by:                                                          *
! P.O. Widmark, M.P. Fuelscher and P. Borowski                         *
! University of Lund, Sweden, 1992                                     *
! modified by M. Schuetz, 1996                                         *
! - traces are recomputed for all iterations between iterLw & iter     *
!                                                                      *
!***********************************************************************
      use InfSCF, only: iDKeep, Iter, nBT, MapDns, iDisk
      use stdalloc, only: mma_allocate, mma_deallocate
      use SCF_Arrays, only: OneHam, TwoHam, Vxc, Dens, TrDh, TrDP, TrDD
      use Constants, only: Zero
      Implicit None
      Integer nD, IterLw
!---- Define local variables
      Integer ii, iPosL, iD, i, iPos
      Real*8, External:: DDot_
      Real*8, Dimension(:,:), Allocatable, Target:: Aux1, Aux2, Aux3
      Real*8, Dimension(:,:), Pointer:: pDens, pTwoHam, pVxc
#ifdef _DEBUGPRINT_
      Integer iR, iC
#endif
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
      If (iDKeep.lt.0) Return
!
!----------------------------------------------------------------------*
!                                                                      *
! Expand the DFT contribution around the reference density             *
! i.e. the last density.                                               *
!                                                                      *
! Note that for ii larger than 1 that the stored density is the        *
! density difference between the two last iterations.                  *
!                                                                      *
!----------------------------------------------------------------------*
!
!
!     Loop over densities to interpolate over
!
      Do ii = iterLw, iter
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!        Get the one-particle density_ii
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!        Get pointer to the density
!
         iPosL=MapDns(ii)
!
         If(iPosL.le.0) Then
!
!           If not in memory pick up from disk
!
            If (.Not.Allocated(Aux1)) Call mma_allocate(Aux1,nBT,nD,Label='Aux1')
!
!           Pick up the density matrix and external potential
!
            Call RWDTG(-iPosL,Aux1,nBT*nD,'R','DENS  ',iDisk,SIZE(iDisk,1))
            pDens => Aux1
         Else
            pDens => Dens(1:nBT,1:nD,iPosL)
         End If
!
         Do iD = 1, nD
!
!           Trace the one-electron density with the one-electron
!           Hamiltonian.
!
            TrDh(ii,ii,iD)=DDot_(nBT,pDens(:,iD),1,OneHam,1)
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
         End Do ! iD
!
         Nullify(pDens)
!
      End Do ! ii
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!
#ifdef _DEBUGPRINT_
      Do iD = 1, nD
         Write(6,'(a)') 'traclc: TrDh'
         Write(6,'(6f16.8)') (TrDh(ii,ii,iD),ii=1,iter)
      End Do
#endif
      If (Allocated(Aux1)) Call mma_deallocate(Aux1)
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
      Do ii = iterLw, iter
!
         Do iD = 1, nD
!
!           Trace the two-electron contribution of the Fock matrix with
!           the density. Diagonal terms.
!
            iPosL=MapDns(ii)
            If (iPosL.gt.0) Then
               TrDP(ii,ii,iD)=DDot_(nBT,Dens(1,iD,iPosL),1,TwoHam(1,iD,iPosL),1)    &
                             +DDot_(nBT,Dens(1,iD,iPosL),1,Vxc   (1,iD,iPosL),1)
               TrDD(ii,ii,iD)=DDot_(nBT,Dens(1,iD,iPosL),1,Dens(1,iD,iPosL),1)
            Else
               Write(6,'(a)') 'traclc: should not happen!!!'
               TrDP(ii,ii,iD)=Zero
               Call Abend()
            End If
         End Do ! iD
!
         Do i = 1, ii - 1
!
!           Get pointer to density
!
            iPos=MapDns(i)
            If(iPos.le.0) Then
               If (.NOT.Allocated(Aux1)) Then
                  Call mma_allocate(Aux1,nBT,nD,Label='Aux1')
                  Call mma_allocate(Aux2,nBT,nD,Label='Aux2')
                  Call mma_allocate(Aux3,nBT,nD,Label='Aux3')
               End If
               Call RWDTG(-iPos,Aux1,nBT*nD,'R','TWOHAM',iDisk,SIZE(iDisk,1))
               Call RWDTG(-iPos,Aux2,nBT*nD,'R','dVxcdR',iDisk,SIZE(iDisk,1))
               Call RWDTG(-iPos,Aux3,nBT*nD,'R','DENS  ',iDisk,SIZE(iDisk,1))
               pTwoHam => Aux1
               pVxc    => Aux2
               pDens   => Aux3
            Else
               pTwoHam => TwoHam(1:nBT,1:nD,iPos)
               pVxc    => Vxc   (1:nBT,1:nD,iPos)
               pDens   => Dens  (1:nBT,1:nD,iPos)
            End If
!
            Do iD = 1, nD
               TrDP(i,ii,iD) = DDot_(nBT,Dens(1,iD,iPosL),1,pTwoHam(:,iD),1)
               TrDP(ii,i,iD) = TrDP(i,ii,iD)
               TrDP(i,ii,iD) = TrDP(i,ii,iD)+ DDot_(nBT,Dens(1,iD,iPosL),1,pVxc   (:,iD),1)
               TrDP(ii,i,iD) = TrDP(ii,i,iD)+ DDot_(nBT,Vxc (1,iD,iPosL),1,pDens(:,iD),1)
               TrDD(i,ii,iD) = DDot_(nBT,Dens(1,iD,iPosl),1,pDens(:,iD),1)
               TrDD(ii,i,iD) = TrDD(i,ii,iD)
!
            End Do ! iD
!
#ifdef _DEBUGPRINT_
            Do iD = 1, nD
               Write(6,*) 'iteration:',ii
               Write(6,'(a)') 'traclc: TrDh'
               Do iR = 1, ii
                  Write (6,'(6f16.8)')TrDh(iR,iR,iD)
               End Do
               Write(6,'(a)') 'traclc: TrDP'
               Do iR = 1, ii
                  Write (6,'(6f16.8)')(TrDP(iR,iC,iD),iC=1,ii)
               End Do
               Write(6,'(a)') 'traclc: TrDD'
               Do iR = 1, ii
                  Write (6,'(6f16.8)')(TrDD(iR,iC,iD),iC=1,ii)
               End Do
            End Do ! iD
#endif
            Nullify(pTwoHam)
            Nullify(pVxc   )
            Nullify(pDens  )
!
         End Do ! i
!
      End Do ! ii
!
      If (Allocated(Aux1)) Then
         Call mma_deallocate(Aux1)
         Call mma_deallocate(Aux2)
         Call mma_deallocate(Aux3)
      End If
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
      Return
      End SubRoutine TraClc_i
