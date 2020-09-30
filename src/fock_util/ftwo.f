************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996, Markus P. Fuelscher                              *
*               1998, Roland Lindh                                     *
************************************************************************
      Subroutine Ftwo(case,ExFac,
     &                iSym,kSym,
     &                iBas,kBas,
     &                off_sqMat,off_ltMat,
     &                D1I,FI,D1A,FA,PQRS)
************************************************************************
*                                                                      *
*     Assemble Fock matrices FI and FA (in AO-basis)                   *
*                                                                      *
*     calling arguments:                                               *
*     case    : input, integer                                         *
*               symmetry case number; (II!II)=1, (II!KK)=2, (IK!IK)=3  *
*     ExFac   : input, real*8                                          *
*               coefficient of "exact exchange"                        *
*     iSym    : input, integer                                         *
*               symmetry species iSym                                  *
*     kSym    : input, integer                                         *
*               symmetry species kSym                                  *
*     iBas    : input, integer                                         *
*               1-st, fixed index of ERIs                              *
*     kBas    : input, integer                                         *
*               2-nd, fixed index of ERIs                              *
*     off_sqMat : input, array of integer                              *
*               offset of one-electron integrals (squared format)      *
*     off_ltMat : input, array of integer                              *
*               offset of one-electron integrals (lower triangular)    *
*     D1I     : input, array of real*8                                 *
*               one-body density matrix (frozen+inactive, AO-basis)    *
*     FI      : output, array of real*8                                *
*               Fock matrix (frozen+inactive, AO-basis)                *
*     D1A     : input, array of real*8                                 *
*               one-body density matrix (active, AO-basis)             *
*     FA      : output, array of real*8                                *
*               Fock matrix (active, AO-basis)                         *
*     PQRS    : input, array of real*8                                 *
*               two-electron integrals (AO-basis)                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
*     Modified to only need the unique symmetry blocks, R. Lindh,      *
*     March '98.                                                       *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Parameter ( Zero=0.0D0, Half=-0.5d0 , One=1.0D0, Two=2.0d0)
#include "rasdim.fh"
#include "general.fh"
*
      Integer   case, off_sqMat, off_ltMat
      Dimension off_sqMat(*), off_ltMat(*)
      Dimension D1I(*), FI(*), D1A(*), FA(*), PQRS(*)
*
      iTri(i)=(i*i-i)/2
      iSqr(i)=(i-1)*nBas(iSym)
      kSqr(i)=(i-1)*nBas(kSym)
*
*
*                                                                      *
************************************************************************
*                                                                      *
*     symmetry case (II|II)
*     Both Coulomb and Exchange contributions
*                                                                      *
************************************************************************
*                                                                      *
      If ( case.eq.1 ) Then
*
*------- Coulombic contribution
*
         iOff=off_ltMat(iSym)+iTri(iBas)+kBas
         kOff=off_sqMat(kSym)+1
         FI(iOff) = FI(iOff)
     &            + dDot_(nBas(kSym)*nBas(kSym),D1I(kOff),1,PQRS,1)
         FA(iOff) = FA(iOff)
     &            + dDot_(nBas(kSym)*nBas(kSym),D1A(kOff),1,PQRS,1)
*
*------- Exchange contribution
*
         If (ExFac.ne.Zero) Then
            iOff=off_sqMat(iSym)+iSqr(iBas)+1
            kOff=off_ltMat(kSym)+iTri(kBas)+1
*            Call DGeMX (kBas,nBas(iSym),Half*ExFac,
*     &                  PQRS,nBas(iSym),
*     &                  D1I(iOff),1,
*     &                  FI(kOff),1)
*            Call DGeMX (kBas,nBas(iSym),Half*ExFac,
*     &                  PQRS,nBas(iSym),
*     &                  D1A(iOff),1,
*     &                  FA(kOff),1)
            CALL DGEMV_('N',kBas,nBas(iSym),(Half*ExFac),
     &                  PQRS,nBas(iSym),
     &                  D1I(iOff),1,1.0D0,FI(kOff),1)
            CALL DGEMV_('N',kBas,nBas(iSym),(Half*ExFac),
     &                  PQRS,nBas(iSym),
     &                  D1A(iOff),1,1.0D0,FA(kOff),1)
            If ( iBas.ne.kBas ) then
               iOff=off_ltMat(iSym)+iTri(iBas)+1
               kOff=off_sqMat(kSym)+kSqr(kBas)+1
*               Call DGeMX (iBas,nBas(iSym),Half*ExFac,
*     &                     PQRS,nBas(iSym),
*     &                     D1I(kOff),1,
*     &                     FI(iOff),1)
*               Call DGeMX (iBas,nBas(iSym),Half*ExFac,
*     &                     PQRS,nBas(iSym),
*     &                     D1A(kOff),1,
*     &                     FA(iOff),1)
               CALL DGEMV_('N',iBas,nBas(iSym),(Half*ExFac),
     &                     PQRS,nBas(iSym),
     &                     D1I(kOff),1,1.0D0,FI(iOff),1)
               CALL DGEMV_('N',iBas,nBas(iSym),(Half*ExFac),
     &                     PQRS,nBas(iSym),
     &                     D1A(kOff),1,1.0D0,FA(iOff),1)
            End If
         End If
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     symmetry case (II|KK) and (KK|II)
*     Only Coulomb contributions
*     This is limited to iBas >= jBas
*     Code modified to only require uniqe symmetry blocks of
*     two-electron integrals. R. Lindh, March '98.
*                                                                      *
************************************************************************
*                                                                      *
      If ( case.eq.2 .and. iSym.gt.kSym) then
         jBas=kBas    ! To avoid confusion.
         iOff_ij=off_ltMat(iSym)+iTri(iBas)+jBas
         kOff=off_sqMat(kSym)+1
*
         FI(iOff_ij) = FI(iOff_ij)
     &               + DDot_(nBas(kSym)**2,D1I(kOff),1,PQRS,1)
         FA(iOff_ij) = FA(iOff_ij)
     &               + DDot_(nBas(kSym)**2,D1A(kOff),1,PQRS,1)
*
         D1I_ij=D1I(off_sqMat(iSym)+(jBas-1)*nBas(iSym)+iBas)
         D1A_ij=D1A(off_sqMat(iSym)+(jBas-1)*nBas(iSym)+iBas)
         If (iBas.ne.jBas) Then
            D1I_ij=D1I_ij*Two
            D1A_ij=D1A_ij*Two
         End If
         Do k = 1, nBas(kSym)
            Do l = 1, k
               kl = (l-1)*nBas(kSym)+k
               iOff_kl = off_ltMat(kSym)+iTri(k)+l
               FI(iOff_kl)=FI(iOff_kl)+D1I_ij*PQRS(kl)
               FA(iOff_kl)=FA(iOff_kl)+D1A_ij*PQRS(kl)
            End Do
         End Do
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
*     symmetry case (IK!IK)
*     Only Exchange contributions
*                                                                      *
************************************************************************
*                                                                      *
      If ( case.eq.3 .and. ExFac.ne.Zero) then
         iOff=off_sqMat(iSym)+iSqr(iBas)+1
         kOff=off_ltMat(kSym)+iTri(kBas)+1
*         Call DGeMX  (kBas,nBas(iSym),Half*ExFac,
*     &                PQRS,nBas(kSym),
*     &                D1I(iOff),1,
*     &                FI(kOff),1)
*         Call DGeMX  (kBas,nBas(iSym),Half*ExFac,
*     &                PQRS,nBas(kSym),
*     &                D1A(iOff),1,
*     &                FA(kOff),1)
         CALL DGEMV_('N',kBas,nBas(iSym),(Half*ExFac),
     &                PQRS,nBas(kSym),
     &                D1I(iOff),1,1.0D0,FI(kOff),1)
         CALL DGEMV_('N',kBas,nBas(iSym),(Half*ExFac),
     &                PQRS,nBas(kSym),
     &                D1A(iOff),1,1.0D0,FA(kOff),1)
         iOff=off_ltMat(iSym)+iTri(iBas)+1
         kOff=off_sqMat(kSym)+kSqr(kBas)+1
*         Call DGeMTX (nBas(kSym),iBas,Half*ExFac,
*     &                PQRS,nBas(kSym),
*     &                D1I(kOff),1,
*     &                FI(iOff),1)
*         Call DGeMTX (nBas(kSym),iBas,Half*ExFac,
*     &                PQRS,nBas(kSym),
*     &                D1A(kOff),1,
*     &                FA(iOff),1)
         CALL DGEMV_('T',nBas(kSym),iBas,(Half*ExFac),
     &                PQRS,nBas(kSym),
     &                D1I(kOff),1,1.0D0,FI(iOff),1)
         CALL DGEMV_('T',nBas(kSym),iBas,(Half*ExFac),
     &                PQRS,nBas(kSym),
     &                D1A(kOff),1,1.0D0,FA(iOff),1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
*
      Return
      End
