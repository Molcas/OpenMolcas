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
* Copyright (C) 1995, Martin Schuetz                                   *
************************************************************************
      SubRoutine ChkOrt(CMO,nCMO,Ovl,nOvl,OffMx)
************************************************************************
*                                                                      *
*     purpose: Check orthogonality of CMOs                             *
*                                                                      *
*     input:                                                           *
*       CMO     : orthonormal (?) vectors                              *
*       Ovl     : overlap matrix in AO basis                           *
*                                                                      *
*     output:                                                          *
*       OffMx   : maximal off diagonal element                         *
*                                                                      *
*     called from: NewOrb                                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. Schuetz                                                       *
*     University of Lund, Sweden, 1995                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
*
      Real*8 CMO(nCMO),Ovl(nOvl)
*
      Parameter (OrtThr = 1.0d-9)
*
*     declaration of local vars
      Integer iOffMx,jOffMx,iDgNo1
      Real*8 OffMx,DgNo1
      Logical termin
      Real*8, Dimension(:), Allocatable:: OvlS, Aux
*
      Call mma_allocate(OvlS,MaxBas**2,Label='OvlS')
      Call mma_allocate(Aux,MaxBas**2,Label='Aux')
*
      ij   = 1
      iCMO = 1
      termin = .FALSE.
      Do iSym = 1, nSym
        OffMx = Zero
        DgNo1 = Zero
        nBs = nBas(iSym)
        nOr = nOrb(iSym)
        If (nOr.gt.0) Then
           Call Square(Ovl(ij),OvlS,1,nBs,nBs)
           Call DGEMM_('N','N',
     &                 nBs,nOr,nBs,
     &                 1.0d0,OvlS,nBs,
     &                       CMO(iCMO),nBs,
     &                 0.0d0,Aux,nBs)
           Call DGEMM_('T','N',
     &                 nOr,nOr,nBs,
     &                 1.0d0,CMO(iCMO),nBs,
     &                       Aux,nBs,
     &                 0.0d0,OvlS,nOr)
*          get largest non zero off diag element
           Do i=1,nOr
             Do j = 1, i-1
               iOff = (j-1)*nOr+i
               OffMx = Max(OffMx,Abs(OvlS(iOff)))
             End Do
           End Do
*          get diag element most different from one
           Do i=1,nOr
             iOff = (i-1)*nOr+i
             DgNo1 = Max(DgNo1,Abs(OvlS(iOff)-One))
           End Do
*          check, if orthogonality violated
           If ((OffMx.gt.OrtThr).OR.(DgNo1.gt.OrtThr)) Then
*            Ooooops...
C            Write(6,*) 'WARNING: reorthonormalizing MOs...',OffMx,DgNo1
*
*            try to re-orthonormalize...
             Call Orthox(OvlS,CMO(iCMO),nOr,nBs)
*
*            and test again...
             OffMx = Zero
             DgNo1 = Zero
             Call Square(Ovl(ij),OvlS,1,nBs,nBs)
             Call DGEMM_('N','N',
     &                   nBs,nOr,nBs,
     &                   1.0d0,OvlS,nBs,
     &                         CMO(iCMO),nBs,
     &                   0.0d0,Aux,nBs)
             Call DGEMM_('T','N',
     &                   nOr,nOr,nBs,
     &                   1.0d0,CMO(iCMO),nBs,
     &                         Aux,nBs,
     &                   0.0d0,OvlS,nOr)
*            get largest non zero off diag element
             Do i=1,nOr
               Do j = 1, i-1
                 iOff = (j-1)*nOr+i
                 OffMx = Max(OffMx,Abs(OvlS(iOff)))
               End Do
             End Do
*            get diag element most different from one
             Do i=1,nOr
               iOff = (i-1)*nOr+i
               DgNo1 = Max(DgNo1,Abs(OvlS(iOff)-One))
             End Do
*            check, if orthogonality violated
             If ((OffMx.gt.OrtThr).OR.(DgNo1.gt.OrtThr)) Then
               termin=.TRUE.
               If (OffMx.gt.OrtThr) Then
*                off diag element too large. Now we have time, since
*                program will terminate anyway -> Go through matrix
*                again, this time we want element indices...
                 Do i=1,nOr
                   Do j = 1, i-1
                     iOff = (j-1)*nOr+i
                     If (Abs(OvlS(iOff)).ge.OffMx) Then
*                      found
                       iOffMx=i
                       jOffMx=j
                       OffMx=OvlS(iOff)
*                      exit loop
                       GoTo 100
                     End If
                   End Do
                 End Do
  100            Continue
                 call WarningMessage(0, 'Orthogonality violated')
c                 Write(6,*)' iSym =',iSym
c                Write(6,*)' largest off diag element:',
c     &                     ' [',iOffMx,',',jOffMx,']',' = ',OffMx
               End If
               If (DgNo1.gt.OrtThr) Then
*                diag element too different from One
                 Do i=1,nOr
                   iOff = (i-1)*nOr+i
                   If (Abs(OvlS(iOff)-One).ge.DgNo1) Then
*                    found
                     iDgNo1=i
                     DgNo1=OvlS(iOff)
*                    exit loop
                     GoTo 110
                   End If
                 End Do
  110            Continue
                 call WarningMessage(0, 'Orthogonality violated')
c                 Write(6,*)
c                 Write(6,*)' ***** Orthogonality violated *****'
c                 Write(6,*)' iSym =',iSym
c                 Write(6,*)' diag element most different from 1.0:',
c     &                     ' [',iDgNo1,',',iDgNo1,']',' = ',DgNo1
               End If
             End If
           End If
        End If
        ij   = ij   + nBs*(nBs + 1)/2
        iCMO = iCMO + nBs*nOr
      End Do
*
      Call mma_deallocate(Aux)
      Call mma_deallocate(OvlS)
      If (termin) Call Abend
*
      Return
      End
