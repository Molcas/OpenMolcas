************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Thouless_T1(CMO,nSym,nBas,nFro,nOcc,nSsh,T1amp)

      Implicit Real*8 (a-h,o-z)
      Integer nSym, nBas(nSym), nFro(nSym), nOcc(nSym), nSsh(nSym)
      Real*8 CMO(*), T1amp(*)
      Character(LEN=40) OrbTit
      Real*8 Dummy(1)
      Integer iDummy(1)
#include "stdalloc.fh"

      Real*8, Allocatable:: S(:), X(:), Scr(:), U(:)
      Real*8, Allocatable:: W(:), Y(:), Z(:), R(:)

C
C     Compute the T1 amplitudes according to Thouless formula
C     --------------------------------------------------------

      l_S = nBas(1)**2
      lScr=nBas(1)*(nFro(1)+nOcc(1))
      l_O=nOcc(1)
      Do iSym = 2,nSym
         l_S = l_S + nBas(iSym)**2
         lScr=Max(lScr,nBas(iSym)*(nFro(iSym)+nOcc(iSym)))
         l_O=Max(l_O,nOcc(iSym))
      End Do
      l_O2=l_O**2

      Call mma_allocate(Scr,lScr,Label='Scr')
      Call mma_allocate(U,lScr,Label='U')

      Call mma_allocate(W,l_O2,Label='W')
      Call mma_allocate(Y,l_O2,Label='Y')
      Call mma_allocate(Z,l_O2,Label='Z')
      Call mma_allocate(R,l_O2,Label='R')

      Call mma_allocate(S,l_S,Label='S')
      Call mma_allocate(X,l_S,Label='X')

      Call GetOvlp_Localisation(S,'Sqr',nBas,nSym)

      Lu=12
      Call RdVec('INPORB',Lu,'C',nSym,nBas,nBas,X,Dummy,Dummy,
     &                    iDummy,OrbTit,1,iErr)

      write(6,*)
      write(6,*) '      Thouless singles amplitudes from: '
      write(6,*) '      '//OrbTit
      write(6,*)

      iOff=0
      kOff=0
      Do iSym=1,nSym
         jp_S = 1 + iOff
         jp_C = 1 + iOff + nBas(iSym)*nFro(iSym)
         jp_X = 1 + iOff + nBas(iSym)*nFro(iSym)
         nOrb=nOcc(iSym)+nSsh(iSym)

         Call GetUmat_T1(U,CMO(jp_C),S(jp_S),X(jp_X),
     &                   Scr,lScr,nBas(iSym),
     &                   nOrb,nOcc(iSym))

         iU=1
         Do j=1,nOcc(iSym)
            ifr=1+nOrb*(j-1)
            ito=1+nOcc(iSym)*(j-1)
            call dcopy_(nOcc(iSym),U(ifr),1,Scr(ito),1)
            jU=ifr+nOcc(iSym)
            Do i=1,nSsh(iSym)
               U(iU)=U(jU)
               iU=iU+1
               jU=jU+1
            End Do
         End Do

c --- SVD of U in the oo space:   U = Y * w * Z'

         Call SVD(nOcc(iSym),nOcc(iSym),nOcc(iSym),Scr,
     &            W,.true.,Y,.true.,Z,ierr,R)

         If (ierr.ne.0) Then
            write(6,*)
            write(6,*) ' *** Warning: SVD failed to get singval: ',ierr
            write(6,*) ' *** Located in Thouless_T1 -- call to SVD .'
            write(6,*)
            write(6,*) ' omega= ',(W(k),k=1,nOcc(iSym))
         EndIf

         Call FZero(R,nOcc(iSym)**2)
         Do k=1,nOcc(iSym)
            omega=W(k)
            kk=nOcc(iSym)*(k-1)+k
            If (omega.gt.1.0d-8) Then
               R(kk)=1.0d0/omega
            EndIf
         End Do

c --- Compute U^-1 = Z * w^-1 * Y'

         Call DGEMM_('N','T',nOcc(iSym),nOcc(iSym),nOcc(iSym),
     &                           1.0d0,R,nOcc(iSym),
     &                                 Y,nOcc(iSym),
     &                           0.0d0,W ,nOcc(iSym))

         Call DGEMM_('N','N',nOcc(iSym),nOcc(iSym),nOcc(iSym),
     &                           1.0d0,Z,nOcc(iSym),
     &                                 W,nOcc(iSym),
     &                           0.0d0,Scr,nOcc(iSym))

         jp_T = 1 + kOff
         Call DGEMM_('T','T',nOcc(iSym),nSsh(iSym),nOcc(iSym),
     &                           1.0d0,Scr,nOcc(iSym),
     &                                 U,nSsh(iSym),
     &                           0.0d0,T1amp(jp_T),nOcc(iSym))

         iOff = iOff + nBas(iSym)**2
         kOff = kOff + nOcc(iSym)*nSsh(iSym)
      End Do

      Call mma_deallocate(X)
      Call mma_deallocate(S)
      Call mma_deallocate(R)
      Call mma_deallocate(Z)
      Call mma_deallocate(Y)
      Call mma_deallocate(W)
      Call mma_deallocate(U)
      Call mma_deallocate(Scr)
      Return
      End


      SubRoutine GetUmat_T1(U,C,S,X,Scr,lScr,nBas,nOrb1,nOrb2)
C
C     Purpose: compute transformation matrix U=C^TSX.
C
      Implicit None
      Real*8  U(*), C(*), S(*), X(*)
      Integer lScr
      Real*8  Scr(lScr)
      Integer nBas, nOrb1, nOrb2

      Character*80 Txt
      Character*10 SecNam
      Parameter (SecNam = 'GetUmat_T1')

      Real*8 d0, d1
      Parameter (d0 = 0.0d0, d1 = 1.0d0)

      Integer Need

      If (nOrb1*nOrb2.lt.1 .or. nBas.lt.1) Return

      Need = nBas*nOrb2
      If (lScr .lt. Need) Then
         Write(Txt,'(A,I9,A,I9)')
     &   'lScr =',lScr,'     Need =',Need
         Call SysAbendMsg(SecNam,
     &                   'Insufficient dimension of scratch array!',Txt)
      End If

      Call DGEMM_('N','N',nBas,nOrb2,nBas,
     &                 d1,S,nBas,
     &                    X,nBas,
     &                 d0,Scr,nBas)
      Call DGEMM_('T','N',nOrb1,nOrb2,nBas,
     &                 d1,C,nBas,
     &                    Scr,nBas,
     &                 d0,U,nOrb1)

      End
