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
      Character*40 OrbTit
#include "WrkSpc.fh"

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

      Call GetMem('Scr','Allo','Real',ip_Scr,2*lScr)
      ip_U=ip_Scr+lScr

      Call GetMem('WYZR','Allo','Real',ip_w,4*l_O2)
      ip_y = ip_w + l_O2
      ip_z = ip_y + l_O2
      ip_r = ip_z + l_O2

      Call GetMem('S','Allo','Real',ip_S,2*l_S)
      ip_X=ip_S+l_S

      Call GetOvlp_Localisation(Work(ip_S),'Sqr',nBas,nSym)

      Lu=12
      Call RdVec('INPORB',Lu,'C',nSym,nBas,nBas,Work(ip_X),Dummy,Dummy,
     &                    iDummy,OrbTit,1,iErr)

      write(6,*)
      write(6,*) '      Thouless singles amplitudes from: '
      write(6,*) '      '//OrbTit
      write(6,*)

      iOff=0
      kOff=0
      Do iSym=1,nSym
         jp_S = ip_S + iOff
         jp_C = 1 + iOff + nBas(iSym)*nFro(iSym)
         jp_X = ip_X + iOff + nBas(iSym)*nFro(iSym)
         nOrb=nOcc(iSym)+nSsh(iSym)

         Call GetUmat_T1(Work(ip_U),CMO(jp_C),Work(jp_S),Work(jp_X),
     &                   Work(ip_Scr),lScr,nBas(iSym),
     &                   nOrb,nOcc(iSym))

         iU=ip_U
         Do j=1,nOcc(iSym)
            ifr=ip_U+nOrb*(j-1)
            ito=ip_Scr+nOcc(iSym)*(j-1)
            call dcopy_(nOcc(iSym),Work(ifr),1,Work(ito),1)
            jU=ifr+nOcc(iSym)
            Do i=1,nSsh(iSym)
               Work(iU)=Work(jU)
               iU=iU+1
               jU=jU+1
            End Do
         End Do

c --- SVD of U in the oo space:   U = Y * w * Z'

         Call SVD(nOcc(iSym),nOcc(iSym),nOcc(iSym),Work(ip_Scr),
     &            Work(ip_w),.true.,Work(ip_y),.true.,Work(ip_z),
     &            ierr,Work(ip_r))

         If (ierr.ne.0) Then
            write(6,*)
            write(6,*) ' *** Warning: SVD failed to get singval: ',ierr
            write(6,*) ' *** Located in Thouless_T1 -- call to SVD .'
            write(6,*)
            write(6,*) ' omega= ',(Work(ip_w+k),k=0,nOcc(iSym)-1)
         EndIf

         Call FZero(Work(ip_r),nOcc(iSym)**2)
         Do k=1,nOcc(iSym)
            omega=Work(ip_w+k-1)
            kk=nOcc(iSym)*(k-1)+k-1
            If (omega.gt.1.0d-8) Then
               Work(ip_r+kk)=1.0d0/omega
            EndIf
         End Do

c --- Compute U^-1 = Z * w^-1 * Y'

         Call DGEMM_('N','T',nOcc(iSym),nOcc(iSym),nOcc(iSym),
     &                           1.0d0,Work(ip_r),nOcc(iSym),
     &                                 Work(ip_y),nOcc(iSym),
     &                           0.0d0,Work(ip_w),nOcc(iSym))

         Call DGEMM_('N','N',nOcc(iSym),nOcc(iSym),nOcc(iSym),
     &                           1.0d0,Work(ip_z),nOcc(iSym),
     &                                 Work(ip_w),nOcc(iSym),
     &                           0.0d0,Work(ip_Scr),nOcc(iSym))

         jp_T = 1 + kOff
         Call DGEMM_('T','T',nOcc(iSym),nSsh(iSym),nOcc(iSym),
     &                           1.0d0,Work(ip_Scr),nOcc(iSym),
     &                                 Work(ip_U),nSsh(iSym),
     &                           0.0d0,T1amp(jp_T),nOcc(iSym))

         iOff = iOff + nBas(iSym)**2
         kOff = kOff + nOcc(iSym)*nSsh(iSym)
      End Do

      Call GetMem('S','Free','Real',ip_S,2*l_S)
      Call GetMem('WYZR','Free','Real',ip_w,4*l_O2)
      Call GetMem('Scr','Free','Real',ip_Scr,2*lScr)
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
