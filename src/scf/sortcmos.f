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
      Subroutine SorbCMOs(CMO,nCMO,nD,EOrb,Occ,nnB,nBas,nOrb,nSym)
      Implicit Real*8 (a-h,o-z)
      Real*8 CMO(nCMO,nD), EOrb(nnB,nD), Occ(nnB,nD)
      Integer nBas(nSym), nOrb(nSym)
*define _DEBUG_
#ifdef _DEBUG_
      Do iD = 1, nD
         iOff=1
         jOff=1
         Do iSym = 1, nSym
            Call RecPrt('Occ','(10F6.2)',Occ(iOff,iD),1,nOrb(iSym))
            Call RecPrt('CMO','(10F6.2)',CMO(jOff,iD),
     &                  nBas(iSym),nOrb(iSym))
            iOff = iOff + nOrb(iSym)
            jOff = jOff + nBas(iSym)*nOrb(iSym)
         End Do
      End Do
#endif
*
*     Sort orbitals into
*     1) occupied and virtual orbitals
*     2) within each block sort according to the orbital energy
*
      Do iD = 1, nD
*        Write (6,*)
*        Write (6,*) 'iD=',iD
*        Write (6,*)
         iOff1= 0
         iOff2= 1
         Do iSym = 1, nSym
*
            nOcc = 0
            If (nOrb(iSym).eq.0) Go To 100
*
*           Sort first the orbitals according to the occupation numbers.
*
            Do iOrb = 1, nOrb(iSym)-1
*
               Occ_i = Occ(iOff1+iOrb,iD)
*              Write (6,*) 'Occ_i,iOrb=',Occ_i,iOrb
               kOrb = 0
               Do jOrb = iOrb + 1, nOrb(iSym)
                  Occ_j = Occ(iOff1+jOrb,iD)
*                 Write (6,*) 'Occ_j,jOrb=',Occ_j,jOrb
                  If (Occ_j.gt.Occ_i) Then
                     Occ_i=Occ_j
                     kOrb = jOrb
                  End If
               End Do
*              Write (6,*) 'kOrb=',kOrb
               If (kOrb.ne.0) Then
                  Occ_i = Occ(iOff1+iOrb,iD)
                  Occ(iOff1+iOrb,iD)=Occ(iOff1+kOrb,iD)
                  Occ(iOff1+kOrb,iD)=Occ_i
                  EOrb_i = EOrb(iOff1+iOrb,iD)
                  EOrb(iOff1+iOrb,iD)=EOrb(iOff1+kOrb,iD)
                  EOrb(iOff1+kOrb,iD)=EOrb_i
                  Call DSwap_(nBas(iSym),
     &                        CMO(iOff2+(iOrb-1)*nBas(iSym),iD),1,
     &                        CMO(iOff2+(kOrb-1)*nBas(iSym),iD),1)
               End If
*
               If (Occ(iOff1+iOrb,iD).ne.0.0D0) nOcc = nOcc + 1
*
            End Do
*
*           Now sort the block with respect to the lowest possible
*           orbital energy.
*
            Do iBlock = 1, 2
*
               If (iBlock.eq.1) Then
                  iStr=1
                  iEnd=nOcc
               Else
                  iStr=nOcc+1
                  iEnd=nOrb(iSym)
               End If
*
               Do iOrb = iStr, iEnd-1
*
                  EOrb_i = EOrb(iOff1+iOrb,iD)
                  kOrb = 0
                  Do jOrb = iOrb + 1, iEnd
                     EOrb_j = EOrb(iOff1+jOrb,iD)
                     If (EOrb_j.lt.EOrb_i) Then
                        EOrb_i=EOrb_j
                        kOrb = jOrb
                     End If
                  End Do
                  If (kOrb.ne.0) Then
                     Occ_i = Occ(iOff1+iOrb,iD)
                     Occ(iOff1+iOrb,iD)=Occ(iOff1+kOrb,iD)
                     Occ(iOff1+kOrb,iD)=Occ_i
                     EOrb_i = EOrb(iOff1+iOrb,iD)
                     EOrb(iOff1+iOrb,iD)=EOrb(iOff1+kOrb,iD)
                     EOrb(iOff1+kOrb,iD)=EOrb_i
                     Call DSwap_(nBas(iSym),
     &                           CMO(iOff2+(iOrb-1)*nBas(iSym),iD),1,
     &                           CMO(iOff2+(kOrb-1)*nBas(iSym),iD),1)
                  End If
*
               End Do
            End Do     ! iBlock
*
 100        Continue
*
            iOff1 = iOff1 + nOrb(iSym)
            iOff2 = iOff2 + nBas(iSym)*nOrb(iSym)
         End Do        ! iSym
      End Do           ! iD
#ifdef _DEBUG_
      Do iD = 1, nD
         iOff=1
         jOff=1
         Do iSym = 1, nSym
            Call RecPrt('Occ','(10F6.2)',Occ(iOff,iD),1,nOrb(iSym))
            Call RecPrt('CMO','(10F6.2)',CMO(jOff,iD),
     &                  nBas(iSym),nOrb(iSym))
            iOff = iOff + nOrb(iSym)
            jOff = jOff + nBas(iSym)*nOrb(iSym)
         End Do
      End Do
#endif
*
      Return
      End
