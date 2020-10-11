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
      Subroutine Sorb2CMOs(CMO,nCMO,nD,Occ,nnB,nBas,nOrb,nSym,OrbType)
      Implicit Real*8 (a-h,o-z)
      Real*8 CMO(nCMO,nD), Occ(nnB,nD)
      Integer nBas(nSym), nOrb(nSym), OrbType(nnB,nD)
*define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
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
                  iTmp=OrbType(iOff1+iOrb,iD)
                  OrbType(iOff1+iOrb,iD)=OrbType(iOff1+kOrb,iD)
                  OrbType(iOff1+kOrb,iD)=iTmp
*
                  Occ_i = Occ(iOff1+iOrb,iD)
                  Occ(iOff1+iOrb,iD)=Occ(iOff1+kOrb,iD)
                  Occ(iOff1+kOrb,iD)=Occ_i
                  Call DSwap_(nBas(iSym),
     &                        CMO(iOff2+(iOrb-1)*nBas(iSym),iD),1,
     &                        CMO(iOff2+(kOrb-1)*nBas(iSym),iD),1)
               End If
*
               If (Occ(iOff1+iOrb,iD).ne.0.0D0) nOcc = nOcc + 1
*
            End Do
*
            iOff1 = iOff1 + nOrb(iSym)
            iOff2 = iOff2 + nBas(iSym)*nOrb(iSym)
         End Do        ! iSym
      End Do           ! iD
#ifdef _DEBUGPRINT_
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
