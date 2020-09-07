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
      Subroutine PotNuc_nad(nSym,nAtoms,ReCharge,ZRE_nad)
************************************************************************
*                                                                      *
*     purpose: Computes NAD part of the Nuclear Repulsion Energy.      *
*              An array of reference charges (ReCharge) is used to     *
*              identify which atoms were alternating as ghosts         *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "angstr.fh"
      Integer nSym, nAtoms
      Real*8  ReCharge(nAtoms)
      Integer iGen(3), iCoSet(0:7,0:7), iStab(0:7), iOper(0:7)
*----------------------------------------------------------------------*
*     Prologue                                                         *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*     Read symm. oper per symm. species                                *
*----------------------------------------------------------------------*
      Call Get_iArray('Symmetry operations',iOper,nSym)
*----------------------------------------------------------------------*
*     Read atom Charges                                                *
*----------------------------------------------------------------------*
      Call GetMem('Charge','ALLO','REAL',lw0,8*nAtoms)
      Call Get_dArray('Effective nuclear Charge',Work(lw0),nAtoms)
*----------------------------------------------------------------------*
*     Read coordinates of atoms                                        *
*----------------------------------------------------------------------*
      Call GetMem('Coor','ALLO','REAL',lw1,3*8*nAtoms)
      Call Get_dArray('Unique Coordinates',Work(lw1),3*nAtoms)
*----------------------------------------------------------------------*
*     Apply the symmetry operations                                    *
*----------------------------------------------------------------------*
      nGen=0
      If (nSym.eq.2) nGen=1
      If (nSym.eq.4) nGen=2
      If (nSym.eq.8) nGen=3
      If (nGen.ge.1) iGen(1)=iOper(1)
      If (nGen.ge.2) iGen(2)=iOper(2)
      If (nGen.ge.3) iGen(3)=iOper(4)
*
      iAll_Atom=0
      MaxDCR=0
      iAll_Atom=nAtoms
      Do iAtom = 1, nAtoms
         iChAtom=iChxyz(Work(lw1+(iAtom-1)*3),iGen,nGen)
         Call Stblz(iChAtom,iOper,nSym,nStab,iStab,MaxDCR,iCoSet)
         nCoSet=nSym/nStab
         XOld=Work(lw1+(iAtom-1)*3 )
         YOld=Work(lw1+(iAtom-1)*3+1)
         ZOld=Work(lw1+(iAtom-1)*3+2)
         Charge_=Work(lw0+iAtom-1)

*
         Do iCo = 1, nCoSet-1
*
            Work(lw0+iAll_Atom)=Charge_
*
            iAll_Atom = iAll_Atom + 1
            Work(lw1+(iAll_Atom-1)*3  )=
     &           XOld*DBLE(iPhase(1,iCoSet(iCo,0)))
            Work(lw1+(iAll_Atom-1)*3+1)=
     &           YOld*DBLE(iPhase(2,iCoSet(iCo,0)))
            Work(lw1+(iAll_Atom-1)*3+2)=
     &           ZOld*DBLE(iPhase(3,iCoSet(iCo,0)))
*
         End Do
*
      End Do
*
*
*----------------------------------------------------------------------*
*     Compute NAD part of the nuclear repulsion energy                 *
*----------------------------------------------------------------------*
      ZRE_nad=0.0d0
*
      If (ReCharge(1).gt.0.0d0) Then

         Do jAt=0,iAll_Atom-1
           pCharge=Work(lw0+jAt)
           If (pCharge.gt.0.0d0) Then
              Do iAt=0,jAt-1  ! loop downwards
                 kAt=mod(iAt+1,nAtoms)
                 If (kAt.eq.0) kAt=nAtoms
                 qCharge=ReCharge(kAt)
                 If (qCharge.gt.0.0d0) Then
                    Xpq = Work(lw1+3*iAt) - Work(lw1+3*jAt)
                    Ypq = Work(lw1+3*iAt+1) - Work(lw1+3*jAt+1)
                    Zpq = Work(lw1+3*iAt+2) - Work(lw1+3*jAt+2)
                    Rpq = sqrt(Xpq**2+Ypq**2+Zpq**2)
                    pq_rep = pCharge*qCharge/Rpq
                    ZRE_nad = ZRE_nad + pq_rep
                 EndIf
              End Do
           EndIf
         End Do

      Else

         Do jAt=0,iAll_Atom-1
           pCharge=Work(lw0+jAt)
           If (pCharge.gt.0.0d0) Then
              Do iAt=jAt+1,iAll_Atom-1   ! loop upwards
                 kAt=mod(iAt+1,nAtoms)
                 If (kAt.eq.0) kAt=nAtoms
                 qCharge=ReCharge(kAt)
                 If (qCharge.gt.0.0d0) Then
                    Xpq = Work(lw1+3*iAt) - Work(lw1+3*jAt)
                    Ypq = Work(lw1+3*iAt+1) - Work(lw1+3*jAt+1)
                    Zpq = Work(lw1+3*iAt+2) - Work(lw1+3*jAt+2)
                    Rpq = sqrt(Xpq**2+Ypq**2+Zpq**2)
                    pq_rep = pCharge*qCharge/Rpq
                    ZRE_nad = ZRE_nad + pq_rep
                 EndIf
              End Do
           EndIf
         End Do

      EndIf

#ifdef _DEBUG_
*----------------------------------------------------------------------*
*     Print coordinates of the system  / ZRE_nad energy                *
*----------------------------------------------------------------------*
      Write(6,*)
      Write(6,'(6X,A)')'Atoms cartesian coordinates in Angstrom:'
      Write(6,'(6X,A)')'-----------------------------------------------'
      Write(6,'(6X,A)')'No.  Charge A/B      X         Y         Z     '
      Write(6,'(6X,A)')'-----------------------------------------------'
      Do iAt=0,iAll_Atom-1
        kAt=mod(iAt+1,nAtoms)
        Write(6,'(4X,I4,2X,F4.0,1X,F4.0,2X,3F10.5)')
     &  iAt+1,Work(lw0+iAt),Recharge(kAt),
     &  Angstr*Work(lw1+3*iAt),
     &  Angstr*Work(lw1+3*iAt+1),
     &  Angstr*Work(lw1+3*iAt+2)
      End Do
      Write(6,'(6X,A)')'-----------------------------------------------'
      Write(6,'(6X,A,F12.6)')'Nuclear repulsion energy (NAD) =',ZRE_nad
      Write(6,*)
      Write(6,*)
#endif
*----------------------------------------------------------------------*
*     Normal exit                                                      *
*----------------------------------------------------------------------*
      Call GetMem('Charge','Free','REAL',lw0,8*nAtoms)
      Call GetMem('Coor','FREE','REAL',lw1,3*8*nAtoms)
      Return
      End
