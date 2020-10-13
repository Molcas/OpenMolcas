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
      Subroutine Init_LoProp(nSym,nBas,nOrb,CoC,nAtoms,ipC,
     &                       ipQ_Nuc,ip_ANr,ip_Type,ip_Center,
     &                       nSize,nBas1,nBas2,nBasMax,ipP,ipPInv)
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Integer Occ, Vir
      Parameter(Occ=1,Vir=0)
      Integer nBas(8), nOrb(8)
      Real*8 CoC(3)
      Logical lOrb
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      Call qpg_iarray('nOrb',lOrb,iDum)
      If (lOrb) Then
         Call Get_iArray('nOrb',nOrb,nSym)
      Else
         Call ICopy(nSym,nBas,1,nOrb,1)
      End If
*                                                                      *
************************************************************************
*                                                                      *
      nSize=0
      nBas1=0
      nBas2=0
      nBasMax=0
      Do iSym = 1, nSym
         nSize = nSize + nBas(iSym)*(nBas(iSym)+1)/2
         nBas1 = nBas1 + nBas(iSym)
         nBas2 = nBas2 + nBas(iSym)**2
         nBasMax = Max(nBasMax,nBas(iSym))
      End Do
      nSize=nSize + 4
*                                                                      *
************************************************************************
*                                                                      *
*---- Center of Charge
      Call Get_dArray('Center of Charge',CoC,3)
*
*---- List coordinates of Coordinates
      Call Get_iScalar('LP_nCenter',nAtoms)
      Call Allocate_Work(ipC,3*nAtoms)
      Call Get_dArray('LP_Coor',Work(ipC),3*nAtoms)
      Call Allocate_Work(ipQ_Nuc,nAtoms)
*
*---- Effective charge at each center
      Call Get_dArray('LP_Q',Work(ipQ_Nuc),nAtoms)
*
*---- Atom number of each center
      Call Allocate_iWork(ip_ANr,nAtoms)
      Call Get_iArray('LP_A',iWork(ip_ANr),nAtoms)
*
*---- Pick up information of orbital type. Occ/Vir
      Call Allocate_iWork(ip_type,nbas1)
      Call Get_iArray('Orbital Type',iWork(ip_type),nBas1)
      Do i = ip_type, ip_type+nBas1-1
         If (iWork(i).ne.Occ.and.iWork(i).ne.Vir) Then
            Write (6,*) 'Orbital type vector is corrupted!'
            Call Abend()
         End If
      End Do
*
*---- Pick up index array of which center a basis function belong.
      Call Allocate_iWork(ip_center,nbas1)
      Call Get_iArray('Center Index',iWork(ip_center),nBas1)
*
#ifdef _DEBUGPRINT_
      Write (6,*) '******* LoProp Debug Info *******'
      Call RecPrt('Coordinates',' ',Work(ipC),3,nAtoms)
      Call RecPrt('Charges',' ',Work(ipQ_Nuc),1,nAtoms)
      Write (6,*) 'Atom Nr:',(iWork(i),i=ip_ANr,ip_ANr+nAtoms-1)
      Do iBas = 1, nBas1
         If (iWork(ip_type+iBas-1).eq.Occ) Then
            Write (6,'(A,I3,A)') 'Basis function ',iBas,
     &            ' is occupied'
         Else
            Write (6,'(A,I3,A)') 'Basis function ',iBas,
     &            ' is virtual'
         End If
      End Do
      Write (6,*)
      Do iBas = 1, nBas1
         Write (6,'(A,I3,A,I3)') 'Basis function ',iBas,' belongs'//
     &         ' to center ',iWork(ip_center+iBas-1)
      End Do
      Write (6,*) '*********************************'
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     In case of symmetry we need the desymmetrization matrix
*
      If (nSym.eq.1) Go To 99
*
      Call Allocate_Work(ipP,nbas1**2)
      Call Allocate_Work(ipPInv,nbas1**2)
      Call Get_dArray('SM',Work(ipP),nbas1**2)
#ifdef _DEBUGPRINT_
      Call RecPrt('SM',' ',Work(ipP),nbas1,nbas1)
#endif
      Call MINV(Work(ipP),Work(ipPInv),ISING,DET,nBas1)
#ifdef _DEBUGPRINT_
      Call RecPrt('SMInv',' ',Work(ipPInv),nbas1,nbas1)
#endif
      Call DGeTMi(Work(ipPInv),nbas1,nbas1)
*                                                                      *
************************************************************************
*                                                                      *
 99   Return
      End
