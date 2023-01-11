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
* Copyright (C) 1992, Martin Schuetz                                   *
*               1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
*               2017,2022, Roland Lindh                                *
************************************************************************
      SubRoutine SOiniH()
************************************************************************
*                                                                      *
*     purpose: generate initial Hessian (diagonal) from                *
*              orbital energies (for second order update)              *
*                                                                      *
************************************************************************
      use Orb_Type, only: OrbType
      use InfSCF, only: nSym, nFro, nOrb, nOcc, RGEK
*     use SCF_Arrays, only: HDiag, FockMO, EOrb, CMO_Ref
      use SCF_Arrays, only: HDiag, FockMO
      use Constants, only: Zero, Four
      Implicit None
*
*     declaration local variables
      Integer iD, nD
      Integer iSym,iOcc,iVir,ioffs,iOff_H,nOccmF,nOrbmF, iOff_F
      Integer jOcc, jVir
      Real*8, Parameter:: Hii_Min=0.05D0
      Real*8, Parameter:: Hii_Max=1.00D0
      Real*8, Pointer:: Fock(:,:)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*     Set the array to silly large values. In the case of UHF these
*     will remain but should not make any difference. They are actully
*     needed to make the rs-rfo code work.
*
*     Compute the diagonal values of the Fock matrix, stored in EOrb.
!     Call Mk_EOrb(CMO_Ref,Size(CMO_Ref,1),Size(CMO_Ref,2))

      nD   =Size(FockMO,2)
      HDiag(:)=1.0D+99
*
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
      Write (6,*) 'nD=',nD
      Do iD = 1, nD
         Write (6,*) 'iD=',iD
         Write (6,'(A,8I3)') 'nOcc',(nOcc(iSym,iD),iSym=1,nSym)
      End Do
      Write (6,'(A,8I3)') 'nFro',(nFro(iSym),iSym=1,nSym)
      Write (6,'(A,8I3)') 'nOrb',(nOrb(iSym),iSym=1,nSym)
#endif
      iOff_H=1
      Do iD = 1, nD
*
         iOffs=0
         iOff_F=0
         Do iSym=1,nSym
*
*            loop over all occ orbitals in sym block
*
             iOffs=iOffs+nFro(iSym)
             ! number of Occupied, excluding frozen
             nOccmF=nOcc(iSym,iD)-nFro(iSym)
             ! number of Orbitals, excluding frozen
             nOrbmF=nOrb(iSym)-nFro(iSym)

             Fock(1:nOrb(iSym),1:nOrb(iSym)) =>
     &            FockMO(iOff_F+1:iOff_F+nOrb(iSym)**2,iD)
*
             Do iOcc= ioffs+1, ioffs+nOccmF
                jOcc = iOcc - iOffs
*
*               loop over all virt orbitals in sym block
*
                Do iVir=ioffs+nOccmF+1,ioffs+nOrbmF
                   jVir = iVir - iOffs
*
                   If (OrbType(iVir,iD).eq.OrbType(iOcc,iD)) Then

                      HDiag(iOff_H)=
     &                  Four*(Fock(jVir,jVir)-Fock(jOcc,jOcc))/DBLE(nD)

                      If (HDiag(iOff_H)<Zero.and..NOT.RGEK) Then
                         Write (6,*) 'SOIniH: Hii<0.0, Hii=',
     &                                HDiag(iOff_H)
                         HDiag(iOff_H)=Max(Hii_Max,Abs(HDiag(iOff_H)))
                      Else If (Abs(HDiag(iOff_H)).lt.Hii_Min) Then
*                        Write (6,*) 'SOIniH: Abs(Hii)<0.05, Hii=',
*    &                                HDiag(iOff_H)
*                        Write (6,*) 'jVir,jOcc=',jVir,jOcc
*                        Write (6,*) 'Fock(jOcc,jOcc)=',Fock(jOcc,jOcc)
*                        Write (6,*) 'Fock(jVir,jVir)=',Fock(jVir,jVir)
                         HDiag(iOff_H)=Hii_Min
                      End If
                   End If
*
                   iOff_H=iOff_H+1
*
                End Do  ! iVir
*
             End Do     ! iOcc
*
             Nullify(Fock)
             iOff_F = iOff_F + nOrb(iSym)**2
             iOffs=iOffs+nOrbmF
*
          End Do ! iSym
      End Do ! iD

#ifdef _DEBUGPRINT_
      Call RecPrt('HDiag',' ',HDiag(:),1,Size(HDiag))
#endif
      Return
      End
