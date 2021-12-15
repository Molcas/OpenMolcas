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
      Subroutine Print_T_Values(T_Values,iT_Sets,iANr,EC,Bond_Threshold,
     &                          nAtoms,nij,Standard,iWarnings,
     &                          Num_Warnings,iPrint)
      Implicit Real*8 (A-H,O-Z)
#include "Molcas.fh"
#include "real.fh"
      Real*8 T_Values(nij),EC(3,nij)
      Integer iT_Sets(nij),iANr(nAtoms),iWarnings(nij)
      Character*(LENIN) AtomLbl(MxAtom)
      Character*(LENIN4) AtomLbl4(MxAtom)
      Character*17 BondLbl
      Parameter (iLength=25)
      Character*(iLength) Warning
      Logical Standard
*
* Print header
*
      Call Get_cArray('LP_L',AtomLbl4,(LENIN4)*nAtoms)
      Do i=1,nAtoms
       AtomLbl(i)(1:LENIN)=AtomLbl4(i)(1:LENIN)
      EndDo
      ij = 0
      Write(6,*)
      If (Num_Warnings .gt. 0) Then
         Write(6,'(A,I3,A)')
     &              'During optimization of the expansion centers ',
     &              Num_Warnings, ' warnings were encountered.'
         Write(6,*)
         Write(6,'(A)') ' iAtom   jAtom   Atom(s)          Factor'
     &               // '   Bragg-Slater      t      Warning'
      Else
         Write(6,'(A)') ' iAtom   jAtom   Atom(s)          Factor'
     &               // '   Bragg-Slater      t'
      End If
*
* Print informations for the atoms
*
      Do iAtom = 1, nAtoms
         ii = iAtom*(iAtom+1)/2
         If ((iT_sets(ii) .eq. 1 .OR. iPrint .ge. 2)) Then
            Write(BondLbl, '(A)') AtomLbl(iAtom)
            Last_NonBlank = 0
            j = 0
            Do i = 1, 17
               If (BondLbl(i:i) .ne. ' ') Then
                  j = j + 1
                  BondLbl(j:j) = BondLbl(i:i)
                  Last_NonBlank = j
               End If
            End Do
            Do i = Last_NonBlank+1,17
               BondLbl(i:i) = ' '
            End Do
*
            If (Num_Warnings .gt. 0) Then
               Call Warnings(iWarnings(ii),Warning,iLength)
               If (iT_sets(ii) .eq. 1) Then
                  Write(6,'(1X,I3,5X,8X,A17,24X,F7.4,3X,A)')
     &               iAtom,BondLbl,T_Values(ij),Warning
               Else If (Standard) Then
                  Write(6,'(1X,I3,5X,8X,A17,24X,A,3X,A)')
     &               iAtom,BondLbl,'Standard',Warning
               Else
                  Write(6,'(1X,I3,5X,8X,A17,24X,A,3X,A)')
     &               iAtom,BondLbl,'Skipped',Warning
               End If
            Else
               If (iT_sets(ii) .eq. 1) Then
                  Write(6,'(1X,I3,5X,8X,A17,24X,F7.4)')
     &               iAtom,BondLbl,T_Values(ij)
               Else If (Standard) Then
                  Write(6,'(1X,I3,5X,8X,A17,24X,A)')
     &               iAtom,BondLbl,'Standard'
               Else
                  Write(6,'(1X,I3,5X,8X,A17,24X,A)')
     &               iAtom,BondLbl,'Skipped'
               End If
            End If
         End If
      End Do
*
      Write(6,'(79A)') ('-',i=1,79)
*
* Print informations for the bonds
*
      Do iAtom = 1, nAtoms
         ii = iAtom*(iAtom+1)/2
         Do jAtom = 1, iAtom-1
            ij = iAtom*(iAtom-1)/2+jAtom
            If ((iT_sets(ij) .eq. 1 .OR. iPrint .ge. 2)) Then
               jj = jAtom*(jAtom+1)/2
               Factor      = Zero
               Write(BondLbl,'(3A)') AtomLbl(iAtom),'-',
     &                               AtomLbl(jAtom)
               Radius_i    = Bragg_Slater(iANr(iAtom))
               Radius_j    = Bragg_Slater(iANr(jAtom))
               Bond_Length = Sqrt((EC(1,ii)-EC(1,jj))**2
     &                           +(EC(2,ii)-EC(2,jj))**2
     &                           +(EC(3,ii)-EC(3,jj))**2)
               Bond_Max    = Bond_Threshold*(Radius_i+Radius_j)
               Factor      = Bond_Length/Bond_Max
               bs_t        = Radius_i/(Radius_i+Radius_j)-Half
*
               Last_NonBlank = 0
               j = 0
               Do i = 1, 17
                  If (BondLbl(i:i) .ne. ' ') Then
                     j = j + 1
                     BondLbl(j:j) = BondLbl(i:i)
                     Last_NonBlank = j
                  End If
               End Do
               Do i = Last_NonBlank+1,17
                  BondLbl(i:i) = ' '
               End Do
*
               If (Num_Warnings .gt. 0) Then
                  Call Warnings(iWarnings(ij),Warning,iLength)
                  If (iT_sets(ij) .eq. 1) Then
               Write(6,'(1X,I3,5X,I3,5X,A17,F6.3,5X,F7.4,6X,F7.4,3X,A)')
     &              iAtom,jAtom,BondLbl,Factor,BS_t,T_Values(ij),Warning
                  Else If (Standard) Then
                  Write(6,'(1X,I3,5X,I3,5X,A17,F6.3,5X,F7.4,6X,A,3X,A)')
     &                iAtom,jAtom,BondLbl,Factor,BS_t,'Standard',Warning
                  Else
                  Write(6,'(1X,I3,5X,I3,5X,A17,F6.3,5X,F7.4,6X,A,3X,A)')
     &                 iAtom,jAtom,BondLbl,Factor,BS_t,'Skipped',Warning
                  End If
               Else
                  If (iT_sets(ij) .eq. 1) Then
                   Write(6,'(1X,I3,5X,I3,5X,A17,F6.3,5X,F7.4,6X,F7.4)')
     &                  iAtom,jAtom,BondLbl,Factor,BS_t,T_Values(ij)
                  Else If (Standard) Then
                     Write(6,'(1X,I3,5X,I3,5X,A17,F6.3,5X,F7.4,6X,A)')
     &                  iAtom,jAtom,BondLbl,Factor,BS_t,'Standard'
                  Else
                     Write(6,'(1X,I3,5X,I3,5X,A17,F6.3,5X,F7.4,6X,A)')
     &                  iAtom,jAtom,BondLbl,Factor,BS_t,'Skipped'
                  End If
               End If
            End If
         End Do
      End Do
*
      Return
      End
