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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
      Subroutine RPA_RdRun()
C
C     Thomas Bondo Pedersen
C
C     Read data from Runfile.
C
      Implicit None
#include "rpa_config.fh"
#include "rpa_data.fh"
      Character*9 SecNam
      Parameter (SecNam='RPA_RdRun')

      Integer  RPA_iUHF
      External RPA_iUHF

      Character*8 Model

      Logical Warn

      Integer iUHF
      Integer iSym
      Integer i

      ! Register entry
      Call qEnter(SecNam)

      ! Set type of SCF reference wave function
      ! Note: in RPA, iUHF=1 means restricted, 2 means unrestricted.
      Call Get_cArray('Relax Method',Model,8)
      Call Get_iScalar('SCF mode',iUHF)
      iUHF=iUHF+1 ! correct for convention in SCF
      If (Model(1:7).eq.'RHF-SCF') Then
         If (iUHF.eq.1) Then
            Reference='RHF'
            Warn=.false.
         Else
            Reference='UHF'
            warn=.true.
         End If
      Else If (Model(1:7).eq.'UHF-SCF') Then
         If (iUHF.eq.1) Then
            Reference='RHF'
            Warn=.true.
         Else
            Reference='UHF'
            Warn=.false.
         End If
      Else If (Model(1:6).eq.'KS-DFT') Then
         Warn=.false.
         If (iUHF.eq.1) Then
            Reference='RKS'
         Else
            Reference='UKS'
         End If
      Else If (Model(1:8).eq.'dRPA@RHF' .or.
     *         Model(1:8).eq.'SOSX@RHF') Then
         If (iUHF.eq.1) Then
            Reference='RHF'
            Warn=.false.
         Else
            Reference='UHF'
            Warn=.true.
         End If
      Else If (Model(1:8).eq.'dRPA@UHF' .or.
     *         Model(1:8).eq.'SOSX@UHF') Then
         If (iUHF.eq.1) Then
            Reference='RHF'
            Warn=.true.
         Else
            Reference='UHF'
            Warn=.false.
         End If
      Else If (Model(1:8).eq.'dRPA@RKS' .or.
     *         Model(1:8).eq.'SOSX@RKS') Then
         If (iUHF.eq.1) Then
            Reference='RKS'
            Warn=.false.
         Else
            Reference='UKS'
            Warn=.true.
         End If
      Else If (Model(1:8).eq.'dRPA@UKS' .or.
     *         Model(1:8).eq.'SOSX@UKS') Then
         If (iUHF.eq.1) Then
            Reference='RKS'
            Warn=.true.
         Else
            Reference='UKS'
            Warn=.false.
         End If
      Else
         Write(6,'(A,A)') 'Reference model from Runfile: ',Model
         Write(6,'(A,I8)') 'iUHF from Runfile:            ',iUHF-1
         Call RPA_Warn(2,'Illegal reference wave function in RPA')
         Reference='Non'
         Warn=.false.
      End If
      If (Warn) Then
         Call RPA_Warn(1,
     *                'Runfile restricted/unrestricted conflict in RPA')
         Write(6,'(A,A)') 'Reference model from Runfile: ',Model
         Write(6,'(A,I8)') 'iUHF from Runfile:            ',iUHF-1
         Write(6,'(A,A,A)') 'Assuming ',Reference,' reference!'
         Call xFlush(6)
      End If

      ! Get nuclear potential energy
      Call Get_dScalar('PotNuc',NuclearRepulsionEnergy(1))

      ! Get DFT functional
      If (Reference(2:3).eq.'KS') Then
         Call Get_cArray('DFT functional',DFTFunctional,16)
      Else
         DFTFunctional='Hartree-Fock'
      End If

      ! Get number of irreps
      Call Get_iScalar('nSym',nSym)
      If (nSym.lt.1 .or. nSym.gt.8) Then
         Call RPA_Warn(3,'nSym out of bonds in RPA')
      End If

      ! Get iUHF (RPA convention: 1->RHF, 2->UHF)
      !          (as opposed to SCF: 0->RHF, 1->UHF)
      iUHF=RPA_iUHF()

      ! Get number of basis functions, orbitals, occupied,
      ! frozen (in SCF), deleted
      Call Get_iArray('nBas',nBas,nSym)
      Call Get_iArray('nOrb',nOrb,nSym)
      Call Get_iArray('nDel',nDel(1,1),nSym)
      Call Get_iArray('nFro',nFro(1,1),nSym)
      Call Get_iArray('nIsh',nOcc(1,1),nSym)
      If (iUHF.eq.2) Then
         ! unrestricted: read data for beta spin
         Call Get_iArray('nIsh_ab',nOcc(1,2),nSym)
      End If

      ! Check for orbitals frozen in SCF and check consistency.
      Do iSym=1,nSym
         If (nFro(iSym,1).ne.0) Then
            Write(6,'(A,8I8)') 'nFro=',(nFro(i,1),i=1,nSym)
            Call RPA_Warn(4,
     &                    SecNam//': Some orbitals were frozen in SCF!')
         End If
         If (nDel(iSym,1).ne.(nBas(iSym)-nOrb(iSym))) Then
            Write(6,'(A,8I8)') 'nBas=     ',(nBas(i),i=1,nSym)
            Write(6,'(A,8I8)') 'nOrb=     ',(nOrb(i),i=1,nSym)
            Write(6,'(A,8I8)') 'nBas-nOrb=',((nBas(i)-nOrb(i)),i=1,nSym)
            Write(6,'(A,8I8)') 'nDel=     ',(nDel(i,1),i=1,nSym)
            Call RPA_Warn(4,SecNam//': nDel != nBas-nOrb')
         End If
      End Do

      ! Set default frozen (core) orbitals
      Call Get_iArray('Non valence orbitals',nFro(1,1),nSym)
      Do iSym=1,nSym
         If (nFro(iSym,1).gt.nOcc(iSym,1)) Then
            If (iUHF.eq.1) Then
               Write(6,'(A,8I8)') 'nOcc=',(nOcc(i,1),i=1,nSym)
               Write(6,'(A,8I8)') 'nFro=',(nFro(i,1),i=1,nSym)
               Call RPA_Warn(4,SecNam//': nFro > nOrb')
            Else
               Write(6,'(A,8I8)') 'nOcc(alpha)=',(nOcc(i,1),i=1,nSym)
               Write(6,'(A,8I8)') 'nFro(alpha)=',(nFro(i,1),i=1,nSym)
               Call RPA_Warn(4,SecNam//': nFro > nOrb [alpha]')
            End If
         End If
      End Do
      If (iUHF.eq.2) Then
         Do iSym=1,nSym
            nFro(iSym,2)=nFro(iSym,1)
            If (nFro(iSym,2).gt.nOcc(iSym,2)) Then
               Write(6,'(A,8I8)') 'nOcc(beta)=',(nOcc(i,2),i=1,nSym)
               Write(6,'(A,8I8)') 'nFro(beta)=',(nFro(i,2),i=1,nSym)
               Call RPA_Warn(4,SecNam//': nFro > nOrb [beta]')
            End If
         End Do
      Else
         Do iSym=1,nSym
            nFro(iSym,2)=0
         End Do
      End If

      ! Compute number of virtual orbitals
      Do i=1,iUHF
         Do iSym=1,nSym
            nVir(iSym,i)=nOrb(iSym)-nOcc(iSym,i)
         End Do
      End Do

      ! Set default deleted (virtual) orbitals
      Call iZero(nDel,16)

      ! Register exit
      Call qExit(SecNam)

      End
