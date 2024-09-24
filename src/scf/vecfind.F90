!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) Per-Olof Widmark                                       *
!               2017, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      Subroutine VecFind(OccSet,FermSet,CharSet,SpinSet)
!***********************************************************************
!                                                                      *
! This routine figure out which set of starting orbitals are used.     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden                                     *
!                                                                      *
!          R. Lindh                                                    *
!          Uppsala University, Sweden, Feb 2017                        *
!          Remove Work                                                 *
!                                                                      *
!***********************************************************************
!                                                                      *
! Input vector mappings                                                *
!                                                                      *
! LstVec              InVec                                            *
! -1 Die               0 Core                                          *
!  0 Old SCF           1 NDDO (not used)                               *
!  1 Guessorb          2 Lumorb                                        *
!  2 Lumorb            3 Old density                                   *
!  3 Old Density       4 Restart (not used)                            *
!  4 Core              8 Old SCF                                       *
!  5 NDDO (use not)    9 Guessorb                                      *
!                                                                      *
!***********************************************************************
#ifdef _HDF5_
      Use mh5, Only: mh5_fetch_attr
      use InfSCF, only: FileOrb_ID
#endif
      use InfSCF, only: iAu_ab, InVec, isHDF5, nD, nSym, nStOpt, SCF_FileOrb, Tot_Charge, Tot_El_Charge, &
                        Tot_Nuc_Charge, nBas, LstVec, nOcc, nAufb
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Half
      Implicit None
!----------------------------------------------------------------------*
! Dummy arguments                                                      *
!----------------------------------------------------------------------*
      Logical OccSet
      Logical FermSet
      Logical CharSet
      Logical SpinSet
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
      Character*80 cList
      Logical Found
      Integer mynSym
      Integer mynBas(8)
      Integer mynOrb(8)
      Character*10 infoLbl
      Real*8, Dimension(:), Allocatable:: EOrb
      Integer nSQRSum, iSym, i, nData, iVer, j, N2, N1, iDSpin, nEle, iTmp, nEle1, nEle2, mTmp, iOff, n, iBas, iRC
      Real*8 GAP, eAlpha, eBeta, tmp
!----------------------------------------------------------------------*
! Setup                                                                *
!----------------------------------------------------------------------*
      nSqrSum=0
      Do iSym = 1, nSym
         nSqrSum=nSqrSum+nBas(iSym)*nBas(iSym)
      End Do
!----------------------------------------------------------------------*
! Check start orbital priority list                                    *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
      Write(6,*) 'VecFind: LstVec',LstVec
#endif
      Do i=1,nStOpt
         Found=.false.
         If(LstVec(i).eq.0) Then
#ifdef _DEBUGPRINT_
            Write(6,'(2i5,a)') i,LstVec(i),' Old scf orbitals'
#endif
            Call qpg_darray('SCF orbitals',Found,nData)
            If(Found .and. nData.eq.nSqrSum) Then
#ifdef _DEBUGPRINT_
               Write(6,*) 'Found orbitals'
#endif
               Call qpg_darray('OrbE',Found,nData)
               If(Found) Then
#ifdef _DEBUGPRINT_
                  Write(6,*) 'Found energies'
#endif
                  InVec=8
                  GoTo 101
               End If
            End If
         Else If(LstVec(i).eq.1) Then
#ifdef _DEBUGPRINT_
            Write(6,'(2i5,a)') i,LstVec(i),' Guessorb orbitals'
#endif
            Call qpg_darray('Guessorb',Found,nData)
            If(Found .and. nData.eq.nSqrSum) Then
               Call qpg_darray('Guessorb energies',Found,nData)
               If(Found) Then
                  InVec=9
                  GoTo 101
               End If
            End If
         Else If(LstVec(i).eq.2) Then
#ifdef _DEBUGPRINT_
            Write(6,'(2i5,a)') i,LstVec(i),' Lumorb orbitals'
#endif
            Call F_Inquire(SCF_FileOrb,Found)
            If(Found) Then
               If (isHDF5) Then
#ifdef _HDF5_
                 iVer=-1
                 Call mh5_fetch_attr(fileorb_id,'NSYM',mynSym)
                 Call mh5_fetch_attr(fileorb_id,'NBAS',mynBas)
#endif
               Else
                 Call ChkVec(SCF_FileOrb,iVer,mynSym,mynBas,mynOrb,InfoLbl,iRc)
               End If
               If(iVer.eq.0) Found=.false.
               If(mynSym.ne.nSym) Found=.false.
               Do j=1,nSym
                 If(mynBas(j).ne.nBas(j)) Found=.false.
               End Do
            End If
            If(Found) Then
#ifdef _DEBUGPRINT_
               Write(6,*) 'Found INPORB'
#endif
               InVec=2
               GoTo 101
            End If
         Else If(LstVec(i).eq.3) Then
#ifdef _DEBUGPRINT_
            Write(6,'(2i5,a)') i,LstVec(i),' Density'
#endif
            InVec=0
         Else If(LstVec(i).eq.4) Then
#ifdef _DEBUGPRINT_
            Write(6,'(2i5,a)') i,LstVec(i),' Core orbitals'
#endif
            InVec=0
            Found=.true.
            GoTo 101
         Else
#ifdef _DEBUGPRINT_
            Write(6,'(2i5,a)') i,LstVec(i),' Die'
#endif
            InVec=0
            Found=.false.
            GoTo 101
         End If
#ifdef _DEBUGPRINT_
         Write(6,*) 'LstVec(i),Found:',LstVec(i),Found
#endif
      End Do
      Found=.false.
101   Continue
#ifdef _DEBUGPRINT_
      Write(6,*) 'VecFind: InVec, Found=',InVec, Found
#endif
!----------------------------------------------------------------------*
! Did we find the requested orbitals?                                  *
!----------------------------------------------------------------------*
      If(.not.Found) Then
         cList=''
         n2=0
         Do i=1,nStOpt
            If(LstVec(i).eq.0) Then
               n1=n2+1
               n2=n1+9
               Write(cList(n1:n2),'(a)') 'Old Scf, '
            Else If (LstVec(i).eq.1) Then
               n1=n2+1
               n2=n1+10
               Write(cList(n1:n2),'(a)') 'Guessorb, '
            Else If (LstVec(i).eq.2) Then
               n1=n2+1
               n2=n1+8
               If (isHDF5) Then
                  Write(cList(n1:n2),'(a)') 'HDF5, '
               Else
                  Write(cList(n1:n2),'(a)') 'Lumorb, '
               End If
            Else If (LstVec(i).eq.3) Then
               n1=n2+1
               n2=n1+13
               Write(cList(n1:n2),'(a)') 'Old density, '
            Else If (LstVec(i).eq.4) Then
               n1=n2+1
               n2=n1+6
               Write(cList(n1:n2),'(a)') 'Core, '
            Else If (LstVec(i).eq.5) Then
               n1=n2+1
               n2=n1+6
               Write(cList(n1:n2),'(a)') 'NDDO, '
            Else If (LstVec(i).eq.-1) Then
               GoTo 290
            Else
               n1=n2+1
               n2=n1+9
               Write(cList(n1:n2),'(a)') 'Unknown, '
            End If
         End Do
290      Continue
         Call SysAbendMsg('SCF:','Cannot find start orbitals according to list:',cList)
      End If
!----------------------------------------------------------------------*
! What are the defaults for the different cases?                       *
!----------------------------------------------------------------------*
      Call Get_dScalar('Total Nuclear Charge',Tot_Nuc_Charge)
      Tot_El_Charge=Tot_Charge-Tot_Nuc_Charge
      If(InVec.eq.0) Then
!
! We will use core diagonalization
!
#ifdef _DEBUGPRINT_
         Write(6,*) 'Using core diagonalization'
#endif
         OccSet=.false.
         FermSet=.true.
      Else If(Invec.eq.1) Then
!
! We will use NDDO orbitals, should not be used!
!
#ifdef _DEBUGPRINT_
         Write(6,*) 'Using NDDO orbitals'
#endif
      Else If(Invec.eq.2) Then
!
! We will use Lumorb orbitals
!
#ifdef _DEBUGPRINT_
         Write(6,*) 'Using Lumorb orbitals'
#endif
         Call ChkLumo(OccSet,FermSet,SpinSet)
      Else If(Invec.eq.3) Then
!
! We will use density as start, does it even work?
!
#ifdef _DEBUGPRINT_
         Write(6,*) 'Using density'
#endif
      Else If(Invec.eq.4) Then
!
! This is a restart case
!
#ifdef _DEBUGPRINT_
         Write(6,*) 'Using Restart'
#endif
      Else If(Invec.eq.8) Then
!
! We will use old SCF orbitals
!
#ifdef _DEBUGPRINT_
         Write(6,*) 'Using SCF orbitals'
         Write(6,*) 'tot_charge',tot_charge
         Write(6,*) 'tot_el_charge',tot_el_charge
         Write(6,*) 'tot_nuc_charge',tot_nuc_charge
#endif
         Call qpg_iarray('SCF nOcc',Found,nData)
         idspin=0
         If(Found) Then
#ifdef _DEBUGPRINT_
            Write(6,*) 'vecfind: Alright, old scf orbitals it is'
#endif
            nEle=0
            If(nD==1) Then
               Call Get_iarray('SCF nOcc',nOcc(1,1),nSym)
               Call Get_iScalar('SCF mode',iTmp)
               If(iTmp.eq.0) Then
#ifdef _DEBUGPRINT_
                  Write(6,*) 'Starting RHF with RHF orbitals'
#endif
                  Do iSym=1,nSym
                     nEle=nEle+2*nOcc(iSym,1)
                  End Do
               Else
#ifdef _DEBUGPRINT_
                  Write(6,*) 'Starting RHF with UHF orbitals'
#endif
                  Call Get_iarray('SCF nOcc_ab',nOcc(1,2),nSym)
                  nEle1=0
                  nEle2=0
                  Do iSym=1,nSym
                     nEle1=nEle1+2*nOcc(iSym,1)
                     nEle2=nEle2+nOcc(iSym,1)+nOcc(iSym,2)
                     idspin=idspin+nOcc(iSym,1)-nOcc(iSym,2)
                  End Do
                  If(nEle1.ne.nEle2) Then
!                    Tot_Charge=Tot_Nuc_Charge-nEle1
!                    Tot_El_Charge=-nEle1
                     CharSet=.true.
                  End If
                  nEle=nEle2
#ifdef _DEBUGPRINT_
                  Write(6,*) 'After strange code'
                  Write(6,*) 'tot_charge',tot_charge
                  Write(6,*) 'tot_el_charge',tot_el_charge
                  Write(6,*) 'tot_nuc_charge',tot_nuc_charge
#endif
               End If
            Else
#ifdef _DEBUGPRINT_
               Write(6,*) 'Starting UHF with RHF/UHF orbitals'
#endif
               Call Get_iarray('SCF nOcc',nOcc(1,1),nSym)
               Call qpg_iarray('SCF nOcc_ab',Found,nData)
               If(Found) Then
                  Call Get_iarray('SCF nOcc_ab',nOcc(1,2),nSym)
               Else
                  Call Get_iarray('SCF nOcc',nOcc(1,2),nSym)
               End If
               Do iSym=1,nSym
                  nEle=nEle+nOcc(iSym,1)+nOcc(iSym,2)
                  idspin=idspin+nOcc(iSym,1)-nOcc(iSym,2)
               End Do
            End If
#ifdef _DEBUGPRINT_
            Write(6,*) 'idspin',idspin
            Write(6,*) 'iAu_ab',iAu_ab
#endif
            idspin=idspin-iAu_ab
            If(Abs(Tot_El_Charge+nEle).gt.Half.or.idspin.ne.0) Then
               If(Abs(Tot_El_Charge+nEle).gt.Half) Then
#ifdef _DEBUGPRINT_
                  Write(6,*) 'System have changed charge!'
#endif
               End If
               If(idspin.ne.0) Then
#ifdef _DEBUGPRINT_
                  Write(6,*) 'System have changed spin!'
#endif
               End If
               If(CharSet.or.idspin.ne.0) Then
                  OccSet=.false.
                  FermSet=.true.
               Else
                  OccSet=.true.
                  FermSet=.false.
!                 Tot_Charge=Tot_Nuc_Charge-nEle
!                 Tot_El_Charge=-nEle
               End If
            Else
#ifdef _DEBUGPRINT_
               Write(6,*) 'System have same spin and charge'
#endif
               OccSet=.true.
               FermSet=.false.
            End If
#ifdef _DEBUGPRINT_
            Write(6,*) 'OccSet  ',OccSet
            Write(6,*) 'FermSet ',FermSet
            Write(6,*) 'CharSet ',CharSet
            Write(6,*) 'SpinSet ',SpinSet
            Write(6,*) 'nOcc',nOcc
#endif
         Else
            OccSet=.false.
            FermSet=.true.
         End If
      Else If(Invec.eq.9) Then
!
! We will use Guessorb orbitals
!
#ifdef _DEBUGPRINT_
         Write(6,*) 'Using Guessorb orbitals'
#endif
         If(OccSet) Then
#ifdef _DEBUGPRINT_
            Write(6,*) 'Occupation is set'
#endif
            !continue
         Else If(FermSet) Then
#ifdef _DEBUGPRINT_
            Write(6,*) 'Fermi is set'
#endif
            !continue
         Else
#ifdef _DEBUGPRINT_
            Write(6,*) 'Must decide if to use Fermi'
#endif
            If(nAufb(1).eq.-1) Then
               mtmp=Int(-Tot_El_Charge+0.1D0)
               If(nD==1) Then
                  If(Mod(mtmp,2).ne.0) Then
                     Write(6,*) 'VecFind: Error in number of electrons'
                     Write(6,*) '         An even number of electrons ','are required by RHF, use UHF'
                     Write(6,*)
                     Call Abend()
                  End If
                  nAufb(1)=mtmp/2
               Else
                  nAufb(2)=(mtmp-iAu_ab)/2
                  nAufb(1)=Int(-Tot_El_Charge-nAufb(2))
               End If
            End If
#ifdef _DEBUGPRINT_
            Write(6,*) 'nAufb',nAufb
            Write(6,*) 'Now figure out homo-lumo gap'
#endif
            Call qpg_darray('Guessorb energies',Found,nData)
            Call mma_allocate(EOrb,nData,Label='EOrb')
            Call get_darray('Guessorb energies',Eorb,nData)
            If(nD==1) Then
               Call GetGap(Eorb,nData,nAufb(1),Gap,Ealpha)
            Else
               Call GetGap(Eorb,nData,nAufb(1),tmp,Ealpha)
               Call GetGap(Eorb,nData,nAufb(2),Gap,Ebeta)
               Gap=Min(tmp,Gap)
            End If
            Call get_darray('Guessorb energies',Eorb,nData)
            If(Gap.ge.Half) Then
               If(nD==1) Then
                  iOff=0
                  Do iSym=1,nSym
                     n=0
                     Do iBas=1,nBas(iSym)
                        If(EOrb(iOff+iBas).lt.Ealpha) n=n+1
                     End Do
                     nOcc(iSym,1)=n
                     iOff=iOff+nBas(iSym)
                  End Do
               Else
                  iOff=0
                  Do iSym=1,nSym
                     n=0
                     Do iBas=1,nBas(iSym)
                        If(EOrb(iOff+iBas).lt.Ealpha) n=n+1
                     End Do
                     nOcc(iSym,1)=n
                     iOff=iOff+nBas(iSym)
                  End Do
                  iOff=0
                  Do iSym=1,nSym
                     n=0
                     Do iBas=1,nBas(iSym)
                        If(EOrb(iOff+iBas).lt.Ebeta) n=n+1
                     End Do
                     nOcc(iSym,2)=n
                     iOff=iOff+nBas(iSym)
                  End Do
               End If
               OccSet=.true.
               FermSet=.false.
#ifdef _DEBUGPRINT_
               Write(6,*) 'Decided on occupation list'
#endif
            Else
               OccSet=.false.
               FermSet=.true.
#ifdef _DEBUGPRINT_
               Write(6,*) 'Decided on Fermi aufbau'
#endif
            End If
            Call mma_deallocate(EOrb)
#ifdef _DEBUGPRINT_
            Write(6,*) 'Gap is',Gap
#endif
         End If
      Else
!
! This case should not appear
!
         Call SysAbendMsg('SCF:','Internal error in VecFind!','InVec have illegal value')
      End If
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
      Write(6,'(a,i2)') 'VecFind: InVec=',InVec
      Write(6,*) 'OccSet=',OccSet
      Write(6,*) 'FermSet=',FermSet
#endif
!----------------------------------------------------------------------*
!                                                                      *
!----------------------------------------------------------------------*
      Return
      End Subroutine VecFind
