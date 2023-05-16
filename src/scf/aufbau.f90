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
! Copyright (C) 1995, Martin Schuetz                                   *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine Aufbau(nAuf,Occup,nOccup,iOK,nD)
!***********************************************************************
!                                                                      *
!     purpose: sets the orbital occupation numbers in the different    *
!              irreps according to the Aufbau scheme...                *
!                                                                      *
!     method:  sets up a map vector, which is sorted with respect to   *
!              the orbital energies, neglecting the boundaries of the  *
!              different irrep blocks. The lowest orbitals then are    *
!              occupied...                                             *
!     input:                                                           *
!       nAuf          : # (doubly) occupied orbitals                   *
!                                                                      *
!     output:                                                          *
!       Occup(nOccup) : orbital occupation numbers                     *
!***********************************************************************
      Use InfSCF, only: nSym, nOcc, TEEE, nFro, nOrb, rTemp
      Use SCF_Arrays, only: EOrb
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero, Half, Two, Three
      Implicit None
#include "Molcas.fh"
!
!     declaration subroutine parameters
      Integer nOccup,nD, iOK
      Real*8 Occup(nOccup,nD)
      Integer nAuf(2)
!
!     declaration of local variables...
      Integer nEOrb, iOrBas, iSym, iOrb, iD, nOrBas, jOrBas, mD, ipOcc, jSym, nElec
      Real*8 UHF_Occ, EFerm, UnlikelyOcc, Fact, Fact2
      Real*8, External:: Fermipop
      Integer, Dimension(:,:), Allocatable:: Map, Irp
      Real*8 Sum_el(2)
      Integer nOrb_AS(2), mOrb_AS(2)
!
!     These occupation number vectors are used to determine if we have
!     convergence.
!
      Integer, Save :: kOccAuf=-1
      Integer, Save :: nOccAuf(MxSym,2,2)
!
      nEOrb = Size(EOrb,1)

      Call mma_allocate(Map,nEOrb,nD,Label='Map')
      Call mma_allocate(Irp,nEOrb,nD,Label='Irp')
!----------------------------------------------------------------------*
! Initialize convergence detection                                     *
!----------------------------------------------------------------------*
      If (kOccAuf.eq.-1) Then
         nOccAuf(:,:,:)=-1
         kOccAuf=1
      End If
!
!---- Set up map...
!
      iOrbAS=1
      Do iSym = 1, nSym
         Do iOrb = 1, nOrb(iSym)-nFro(iSym)
            Do iD = 1, nD
               Irp(iOrbAS,iD)=iSym
               Map(iOrbAS,iD)=iOrbAS
            End Do
            iOrbAS=iOrbAS+1
         End Do
      End Do
      nOrbAS=iOrbAS-1
!
!---- Now sort map with respect to orbital energies (bubblesort)
!
      Do iOrbAS = 1, nOrbAS-1
         Do jOrbAS = nOrbAS-1, iOrbAS, -1
            Do iD = 1, nD
               If (EOrb(Map(  jOrbAS,iD),iD).gt.EOrb(Map(1+jOrbAS,iD),iD))   &
                 Call Swap_Seward(Map(  jOrbAS,iD),Map(1+jOrbAS,iD))
            End Do
         End Do
      End Do
!
!---- and fill up the orbitals...
!
      nOcc(:,:)=0
      call dcopy_(nOccup*nD,[Zero],0,Occup,1)
!
      If (Teee) then
!
         UHF_occ=Three-DBLE(nD)
         mD = 2/nD
         Do iD = 1, nD
            eferm=FermiPop(EOrb(1,iD),Occup(1,iD),nOrbAS,RTemp,nAuf(iD)*mD,UHF_occ)
#ifdef _DEBUGPRINT_
            Write (6,'(A,G20.10)')'         E(Fermi)=',eferm
#endif
         End Do
!
         iOrbAS=0
         Do iSym = 1, nSym
            Do iD = 1, nD
               nOrb_AS(iD)=0
               mOrb_AS(iD)=0
               sum_el(iD)=Zero
            End Do

            jOrbAS = iOrbAS
            unlikelyOcc=0.19D0
            Do iOrb = 1, nOrb(iSym)-nFro(iSym)
               iOrbAS = iOrbAS + 1
               Do iD = 1, nD
                  If (Occup(iOrbAS,iD).ge.UHF_occ-unlikelyOcc) nOrb_AS(iD) = nOrb_AS(iD) + 1
                  If (Occup(iOrbAS,iD).lt.unlikelyOcc) mOrb_AS(iD) = mOrb_AS(iD) + 1
                  sum_el(iD)=sum_el(iD)+Occup(iOrbAS,iD)
               End Do
            End Do
            Fact=nD*half
            Fact2=0.99D0 + DBLE(2-nD)
            Do iD = 1, nD
               nOccAuf(iSym,kOccAuf,iD)=nOrb_AS(iD)
               nOcc(iSym,iD)=INT(Fact*(sum_el(iD)+Fact2/nSym))
            End Do
         End Do
         kOccAuf=3-kOccAuf
!
      Else
!
         Fact=Two/DBLE(nD)
         Do iD = 1, nD
            Do iOrbAS = 1, nAuf(iD)
               iSym=Irp(Map(iOrbAS,iD),iD)
               nOcc(iSym,iD)=nOcc(iSym,iD)+1
               ipOcc=0
               Do jSym = 1, iSym-1
                  ipOcc=ipOcc+nOrb(jSym)
               End Do
               Occup(ipOcc+nOcc(iSym,iD),iD)=Fact
            End Do
         End Do
!
      End If
!
      iOK=1
      Do iD = 1, nD
         nElec=0
         Do iSym=1,nSym
            If (nOccAuf(iSym,1,iD).ne.nOccAuf(iSym,2,iD)) iOK=0
            nElec=nElec+nOccAuf(iSym,1,iD)
         End Do
         If(nElec.ne.nAuf(iD)) iOK=0
      End Do
      If (iOK.eq.1) Then
         Do iSym=1,nSym
            Do iD = 1, nD
               nOcc(iSym,iD)=nOccAuf(iSym,1,iD)
            End Do
         End Do
      End If
!
!---- Write new occupation on the RUNFILE
!
      Call Put_iArray('nIsh',nOcc(1,1),nSym)
      if(nD==2) Call Put_iArray('nIsh beta',nOcc(1,2),nSym)
!
      Call mma_deallocate(Irp)
      Call mma_deallocate(Map)
!
      Return
#ifdef _WARNING_WORKAROUND_
      IF (.False.) Call Unused_real(eferm)
#endif
      End SubRoutine Aufbau
!***********************************************************************
!                                                                      *
! This function computes the Fermi energy level for a number of        *
! energy levels and populates them. Each level is populated with up    *
! to 2 electrons.                                                      *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University                                             *
!          Sweden                                                      *
! Written: June 1999                                                   *
! History:                                                             *
!                                                                      *
!***********************************************************************
      Real*8 Function FermiPop(e,o,n,T,nEle,UHF_occ)
      use Constants, only: Zero, Half, One, Three, Ten
      Implicit None
!----------------------------------------------------------------------*
! Dummy arguments:                                                     *
! e(*) -- Orbital energies, input.                                     *
! o(*) -- Orbital occupations, output.                                 *
! n    -- Number of orbitals, input.                                   *
! T    -- Temperature, input.                                          *
! nEle -- Number of electrons.                                         *
!----------------------------------------------------------------------*
      Integer n,nEle
      Real*8  e(n),o(n),T,UHF_occ
!----------------------------------------------------------------------*
! Local variables:                                                     *
!----------------------------------------------------------------------*
      Real*8  ef,beta,f,f_old,Step,ff
      Real*8  x0,x1,x2,y0,y2,z
      Integer i,iter
#ifdef _DEBUGPRINT_
      Real*8  y1
#endif
!----------------------------------------------------------------------*
! Initialize                                                           *
!----------------------------------------------------------------------*
      ef=Zero
      If (T.le.Zero) Then
         beta=1.0D99
      Else
         beta=One/T
      End If
!----------------------------------------------------------------------*
! Scan for Fermi level                                                 *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
      Write(6,'(a)') 'Scan for Fermi energy level'
      Write(6,'(a)') '       ef             y       '
      Write(6,'(a)') ' -------------- --------------'
#endif

      f=-nEle
      f_old=f
      Do i=1,n
!        Write (6,'(A,G20.10)') 'e(i)=',e(i)
         z=beta*(e(i)-ef)
         z=Min(z,30.d0)
         f=f+UHF_occ/(One+exp(z))
      End Do
      If(f.gt.Zero) Then
         Step=-One
      Else
         Step=One
      End If
      Iter=0
100   Continue
         Iter=Iter+1
         If(Iter.gt.100000) GoTo 101
         f_old=f
         ef=ef+Step
!         f=-nEle
         ff=Zero
!vv overoptimization with Intel compiler
         i=1
300      continue
!         Do i=1,n
            z=beta*(e(i)-ef)
            z=Min(z,30.d0)
            ff=ff+1/(One+exp(z))
            i=i+1
            if(i.le.n) goto 300
!         End Do
         f=-nEle+ff*UHF_occ
#ifdef _DEBUGPRINT_
         Write(6,'(2G20.10)') ef,f
#endif
         If(f*f_old.gt.Zero) GoTo 100
101   Continue
!----------------------------------------------------------------------*
! Refine with interval halving.                                       *
!----------------------------------------------------------------------*
#ifdef _DEBUGPRINT_
      Write(6,'(a)') 'Refine Fermi level with interval halving'
      Write(6,'(a)') '       y0            y2             y1       '
      Write(6,'(a)') ' -------------- -------------- --------------'
      y1=f
#endif
      x0=ef-Step
      x1=ef
      y0=f_old
      x2=half*(x0+x1)
      Iter=0
200   Continue
         Iter=iter+1
         If(Iter.gt.1000) GoTo 201
         ef=x2
         f=-nEle
         Do i=1,n
            z=beta*(e(i)-ef)
            z=Min(z,three*ten)
            f=f+UHF_occ/(One+exp(z))
         End Do
         y2=f
#ifdef _DEBUGPRINT_
         Write(6,'(3f15.8)') y0,y2,y1
#endif
         If(abs(y2).lt.1.0d-9) GoTo 201
         If(y0*y2.le.Zero) Then
            x1=x2
#ifdef _DEBUGPRINT_
            y1=y2
#endif
         Else
            x0=x2
            y0=y2
         End If
         x2=half*(x0+x1)
         GoTo 200
201   Continue

!----------------------------------------------------------------------*
! Populate occupation number vector.                                   *
!----------------------------------------------------------------------*
!     Write (*,*)
      f=Zero
      Do i=1,n
!        Write(*,'(2G20.10)') e(i),ef
         z=beta*(e(i)-ef)
!        Write (*,*) 'z,Beta=',z,Beta
         z=Min(z,Three*Ten)
!        Write (*,*) 'z=',z
         o(i)=UHF_occ/(One+exp(z))
!        Write(*,'(1G20.10)') o(i)
         f=f+o(i)
      End Do
      f=nEle/f
!     Write (*,*)
!     Write(*,'(1G20.10)') f
!     Write (*,*)
      Do i=1,n
         o(i)=f*o(i)
!        Write(*,'(1G20.10)') o(i)
      End Do
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*

      FermiPop=ef
      Return
      End Function FermiPop
