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
      SubRoutine ReadBas(Lhigh,makemean,bonn,breit,
     &                   symmetry,sameorb,AIMP,oneonly,ncont4,
     &                   numballcart,IN,ifinite)
*
*     Suposed to read the maximum of l-values, the number of primitive and
*     contracted functions, the exponents and contraction coefficients
*
      Implicit Real*8 (a-h,o-z)
#include "para.fh"
#include "param.fh"
#include "ired.fh"
#include "Molcas.fh"
#include "stdalloc.fh"
      Integer, Allocatable:: nOff(:,:)
      Character*4 Word
      Character*4 Symmetry
#ifdef _DEBUGPRINT_
      Character*21 chCharge
#endif
      Character*54 Stars
      Logical MakeMean, Bonn, Breit, SameOrb, AIMP, OneOnly, IfTest
      Common /Nucleus/ Charge, Exp_Finite
      Integer OUT, iBeginIRed(8), iDelperSym(8)
      Data IfTest/.False./
*
#ifdef _DEBUGPRINT_
      IfTest=.True.
      chCharge='  Charge of nucleus: '
#endif
      OUT=6
      Stars   ='******************************************************'
      Bonn    =.False.
      Breit   =.False.
      SameOrb =.False.
      AIMP    =.False.
      OneOnly =.False.
      MakeMean=.True.
*
      If (IfTest) Then
         Write(OUT,'(/,/,/,24X,A)') Stars
         Write(OUT,'(24X,2A)') '******** Starting Atomic Spin-Orbit MF',
     &                         ' code ********'
         Write(OUT,'(24X,A,/,/)') Stars
      End If
*
      Do I = 0, lMax
         iCore(I) = 0
      End Do
      Call RdNLst(IN,'AMFI')
 123  Read(IN,'(A4)') Word
      If (IfTest) Write(OUT,'(A4)') Word
      Call UpCase(Word)
      If (WORD.eq.'BONN') Then
         Bonn=.True.
         GoTo 123
      ElseIf (WORD.eq.'BREI') Then
         Breit=.True.
         GoTo 123
      ElseIf (WORD.eq.'FINI') Then
         iFinite=1
         Read(IN,*) Exp_finite
         GoTo 123
      ElseIf (WORD.eq.'SAME') then
         SameOrb=.True.
         GoTo 123
      ElseIf (WORD.eq.'AIMP') then
         AIMP=.True.
                Read(IN,*) lDel,(iCore(I),I=0,ldel)
         If (IfTest) Then
            Write(OUT,*)
            Write(OUT,*) 'CORE to be deleted '
                   Write(OUT,*) '   L   #orbs.  '
            Write(OUT,*)
                   Do I = 0, lDel
                      Write(OUT,'(2I5)') I,iCore(I)
                   End Do
         End If
         GoTo 124
      ElseIf (Word.eq.'ONEO') Then
         OneOnly=.True.
         Write(OUT,*) ' Only one-electron integrals!!'
         Write(OUT,*) ' Probably useful for test-purposes only'
         GoTo 123
      End If
*
 124  Continue
      If (IfTest) Then
         Write(OUT,*) ' AMFI: '
         If (BONN) Then
            Write(OUT,*) ' Bonn-approach for spin-other-orbit part'
         End If
         If (BREIT) Then
            Write(OUT,*) ' Breit-Pauli type of the SO operator'
         Else
            Write(OUT,*) ' Douglas-Kroll type of the SO operator'
         End If
         If (iFinite.eq.0) Then
            Write(OUT,*) ' Point nucleus '
         Else
            Write(OUT,*) ' Finite nucleus'
         End If
      End If
*
      Symmetry='D2H'
      NumbofSym=8
      If (IfTest) Then
         Write(OUT,*) ' Symmetry is D2H'
         If (SameOrb) then
            Write(OUT,*) ' Same-Orbit only'
         Else
            Write(OUT,*) ' Other-Orbit included'
         End If
      End If
      Read(IN,*) Charge,Lhigh
      If (Lhigh.gt.Lmax) Then
          Write(OUT,*) ' Sorry, so far this code deals only ',
     &                 '    with maximum l-values of ',Lmax
          Call Abend()
      End If
#ifdef _DEBUGPRINT_
      Write(OUT,'(A21,F5.2)') chCharge, Charge
#endif
      Call InitiRed
      Do iredrun=1,numbofsym
         Do Lrun=0,Lhigh
            nmbMperIRL(iredrun,Lrun)=0
         End Do
      End Do
      If (IfTest) Write(OUT,'(/,A)') '  Used SOC basis set: '
      Do Lrun=0,Lhigh
         Read(IN,*) nprimit(Lrun),ncontrac(Lrun)
         If (IfTest) Then
            Write(OUT,'(/,A,I2,A,I2)') '  nExp: ', nprimit(Lrun),
     &                                 ' lAng: ',lRun
            Write(OUT,'(I3,I3)') nprimit(Lrun),ncontrac(Lrun)
         End If
         If (nprimit(Lrun).gt.MxprimL) Then
            Write(OUT,*) 'To many primitives for L=',Lrun,
     &        ' increase MxprimL in para.fh or reduce ',
     &        ' the number of primitives to at least ',MxprimL
            Call Abend()
         End If
         If (ncontrac(Lrun).gt.MxcontL) Then
            Write(OUT,*) ' To many contracted fncts for L=',Lrun,
     &                   ' increase MxcontL in para.fh or ',
     &                   ' reduce the number of contracted functions',
     &                   ' to at most ',MxcontL
            Call Abend()
            End If
         If (ncontrac(Lrun).gt.nprimit(Lrun)) Then
            Write(OUT,*) ' You have more contracted than ',
     &                   ' uncontracted functions, I don''t believe ',
     &                   ' that. Sorry!! '
            Call Abend()
      End If
*
*     Read input in MOLCAS-style
*
              Read(IN,*) (exponents(ILINE,Lrun),ILINE=1,nprimit(Lrun))
              Do ILINE=1,nprimit(Lrun)
                 Read(IN,*) (cntscrtch(ILINE,JRUN,Lrun),Jrun=1,
     &                             ncontrac(Lrun))
              End Do
*
*     End of reading for the current L-value
*
      If (IfTest) Then
         Write(OUT,'(5E18.8)') (exponents(ILINE,Lrun),
     &                          ILINE=1,nprimit(Lrun))
         Do Irun = 1, ncontrac(Lrun)
            Write(OUT,*) ' orbital : ',irun
            Write(OUT,'(6(1X,F12.7))')
     &           (cntscrtch(I,Irun,Lrun),I=1,nprimit(Lrun))
         End Do
      End If
*
*     Setting the numbers of cartesians per IR
*
         Do iRedRun = 1, NumbofSym
            nFunctions(iRedRun,Lrun)=0
         End Do
         Do mRun=-Lrun,Lrun
            nfunctions(ipow2ired(ipowxyz(1,mrun,Lrun),
     &      ipowxyz(2,mrun,Lrun),Ipowxyz(3,mrun,Lrun)),Lrun)=
     &      nfunctions(ipow2ired(ipowxyz(1,mrun,Lrun),
     &      ipowxyz(2,mrun,Lrun),ipowxyz(3,mrun,Lrun)),Lrun)+
     &      ncontrac(Lrun)
         End Do
         Do mRun=-Lrun,Lrun
            nmbMperIRL(ipow2ired(ipowxyz(1,mrun,Lrun),
     &      ipowxyz(2,mrun,Lrun),Ipowxyz(3,mrun,Lrun)),lruN)=
     &      nmbMperIRL(ipOw2ired(ipowxyz(1,mrun,Lrun),
     &      ipowxyz(2,mrun,Lrun),IpowxYz(3,mrun,Lrun)),lruN)+1
         End Do
         If (IfTest) Then
            Write(OUT,'(A,8I4)')
     &      ' Number of functions per IR: ',(nfunctions(iredrun,Lrun),
     &                                           iredrun=1,numbofsym)
         End If
      End Do   ! End Do for loop over L-values
*
      If (IfTest) Then
          Write(OUT,*) ' Distribution of M-values'
          Do Lrun=0,Lhigh
             Write(OUT,*) (nmbMperIRL(nsym,Lrun),nsym=1,numbofsym)
          End Do
      End If
*
      numbofcart=0
      Do lrun=0,Lhigh
         numbofcart=numbofcart+(Lrun+Lrun+1)*
     &              ncontrac(Lrun)
      End Do
*
      Call mma_allocate(nOff,numbofcart,2,Label='nOff')
*
      Do iredrun=1,numbofsym
         nfunctperIRED(iredrun)=0
      End Do
      Do Lrun=0,Lhigh
         Do iredrun=1,numbofsym
         nfunctperIRED(iredrun)=nfunctperIRED(iredrun)+
     &                          nfunctions(iredrun,Lrun)
         End Do
      End Do
      If (IfTest) Then
         Write(OUT,'(A,8I3)')
     &        ' Total number of atomic functions per IRED ',
     &        (nfunctperIRED(iredrun),iredrun=1,numbofsym)
      End If
      isum=0
      Do iredrun=1,numbofsym
         itotalperIR(iredrun)=nfunctperIRED(iredrun)
         isum=isum+itotalperIR(iredrun)
      End Do
      numballcart=isum
      iorbrun=0
      Do iredrun=1,numbofsym
         Do inired=1,itotalperIR(iredrun)
            iorbrun=iorbrun+1
            IREDoffunctnew(Iorbrun)=iredrun
         End Do
      End Do
      If (IfTest) Then
         Write(OUT,'(A,8I3)')
     &        'including additional functions per IRED ',
     &        (itotalperIR(iredrun),iredrun=1,numbofsym)
      End If
      Do iredrun=1,numbofsym
         ibeginIRED(iredrun)=0
      End Do
      Do lrun=0,Lhigh
         Do mrun=-lrun,lrun
            iredLM(mrun,lrun)=ipow2ired(ipowxyz(1,mrun,Lrun),
     &      ipowxyz(2,mrun,Lrun),
     &      ipowxyz(3,mrun,Lrun))
            incrLM(mrun,lrun)=ibeginIRED(iredLM(mrun,lrun))
            ibeginIRED(iredLM(mrun,lrUn))=
     &      ibeginIRED(iredLM(mrun,lrun))+ncontrac(lrun)
         EndDo
      EndDo
      If (IfTest) Then
         Do lrun=0,Lhigh
            Write(OUT,'(A,I4,A,21I3)') 'L= ',lrun,
     &                                 ' shifts inside the IRED',
     &                       (incrLM(mrun,lrun),mrun=-lrun,lrun)
         End Do
      End If
      shiftIRED(1)=0
      Do iredrun=2,numbofsym
         shiftIRED(iredrun)=shiftIRED(iredrun-1)
     &                   +itotalperIR(iredrun-1)
      End Do
      If (IfTest) Then
         Write(OUT,'(A,8I4)') 'shifts for the IREDs ',
     &       (shiftIRED(iredrun),iredrun=1,numbofsym)
         Do lrun=0,Lhigh
            Do mrun=-Lrun,Lrun
               Do irun=1,ncontrac(lrun)
                  Write(OUT,*) 'L,M,contr funct, absolute number ',
     &            lrun,mrun,irun,shiftired(iredLM(mrun,lrun))+
     &            incrLM(mrun,Lrun)+irun
               End Do
            End Do
         End Do
      End If
      shiftIRIR(1)=0
      irun=1
      Do ired1=2,numbofsym
         Do ired2=1,ired1
            irun=irun+1
            If (ired2.eq.1) Then
               shiftIRIR(irun)=shiftIRIR(irun-1)+
     &         (itotalperIR(ired1-1)*itotalperIR(ired1-1)+
     &         itotalperIR(ired1-1))/2
            Else
               shiftIRIR(irun)=shiftIRIR(irun-1)+
     &         itotalperIR(ired1)*itotalperIR(ired2-1)
            End If
         End Do
      End Do
      Do lrun=0,Lhigh
         Do Mrun=-Lrun,Lrun
            ired=iredLM(Mrun,Lrun)
            ishifter=shiftIRED(ired)+incrLM(mrun,lrun)
            Do icart=1,ncontrac(Lrun)
               moffunction(ishifter+icart)=Mrun
               Loffunction(ishifter+icart)=Lrun
               IREDoffunction(ishifter+Icart)=ired
               nOff(ishifter+Icart,2)=icart
            End Do
         End Do
      End Do
      Do irun = 1, numbofcart
         nOff(irun,1)=irun
      End Do
      Do nsymrun=1,numbofsym
         idelpersym(nsymrun)=0
      End Do
      Do nsymrun=1,numbofsym
         nrtofiperIR(nsymrun)=itotalperIR(nsymrun)
      End Do
      If (AIMP) Then
*
*     Generate list of orbitals to be removed
*
         If (IfTest)
     &      Write(OUT,'(/,A)') '  Core removed for use with AIMP'
         ikeeporb=0
         numbprev=0
         Do irun=1,numbofcart
4712        If (irun.eq.1.or.(irun.ge.2.and.noff(irun,1).eq.
     &         numbprev+1)) Then
               Lval=Loffunction(irun)
               number=nOff(irun,1)
               itype=nOff(irun,2)
               If (itype.le.icore(lval)) then
                  Write(OUT,777) number,itype,lval
                  idelpersym(IREDoffunction(irun))=
     &            idelpersym(IREDoffunction(irun))+1
                  numbprev=number
               Else
                  ikeeporb=ikeeporb+1
                  ikeeplist(ikeeporb)=number
                  numbprev=number
               End If
            Else
               ikeeporb=ikeeporb+1
               ikeeplist(ikeeporb)=numbprev+1
               numbprev=numbprev+1
               GoTo 4712
            End If
         End Do
         ikeeporb=0
         Do nsymrun=1,numbofsym
            nrtofiperIR(nsymrun)=
     &      itotalperIR(nsymrun)-idelpersym(nsymrun)
         End Do
         Do nsymrun=1,numbofsym
            ikeeporb=ikeeporb+nrtofiperIR(nsymrun)
         End Do
         If (IfTest) Then
            Write(OUT,'(A,8I3)')
     &      '  Number of funct. per IRED after removing core: ',
     &           (nrtofiperIR(iredrun),iredrun=1,numbofsym)
            Write(OUT,'(I4,A)') ikeeporb,
     &      ' orbitals left after deleting core'
         End If
      End If
      nmax=max(6,ncontrac(0))
      Do lrun=1,Lhigh
         nmax=max(nmax,ncontrac(lrun))
      End Do
      ncont4=nmax*nmax*nmax*nmax
*
      Call mma_deallocate(nOff)
*
      Return
777   Format('  Orbital number ',I4,' is the ',I3,'th of L-value ',I3,
     &       ' it will be removed !!!')
      End
