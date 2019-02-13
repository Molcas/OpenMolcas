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
      SUBROUTINE PRIMO(Header,PrOcc,PrEne,ThrOcc,ThrEne,
     *                 nSym,nBas,nOrb,Name,Ene,Occ,CMO,iPrForm)
*
************************************************************************
*                                                                      *
*     Purpose: print a set of orbitals                                 *
*                                                                      *
*     Calling parameters:                                              *
*     Header : Header line which is printed together with the          *
*              section title                                           *
*     nSym   : number of irreducible representations                   *
*     nOrb   : Total number of orbitals per irred. rep.                *
*     nBas   : Number of basis functions per irred. rep.               *
*     Name   : Center and function type label per basis function       *
*     Ene    : Orbital Energies                                        *
*     CMO    : Orbital coefficients                                    *
*     Occ    : Orbital Occupations                                     *
*     PrEne  : Switch to turn printing of orbital Energies ON/OFF      *
*     PrOcc  : Switch to turn printing of orbital Occupatios ON/OFF    *
*     ThrOcc : Threshold for Occupation number to be printed           *
*              Orbitals with Occupation number less than ThrOcc are    *
*              not printed                                             *
*     ThrEne : Threshold for orbitals Energies to be printed           *
*              Orbitals with Energies larger than ThrEne are           *
*              not printed                                             *
*     iPrForm: print level: -1 - Old way, 0- None, 1-list, 2-Short,    *
*              3-long                                                  *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "real.fh"
*
#include "Molcas.fh"
      DIMENSION nBas(*),nOrb(*),Ene(*),CMO(*),Occ(*)
      CHARACTER*(LENIN8) Name(*)
      CHARACTER*(*) Header
      CHARACTER*24 FMT0,FMT1,FMT2,LABEL0,LABEL1,LABEL2
      Character*3 IrrepName(8)
      Character*20 Fmt_s, Fmt_l, Fmt_f, Fmt
      Character*180 Line
      Parameter (Magic=5+1+LENIN8+1+6+3)
      Character*(Magic) ChCMO
      Character*4 Star(10)
      Real*8 Shift(10)
      LOGICAL PrEne,PrOcc, Large, Go, Debug, Header_Done
      Logical Reduce_Prt
      External Reduce_Prt
      Character*(LENIN8) Clean_BName !14
      External Clean_BName
      Debug=.false.
#ifdef _DEBUG_
      Debug=.true.
#endif
*
*----------------------------------------------------------------------*
*                                                                      *
      Do i=1,10
         Star(i)=' '
      End Do
      iPL=iPrintLevel(-1)
      If (Reduce_Prt().and.iPL.lt.3) iPL=0
      If (iPL.le.1) Return
*
      If (iPrForm.eq.0) Return
*
*----------------------------------------------------------------------*
*     Print Header                                                     *
*----------------------------------------------------------------------*
*
*
      Write(6,*)
      Call CollapseOutput(1,'   Molecular orbitals:')
      Write(6,'(6X,A)') '-------------------'
      Write(6,*)
      Write(6,'(6X,2A)') 'Title: ',Header(:mylen(Header))
c     Write(6,*)
c     Write(6,*) 'test print out'
*----------------------------------------------------------------------*
*     Legend (general)                                                 *
*----------------------------------------------------------------------*
      If(iPrForm.ge.4) Then
         Write (6,'(6X,A)') 'LEGEND'
         Write (6,'(6X,A)') '============================'
     &                    //'============================='
     &                    //'============================'
     &                    //'============================'
         Write (6,'(6X,A)') 'A basis set label has the st'
     &                    //'ructure "nl+m" or "nl-m" wher'
     &                    //'e n is the principle quantum'
     &                    //'number for atomic basis'
         Write (6,'(6X,A)') 'functions (a star, "*", deno'
     &                    //'tes a primitive, diffuse or a'
     &                    //' polarization function), l d'
     &                    //'enotes the angular quantum, '
         Write (6,'(6X,A)') 'and m the magnetic quantum n'
     &                    //'umber. The basis functions ar'
     &                    //'e normally real Spherical ha'
     &                    //'rmonics. For Cartesian type '
         Write (6,'(6X,A)') 'basis functions we use the n'
     &                    //'otation "ijk" to denote the '
     &                    //'power of the x, y, and z com'
     &                    //'ponent in the radial part.'
         Write (6,'(6X,A)') 'For p-functions we always use'
     &                    //' the notation px, py, and pz.'
      End If
*----------------------------------------------------------------------*
*     Define standard output format                                    *
*----------------------------------------------------------------------*
      NCOLS=10
      LABEL0='Orbital '
      LABEL1='Energy  '
      LABEL2='Occ. No.'
      FMT0='(10X,A12,3X,10(I5,A,1X))'
      FMT1='(10X,A12,2X,10F10.4)'
      FMT2='(5X,I4,1X,A,10F10.4)'
*
*----------------------------------------------------------------------*
*     Set up list of pointers                                          *
*----------------------------------------------------------------------*
*
      ISX=1
      ISB=1
      ISCMO=1
*
*----------------------------------------------------------------------*
*     Start loop over all irreducible representations                  *
*----------------------------------------------------------------------*
*
      nTot=0
      Do iSym = 1, nSym
         nTot=nTot+nBas(iSym)
      End Do
      Large=.false.
      if (iPrForm.eq.-1) Large=nTot.gt.256
      if (iPrForm.eq.1)  Large=.false.
      if (iPrForm.eq.2)  Large=.true.
      if (iPrForm.eq.3)  Large=.false.
      if (iPrForm.eq.4)  Large=.false.
      Call Get_cArray('Irreps',IrrepName,24)
*
      Do 100 iSym=1,nSym
        nB=nBas(iSym)
        nO=nOrb(iSym)
        If( nO.EQ.0 ) Go To 100
*
        Header_Done=.False.
*
*----------------------------------------------------------------------*
*     Start loop over columns                                          *
*----------------------------------------------------------------------*
*
        If (Large) Then
           FMT_s='(I5,1X,A,A,F6.3,A)'
           FMT_l='(I5,1X,A,A,F6.2,A)'
           FMT_f='(I5,1X,A,A,F6.1,A)'
           ThrCff=0.10d00
*
           Do iSO = 0, nO-1
            If (PrEne.and.PrOcc) Then
               Go = Ene(ISX+IsO).LE.ThrEne .AND. Occ(ISX+IsO).GE.ThrOcc
            Else If (PrEne) Then
               Go = Ene(ISX+IsO).LE.ThrEne
            Else If (PrOcc) Then
               Go = Occ(ISX+IsO).GE.ThrOcc
            Else
               Go = .False.
            End If
            If (Debug) Write (6,*) 'Go=',Go
            If (Go) Then
               If (.Not.Header_Done) Then
                  Write(6,'(/6X,A,I2,A,A)')
     &                 'Molecular orbitals for symmetry species',iSym,
     &                 ': ',IrrepName(iSym)
                  Write (6,*)
*
*----------------------------------------------------------------------*
*                 Start by listing the basis of this irrep             *
*----------------------------------------------------------------------*
*
                  Write (6,*)
                  Write (6,'(6X,A)') 'Basis set list:'
                  nCol=Min((nB+9)/10,5)
                  Inc=(nB+nCol-1)/nCol
                  jSB=iSB-1
                  Do iB = 1, Inc
ct.t.;
c                    Write (6,'(6X,5(I3,1X,A,15X))')
                     Write (6,'(4X,5(I5,1X,A,9X))')
ct.t.; end
     &               (jB,Clean_BName(Name(jSB+jB),LENIN),jB=iB, nB, Inc)
                  End Do
                  Write (6,*)
*
                  Write (6,*)
                  Write (6,*)
     &                 ' Index Energy  Occupation Coefficients ...'
                  Write (6,*)
                  Header_Done=.True.
               End If
*............. This will occupy the first 25 positions
               If (PrEne) Then
                  Energy=Ene(ISX+iSO)
               Else
                  Energy=Zero
               End If
               If (PrOcc) Then
                  Occupation=Occ(ISX+iSO)
               Else
                  Occupation=Zero
               End If
               Write(Line,'(I5,2F10.4)') iSO+1,Energy,Occupation
               If (Debug) Write (6,'(A,A)') 'Line=',Line
               iPos=1
               Cff_Mx=0.0D0
               Do iB = 0, NB-1
                  If (Abs(CMO(isCMO+iSO*nB+iB)).gt.Cff_Mx)
     &               Cff_Mx=Abs(CMO(isCMO+iSO*nB+iB))
               End Do
               If (Debug) Write (6,*) 'Cff_Mx=',Cff_Mx
               Do iB = 0, NB-1
                  If (Debug) Write (6,*) ' iB=',iB
                  If (Abs(CMO(isCMO+iSO*nB+iB)).ge.Cff_Mx*0.5D0) Then
                     If (Debug) Write (6,*) ' Process iB=',iB
                     If (Abs(CMO(isCMO+iSO*nB+iB)).ge.100.0D0) Then
                        Fmt=Fmt_f
                     Else If (Abs(CMO(isCMO+iSO*nB+iB)).ge.10.0D0) Then
                        Fmt=Fmt_l
                     Else
                        Fmt=Fmt_s
                     End If
                     If (Debug) Write (6,*) ' Fmt=',Fmt
                     Write (ChCMO,Fmt)
     &                  iB+1,Clean_BName(Name(ISB+IB),LENIN),
     &                  '(',CMO(isCMO+iSO*nB+iB),'), '
                     If (Debug) Write (6,'(A)') ChCMO
                     If (iPos.eq.1) Then
                        Line(2+Magic:2+Magic*2-1)=ChCMO
                        iPos=2
                     Else If (iPos.eq.2) Then
                        Line(2+Magic*2:2+Magic*3-1)=ChCMO
                        iPos=3
                     Else If (iPos.eq.3) Then
                        Line(2+Magic*3:2+Magic*4-1)=ChCMO
                        iPos=4
                     Else If (iPos.eq.4) Then
                        Line(2+Magic*4:2+Magic*5-1)=ChCMO
                        iPos=1
                        Write (6,'(A)') Line
                        Line=' '
                     End If
                  End If
               End Do
               If (iPos.ne.1) Write (6,'(A)') Line
               Write (6,*)
            End If
           End Do
*
        Else
*
        Do 110 ISO=0,NO-1,NCOLS
          IEO=ISO+NCOLS-1
          IEO=MIN((NO-1),IEO)
          IEO2=IEO

cvv NAG
c          Do IO=ISO,IEO2
c            If (PrEne.and.PrOcc) Then
c               If ( Ene(ISX+IO).GT.ThrEne  .OR.
c     &              Occ(ISX+IO).LT.ThrOcc ) IEO=IEO-1
c            Else If (PrEne) Then
c               If ( Ene(ISX+IO).GT.ThrEne ) IEO=IEO-1
c            Else If (PrOcc) Then
c               If ( Occ(ISX+IO).LT.ThrOcc ) IEO=IEO-1
c            End If
c          End Do

            If (PrEne.and.PrOcc) Then
             Do IO=ISO,IEO2
               If ( Ene(ISX+IO).GT.ThrEne  .OR.
     &              Occ(ISX+IO).LT.ThrOcc ) IEO=IEO-1
             enddo
            Else If (PrEne) Then
             Do IO=ISO,IEO2
               If ( Ene(ISX+IO).GT.ThrEne ) IEO=IEO-1
             enddo
            Else If (PrOcc) Then
             Do IO=ISO,IEO2
               If ( Occ(ISX+IO).LT.ThrOcc ) IEO=IEO-1
             enddo
            End If
          If ( IEO.GE.ISO ) Then
             If (.Not.Header_Done) Then
                Write(6,'(/6X,A,I2,A,A)')
     &               'Molecular orbitals for symmetry species',iSym,
     &               ': ',IrrepName(iSym)
                Header_Done=.True.
             End If
             If (PrEne) Then
                tmp=0.0D0
                Do IO=ISO,IEO
                   Shift(IO-ISO+1)=1.0D0
                   tmp=Max(Abs(Ene(ISX+IO)),tmp)
                End Do
                If (tmp.gt.1.0D3) Then
                   Write (6,*)
                   Write (6,'(10X,A)') 'Some orbital energies have '//
     &                   'been scaled by powers of 10, the power is '//
     &                   'written right after the orbital index'
                   Write (6,*)
                   Do IO=ISO,IEO
                      tmp=Abs(Ene(ISX+IO))
                      If (tmp.gt.1.0D3) Then
                         tmp=tmp*1.0D-2
                         itmp=INT(LOG10(tmp))
                         Shift(IO-ISO+1)=10.0D0**itmp
                         Write (Star(IO-ISO+1),'(A,I1,A)')
     &                         ' (',itmp,')'
                      Else
                         Star(IO-ISO+1)='    '
                      End If
                   End Do
                Else
                   Do IO=ISO,IEO
                      Star(IO-ISO+1)='    '
                   End Do
                End If
             End If
             Write(6,*)
             Write(6,FMT0) LABEL0,(1+IO,Star(IO-ISO+1),IO=ISO,IEO)
             IF ( PrEne ) Write(6,FMT1) LABEL1,
     &          (Ene(ISX+IO)/Shift(IO-ISO+1),IO=ISO,IEO)
             IF ( PrOcc ) Write(6,FMT1) LABEL2,(Occ(ISX+IO),IO=ISO,IEO)
             Write(6,*)
             if(iPrForm.ne.1) then
              Do IB=0,NB-1
                Write(6,FMT2)
     &            IB+1,Clean_BName(Name(ISB+IB),LENIN),
     &            (CMO(ISCMO+IO*NB+IB),IO=ISO,IEO)
              End Do
             endif
           End If
*
*----------------------------------------------------------------------*
*     End of loop over columns                                         *
*----------------------------------------------------------------------*
*
110     CONTINUE
        End If
*
*----------------------------------------------------------------------*
*     Update list of pointers                                          *
*----------------------------------------------------------------------*
*
        ISX=ISX+NO
        ISB=ISB+NB
        ISCMO=ISCMO+NO*NB
*
*----------------------------------------------------------------------*
*     End of loop over all irreducible representations                 *
*----------------------------------------------------------------------*
*
100   CONTINUE
*
      Call CollapseOutput(0,'   Molecular orbitals:')
      Write (6,*)
*
      Return
      End
