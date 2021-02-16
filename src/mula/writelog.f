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
* Copyright (C) 1995, Niclas Forsberg                                  *
************************************************************************
C!-----------------------------------------------------------------------!
C!
c       Module InOutMod
C!
C!  Contains:
C!    ReadInp
C!    WriteLog
C!    WriteFC
C!    WriteHeader
C!    WrMold
C!
C!  Uses:
C!    Constants
C!    InputData
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
C!-----------------------------------------------------------------------!
C!
C!
c       Contains




C!-----------------------------------------------------------------------!
C!
      Subroutine WriteLog(PotCoef,AtomLbl,AtCoord,Mass,InterVec,
     &  stand_dev,
     &  max_err,energy,Hess,G,V,harmfreq,qMat,Bond,nBond,
     &  xvec,D3,D4,PED,x_anharm,anharmfreq,max_term,nState,
     &  ForceField,NumOfAt,nOsc)
C!
C!  Purpose:
C!    Write results to log file.
C!
C!  Input:
C!    AtomLbl    : Array of character - contains the labels for the
C!                 atoms.
C!    AtCoord    : Two dimensional Real*8 array - contains
C!                 the cartesian coordinates of the atoms.
C!    Mass       : Real*8 array - contains the mass of the
C!                 atoms.
C!    InterVec   : Integer array - contains the atoms that are used
C!                 in the calculations of each internal coordinate.
C!    ipow       : Two dimesional integer array - terms of polynomial.
C!    PotCoef    : Real*8 two dimensional array - coefficients
C!                 of terms given in ipow.
C!    stand_dev  : Real*8 variable - standard deviation of fitted
C!                 values of polynomial.
C!    max_err    : Real*8 variable - maximum error of fitted
C!                 values of polynomial.
C!    energy     : Real*8 variable - energy in minimum.
C!    Hess       : Real*8 array -  contains the force
C!                 constants expressed in internal coordinates.
C!    G          : Real*8 array.
C!    V          : Real*8 array - contains the eigenvectors
C!                 of G*F as columns.
C!    harmfreq   : Real*8 array - harmonic frequencies.
C!    qMat       : Real*8 array - contains the cartesian
C!                 displacements of the atoms.
C!    Bond       : Integer array - contains atom pairs that are to be
C!                 bonded together in a plot.
C!    nBond      : Integer - dim of Bond, i.e. 2*(number of bonds).
C!    xvec       : Real*8 array - contains the internal coordinates
C!                 at equilibrium.
C!    D3         : Real*8 array - contains the third derivatives.
C!    D4         : Real*8 array - contains the fourth derivatives.
C!    PED        : Real*8 three dimensional array - potential
C!                 energy distribution.
C!    x_anharm   : Real*8 two dimensional array - anharmonicity
C!                 constants.
C!    anharmfreq : Real*8 array - anharmonic frequencies.
C!    max_term   : Integer - highest power of term in fitted polynomial.
C!    nState     : Integer - 1 or 2, depending upon which state it is.
C!
C!  Output:
C!    Log file
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
      Implicit Real*8 ( a-h,o-z )
#include "Constants_mula.fh"
#include "dims.fh"
#include "inputdata.fh"
#include "indims.fh"
#include "inout.fh"
      Integer   VibPlotUnit
      Real*8   energy
      Real*8 Hess(nOsc,nOsc),G(nOsc,nOsc), V(nOsc,nOsc)
      Real*8 D3(nOsc,nOsc,nOsc)
      Real*8 D4(nOsc,nOsc,nOsc,nOsc)
      Real*8 qMat(3*NumOfAt,nOsc)
      Real*8 harmfreq(nOsc)
      Real*8 xvec(nOsc)
      Integer InterVec(MaxNumAt*15)
      Real*8 Mass(MaxNumAt)
      Real*8 AtCoord(3,NumOfAt)
      Character*4 AtomLbl(MaxNumAt)
      Real*8 PED(nOsc,nOsc,nOsc)
      Real*8 x_anharm(nOsc,nOsc)
      Real*8 anharmfreq(nOsc)
      Integer Bond(2*MaxNumAt)
c       Logical    Plot
      Logical    ForceField
      Real*8     stand_dev,max_err
      Character*11  VibPlotFile
      Character*8   BondString
C!
      Real*8 PotCoef(nPolyTerm,1)
#include "WrkSpc.fh"
C!
C!---- Format declarations.
      Character*32 F1
      Character*32 F2
      F1='(a2,a)'
      F2='(a2,a4,3f14.8,f20.8)'
C!
      NumInt  = nOsc
      l_aNormModes=NumInt
      Call GetMem('aNormModes','Allo','Inte',
     &  ipaNormModes,l_aNormModes)
      Do i = 1,NumInt
      iWork(ipaNormModes+i-1) = i
      End Do
C!
      If(nState.eq.1) Then
      Write(6,*)
      Write(6,*)
      Write(6,'(a27,a)') ' ',
     &       ' ================================================='
      Write(6,'(a27,a)') ' ',
     &       '|                                                 |'
      Write(6,'(a27,a)') ' ',
     &       '|           Results for the first state           |'
      Write(6,'(a27,a)') ' ',
     &       '|                                                 |'
      Write(6,'(a27,a)') ' ',
     &       ' ================================================='
      Write(6,*)
      Endif
      If(nState.eq.2) Then
      Write(6,*)
      Write(6,*)
      Write(6,*)
      Write(6,'(a27,a)') ' ',
     &       ' ================================================='
      Write(6,'(a27,a)') ' ',
     &       '|                                                 |'
      Write(6,'(a27,a)') ' ',
     &       '|           Results for the second state          |'
      Write(6,'(a27,a)') ' ',
     &       '|                                                 |'
      Write(6,'(a27,a)') ' ',
     &       ' ================================================='
      Write(6,*)
      End If
C!
C!---- Write coefficients of polynomial fit, standard deviation and
C!     maximum error to log.
      If ( .not.ForceField ) Then
      Write(6,*)
      Write(6,'(A)')
     &    '  Coefficients of different terms in polynomial:'
      Write(6,'(A)')
     &    '  ----------------------------------------------'
      Do i = 1,nPolyTerm
      Write(6,'(1x,f18.12,10(3x,i2))') PotCoef(i,1),
     &       (iWork(ipipow+i+nPolyTerm*(ivar-1)-1),ivar=1,nvar)
      End Do
      Write(6,*)
      Write(6,*)
      Write(6,'(A,es14.6)')
     &    '  Standard deviation of fitted values:',stand_dev
      Write(6,'(A)')
     &    '  ------------------------------------'
      Write(6,*)
      Write(6,'(A,es14.6)')
     &    '  Maximum error of fitted values:     ',max_err
      Write(6,'(A)') '  -------------------------------'
C!
C!---- Write energy in minimum.
      Write(6,*)
      Write(6,*)
      Write(6,'(A,F15.8,A)') '  Energy in minimum:',energy
      Write(6,'(A)') '  ------------------'
      End If
C!
C!---- Write force constant matrix to log file.
      If ( Huge_Print ) Then
      Write(6,*)
      Write(6,*)
      Write(6,*) ' ','Force Constant Matrix:'
      Write(6,*) ' ','======================'
      Write(6,*)
      maxCol = 10
      nRow = NumInt/maxCol+1
      nCol = Mod(NumInt,maxCol)
      m1 = 1
      m2 = 1
      cont = .true.
      Do i = 1,nRow
      If (( i.eq.nRow ).or.(( i.eq.1 ).and.
     &       ( NumInt.lt.maxCol ))) Then
      If ( i.eq.1 ) Then
      m2 = m2+nCol-1
      Else
      m2 = m2+nCol
      End If
      Else
      If ( i.eq.1 ) Then
      m2 = m2+maxCol-1
      Else
      m2 = m2+maxCol
      End If
      End If
      If (( nCol.eq.0 ).and.( i.eq.nRow )) cont = .false.
      If ( cont ) Then
      Write(6,'(10I12)') (nInt,nInt=m1,m2)
      Do mInt = 1,NumInt
      Write(6,'(10F12.6)') (Hess(mInt,nInt),nInt=m1,m2)
      End Do
      Write(6,*)
      Write(6,*)
      m1 = m2+1
      End If
      End Do
C!
C!---- Write matrix G to log file.
      Write(6,*)
      Write(6,*)
      Write(6,*) ' ','Inverse Mass Tensor :'
      Write(6,*) ' ','====================='
      Write(6,*)
      maxCol = 10
      nRow = NumInt/maxCol+1
      nCol = Mod(NumInt,maxCol)
      m1 = 1
      m2 = 1
      cont = .true.
      Do i = 1,nRow
      If (( i.eq.nRow ).or.(( i.eq.1 ).and.
     &       ( NumInt.lt.maxCol ))) Then
      If ( i.eq.1 ) Then
      m2 = m2+nCol-1
      Else
      m2 = m2+nCol
      End If
      Else
      If ( i.eq.1 ) Then
      m2 = m2+maxCol-1
      Else
      m2 = m2+maxCol
      End If
      End If
      If (( nCol.eq.0 ).and.( i.eq.nRow )) cont = .false.
      If ( cont ) Then
      Write(6,'(10I12)') (nInt,nInt=m1,m2)
      Do mInt = 1,NumInt
      Write(6,'(10F12.6)') (G(mInt,nInt),nInt=m1,m2)
      End Do
      Write(6,*)
      Write(6,*)
      m1 = m2+1
      End If
      End Do
      Write(6,*)
C!
C!---- Eigenvectors.
      Write(6,*)
      Write(6,*)
      Write(6,*) ' ','Eigenvectors :'
      Write(6,*) ' ','=============='
      Write(6,*)
      maxCol = 10
      nRow = NumInt/maxCol+1
      nCol = Mod(NumInt,maxCol)
      m1 = 1
      m2 = 1
      cont = .true.
      Do i = 1,nRow
      If (( i.eq.nRow ).or.(( i.eq.1 ).and.
     &       ( NumInt.lt.maxCol ))) Then
      If ( i.eq.1 ) Then
      m2 = m2+nCol-1
      Else
      m2 = m2+nCol
      End If
      Else
      If ( i.eq.1 ) Then
      m2 = m2+maxCol-1
      Else
      m2 = m2+maxCol
      End If
      End If
      If (( nCol.eq.0 ).and.( i.eq.nRow )) cont = .false.
      If ( cont ) Then
      Write(6,'(10I12)') (nInt,nInt=m1,m2)
      Do mInt = 1,NumInt
      Write(6,'(10F12.6)') (V(mInt,nInt),nInt=m1,m2)
      End Do
      Write(6,*)
      Write(6,*)
      m1 = m2+1
      End If
      End Do
      Write(6,*)
      End If
C!
C!---- Third derivatives.
      If ( max_term.gt.2 ) Then
      Write(6,*)
      Write(6,*)
      Write(6,'(A)') '  Cubic force constants:'
      Write(6,'(A)') '  ----------------------'
      Do i = 1,NumInt
      Do j = i,NumInt
      Do k = j,NumInt
      Write(6,'(a,i2,i2,i2,f15.8)') '  ',i,j,k,D3(i,j,k)
      End Do
      End Do
      End Do
      End If
C!
C!---- Fourth derivatives.
      If ( max_term.gt.3 ) Then
      Write(6,*)
      Write(6,*)
      Write(6,'(A)') '  Quartic force constants:'
      Write(6,'(A)') '  ------------------------'
      Do i = 1,NumInt
      Do j = i,NumInt
      Do k = j,NumInt
      Do l = k,NumInt
      Write(6,'(a,i2,i2,i2,i2,f15.8)')
     &                '  ',i,j,k,l,D4(i,j,k,l)
      End Do
      End Do
      End Do
      End Do
      End If
C!
C!---- Anharmonicity constants.
      If ( max_term.gt.2 ) Then
      Write(6,*)
      Write(6,*)
      Write(6,'(A)') '  Anharmonicity constants:'
      Write(6,'(A)') '  ------------------------'
      Do i = 1,NumInt
      Write(6,'(a,20f15.8)') ' ',(HarToRcm*x_anharm(i,j),j=1,i)
      End Do
      End If
C!
C!---- Write labels, coordinates and masses to log file.
      Call WriteCartCoord(AtomLbl,AtCoord,Mass,NumOfAt)
C!
C!---- Write internal coordinates to log file.
      Call WriteIntCoord(InterVec,AtomLbl,xvec,NumInt)
C!
C!---- Harmonic frequencies in reciprocal cm, GHz and hartrees.
      Call WriteFreq(harmfreq,iWork(ipaNormModes),
     &  l_aNormModes,'Harmonic frequencies')
C!
C!---- Fundamental frequencies in reciprocal cm, GHz and hartrees.
      If ( max_term.gt.2 )
     &       Call WriteFreq(anharmfreq,iWork(ipaNormModes),
     &   l_aNormModes,
     &     'Anharmonic frequencies')
C!
C!---- Write Potential Energy Distribution for each mode to log
C!     if the molecule is not too large.
      If ( Huge_Print .AND. NumInt.le.10 ) Then
      Write(6,*)
      Write(6,*)
      Write(6,*) ' ','Potential Energy Distribution :'
      Write(6,*) ' ','==============================='
      Do imode = 1,NumInt
      Write(6,*)
      Write(6,*) ' ','Mode',imode
      Write(6,*) ' ','-------'
      Do k = 1,NumInt
      If ( NumInt.lt.21 ) Then
      Write(6,'(60F6.2)') (PED(k,l,imode),l=1,NumInt)
      Else
      Write(6,'(60F5.2)') (PED(k,l,imode),l=1,NumInt)
      End If
      End Do
      End Do
      End If
C!
C!---- Find length of longest displacement vector.
      vMax = Dnrm2_(3,qMat(1,1),1)
      Do nAtom = 1,NumOfAt
      Do j = 1,NumInt
      vLength = Dnrm2_(3,qMat((nAtom-1)*3+1,j),1)
      If ( vLength.gt.vMax ) Then
      vMax = vLength
      End If
      End Do
      End Do
C!
C!---- Print displacement vectors scaled with length of longest vector
C!     if the molecule is not too large.
      If ( Huge_Print .AND. NumInt.le.10 ) Then
      Write(6,*)
      Write(6,*)
      Write(6,*) ' ','Scaled displacement vectors :'
      Write(6,*) ' ','============================='
      Write(6,*) ' ','(length of longest vector set equal to one)'
      Write(6,*)
      maxCol = 10
      Do nAtom = 1,NumOfAt
      Write(6,*) ' ',' Atom: ',AtomLbl(nAtom)
      Write(6,*) ' ',' ---------'
      nRow = NumInt/maxCol
      nCol = Mod(NumInt,maxCol)
      m1 = 1
      m2 = 0
      Do i = 0,nRow
      If ( i.eq.nRow ) Then
      m2 = m2+nCol
      Else
      m2 = m2+maxCol
      End If
      If ( m2.ne.0 ) Then
      Write(6,'(10I12)') (j,j=m1,m2)
      Write(6,'(A4,10F12.8)')
     &                    '  x ',(qMat((nAtom-1)*3+1,j)/vMax,j=m1,m2)
      Write(6,'(A4,10F12.8)')
     &                    '  y ',(qMat((nAtom-1)*3+2,j)/vMax,j=m1,m2)
      Write(6,'(A4,10F12.8)')
     &                    '  z ',(qMat((nAtom-1)*3+3,j)/vMax,j=m1,m2)
      Write(6,*)
      m1 = m2+1
      End If
      End Do
      End Do
      Write(6,*)
      Write(6,*)
      End If
C!
C! Generate MOLDEN input:
C!
      If( nState.eq.1) Then
      Call WrMold('MOLDEN-1',NumOfAt,AtomLbl,AtCoord,
     &   NumInt,HarmFreq,QMat)
      Else
      Call WrMold('MOLDEN-2',NumOfAt,AtomLbl,AtCoord,
     &   NumInt,HarmFreq,QMat)
      End If
C!
C!---- Print geometry of molecule, together with cartesian displacements,
C!     in a file.
      VibPlotUnit=77
      If ( VibModPlot ) Then
      If (nState.eq.1) Then
      VibPlotUnit = VibPlotUnit1
      VibPlotFile = 'plot.modes1'
      Endif
      If (nState.eq.2) Then
      VibPlotUnit = VibPlotUnit2
      VibPlotFile = 'plot.modes2'
      End If
      call molcas_Open(vibplotunit,vibplotfile)
c      Open (Unit=VibPlotUnit,file=VibPlotFile)
      Write(VibPlotUnit,*) 2*NumOfAt
      Write(VibPlotUnit,*) nBond/2
      Write(VibPlotUnit,*) '1'
      Write(VibPlotUnit,*) NumInt
      iBond = 1
      Do While ( iBond+1.le.nBond )
      Write(BondString,'(i4,i4)') Bond(iBond),Bond(iBond+1)
      Write(VibPlotUnit,'(a8)') BondString
      iBond = iBond+2
      End Do
      Do nAtom = 1,NumOfAt
      x1 = AtCoord(1,nAtom)
      x2 = AtCoord(2,nAtom)
      x3 = AtCoord(3,nAtom)
      Write(VibPlotUnit,'(a,3f12.3)') '>',x2,x3,x1
      End Do
      Do iInt = 1,NumInt
      Do nAtom = 1,NumOfAt
      x1 = AtCoord(1,nAtom)
      x2 = AtCoord(2,nAtom)
      x3 = AtCoord(3,nAtom)
      Write(VibPlotUnit,'(a,3f12.3)') AtomLbl(nAtom),
     &                     x2+qMat((nAtom-1)*3+2,iInt)/vMax,
     &                     x3+qMat((nAtom-1)*3+3,iInt)/vMax,
     &                     x1+qMat((nAtom-1)*3+1,iInt)/vMax
      End Do
      End Do
      Close (VibPlotUnit)
      End If
C!
      Call GetMem('aNormModes','Free','Inte',
     &  ipaNormModes,l_aNormModes)
C!
      End
