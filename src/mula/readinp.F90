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
! Copyright (C) 1995, Niclas Forsberg                                  *
!               2017, Ignacio Fdez. Galvan                             *
!***********************************************************************
!!-----------------------------------------------------------------------!
!!
      Subroutine ReadInp(Title,AtomLbl,Mass,InterVec,Bond,nBond,        &
     &   NumInt,NumOfAt,                                                &
     &   trfName1,trfName2,m_max,n_max,max_dip,max_term,MatEl,          &
     &   ForceField,Cartesian,lExpan,lISC,iCode,dMinWind)
!!
!!  Purpose:
!!    Read the input file.
!!
!!  Output:
!!    Title    : String - title of the project.
!!    AtomLbl  : Array of character - contains the labels for the
!!               atoms.
!!    Mass     : Real*8 array - contains the mass of the
!!               atoms.
!!    InterVec : Integer array - containis the atoms that are used
!!               in the calculations of each internal coordinate.
!!    Bond     : Integer array - contains atom pairs that are to be
!!               bonded together in a plot.
!!    nBond    : Integer - dim of Bond, i.e. 2*(number of bonds).
!!    NumInt   : Integer - the total number of internal coordinates.
!!    NumOfAt  : Integer - the number of atoms.
!!    trfName  : Character array - type of transformation of variables.
!!    m_max    : Integer - maximum level for the first state.
!!    n_max    : Integer - maximum level for the second state.
!!    m_plot   : Integer array - level(s) to plot for the first state.
!!    n_plot   : Integer array - level(s) to plot for the second state.
!!    max_dip  : Integer - highest order of term in transition dipole.
!!    max_term : Integer - highest power of a term in polynomial fitted
!!               to energy values.
!!    MatEl    : Logical
!!    lISC     : Logical to calculate InterSystem Crossing
!!
!!  Uses:
!!    IOTools
!!
!!  Written by:
!!    Niclas Forsberg,
!!    Dept. of Theoretical Chemistry, Lund University, 1995.
!!
!!  Modified to use isotopes module
!!    Ignacio Fdez. Galvan, 2017
!!
      Use Isotopes
      Implicit Real*8 ( a-h,o-z )
#include "inout.fh"
#include "Constants_mula.fh"
#include "inputdata.fh"
#include "indims.fh"
      Integer InterVec(MaxNumAt*15)
      Character*80  trfName1(MaxNumAt),trfName2(MaxNumAt)
      Real*8  Mass(MaxNumAt)
      Character*4  AtomLbl(MaxNumAt)
      Integer Bond(MaxNumAt*2)
      Character*80   Title
!       Logical   Plot
      Character*9  CoordType
      Character*1  c
      Character*7  AtInp
      Character*2  AtName
      Character*2  Eunit
      Character*3  DegOrRad
      Character*2  AAOrAu
      Character*4  Atom
      Character*6  PlUnit
      Character*4  Atom1,Atom2, Atom3,Atom4
!      Character*2 AtList( 0:120 )
!      Integer AtData( 120,4 )
!      Real*8 MassData (400)
      Real*8 coord(1000,10)
      Real*8 GrdVal(1000,10)
      Logical  MatEl,ForceField,lISC
      Logical  Cartesian
      Logical  exist
      Logical  lExpan
      Real*8   m1,m2,m3,m12,m23
      Real*8 TmpCoord(3)
      parameter (max_len_xvec=150)
      Real*8 xvec(max_len_xvec)
      parameter (max_plot_temp=100)
      Integer plot_temp(max_plot_temp), iCode

      Logical  XnotFound, YnotFound,ZnotFound
!!- Format declarations.
      Character*8 Format
!!- User-defined functions called:
      External StrToDble, iStrToInt
!!
#include "WrkSpc.fh"
!!
      Format='(A80)'
!!
!!---- Read from MassFile:
!!     - name of atoms into array of string - AtList.
!!     - atomic number, default mass number, smallest mass number and
!!       table offset into two dimensional array - AtData.
!!     - masses of all isotopes into array - MassData.
!!
!     i = 0
!     Call Molcas_Open(massUnit,'MASSUNIT')
!     Read(massUnit,Format) InLine
! Read sequential lines from atomic data file.
!     Call Normalize(InLine,OutLine)
!     Do While ( OutLine(1:4) .ne. 'END ')
!     l = Index(OutLine,'*')
!     If ( l .eq. 0 ) Then
!     i = i+1
!     AtList(i)=OutLine(1:2)
!     Read(OutLine(3:80),*) (AtData(i,j),j=1,4)
!     End If
!     Read(massUnit,Format) InLine
!     Call Normalize(InLine,OutLine)
!     End Do
!     iAtList=i
!!
!     i = 0
!     Read(massUnit,Format) InLine
!     Call Normalize(InLine,OutLine)
!     Do While ( OutLine(1:4) .ne. 'END ')
!     l = Index(OutLine,'*')
!     If ( l .eq. 0 ) Then
!     i = i+1
!     Read(OutLine,*) MassData(i)
!     End If
!     Read(massUnit,Format) InLine
!     Call Normalize(InLine,OutLine)
!     End Do
!     close(massUnit)
      iPrint = iPrintLevel(-1)

! ----------------------------------------------------------------------
! --- TITLe: Read the title into a string - Title
!
      Title=' '
      Call KeyWord(inpUnit,'TITL',.true.,exist)
      If (exist) Read(inpUnit,Format) Title

! ----------------------------------------------------------------------
! --- ATOMs:  Process ATOMS
!
      Call KeyWord(inpUnit,'ATOM',.true.,exist)
      NumOfAt = 0
      Read(inpUnit,Format) InLine
      Call Normalize(InLine,OutLine)
      Do While ( OutLine(1:4) .ne. 'END ')
        NumOfAt = NumOfAt+1
        Read(inpUnit,Format) InLine
        Call Normalize(InLine,OutLine)
      End Do
!!
!! Read the labels of atoms into an array - AtomLbl.
!!
      Call KeyWord(inpUnit,'ATOM',.true.,exist)
      Do nAtom = 1,NumOfAt
        Read(inpUnit,Format) InLine
        Call Normalize(InLine,OutLine)
        k = 1
        Call WordPos(k,OutLine,iStart,iStop)
        k = iStop+1
        AtInp = '       '
        AtInp = OutLine(iStart:iStop)
!! ---- Check if a massnumber is given.
        Num = 0
        Do i=1,7
          c = AtInp(i:i)
          IntVal = Index('0123456789',c)-1
          If( IntVal .ge. 0 ) then
            Num = 10*Num+IntVal
          Else
            goto 600
          Endif
        End Do
600     If ( (iStart+i-1) .gt. iStop ) Then
          Call WordPos(k,OutLine,iStart,iStop)
          AtInp = '      '
          AtInp = OutLine(iStart:iStop)
          i = 1
          c = AtInp(i:i)
        End If
!! ---- Extract atomic name and possible label.
        AtName = '  '
        AtomLbl(nAtom) = '    '
        j = 1
        do ii=i,7
          IntVal = Index('0123456789',c)-1
          if (( IntVal .lt. 0 ).and.( c .ne. ' ' )) then
            AtName(j:j) = c
            AtomLbl(nAtom)(j:j) = c
            j = j+1
            i = i+1
            c = AtInp(i:i)
          endif
        End Do
        do ii=i,7
          IntVal = Index('0123456789',c)-1
          if (( (iStart+i-1) .le. iStop ).and.( IntVal .ge. 0 )) then
            AtomLbl(nAtom)(j:j) = c
            j = j+1
            i = i+1
            c = AtInp(i:i)
          endif
        End Do
!! ---- Remove massnumber, atom and label from OutLine.
        Length = Len(OutLine)
        OutLine = OutLine(iStop+1:Length)
!!
!! ---- If mass is given in input, then read mass, else use MassData.
        k = 1
        Call WordPos(k,OutLine,iStart,iStop)
        k = iStop
        If ( k .ne. Length ) Then
          Read(OutLine,*) Mass(nAtom)
        Else
!          Do i=1,iAtList
!            if (AtList(i) .eq. AtName ) goto 1111
!          End Do
!1111      continue
!          If ( Num .eq. 0 ) Then
!            Num = AtData(i,2)
!          End If
!          Mass(nAtom) = MassData(AtData(i,4)+Num+1-AtData(i,3))
          Call Isotope(Num,AtName,Mass(nAtom))
          Mass(nAtom)=Mass(nAtom)/uToau
!!  ----  Print error message if the isotope is not listed in
!!        MassFile.
!          nOffset = AtData(i+1,4)
!          If (( Num .gt. nOffset ).or.(Mass(nAtom) .eq. (-1.0d0))) Then
!            Write(6,*)
!            Write(6,*) ' *************** ERROR *****************'
!            Write(6,'(A,A2,A,I3,A)') ' The isotope of ',AtName,
!     &                               ' with mass number ',Num
!            Write(6,*) ' is not listed in MassFile.             '
!            Write(6,*) '****************************************'
!            Call Quit_OnUserError()
!          End If
        End If
      End Do
!
! --- ATOMs: End -------------------------------------------------------

! ----------------------------------------------------------------------
! --- INTErnal: Resolve the different internal coordinates specified
!               in the input.
!
      Call KeyWord(inpUnit,'INTE',.true.,exist)
      j = 1
      NumInt = 0
      nBond = 0
      Read(inpUnit,Format) InLine
      Call Normalize(InLine,OutLine)

      DO WHILE ( OutLine(1:4) .ne. 'END ')

!!---- Bond Stretching.
      l = Index(OutLine,'BOND')
      If ( l .ne. 0 ) Then
      k = l+Len('BOND')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      InterVec(j) = 1
      Do i = 1,NumOfAt
        If ( Atom1 .eq. AtomLbl(i) ) Then
          InterVec(j+1) = i
        Else If ( Atom2 .eq. AtomLbl(i) ) Then
          InterVec(j+2) = i
        End If
      End Do
      Bond(nBond+1) = InterVec(j+1)
      Bond(nBond+2) = InterVec(j+2)
      nBond = nBond+2
      j = j+3
      NumInt = NumInt+1
      End If
!!---- Valence Angle Bending.
      l = Index(OutLine,'ANGLE')
      If ( l .ne. 0 ) Then
      k = l+Len('ANGLE')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom3 = OutLine(iStart:iStop)
      InterVec(j) = 2
      Do i = 1,NumOfAt
        If ( Atom1 .eq. AtomLbl(i) ) Then
          InterVec(j+1) = i
        Else If ( Atom2 .eq. AtomLbl(i) ) Then
          InterVec(j+2) = i
        Else If ( Atom3 .eq. AtomLbl(i) ) Then
          InterVec(j+3) = i
        End If
      End Do
      j = j+4
      NumInt = NumInt+1
      End If
!!---- Linear Valence Angle.
      l = Index(OutLine,'LINANG')
      If ( l .ne. 0 ) Then
      k = l+Len('LINANG')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom3 = OutLine(iStart:iStop)
      InterVec(j) = 3
      Do i = 1,NumOfAt
        If ( Atom1 .eq. AtomLbl(i) ) Then
          InterVec(j+1) = i
        Else If ( Atom2 .eq. AtomLbl(i) ) Then
          InterVec(j+2) = i
        Else If ( Atom3 .eq. AtomLbl(i) ) Then
          InterVec(j+3) = i
        End If
      End Do
      j = j+4
      NumInt = NumInt+2
      End If
!!---- Torsion.
      l = Index(OutLine,'TORSION')
      If ( l .ne. 0 ) Then
      k = l+Len('TORSION')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom3 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom4 = OutLine(iStart:iStop)
      InterVec(j) = 4
      Do i = 1,NumOfAt
        If ( Atom1 .eq. AtomLbl(i) ) Then
          InterVec(j+1) = i
        Else If ( Atom2 .eq. AtomLbl(i) ) Then
          InterVec(j+2) = i
        Else If ( Atom3 .eq. AtomLbl(i) ) Then
          InterVec(j+3) = i
        Else If ( Atom4 .eq. AtomLbl(i) ) Then
          InterVec(j+4) = i
        End If
      End Do
      j = j+5
      NumInt = NumInt+1
      End If
!!---- Out of Plane Angle.
      l = Index(OutLine,'OUTOFPL')
      If ( l .ne. 0 ) Then
      k = l+Len('OUTOFPL')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom3 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom4 = OutLine(iStart:iStop)
      InterVec(j) = 5
      Do i = 1,NumOfAt
        If ( Atom1 .eq. AtomLbl(i) ) Then
          InterVec(j+1) = i
        Else If ( Atom2 .eq. AtomLbl(i) ) Then
          InterVec(j+2) = i
        Else If ( Atom3 .eq. AtomLbl(i) ) Then
          InterVec(j+3) = i
        Else If ( Atom4 .eq. AtomLbl(i) ) Then
          InterVec(j+4) = i
        End If
      End Do
      j = j+5
      NumInt = NumInt+1
      End If

      Read(inpUnit,Format) InLine
      Call Normalize(InLine,OutLine)
      END DO

      Call GetMem('AtCoord1','Allo','Real',ipAtCoord1,3*NumOfAt)
      Call GetMem('AtCoord2','Allo','Real',ipAtCoord2,3*NumOfAt)
      l_Hess1=NumInt
      Call GetMem('Hess1','Allo','Real',ipHess1,l_Hess1*l_Hess1)
      Call GetMem('Hess2','Allo','Real',ipHess2,l_Hess1*l_Hess1)
      n = 3*NumOfAt
      Call GetMem('TranDipGrad','Allo','Real',ipTranDipGrad,3*n)
!
! --- INTErnal: End ----------------------------------------------------

! ----------------------------------------------------------------------
! --- ISC: InterSystem Crossing. GG-30-Dec-08                     ! CGGn
!
      lISC=.False.
      dMinWind = 0.0d0
      Call KeyWord(inpUnit,'ISC ',.true.,lISC)
      If (lISC) then
         Call KeyWord(inpUnit,'EXPF',.true.,exist)
         If (exist) then
           Read(InpUnit,Format) InLine
           Call Normalize(InLine,Outline)
           Read(OutLine,*) dMinWind
         EndIf
         Call KeyWord(inpUnit,'OLDC',.true.,exist)
         If (exist) iCode = iCode  +  1
         Call KeyWord(inpUnit,'DISK',.true.,exist)
         If (exist) iCode = iCode  + 10
      EndIf
! --- ISC: End ---------------------------------------------------------

! ----------------------------------------------------------------------
! --- MODEs: Chose which modes to use in intensity calculations.
!
      Call KeyWord(inpUnit,'MODE',.true.,exist)
      If ( exist ) Then
      Call GetMem('tempmodes','Allo','Inte',iptempModes,NumInt)
      Read(inpUnit,Format) InLine
      Call Normalize(InLine,OutLine)
      i = 1
      Do While ( OutLine(1:4) .ne. 'END ')
        k = 1
        Call WordPos(k,OutLine,iStart,iStop)
        Do While ( iStop .lt. len(OutLine) )
          iWork(iptempModes+i-1) =                                      &
     &          iStrToInt(OutLine(iStart:iStop))
          If ( iWork(iptempModes+i-1) .gt. NumInt ) Then
            Write(6,*)
            Write(6,*) ' *********** ERROR *************'
            Write(6,*) ' Too many normal modes specified'
            Write(6,*) ' *******************************'
            Call Quit_OnUserError()
          End If
          i = i+1
          Call WordPos(k,OutLine,iStart,iStop)
        End Do
        Read(inpUnit,Format) InLine
        Call Normalize(InLine,OutLine)
      End Do
      n = i-1
      l_NormModes=n
      Call GetMem('NormModes','Allo','Inte',                            &
     &    ipNormModes,l_NormModes)
      do iv=1,n
        iWork(ipNormModes+iv-1) = iWork(iptempModes+iv-1)
      enddo
      Call GetMem('tempmodes','Free','Inte',iptempModes,NumInt)
      Do i = 1,n-1
        Do j = i+1,n
          If ( iWork(ipNormModes+j-1) .lt.                              &
     &          iWork(ipNormModes+i-1) ) Then
            modeTemp     = iWork(ipNormModes+i-1)
            iWork(ipNormModes+i-1) = iWork(ipNormModes+j-1)
            iWork(ipNormModes+j-1) = modeTemp
          End If
        End Do
      End Do
      Else
! --- No MODEs specified
        l_NormModes=NumInt
        Call GetMem('NormModes','Allo','Inte',                          &
     &    ipNormModes,l_NormModes)
        Do i = 1,NumInt
          iWork(ipNormModes+i-1) = i
        End Do
        If (iPrint.GE.1) then
          Write(6,*) ' *** Warning: MODEs not specified !'
          Write(6,*) '     All ',l_NormModes,' will be calculated.'
          Write(6,*)
        EndIf
      End If
!
! --- MODEs: End -------------------------------------------------------

! ----------------------------------------------------------------------
! --- MXLEvels: Maximun number of vibr. quanta in first and second state
!
      Call KeyWord(inpUnit,'MXLE',.true.,exist)
      m_max = 0 ! Defaul
      n_max = 1 ! Defaul
      If(.NOT.exist) then
        If (iPrint.GE.1) then
          Write(6,*) ' *** Warning: MXLEvels not specified !'
          Write(6,*) '     Default values are 0 and 1.      '
          Write(6,*)
        EndIf
      else
        Read(inpUnit,*) m_max,n_max
      Endif
!
! --- MXLEvels: End ----------------------------------------------------

! ----------------------------------------------------------------------
! --- VARIational: Check if a full matrix element calculation or a
!                  simpler harmonic approximation is wanted.
!
      MatEl = .false.
      Call KeyWord(inpUnit,'VARI',.true.,MatEl)
!
! --- VARIational: End -------------------------------------------------

! ----------------------------------------------------------------------
! --- TRANsitions: Check which transitions that are wanted in the
!                  output.
      Call KeyWord(inpUnit,'TRAN',.true.,exist)
      If(.NOT.exist) then
        If (.NOT.lISC) then
          Write(6,*) ' *** Warning: TRANsitions not specified!'
          Write(6,*) '     All Levels will be printed.        '
          Write(6,*)
        EndIf
        n = m_max + 1
        Do i = 1, n
          plot_temp(i) = i-1
        EndDo
      else
        Call KeyWord(inpUnit,'FIRS',.false.,exist)
        Read(inpUnit,Format) InLine
        Call Normalize(InLine,OutLine)
        k = 1
        i = 1
        Call WordPos(k,OutLine,iStart,iStop)
        Do While ( iStop .lt. len(OutLine) )
          if(i .gt. max_plot_temp) then
            Write(6,*)
            Write(6,*) ' ************* ERROR **************'
            write(6,*) ' TRANsitions: i .gt. max_plot_temp '
            Write(6,*) ' **********************************'
            Call Quit_OnUserError()
          endif
          plot_temp(i) = iStrToInt(OutLine(iStart:iStop))
          If ( plot_temp(i) .gt. m_max ) Then
            Write(6,*)
            Write(6,*) ' ************* ERROR **************'
            Write(6,*) ' A quantum specified for the output'
            Write(6,*) ' is larger than max quantum.       '
            Write(6,*) ' **********************************'
            Call Quit_OnUserError()
          End If
          i = i+1
          Call WordPos(k,OutLine,iStart,iStop)
        End Do
        n = i-1
      EndIf
      l_m_plot=n
      Call GetMem('m_plot','Allo','Inte',ipm_plot,l_m_plot)
      do iv=1,n
        iWork(ipm_plot+iv-1) = plot_temp(iv)
      enddo

      Call KeyWord(inpUnit,'TRAN',.true.,exist)
      If(.NOT.exist) then
        n = n_max + 1
        Do i = 1, n
          plot_temp(i) = i-1
        EndDo
      else
        Call KeyWord(inpUnit,'SECO',.false.,exist)
        Read(inpUnit,Format) InLine
        Call Normalize(InLine,OutLine)
        k = 1
        i = 1
        Call WordPos(k,OutLine,iStart,iStop)
        Do While ( iStop .lt. len(OutLine) )
          If(i .gt. max_plot_temp) then
            Write(6,*)
            Write(6,*) ' ************* ERROR **************'
            write(6,*) ' TRANsitions: i .gt. max_plot_temp '
            Write(6,*) ' **********************************'
            Call Quit_OnUserError
          EndIf
          plot_temp(i) = iStrToInt(OutLine(iStart:iStop))
          If ( plot_temp(i) .gt. n_max ) Then
            Write(6,*)
            Write(6,*) ' ************* ERROR **************'
            Write(6,*) ' A quantum specified for the output'
            Write(6,*) ' is larger than max quantum.'
            Write(6,*) ' **********************************'
            Call Quit_OnUserError()
          EndIf
          i = i+1
          Call WordPos(k,OutLine,iStart,iStop)
        End Do
        n = i-1
      EndIf
      l_n_plot=n
      Call GetMem('n_plot','Allo','Inte',ipn_plot,l_n_plot)
      do iv=1,n
        iWork(ipn_plot+iv-1) = plot_temp(iv)
      enddo
!
! --- TRANsitions: End -------------------------------------------------

! -----------------------------------------------------------------------!
!
!                            Forcefield part
!
! -----------------------------------------------------------------------!
!
      Call KeyWord(inpUnit,'FORC',.true.,Forcefield)
      If ( .not.Forcefield ) then
        Write(6,*) ' Forcefield not found. Will use Energy Surface.'
        Goto 10
      EndIf

! ----------------------------------------------------------------------
! --- ENERgy: Read minimum energies for the two surfaces.
!
      Call KeyWord(inpUnit,'ENER',.true.,exist)
      If(.NOT.exist) then
        Write(6,*)
        Write(6,*) ' ************** ERROR *************'
        Write(6,*) ' Electronic T_e energies not given.'
        Write(6,*) ' **********************************'
        Call Quit_OnUserError()
      EndIf
      Call KeyWord(inpUnit,'FIRS',.false.,exist)
      If(.NOT.exist) then
        Write(6,*)
        Write(6,*) ' ******************** ERROR *********************'
        Write(6,*) ' Electronic T_e energy for first state not given.'
        Write(6,*) ' ************************************************'
        Call Quit_OnUserError()
      EndIf
!      Read(inpUnit,*) energy1,Eunit
!      Call Upcase(Eunit)
! Replace above two lines with:
      Read(InpUnit,Format) InLine
      Call Normalize(InLine,Outline)
      Read(OutLine,*) Energy1
      EUnit='AU'
      If(Index(OutLine,'EV') .gt. 0) EUnit='EV'
      If(Index(OutLine,'CM') .gt. 0) EUnit='CM'
! End of replacement.

      rfact=1.0d0
      If ( Eunit .eq. 'AU' ) rfact = 1.0d0
      If ( Eunit .eq. 'EV' ) rfact = 1.0d0/27.2114d0
      If ( Eunit .eq. 'CM' ) rfact = 1.0d0/HarToRcm
      energy1 = energy1*rfact
      Call KeyWord(inpUnit,'SECO',.false.,exist)
      If(.NOT.exist) then
        Write(6,*)
        Write(6,*) ' ******************** ERROR **********************'
        Write(6,*) ' Electronic T_e energy for second state not given.'
        Write(6,*) ' *************************************************'
        Call Quit_OnUserError()
      EndIf
!      Read(inpUnit,*) energy2,Eunit
!      Call Upcase(Eunit)
! Replace above two lines with:
      Read(InpUnit,Format) InLine
      Call Normalize(InLine,Outline)
      Read(OutLine,*) Energy2
      EUnit='AU'
      If(Index(OutLine,'EV') .gt. 0) EUnit='EV'
      If(Index(OutLine,'CM') .gt. 0) EUnit='CM'
! End of replacement.
      If ( Eunit .eq. 'AU' ) rfact = 1.0d0
      If ( Eunit .eq. 'EV' ) rfact = 1.0d0/27.2114d0
      If ( Eunit .eq. 'CM' ) rfact = 1.0d0/HarToRcm
      energy2 = energy2*rfact
!
! --- ENERgy: : End ----------------------------------------------------

! ----------------------------------------------------------------------
! --- GEOMetries: Read geometries for the two surfaces.
!
      Call KeyWord(inpUnit,'GEOM',.true.,exist)
      If(.NOT.exist) then
        Write(6,*)
        Write(6,*) ' ***** ERROR *******'
        Write(6,*) ' GEOMetry not found.'
        Write(6,*) ' *******************'
        Call Quit_OnUserError()
      EndIf
!!D Write(6,*)' READINP: Read geometry (''GEOMETRY'').'
!      Read(inpUnit,*) CoordType
!      Call UpCase(CoordType)
! Replace above two lines with:
      Read(InpUnit,Format) InLine
      Call Normalize(InLine,Outline)
      iend=Index(OutLine,' ')
      CoordType=OutLine(1:iend-1)
! End of replacement.

      If (CoordType .eq. 'FILE') Then
        Coordtype='CARTESIAN'
        inpUnit1 = 31
        inpUnit2 = 32
        call molcas_open(31,'UNSYM1')
!      Open (unit=31,file='UNSYM1')
        Call KeyWord(inpUnit1,'*LABEL COORDINATES CHARGE',              &
     &    .false.,exist)
        call molcas_open(32,'UNSYM2')
!      Open (unit=32,file='UNSYM2')
        Call KeyWord(inpUnit2,'*LABEL COORDINATES CHARGE',              &
     &    .false.,exist)
      Else
        inpUnit1 = inpunit
        inpUnit2 = inpunit
      End If
!!
      If ( CoordType .eq. 'CARTESIAN' ) Then
        Cartesian = .true.
!!D Write(6,*)' READINP: Read cartesian coordinates for first state.'
        Call KeyWord(inpUnit,'FIRS',.false.,exist)
        If(.NOT.exist) then
          Write(6,*)
          Write(6,*) ' ****** ERROR ********'
          Write(6,*) ' Wrong geometry input.'
          Write(6,*) ' *********************'
          Call Quit_OnUserError()
        EndIf
        ierror=0
        Do iAtom = 1,NumOfAt
          Read(inpUnit1,Format) InLine
          Call UpCase(InLine)
          Atom=InLine(1:4)
          Read(InLine(5:80),*) (TmpCoord(i),i=1,3)
          j=10000000
          Do i=1,NumOfAt
            If (Atom .eq. AtomLbl(i)) j=i
          End Do
          If(j .gt. NumOfAt) then
            Write(6,*)
            Write(6,*)'*********** ERROR ***********'
            Write(6,*)' Unrecognized atom label '''//Atom//''''
            Write(6,*)'*****************************'
            ierror=ierror+1
          else
            do iv=1,3
              Work(ipAtCoord1+iv+3*(j-1)-1) = TmpCoord(iv)
            enddo
          End if
        End Do
!!D Write(6,*)' READINP: Read cartesian coordinates for second state.'
        Call KeyWord(inpUnit,'SECO',.false.,exist)
        If(.NOT.exist) then
          Write(6,*)
          Write(6,*) ' ****** ERROR ********'
          Write(6,*) ' Wrong geometry input.'
          Write(6,*) ' *********************'
          Call Quit_OnUserError()
        EndIf
        Do iAtom = 1,NumOfAt
          Read(inpUnit2,Format) InLine
          Call UpCase(InLine)
          Atom=InLine(1:4)
          Read(InLine(5:80),*) (TmpCoord(i),i=1,3)
          j = 1
          Do While ( Atom .ne. AtomLbl(j) )
            j = j+1
          End Do
          If(j .gt. NumOfAt) then
            Write(6,*)
            Write(6,*)'*********** ERROR ***********'
            Write(6,*)' Unrecognized atom label '''//Atom//''''
            Write(6,*)'*****************************'
            ierror=ierror+1
          else
            do iv=1,3
              Work(ipAtCoord2+iv+3*(j-1)-1) = TmpCoord(iv)
            enddo
          end if
        End Do

      Else If ( CoordType .eq. 'INTERNAL' ) Then

      Cartesian = .false.
      Call KeyWord(inpUnit,'FIRS',.false.,exist)
!!D Write(6,*)' READINP: Read internal coordinates for first state.'
      n = 5*(3*NumOfAt-5)
!         ! Largest possible necessary GeoVec
      Call GetMem('GeoVec','Allo','Inte',ipGeoVec,5*(3*NumOfAt-5))
      iInt = 1
      j = 1
      Read(inpUnit,Format) InLine
      Call Normalize(InLine,OutLine)
      Do While ( OutLine(1:4) .ne. 'END ')
      !---- Bond distance.
      l = Index(OutLine,'BOND')
      If ( l .ne. 0 ) Then
      k = l+Len('BOND')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      iwork(ipGeoVec+j-1) = 1
      Do i = 1,NumOfAt
      If ( Atom1 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j) = i
      Else If ( Atom2 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+1) = i
      End If
      End Do
      Call WordPos(k,OutLine,iStart,iStop)
      xvec(iInt) = StrToDble(OutLine(iStart:iStop))
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      AAOrAu = OutLine(iStart:iStop)
      If ( AAOrAu .eq. 'AA' ) Then
      xvec(iInt) = xvec(iInt)/Angstrom
      End If
      j = j+3
      End If
      !---- Valence Angle.
      l = Index(OutLine,'ANGLE')
      If ( l .ne. 0 ) Then
      k = l+Len('ANGLE')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom3 = OutLine(iStart:iStop)
      iWork(ipGeoVec+j-1) = 2
      Do i = 1,NumOfAt
      If ( Atom1 .eq. AtomLbl(i) ) Then
      iWOrk(ipGeoVec+j) = i
      Else If ( Atom2 .eq. AtomLbl(i) ) Then
      iWOrk(ipGeoVec+j+1) = i
      Else If ( Atom3 .eq. AtomLbl(i) ) Then
      iWOrk(ipGeoVec+j+2) = i
      End If
      End Do
      Call WordPos(k,OutLine,iStart,iStop)
      xvec(iInt) = StrToDble(OutLine(iStart:iStop))
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      DegOrRad = OutLine(iStart:iStop)
      If ( DegOrRad .eq. 'DEG' ) Then
      xvec(iInt) = xvec(iInt)*rpi/180.0d0
      End If
      j = j+4
      End If
      !---- Linear Valence Angle.
      l = Index(OutLine,'LINANG')
      If ( l .ne. 0 ) Then
      k = l+Len('LINANG')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom3 = OutLine(iStart:iStop)
      iWork(ipGeoVec+j-1) = 3
      Do i = 1,NumOfAt
      If ( Atom1 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j) = i
      Else If ( Atom2 .eq. AtomLbl(i) ) Then
      iWOrk(ipGeoVec+j+1) = i
      Else If ( Atom3 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+2) = i
      End If
      End Do
      Call WordPos(k,OutLine,iStart,iStop)
      xvec(iInt) = StrToDble(OutLine(iStart:iStop))
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      DegOrRad = OutLine(iStart:iStop)
      If ( DegOrRad .eq. 'DEG' ) Then
      xvec(iInt) = xvec(iInt)*rpi/180.0d0
      End If
      j = j+4
      End If
      !---- Dihedral angle.
      l = Index(OutLine,'TORSION')
      If ( l .ne. 0 ) Then
      k = l+Len('TORSION')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom3 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom4 = OutLine(iStart:iStop)
      iWork(ipGeoVec+j-1) = 4
      Do i = 1,NumOfAt
      If ( Atom1 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j) = i
      Else If ( Atom2 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+1) = i
      Else If ( Atom3 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+2) = i
      Else If ( Atom4 .eq. AtomLbl(i) ) Then
      iWOrk(ipGeoVec+j+3) = i
      End If
      End Do
      Call WordPos(k,OutLine,iStart,iStop)
      xvec(iInt) = StrToDble(OutLine(iStart:iStop))
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      DegOrRad = OutLine(iStart:iStop)
      If ( DegOrRad .eq. 'DEG' ) Then
      xvec(iInt) = xvec(iInt)*rpi/180.0d0
      End If
      j = j+5
      End If
      !---- Out of Plane Angle.
      l = Index(OutLine,'OUTOFPL')
      If ( l .ne. 0 ) Then
      k = l+Len('OUTOFPL')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom3 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom4 = OutLine(iStart:iStop)
      iWork(ipGeoVec+j-1) = 5
      Do i = 1,NumOfAt
      If ( Atom1 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j) = i
      Else If ( Atom2 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+1) = i
      Else If ( Atom3 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+2) = i
      Else If ( Atom4 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+3) = i
      End If
      End Do
      Call WordPos(k,OutLine,iStart,iStop)
      xvec(iInt) = StrToDble(OutLine(iStart:iStop))
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      DegOrRad = OutLine(iStart:iStop)
      If ( DegOrRad .eq. 'DEG' ) Then
      xvec(iInt) = xvec(iInt)*rpi/180.0d0
      End If
      j = j+5
      End If
      iInt = iInt+1
      Read(inpUnit,Format) InLine
      Call Normalize(InLine,OutLine)
      End Do
      iInt = iInt-1
!!    Call Int_to_Cart(GeoVec,xvec,AtCoord1,NumOfAt,iInt,Mass)
      l_a=NumOfAt
      l_n=max_len_xvec
      Call Int_to_Cart1(iWork(ipGeoVec),xvec,Work(ipAtCoord1),          &
     &    l_a,l_n)

      Call GetMem('GeoVec','Free','Inte',ipGeoVec,5*(3*NumOfAt-5))
!!
      Call KeyWord(inpUnit,'SECO',.false.,exist)
!!D Write(6,*)' READINP: Read internal coordinates for second state.'
      n = 5*(3*NumOfAt-5)
      Call GetMem('GeoVec','Allo','Inte',ipGeoVec,5*(3*NumOfAt-5))
      iInt = 1
      j = 1
      Read(inpUnit,Format) InLine
      Call Normalize(InLine,OutLine)
      Do While ( OutLine(1:4) .ne. 'END ')
      !---- Bond distance.
      l = Index(OutLine,'BOND')
      If ( l .ne. 0 ) Then
      k = l+Len('BOND')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      iWork(ipGeoVec+j-1) = 1
      Do i = 1,NumOfAt
      If ( Atom1 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j) = i
      Else If ( Atom2 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+1) = i
      End If
      End Do
      Call WordPos(k,OutLine,iStart,iStop)
      xvec(iInt) = StrToDble(OutLine(iStart:iStop))
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      AAOrAu = OutLine(iStart:iStop)
      If ( AAOrAu .eq. 'AA' ) Then
      xvec(iInt) = xvec(iInt)/Angstrom
      End If
      j = j+3
      End If
      !---- Valence Angle.
      l = Index(OutLine,'ANGLE')
      If ( l .ne. 0 ) Then
      k = l+Len('ANGLE')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom3 = OutLine(iStart:iStop)
      iWork(ipGeoVec+j-1) = 2
      Do i = 1,NumOfAt
      If ( Atom1 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j) = i
      Else If ( Atom2 .eq. AtomLbl(i) ) Then
      iWOrk(ipGeoVec+j+1) = i
      Else If ( Atom3 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+2) = i
      End If
      End Do
      Call WordPos(k,OutLine,iStart,iStop)
      xvec(iInt) = StrToDble(OutLine(iStart:iStop))
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      DegOrRad = OutLine(iStart:iStop)
      If ( DegOrRad .eq. 'DEG' ) Then
      xvec(iInt) = xvec(iInt)*rpi/180.0d0
      End If
      j = j+4
      End If
      !---- Linear Valence Angle.
      l = Index(OutLine,'LINANG')
      If ( l .ne. 0 ) Then
      k = l+Len('LINANG')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      Atom3 = OutLine(iStart:iStop)
      iWork(ipGeoVec+j-1) = 3
      Do i = 1,NumOfAt
      If ( Atom1 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j) = i
      Else If ( Atom2 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+1) = i
      Else If ( Atom3 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+2) = i
      End If
      End Do
      Call WordPos(k,OutLine,iStart,iStop)
      xvec(iInt) = StrToDble(OutLine(iStart:iStop))
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      DegOrRad = OutLine(iStart:iStop)
      If ( DegOrRad .eq. 'DEG' ) Then
      xvec(iInt) = xvec(iInt)*rpi/180.0d0
      End If
      j = j+4
      End If
      !---- Dihedral angle.
      l = Index(OutLine,'TORSION')
      If ( l .ne. 0 ) Then
      k = l+Len('TORSION')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom3 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom4 = OutLine(iStart:iStop)
      iWork(ipGeoVec+j-1) = 4
      Do i = 1,NumOfAt
      If ( Atom1 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j) = i
      Else If ( Atom2 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+1) = i
      Else If ( Atom3 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+2) = i
      Else If ( Atom4 .eq. AtomLbl(i) ) Then
      iWOrk(ipGeoVec+j+3) = i
      End If
      End Do
      Call WordPos(k,OutLine,iStart,iStop)
      xvec(iInt) = StrToDble(OutLine(iStart:iStop))
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      DegOrRad = OutLine(iStart:iStop)
      If ( DegOrRad .eq. 'DEG' ) Then
      xvec(iInt) = xvec(iInt)*rpi/180.0d0
      End If
      j = j+5
      End If
      !---- Out of Plane Angle.
      l = Index(OutLine,'OUTOFPL')
      If ( l .ne. 0 ) Then
      k = l+Len('OUTOFPL')
      Call WordPos(k,OutLine,iStart,iStop)
      Atom1 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom2 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom3 = OutLine(iStart:iStop)
      Call WordPos(k,OutLine,iStart,iStop)
      Atom4 = OutLine(iStart:iStop)
      iWOrk(ipGeoVec+j-1) = 5
      Do i = 1,NumOfAt
      If ( Atom1 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j) = i
      Else If ( Atom2 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+1) = i
      Else If ( Atom3 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+2) = i
      Else If ( Atom4 .eq. AtomLbl(i) ) Then
      iWork(ipGeoVec+j+3) = i
      End If
      End Do
      Call WordPos(k,OutLine,iStart,iStop)
      xvec(iInt) = StrToDble(OutLine(iStart:iStop))
      k = iStop+1
      Call WordPos(k,OutLine,iStart,iStop)
      DegOrRad = OutLine(iStart:iStop)
      If ( DegOrRad .eq. 'DEG' ) Then
      xvec(iInt) = xvec(iInt)*rpi/180.0d0
      End If
      j = j+5
      End If
      iInt = iInt+1
      Read(inpUnit,Format) InLine
      Call Normalize(InLine,OutLine)
      End Do
      iInt = iInt-1
!!    Call Int_to_Cart(GeoVec,xvec,AtCoord2,NumOfAt,iInt,Mass)
      l_a=NumOfAt
      l_n=max_len_xvec
      Call Int_to_Cart1(iWork(ipGeoVec),xvec,Work(ipAtCoord2),          &
     &    l_a,l_n)
      Call GetMem('GeoVec','Free','Inte',ipGeoVec,5*(3*NumOfAt-5))
      End If

      If ( inpunit1 .eq. 31 ) Then
        close(inpUnit1)
        close(inpUnit2)
      End If
!
! --- GEOMetries: End --------------------------------------------------

! ----------------------------------------------------------------------
! --- MXORder: Read maximum order for the transition dipole moment.
!
      Call KeyWord(inpUnit,'MXOR',.true.,exist)
      If (.NOT.exist) then
        If (.NOT.lISC) then
          Write(6,*) ' *** Warning: MXORder not specified !'
          Write(6,*) '     Transition Dipole is assumed as constant.'
          Write(6,*)
        EndIf
        max_dip = 0
      else
        Read(inpUnit,*) max_dip
      EndIf
!
! --- MXORder: End -----------------------------------------------------

! ----------------------------------------------------------------------
! --- OSCStr: Determines if you get the oscillator strength
!     or the transition intensity in the output.
!
      OscStr=.False.
      If (.NOT.lISC) then
        Call KeyWord(inpUnit,'OSCS',.true.,OscStr)
        If (OscStr) then
          Write(6,*) '     The Oscillator Strength will be calculated.'
          Write(6,*)
        else
          Write(6,*) '     The Intensity will be calculated.'
          Write(6,*)
        EndIf
      EndIf
!
! --- OSCStr: End ------------------------------------------------------

! ----------------------------------------------------------------------
! --- NANOmeters: Generate plot file in nanometers.
      Use_nm=.False.
      If (.NOT.lISC) Call KeyWord(inpUnit,'NANO',.true.,Use_nm)
!
! ----------------------------------------------------------------------
! --- CM-1: Generate plot file in cm-1.
      Use_cm=.False.
      If (.NOT.lISC) Call KeyWord(inpUnit,'CM-1',.true.,Use_cm)
!
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! --- PLOT: Enter the limits (in eV, cm-1 or in nm) for the plot file.
!
      cmstart=0.0d0
      cmend  =0.0d0
      plotwindow=.False.
      If (.NOT.lISC) then
        Call KeyWord(inpUnit,'PLOT',.true.,plotwindow)
        If(plotwindow) then
          read(inpunit,*) cmstart,cmend
        end if
      EndIf
!
! --- PLOT : End -------------------------------------------------------

      PlUnit='      '
      broadplot=.False.
      LifeTime = 130.0d-15
      If (lISC) GoTo 900
      If (Use_nm .AND. Use_cm) then
        Write(6,*)
        Write(6,*) ' ***************** ERROR *******************'
        Write(6,*) ' NANOmeters and CM-1 are mutually exclusive.'
        Write(6,*) ' *******************************************'
        Call Quit_OnUserError()
      EndIf
      If (Use_nm) then
        PlUnit=' nm.  '
        Write(6,*) '     The plot file will be in nanometers and the '
      else If (Use_cm) then
        PlUnit=' cm-1.'
        Write(6,*) '     The plot file will be in cm-1 and the'
      else
        PlUnit=' eV.  '
        Write(6,*) '     The plot file will be in eV and the'
      EndIf
      If(plotwindow) then
        Write(6,*) '     interval is from ',cmstart,' to ',cmend,PlUnit
      else
        Write(6,*) '     interval automatically defined.'
      EndIf
      Write(6,*)
! ----------------------------------------------------------------------
! --- BROAplot: Gives the peaks in the spectrum an artificial halfwidth.
!
      Call KeyWord(inpUnit,'BROA',.true.,broadplot)
      Call KeyWord(inpUnit,'LIFE',.true.,exist)
      If (exist) Read(inpunit,*) LifeTime
      Write(6,*) '     Life time : ',LifeTime,' sec'
 900  FWHM = hbarcm/(2.0d0*LifeTime)
      If (.NOT.lISC) Write(6,'(1X,A,F8.1,A,F9.6,A)')                    &
     &    '     FWHM : ',FWHM,' cm-1 /',FWHM/8065.5D0,' eV'

!
! --- BROAplot: End ----------------------------------------------------

! ----------------------------------------------------------------------
! --- VIBWrite: Writes vibrational levels to logfile.
!
      WriteVibLevels=.False.
      Call KeyWord(inpUnit,'VIBW',.true.,WriteVibLevels)
!
! ----------------------------------------------------------------------
! --- VibPlot: Check for keyword VibPlot.
!
      VibModPlot=.False.
      Call KeyWord(inpUnit,'VIBP',.true.,VibModPlot)
      If ( VibModPlot ) Then
        Call KeyWord(inpUnit,'CYCL',.false.,exist)
        If ( exist ) Then
          Call KeyWord(inpUnit,'VIBP',.true.,exist)
          Read(inpUnit,Format) InLine
          Call Normalize(InLine,OutLine)
          k = 1
          Call WordPos(k,OutLine,iStart,iStop)
          Call WordPos(k,OutLine,iStart,iStop)
          Bond(nBond+1) = iStrToInt(OutLine(iStart:iStop))
          Call WordPos(k,OutLine,iStart,iStop)
          Bond(nBond+2) = iStrToInt(OutLine(iStart:iStop))
          nBond = nBond+2
        End If
      End If
!
! --- VibPlot: End -----------------------------------------------------

! ----------------------------------------------------------------------

! --- HUGElog: If exists, then the logfile will be more detailed.
!
      Huge_Print=.False.
      Call KeyWord(inpUnit,'HUGE',.true.,Huge_Print)
!
! --- HUGElog: End -----------------------------------------------------

! ----------------------------------------------------------------------
! --- EXPAnsion: Program will be aborted after expansion point geometry
!                is calculated.
      lExpan=.False.
      Call KeyWord(inpUnit,'EXPA',.true.,lExpan)
!
! --- HUGElog: End -----------------------------------------------------

! ----------------------------------------------------------------------
! --- FORCe: Here really read force constants.
      Call KeyWord(inpUnit,'FORC',.true.,exist)
      max_term = 2
      Call KeyWord(inpUnit,'FIRS',.false.,exist)
!!D Write(6,*)' READINP: Read force constants for first state.'
      Read(InpUnit,Format) InLine
      Call Normalize(InLine,Outline)
      iend=Index(OutLine,' ')
      CoordType=OutLine(1:iend-1)
!!D Write(6,*)' READINP: CoordType=',CoordType
!!
!!---- Read Hessian. If Hessian is given in cartesian coordinates,
!!     first remove total translation and total rotation and then
!!     transform it to internal coordinates.
      If ( CoordType .eq. 'FILE' ) Then
        Coordtype = 'CARTESIAN'
        inpUnit1 = 31
        call molcas_open(31,'UNSYM1')
!        Open (unit=31,file='UNSYM1')
        Call KeyWord(inpUnit1,'UNSYMMETRIZED HESSIAN',.false.,exist)
        Read(inpUnit1,'(a17)') Inline
        Read(inpUnit1,'(a17)') Inline
      Else
        inpUnit1 = inpUnit
      End If
!!
      If ( CoordType .eq. 'INTERNAL' ) Then
        Do i = 1,NumInt
          Read(inpUnit,*)                                               &
     &       (Work(ipHess1+i+l_Hess1*(j-1)-1),j=1,NumInt)
        End Do
      Else If ( CoordType .eq. 'CARTESIAN' ) Then
        Call GetMem('SS','Allo','Real',ipSS,3*NumOfAt*NumInt)
        call dcopy_(3*NumOfAt*NumInt,[0.0d0],0,Work(ipSS),1)
!          SS = 0.0d0
        Call CalcS(Work(ipAtCoord1),InterVec,Work(ipSS),NumInt,NumOfAt)
!!---- Invert S matrix and remove total translation and
!!     total rotation.
        Call GetMem('Sinv','Allo','Real',ipSinv,3*NumOfAt*NumInt)
        Call RotTranRem(Work(ipSinv),Work(ipSS),Mass,                   &
     &    Work(ipAtCoord1),NumOfAt,NumInt)
!!
!!---- Read cartesian force constant matrix.
        n = 3*NumOfAt
        nFcart=n
        Call GetMem('Fcart','Allo','Real',ipFcart,nFcart*nFcart)
        Do m = 1,(3*NumOfAt)
          Read(inpUnit1,'(a17)',err=999) Inline
          Read(inpUnit1,*,err=999)                                      &
     &       (Work(ipFcart+m+nFcart*(n-1)-1),                           &
     &                 n=1,(3*NumOfAt))
          GoTo 998
  999     Continue
          Write(6,*)
          Write(6,*) ' ***************** ERROR *******************'
          Write(6,*) ' Error trying to read cartesian force cnsts.'
          Write(6,*) ' for the first state. Probably wrong input  '
          Write(6,*) ' structure, perhaps too few values.         '
          Write(6,*) ' *******************************************'
          Call Quit_OnUserError()
  998     Continue
        End Do
!!
!!---- Check if Hessian is symmetric.
        error = 0.0d0
        Do m = 1,(3*NumOfAt)
          Do n = 1,(3*NumOfAt)
            error = error+(Work(ipFcart+m+nFcart*(n-1)-1)-              &
     &          Work(ipFcart+n+nFcart*(m-1)-1))**2
          End Do
        End Do
        If ( error .gt. 1.0d-10) Then
          Write(6,*)
          Write(6,*)' ******************** ERROR **********************'
          Write(6,*) ' Non-symmetric error in cartesian Hessian:',error
          Write(6,*)' *************************************************'
          Call Quit_OnUserError()
        End If
!!
!!---- Transform cartesian F matrix to internal coordinates.
!!PAM01: Here follows an n**4-scaling 2-index transformation
!!It has been replaced by two n**3 matrix multiplies
!!    Do j = 1,NumInt
!!       Do i = 1,NumInt
!!          Fsum = 0.0d0
!!          Do n = 1,(3*NumOfAt)
!!             Do m = 1,(3*NumOfAt)
!!                Fsum = Fsum+Fcart(m,n)* &
!!                   Sinv(1+Mod((m+2),3),Int((m+2)/3),i)*  &
!!                   Sinv(1+Mod((n+2),3),Int((n+2)/3),j)
!!             End Do
!!          End Do
!!          Hess1(i,j) = Fsum
!!       End Do
!!    End Do
!!PAM01 Here follows replacement code:
        Call GetMem('Temp','Allo','Real',ipTemp,3*NumOfAt*NumInt)
        Call DGEMM_('N','N',                                            &
     &            3*NumOfAt,NumInt,3*NumOfAt,                           &
     &            1.0d0,Work(ipFCart),3*NumOfAt,                        &
     &            work(ipSInv),3*NumOfAt,                               &
     &            0.0d0,Work(ipTemp),3*NumOfAt)
        Call DGEMM_('T','N',                                            &
     &            NumInt,NumInt,3*NumOfAt,                              &
     &            1.0d0,Work(ipSInv),3*NumOfAt,                         &
     &            Work(ipTemp),3*NumOfAt,                               &
     &            0.0d0,Work(ipHess1),NumInt)
        Call GetMem('Temp','Free','Real',ipTemp,3*NumOfAt*NumInt)
        Call GetMem('Sinv','Free','Real',ipSinv,3*NumOfAt*NumInt)
        Call GetMem('SS','Free','Real',ipSS,3*NumOfAt*NumInt)
        call GetMem('Fcart','Free','Real',ipFcart,nFcart*nFcart)

      End If
!!
!!---- Scale Hessian if scaling factors were given.
      Call KeyWord(inpUnit,'SCAL',.true.,exist)
      If ( exist ) Then
      Call GetMem('Scale1','Allo','Real',ipScaleParam1,NumInt)

      Call KeyWord(inpUnit,'FIRS',.false.,exist)
      Do i = 1,NumInt
      Read(inpUnit,*) Work(ipScaleParam1+i-1)
      Work(ipScaleParam1+i-1) = sqrt(Work(ipScaleParam1+i-1))
      End Do
      Do j = 1,Numint
      Do i = 1,NumInt
      Work(ipHess1+i+l_Hess1*(j-1)-1) =                                 &
     &          Work(ipHess1+i+l_Hess1*(j-1)-1)*                        &
     &          Work(ipScaleParam1+i-1)*                                &
     &          Work(ipScaleParam1+j-1)
      End Do
      End Do
      Call GetMem('Scale1','Free','Real',ipScaleParam1,NumInt)
      End if
!!D Write(6,*)' Scaled.'
!!
!!---- Hessian for second surface.
      Call KeyWord(inpUnit,'FORC',.true.,exist)
      Call KeyWord(inpUnit,'SECO',.false.,exist)
!!D Write(6,*)' READINP: Read force constants for second state.'
!      Read(inpUnit,*) CoordType
!      Call UpCase(CoordType)
! Replace above two lines with:
      Read(InpUnit,Format) InLine
      Call Normalize(InLine,Outline)
      iend=Index(OutLine,' ')
      CoordType=OutLine(1:iend-1)
! End of replacement
      If ( CoordType .eq. 'FILE' ) Then
        Coordtype = 'CARTESIAN'
        inpUnit2 = 32
        Call molcas_open(32,'UNSYM2')
!      Open (unit=32,file='UNSYM2')
        Call KeyWord(inpUnit2,'UNSYMMETRIZED HESSIAN',.false.,exist)
        Read(inpUnit2,'(a17)') Inline
        Read(inpUnit2,'(a17)') Inline
      Else
        inpUnit2 = inpUnit
      End If
!!
      If ( CoordType .eq. 'INTERNAL' ) Then
        Do i = 1,NumInt
        Read(inpUnit,*)                                                 &
     &       (Work(ipHess2+i+l_Hess1*(j-1)-1),j=1,NumInt)
        End Do
      Else If ( CoordType .eq. 'CARTESIAN' ) Then
        Call GetMem('SS','Allo','Real',ipSS,3*NumOfAt*NumInt)

        call dcopy_(3*NumOfAt*NumInt,[0.0d0],0,Work(ipSS),1)
!          SS = 0.0d0
        Call CalcS(Work(ipAtCoord2),InterVec,Work(ipSS),NumInt,NumofAt)
!!
!!---- Invert S matrix and remove total translation and
!!     total rotation.
        Call GetMem('Sinv','Allo','Real',ipSinv,3*NumOfAt*NumInt)

        Call RotTranRem(Work(ipSinv),Work(ipSS),Mass,                   &
     &    Work(ipAtCoord2),NumOfAt,NumInt)
!!
!!---- Read cartesian force constant matrix.
        n = 3*NumOfAt
        nFcart=n
        call GetMem('Fcart','Allo','Real',ipFcart,nFcart*nFcart)

        Do m = 1,(3*NumOfAt)
          Read(inpUnit2,'(a17)',err=997) Inline
          Read(inpUnit2,*,err=997) (Work(ipFcart+m+nFcart*(n-1)-1),     &
     &              n=1,(3*NumOfAt))
          GoTo 996
 997      Continue
          Write(6,*)
          Write(6,*) ' ***************** ERROR *******************'
          Write(6,*) ' Error trying to read cartesian force cnsts.'
          Write(6,*) ' for the second state. Probably wrong input '
          Write(6,*) ' structure, perhaps too few values.         '
          Write(6,*) ' *******************************************'
          Call Quit_OnUserError()
 996      Continue
        End Do
!!
!!---- Check if Hessian is symmetric.
        error = 0.0d0
        Do n = 1,(3*NumOfAt)
          Do m = 1,(3*NumOfAt)
          error = error+(Work(ipFcart+m+nFcart*(n-1)-1)-                &
     &          Work(ipFcart+n+nFcart*(m-1)-1))**2
          End Do
        End Do
        If ( error .gt. 1.0d-2 ) Then
          Write(6,*)
          Write(6,*)' ******************** ERROR **********************'
          Write(6,*)' Non-symmetric error in cartesian Hessian:',error
          Write(6,*)' *************************************************'
          Call Quit_OnUserError()
        End If
!!
!!---- Transform cartesian F matrix to internal coordinates.
!!PAM01: Here follows an n**4-scaling 2-index transformation
!!It has been replaced by two n**3 matrix multiplies
!!    Do j = 1,NumInt
!!       Do i = 1,NumInt
!!          Fsum = 0.0d0
!!          Do n = 1,(3*NumOfAt)
!!             Do m = 1,(3*NumOfAt)
!!                Fsum = Fsum+Fcart(m,n)* &
!!                    Sinv(1+Mod((m+2),3),Int((m+2)/3),i)*  &
!!                    Sinv(1+Mod((n+2),3),Int((n+2)/3),j)
!!             End Do
!!          End Do
!!          Hess2(i,j) = Fsum
!!       End Do
!!    End Do
!!PAM01 Here follows replacement code:
      Call GetMem('Temp','Allo','Real',ipTemp,3*NumOfAt*NumInt)
      Call DGEMM_('N','N',                                              &
     &            3*NumOfAt,NumInt,3*NumOfAt,                           &
     &            1.0d0,Work(ipFCart),3*NumOfAt,                        &
     &            work(ipSInv),3*NumOfAt,                               &
     &            0.0d0,Work(ipTemp),3*NumOfAt)
      Call DGEMM_('T','N',                                              &
     &            NumInt,NumInt,3*NumOfAt,                              &
     &            1.0d0,Work(ipSInv),3*NumOfAt,                         &
     &            Work(ipTemp),3*NumOfAt,                               &
     &            0.0d0,Work(ipHess2),NumInt)
      Call GetMem('Temp','Free','Real',ipTemp,3*NumOfAt*NumInt)
      Call GetMem('Sinv','Free','Real',ipSinv,3*NumOfAt*NumInt)
      Call GetMem('SS','Free','Real',ipSS,3*NumOfAt*NumInt)
      call GetMem('Fcart','Free','Real',ipFcart,nFcart*nFcart)

      End If
!!
      If ( inpunit1 .eq. 31) Then
        Close(inpunit1)
        Close(inpunit2)
      End If
!
! --- FORCe: End -------------------------------------------------------

! ----------------------------------------------------------------------
! --- SCALe: Scale Hessian if scaling factors were given.
      Call KeyWord(inpUnit,'SCAL',.true.,exist)
      If ( exist ) Then
        Call GetMem('Scale2','Allo','Real',ipScaleParam2,NumInt)

        Call KeyWord(inpUnit,'SECO',.false.,exist)
        Do i = 1,NumInt
          Read(inpUnit,*) Work(ipScaleParam2+i-1)
          Work(ipScaleParam2+i-1) = sqrt(Work(ipScaleParam2+i-1))
        End Do
        Do j = 1,NumInt
          Do i = 1,NumInt
            Work(ipHess2+i+l_Hess1*(j-1)-1) =                           &
     &      Work(ipHess2+i+l_Hess1*(j-1)-1)*                            &
     &      Work(ipScaleParam2+i-1)*                                    &
     &      Work(ipScaleParam2+j-1)
          End Do
        End Do
        Call GetMem('Scale2','Allo','Real',ipScaleParam2,NumInt)
      End If
!
! --- SCALe: End -------------------------------------------------------

! ----------------------------------------------------------------------
! --- DIPOles: Read Transition Dipoles or Spin-Orbit Coupling.
      Call KeyWord(inpUnit,'DIPO',.true.,exist)
      If (.NOT.exist) Call KeyWord(inpUnit,'SOC ',.true.,exist)
      If(.not.exist) then
        Write(6,*)
        Write(6,*) ' ************ ERROR ****************'
        Write(6,*) ' The DIPOLES/SOC keyword is missing.'
        Write(6,*) ' ***********************************'
        Call Quit_OnUserError()
      End if
!!D Write(6,*)' READINP: Read dipoles (''DIPOLES'').'
      Read(inpUnit,Format) Inline
      Call Normalize(InLine,OutLine)
      If ( OutLine(1:4) .eq. 'FILE' ) Then
!!D  Write(6,*)' Read trans dips from file UNSYM21'
        call molcas_open(31,'UNSYM21')
!      Open (unit=31,file='UNSYM21')
        If ( max_dip .eq. 1 ) then
          Call KeyWord(31,'*BEGIN TRANSDIPDER FOR COMPONENT X',         &
     &      .false.,exist)
          Read(31,*) (Work(ipTranDipGrad+3*(i-1)),i=1,3*NumOfAt)
          Call KeyWord(31,'*BEGIN TRANSDIPDER FOR COMPONENT Y',         &
     &      .false.,exist)
          Read(31,*) (Work(ipTranDipGrad+1+3*(i-1)),i=1,3*NumOfAt)
          Call KeyWord(31,'*BEGIN TRANSDIPDER FOR COMPONENT Z',         &
     &      .false.,exist)
          Read(31,*) (Work(ipTranDipGrad+2+3*(i-1)),i=1,3*NumOfAt)
        End If
        XnotFound=.true.
        YnotFound=.true.
        ZnotFound=.true.
        Call KeyWord(31,'*BEGIN TRANSITION PROPERTIES',.false.,exist)
        Do While(XnotFound.or.YnotFound.or.ZnotFound)
          Read(31,'(A37)') InLine
          If ( InLine(1:20) .eq. 'MLTPL  1 COMPONENT1 ' ) Then
            XnotFound = .false.
            Read(InLine(22:37),*) TranDip(1)
          End If
          If ( InLine(1:20) .eq. 'MLTPL  1 COMPONENT2 ' ) Then
            YnotFound = .false.
            Read(InLine(22:37),*) TranDip(2)
          End If
          If ( InLine(1:20) .eq. 'MLTPL  1 COMPONENT3 ' ) Then
            ZnotFound = .false.
            Read(InLine(22:37),*) TranDip(3)
          End If
        End Do
        Close(31)
      Else
        Call KeyWord(inpUnit,'DIPO',.true.,exist)
!!D Write(6,*)' READINP: Read dipoles (''DIPOLES'').'
        Read(inpUnit,*) (TranDip(i),i=1,3)
        If ( max_dip .eq. 1 ) Then
          Read(inpUnit,*)                                               &
     &       (Work(ipTranDipGrad+3*(i-1)),i=1,3*NumOfAt)
          Read(inpUnit,*)                                               &
     &       (Work(ipTranDipGrad+1+3*(i-1)),i=1,3*NumOfAt)
          Read(inpUnit,*)                                               &
     &       (Work(ipTranDipGrad+2+3*(i-1)),i=1,3*NumOfAt)
        End If
      End If

      Goto 20
!
!                         End of Forcefield part
! ----------------------------------------------------------------------

 10    Continue

!!-----------------------------------------------------------------------!
!!
!!                       Energy surface input part
!!
!!-----------------------------------------------------------------------!
!!
      ForceField = .false.
!!
!!---- Read transformations.
      Call KeyWord(inpUnit,'NONL',.true.,exist)
      If(.not.Exist) Then
        Write(6,*)
        Write(6,*) ' *********** ERROR **************'
        Write(6,*) ' Cannot find keyword  NONLINEAR !'
        Write(6,*) ' ********************************'
        Call Quit_OnUserError()
      End If
      Do icoord = 1,NumInt
        Read(inpUnit,Format) InLine
        Call Normalize(InLine,OutLine)
        trfName1(icoord) = OutLine
      End Do
      Do icoord = 1,NumInt
        Read(inpUnit,Format) InLine
        Call Normalize(InLine,OutLine)
        trfName2(icoord) = OutLine
      End Do
!!
!!---- Read terms in polynomial.
      Call KeyWord(inpUnit,'POLY',.true.,exist)
      If(.not.Exist) Then
        Write(6,*)
        Write(6,*) ' ************ ERROR **************'
        Write(6,*) ' Cannot find keyword  POLYNOMIAL !'
        Write(6,*) ' *********************************'
        Call Quit_OnUserError()
      End If
!!D Write(6,*)' READINP: Read polynomial.'
      k = 1
      nvar = 0
      Read(inpUnit,Format) InLine
      Call WordPos(k,InLine,iStart,iStop)
      Do While ( iStop .lt. Len(Inline) )
        k = iStop+1
        nvar = nvar+1
        Call WordPos(k,InLine,iStart,iStop)
      End Do
      nPolyTerm = 1
      Read(inpUnit,Format) InLine
      Call Normalize(InLine,OutLine)
      Do While ( OutLine(1:4) .ne. 'END ')
        nPolyTerm = nPolyTerm+1
        Read(inpUnit,Format) InLine
        Call Normalize(InLine,OutLine)
      End Do
      Call KeyWord(inpUnit,'POLY',.true.,exist)
      If(.not.Exist) Then
        Write(6,*)
        Write(6,*) ' ************ ERROR **************'
        Write(6,*) ' Cannot find keyword  POLYNOMIAL !'
        Write(6,*) ' *********************************'
        Call Quit_OnUserError()
      End If

      Call GetMem('ipow','Allo','Inte',ipipow,nPolyTerm*nvar)
      Do iterm = 1,nPolyTerm
      Read(inpUnit,*)                                                   &
     &    (iWork(ipipow+iterm+nPolyTerm*(ivar-1)-1),ivar=1,nvar)
      End Do
!!
!!---- Find highest power of term in polynomial.
      max_term = 0
      Do iterm = 1,nPolyTerm
        nsum = 0
        Do jvar = 1,nvar
          nsum = nsum+iWork(ipipow+iterm+nPolyTerm*(jvar-1)-1)
          If ( nsum .gt. max_term ) max_term = nsum
        End Do
      End Do
!!
!!---- Read grid points and energies.
      CoordType = 'INTERNAL'
      Call KeyWord(inpUnit,'DATA',.true.,exist)
      If(.not.Exist) Then
        Write(6,*)
        Write(6,*) ' ********* ERROR ***********'
        Write(6,*) ' Cannot find keyword  DATA !'
        Write(6,*) ' ***************************'
        Call Quit_OnUserError()
      End If
!!D Write(6,*)' READINP: Read grid points (Key word ''DATA'').'
      idata = 0
      Do while(.true.)
      idata = idata+1
      Read(inpUnit,*,end=30) (coord(idata,jvar),jvar=1,nvar),           &
     &   (GrdVal(idata,j),j=1,10)
      End Do
 30    Continue
!!
      ndata = idata-1
      Call GetMem('var','Allo','Real',ipvar,ndata*nvar)
      Do idata = 1,ndata
      Do jvar = 1,nvar
      Work(ipvar+idata+ndata*(jvar-1)-1) = coord(idata,jvar)
      End Do
      End Do
!!
!!---- CASPT2 energies for ground state.
      Call GetMem('yin1','Allo','Real',ipyin1,ndata)

      Do idata = 1,ndata
      Work(ipyin1+idata-1) = GrdVal(idata,2)
      End Do
!!
!!---- CASPT2 energies for excited state.
      Call GetMem('yin2','Allo','Real',ipyin2,ndata)
      Do idata = 1,ndata
      Work(ipyin2+idata-1) = GrdVal(idata,6)
      End Do
!!
!!---- If three atomic molecule, add Watson correction to energy values.
      If ( NumOfAt .eq. 3 ) Then
      m1 = Mass(1)*uToAu
      m2 = Mass(2)*uToAu
      m3 = Mass(3)*uToAu
      m12 = (m1*m2)/(m1+m2)
      m23 = (m2*m3)/(m2+m3)
      Do i = 1,ndata
      r1 = Work(ipvar+i-1)
      r2 = Work(ipvar+i+ndata -1)
      theta = Work(ipvar+i+ndata*2 -1)*(rpi/180.0d0)
      const = -(1.0d0/8.0d0)*(1.0d0/(m12*r1**2)+1.0d0/(m23*r2**2))*     &
     &       (1.0d0+(1.0d0/(sin(theta)**2)))
      const = const+0.25d0*(cos(theta)/(m2*r1*r2))*                     &
     &       (cos(theta)/sin(theta))**2
      Work(ipyin1+i-1) = Work(ipyin1+i-1)+const
      Work(ipyin2+i-1) = Work(ipyin2+i-1)+const
      End Do
      End If
!!
      Call GetMem('t_dipin1','Allo','Real',ipt_dipin1,ndata)
      Do idata = 1,ndata
      Work(ipt_dipin1+idata-1) = GrdVal(idata,9)
      End Do
!!
      Call GetMem('t_dipin2','Allo','Real',ipt_dipin2,ndata)
      Do idata = 1,ndata
      Work(ipt_dipin2+idata-1) = GrdVal(idata,10)
      End Do
      Call GetMem('t_dipin3','Allo','Real',ipt_dipin3,ndata)
      call dcopy_(ndata,[0.0d0],0,Work(ipt_dipin3),1)
!       t_dipin3 = 0.0d0
!!
!!---- Make sure that all transition dipole values have the same sign.
      sign1_1 = Work(ipt_dipin1)/abs(Work(ipt_dipin1))
      sign2_1 = Work(ipt_dipin2)/abs(Work(ipt_dipin2))
      Do idata = 2,ndata
      sign1_2 =                                                         &
     &    Work(ipt_dipin1+idata-1)/abs(Work(ipt_dipin1+idata-1))
      sign2_2 =                                                         &
     &    Work(ipt_dipin2+idata-1)/abs(Work(ipt_dipin2+idata-1))
      If ( int(sign1_1*sign2_1) .ne. int(sign1_2*sign2_2) ) Then
      Write(6,*) 'Sign shift in transition dipole data:',idata
      End If
      Work(ipt_dipin1+idata-1) =                                        &
     &    sign1_1*sign1_2*Work(ipt_dipin1+idata-1)
      Work(ipt_dipin2+idata-1) =                                        &
     &    sign2_1*sign2_2*Work(ipt_dipin2+idata-1)
      End Do
!!D Write(6,*)' READINP: Transition dipoles finished.'

 20     Continue
!
!                         End of Energy surface input part
! ----------------------------------------------------------------------

!!D Write(6,*)' Ending READINP.'
!!
      End
