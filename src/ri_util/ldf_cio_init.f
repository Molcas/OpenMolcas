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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_CIO_Init(BufFrac,irc)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Initialize coefficient I/O.
C
C     BufFrac, which must be between 0.0d0 and 1.0d0, specifies the
C     fraction of available memory that will at most be allocated as
C     coefficient buffer. The buffer will never exceed the max needed
C     memory.
C
      Implicit None
      Real*8 BufFrac
      Integer irc
#include "WrkSpc.fh"
#include "ldf_cio.fh"
#include "ldf_atom_pair_info.fh"

#if defined (_DEBUG_)
      Character*12 SecNam
      Parameter (SecNam='LDF_CIO_Init')
      Real*8 Byte
      Character*2 Unt
#endif

      Integer  LDF_X_OpenCoefficientFile
      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_X_OpenCoefficientFile
      External LDF_nBas_Atom, LDF_nBasAux_Pair

      Logical CFileExists

      Real*8 Frac

      Integer ip, l
      Integer lBuf_Max
      Integer iAtomPair
      Integer iAtom, jAtom
      Integer iAddr

      Integer i, j
      Integer AP_DiskC
      Integer AP_Atoms
      AP_DiskC(i)=iWork(ip_AP_DiskC-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Init return code
      irc=0

      ! Init variables in ldf_cio.fh
      Lu_LDFC=0
      ip_LDFC_Buffer=0
      l_LDFC_Buffer=0
      ip_LDFC_Blocks=0
      l_LDFC_Blocks=0
      LastAtomPair=0

      ! Check that the coefficient file exists
      Call f_Inquire('LDFC',CFileExists)
      If (.not.CFileExists) Then
         irc=-1
#if defined (_DEBUG_)
         Write(6,'(/,A,A)') SecNam,': Coefficient file not found!!'
         Call xFlush(6)
#endif
         Return
      End If

      ! Open coefficient file
      Lu_LDFC=LDF_X_OpenCoefficientFile()

      ! Normalize BufFrac
      Frac=min(max(BufFrac,0.0d0),1.0d0)

      ! Return if no buffer requested
      If (Frac.lt.1.0d-14) Then
#if defined (_DEBUG_)
         ! Set lBuf_Max for the sake of printing
         lBuf_Max=0
#endif
         Go To 1  ! return
      End If

      ! Allocate AP block pointer array
      Call GetMem('MaxiMem','Max ','Inte',ip,l)
      If (l.lt.NumberOfAtomPairs) Return ! insufficient memory
      l_LDFC_Blocks=NumberOfAtomPairs
      Call GetMem('LDFC_Blk','Allo','Inte',ip_LDFC_Blocks,l_LDFC_Blocks)

      ! Compute max buffer size
      Call GetMem('MaxMem','Max ','Real',ip,l)
      lBuf_Max=int(Frac*dble(l))
      If (lBuf_Max.lt.1) Then ! no memory available
         Call GetMem('LDFC_Blk','Free','Inte',
     &               ip_LDFC_Blocks,l_LDFC_Blocks)
         ip_LDFC_Blocks=0
         l_LDFC_Blocks=0
         Go To 1  ! return
      End If

      ! Compute buffer size
      l=0
      iAtomPair=0
      Do While (l.lt.lBuf_Max .and. iAtomPair.lt.NumberOfAtomPairs)
         iAtomPair=iAtomPair+1
         iAtom=AP_Atoms(1,iAtomPair)
         jAtom=AP_Atoms(2,iAtomPair)
         l=l
     &    +LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
     &    *LDF_nBasAux_Pair(iAtomPair)
         If (l.le.lBuf_Max) Then
            LastAtomPair=iAtomPair
            l_LDFC_Buffer=l
         End If
      End Do

      If (l_LDFC_Buffer.gt.0) Then
         ! Reallocate AP block pointer array
         If (LastAtomPair.lt.NumberOfAtomPairs) Then
            Call GetMem('LDFC_Blk','Free','Inte',
     &                  ip_LDFC_Blocks,l_LDFC_Blocks)
            l_LDFC_Blocks=LastAtomPair
            Call GetMem('LDFC_Blk','Allo','Inte',
     &                  ip_LDFC_Blocks,l_LDFC_Blocks)
         End If
         ! Allocate buffer
         Call GetMem('CBuffer','Allo','Real',
     &               ip_LDFC_Buffer,l_LDFC_Buffer)
         ! Read coefficients
         ip=ip_LDFC_Buffer
         Do iAtomPair=1,LastAtomPair
            iAtom=AP_Atoms(1,iAtomPair)
            jAtom=AP_Atoms(2,iAtomPair)
            l=LDF_nBas_Atom(iAtom)*LDF_nBas_Atom(jAtom)
     &       *LDF_nBasAux_Pair(iAtomPair)
            iAddr=AP_DiskC(iAtomPair)
            Call dDAFile(Lu_LDFC,2,Work(ip),l,iAddr)
            iWork(ip_LDFC_Blocks-1+iAtomPair)=ip
            ip=ip+l
         End Do
#if defined (_DEBUG_)
         If (ip.ne.(ip_LDFC_Buffer+l_LDFC_Buffer)) Then
            Call WarningMessage(2,SecNam//': buffer error!')
            Call LDF_Quit(1)
         End If
#endif
      End If

    1 Continue
#if defined (_DEBUG_)
      Write(6,'(/,A,A)')
     & SecNam,': coefficient I/O has been initialized!'
      Call Cho_Head('LDF coefficient buffer info','=',80,6)
      Write(6,'(A,F9.4)')
     & 'Max memory fraction (input).......',BufFrac
      Write(6,'(A,F9.4)')
     & 'Max memory fraction (normalized)..',Frac
      Call Cho_Word2Byte(lBuf_Max,8,Byte,Unt)
      Write(6,'(A,I9,A,F7.2,1X,A,A)')
     & 'Max possible buffer size..........',lBuf_Max,' words (',
     & Byte,Unt,')'
      Call Cho_Word2Byte(l_LDFC_Buffer,8,Byte,Unt)
      Write(6,'(A,I9,A,F7.2,1X,A,A)')
     & 'Allocated buffer..................',l_LDFC_Buffer,' words (',
     & Byte,Unt,')'
      Write(6,'(A,I9)')
     & 'Total number of atom pairs........',NumberOfAtomPairs
      If (NumberOfAtomPairs.gt.0) Then
         Write(6,'(A,I9,A,F7.2,A)')
     &   'Number of atom pairs in buffer....',LastAtomPair,
     &   ' (',1.0d2*dble(LastAtomPair)/dble(NumberOfAtomPairs),'%)'
      Else
         Write(6,'(A,I9)')
     &   'Number of atom pairs in buffer....',LastAtomPair
      End If
      If (LastAtomPair.gt.NumberOfAtomPairs) Then
         Call WarningMessage(2,
     &                     SecNam//': LastAtomPair > NumberOfAtomPairs')
         Call LDF_Quit(1)
      End If
      Write(6,*)
      Call xFlush(6)
#endif

      End
