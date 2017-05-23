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
      Subroutine LDF_PrintAtomPairDiagonal(iAtomPair)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Purpose: print information about the atom pair diagonal
C
      Implicit None
      Integer iAtomPair
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"

      Character*25 SecNam
      Parameter (SecNam='LDF_PrintAtomPairDiagonal')

      Integer  LDF_AtomPair_DiagDim
      External LDF_AtomPair_DiagDim

      real*8 ddot_

      Integer l
      Integer ipD, ipDB
      Integer nNegD, nNegDB

      Real*8 NormD, NormDB
      Real*8 SumD, SumDB
      Real*8 AvgD, AvgDB
      Real*8 StdDevD, StdDevDB
      Real*8 MinD, MinDB
      Real*8 MaxD, MaxDB

      Integer i, j
      Integer ip_DB, ip_D, AP_Atoms
      ip_DB(i)=iWork(ip_AP_DiagBak-1+i)
      ip_D(i)=iWork(ip_AP_Diag-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      l=LDF_AtomPair_DiagDim(iAtomPair)
      If (l.lt.1) Then
         Call WarningMessage(2,SecNam//': l < 1')
         Call LDF_Quit(1)
      End If

      ipDB=ip_DB(iAtomPair)
      NormDB=sqrt(dDot_(l,Work(ipDB),1,Work(ipDB),1))
      ipD=ip_D(iAtomPair)
      NormD=sqrt(dDot_(l,Work(ipD),1,Work(ipD),1))

      SumDB=Work(ipDB)
      SumD=Work(ipD)
      Do i=1,l-1
         SumDB=SumDB+Work(ipDB+i)
         SumD=SumD+Work(ipD+i)
      End Do
      AvgDB=SumDB/dble(l)
      AvgD=SumD/dble(l)

      StdDevDB=(Work(ipDB)-AvgDB)**2
      StdDevD=(Work(ipD)-AvgD)**2
      Do i=1,l-1
         StdDevDB=(Work(ipDB+i)-AvgDB)**2
         StdDevD=(Work(ipD+i)-AvgD)**2
      End Do
      StdDevDB=sqrt(StdDevDB/dble(l))
      StdDevD=sqrt(StdDevD/dble(l))

      MinDB=Work(ipDB)
      MinD=Work(ipD)
      MaxDB=Work(ipDB)
      MaxD=Work(ipD)
      Do i=1,l-1
         MinDB=min(MinDB,Work(ipDB+i))
         MinD=min(MinD,Work(ipD+i))
         MaxDB=max(MaxDB,Work(ipDB+i))
         MaxD=max(MaxD,Work(ipD+i))
      End Do

      nNegDB=0
      nNegD=0
      Do i=0,l-1
         If (Work(ipDB+i).lt.0.0d0) nNegDB=nNegDB+1
         If (Work(ipD+i).lt.0.0d0) nNegD=nNegD+1
      End Do

      Write(6,'(/,A,I10)') 'Atom Pair............',iAtomPair
      Write(6,'(A,2I10)')  'Atoms................',
     & AP_Atoms(1,iAtomPair),AP_Atoms(2,iAtomPair)
      Write(6,'(A,I10)')   'Diagonal dimension...',l
      Write(6,'(/,17X,A,10X,A)') 'Original','Current'
      Write(6,'(A,1P,2(1X,D16.6))') 'Norm    ',NormDB,NormD
      Write(6,'(A,1P,2(1X,D16.6))') 'Sum     ',SumDB,SumD
      Write(6,'(A,1P,2(1X,D16.6))') 'Average ',AvgDB,AvgD
      Write(6,'(A,1P,2(1X,D16.6))') 'Std Dev ',StdDevDB,StdDevD
      Write(6,'(A,1P,2(1X,D16.6))') 'Min     ',MinDB,MinD
      Write(6,'(A,1P,2(1X,D16.6))') 'Max     ',MaxDB,MaxD
      Write(6,'(A,7X,I10,7X,I10)')  'Negative',nNegDB,nNegD
      Call xFlush(6)

      End
