      module opacity_switches

!.... 2007 SEP - CHANGED DEFAULT FOR ifop(17) TO .true.

      implicit none

!.... OPACITY

      logical :: if_op(20) = [
!....         H1         H2+        H-         Hray       He1     
     &       .true.,    .true.,    .true.,    .true.,    .true.,  
!....         He2        He-        Heray      Cool       Luke    
     &       .true.,    .true.,    .true.,    .true.,    .true.,
!....         Hot        Elec       H2ray      Hline      Lines   
     &       .true.,    .true.,    .true.,    .false.,   .true.,
!....         Lscat      Xline      Xlsct      Xcont      Xscat   
     &       .false.,   .true.,    .false.,   .false.,   .false. ]

      end module opacity_switches

!***********************************************************************
