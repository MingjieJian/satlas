      module opacity_switches

      implicit none

      logical :: if_op(20) = [
!....         H1         H2+        H-         Hray       He1     
     &       .true.,    .true.,    .true.,    .true.,    .true.,  
!....         He2        He-        Heray      Cool       Luke    
     &       .true.,    .true.,    .true.,    .true.,    .true.,
!....         Hot        Elec       H2ray      Hline      Lines   
     &       .true.,    .true.,    .true.,    .false.,   .true.,
!....         Lscat      Xline      Xlsct      Xcont      Xscat   
     &       .false.,   .false.,   .false.,   .false.,   .false. ]

      end module opacity_switches

!***********************************************************************
