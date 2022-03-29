!*****************************************************************************
module MD_precision_module
!*****************************************************************************
!
!  Project: MolecularDynamics   
!
!  Created on: Wed  7 Mar 18:30:26 2018  by Eduardo R. Hernandez
!
!  Last Modified: Fri 11 Jan 12:38:43 2019
!
!---
! Defines some precision-related parameters
!---
!
!*****************************************************************************
!  modules used

!*****************************************************************************
!  No implicits please!
   
   implicit none

!*****************************************************************************
!  variables

   integer, parameter :: i8b = selected_int_kind(18)

   integer, parameter :: sp = selected_real_kind(6,30)
   integer, parameter :: dp = selected_real_kind(14,100)

end module MD_precision_module
