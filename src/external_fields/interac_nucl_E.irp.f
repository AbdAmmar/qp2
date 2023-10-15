
! ---

BEGIN_PROVIDER [double precision, Int_nuclE]

  implicit none
  integer :: l

  BEGIN_DOC
  !
  ! Nuclear interaction with the external Field E
  !
  END_DOC

  PROVIDE nucl_coord nucl_charge nucl_num 
  PROVIDE elec_field

  Int_nuclE = 0.d0
  do l = 1, nucl_num
    Int_nuclE = Int_nuclE + nucl_charge(l) * nucl_coord(l,3)
  enddo
  Int_nuclE = Int_nuclE * elec_field

  call write_double(6, Int_nuclE, 'Nuclear interaction with E')

END_PROVIDER

! ---

