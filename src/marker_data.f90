MODULE marker_data
  integer, parameter :: max_markers = 10000000
  SAVE
  type marker
     sequence
     real(8) :: a1,a2         ! baricentric coordinates
     real(8) :: x,y           ! Euler coordinates
     real(8) :: age           ! creation time
     integer(4) :: dead
     integer(4) :: ntriag     ! number of FE-triangle
     integer(4) :: phase
     integer(4) :: ID         ! unique ID-number
  end type marker

  type(marker) :: mark(max_markers)

END MODULE marker_data
