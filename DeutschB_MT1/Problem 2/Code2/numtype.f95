
module NumType

	save
	integer, parameter :: dp = kind(1.d0)
	real(dp), parameter :: pi = 4*atan(1._dp) ! dp is double precision; defining pi
	complex(dp), parameter :: iic = (0._dp,1._dp)
	
end module NumType
