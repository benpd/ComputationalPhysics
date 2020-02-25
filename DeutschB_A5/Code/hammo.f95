program hammo
	use numtype
	implicit none
	
	integer :: i, j, n, info
!gotta set that N manuelly, superneat! 
	integer, parameter :: ndim = 1002, lwork=16*ndim
	real(dp), dimension(ndim,ndim) :: x, x2, p, p2, h, h2
	real(dp), dimension(ndim) :: eigen
	real(dp) :: work(lwork), w2
	
	do  i = 1, ndim
		do  j = 1, ndim 

			if ( j-i == 1 ) then                                                                                 
                x(i,j) = sqrt(float(i))
                p(i,j) = -sqrt(float(i))  ! C*C^ is real, no iic
            else if ( i-j == 1 ) then 
				x(i,j) = sqrt(float(j))
                p(i,j) = sqrt(float(j))   ! C*C^ is real, no iic                                                                        
            else
                x(i,j) = 0._dp
                p(i,j) = 0._dp
			end if 
		
		end do 
	end do
x2(1:ndim, 1:ndim)=matmul(x(1:ndim, 1:ndim), x(1:ndim, 1:ndim))
p2(1:ndim, 1:ndim)=matmul(p(1:ndim, 1:ndim), p(1:ndim, 1:ndim))

!omegaguys .5, 1, 2, 3 
w2 = 1
h(1:ndim-1,1:ndim-1) = (1/4._dp)*(-(p2(1:ndim-1,1:ndim-1)) + w2*(x2(1:ndim-1,1:ndim-1))) 
h2(1:ndim-1,1:ndim-1) = h(1:ndim-1,1:ndim-1)

!cooleigenguys
	call dsyev ('n', 'u', ndim-1, h2, ndim, eigen, work, lwork, info)
	print '(f10.2)', eigen(1: ndim-1)
	print *,

!graphthoseguys
		do j = 0, ndim-2                                                                                            
            write(10,*) j, eigen(j+1)
			write(20,*) (j + 0.5_dp)*sqrt(w2)                                                                                                                                
        end do
end program hammo