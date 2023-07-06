module ASenD_r_overloadFunctions
    use ASenD_constantVals

    contains
	
	subroutine r_greater(gt,inp1,inp2)
	    implicit none
		
		integer, intent(out) :: gt
		real*8, intent(in) :: inp1, inp2
		
		real*8 :: i2shft
		
		i2shft = inp2 + 1e-1*abs(inp1-inp2)*r_1
		
		if(abs(inp1-i2shft) .lt. abs(inp1-inp2)) then
		    gt = 1
		else
		    gt = 0
		endif
		
	end subroutine r_greater

    subroutine r_sqrt(sqt,inpt)
	    implicit none
		
		real*8, intent(out) :: sqt
		real*8, intent(in) :: inpt
		
		real*8 :: del, fVal, dF
		integer :: ct
		
		sqt = r_1
		del = r_1
		ct = 0
		do while(abs(del) .gt. 1e-13 .and. ct .lt. 50)
		    fVal = sqt*sqt - inpt
			dF = r_2*sqt
			del = -fVal/dF
			sqt = sqt + del
			ct = ct + 1
		enddo
		
	end subroutine r_sqrt

    subroutine r_sin0_90(sT,theta)
	    implicit none
		
		real*8, intent(out) :: sT
		real*8, intent(in) :: theta
		
		real*8 :: coeffs(7), thetSq
		integer :: i1
		
		
		if(abs(theta) .lt. 0.7853981633974483) then
		    thetSq = theta*theta
		    coeffs(1) = r_1
			coeffs(2) = -r_1*1.666666666666666e-1
			coeffs(3) = r_1*8.333333333333333e-3
			coeffs(4) = -r_1*1.984126984126984e-4
			coeffs(5) = r_1*2.755731922398589e-6
			coeffs(6) = -r_1*2.505210838544171e-8
			coeffs(7) = r_1*1.6059043836821614e-10
			sT = coeffs(7)
			do i1 = 6, 1, -1
			    sT = sT*thetSq + coeffs(i1)
			enddo
			sT = sT*theta
		else
		    thetSq = (theta-0.5d0*r_pi)*(theta-0.5d0*r_pi)
		    coeffs(1) = r_1
			coeffs(2) = -r_p5
			coeffs(3) = r_1*4.1666666666666666e-2
			coeffs(4) = -r_1*1.388888888888889e-3
			coeffs(5) = r_1*2.4801587301587301e-5
			coeffs(6) = -r_1*2.755731922398589e-7
			coeffs(7) = r_1*2.0876756987868098e-9
			sT = coeffs(7)
			do i1 = 6, 1, -1
			    sT = sT*thetSq + coeffs(i1)
			enddo
		endif
		
	end subroutine r_sin0_90
	
	subroutine r_sin(sT,theta)
	    implicit none
		
		real*8, intent(out) :: sT
		real*8, intent(in) :: theta
		
		integer :: gt, i1
		real*8 :: mapTheta1, mapTheta2

        mapTheta1 = theta	
	    call r_greater(gt,theta,2d0*r_pi)
	    if(gt .eq. 1) then
		    i1 = abs(r_p5*theta/r_pi)
			mapTheta1 = theta-i1*2d0*r_pi
		endif
		
		call r_greater(gt,r_0,theta)
		if(gt .eq. 1) then
		    i1 = abs(r_p5*theta/r_pi)
			mapTheta1 = theta+(i1+1)*2d0*r_pi
		endif
		
		call r_greater(gt,mapTheta1,1.5d0*r_pi)
		if(gt .eq. 1) then
			mapTheta2 = 2d0*r_pi - mapTheta1
			call r_sin0_90(sT,mapTheta2)
			sT = -sT
			return
		endif
		
		call r_greater(gt,mapTheta1,r_pi)
		if(gt .eq. 1) then
			mapTheta2 = mapTheta1 - r_pi
			call r_sin0_90(sT,mapTheta2)
			sT = -sT
			return
		endif
		
		call r_greater(gt,mapTheta1,0.5d0*r_pi)
		if(gt .eq. 1) then
			mapTheta2 = r_pi - mapTheta1
			call r_sin0_90(sT,mapTheta2)
			return
		endif
		
		call r_sin0_90(sT,mapTheta1)
		
		
	end subroutine r_sin
	
	subroutine r_cos(cT,theta)
	    implicit none
		
		real*8, intent(out) :: cT
		real*8, intent(in) :: theta
		
		real*8 :: thetaAdj
		
		thetaAdj = theta + 0.5d0*r_pi
		
		call r_sin(cT,thetaAdj)
		
	end subroutine r_cos
	
	subroutine r_asin(aS,sn)
	    implicit none
		
		real*8, intent(in) :: sn
		real*8, intent(out) :: aS
        
        real*8 :: dt, fVal, dF
		integer :: i1
		
		aS = r_0
		dt = r_1
		i1 = 0
		do while(abs(dt) .gt. 1e-13 .and. i1 .lt. 20)
            call r_sin(fVal,aS)
			fVal = fVal - sn
			call r_cos(dF,aS)
			dt = -fVal/dF
			aS = aS + dt
			i1 = i1 + 1
		enddo
		
	end subroutine r_asin
	
	subroutine r_acos(aC,cs)
	    implicit none
		
		real*8, intent(in) :: cs
		real*8, intent(out) :: aC
        
        real*8 :: dt, fVal, dF
		integer :: i1
		
		aC = 0.5d0*r_pi
		dt = r_1
		i1 = 0
		do while(abs(dt) .gt. 1e-13 .and. i1 .lt. 20)
            call r_cos(fVal,aC)
			fVal = fVal - cs
			call r_sin(dF,aC)
			dF = -dF
			dt = -fVal/dF
			aC = aC + dt
			i1 = i1 + 1
		enddo
		
	end subroutine r_acos

end module ASenD_r_overloadFunctions