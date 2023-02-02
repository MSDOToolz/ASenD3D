module ASenD_r_elementEqns
    use ASenD_globalData
	use ASenD_constantVals
	use ASenD_r_overloadFunctions
	use ASenD_r_designPropertyFunctions

    contains
	
	!Evaluates the array of basis function values N, at a certain point for a certain type of 
    !element.  Here the element natural coordinats are denoted s1, s2, and s3.
    subroutine r_EvalN(nVec, sPt, eType)
        implicit none
        
        real*8, intent(in) :: sPt(3)
        real*8, intent(out) :: nVec(11)
        integer, intent(in) :: eType
        
        real*8 :: s1, s2, s3
        integer :: i1, i2, i3, i4
        
        s1 = sPt(1)
        s2 = sPt(2)
        s3 = sPt(3)

	    nVec(:) = r_0
        
	    if(eType .eq. 8 .or. eType .eq. 81) then
            nVec(1) = r_p125*(r_1-s1)*(r_1-s2)*(r_1-s3)
            nVec(2) = r_p125*(r_1+s1)*(r_1-s2)*(r_1-s3)
            nVec(3) = r_p125*(r_1+s1)*(r_1+s2)*(r_1-s3)
            nVec(4) = r_p125*(r_1-s1)*(r_1+s2)*(r_1-s3)
            nVec(5) = r_p125*(r_1-s1)*(r_1-s2)*(r_1+s3)
            nVec(6) = r_p125*(r_1+s1)*(r_1-s2)*(r_1+s3)
            nVec(7) = r_p125*(r_1+s1)*(r_1+s2)*(r_1+s3)
            nVec(8) = r_p125*(r_1-s1)*(r_1+s2)*(r_1+s3)
	    endif

	    if(eType .eq. 81) then
            nVec(9) = r_1-s1*s1
            nVec(10) = r_1-s2*s2
            nVec(11) = r_1-s3*s3
	    endif

	    if(eType .eq. 41) then
	        nVec(1) = r_p25*(r_1-s1)*(r_1-s2)
            nVec(2) = r_p25*(r_1+s1)*(r_1-s2)
            nVec(3) = r_p25*(r_1+s1)*(r_1+s2)
            nVec(4) = r_p25*(r_1-s1)*(r_1+s2)
			nVec(5) = r_1-s1*s1
            nVec(6) = r_1-s2*s2
	        nVec(7) = (r_1-s1*s1)*(r_1-s2)
	        nVec(8) = (r_1-s2*s2)*(r_1+s1)
	        nVec(9) = (r_1-s1*s1)*(r_1+s2)
	        nVec(10) = (r_1-s2*s2)*(r_1-s1)
	    endif
        
        if(eType .eq. 4) then
            nVec(1) = r_1-s1-s2-s3
            nVec(2) = s1
            nVec(3) = s2
            nVec(4) = s3
        endif
        
        if(eType .eq. 6) then
            nVec(1) = r_p5*(r_1-s1-s2)*(r_1-s3)
	        nVec(2) = r_p5*s1*(r_1-s3)
	        nVec(3) = r_p5*s2*(r_1-s3)
	        nVec(4) = r_p5*(r_1-s1-s2)*(r_1+s3)
	        nVec(5) = r_p5*s1*(r_1+s3)
	        nVec(6) = r_p5*s2*(r_1+s3)
        endif

		if(eType .eq. 3) then
			nVec(1) = r_1-s1-s2
			nVec(2) = s1
			nVec(3) = s2
			nVec(4) = s1*(r_1-s1-s2)
			nVec(5) = s1*s2
			nVec(6) = s2*(r_1-s1-s2)
		endif

        if(eType .eq. 2) then
            nVec(1) = r_p5*(r_1 - s1)
            nVec(2) = r_p5*(r_1 + s1)
            nVec(3) = r_1 - s1*s1
        endif
        
        return
    end subroutine r_EvalN
    
    !Evaluate the matrix of partial derivatives of the basis functions with respect to natural coordinates
    !at a certain point, for a certain element type.  Here the element natural coordinates are denoted 
    !s1, s2, and s3.
    subroutine r_EvalNs(Ns, sPt, eType)
        implicit none
        
        real*8, intent(in) :: sPt(3)
        real*8, intent(out) :: Ns(11,3)
        integer, intent(in) :: eType
        
        real*8 :: s1, s2, s3
        integer :: i1, i2, i3, i4, i5

	    s1 = sPt(1)
	    s2 = sPt(2)
	    s3 = sPt(3)
	
	    Ns(:,:) = r_0

	    if(eType .eq. 8 .or. eType .eq. 81) then
            Ns(1,1) = -r_p125*(r_1-s2)*(r_1-s3)
            Ns(2,1) = r_p125*(r_1-s2)*(r_1-s3)
            Ns(3,1) = r_p125*(r_1+s2)*(r_1-s3)
            Ns(4,1) = -r_p125*(r_1+s2)*(r_1-s3)
            Ns(5,1) = -r_p125*(r_1-s2)*(r_1+s3)
            Ns(6,1) = r_p125*(r_1-s2)*(r_1+s3)
            Ns(7,1) = r_p125*(r_1+s2)*(r_1+s3)
            Ns(8,1) = -r_p125*(r_1+s2)*(r_1+s3)
	        Ns(1,2) = -r_p125*(r_1-s1)*(r_1-s3)
            Ns(2,2) = -r_p125*(r_1+s1)*(r_1-s3)
            Ns(3,2) = r_p125*(r_1+s1)*(r_1-s3)
            Ns(4,2) = r_p125*(r_1-s1)*(r_1-s3)
            Ns(5,2) = -r_p125*(r_1-s1)*(r_1+s3)
            Ns(6,2) = -r_p125*(r_1+s1)*(r_1+s3)
            Ns(7,2) = r_p125*(r_1+s1)*(r_1+s3)
            Ns(8,2) = r_p125*(r_1-s1)*(r_1+s3)
	        Ns(1,3) = -r_p125*(r_1-s1)*(r_1-s2)
            Ns(2,3) = -r_p125*(r_1+s1)*(r_1-s2)
            Ns(3,3) = -r_p125*(r_1+s1)*(r_1+s2)
            Ns(4,3) = -r_p125*(r_1-s1)*(r_1+s2)
            Ns(5,3) = r_p125*(r_1-s1)*(r_1-s2)
            Ns(6,3) = r_p125*(r_1+s1)*(r_1-s2)
            Ns(7,3) = r_p125*(r_1+s1)*(r_1+s2)
            Ns(8,3) = r_p125*(r_1-s1)*(r_1+s2)
	    endif

	    if(eType .eq. 81) then
            Ns(9,1) = -s1 - s1
            Ns(10,2) = -s2 - s2
            Ns(11,3) = -s3 - s3
	    endif

	    if(eType .eq. 41) then
			Ns(1,1) = -r_p25*(r_1-s2)
			Ns(1,2) = -r_p25*(r_1-s1)
			Ns(2,1) = r_p25*(r_1-s2)
			Ns(2,2) = -r_p25*(r_1+s1)
			Ns(3,1) = r_p25*(r_1+s2)
			Ns(3,2) = r_p25*(r_1+s1)
			Ns(4,1) = -r_p25*(r_1+s2)
			Ns(4,2) = r_p25*(r_1-s1)
			Ns(5,1) = -s1 - s1
			Ns(6,2) = -s2 - s2
			Ns(7,1) = (-s1 - s1)*(r_1-s2)
			Ns(7,2) = -(r_1-s1*s1)
			Ns(8,1) = (r_1-s2*s2)
			Ns(8,2) = (-s2 - s2)*(r_1+s1)
			Ns(9,1) = (-s1 - s1)*(r_1+s2)
			Ns(9,2) = (r_1-s1*s1)
			Ns(10,1) = -(r_1-s2*s2)
			Ns(10,2) = (-s2 - s2)*(r_1-s1)
	    endif

	    if(eType .eq. 4) then
            Ns(1,1) = -r_1
            Ns(1,2) = -r_1
			Ns(1,3) = -r_1
            Ns(2,1) = r_1
            Ns(3,2) = r_1
            Ns(4,3) = r_1
        endif

	    if(eType .eq. 6) then
            Ns(1,1) = -r_p5*(r_1-s3)
			Ns(2,1) = r_p5*(r_1-s3)
			Ns(4,1) = -r_p5*(r_1+s3)
			Ns(5,1) = r_p5*(r_1+s3)
			Ns(1,2) = -r_p5*(r_1-s3)
			Ns(3,2) = r_p5*(r_1-s3)
			Ns(4,2) = -r_p5*(r_1+s3)
			Ns(6,2) = r_p5*(r_1+s3)
			Ns(1,3) = -r_p5*(r_1-s1-s2)
			Ns(2,3) = -r_p5*s1
			Ns(3,3) = -r_p5*s2
			Ns(4,3) = r_p5*(r_1-s1-s2)
			Ns(5,3) = r_p5*s1
			Ns(6,3) = r_p5*s2
        endif

		if(eType .eq. 3) then	
			Ns(1,1) = -r_1
			Ns(1,2) = -r_1
			Ns(2,1) = r_1
			Ns(3,2) = r_1
			Ns(4,1) = r_1 - s1 - s2 - s1
			Ns(4,2) = -s1
			Ns(5,1) = s2
			Ns(5,2) = s1
			Ns(6,1) = -s2
			Ns(6,2) = r_1 - s1 - s2 - s2
		endif

        if(eType .eq. 2) then
            Ns(1,1) = -r_p5
            Ns(2,1) = r_p5
            Ns(3,1) = -s1 - s1
        endif
        
        return
    end subroutine r_EvalNs
	
	subroutine r_getDetInv(D,I,A)
        implicit none
    
        real*8, intent(out) :: I(3,3)
        real*8, intent(in) :: A(3,3)
        real*8, intent(out) :: D
        
        real*8 :: Adj11, Adj12, Adj13, Dinv
        
        Adj11 = A(2,2)*A(3,3) - A(2,3)*A(3,2)
        Adj12 = A(2,3)*A(3,1) - A(2,1)*A(3,3)
        Adj13 = A(2,1)*A(3,2) - A(2,2)*A(3,1)
        
        D = A(1,1)*Adj11 + A(1,2)*Adj12 + A(1,3)*Adj13
        Dinv = r_1/D
        
        I(1,1) = Dinv*Adj11
        I(2,1) = Dinv*Adj12
        I(3,1) = Dinv*Adj13
        I(1,2) = Dinv*(A(1,3)*A(3,2) - A(1,2)*A(3,3))
        I(2,2) = Dinv*(A(1,1)*A(3,3) - A(1,3)*A(3,1))
        I(3,2) = Dinv*(A(1,2)*A(3,1) - A(1,1)*A(3,2))
        I(1,3) = Dinv*(A(1,2)*A(2,3) - A(1,3)*A(2,2))
        I(2,3) = Dinv*(A(1,3)*A(2,1) - A(1,1)*A(2,3))
        I(3,3) = Dinv*(A(1,1)*A(2,2) - A(1,2)*A(2,1))
        return
    end subroutine r_getDetInv
	
	subroutine correctAlpha(alpha,nodes,eType)
        implicit none

        integer, intent(in) :: eType
        real*8, intent(in) :: nodes(3,10)
        real*8, intent(out) :: alpha(3,3)

        real*8 :: dp, mag, theta, v1(3), v2(3), v3(3), alpha0(3,3)
		integer :: check

        select case(eType)
            case(41,3)
                select case(eType)
                   case(41)
                       v1 = nodes(:,3) - nodes(:,1)
                       v2 = nodes(:,4) - nodes(:,2)
                   case(3)
                       v1 = nodes(:,2) - nodes(:,1)
                       v2 = nodes(:,3) - nodes(:,1)
                end select
                v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
                v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
                v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
				v1 = alpha(3,:)
				dp = v1(1)*v3(1) + v1(2)*v3(2) + v1(3)*v3(3)
				call r_greater(check,r_0,dp)
				if(check .eq. 1) then
				    v3 = -v3
				endif
				dp = v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3)
				call r_sqrt(mag,dp)
                v2 = (r_1/mag)*v3   !! Unit vector normal to shell
                v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
                v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
                v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
				dp = v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3)
				call r_sqrt(mag,dp)
                if(abs(mag) .gt. 1e-5) then
				    call r_asin(theta,mag)
                    v3 = (theta/mag)*v3
                    alpha0(:,:) = alpha
					call r_rotateAlpha(alpha,alpha0,v3,0,0)
                endif
			case(2)
			    v1 = alpha(1,:)
			    v2 = nodes(:,2) - nodes(:,1)
				dp = v2(1)+v2(1) + v2(2)*v2(2) + v2(3)*v2(3)
				call r_sqrt(mag,dp)
				v2 = (r_1/mag)*v2
				v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
                v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
                v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
				dp = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
				call r_greater(check,r_0,dp)
				if(check .eq. 1) then
				    v3 = -v3
				endif
				dp = v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3)
				call r_sqrt(mag,dp)
				call r_asin(theta,mag)
				v3 = (theta/mag)*v3
				alpha0(:,:) = alpha
			    call r_rotateAlpha(alpha,alpha0,v3,0,0)
        end select

        return
    end subroutine correctAlpha
	
	subroutine r_getElementFaces(faces,numFaces,elNum)	
	    implicit none
		
		integer, intent(in) :: elNum
		integer, intent(out) :: faces(6,6), numFaces
		
		integer :: eType
		
		eType = elementType(elNum)
		faces(:,:) = 0
		if(eType .eq. 4) then
		    faces(1:3,1) = (/1,3,2/)
			faces(1:3,2) = (/1,2,4/)
			faces(1:3,3) = (/2,3,4/)
			faces(1:3,4) = (/1,4,3/)
			numFaces = 4
		elseif(eType .eq. 6) then
		    faces(1:3,1) = (/1,3,2/)
            faces(1:3,2) = (/4,5,6/)
            faces(1:4,3) = (/1,2,5,4/)
            faces(1:4,4) = (/2,3,6,5/)
            faces(1:4,5) = (/1,4,6,3/)
			numFaces = 5
		elseif(eType .eq. 8 .or. eType .eq. 81) then
		    faces(1:4,1) = (/4,3,2,1/)
			faces(1:4,2) = (/5,6,7,8/)
			faces(1:4,3) = (/1,2,6,5/)
			faces(1:4,4) = (/2,3,7,6/)
			faces(1:4,5) = (/3,4,8,7/)
			faces(1:4,6) = (/4,1,5,8/)
			numFaces = 6
	    elseif(eType .eq. 3) then
		    faces(1:3,1) = (/1,2,3/)
			faces(1:3,2) = (/1,3,2/)
			numFaces = 2
		elseif(eType .eq. 41) then
		    faces(1:4,1) = (/1,2,3,4/)
			faces(1:4,2) = (/1,4,3,2/)
			numFaces = 2
		elseif(eType .eq. 2) then
		    faces(1:2,1) = (/1,2/)
			numFaces = 1
		endif
		
	end subroutine r_getElementFaces
	
	subroutine r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,elNum)
	    implicit none
		
		integer, intent(in) :: elNum
		integer, intent(out) :: numNds, dofPerNd, numIntDof, numIntPts
		integer, intent(out) :: dofTable(2,33)
		real*8, intent(out) :: intPts(3,8), ipWt(8)
		
		integer :: i1, i2, i3, eType
		
		eType = elementType(elNum)
		
		dofTable(:,:) = 0
		intPts(:,:) = r_0
		ipWt(:) = r_0
		if(eType .eq. 4) then
		    numNds = 4
			dofPerNd = 3
			numIntDof = 0
			numIntPts = 1
			intPts(:,1) = r_p25*(/r_1,r_1,r_1/)
			ipWt(1) = r_1o6
		elseif(eType .eq. 6) then
		    numNds = 6
			dofPerNd = 3
			numIntDof = 0
			numIntPts = 2
			intPts(1,1:2) = r_1o3*(/r_1,r_1/)
			intPts(2,1:2) = r_1o3*(/r_1,r_1/)
			intPts(3,1:2) = r_1rt3*(/-r_1,r_1/)
			ipWt(1:2) = r_p5
		elseif(eType .eq. 8) then
		    numNds = 8
			dofPerNd = 3
			numIntDof = 0
			numIntPts = 8
			intPts(:,1) = r_1rt3*(/-r_1,-r_1,-r_1/)
			intPts(:,2) = r_1rt3*(/r_1,-r_1,-r_1/)
			intPts(:,3) = r_1rt3*(/-r_1,r_1,-r_1/)
			intPts(:,4) = r_1rt3*(/r_1,r_1,-r_1/)
			intPts(:,5) = r_1rt3*(/-r_1,-r_1,r_1/)
			intPts(:,6) = r_1rt3*(/r_1,-r_1,r_1/)
			intPts(:,7) = r_1rt3*(/-r_1,r_1,r_1/)
			intPts(:,8) = r_1rt3*(/r_1,r_1,r_1/)
			ipWt(1:8) = r_1
		elseif(eType .eq. 81) then
		    numNds = 8
			dofPerNd = 3
			numIntDof = 9
			numIntPts = 8
			intPts(:,1) = r_1rt3*(/-r_1,-r_1,-r_1/)
			intPts(:,2) = r_1rt3*(/r_1,-r_1,-r_1/)
			intPts(:,3) = r_1rt3*(/-r_1,r_1,-r_1/)
			intPts(:,4) = r_1rt3*(/r_1,r_1,-r_1/)
			intPts(:,5) = r_1rt3*(/-r_1,-r_1,r_1/)
			intPts(:,6) = r_1rt3*(/r_1,-r_1,r_1/)
			intPts(:,7) = r_1rt3*(/-r_1,r_1,r_1/)
			intPts(:,8) = r_1rt3*(/r_1,r_1,r_1/)
			i3 = 24
			do i1 = 1, 3
			    do i2 = 9, 11
				    i3 = i3 + 1
					dofTable(1,i3) = i1
					dofTable(2,i3) = i2
				enddo
			enddo
			ipWt(1:8) = r_1
		elseif(eType .eq. 41) then
		    numNds = 4
			dofPerNd = 6
			numIntDof = 8
			numIntPts = 4
			intPts(:,1) = r_1rt3*(/-r_1,-r_1,r_0/)
			intPts(:,2) = r_1rt3*(/r_1,-r_1,r_0/)
			intPts(:,3) = r_1rt3*(/-r_1,r_1,r_0/)
			intPts(:,4) = r_1rt3*(/r_1,r_1,r_0/)
			! dofTable(:,25) = (/3,7/)
			! dofTable(:,26) = (/3,8/)
			! dofTable(:,27) = (/3,9/)
			! dofTable(:,28) = (/3,10/)
			dofTable(:,25) = (/1,5/)
			dofTable(:,26) = (/1,6/)
			dofTable(:,27) = (/2,5/)
			dofTable(:,28) = (/2,6/)
			dofTable(:,29) = (/3,7/)
			dofTable(:,30) = (/3,8/)
			dofTable(:,31) = (/3,9/)
			dofTable(:,32) = (/3,10/)
			ipWt(1:4) = r_1
		elseif(eType .eq. 3) then
		    numNds = 3
			dofPerNd = 6
			numIntDof = 3
			numIntPts = 3
			intPts(:,1) = r_1o6*(/r_1,r_1,r_0/)
			intPts(:,2) = r_1o6*(/4d0,1d0,0d0/)
			intPts(:,3) = r_1o6*(/1d0,4d0,0d0/)
			dofTable(:,19) = (/3,4/)
			dofTable(:,20) = (/3,5/)
			dofTable(:,21) = (/3,6/)
			ipWt(1:3) = r_1o6
		elseif(eType .eq. 2) then
		    numNds = 2
			dofPerNd = 6
			numIntDof = 2
			numIntPts = 2
			intPts(1,1:2) = r_1rt3*(/-r_1,r_1/)
			dofTable(:,13) = (/2,3/)
			dofTable(:,14) = (/3,3/)
			ipWt(1:2) = r_1
		endif
		
		i3 = 0
		do i1 = 1, dofPerNd
		    do i2 = 1, numNds
			    i3 = i3 + 1
				dofTable(1,i3) = i1
				dofTable(2,i3) = i2
			enddo
		enddo
		
	end subroutine r_getElementProfile
	
	subroutine r_getElementProperties(locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, &
	             ABD, sMass, sTELd, stCond, ssHeat, bStiff, bMass, bTELd, btCond, bsHeat, elNum)
	    implicit none
		
		integer, intent(in) :: elNum
		real*8, intent(out) :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8, intent(out) :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		real*8, intent(out) :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		
		integer :: i1, i2, i3, eType
		
		eType = elementType(elNum)
		
		call r_getNodeCoord(globNds,elNum)
		call r_getOrientation(orient,elNum)
		
		if(eType .eq. 41 .or. eType .eq. 3) then
		    call correctAlpha(orient,globNds,eType)
		endif
		
		locNds(:,:) = r_0
		do i1 = 1, 3
		    do i2 = 1, 10
			    do i3 = 1, 3
				    locNds(i1,i2) = locNds(i1,i2) + orient(i1,i3)*globNds(i3,i2)
				enddo
			enddo
		enddo
		
		if(eType .eq. 2) then
			call r_getBeamStiffness(bStiff,elNum)
			call r_getBeamMass(bMass,elNum)
			call r_getBeamExpLoad(bTELd,elNum)
			call r_getBeamTCond(btCond,elNum)
			call r_getBeamSpecHeat(bsHeat,elNum)
		elseif(eType .eq. 41 .or. eType .eq. 3) then
			call r_getShellStiffness(ABD,elNum)
			call r_getShellMass(sMass,elNum)
			call r_getShellExpLoad(sTELd,elNum)
			call r_getShellThermCond(stCond,elNum)
			call r_getShellSpecHeat(ssHeat,elNum)
        else
			call r_getMaterialStiffness(cMat,elNum)
			call r_getMaterialDensity(den,elNum)
			call r_getMaterialThermExp(tExp,elNum)
			call r_getMaterialThermCond(tCond,elNum)
			call r_getMaterialSpecHeat(sHeat,elNum)
		endif
		
		
	end subroutine r_getElementProperties
	
	subroutine r_getElementSolution(temp, Tdot, pTemp, pTdot, disp, vel, acc, pDisp, pVel, pAcc, &
	                                numNds, dofPerNd, numIntDof, dofTable, elNum)
	    implicit none
		
		integer, intent(in) :: elNum, numNds, dofPerNd, numIntDof, dofTable(2,33)
		real*8, intent(out) :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8, intent(out) :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
		
		integer :: i1, i2, i3, i4, i5, i6
		
		temp(:) = r_0
		Tdot(:) = r_0
		pTemp(:) = r_0
		pTdot(:) = r_0
		disp(:,:) = r_0
		vel(:,:) = r_0
		acc(:,:) = r_0
		pDisp(:,:) = r_0
		pVel(:,:) = r_0
		pAcc(:,:) = r_0
		
		do i1 = 1, numNds
			i3 = currentRank(elementList(i1,elNum))
			pTemp(i1) = r_1*prevTemp(i3)
			pTdot(i1) = r_1*prevTdot(i3)
			temp(i1) = r_1*nodeTemp(i3)
			Tdot(i1) = r_1*nodeTdot(i3)
			do i2 = 1, dofPerNd
			    i4 = nDofIndex(i2,elementList(i1,elNum))
				pDisp(i2,i1) = r_1*prevDisp(i4)
				pAcc(i2,i1) = r_1*prevAcc(i4)
				pVel(i2,i1) = r_1*prevVel(i4)
				disp(i2,i1) = r_1*nodeDisp(i4)
				acc(i2,i1) = r_1*nodeAcc(i4)
				vel(i2,i1) = r_1*nodeVel(i4)
			enddo
		enddo
		
		if(numIntDof .ne. 0) then
		    i1 = intVecRange(elNum-1)
		    i2 = numNds*dofPerNd
		    do i3 = 1, numIntDof
			    i4 = i2 + i3
				i5 = dofTable(1,i4)
				i6 = dofTable(2,i4)
				disp(i5,i6) = r_1*internalDisp(i1+i3)
				pDisp(i5,i6) = r_1*prevIntDisp(i1+i3)
			enddo
		endif
		
	end subroutine r_getElementSolution
	
	subroutine r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
	    locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
		ABD, sMass, sTELd, stCond, ssHeat, &
		bStiff, bMass, bTELd, btCond, bsHeat, &
	    temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,elNum)
	    implicit none
		
		integer, intent(in) :: elNum
		integer, intent(out) :: numNds, dofPerNd, numIntDof, numIntPts
		integer, intent(out) :: dofTable(2,33)
		real*8, intent(out) :: intPts(3,8), ipWt(8)
		real*8, intent(out) :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8, intent(out) :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		real*8, intent(out) :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		real*8, intent(out) :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8, intent(out) :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
		
		integer :: eType
		integer :: i1, i2, i3, i4, i5, i6
		
        call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,elNum)
		call r_getElementProperties(locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, &
	             ABD, sMass, sTELd, stCond, ssHeat, bStiff, bMass, bTELd, btCond, bsHeat, elNum)
		call r_getElementSolution(temp, Tdot, pTemp, pTdot, disp, vel, acc, pDisp, pVel, pAcc, &
	                                numNds, dofPerNd, numIntDof, dofTable, elNum)
		
	end subroutine r_getElementData
	
	subroutine r_getIntPtData(Nvec,Nx,detJ,locNds,sPt,eType,numNds)
	    implicit none
		
		integer, intent(in) :: eType, numNds
		real*8, intent(out) :: Nvec(11), Nx(11,3), detJ
		real*8, intent(in) :: locNds(3,10), sPt(3)
		
		real*8 :: Ns(11,3), dxds(3,3), dsdx(3,3), v3(3), mag, dp
		integer :: i1, i2, i3, gt
		
		call r_EvalN(Nvec, sPt, eType)
		call r_EvalNs(Ns, sPt, eType)
		dxds(:,:) = r_0
		do i1 = 1, 3
		    do i2 = 1, 3
			    do i3 = 1, numNds
				    dxds(i1,i2) = dxds(i1,i2) + locNds(i1,i3)*Ns(i3,i2)
				enddo
			enddo
		enddo
		
		if(eType .eq. 41 .or. eType .eq. 3) then
		    v3(1) = dxds(2,1)*dxds(3,2) - dxds(3,1)*dxds(2,2)
			v3(2) = dxds(3,1)*dxds(1,2) - dxds(1,1)*dxds(3,2)
			v3(3) = dxds(1,1)*dxds(2,2) - dxds(2,1)*dxds(1,2)
			dp = v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3)
			call r_sqrt(mag,dp)
			dxds(:,3) = (r_1/mag)*v3(:)
		elseif(eType .eq. 2) then
		    dxds(2,2) = r_1
			call r_greater(gt,dxds(1,1),r_0)
			if(gt .eq. 1) then
			    dxds(3,3) = r_1
			else
			    dxds(3,3) = -r_1
			endif
		endif
		
		call r_getDetInv(detJ,dsdx,dxds)
		
		Nx(:,:) = r_0
		
		do i1 = 1, 11
		    do i2 = 1, 3
			    do i3 = 1, 3
				    Nx(i1,i2) = Nx(i1,i2) + Ns(i1,i3)*dsdx(i3,i2)
				enddo
			enddo
		enddo
	
	end subroutine r_getIntPtData
	
	subroutine r_getNLStrain(strain,ux,Nx,orient,dv,dv2,dofTable)
        implicit none

        integer, intent(in) :: dv, dv2
		integer, intent(in) :: dofTable(2,33)
        real*8, intent(in) :: ux(3,3), Nx(11,3), orient(3,3)
        real*8, intent(out) :: strain(6)

        real*8 :: uxL(3,3)
        integer :: i1, i2, i3, nd, var, nd2, var2

        strain(:) = r_0
		if(dv + dv2 .eq. 0) then
		    uxL(:,:) = r_0
		    do i1 = 1, 3
			    do i2 = 1, 3
				    do i3 = 1, 3
					    uxL(i1,i2) = uxL(i1,i2) + orient(i1,i3)*ux(i3,i2)
					enddo
				enddo
			enddo
			strain(1) = uxL(1,1)
			strain(2) = uxL(2,2)
			strain(3) = uxL(3,3)
			strain(4) = uxL(1,2) + uxL(2,1)
			strain(5) = uxL(1,3) + uxL(3,1)
			strain(6) = uxL(2,3) + uxL(3,2)
			if(nLGeom .ne. 0) then
				strain(1) = strain(1) + r_p5*(ux(1,1)*ux(1,1) + ux(2,1)*ux(2,1) + ux(3,1)*ux(3,1))
				strain(2) = strain(2) + r_p5*(ux(1,2)*ux(1,2) + ux(2,2)*ux(2,2) + ux(3,2)*ux(3,2))
				strain(3) = strain(3) + r_p5*(ux(1,3)*ux(1,3) + ux(2,3)*ux(2,3) + ux(3,3)*ux(3,3))
				strain(4) = strain(4) + ux(1,1)*ux(1,2) + ux(2,1)*ux(2,2) + ux(3,1)*ux(3,2)
				strain(5) = strain(5) + ux(1,1)*ux(1,3) + ux(2,1)*ux(2,3) + ux(3,1)*ux(3,3)
				strain(6) = strain(6) + ux(1,3)*ux(1,2) + ux(2,3)*ux(2,2) + ux(3,3)*ux(3,2)
			endif
		elseif(dv*dv2 .eq. 0) then
		    if(dv .ne. 0) then
			    var = dofTable(1,dv)
			    nd = dofTable(2,dv)
			else
			    var = dofTable(1,dv2)
				nd = dofTable(2,dv2)
			endif
			strain(1) = orient(1,var)*Nx(nd,1)
			strain(2) = orient(2,var)*Nx(nd,2)
			strain(3) = orient(3,var)*Nx(nd,3)
			strain(4) = orient(1,var)*Nx(nd,2) + orient(2,var)*Nx(nd,1)
			strain(5) = orient(1,var)*Nx(nd,3) + orient(3,var)*Nx(nd,1)
			strain(6) = orient(2,var)*Nx(nd,3) + orient(3,var)*Nx(nd,2)
			if(nLGeom .ne. 0) then
				strain(1) = strain(1) + Nx(nd,1)*ux(var,1)
				strain(2) = strain(2) + Nx(nd,2)*ux(var,2)
				strain(3) = strain(3) + Nx(nd,3)*ux(var,3)
				strain(4) = strain(4) + Nx(nd,1)*ux(var,2) + Nx(nd,2)*ux(var,1)
				strain(5) = strain(5) + Nx(nd,1)*ux(var,3) + Nx(nd,3)*ux(var,1)
				strain(6) = strain(6) + Nx(nd,2)*ux(var,3) + Nx(nd,3)*ux(var,2)
			endif
		else
		    if(nLGeom .ne. 0) then
				var = dofTable(1,dv)
				nd = dofTable(2,dv)
				var2 = dofTable(1,dv2)
				nd2 = dofTable(2,dv2)
				if(var .eq. var2) then
					strain(1) = Nx(nd,1)*Nx(nd2,1)
					strain(2) = Nx(nd,2)*Nx(nd2,2)
					strain(3) = Nx(nd,3)*Nx(nd2,3)
					strain(4) = Nx(nd,1)*Nx(nd2,2) + Nx(nd2,1)*Nx(nd,2)
					strain(5) = Nx(nd,1)*Nx(nd2,3) + Nx(nd2,1)*Nx(nd,3)
					strain(6) = Nx(nd,2)*Nx(nd2,3) + Nx(nd2,2)*Nx(nd,3)
				endif
			endif
		endif

        return
    end subroutine r_getNLStrain
	
	subroutine r_getShellDef(def,ux,rot,rx)
	    implicit none
		
		real*8, intent(out) :: def(9)
		real*8, intent(in) :: ux(3,3), rot(3), rx(3,3)
		
		def(1) = ux(1,1)
		def(2) = ux(2,2)
		def(3) = ux(1,2) + ux(2,1)
		def(4) = rx(2,1)
		def(5) = -rx(1,2)
		def(6) = rx(2,2) - rx(1,1)
		def(7) = ux(3,1) + rot(2)
		def(8) = ux(3,2) - rot(1)
		def(9) = r_2*rot(3) - ux(2,1) + ux(1,2)
		
	end subroutine r_getShellDef
	
	subroutine r_getBeamDef(def,ux,rot,rx)
	    implicit none
		
		real*8, intent(out) :: def(6)
		real*8, intent(in) :: ux(3,3), rot(3), rx(3,3)
		
		def(1) = ux(1,1)
		def(2) = ux(2,1) - rot(3)
		def(3) = ux(3,1) + rot(2)
		def(4) = rx(1,1)
		def(5) = rx(2,1)
		def(6) = rx(3,1)
		
	end subroutine r_getBeamDef
	
	subroutine r_getInstOrient(instOri,locOri,disp,numNds,static)
	    implicit none

        integer, intent(in) :: numNds, static
		real*8, intent(in) :: locOri(3,3), disp(6,11)
		real*8, intent(out) :: instOri(3,240)
		
		real*8 :: prod(3,3), tempAl(3,3), tempAl2(3,3), elRot(3)
		
		integer :: i1, i2, i3, i4, i5, stCol
		
		instOri(:,:) = r_0
		
		if(nLGeom .eq. 0) then
		    instOri(:,1:3) = locOri(:,:)
			
			do i1 = 1, numNds
			    stCol = 48*i1
				instOri(:,stCol+1:stCol+3) = locOri(:,:)
			enddo
			
			do i1 = 1, 3
			    elRot(:) = r_0
			    call r_rotateAlpha(tempAl,locOri,elRot,i1,0)
				do i2 = 1, numNds
				    stCol = 48*i2 + 12*i1
					instOri(:,stCol+1:stCol+3) = tempAl(:,:)
					stCol = 48*i2 + 3*i1
					instOri(:,stCol+1:stCol+3) = tempAl(:,:)
					stCol = 48*i2
				    instOri(:,stCol+1:stCol+3) = instOri(:,stCol+1:stCol+3) + tempAl(:,:)*disp(i1+3,i2)
				enddo
			enddo
		else
			elRot(:) = r_0
			do i1 = 1, numNds
				elRot(:) = elRot(:) + disp(4:6,i1)
			enddo
			elRot(:) = (r_1/numNds)*elRot(:)
			
			call r_rotateAlpha(instOri(:,1:3),locOri,elRot,0,0)
			if(static .eq. 0) then
				do i1 = 1, 3
					call r_rotateAlpha(tempAl,locOri,elRot,i1,0)
					stCol = 12*i1
					instOri(:,stCol+1:stCol+3) = tempAl(:,:)
					stCol = 3*i1
					instOri(:,stCol+1:stCol+3) = tempAl(:,:)
					do i2 = i1, 3
						call r_rotateAlpha(tempAl,locOri,elRot,i1,i2)
						stCol = 12*i1 + 3*i2
						instOri(:,stCol+1:stCol+3) = tempAl(:,:)
						stCol = 3*i1 + 12*i2
						instOri(:,stCol+1:stCol+3) = tempAl(:,:)
					enddo
				enddo
			endif
			
			do i3 = 1, numNds
				stCol = 48*i3
				call r_rotateAlpha(instOri(:,stCol+1:stCol+3),locOri,disp(4:6,i3),0,0)
				do i1 = 1, 3
					call r_rotateAlpha(tempAl,locOri,disp(4:6,i3),i1,0)
					stCol = 48*i3 + 12*i1
					instOri(:,stCol+1:stCol+3) = tempAl(:,:)
					stCol = 48*i3 + 3*i1
					instOri(:,stCol+1:stCol+3) = tempAl(:,:)
					do i2 = i1, 3
						call r_rotateAlpha(tempAl,locOri,disp(4:6,i3),i1,i2)
						stCol = 48*i3 + 12*i1 + 3*i2
						instOri(:,stCol+1:stCol+3) = tempAl(:,:)
						stCol = 48*i3 + 3*i1 + 12*i2
						instOri(:,stCol+1:stCol+3) = tempAl(:,:)
					enddo
				enddo
			enddo
		endif
		

		
	end subroutine r_getInstOrient
	
	subroutine r_getInstDof(instU,globU,globNds,instOri,locOri,numNds,dofTable,dv,dv2)
	    implicit none
		
		integer, intent(in) :: numNds, dv, dv2
		integer, intent(in) :: dofTable(2,33)
		real*8, intent(in) :: globU(6,11), globNds(3,10), instOri(3,240), locOri(3,3)
		real*8, intent(out) :: instU(6,11)

        real*8 :: nnInv, nnInv2, delRot(3), ddelRot(3)
        integer :: i1, i2, i3, i4, i5, var, nd, var2, nd2, stCol1, stCol2
		
		if(dv + dv2 .eq. 0) then
		    instU(:,:) = r_0
			do i1 = 1, 3
			    do i2 = 1, numNds
				    do i3 = 1, 3
					    instU(i1,i2) = instU(i1,i2) + instOri(i1,i3)*(globU(i3,i2) + globNds(i3,i2)) - locOri(i1,i3)*globNds(i3,i2)
					enddo
				enddo
			enddo
			do i1 = 1, 33
			    if(dofTable(2,i1) .gt. numNds) then
				    i2 = dofTable(1,i1)
					i3 = dofTable(2,i1)
				    instU(i2,i3) = globU(i2,i3)
				endif
			enddo
			do i1 = 1, numNds
			    instU(4:6,i1) = r_0
				stCol1 = 48*i1
				stCol2 = 0
				do i2 = 1, 3
				    instU(4,i1) = instU(4,i1) + instOri(2,stCol1+i2)*instOri(3,stCol2+i2)
					instU(5,i1) = instU(5,i1) + instOri(3,stCol1+i2)*instOri(1,stCol2+i2)
					instU(6,i1) = instU(6,i1) + instOri(1,stCol1+i2)*instOri(2,stCol2+i2)
				enddo
			enddo
		elseif(dv*dv2 .eq. 0) then
		    if(dv .ne. 0) then
			    nd = dofTable(2,dv)
				var = dofTable(1,dv)
			else
			    nd = dofTable(2,dv2)
				var = dofTable(1,dv2)
			endif
		    if(var .gt. 0 .and. var .lt. 4) then  !! if we're differentiating wrt a global displacement
			    if(nd .le. numNds) then  !! if it's a nodal displacement
					instU(:,:) = r_0
					instU(1:3,nd) = instOri(:,var)
				else
				    instU(:,:) = r_0
					instU(var,nd) = r_1
				endif
			else
				nnInv = numNds
				nnInv = r_1/nnInv
				instU(:,:) = r_0
				stCol1 = 12*(var-3)
				do i1 = 1, 3
					do i2 = 1, numNds
						do i3 = 1, 3
							instU(i1,i2) = instU(i1,i2) + nnInv*instOri(i1,stCol1+i3)*(globU(i3,i2) + globNds(i3,i2))
						enddo
					enddo
				enddo
				instU(4:6,i1) = r_0
				do i1 = 1, numNds
					stCol1 = 48*i1
					stCol2 = 12*(var-3)
					do i2 = 1, 3
						instU(4,i1) = instU(4,i1) + nnInv*instOri(2,stCol1+i2)*instOri(3,stCol2+i2)
						instU(5,i1) = instU(5,i1) + nnInv*instOri(3,stCol1+i2)*instOri(1,stCol2+i2)
						instU(6,i1) = instU(6,i1) + nnInv*instOri(1,stCol1+i2)*instOri(2,stCol2+i2)
					enddo
				enddo
				i1 = nd
				stCol1 = 48*i1 + 12*(var-3)
				stCol2 = 0
				do i2 = 1, 3
					instU(4,i1) = instU(4,i1) + instOri(2,stCol1+i2)*instOri(3,stCol2+i2)
					instU(5,i1) = instU(5,i1) + instOri(3,stCol1+i2)*instOri(1,stCol2+i2)
					instU(6,i1) = instU(6,i1) + instOri(1,stCol1+i2)*instOri(2,stCol2+i2)
				enddo
			endif
		else
		    if(dofTable(1,dv) .lt. dofTable(1,dv2)) then
				nd = dofTable(2,dv)
				var = dofTable(1,dv)
				nd2 = dofTable(2,dv2)
				var2 = dofTable(1,dv2)
			else
			    nd = dofTable(2,dv2)
				var = dofTable(1,dv2)
				nd2 = dofTable(2,dv)
				var2 = dofTable(1,dv)
			endif
			nnInv = numNds
			nnInv = r_1/nnInv
			instU(:,:) = r_0
			if(var .lt. 4 .and. var2 .ge. 4 .and. nd .le. numNds) then
			    stCol1 = 12*(var2-3)
				instU(1:3,nd) = nnInv*instOri(:,stCol1+var)
			elseif(var .ge. 4 .and. var2 .ge. 4) then
			    nnInv2 = nnInv*nnInv
			    stCol1 = 12*(var-3) + 3*(var2-3)
				do i2 = 1, 3
				    do i3 = 1, 3
					    do i4 = 1, numNds
						    instU(i2,i4) = instU(i2,i4) + instOri(i2,stCol1+i3)*(globU(i3,i4) + globNds(i3,i4))
						enddo
					enddo
				enddo
				instU(1:3,1:numNds) = nnInv2*instU(1:3,1:numNds)
				do i1 = 1, numNds
					stCol1 = 48*i1
					stCol2 = 12*(var-3) + 3*(var2-3)
					do i2 = 1, 3
						instU(4,i1) = instU(4,i1) + nnInv2*instOri(2,stCol1+i2)*instOri(3,stCol2+i2)
						instU(5,i1) = instU(5,i1) + nnInv2*instOri(3,stCol1+i2)*instOri(1,stCol2+i2)
						instU(6,i1) = instU(6,i1) + nnInv2*instOri(1,stCol1+i2)*instOri(2,stCol2+i2)
					enddo
				enddo
				i1 = nd
				stCol1 = 48*i1 + 12*(var-3)
				stCol2 = 12*(var2-3)
				do i2 = 1, 3
					instU(4,i1) = instU(4,i1) + nnInv*instOri(2,stCol1+i2)*instOri(3,stCol2+i2)
					instU(5,i1) = instU(5,i1) + nnInv*instOri(3,stCol1+i2)*instOri(1,stCol2+i2)
					instU(6,i1) = instU(6,i1) + nnInv*instOri(1,stCol1+i2)*instOri(2,stCol2+i2)
				enddo
				i1 = nd2
				stCol1 = 48*i1 + 12*(var2-3)
				stCol2 = 12*(var-3)
				do i2 = 1, 3
					instU(4,i1) = instU(4,i1) + nnInv*instOri(2,stCol1+i2)*instOri(3,stCol2+i2)
					instU(5,i1) = instU(5,i1) + nnInv*instOri(3,stCol1+i2)*instOri(1,stCol2+i2)
					instU(6,i1) = instU(6,i1) + nnInv*instOri(1,stCol1+i2)*instOri(2,stCol2+i2)
				enddo
				if(nd .eq. nd2) then
				    i1 = nd
					stCol1 = 48*i1 + 12*(var-3) + 3*(var2-3)
					stCol2 = 0
					do i2 = 1, 3
						instU(4,i1) = instU(4,i1) + instOri(2,stCol1+i2)*instOri(3,stCol2+i2)
						instU(5,i1) = instU(5,i1) + instOri(3,stCol1+i2)*instOri(1,stCol2+i2)
						instU(6,i1) = instU(6,i1) + instOri(1,stCol1+i2)*instOri(2,stCol2+i2)
					enddo
				endif
			endif
		endif
		
	end subroutine r_getInstDof
	
	subroutine r_getElementDisp(disp,ddispdU,elNum,sPt)
	    implicit none
		
		integer, intent(in) :: elNum
		real*8, intent(in) :: spt(3)
		real*8, intent(out) :: disp(6), ddispdU(6,33)
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		
		real*8 :: Nvec(11)
        integer :: i1, i2, i3, i4
		
		call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,elNum)
		
		eType = elementType(elNum)
		call r_EvalN(Nvec, sPt, eType)
		disp(:) = r_0
		ddispdU(:,:) = r_0
		do i1 = 1, numNds*dofPerNd
		    i2 = dofTable(1,i1)
			i3 = dofTable(2,i1)
			i4 = nDofIndex(i2,elementList(i3,elNum))
			disp(i2) = disp(i2) + r_1*nodeDisp(i4)*Nvec(i3)
			ddispdU(i2,i1) = ddispdU(i2,i1) + Nvec(i3)
		enddo
	
	end subroutine r_getElementDisp
	
	subroutine r_getSolidStress(stress,strain,elNum,sPt)
	    implicit none
		
		integer, intent(in) :: elNum
		real*8, intent(in) :: sPt(3)
		real*8, intent(out) :: stress(6), strain(6)
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
        
        real*8 :: ux(3,3), Nx(11,3), Nvec(11), detJ, PtT
		integer :: i1, i2, i3
		
		eType = elementType(elNum)
		call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
			locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
			ABD, sMass, sTELd, stCond, ssHeat, &
			bStiff, bMass, bTELd, btCond, bsHeat, &
			temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,elNum)
		
		call r_getIntPtData(Nvec,Nx,detJ,locNds,sPt,eType,numNds)
		
		ux(:,:) = r_0
		do i1 = 1, 3
		    do i2 = 1, 3
			    do i3 = 1, 11
				    ux(i1,i2) = ux(i1,i2) + disp(i1,i3)*Nx(i3,i2)
				enddo
			enddo
		enddo
		
		call r_getNLStrain(strain,ux,Nx,orient,0,0,dofTable)
		
		PtT = r_0
		do i1 = 1, numNds
		    PtT = PtT + temp(i1)*Nvec(i1)
		enddo
		
		stress(:) = r_0
		do i1 = 1, 6
		    do i2 = 1, 6
			    stress(i1) = stress(i1) + cMat(i1,i2)*(strain(i2) - PtT*tExp(i2))
			enddo
		enddo
		
	end subroutine r_getSolidStress
	
	subroutine r_getdSolidStressdU(dsdU,dsdT,dedU,elNum,sPt)
	    implicit none
		
		integer, intent(in) :: elNum
		real*8, intent(in) :: sPt(3)
		real*8, intent(out) :: dsdU(6,33), dsdT(6,10), dedU(6,33)
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
        
        real*8 :: ux(3,3), Nx(11,3), Nvec(11), detJ, PtT
		real*8 :: dPtTdT(10)
		integer :: i1, i2, i3, totDof
		
		eType = elementType(elNum)
		call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
			locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
			ABD, sMass, sTELd, stCond, ssHeat, &
			bStiff, bMass, bTELd, btCond, bsHeat, &
			temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,elNum)
		
		call r_getIntPtData(Nvec,Nx,detJ,locNds,sPt,eType,numNds)
		
		ux(:,:) = r_0
		do i1 = 1, 3
		    do i2 = 1, 3
			    do i3 = 1, 11
				    ux(i1,i2) = ux(i1,i2) + disp(i1,i3)*Nx(i3,i2)
				enddo
			enddo
		enddo
		
		totDof = numNds*dofPerNd + numIntDof
		do i2 = 1, totDof
		    call r_getNLStrain(dedU(:,i2),ux,Nx,orient,i2,0,dofTable)
		enddo
		
		dPtTdT(:) = r_0
		do i1 = 1, numNds
		    dPtTdT(i1) = Nvec(i1)
		enddo
		
		dsdU(:,:) = r_0
		do i1 = 1, 6
		    do i2 = 1, 6
			    do i3 = 1, totDof
			        dsdU(i1,i3) = dsdU(i1,i3) + cMat(i1,i2)*dedU(i2,i3)
				enddo
			enddo
		enddo
		
		dsdT(:,:) = r_0
		do i1 = 1, 6
		    do i2 = 1, 6
			    do i3 = 1, numNds
			        dsdT(i1,i3) = dsdT(i1,i3) - cMat(i1,i2)*dPtTdT(i3)*tExp(i2)
				enddo
			enddo
		enddo
		
	end subroutine r_getdSolidStressdU
	
	subroutine r_getShellFrcMom(frcMom,shDef,ptT,elNum,sPt)
	    implicit none
		
		integer, intent(in) :: elNum
		real*8, intent(in) :: sPt(3)
		real*8, intent(out) :: frcMom(9), shDef(9), ptT
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
        
        real*8 :: ux(3,3), rx(3,3), Nx(11,3), Nvec(11), detJ
		real*8 :: instOri(3,240), drIdrG(3,16), rot(3), ptRot(3)
		real*8 :: instU(6,11)
		real*8 :: rnNds
		integer :: i1, i2, i3
		
		eType = elementType(elNum)
		if(eType .ne. 41 .and. eType .ne. 3) then
		    return
		endif
		call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
			locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
			ABD, sMass, sTELd, stCond, ssHeat, &
			bStiff, bMass, bTELd, btCond, bsHeat, &
			temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,elNum)
		
		call r_getIntPtData(Nvec,Nx,detJ,locNds,sPt,eType,numNds)
		
		call r_getInstOrient(instOri,orient,disp,numNds,0)
		call r_getInstDof(instU,disp,globNds,instOri,orient,numNds,dofTable,0,0)
		
		ux(:,:) = r_0
		rx(:,:) = r_0
		do i1 = 1, 3
		    do i2 = 1, 3
			    do i3 = 1, 10
				    ux(i1,i2) = ux(i1,i2) + instU(i1,i3)*Nx(i3,i2)
					rx(i1,i2) = rx(i1,i2) + instU(i1+3,i3)*Nx(i3,i2)
				enddo
			enddo
		enddo
		
		ptRot(:) = r_0
		do i1 = 1, 3
		    do i2 = 1, numNds
			    ptRot(i1) = ptRot(i1) + instU(i1+3,i2)*Nvec(i2)
			enddo
		enddo
		
		call r_getShellDef(shDef,ux,ptRot,rx)
		ptT = r_0
		do i2 = 1, numNds
			ptT = ptT + Nvec(i2)*temp(i2)
		enddo
		frcMom(:) = r_0
		do i2 = 1, 9
			do i3 = 1, 9
				frcMom(i2) = frcMom(i2) + ABD(i2,i3)*shDef(i3)
			enddo
		enddo
		frcMom(1:6) = frcMom(1:6) - ptT*sTELd(1:6)
	
	end subroutine r_getShellFrcMom
	
	subroutine r_getdShellFrcMomdU(dfMdU,dfMdT,dsDdU,dPtTdT,elNum,sPt)
	    implicit none
		
		integer, intent(in) :: elNum
		real*8, intent(in) :: sPt(3)
		real*8, intent(out) :: dfMdU(9,33), dfMdT(9,10), dsDdU(9,33), dPtTdT(10)
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
        
        real*8 :: ux(3,3), rx(3,3), Nx(11,3), Nvec(11), detJ
		real*8 :: instOri(3,240), drIdrG(3,16), rot(3), ptRot(3)
		real*8 :: instU(6,11)
		real*8 :: rnNds
		integer :: i1, i2, i3, i4, i5
		
		eType = elementType(elNum)
		if(eType .ne. 41 .and. eType .ne. 3) then
		    return
		endif
		call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
			locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
			ABD, sMass, sTELd, stCond, ssHeat, &
			bStiff, bMass, bTELd, btCond, bsHeat, &
			temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,elNum)
		
		call r_getIntPtData(Nvec,Nx,detJ,locNds,sPt,eType,numNds)
		
		call r_getInstOrient(instOri,orient,disp,numNds,0)
		i1 = numNds*dofPerNd + numIntDof
		do i5 = 1, i1
			call r_getInstDof(instU,disp,globNds,instOri,orient,numNds,dofTable,i5,0)
			
			ux(:,:) = r_0
			rx(:,:) = r_0
			if(dofTable(1,i5) .le. 3) then
			    i3 = dofTable(2,i5)
				do i1 = 1, 3
					do i2 = 1, 3
						ux(i1,i2) = ux(i1,i2) + instU(i1,i3)*Nx(i3,i2)
					enddo
				enddo
			else
			    i4 = dofTable(2,i5)
			    do i1 = 1, 3
					do i2 = 1, 3
						do i3 = 1, 10
							ux(i1,i2) = ux(i1,i2) + instU(i1,i3)*Nx(i3,i2)
						enddo
						rx(i1,i2) = rx(i1,i2) + instU(i1+3,i4)*Nx(i4,i2)
					enddo
				enddo
				i4 = dofTable(1,i5) - 3
			endif
			
			ptRot(:) = r_0
			do i1 = 1, 3
				do i2 = 1, numNds
					ptRot(i1) = ptRot(i1) + instU(i1+3,i2)*Nvec(i2)
				enddo
			enddo
			
			call r_getShellDef(dsDdU(:,i5),ux,ptRot,rx)
			dfMdU(:,i5) = r_0
			do i2 = 1, 9
				do i3 = 1, 9
					dfMdU(i2,i5) = dfMdU(i2,i5) + ABD(i2,i3)*dsDdU(i3,i5)
				enddo
			enddo	
		enddo
        
        dfMdT(:,:) = r_0
		dPtTdT(:) = r_0
		do i1 = 1, numNds
		    dfMdT(1:6,i1) = -Nvec(i1)*sTELd(1:6)
			dPtTdT(i1) = Nvec(i1)
		enddo
		
	end subroutine r_getdShellFrcMomdU
	
	subroutine r_getBeamFrcMom(frcMom,bmDef,elNum,sPt)
	    implicit none
		
		integer, intent(in) :: elNum
		real*8, intent(in) :: sPt(3)
		real*8, intent(out) :: frcMom(6), bmDef(6)
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
        
        real*8 :: ux(3,3), rx(3,3), Nx(11,3), Nvec(11), detJ, PtT
		real*8 :: instOri(3,240), drIdrG(3,16), rot(3), ptRot(3)
		real*8 :: instU(6,11)
		real*8 :: rnNds
		integer :: i1, i2, i3
		
		eType = elementType(elNum)
		if(eType .ne. 2) then
		    return
		endif
		call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
			locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
			ABD, sMass, sTELd, stCond, ssHeat, &
			bStiff, bMass, bTELd, btCond, bsHeat, &
			temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,elNum)
		
		call r_getIntPtData(Nvec,Nx,detJ,locNds,sPt,eType,numNds)
		
		call r_getInstOrient(instOri,orient,disp,numNds,0)
		call r_getInstDof(instU,disp,globNds,instOri,orient,numNds,dofTable,0,0)
		
		ux(:,:) = r_0
		rx(:,:) = r_0
		do i1 = 1, 3
		    do i2 = 1, 3
			    do i3 = 1, 10
				    ux(i1,i2) = ux(i1,i2) + instU(i1,i3)*Nx(i3,i2)
					rx(i1,i2) = rx(i1,i2) + instU(i1+3,i3)*Nx(i3,i2)
				enddo
			enddo
		enddo
		
		ptRot(:) = r_0
		do i1 = 1, 3
		    do i2 = 1, numNds
			    ptRot(i1) = ptRot(i1) + instU(i1+3,i2)*Nvec(i2)
			enddo
		enddo
		
		call r_getBeamDef(bmDef,ux,ptRot,rx)
		ptT = r_0
		do i2 = 1, numNds
			ptT = ptT + Nvec(i2)*temp(i2)
		enddo
		frcMom(:) = r_0
		do i2 = 1, 6
			do i3 = 1, 6
				frcMom(i2) = frcMom(i2) + bStiff(i2,i3)*bmDef(i3)
			enddo
		enddo
		frcMom(1:6) = frcMom(1:6) - ptT*bTELd(1:6)
	
	end subroutine r_getBeamFrcMom
	
	subroutine r_getdBeamFrcMomdU(dfMdU,dfMdT,dbDdU,elNum,sPt)
	    implicit none
		
		integer, intent(in) :: elNum
		real*8, intent(in) :: sPt(3)
		real*8, intent(out) :: dfMdU(6,33), dfMdT(6,10), dbDdU(6,33)
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
        
        real*8 :: ux(3,3), rx(3,3), Nx(11,3), Nvec(11), detJ, dPtTdT(10)
		real*8 :: instOri(3,240), drIdrG(3,16), rot(3), ptRot(3)
		real*8 :: instU(6,11)
		real*8 :: rnNds
		integer :: i1, i2, i3, i4, i5
		
		eType = elementType(elNum)
		if(eType .ne. 2) then
		    return
		endif
		call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
			locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
			ABD, sMass, sTELd, stCond, ssHeat, &
			bStiff, bMass, bTELd, btCond, bsHeat, &
			temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,elNum)
		
		call r_getIntPtData(Nvec,Nx,detJ,locNds,sPt,eType,numNds)
		
		call r_getInstOrient(instOri,orient,disp,numNds,0)
		i1 = numNds*dofPerNd + numIntDof
		do i5 = 1, i1
			call r_getInstDof(instU,disp,globNds,instOri,orient,numNds,dofTable,i5,0)
			
			ux(:,:) = r_0
			rx(:,:) = r_0
			if(dofTable(1,i5) .le. 3) then
			    i3 = dofTable(2,i5)
				do i1 = 1, 3
					do i2 = 1, 3
						ux(i1,i2) = ux(i1,i2) + instU(i1,i3)*Nx(i3,i2)
					enddo
				enddo
			else
			    i4 = dofTable(2,i5)
			    do i1 = 1, 3
					do i2 = 1, 3
						do i3 = 1, 10
							ux(i1,i2) = ux(i1,i2) + instU(i1,i3)*Nx(i3,i2)
						enddo
						rx(i1,i2) = rx(i1,i2) + instU(i1+3,i4)*Nx(i4,i2)
					enddo
				enddo
				i4 = dofTable(1,i5) - 3
			endif
			
			ptRot(:) = r_0
			do i1 = 1, 3
				do i2 = 1, numNds
					ptRot(i1) = ptRot(i1) + instU(i1+3,i2)*Nvec(i2)
				enddo
			enddo
			
			call r_getBeamDef(dbDdU(:,i5),ux,ptRot,rx)
			dfMdU(:,i5) = r_0
			do i2 = 1, 6
				do i3 = 1, 6
					dfMdU(i2,i5) = dfMdU(i2,i5) + bStiff(i2,i3)*dbDdU(i3,i5)
				enddo
			enddo	
		enddo
        
        dfMdT(:,:) = r_0	
		do i1 = 1, numNds
		    dfMdT(1:6,i1) = -Nvec(i1)*bTELd(1:6)
		enddo
		
	end subroutine r_getdBeamFrcMomdU
	
	subroutine r_getElMassVol(mass,vol,elNum)
	    implicit none
		
		integer, intent(in) :: elNum
		real*8, intent(out) :: mass, vol
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat

        real*8 :: Nvec(11), Nx(11,3), detJ, thick, zOff, area
        integer :: i1, i2, i3, i4

        eType = elementType(elNum)
        call r_getElementProfile(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt,elNum)
		call r_getElementProperties(locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, &
	             ABD, sMass, sTELd, stCond, ssHeat, bStiff, bMass, bTELd, btCond, bsHeat, elNum)
	    
		if(eType .eq. 41 .or. eType .eq. 3) then
		    call r_getElThickOffset(thick,zOff,elNum)
		endif
		
		if(eType .eq. 2) then
		    i4 = elementSection(elNum)
		    area = beamProperties(1,i4)
		    do i1 = elToDRange(elNum-1)+1, elToDRange(elNum)
			    i2 = elToD(i1)
				if(dCategory(i2) .eq. 'area') then
					area = area + elToCoef(i1)*r_dVec(i2)
				endif
			enddo
		endif
		
		mass = r_0
		vol = r_0
		do i1 = 1, numIntPts
		    call r_getIntPtData(Nvec,Nx,detJ,locNds,intPts(:,i1),eType,numNds)
			if(eType .eq. 41 .or. eType .eq. 3) then
			    vol = vol + thick*detJ*ipWt(i1)
			    mass = mass + sMass(1,1)*detJ*ipWt(i1)
			elseif(eType .eq. 2) then
			    vol = vol + area*detJ*ipWt(i1)
			    mass = mass + bMass(1,1)*detJ*ipWt(i1)
			else
				vol = vol + detJ*ipWt(i1)
				mass = mass + den*detJ*ipWt(i1)
			endif
		enddo
		
	end subroutine r_getElMassVol
	
	subroutine r_getLayerStress(stress,strain,shDef,ptT,elNum,layerNum)
	    implicit none

        integer, intent(in) :: elNum, layerNum
		real*8, intent(in) :: shDef(9), ptT
		real*8, intent(out) :: stress(6), strain(6)

        real*8 :: secThick, zOff, zCrd, layerThick, zNext, zMid
		real*8 :: globE(3), locE(3), locS(3)
		real*8 :: layerProps(6), tExp(6), SMat(3,3), QMat(3,3), fVec(3), zVec(3)
		real*8 :: thetaRad, a11, a22, a12, a21, Te(3,3), Ts(3,3)
		integer :: secNum, i1, i2, i3, i4, matID
		
		secNum = elementSection(elNum)
		
		call r_getElThickOffset(secThick,zOff,elNum)
		
		zCrd = -r_p5*secThick*(zOff+r_1)  !! initialize to zMin
		i4 = 0
		do i1 = secLayupRange(secNum-1)+1, secLayupRange(secNum-1)+layerNum
		    i4 = i4 + 1
			matID = layupMatId(i1)
			layerThick = r_1*layupThickness(i1)
			do i2 = elToDRange(elNum-1)+1, elToDRange(elNum)
			    i3 = elToD(i2)
			    if(dCategory(i3) .eq. 'thickness' .and. dLayer(i3) .eq. i4) then
					layerThick = layerThick + elToCoef(i2)*r_dVec(i3)
				endif
			enddo
			zNext = zCrd + layerThick
			zMid = r_p5*(zCrd + zNext)
			zCrd = zNext
		enddo

        globE(1) = shDef(1) + zMid*shDef(4)
		globE(2) = shDef(2) + zMid*shDef(5)
		globE(3) = shDef(3) + zMid*shDef(6)
		
		i1 = secLayupRange(secNum-1)+layerNum
		layerProps(1) = r_1*materialElastic(1,matID)
		layerProps(2) = r_1*materialElastic(2,matID)
		layerProps(3) = r_1*materialElastic(4,matID)
		layerProps(4) = r_1*materialElastic(7,matID)
		tExp(:) = r_1*materialThermExp(:,matID)
		layerProps(6) = r_1*layupAngle(i1)
		do i2 = elToDRange(elNum-1)+1, elToDRange(elNum)
			i3 = elToD(i2)
			if(dLayer(i3) .eq. layerNum) then
				if(dCategory(i3) .eq. 'modulus') then
					i4 = dComponent(i3)
					layerProps(i4) = layerProps(i4) + elToCoef(i2)*r_dVec(i3)
				elseif(dCategory(i3) .eq. 'poissonRatio') then
					layerProps(3) = layerProps(3) + elToCoef(i2)*r_dVec(i3)
				elseif(dCategory(i3) .eq. 'shearModulus') then
					layerProps(4) = layerProps(4) + elToCoef(i2)*r_dVec(i3)
				elseif(dCategory(i3) .eq. 'thermalExp') then
				    i4 = dComponent(i3)
					tExp(i4) = tExp(i4) + elToCoef(i2)*r_dVec(i3)
				elseif(dCategory(i3) .eq. 'angle') then
					layerProps(6) = layerProps(6) + elToCoef(i2)*r_dVec(i3)
				endif
			endif
		enddo
		SMat(:,:) = r_0
		SMat(1,1) = r_1/layerProps(1)
		SMat(1,2) = -layerProps(3)/layerProps(1)
		SMat(2,2) = r_1/layerProps(2)
		SMat(3,3) = r_1/layerProps(4)
		SMat(2,1) = SMat(1,2)
		call r_GetInvLU(QMat,SMat,fVec,zVec,3,1,3)
		thetaRad = r_pi180*layerProps(6)
		call r_cos(a11,thetaRad)
		call r_sin(a12,thetaRad)
		a21 = -a12
		a22 = a11
		Ts(1,1) = a11*a11
		Ts(1,2) = a12*a12
		Ts(1,3) = r_2*a11*a12
		Ts(2,1) = a21*a21
		Ts(2,2) = a22*a22
		Ts(2,3) = r_2*a22*a21
		Ts(3,1) = a11*a21
		Ts(3,2) = a12*a22
		Ts(3,3) = a11*a22 + a12*a21
		Te(:,:) = Ts(:,:)
		Te(:,3) = r_p5*Te(:,3)
		Te(3,:) = r_2*Te(3,:)

        locE(:) = r_0
        do i1 = 1, 3
		    do i2 = 1, 3
			    locE(i1) = locE(i1) + Te(i1,i2)*globE(i2)
			enddo
        enddo
		
		tExp(3) = tExp(4)
		locS(:) = r_0
		do i1 = 1, 3
		    do i2 = 1, 3
			    locS(i1) = locS(i1) + QMat(i1,i2)*(locE(i2) - ptT*tExp(i2))
			enddo
        enddo
		
		stress(:) = r_0
		stress(1) = locS(1)
		stress(2) = locS(2)
		stress(4) = locS(3)
		
		strain(:) = r_0
		strain(1) = locE(1)
		strain(2) = locE(2)
		strain(4) = locE(3)
		
		
	end subroutine r_getLayerStress
	
	subroutine r_getdLayerStressdU(dsdU,dsdT,dedU,dsDdU,dPtTdT,elNum,layer)
	    implicit none
		
	    integer, intent(in) :: elNum, layer
		real*8, intent(in) :: dsDdU(9,33), dPtTdT(10)
		real*8, intent(out) :: dsdU(6,33), dsdT(6,10), dedU(6,33)
		
		real*8 :: secThick, zOff, zCrd, layerThick, zNext, zMid
		real*8 :: dglobEdU(3,33), dlocEdU(3,33), dlocSdU(3,33), dlocSdT(3,10)
		real*8 :: layerProps(6), tExp(6), SMat(3,3), QMat(3,3), fVec(3), zVec(3)
		real*8 :: thetaRad, a11, a22, a12, a21, Te(3,3), Ts(3,3)
		integer :: secNum, i1, i2, i3, i4, matID
		
		secNum = elementSection(elNum)
		
		call r_getElThickOffset(secThick,zOff,elNum)
		
		zCrd = -r_p5*secThick*(zOff+r_1)  !! initialize to zMin
		i4 = 0
		do i1 = secLayupRange(secNum-1)+1, secLayupRange(secNum-1)+layer
		    i4 = i4 + 1
			matID = layupMatId(i1)
			layerThick = r_1*layupThickness(i1)
			do i2 = elToDRange(elNum-1)+1, elToDRange(elNum)
			    i3 = elToD(i2)
			    if(dLayer(i3) .eq. i4) then
					if(dCategory(i3) .eq. 'thickness') then
					    layerThick = layerThick + elToCoef(i2)*r_dVec(i3)
					endif
				endif
			enddo
			zNext = zCrd + layerThick
			zMid = r_p5*(zCrd + zNext)
			zCrd = zNext
		enddo

        dglobEdU(1,:) = dsDdU(1,:) + zMid*dsDdU(4,:)
		dglobEdU(2,:) = dsDdU(2,:) + zMid*dsDdU(5,:)
		dglobEdU(3,:) = dsDdU(3,:) + zMid*dsDdU(6,:)
		
		i1 = secLayupRange(secNum-1)+layer
		layerProps(1) = r_1*materialElastic(1,matID)
		layerProps(2) = r_1*materialElastic(2,matID)
		layerProps(3) = r_1*materialElastic(4,matID)
		layerProps(4) = r_1*materialElastic(7,matID)
		tExp(:) = r_1*materialThermExp(:,matID)
		layerProps(6) = r_1*layupAngle(i1)
		do i2 = elToDRange(elNum-1)+1, elToDRange(elNum)
			i3 = elToD(i2)
			if(dLayer(i3) .eq. layer) then
				if(dCategory(i3) .eq. 'modulus') then
					i4 = dComponent(i3)
					layerProps(i4) = layerProps(i4) + elToCoef(i2)*r_dVec(i3)
				elseif(dCategory(i3) .eq. 'poissonRatio') then
					layerProps(3) = layerProps(3) + elToCoef(i2)*r_dVec(i3)
				elseif(dCategory(i3) .eq. 'shearModulus') then
					layerProps(4) = layerProps(4) + elToCoef(i2)*r_dVec(i3)
				elseif(dCategory(i3) .eq. 'thermalExp') then
				    i4 = dComponent(i3)
					tExp(i4) = tExp(i4) + elToCoef(i2)*r_dVec(i3)
				elseif(dCategory(i3) .eq. 'angle') then
					layerProps(6) = layerProps(6) + elToCoef(i2)*r_dVec(i3)
				endif
			endif
		enddo
		SMat(:,:) = r_0
		SMat(1,1) = r_1/layerProps(1)
		SMat(1,2) = -layerProps(3)/layerProps(1)
		SMat(2,2) = r_1/layerProps(2)
		SMat(3,3) = r_1/layerProps(4)
		SMat(2,1) = SMat(1,2)
		call r_GetInvLU(QMat,SMat,fVec,zVec,3,1,3)
		thetaRad = r_pi180*layerProps(6)
		call r_cos(a11,thetaRad)
		call r_sin(a12,thetaRad)
		a21 = -a12
		a22 = a11
		Ts(1,1) = a11*a11
		Ts(1,2) = a12*a12
		Ts(1,3) = r_2*a11*a12
		Ts(2,1) = a21*a21
		Ts(2,2) = a22*a22
		Ts(2,3) = r_2*a22*a21
		Ts(3,1) = a11*a21
		Ts(3,2) = a12*a22
		Ts(3,3) = a11*a22 + a12*a21
		Te(:,:) = Ts(:,:)
		Te(:,3) = r_p5*Te(:,3)
		Te(3,:) = r_2*Te(3,:)

        dlocEdU(:,:) = r_0
        do i1 = 1, 3
		    do i2 = 1, 3
			    dlocEdU(i1,:) = dlocEdU(i1,:) + Te(i1,i2)*dglobEdU(i2,:)
			enddo
        enddo
		
		tExp(3) = tExp(4)
		dlocSdU(:,:) = r_0
		do i1 = 1, 3
		    do i2 = 1, 3
			    dlocSdU(i1,:) = dlocSdU(i1,:) + QMat(i1,i2)*dlocEdU(i2,:)
			enddo
        enddo
		
		dlocSdT(:,:) = r_0
		do i1 = 1, 3
		    do i2 = 1, 3
			    dlocSdT(i1,:) = dlocSdT(i1,:) - QMat(i1,i2)*dPtTdT(:)*tExp(i2)
			enddo
        enddo
		
		dsdU(:,:) = r_0
		dsdU(1,:) = dlocSdU(1,:)
		dsdU(2,:) = dlocSdU(2,:)
		dsdU(4,:) = dlocSdU(3,:)
		
		dsdT(:,:) = r_0
		dsdT(1,:) = dlocSdT(1,:)
		dsdT(2,:) = dlocSdT(2,:)
		dsdT(4,:) = dlocSdT(3,:)
		
		dedU(:,:) = r_0
		dedU(1,:) = dlocEdU(1,:)
		dedU(2,:) = dlocEdU(2,:)
		dedU(4,:) = dlocEdU(3,:)
		
	end subroutine r_getdLayerStressdU
	
	subroutine r_getElStress(stress,strain,elNum,layer,sPt)
	    implicit none
		
		integer, intent(in) :: elNum, layer
		real*8, intent(in) :: sPt(3)
		real*8, intent(out) :: stress(6), strain(6)
		
		real*8 :: shDef(9), frcMom(9), ptT
		
		if(layer .ne. 0) then
		    call r_getShellFrcMom(frcMom,shDef,ptT,elNum,sPt)
		    call r_getLayerStress(stress,strain,shDef,ptT,elNum,layer)
		else
		    call r_getSolidStress(stress,strain,elNum,sPt)
		endif
		
	end subroutine r_getElStress
	
	subroutine r_getdElStressdU(dsdU,dsdT,dedU,elNum,layer,sPt)
	    implicit none
		
		integer, intent(in) :: elNum, layer
		real*8, intent(in) :: sPt(3)
		real*8, intent(out) :: dsdU(6,33), dsdT(6,10), dedU(6,33)
		
		real*8 :: dfMdU(9,33), dfMdT(9,10), dsDdU(9,33), dPtTdT(10)
		
		if(layer .ne. 0) then
		    call r_getdShellFrcMomdU(dfMdU,dfMdT,dsDdU,dPtTdT,elNum,sPt)
			call r_getdLayerStressdU(dsdU,dsdT,dedU,dsDdU,dPtTdT,elNum,layer)
		else
		    call r_getdSolidStressdU(dsdU,dsdT,dedU,elNum,sPt)
		endif
		
	end subroutine r_getdElStressdU
	
	subroutine r_getElStrainEnergy(stEn,dSEdU,dSEdT,elNum)
	    implicit none
		
		integer, intent(in) :: elNum
		real*8, intent(out) :: stEn, dSEdU(33), dSEdT(10)
		
		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), stExp(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), btExp(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)

        real*8 :: dRdU(33,33), dRdT(33,10)
        integer :: i1, i2, i3, i4
		
		call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
	    locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
		ABD, sMass, stExp, stCond, ssHeat, &
		bStiff, bMass, btExp, btCond, bsHeat, &
	    temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,elNum)

        eType = elementType(elNum)		
		call r_getElRuk(dSEdU,dRdU,dRdT,stEn,dSEdT,0,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
	    locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,disp)
		
	end subroutine r_getElStrainEnergy
	
	subroutine r_getElFlux(flux,dfluxdT,Tx,dTxdT,elNum,sPt)
	    implicit none
		
		integer, intent(in) :: elNum
		real*8, intent(in) :: sPt(3)
		real*8, intent(out) :: flux(3), dfluxdT(3,10), Tx(3), dTxdT(3,10)

		integer :: eType
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		integer :: dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), sTELd(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), bTELd(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)

        real*8 :: thisCond(3,3), Nvec(11), Nx(11,3), detJ
        integer :: i1, i2, i3, i4
		
		call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
	    locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
		ABD, sMass, sTELd, stCond, ssHeat, &
		bStiff, bMass, bTELd, btCond, bsHeat, &
	    temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,elNum)        

        eType = elementType(elNum)
		if(eType .eq. 2) then
		    thisCond = btCond
		elseif(eType .eq. 41 .or. eType .eq. 3) then
		    thisCond = stCond
		else
		    thisCond = tCond
		endif
		
		call r_getIntPtData(Nvec,Nx,detJ,locNds,sPt,eType,numNds)
		Tx(:) = r_0
		dTxdT(:,:) = r_0
		do i2 = 1, numNds
			do i3 = 1, 3
				Tx(i3) = Tx(i3) + Nx(i2,i3)*temp(i2)
				dTxdT(i3,i2) = Nx(i2,i3)
			enddo
		enddo
		flux = r_0
		dfluxdT(:,:) = r_0
		do i2 = 1, 3
			do i3 = 1, 3
				flux(i2) = flux(i2) - thisCond(i2,i3)*Tx(i3)
				dfluxdT(i2,:) = dfluxdT(i2,:) - thisCond(i2,i3)*dTxdT(i3,:)
			enddo
		enddo
		
	end subroutine r_getElFlux
	
	subroutine r_getElBodyForce(ndFrc,dofTable,ndDof,bdyFrc,isGrav,elNum)
	    implicit none
		
		integer, intent(in) :: isGrav, elNum
		real*8, intent(in) :: bdyFrc(6)
		real*8, intent(out) :: ndFrc(33)
		integer, intent(out) :: dofTable(2,33), ndDof
		
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), stExp(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), btExp(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
		real*8 :: statInOri(3,240)		

        real*8 :: dRdA(33,33), fakeAcc(6,11)
		real*8 :: Nvec(11), Nx(11,3), detJ, nnInv
		real*8 :: thick, zOff, props(7)
		integer :: i1, i2, i3, i4, i5, eType
		
		eType = elementType(elNum)
		
		call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, stExp, stCond, ssHeat, &
				bStiff, bMass, btExp, btCond, bsHeat, &
				temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,elNum)
		
		ndDof = numNds*dofPerNd
		if(isGrav .eq. 1) then
			if(eType .eq. 41 .or. eType .eq. 3 .or. eType .eq. 2) then
			    call r_getInstOrient(statInOri,orient,disp,numNds,1)
			endif
		    do i1 = 1, numNds
			    fakeAcc(:,i1) = bdyFrc(:)
			enddo
			call r_getElRum(ndFrc,dRdA,0,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
	            locNds, globNds, orient, statInOri, den, sMass, bMass, disp, fakeAcc)
		else
		    ndFrc(:) = r_0
			do i1 = 1, numIntPts
			    call r_getIntPtData(Nvec,Nx,detJ,locNds,intPts(:,i1),eType,numNds)
				i3 = numNds*dofPerNd
				detJ = detJ*ipWt(i1)
				do i2 = 1, i3
				    i4 = dofTable(1,i2)
					i5 = dofTable(2,i2)
					ndFrc(i2) = ndFrc(i2) + bdyFrc(i4)*Nvec(i5)*detJ
				enddo
			enddo
			if(eType .eq. 41 .or. eType .eq. 3) then
			    call r_getElThickOffset(thick,zOff,elNum)
				ndFrc(:) = thick*ndFrc(:)
			elseif(eType .eq. 2) then
			    call r_getBeamSecProps(props,elNum)
			    ndFrc(:) = props(1)*ndFrc(:)
			endif
		endif
	
	end subroutine r_getElBodyForce
	
	subroutine r_getElCentAcc(centAcc,center,axis,angVel,elNum)
	    implicit none
		
		integer, intent(in) :: elNum
		real*8, intent(in) :: center(3), axis(3), angVel
		real*8, intent(out) :: centAcc(3)
		
		integer :: numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), stExp(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), btExp(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
		real*8 :: statInOri(3,240), statdrIdrG(3,12), statRot(3)		

        real*8 :: elCent(3), nnInv, posVec(3), dp, mag, unitAxis(3)
		integer :: i1, i2, i3, i4, i5, eType
		
		eType = elementType(elNum)
		
		call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, stExp, stCond, ssHeat, &
				bStiff, bMass, btExp, btCond, bsHeat, &
				temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,elNum)

        elCent(:) = r_0
		
		do i1 = 1, numNds
		    do i2 = 1, 3
		        elCent(i2) = elCent(i2) + globNds(i2,i1) + disp(i2,i1)
			enddo
		enddo
		nnInv = numNds
		nnInv = r_1/nnInv
		elCent(:) = nnInv*elCent(:)
		
		posVec(:) = elCent(:) - center(:)
		
		dp = axis(1)*axis(1) + axis(2)*axis(2) + axis(3)*axis(3)
		call r_sqrt(mag,dp)
		unitAxis(:) = (r_1/mag)*axis(:)
		
		dp = posVec(1)*unitAxis(1) + posVec(2)*unitAxis(2) + posVec(3)*unitAxis(3)
		
		posVec(:) = posVec(:) - dp*unitAxis(:)
		
		centAcc(:) = angVel*angVel*posVec(:)
	
	end subroutine r_getElCentAcc
	
	subroutine r_getElSurfaceTraction(ndFrc,dofTable,ndDof,trac,pressure,elNum,normDir,thetaTol)
	    implicit none
		
		integer, intent(in) :: elNum
        real*8, intent(in) :: trac(6), pressure, normDir(3), thetaTol
		integer, intent(out) :: dofTable(2,33), ndDof
		real*8, intent(out) :: ndFrc(33)
		
		integer :: numNds, dofPerNd, numIntDof, numIntPts
		real*8 :: intPts(3,8), ipWt(8)
		real*8 :: locNds(3,10), globNds(3,10), orient(3,3), cMat(6,6), den, tExp(6), tCond(3,3), sHeat
		real*8 :: ABD(9,9), sMass(6,6), stExp(6), stCond(3,3), ssHeat
		real*8 :: bStiff(6,6), bMass(6,6), btExp(6), btCond(3,3), bsHeat
		real*8 :: temp(10), Tdot(10), pTemp(10), pTdot(10)
		real*8 :: disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)		

        real*8 :: mag, dp, thetaRad, cosTR, tracFinal(6)
        real*8 :: Nvec(11), Nx(11,3), detJ, faceNds(3,10), elNorm(3), unitND(3)
		integer :: i1, i2, i3, i4, i5, i6, i7, eType, faces(6,6), numFaces, faceNum, faceNumNds, gt
		
		ndFrc(:) = r_0
		call r_getElementData(numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds, globNds, orient, cMat, den, tExp, tCond, sHeat, & 
				ABD, sMass, stExp, stCond, ssHeat, &
				bStiff, bMass, btExp, btCond, bsHeat, &
				temp,Tdot,disp,vel,acc,pTemp,pTdot,pDisp,pVel,pAcc,elNum)
		
		ndDof = numNds*dofPerNd
		
		dp = normDir(1)*normDir(1) + normDir(2)*normDir(2) + normDir(3)*normDir(3)
		call r_sqrt(mag,dp)
		unitND(:) = (r_1/mag)*normDir(:)
		
		thetaRad = r_pi180*thetaTol
		call r_cos(cosTR,thetaRad)
		
		call r_getElementFaces(faces,numFaces,elNum)
		
		do faceNum = 1, numFaces
		    if(elementSurfaces(faceNum,elNum) .eq. 1) then
				call r_getElSurfaceNormal(elNorm,elNum,faceNum,globNds)
				
				dp = elNorm(1)*unitND(1) + elNorm(2)*unitND(2) + elNorm(3)*unitND(3)
				
				call r_greater(gt,dp,cosTR)
				if(gt .eq. 1) then
					faceNumNds = 0
					faceNds(:,:) = r_0
					do i1 = 1, 6
						i2 = faces(i1,faceNum)
						if(i2 .ne. 0) then
							faceNds(:,i1) = locNds(:,i2)
							faceNumNds = faceNumNds + 1
						endif
					enddo
					
					if(abs(pressure) .gt. 0d0) then
						tracFinal(1:3) = -pressure*elNorm(1:3)
						tracFinal(4:6) = r_0
					else
						tracFinal(:) = trac(:)
					endif
					
					if(faceNumNds .eq. 3) then
						intPts(:,1) = r_1o6*(/r_1,r_1,r_0/)
						intPts(:,2) = r_1o6*(/4d0,1d0,0d0/)
						intPts(:,3) = r_1o6*(/1d0,4d0,0d0/)
						ipWt(1:3) = r_1o6
						do i1 = 1, 3
							call r_getIntPtData(Nvec,Nx,detJ,faceNds,intPts(:,i1),3,3)
							i3 = numNds*dofPerNd
							detJ = detJ*ipWt(i1)
							do i2 = 1, i3
								i4 = dofTable(1,i2)
								i5 = dofTable(2,i2)
								do i6 = 1, 3
								    i7 = faces(i6,faceNum)
								    if(i5 .eq. i7) then
								        ndFrc(i2) = ndFrc(i2) + tracFinal(i4)*Nvec(i6)*detJ
									endif
								enddo
							enddo
						enddo
					elseif(faceNumNds .eq. 4) then
						intPts(:,1) = r_1rt3*(/-r_1,-r_1,r_0/)
						intPts(:,2) = r_1rt3*(/r_1,-r_1,r_0/)
						intPts(:,3) = r_1rt3*(/-r_1,r_1,r_0/)
						intPts(:,4) = r_1rt3*(/r_1,r_1,r_0/)
						ipWt(1:4) = r_1
						do i1 = 1, 4
							call r_getIntPtData(Nvec,Nx,detJ,faceNds,intPts(:,i1),41,4)
							i3 = numNds*dofPerNd
							detJ = detJ*ipWt(i1)
							do i2 = 1, i3
								i4 = dofTable(1,i2)
								i5 = dofTable(2,i2)
								do i6 = 1, 4
								    i7 = faces(i6,faceNum)
								    if(i5 .eq. i7) then
								        ndFrc(i2) = ndFrc(i2) + tracFinal(i4)*Nvec(i6)*detJ
									endif
								enddo
							enddo
						enddo	
					endif
					
				endif
			endif
		enddo
	
	end subroutine r_getElSurfaceTraction
	
	subroutine r_getElSurfaceNormal(nVec,elNum,faceNum,globNds)
	    implicit none
		
		integer, intent(in) :: elNum, faceNum
		real*8, intent(out) :: globNds(3,10)
		real*8, intent(out) :: nVec(3)
		
		real*8 :: v1(3), v2(3), v3(3), mag, mag2
		integer :: faces(6,6), numFaces
		integer :: n1, n2, n3, n4
		
		call r_getElementFaces(faces,numFaces,elNum)
		
		n1 = faces(1,faceNum)
		n2 = faces(2,faceNum)
		n3 = faces(3,faceNum)
		n4 = faces(4,faceNum)
		
		if(n4 .ne. 0) then
		    v1(:) = globNds(:,n3) - globNds(:,n1)
			v2(:) = globNds(:,n4) - globNds(:,n2)
			v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
			v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
			v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
			mag2 = v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3)
			call r_sqrt(mag,mag2)
			nVec(:) = (r_1/mag)*v3(:)
		elseif(n3 .ne. 0) then
			v1(:) = globNds(:,n2) - globNds(:,n1)
			v2(:) = globNds(:,n3) - globNds(:,n2)
			v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
			v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
			v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
			mag2 = v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3)
			call r_sqrt(mag,mag2)
			nVec(:) = (r_1/mag)*v3(:)
		else
		    nVec(:) = r_0
		endif
		
	end subroutine r_getElSurfaceNormal
	
	subroutine r_getElRtk(Rtk,dRdt,buildMat,numNds,numIntPts,intPts,ipWt, &
	    locNds,orient,tCond,stCond,btCond,temp,eType)
	    implicit none
		
		integer, intent(in) :: buildMat, numNds, numIntPts, eType
		real*8, intent(in) :: intPts(3,8), ipWt(8), locNds(3,10), orient(3,3)
		real*8, intent(in) :: tCond(3,3), stCond(3,3), btCond(3,3), temp(10)
		real*8, intent(out) :: Rtk(10), dRdt(10,10)

        real*8 :: Nvec(11), Nx(11,3), kNx(11,3), detJ, qVec(3), Tx(3), thisCond(3,3)
        integer :: i1, i2, i3, i4
		
		if(eType .eq. 2) then
		    thisCond = btCond
		elseif(eType .eq. 41 .or. eType .eq. 3) then
		    thisCond = stCond
		else
		    thisCond = tCond
		endif
		
		Rtk(:) = r_0
		dRdt(:,:) = r_0
		do i1 = 1, numIntPts
		    call r_getIntPtData(Nvec,Nx,detJ,locNds,intPts(:,i1),eType,numNds)
			Tx(:) = r_0
			do i2 = 1, numNds
			    do i3 = 1, 3
				    Tx(i3) = Tx(i3) + Nx(i2,i3)*temp(i2)
				enddo
			enddo
			qVec = r_0
			do i2 = 1, 3
			    do i3 = 1, 3
				    qVec(i2) = qVec(i2) + thisCond(i2,i3)*Tx(i3)
				enddo
			enddo
			qVec(:) = (detJ*ipWt(i1))*qVec(:)
			do i2 = 1, numNds
			    do i3 = 1, 3
				    Rtk(i2) = Rtk(i2) + Nx(i2,i3)*qVec(i3)
				enddo
			enddo
			if(buildMat .eq. 1) then
			    kNx(:,:) = r_0
				do i2 = 1, numNds
				    do i3 = 1, 3
					    do i4 = 1, 3
						    kNx(i2,i3) = kNx(i2,i3) + Nx(i2,i4)*thisCond(i4,i3)
						enddo
					enddo
				enddo
				kNx = (detJ*ipWt(i1))*kNx
				do i2 = 1, numNds
				    do i3 = 1, numNds
					    do i4 = 1, 3
						    dRdt(i2,i3) = dRdt(i2,i3) + Nx(i2,i4)*kNx(i3,i4)
						enddo
					enddo
				enddo
			endif
		enddo
		
	end subroutine r_getElRtk
	
	subroutine r_getElRtm(Rtm,dRdTdot,buildMat,numNds,numIntPts,intPts,ipWt, &
	    locNds, den, sHeat, ssHeat, bsHeat, Tdot, eType)
	    implicit none
		
		integer, intent(in) :: buildMat, numNds, numIntPts, eType
		real*8, intent(in) :: intPts(3,8), ipWt(8), locNds(3,10), den, sHeat, ssHeat, bsHeat, Tdot(10)
		real*8, intent(out) :: Rtm(10), dRdTdot(10,10)

        real*8 :: thisSH, ptTdot, Nvec(11), Nx(11,3), cpN(11), detJ		
		integer :: i1, i2, i3
		
		if(eType .eq. 2) then
		    thisSH = bsHeat
		elseif(eType .eq. 41 .or. eType .eq. 3) then
		    thisSH = ssHeat
		else
		    thisSH = den*sHeat
		endif
		
		Rtm(:) = r_0
		dRdTdot(:,:) = r_0
		do i1 = 1, numIntPts
		    call r_getIntPtData(Nvec,Nx,detJ,locNds,intPts(:,i1),eType,numNds)
			ptTdot = r_0
			do i2 = 1, numNds
			    ptTdot = ptTdot + Tdot(i2)*Nvec(i2)
			enddo
			ptTdot = detJ*ipWt(i1)*thisSH*ptTdot
			Rtm(1:numNds) = Rtm(1:numNds) + ptTdot*Nvec(1:numNds)
			if(buildMat .eq. 1) then
			    cpN = detJ*ipWt(i1)*thisSH*Nvec
				do i2 = 1, numNds
				    do i3 = 1, numNds
					    dRdTdot(i2,i3) = dRdTdot(i2,i3) + Nvec(i2)*cpN(i3)
					enddo
				enddo
			endif
		enddo
		
	end subroutine r_getElRtm
	
	subroutine r_getElRt(Rt,dRdt,buildMat,numNds,numIntPts,intPts,ipWt, &
	    locNds,orient,tCond,stCond,btCond,den,sHeat,ssHeat,bsHeat, &
		temp,pTemp,pTDot,eType)
	    implicit none
		
		integer, intent(in) :: buildMat, numNds, numIntPts, eType
		real*8, intent(in) :: intPts(3,8), ipWt(8), locNds(3,10), orient(3,3)
		real*8, intent(in) :: tCond(3,3), stCond(3,3), btCond(3,3)
		real*8, intent(in) :: den, sHeat, ssHeat, bsHeat
		real*8, intent(in) :: temp(10), pTemp(10), pTDot(10)
		real*8, intent(out) :: Rt(10), dRdt(10,10)
		
		real*8 :: Rtk(10), Rtm(10), dRdVar(10,10), Tdot(10)
		integer :: gt1, i1

        Rtk(:) = r_0
        call r_getElRtk(Rtk,dRdVar,buildMat,numNds,numIntPts,intPts,ipWt, &
	    locNds,orient,tCond,stCond,btCond,temp,eType)
		
		Rt(:) = Rtk(:)
		if(buildMat .eq. 1) then
		    dRdt(:,:) = dRdVar(:,:)
		endif
		
		if(dynamic .eq. 1) then
		    do i1 = 1, numNds
		        Tdot(i1) = (r_1/nMGamma)*((r_1/delT)*(temp(i1)-pTemp(i1)) - r_1*(1d0-nMGamma)*pTdot(i1))
			enddo
			Rtm(:) = r_0
		    call r_getElRtm(Rtm,dRdVar,buildMat,numNds,numIntPts,intPts,ipWt, &
	            locNds,den,sHeat,ssHeat,bsHeat,Tdot,eType)
			Rt(:) = Rt(:) + Rtm(:)
			if(buildMat .eq. 1) then
			    dRdt(:,:) = dRdt(:,:) + (r_1/(nMGamma*delT))*dRdVar(:,:)
			endif
		endif
		
	end subroutine r_getElRt
	
	subroutine r_getElRuk(Ruk,dRdU,dRdT,stEn,dSEdT,buildMat,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable, &
	    intPts,ipWt,locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,disp)
	    implicit none
		
		real*8, intent(out) :: Ruk(33), dRdU(33,33), dRdT(33,10), stEn, dSEdT(10)
		integer, intent(in) :: buildMat, eType, numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		real*8, intent(in) :: intPts(3,8), ipWt(8), locNds(3,10), globNds(3,10)
		real*8, intent(in) :: orient(3,3), cMat(6,6), tExp(6), ABD(9,9), stExp(6), bStiff(6,6), btExp(6)
		real*8, intent(in) :: temp(10), disp(6,11)
		
		real*8 :: dedU(9,33), cdedU(9,33)
		real*8 :: instOri(3,240), drIdrG(3,16), elRot(3), ptRot(3)
		real*8 :: instDofMat(6,11)
		real*8 :: nnReal
		real*8 :: Nvec(11), Nx(11,3), detJ
		real*8 :: strain(6), stress(6)
		real*8 :: shDef(9), shFrcMom(9), beamDef(6), beamFrcMom(6), ux(3,3), rx(3,3)
		real*8 :: dfmdT(9,10), dsdT(6,10)
		real*8 :: duxdR(3,9)
		real*8 :: ptTemp, teProd, r3Fact
		integer :: i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, nd1, var1, nd2, var2
		
		if(eType .eq. 41 .or. eType .eq. 3 .or. eType .eq. 2) then
			call r_getInstOrient(instOri,orient,disp,numNds,0)
		endif
		
		Ruk(:) = r_0
		dRdU(:,:) = r_0
		dRdT(:,:) = r_0
		stEn = r_0
		dSEdT(:) = r_0
		do i1 = 1, numIntPts
		    call r_getIntPtData(Nvec,Nx,detJ,locNds,intPts(:,i1),eType,numNds)
			if(eType .eq. 41 .or. eType .eq. 3 .or. eType .eq. 2) then
			    call r_getInstDof(instDofMat,disp,globNds,instOri,orient,numNds,dofTable,0,0)
			    ptRot(:) = r_0
				do i2 = 1, numNds
				    do i3 = 1, 3
				        ptRot(i3) = ptRot(i3) + instDofMat(3+i3,i2)*Nvec(i2)
					enddo
				enddo
			    ux(:,:) = r_0
				rx(:,:) = r_0
			    do i2 = 1, 3
				    do i3 = 1, 3
					    do i4 = 1, 10
						    ux(i2,i3) = ux(i2,i3) + instDofMat(i2,i4)*Nx(i4,i3)
							rx(i2,i3) = rx(i2,i3) + instDofMat(i2+3,i4)*Nx(i4,i3)
						enddo
					enddo
				enddo
				if(eType .eq. 41 .or. eType .eq. 3) then
					call r_getShellDef(shDef,ux,ptRot,rx)
					ptTemp = r_0
					do i2 = 1, numNds
						ptTemp = ptTemp + Nvec(i2)*temp(i2)
					enddo
					shFrcMom(:) = r_0
					do i2 = 1, 9
						do i3 = 1, 9
							shFrcMom(i2) = shFrcMom(i2) + ABD(i2,i3)*shDef(i3)
						enddo
					enddo
					shFrcMom(1:6) = shFrcMom(1:6) - ptTemp*stExp(1:6)
					do i2 = 1, 6
					    stEn = stEn + r_p5*shFrcMom(i2)*shDef(i2)*detJ*ipWt(i1)
					    dfMdT(i2,1:numNds) = -Nvec(1:numNds)*stExp(i2)
						dSEdT(1:numNds) = dSEdT(1:numNds) + r_p5*dfMdT(i2,1:numNds)*shDef(i2)*detJ*ipWt(i1)
					enddo
				else
				    call r_getBeamDef(beamDef,ux,ptRot,rx)
					ptTemp = r_0
					do i2 = 1, numNds
						ptTemp = ptTemp + Nvec(i2)*temp(i2)
					enddo
					beamFrcMom(:) = r_0
					do i2 = 1, 6
						do i3 = 1, 6
							beamFrcMom(i2) = beamFrcMom(i2) + bStiff(i2,i3)*beamDef(i3)
						enddo
					enddo
					beamFrcMom(1:6) = beamFrcMom(1:6) - ptTemp*btExp(1:6)
					do i2 = 1, 6
					    stEn = stEn + r_p5*beamFrcMom(i2)*beamDef(i2)*detJ*ipWt(i1)
					    dfMdT(i2,1:numNds) = -Nvec(1:numNds)*btExp(i2)
						dSEdT(1:numNds) = dSEdT(1:numNds) + r_p5*dfMdT(i2,1:numNds)*beamDef(i2)*detJ*ipWt(i1)
					enddo
				endif
				i2 = numNds*dofPerNd + numIntDof
				do i4 = 1, i2
				    call r_getInstDof(instDofMat,disp,globNds,instOri,orient,numNds,dofTable,i4,0)
					ptRot(:) = r_0
					do i5 = 1, numNds
						ptRot(:) = ptRot(:) + instDofMat(4:6,i5)*Nvec(i5)
					enddo
					ux(:,:) = r_0
					rx(:,:) = r_0
					if(dofTable(1,i4) .lt. 4) then  !! differentiating wrt displacement dof
					    i7 = dofTable(2,i4)
						do i5 = 1, 3
							do i6 = 1, 3
								ux(i5,i6) = ux(i5,i6) + instDofMat(i5,i7)*Nx(i7,i6)
							enddo
						enddo
					else  !! wrt rotation dof
						do i6 = 1, 3
							do i7 = 1, 3
								do i8 = 1, numNds
									ux(i6,i7) = ux(i6,i7) + instDofMat(i6,i8)*Nx(i8,i7)
									rx(i6,i7) = rx(i6,i7) + instDofMat(i6+3,i8)*Nx(i8,i7)
								enddo
							enddo
						enddo
					endif
					if(eType .eq. 41 .or. eType .eq. 3) then
						call r_getShellDef(shDef,ux,ptRot,rx)
						if(buildMat .eq. 1) then
						    dedU(:,i4) = shDef(:)
							teProd = r_0
							do i5 = 1, 6
							    teProd = teProd - stExp(i5)*dedU(i5,i4) 
							enddo
							teProd = (detJ*ipWt(i1))*teProd
							dRdT(i4,1:numNds) = dRdT(i4,1:numNds) + teProd*Nvec(1:numNds)
						endif
						shDef = (detJ*ipWt(i1))*shDef
						do i5 = 1, 9
							Ruk(i4) = Ruk(i4) + shFrcMom(i5)*shDef(i5)
						enddo
					else
					    call r_getBeamDef(beamDef,ux,ptRot,rx)
						if(buildMat .eq. 1) then
						    dedU(1:6,i4) = beamDef(:)
							teProd = r_0
							do i5 = 1, 6
							    teProd = teProd - btExp(i5)*dedU(i5,i4) 
							enddo
							teProd = (detJ*ipWt(i1))*teProd
							dRdT(i4,1:numNds) = dRdT(i4,1:numNds) + teProd*Nvec(1:numNds)
						endif
						beamDef = (detJ*ipWt(i1))*beamDef
						do i5 = 1, 6
						    Ruk(i4) = Ruk(i4) + beamFrcMom(i5)*beamDef(i5)
						enddo
					endif
					if(buildMat .eq. 1 .and. nLGeom .ne. 0) then
					    do i10 = i4, i2
							call r_getInstDof(instDofMat,disp,globNds,instOri,orient,numNds,dofTable,i4,i10)
							ptRot(:) = r_0
							do i5 = 1, numNds
								ptRot(:) = ptRot(:) + instDofMat(4:6,i5)*Nvec(i5)
							enddo
							if(dofTable(1,i4) .lt. dofTable(1,i10)) then
							    var1 = dofTable(1,i4)
								nd1 = dofTable(2,i4)
								var2 = dofTable(1,i10)
								nd2 = dofTable(2,i10)
							else
							    var1 = dofTable(1,i10)
								nd1 = dofTable(2,i10)
								var2 = dofTable(1,i4)
								nd2 = dofTable(2,i4)
							endif
							ux(:,:) = r_0
							rx(:,:) = r_0
							if(var1 .lt. 4 .and. var2 .ge. 4) then  !! differentiating wrt displacement dof
								do i5 = 1, 3
									do i6 = 1, 3
										ux(i5,i6) = ux(i5,i6) + instDofMat(i5,nd1)*Nx(nd1,i6)
										! rx(i5,i6) = rx(i5,i6) + instDofMat(i5+3,nd2)*Nx(nd2,i6)
									enddo
								enddo
							elseif(var1 .ge. 4 .and. var2 .ge. 4) then !! wrt rotation dof
								do i6 = 1, 3
									do i7 = 1, 3
										do i8 = 1, numNds
											ux(i6,i7) = ux(i6,i7) + instDofMat(i6,i8)*Nx(i8,i7)
											rx(i6,i7) = rx(i6,i7) + instDofMat(i6+3,i8)*Nx(i8,i7)
										enddo
									enddo
								enddo
							endif
							if(eType .eq. 41 .or. eType .eq. 3) then
								call r_getShellDef(shDef,ux,ptRot,rx)
								shDef = (detJ*ipWt(i1))*shDef
								do i5 = 1, 9
									dRdU(i4,i10) = dRdU(i4,i10) + shFrcMom(i5)*shDef(i5)
								enddo
								dRdU(i10,i4) = dRdU(i4,i10)
							else
								call r_getBeamDef(beamDef,ux,ptRot,rx)
								beamDef = (detJ*ipWt(i1))*beamDef
								do i5 = 1, 6
									dRdU(i4,i10) = dRdU(i4,i10) + beamFrcMom(i5)*beamDef(i5)
								enddo
								dRdU(i10,i4) = dRdU(i4,i10)
							endif					    
						enddo
					endif
				enddo
			else
			    ux(:,:) = r_0
				do i2 = 1, 3
				    do i3 = 1, 3
					    do i4 = 1, 11
						    ux(i2,i3) = ux(i2,i3) + disp(i2,i4)*Nx(i4,i3)
						enddo
					enddo
				enddo
			    call r_getNLStrain(strain,ux,Nx,orient,0,0,dofTable)
				ptTemp = r_0
				do i2 = 1, numNds
				    ptTemp = ptTemp + Nvec(i2)*temp(i2)
				enddo
				stress(:) = r_0
				dsdT(:,:) = r_0
				do i2 = 1, 6
				    do i3 = 1, 6
					    stress(i2) = stress(i2) + cMat(i2,i3)*(strain(i3) - ptTemp*tExp(i3))
						dsdT(i2,1:numNds) = dsdT(i2,1:numNds) - cMat(i2,i3)*Nvec(1:numNds)*tExp(i3)
					enddo
					stEn = stEn + r_p5*stress(i2)*strain(i2)*detJ*ipWt(i1)
					dSEdT(1:numNds) = dSEdT(1:numNds) + r_p5*dsdT(i2,1:numNds)*strain(i2)*detJ*ipWt(i1)
				enddo
				i2 = numNds*dofPerNd + numIntDof
				do i3 = 1, i2
				    call r_getNLStrain(strain,ux,Nx,orient,i3,0,dofTable)
					if(buildMat .eq. 1) then
					    dedU(1:6,i3) = strain(:)
						teProd = r_0
						do i5 = 1, 6
						    do i6 = 1, 6
							    teProd = teProd - tExp(i5)*cMat(i5,i6)*dedU(i6,i3)
                            enddo							
						enddo
						teProd = (detJ*ipWt(i1))*teProd
						dRdT(i4,1:numNds) = dRdT(i4,1:numNds) + teProd*Nvec(1:numNds)
					endif
					strain = (detJ*ipWt(i1))*strain
					do i4 = 1, 6
					    Ruk(i3) = Ruk(i3) + stress(i4)*strain(i4)
					enddo
					if(buildMat .eq. 1 .and. nLGeom .ne. 0) then
						do i5 = i3, i2
							call r_getNLStrain(strain,ux,Nx,orient,i3,i5,dofTable)
							strain = (detJ*ipWt(i1))*strain
							do i4 = 1, 6
								dRdU(i3,i5) = dRdU(i3,i5) + stress(i4)*strain(i4)
							enddo
							dRdU(i5,i3) = dRdU(i3,i5)
						enddo
					endif
				enddo
			endif
			if(buildMat .eq. 1) then
				cdedU(:,:) = r_0
				if(eType .eq. 41 .or. eType .eq. 3) then
					i5 = numNds*dofPerNd + numIntDof
					do i2 = 1, 9
						do i3 = 1, i5
							do i4 = 1, 9
								cdedU(i2,i3) = cdedU(i2,i3) + ABD(i2,i4)*dedU(i4,i3)
							enddo
						enddo
					enddo
					cdedU(:,:) = (detJ*ipWt(i1))*cdedU(:,:)
					do i2 = 1, i5
						do i3 = 1, i5
							do i4 = 1, 9
								dRdU(i2,i3) = dRdU(i2,i3) + dedU(i4,i2)*cdedU(i4,i3)
							enddo
						enddo
					enddo
				elseif(eType .eq. 2) then
					i5 = numNds*dofPerNd + numIntDof
					do i2 = 1, 6
						do i3 = 1, i5
							do i4 = 1, 6
								cdedU(i2,i3) = cdedU(i2,i3) + bStiff(i2,i4)*dedU(i4,i3)
							enddo
						enddo
					enddo
					cdedU(:,:) = (detJ*ipWt(i1))*cdedU(:,:)
					do i2 = 1, i5
						do i3 = 1, i5
							do i4 = 1, 6
								dRdU(i2,i3) = dRdU(i2,i3) + dedU(i4,i2)*cdedU(i4,i3)
							enddo
						enddo
					enddo
				else
				    i5 = numNds*dofPerNd + numIntDof
					do i2 = 1, 6
						do i3 = 1, i5
							do i4 = 1, 6
								cdedU(i2,i3) = cdedU(i2,i3) + cMat(i2,i4)*dedU(i4,i3)
							enddo
						enddo
					enddo
					cdedU(:,:) = (detJ*ipWt(i1))*cdedU(:,:)
					do i2 = 1, i5
						do i3 = 1, i5
							do i4 = 1, 6
								dRdU(i2,i3) = dRdU(i2,i3) + dedU(i4,i2)*cdedU(i4,i3)
							enddo
						enddo
					enddo
				endif
			endif
		enddo
		
	end subroutine r_getElRuk
	
	subroutine r_getElRum(Rum,dRdA,buildMat,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
	    locNds, globNds, orient, statInOri, den, sMass, bMass, disp, acc)
	    implicit none
		
		real*8, intent(out) :: Rum(33), dRdA(33,33)
		integer, intent(in) :: buildMat, eType, numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		real*8, intent(in) :: intPts(3,8), ipWt(8)
		real*8, intent(in) :: locNds(3,10), globNds(3,10), orient(3,3)
		real*8, intent(in) :: statInOri(3,240)
        real*8, intent(in) :: den, sMass(6,6), bMass(6,6)
        real*8, intent(in) :: disp(6,11), acc(6,11)
		
		real*8 :: instDofMat(6,11), instNodalAcc(6,4), iNAM(6,4), inerFrc(6), delU(6), accI
		real*8 :: dudU(6,33), MdudU(6,33), thisMass(6,6)
		real*8 :: Nvec(11), Nx(11,3), detJ
		
		integer :: i1, i2, i3, i4, i5, i6
		
		if(eType .eq. 41 .or. eType .eq. 3 .or. eType .eq. 2) then
		    instDofMat(:,:) = r_0
			call r_getInstDof(instDofMat,disp,globNds,statInOri,orient,numNds,dofTable,0,0)
			i1 = numNds*dofPerNd
			instNodalAcc(:,:) = r_0
			do i2 = 1, i1
			    call r_getInstDof(instDofMat,disp,globNds,statInOri,orient,numNds,dofTable,i2,0)
				i4 = dofTable(1,i2)
				i5 = dofTable(2,i2)
				instNodalAcc(:,:) = instNodalAcc(:,:) + acc(i4,i5)*instDofMat(:,1:4)
			enddo
			if(eType .eq. 41 .or. eType .eq. 3) then
			    thisMass = sMass
			else
			    thisMass = bMass
			endif
			iNAM(:,:) = r_0
			do i1 = 1, 6
				do i2 = 1, numNds
					do i3 = 1, 6
						iNAM(i1,i2) = iNAM(i1,i2) + thisMass(i1,i3)*instNodalAcc(i3,i2)
					enddo
				enddo
			enddo
		endif
		
		Rum(:) = r_0
		dRdA(:,:) = r_0
		do i1 = 1, numIntPts
		    call r_getIntPtData(Nvec,Nx,detJ,locNds,intPts(:,i1),eType,numNds)
			if(eType .eq. 41 .or. eType .eq. 3 .or. eType .eq. 2) then
				inerFrc(:) = r_0
				do i2 = 1, 6
					do i3 = 1, numNds
						inerFrc(i2) = inerFrc(i2) + iNAM(i2,i3)*Nvec(i3)
					enddo
				enddo
				i3 = numNds*dofPerNd
				do i4 = 1, i3
				    delU(:) = r_0
					call r_getInstDof(instDofMat,disp,globNds,statInOri,orient,numNds,dofTable,i4,0)
					do i5 = 1, 6
					    do i6 = 1, numNds
						    delU(i5) = delU(i5) + instDofMat(i5,i6)*Nvec(i6)
						enddo
					enddo
					if(buildMat .eq. 1) then
					    dudU(:,i4) = delU(:)
					endif
					delU = (detJ*ipWt(i1))*delU
					do i5 = 1, 6
					    Rum(i4) = Rum(i4) + inerFrc(i5)*delU(i5)
					enddo
					i2 = i2 + 11
				enddo
				if(buildMat .eq. 1) then
				    MdudU(:,:) = r_0
					do i4 = 1, 6
					    do i5 = 1, i3
						    do i6 = 1, 6
							    MdudU(i4,i5) = MdudU(i4,i5) + thisMass(i4,i6)*dudU(i6,i5)
							enddo
						enddo
					enddo
					MdudU = (detJ*ipWt(i1))*MdudU
					do i4 = 1, i3
					    do i5 = 1, i3
						    do i6 = 1, 6
							    dRdA(i4,i5) = dRdA(i4,i5) + dudU(i6,i4)*MdudU(i6,i5)
							enddo
						enddo
					enddo
				endif
			else
			    dudU(:,:) = r_0
			    i2 = numNds*dofPerNd
				do i3 = 1, i2
				    i4 = dofTable(1,i3)
					i5 = dofTable(2,i3)
					accI = r_0
					do i6 = 1, 11
					    accI = accI + acc(i4,i6)*Nvec(i6)
					enddo
					accI = (detJ*ipWt(i1))*accI
					Rum(i3) = Rum(i3) + den*accI*Nvec(i5)
					if(buildMat .eq. 1) then
					    dudU(i4,i3) = Nvec(i5)
					endif
				enddo
				if(buildMat .eq. 1) then
					MdudU = (den*detJ*ipWt(i1))*dudU
					do i4 = 1, i2
					    do i5 = 1, i2
						    do i6 = 1, 3
							    dRdA(i4,i5) = dRdA(i4,i5) + dudU(i6,i4)*MdudU(i6,i5)
							enddo
						enddo
					enddo
				endif
            endif
		enddo
		
	end subroutine r_getElRum
	
	subroutine r_getElRu(Ru,dRdU,buildMat,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
	    locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,disp,vel,acc,pDisp,pVel,pAcc, &
		statInOri,den,sMass,bMass)
	    implicit none
		
		real*8, intent(out) :: Ru(33), dRdU(33,33)
		integer, intent(in) :: buildMat, eType, numNds, dofPerNd, numIntDof, numIntPts, dofTable(2,33)
		real*8, intent(in) :: intPts(3,8), ipWt(8), locNds(3,10), globNds(3,10), orient(3,3)
		real*8, intent(in) :: cMat(6,6), tExp(6), ABD(9,9), stExp(6), bStiff(6,6), btExp(6)
        real*8, intent(in) :: temp(10), disp(6,11), vel(6,11), acc(6,11), pDisp(6,11), pVel(6,11), pAcc(6,11)
		real*8, intent(in) :: statInOri(3,240), den, sMass(6,6), bMass(6,6)
		
		real*8 :: Ruk(33), Rum(33), Ruc(33), dRdVar(33,33), dRdT(33,10), stEn, dSEdT(10)
		real*8 :: velNext(6,11), accNext(6,11)
		real*8 :: c1, c2
		integer :: gt1, gt2
		
		call r_getElRuk(Ruk,dRdVar,dRdT,stEn,dSEdT,buildMat,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
	    locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,disp)
		
		Ru(:) = Ruk(:)
		if(buildMat .eq. 1) then
		    dRdU(:,:) = dRdVar(:,:)
		endif
		
		if(dynamic .eq. 1) then
		    c1 = r_p5/nMBeta
			c2 = r_2/(delT*delT)
		    accNext(:,:) = c1*(c2*(disp(:,:)-pDisp(:,:)-delT*pVel(:,:)) - (r_1 - r_2*nMBeta)*pAcc(:,:))
			velNext(:,:) = pVel(:,:) + delT*(r_1*(1d0-nMGamma)*pAcc(:,:) + nMGamma*accNext(:,:))
			
			call r_getElRum(Rum,dRdVar,buildMat,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
			locNds,globNds,orient,statInOri,den,sMass,bMass,disp,accNext)
			
			Ru(:) = Ru(:) + Rum(:)
			if(buildMat .eq. 1) then
			    dRdU(:,:) = dRdU(:,:) + (c1*c2)*dRdVar(:,:)
			endif
			
			call r_greater(gt1,r_1*rayCoefK,r_0)
		    call r_greater(gt2,r_1*rayCoefM,r_0)
			
			if(gt1 .eq. 1) then
				call r_getElRuk(Ruc,dRdVar,dRdT,stEn,dSEdT,buildMat,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
				locNds,globNds,orient,cMat,tExp,ABD,stExp,bStiff,btExp,temp,rayCoefK*velNext)
                Ru(:) = Ru(:) + Ruc
				if(buildMat .eq. 1) then
			        dRdU(:,:) = dRdU(:,:) + rayCoefK*(delT*nMGamma)*(c1*c2)*dRdVar(:,:)
			    endif
			endif
			
			if(gt2 .eq. 1) then
				call r_getElRum(Ruc,dRdVar,buildMat,eType,numNds,dofPerNd,numIntDof,numIntPts,dofTable,intPts,ipWt, &
			    locNds,globNds,orient,statInOri,den,sMass,bMass,disp,rayCoefM*velNext)
                Ru(:) = Ru(:) + Ruc
				if(buildMat .eq. 1) then
				    dRdU(:,:) = dRdU(:,:) + rayCoefM*(delT*nMGamma)*(c1*c2)*dRdVar(:,:)
				endif
			endif
		endif
	
	end subroutine r_getElRu
	
	subroutine r_condenseMat(mat,mDim,resVec,bVec,zVec,stRow,endRow)
	    implicit none
		
		integer, intent(in) :: mDim, stRow, endRow
		real*8, intent(out) :: mat(mDim,mDim), resVec(mDim), bVec(mDim), zVec(mDim)
		
		integer :: i1, i2, i3
		
		call r_LUFactor(mat,mDim,stRow,endRow)
		
		do i1 = 1, stRow-1
		    bVec(:) = mat(:,i1)
			resVec(:) = r_0
		    call r_SolveLUxb(resVec, mat, bVec, zVec, mDim, stRow, endRow)
			do i2 = 1, stRow-1
			    do i3 = stRow, endRow
				    mat(i2,i1) = mat(i2,i1) - mat(i2,i3)*resVec(i3)
				enddo
			enddo
		enddo
		
	end subroutine r_condenseMat

end module ASenD_r_elementEqns