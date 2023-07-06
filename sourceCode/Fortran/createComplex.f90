program createComplex
    implicit none
	
	character(len=128) :: sourceFile, destFile
	character(len=256) :: readLine, writeLine, trailingStr
    integer :: iosVal, found, i1

    call get_command_argument(1,sourceFile)
	
	destFile = sourceFile
	i1 = index(sourceFile,'r_')
	destFile(i1:i1+1) = 'c_'
	
	open(unit=9, file=sourceFile, blank='NULL', action='read')
	open(unit=15, file=destFile, action='write', status='replace')
	
	read(9,'(A)',iostat=iosVal) readLine
	do while(iosVal .eq. 0)
	    found = 1
		do while(found .eq. 1)
		    found = 0
			i1 = index(readLine,'real*8')
			if(i1 .gt. 0) then
			    trailingStr = readLine(i1+6:250)
			    writeLine = readLine(1:i1-1) // 'complex*16' // trim(trailingStr)
				readLine = writeLine
				found = 1
			endif
			i1 = index(readLine,'r_')
			if(i1 .gt. 0) then
			    readLine(i1:i1+1) = 'c_'
				found = 1
			endif
		enddo
		write(15,*) trim(readLine)
		read(9,'(A)',iostat=iosVal) readLine
	enddo
	
	close(9)
	close(15)
	
end program