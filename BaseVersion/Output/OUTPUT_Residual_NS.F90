SUBROUTINE OUTPUT_Residual_NS
    USE Global
    IMPLICIT NONE
    integer::i,j,k,NJ1,NK1
    real::msum,rv,ru,rw,wei1,wei2,wei

    !INTEGER,PARAMETER:: myid=0 !Not used until MPI


    OPEN(UNIT=16,FILE=ThisBlock%FLHIS, MODE='WRITE', STATUS='UNKNOWN',ACCESS='APPEND') 

	   write(*,*) 
	   write(*,*) " Residuals of NS equ in blk: ",ThisBlock%id_present_blk

        WRITE(16,40) myid,KI,Itr_Chld,ThisBlock%DRSa,ThisBlock%IMa,ThisBlock%JMa,ThisBlock%KMa,ThisBlock%DUa
        WRITE(*,40) myid,KI,Itr_Chld,ThisBlock%DRSa,ThisBlock%IMa,ThisBlock%JMa,ThisBlock%KMa,ThisBlock%DUa
	   
!        WRITE(16,40) myid,KI,Itr_Chld,ThisBlock%DUa1,ThisBlock%IMa1,ThisBlock%JMa1,ThisBlock%KMa1,ThisBlock%DRSa1
!        WRITE(*,40) myid,KI,Itr_Chld,ThisBlock%DUa1,ThisBlock%IMa1,ThisBlock%JMa1,ThisBlock%KMa1,ThisBlock%DRSa1
!        
!	  WRITE(16,40) myid,KI,Itr_Chld,ThisBlock%DUa2,ThisBlock%IMa2,ThisBlock%JMa2,ThisBlock%KMa2,ThisBlock%DRSa2
!        WRITE(*,40) myid,KI,Itr_Chld,ThisBlock%DUa2,ThisBlock%IMa2,ThisBlock%JMa2,ThisBlock%KMa2,ThisBlock%DRSa2
!        
!	  WRITE(16,40) myid,KI,Itr_Chld,ThisBlock%DUa3,ThisBlock%IMa3,ThisBlock%JMa3,ThisBlock%KMa3,ThisBlock%DRSa3
!        WRITE(*,40) myid,KI,Itr_Chld,ThisBlock%DUa3,ThisBlock%IMa3,ThisBlock%JMa3,ThisBlock%KMa3,ThisBlock%DRSa3
!        
!	  WRITE(16,40) myid,KI,Itr_Chld,ThisBlock%DUa4,ThisBlock%IMa4,ThisBlock%JMa4,ThisBlock%KMa4,ThisBlock%DRSa4
!        WRITE(*,40) myid,KI,Itr_Chld,ThisBlock%DUa4,ThisBlock%IMa4,ThisBlock%JMa4,ThisBlock%KMa4,ThisBlock%DRSa4
!        
!	  WRITE(16,40) myid,KI,Itr_Chld,ThisBlock%DUa5,ThisBlock%IMa5,ThisBlock%JMa5,ThisBlock%KMa5,ThisBlock%DRSa5
!        WRITE(*,40) myid,KI,Itr_Chld,ThisBlock%DUa5,ThisBlock%IMa5,ThisBlock%JMa5,ThisBlock%KMa5,ThisBlock%DRSa5


40      FORMAT(2X,I5,2X,I8,2X,I5,2x,E12.6,2X,3I5,2X,E12.6)
        close(16)
!
    if(ThisBlock%ID_Present_blk==1)then
        NJ1=ThisBlock%NJ1
        NK1=ThisBlock%NK1
            i=4
        do k=1,NK1
        do j=1,NJ1
            wei1=Vol(i-1,j,k)
            wei2=Vol(i,j,k)
            wei=wei1+wei2
            ru=(wei1*V(2,i,j,k)+wei2*V(2,i-1,j,k))/wei
            rv=(wei1*V(3,i,j,k)+wei2*V(3,i-1,j,k))/wei
            rw=(wei1*V(4,i,j,k)+wei2*V(4,i-1,j,k))/wei
            msum=msum+ru*SD(1,1,i,j,k)+rv*SD(1,2,i,j,k)+rw*SD(1,3,i,j,k)
        enddo
        enddo
        msum=msum*Roref*Vref
        OPEN(UNIT=16,FILE="resu/MassIn.dat",MODE='WRITE', STATUS='UNKNOWN',ACCESS='APPEND') 
            write(16,*)KI,msum
        CLOSE(16)
    else if(ThisBlock%ID_Present_blk==2)then
        NJ1=ThisBlock%NJ1
        NK1=ThisBlock%NK1
            i=ThisBlock%NI-4    !5
        do k=1,NK1
        do j=1,NJ1
            wei1=Vol(i-1,j,k)
            wei2=Vol(i,j,k)
            wei=wei1+wei2
            ru=(wei1*V(2,i,j,k)+wei2*V(2,i-1,j,k))/wei
            rv=(wei1*V(3,i,j,k)+wei2*V(3,i-1,j,k))/wei
            rw=(wei1*V(4,i,j,k)+wei2*V(4,i-1,j,k))/wei
            msum=msum+ru*SD(1,1,i,j,k)+rv*SD(1,2,i,j,k)+rw*SD(1,3,i,j,k)
        enddo
        enddo
        msum=msum*Roref*Vref
        OPEN(UNIT=16,FILE="resu/MassOut.dat",MODE='WRITE', STATUS='UNKNOWN',ACCESS='APPEND') 
            write(16,*)KI,msum
        CLOSE(16)
    endif

!    do iblock=1,Max_Block
!        call Global_SetPointers(iblock)
!    enddo

END SUBROUTINE OUTPUT_Residual_NS
