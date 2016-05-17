SUBROUTINE TIME_Update_U
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L
    INTEGER:: NI1, NJ1, NK1
    INTEGER:: MLtemp

    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1




	  if(Kind_model == 4)then
	     MLtemp=6 !ML-1 by xu
	  else
	     MLtemp=ML
	  endif
!
	  if(     Kind_dual == 1 )then
           DO L=1,MLtemp
           DO J=1,NJ1
           DO K=1,NK1
           DO I=1,NI1
	        Wt0(L,I,J,K)=  V(L,I,J,K)
	     ENDDO
	     ENDDO
	     ENDDO
	     ENDDO

	  elseif( Kind_dual == 2)then
           DO L=1,MLtemp
           DO J=1,NJ1
           DO K=1,NK1
           DO I=1,NI1
	           Wt1(L,I,J,K)=Wt0(L,I,J,K)
	           Wt0(L,I,J,K)=  V(L,I,J,K)
	     ENDDO
	     ENDDO
	     ENDDO
	     ENDDO
	  endif

END SUBROUTINE TIME_Update_U