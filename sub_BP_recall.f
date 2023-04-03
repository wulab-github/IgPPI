      subroutine BPrecall(input2,output2,testout,N,w,L,nsample)

      integer L,nsample
      real*8 input2(nsample,101),output2(nsample,101)
      real*8 w(L,101,101)
      integer N(L)
      real*8 testout(nsample,101)

      real*8 y(L,101),T(101)
      real*8 temp
      integer rn,i,j,k

cccccccccccccccccccccccccccc

      do rn=1,nsample

         do j=1,N(1)-1
            y(1,j)=input2(rn,j)
         enddo
         do j=1,N(L)
            T(j)=output2(rn,j)
         enddo

         do i=2,L-1
            do j=1,N(i)-1
               temp=0
               do k=1,N(i-1)-1
                  temp=temp+w(i,j,k)*y(i-1,k)
               enddo
               temp=temp+w(i,j,N(i-1))
               y(i,j)=(exp(temp)-exp(-temp))/(exp(temp)+exp(-temp))
            enddo
         enddo

         do i=L,L
            do j=1,N(i)
               temp=0
               do k=1,N(i-1)-1
                  temp=temp+w(i,j,k)*y(i-1,k)
               enddo
               temp=temp+w(i,j,N(i-1))
               y(i,j)=(exp(temp)-exp(-temp))/(exp(temp)+exp(-temp))
            enddo
         enddo

         do i=1,N(L)
            testout(rn,i)=y(L,i)
         enddo

      enddo

ccccccccccccccccccccccc
      
      return
      end
