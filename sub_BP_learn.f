      subroutine BPlearn(input,output,N,w,L,nsample)

      integer L,nsample
      real*8 input(nsample,101),output(nsample,101)
      real*8 w(L,101,101)
      integer N(L)

      real*8 alpha
      parameter (alpha=0.005)

      real*8 tol
      parameter (tol=0.0001)

      real*8 dw(L,101,101),ddw(L,101,101,2)
      real*8 y(L,101),T(101)
      real*8 d(L,101),temp,temp2,temp3
      integer i,j,k,nloop,rn
      real*8 error(500000)

      real rand2
      double precision r2
      r2=5.0


ccccccccccccccccccccccccccccccccccc

      do nloop=1,500000

         rn=int(rand2(r2)*nsample)+1
         do j=1,N(1)-1
            y(1,j)=input(rn,j)
         enddo
         do j=1,N(L)
            T(j)=output(rn,j)
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



         do j=1,N(L)
            d(L,j)=(T(j)-y(L,j))*(1-y(L,j)**2)
         enddo

         do i=L-1,2,-1
            do j=1,N(i)-1
               temp3=0
               do k=1,N(i+1)
                  temp3=temp3+w(i+1,k,j)*d(i+1,k)
               enddo
               d(i,j)=(1-y(i,j)**2)*temp3
             
            enddo
         enddo

         do i=2,L-1
            do j=1,N(i)-1
               do k=1,N(i-1)-1
                  dw(i,j,k)=0
                  dw(i,j,k)=dw(i,j,k)+alpha*d(i,j)*y(i-1,k)
               enddo
               dw(i,j,N(i-1))=alpha*d(i,j)*1
            enddo
         enddo


         do i=2,L-1
            do j=1,N(i)-1
               do k=1,N(i-1)-1
                  ddw(i,j,k,2)=dw(i,j,k)
               enddo
               ddw(i,j,N(i-1),2)=dw(i,j,N(i-1))
            enddo
         enddo


         do i=L,L
            do j=1,N(i)
               do k=1,N(i-1)-1
                  dw(i,j,k)=0
                  dw(i,j,k)=dw(i,j,k)+alpha*d(i,j)*y(i-1,k)
               enddo
               dw(i,j,N(i-1))=alpha*d(i,j)*1
            enddo
         enddo

         do i=L,L
            do j=1,N(i)
               do k=1,N(i-1)-1
                  ddw(i,j,k,2)=dw(i,j,k)
               enddo
               ddw(i,j,N(i-1),2)=dw(i,j,N(i-1))
            enddo
         enddo

    
         do i=2,L-1
            do j=1,N(i)-1
               do k=1,N(i-1)
                  w(i,j,k)=w(i,j,k)+dw(i,j,k)+0.0*ddw(i,j,k,1)
               enddo
            enddo
         enddo     

         do i=L,L
            do j=1,N(i)
               do k=1,N(i-1)
                  w(i,j,k)=w(i,j,k)+dw(i,j,k)+0.75*ddw(i,j,k,1)
               enddo
            enddo
         enddo     

         do i=2,L
            do j=1,N(i)
               do k=1,N(i-1)
                  ddw(i,j,k,1)=ddw(i,j,k,2)
               enddo
            enddo
         enddo
         
         error(nloop)=0
         do j=1,N(L)
               error(nloop)=error(nloop)+abs(T(j)-y(L,j))
         enddo

         if(sqrt(error(nloop)).le.tol)then
            goto 200
         endif

      enddo

 200  continue

ccccccccccccccccccccc

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc
c>> random number generator


        real  function rand2(r2)
        double precision s,u,v,r2
        s=65536.0
        u=2053.0
        v=13849.0
        m=r2/s
        r2=r2-m*s
        r2=u*r2+v
        m=r2/s
        r2=r2-m*s
        rand2=r2/s
        return
        end
        
ccccccccccccccccccccccccccccccccccccccccccccccccc
