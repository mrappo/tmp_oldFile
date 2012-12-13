      Subroutine jetbin(p,maxparts,passed)
c--- given the event momenta in p (with maxparts entries),
c--- performs some generic WBF cuts (a la Del Duca et al.)
c---  a point that fails the cuts returns passed=.false.
c--- Note: implements Eq. (3.2) of CEZ paper
      implicit none
      include 'constants.f'
      logical passed, first 
      integer j,maxparts,found,j1,j2,j3,NumJets
      double precision pt,etarap,p(mxpart,4),ptj,ptj1,ptj2,
     . etaj1,etaj2,mj1j2
      double precision ptj3,etaj3
      character*2 plabel(mxpart)
      common/plabel/plabel
      data first/.true./
      double precision JetPtMin, JetEtaCut, JetEtaMax,JetDEtaMin, 
     . JetMjjMin, JetMjjMax
      logical jetveto

************************* Selected CUTS *******************************
     
      JetPtMin   = 30d0
      JetEtaMax  = 4.7d0
      NumJets    = 2
      passed     = .false.

***************************** START CUTS *******************************
      found=0
      ptj1=0d0
      ptj2=0d0     
      

c--- check for 2 tagging jets (highest pt) 
      do j=3,maxparts
       if (  (plabel(j) .eq. 'pp') .or. (plabel(j) .eq. 'qj')
     .  .or. (plabel(j) .eq. 'bq') .or. (plabel(j) .eq. 'ba')) then
          ptj=pt(j,p)
c          write(*,*) 'j', j, ptj 
          if (ptj .gt. min(ptj1,ptj2)) then
            found=found+1
            if (found .eq. 1) then
              ptj1=ptj
              j1=j
            endif
            if (found .ge. 2) then
              if (ptj .gt. ptj1) then
                 ptj2=ptj1
                 ptj1=ptj
                 j2=j1
                 j1=j
              else
                 ptj2=ptj
                 j2=j
              endif
            endif
          endif
        endif
      enddo

c      write(*,*) 'found ',found 

      if (found .lt. NumJets ) goto 999
      if (ptj1  .lt. JetPtMin) goto 999
      if (ptj2  .lt. JetPtMin) goto 999

c      write(*,*) 'ptj1 ',ptj1, 'ptj2 ',ptj2 

      etaj1= 0d0
      etaj2= 0d0
      etaj1=etarap(j1,p)
      etaj2=etarap(j2,p)

c     eta acceptance selection                                                                                                                                    

      if(abs(etaj1) .gt. JetEtaMax) goto 999
      if(abs(etaj2) .gt. JetEtaMax) goto 999



********************** END OF CUT CROSS SECTIONS ***********************

      passed=.true.

c--- this is the alternate return if the cuts are failed    
  999 continue
      
      return
      
      end
      
     
