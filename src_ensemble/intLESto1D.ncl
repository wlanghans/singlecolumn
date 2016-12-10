;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get the mean profiles from BOMEX
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
bomexdir="/global/home/users/langhans/BOMEX_LES/"  
fin=addfile(bomexdir+"/BOMEXmean.nc","r")
zmean=doubletofloat(fin->height)
nzmean=dimsizes(zmean)
thl1=doubletofloat(fin->THL(0,:))
qc1=doubletofloat(fin->QC(0,:))
qt1=doubletofloat(fin->QT(0,:))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get the flux profiles from BOMEX
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fin=addfile(bomexdir+"/BOMEXflx.nc","r")
zflx=doubletofloat(fin->height)
nzflx=dimsizes(zflx)
thlflx=doubletofloat(fin->THL(0,:))
qtflx=doubletofloat(fin->QT(0,:))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; get the mean profiles from CPBL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
lesdir=(/"/global/scratch/langhans/witek_les/run12/","/global/scratch/langhans/witek_les/run14/"/)
ncpbl=dimsizes(lesdir)
;loop over cases
nzcpbl1  = floattointeger((3010.-10.)/20.) + 1
zcpbl1=fspan(10.,3010.,nzcpbl1)

do n=0,ncpbl-1

t_ind=(/3.,5./)*3600.
; find right LES file
lesfiles = systemfunc("ls " + lesdir(n) + "dcblm_*.nc")
num_lesfiles=dimsizes(lesfiles)

ind_les=new(dimsizes(t_ind),integer)

do t=0,dimsizes(t_ind)-1
mindt=10000000000.
print("Finding LES file # "+ t)
do i=0,num_lesfiles-1
  f=addfile(lesfiles(i),"r")
  tles=f->time
  dt=abs(tles-t_ind(t))
  if (dt.lt.mindt) then
    ind_les(t) = i
    mindt= dt
  end if
end do
f=addfile(lesfiles(ind_les(t)),"r")
tles=f->time
print("Using Witek CPBL LES results at time " + tles + " sec")
print("Absolute time difference to singlecolumn is:" + mindt + " sec")

if (n.eq.0.and.t.eq.0) then
  thlcpbl = new((/ncpbl*dimsizes(t_ind),nzcpbl1/),float)
  qtcpbl  = new((/ncpbl*dimsizes(t_ind),nzcpbl1/),float)
end if

; get average thetal and qt
thlcpbl(dimsizes(t_ind)*n+t,:) = dim_avg_n(f->t(:,:,0:nzcpbl1-1),(/0,1/))
qtcpbl (dimsizes(t_ind)*n+t,:) = dim_avg_n(f->q(:,:,0:nzcpbl1-1),(/0,1/))


end do ; loop over time
end do ; loop over CPBL cases

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;interpolate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;create target grid
nzmean2  = floattointeger((3010.-30.)/20.) + 1
zmean2=fspan(30.,3010.,nzmean2)
nzflx2  = floattointeger((3000.-20.)/20.) + 1
zflx2=fspan(20.,3000.,nzflx2)

; create target variables
; 3 for BOMEX and 4 each (theta and qt at 3 and 5 hours) of CPBL2 and CPBL4
mean=new((/7,nzmean2/),float)


do k=0,nzmean2-1

  lfound=0
  do l=0,nzmean-2
    if (zmean(l+1).ge.zmean2(k).and.zmean(l).le.zmean2(k)) then
      mean(0,k) = thl1(l) + (zmean2(k)-zmean(l))*(thl1(l+1)-thl1(l))/(zmean(l+1)-zmean(l))
      mean(1,k) = qt1(l) + (zmean2(k)-zmean(l))*(qt1(l+1)-qt1(l))/(zmean(l+1)-zmean(l))
      mean(2,k) = qc1(l) + (zmean2(k)-zmean(l))*(qc1(l+1)-qc1(l))/(zmean(l+1)-zmean(l))
      lfound=1
    end if
  end do
  if (lfound.eq.0) then
      mean(0,k) = -999.
      mean(1,k) = -999.
      mean(2,k) = -999.
  end if

  lfound=0
  do l=0,nzcpbl1-2
    if (zcpbl1(l+1).ge.zmean2(k).and.zcpbl1(l).le.zmean2(k)) then
      mean(3,k) = thlcpbl(0,l) + (zmean2(k)-zcpbl1(l))*(thlcpbl(0,l+1)-thlcpbl(0,l))/(zcpbl1(l+1)-zcpbl1(l))
      mean(4,k) = qtcpbl(1,l)  + (zmean2(k)-zcpbl1(l))*(qtcpbl(1,l+1)-qtcpbl(1,l))  /(zcpbl1(l+1)-zcpbl1(l))
      mean(5,k) = thlcpbl(2,l) + (zmean2(k)-zcpbl1(l))*(thlcpbl(2,l+1)-thlcpbl(2,l))/(zcpbl1(l+1)-zcpbl1(l))
      mean(6,k) = qtcpbl(3,l)  + (zmean2(k)-zcpbl1(l))*(qtcpbl(3,l+1)-qtcpbl(3,l))  /(zcpbl1(l+1)-zcpbl1(l))
      lfound=1
    end if
  end do
  if (lfound.eq.0) then
      mean(3,k) = -999.
      mean(4,k) = -999.
      mean(5,k) = -999.
      mean(6,k) = -999.
  end if
  
  
end do


flx=new((/2,nzflx2/),float)

do k=0,nzflx2-1

  lfound=0
  do l=0,nzflx-2
    if (zflx(l+1).ge.zflx2(k).and.zflx(l).le.zflx2(k)) then
      flx(0,k) = thlflx(l) + (zflx2(k)-zflx(l))*(thlflx(l+1)-thlflx(l))/(zflx(l+1)-zflx(l))
      flx(1,k) = qtflx(l) + (zflx2(k)-zflx(l))*(qtflx(l+1)-qtflx(l))/(zflx(l+1)-zflx(l))
      lfound=1
    end if
  end do
  if (lfound.eq.0) then
      flx(0,k)=-999.
      flx(1,k)=-999.
  end if

end do

       system("/bin/rm -f benchmark.nc")   ; remove any pre-existing file
       ncdf = addfile("benchmark.nc" ,"c")  ; open output netCDF file

       fAtt               = True    
       fAtt@title         = "Benchmark profiles from LES"
       fAtt@source_file   =  "Witek 2011 JAS (CPBL); Suselj 2012 JAS (BOMEX)"
       fAtt@Conventions   = "None"
       fAtt@creation_date = systemfunc ("date")
       fileattdef( ncdf, fAtt )            ; copy file attributes

       outmean=(/mean((/0,3,5/),:)/)

       outmean!0 = "case"
       outmean!1 = "height1"
       outmean&case = (/0,1,2/)
       outmean&height1  = zmean2
       outmean@long_name = "mean profiles"           ; assign attributes 
       outmean@_FillValue = -999.           ; assign attributes 
       outmean2=(/mean(2,:)/)
       outmean2!0 = "height1"
       outmean2&height1  = zmean2
       outmean2@long_name = "qc mean BOMEX profile"           ; assign attributes 
       outmean2@_FillValue = -999.           ; assign attributes 

       ncdf->THL  = outmean                                         
       outmean=(/mean((/1,4,6/),:)/)
       outmean!0 = "case"
       outmean!1 = "height1"
       outmean&case = (/0,1,2/)
       outmean&height1  = zmean2
       outmean@long_name = "mean profiles"           ; assign attributes 
       outmean@_FillValue = -999.           ; assign attributes 

       ncdf->QT    =  outmean                        
       ncdf->QC    = outmean2                         
       outmean3 = (/flx(1,:)/)
       outmean3!0 = "height2"
       outmean3&height2  = zflx2
       outmean3@long_name = "qt flux profile"           
       outmean3@_FillValue = -999.           ; assign attributes 
       ncdf->QTflx =  outmean3                       
       outmean3 = (/flx(0,:)/)
       outmean3!0 = "height2"
       outmean3&height2  = zflx2
       outmean3@long_name = "theta_l flux profile"           
       outmean3@_FillValue = -999.           ; assign attributes 
       ncdf->THLflx = outmean3                         
       zmean2!0 = "height1"
       ncdf->height1 = (/zmean2/)
       zflx2!0 = "height2"
       ncdf->height2 = (/zflx2/)
