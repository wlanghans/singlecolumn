numfiles=1000
fdir="/global/scratch/langhans/metamodel_mm3/"
nruns=1000
matb=new((/7,nruns*numfiles/),float)
pref=(/2.,2.,1000.,1.,550.,1./)
pmin=(/0.,0.,0.,0.,100.,0./)
pmax=(/4.,4.,2000.,2.,1000.,2./)


prefm=new((/6,nruns/),float)
do i=0,nruns-1
  prefm(:,i)=pref
end do

colors=(/"blue","red","orange"/)
colorsp=(/"blue","red","orange","darkgreen"/)

si=0
do i=0,numfiles-1
  filename=fdir+"RMSE_"+i
  matb(0:6,si:si+nruns-1)=asciiread(filename,(/7,nruns/),"float")
  matb(1:6,si:si+nruns-1)=matb(1:6,si:si+nruns-1) + prefm
  si=si+nruns
end do

opt=False
stat=stat_dispersion(matb(0,:),opt)
print("Lower octile: " + stat(4))
print("Lower dectile: " + stat(3))
print("Lower 5%: " + stat(24))
print("Lower 1%: " + stat(23))
best1=stat(23)
best5=stat(24)
best10=stat(3)
best01=stat(22)
loct = stat(23)

num_min=40
p_min=new((/num_min,6/),float)

val_min=new((/num_min/),float)
ind_max=new((/num_min,6/),integer)
val_max=new((/num_min/),float)

;get six minima
tmp=matb(0,:)
tmp@_FillValue=-999.

do i=0,num_min-1

ind1D=ind(tmp.eq.min(tmp))
if (dimsizes(ind1D).gt.1) then
 print(i +" More than one min found, stopping")
 tmpi=ind1D(0)
 delete(ind1D)
 ind1D=tmpi
 delete(tmpi)
 ;exit
end if

p_min(i,0)=matb(1,ind1D)
p_min(i,1)=matb(2,ind1D)
p_min(i,2)=matb(3,ind1D)
p_min(i,3)=matb(4,ind1D)
p_min(i,4)=matb(5,ind1D)
p_min(i,5)=matb(6,ind1D)

val_min(i)=matb(0,ind1D)
tmp(ind1D) = -999.
delete(ind1D)

end do

do i=0,num_min-1
   print("RMSE min "+ sprinti("%i",i+1) + " RMSE="+val_min(i)+" for a_thv="+p_min(i,0)+ " a_qt="+p_min(i,1)+" ctketau="+p_min(i,2)  + " p_min="+p_min(i,3) + " tauneggers="+p_min(i,4) + " Wc="+p_min(i,5) )
end do

rmsedir="/global/scratch/langhans/metamodel_large3/"
f=addfile(rmsedir+"/RMSE_mmthvdqtdrmse4_paperreview.nc","r")


p1=f->alphathv
p2=f->alphaqt
p3=f->ctketau
p4=f->pwmin
p5=f->tauneggers
p6=f->Wc

pnames=(/"a!BF8!q!F!v!N!","a!B!qt!N!","!F8!t!FB!l!NF!","w!B!min!N!","!F8!t!F!!BF8!e!NF!","b!BF8!e!NF!"/)

p=new((/6,dimsizes(p1)/),float)
do i=0,5
  p(i,:) = fspan(pmin(i),pmax(i),dimsizes(p1))
end do


matc=matb
do i=1,6
matc(i,:)=where(matb(0,:).lt.best01,matb(i,:),-999.)
end do
matc@_FillValue=-999.

pdf = new((/6,dimsizes(p1)-1/),float)
pdfcenters = new((/6,dimsizes(p1)-1/),float)
optpdf=True

do i=0,5
  optpdf@bin_min=pmin(i)
  optpdf@bin_max=pmax(i)
  zpdf = pdfx(matc(i+1,:),dimsizes(p1)-1,optpdf)
  pdfcenters(i,:)=doubletofloat(zpdf@bin_center)
  zpdf = zpdf / (pdfcenters(i,1)-pdfcenters(i,0))
  pdf(i,:) = doubletofloat(zpdf)
  delete(zpdf)
end do

pdf=pdf/100.


plotpdf=new(6,graphic)

  wks = gsn_open_wks("eps","pdf01%_mmstoch")               ; send graphics to PNG file
;**********************************************
  res=True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@tiYAxisString=" "
  res@xyLineThicknessF=3.
;  ;res@gsnYRefLine=(/best1,best5/)
  res@xyMarkLineMode="MarkLines"
  res@xyLineColor= colors(1)
  res@xyMarkerColor=colors(1)
  res@xyMarker = 1
  res@xyMarkerSizeF=0.04
  res@tiXAxisFontHeightF=0.03
  res@trYMinF = 0.
  res@gsnXRefLineColor="grey40"
  res@gsnXRefLineThicknessF=4.

  opt1=(/0.2,0.6,0.0,0.21,449.,0.55/)
  opt2=(/0.6,0.88,1187.,0.72,407.,0.48/)

  ymaxvec=(/2.,1.2,0.002,2.,0.0085,4.0/)
  do i=0,5
;  res@gsnXRefLine=(/opt1(i),opt2(i)/)
  res@tiXAxisString=pnames(i)
  res@tiYAxisString="f("+pnames(i)+")"+"/!F8!D!F!"+pnames(i)
  res@tmXBLabels =  p(i,:) ; labels for each box
  res@trXMinF =p(i,0)
  res@trXMaxF = p(i,dimsizes(p1)-1)
  res@trYMaxF = ymaxvec(i)  ;1./(2.*pdfcenters(i,0))

  plotpdf(i)  = gsn_csm_xy (wks,pdfcenters(i,:),pdf(i,:),res)
  end do

  resp=True
  resp@gsnPanelXWhiteSpacePercent=5.

  gsn_panel(wks,plotpdf,(/2,3/),resp)









