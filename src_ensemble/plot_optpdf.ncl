pmin=(/0.,0.,0.,0.,100.,0./)
pmax=(/4.,4.,2000.,2.,1000.,2./)


colors=(/"blue","red","blue","red"/)


filebase=(/"x_srbf_last21","x_srbf_62_82","x_dy_last21","x_dy_62_82"/)

npara=6
data=new((/dimsizes(filebase),npara,420/),float)
do j=0,dimsizes(filebase)-1
do i=0,npara-1
  filename="./"+filebase(j) +"_"+i+".txt"
  tmp= asciiread(filename,(/21,20/),"float")
  data(j,i,:) = pmin(i) + ndtooned(tmp) * (pmax(i)-pmin(i))
end do
end do

rmsedir="/global/scratch/langhans/metamodel_large2/"
f=addfile(rmsedir+"/RMSE_mmthvdqtdrmse.nc","r")


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


pdf = new((/dimsizes(filebase),npara,(dimsizes(p1)-1)/),float)
pdfcenters = new((/6,(dimsizes(p1)-1)/),float)
optpdf=True

do j=0,dimsizes(filebase)-1
do i=0,5
  optpdf@bin_min=pmin(i)
  optpdf@bin_max=pmax(i)
  zpdf = pdfx(data(j,i,:),(dimsizes(p1)-1),optpdf)
  pdfcenters(i,:)=doubletofloat(zpdf@bin_center)
  zpdf = zpdf / (pdfcenters(i,1)-pdfcenters(i,0))
  pdf(j,i,:) = doubletofloat(zpdf)
  delete(zpdf)
end do
end do

pdf=pdf/100.


plotpdf=new(6,graphic)

  wks = gsn_open_wks("eps","pdf_opti")               ; send graphics to PNG file
;**********************************************
  res=True
  res@gsnDraw     = False
  res@gsnFrame    = False
  res@tiYAxisString=" "
  res@xyLineThicknessF=2.
;  ;res@gsnYRefLine=(/best1,best5/)
  res@xyMarkLineMode="MarkLines"
  res@xyLineColors= colors
  res@xyMonoDashPattern = False
  res@xyDashPatterns=(/0,0,2,2/)
  res@xyMarkerColors=colors
  res@xyMarker = 1
  res@xyMarkerSizeF=0.02
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

  plotpdf(i)  = gsn_csm_xy (wks,pdfcenters(i,:),pdf(:,i,:),res)
  end do

  resp=True
  resp@gsnPanelXWhiteSpacePercent=5.

  gsn_panel(wks,plotpdf,(/2,3/),resp)









