function Output=meanGIG2(rho,a,b,j)




Output= (sqrt(b)./sqrt( a)).^j.*besselk(rho+j,sqrt(a.*b),1)./besselk(rho,sqrt(a.*b),1);