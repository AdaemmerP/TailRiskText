function Output=meanGIG(rho,a,b)



Output.mean= sqrt(b)./sqrt( a).*besselk(rho+1,sqrt(a.*b),1)./besselk(rho,sqrt(a.*b),1);
  Output.meaninv=sqrt(a)./sqrt( b).*besselk(rho+1,sqrt(a.*b),1)./besselk(rho,sqrt(a.*b),1)-2*(rho)./b;

Output.mean2= (sqrt(b)./sqrt( a)).^-1.*besselk(rho-1,sqrt(a.*b),1)./besselk(rho,sqrt(a.*b),1);

